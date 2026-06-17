/****************************************************************************
  PackageName  [ tools ]
  Synopsis     [ Block-level benchmark: PhasePoly vs all qsyn rotation strategies ]
  Author       [ Design Verification Lab ]
  Copyright    [ Copyright(c) 2025 DVLab, GIEE, NTU, Taiwan ]

  Usage:
    ./build/benchmark-phasepoly [options] <file.qc> ...
    ./build/benchmark-phasepoly [options] benchmark/qc/optimized/*.qc

  Options:
    --max-rz N        Skip PhasePoly A* for blocks with > N Rz terms (default: 30)
    --max-queue N     A* open-set size cap (default: 5000)
    --max-exp N       A* expansion cap (default: 100000)
    --no-todd         Skip full-circuit Todd comparison

  For each .qc circuit the tool runs synthesis methods on the same input:

  Block-level methods (per CNOT+Rz block):
    pp(k)       PhasePoly A*    SSA-merge group synthesis (k=group size, paper §3.3)
                k=1: per-block independent; k>1: SSA-merged group → single A*
    mst+P       MST             Vandaele MST on P, then PMH(O)
    gstair+P    GStair          GraySynth(staircase) on P + PMH(O)
    gray+P      GraySynth       GraySynth(star) on P + PMH(O)
    naive+P     Naive           naive ladder synthesis on P + PMH(O)

  Full-circuit methods (via Tableau pipeline, different level):
    todd+naive  Todd            T-count opt via Todd, then naive rotation synthesis

  "+P" denotes an added PMH pass to synthesize the block's output matrix O.
  Todd operates on the full Tableau (includes H-gate Clifford parts) and is
  included for completeness; it optimizes T-count, not CNOT-count directly.
****************************************************************************/

// The library uses stop_requested() as a cooperative cancellation hook.
bool stop_requested() { return false; }

#include <fmt/core.h>
#include <spdlog/spdlog.h>

#include <algorithm>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <optional>
#include <string>
#include <vector>

#include "convert/qcir_to_tableau.hpp"
#include "convert/tableau_to_qcir.hpp"
#include "qcir/basic_gate_type.hpp"
#include "qcir/qcir.hpp"
#include "qcir/qcir_io.hpp"
#include "tableau/pauli_rotation.hpp"
#include "tableau/phasepoly/config.hpp"
#include "tableau/phasepoly/extractor.hpp"
#include "tableau/phasepoly/gaussian.hpp"
#include "tableau/phasepoly/phase_block.hpp"
#include "tableau/phasepoly/phase_poly_problem.hpp"
#include "tableau/phasepoly/multiblock.hpp"
#include "tableau/phasepoly/search.hpp"
#include "tableau/stabilizer_tableau.hpp"
#include "tableau/tableau_optimization.hpp"

namespace fs = std::filesystem;

using namespace qsyn;
using namespace qsyn::qcir;
using namespace qsyn::experimental;
using namespace qsyn::experimental::phasepoly;

// ---------------------------------------------------------------------------
// Predefined circuit sets
// ---------------------------------------------------------------------------

// Table 1 from Chen et al. 2025 (19 circuits).
static std::vector<std::string> const TABLE1_CIRCUITS = {
    "tof_3_pyzx", "tof_4_pyzx", "tof_5_pyzx", "tof_10_pyzx",
    "barenco_tof_3_pyzx", "barenco_tof_4_pyzx", "barenco_tof_5_pyzx", "barenco_tof_10_pyzx",
    "grover_5_pyzx",
    "ham15-low_pyzx", "ham15-med_pyzx", "ham15-high_pyzx",
    "mod5_4_pyzx", "mod_mult_55_pyzx", "mod_red_21_pyzx",
    "qcla_com_7_pyzx", "qcla_mod_7_pyzx",
    "rc_adder_6_pyzx", "vbe_adder_3_pyzx",
};

// Extended set: Table 1 + 9 more (28 circuits total).
static std::vector<std::string> const EXTENDED_EXTRA_CIRCUITS = {
    "Adder8_pyzx", "adder_8_pyzx",
    "csla_mux_3_original_pyzx", "csum_mux_9_corrected_pyzx",
    "hwb6_pyzx",
    "mod_adder_1024_pyzx",
    "nth_prime6_pyzx", "qcla_adder_10_pyzx",
    "qft_4_pyzx",
};

// ---------------------------------------------------------------------------
// helpers
// ---------------------------------------------------------------------------

static size_t count_cx(QCir const& qcir) {
    size_t n = 0;
    for (auto* g : qcir.get_gates())
        if (g->get_operation() == CXGate()) ++n;
    return n;
}

static size_t count_rz(QCir const& qcir) {
    size_t n = 0;
    for (auto* g : qcir.get_gates()) {
        auto const& op = g->get_operation();
        if (op.get_underlying_if<qcir::PZGate>() ||
            op.get_underlying_if<qcir::RZGate>()) ++n;
    }
    return n;
}

static std::vector<PauliRotation>
problem_to_pauli_rotations(PhasePolyProblem const& problem) {
    size_t const n = problem.n_qubits;
    size_t const m = problem.num_phase_terms();
    std::vector<PauliRotation> rotations;
    rotations.reserve(m);
    for (size_t j = 0; j < m; ++j) {
        std::vector<Pauli> paulis(n, Pauli::i);
        for (size_t q = 0; q < n; ++q)
            if (problem.phase_matrix.get(q, j)) paulis[q] = Pauli::z;
        rotations.emplace_back(paulis.begin(), paulis.end(),
                               problem.phase_angles[j]);
    }
    return rotations;
}

static constexpr size_t FAIL = SIZE_MAX / 2;

// ---------------------------------------------------------------------------
// per-block result
// ---------------------------------------------------------------------------

struct BlockResult {
    size_t rz;
    size_t cx_phasepoly;
    size_t rz_phasepoly;  // Rz from A* result (may differ from input rz after SSA cancellation)
    size_t cx_mst;
    size_t cx_gstair;
    size_t cx_gray;
    size_t cx_naive;
};

static BlockResult benchmark_block(PhasePolyProblem const& problem,
                                   PhasePolyConfig const& cfg,
                                   size_t max_rz) {
    size_t const m    = problem.num_phase_terms();
    size_t const pmh_o = linear_reversible_cnot_cost(
        problem.output_matrix, LinearSynthesisMode::patel_markov_hayes);

    // PhasePoly A* — skip if block is too large (avoids OOM)
    size_t cx_pp = 0, rz_pp = 0;
    if (m == 0) {
        cx_pp = pmh_o;
        rz_pp = 0;
    } else if (m > max_rz) {
        cx_pp = FAIL;
        rz_pp = FAIL;
    } else {
        auto const sr = synthesize_phasepoly(problem, cfg);
        cx_pp = sr.num_cx;
        rz_pp = sr.num_rz;
    }

    if (m == 0) {
        return {m, cx_pp, rz_pp, pmh_o, pmh_o, pmh_o, pmh_o};
    }

    auto const rotations = problem_to_pauli_rotations(problem);

    auto run_strategy = [&](PauliRotationsSynthesisStrategy const& strategy) -> size_t {
        auto const circ = strategy.synthesize(rotations);
        return circ ? count_cx(*circ) + pmh_o : FAIL;
    };

    size_t cx_mst    = run_strategy(MstSynthesisStrategy{});
    size_t cx_gstair = run_strategy(GraySynthPauliRotationsSynthesisStrategy{
        GraySynthPauliRotationsSynthesisStrategy::Mode::staircase});
    size_t cx_gray   = run_strategy(GraySynthPauliRotationsSynthesisStrategy{});
    size_t cx_naive  = run_strategy(NaivePauliRotationsSynthesisStrategy{});

    return {m, cx_pp, rz_pp, cx_mst, cx_gstair, cx_gray, cx_naive};
}

// ---------------------------------------------------------------------------
// per-circuit result
// ---------------------------------------------------------------------------

struct CircuitResult {
    std::string name;
    size_t n_qubits        = 0;
    size_t n_blocks        = 0;
    size_t orig_cx         = 0;   // CX count of the original input circuit
    size_t orig_rz         = 0;   // Rz count of the original input circuit
    size_t input_rz        = 0;   // total Rz in input phase blocks
    // PhasePoly A* — both CX and Rz from synthesis output
    size_t cx_pp_k1        = 0;
    size_t rz_pp_k1        = 0;
    size_t cx_pp_k2        = 0;
    size_t rz_pp_k2        = 0;
    size_t cx_pp_k3        = 0;
    size_t rz_pp_k3        = 0;
    size_t cx_pp_k5        = 0;
    size_t rz_pp_k5        = 0;
    // Reference methods: CX only (Rz = input_rz, same for all)
    size_t cx_mst          = 0;
    size_t cx_gstair       = 0;
    size_t cx_gray         = 0;
    size_t cx_naive        = 0;
    size_t cx_todd_naive   = FAIL;
    size_t rz_todd_naive   = FAIL;
    bool   pp_has_skipped  = false;
    bool   ok              = false;
};

static CircuitResult benchmark_circuit(fs::path const& path,
                                       PhasePolyConfig const& cfg,
                                       size_t max_rz,
                                       bool   run_todd) {
    CircuitResult r;
    r.name = path.stem().string();

    auto qcir_opt = from_qc(path);
    if (!qcir_opt) return r;
    auto& qcir = *qcir_opt;

    r.n_qubits = qcir.get_num_qubits();
    r.orig_cx  = count_cx(qcir);
    r.orig_rz  = count_rz(qcir);
    auto const extracted = extract_phase_blocks_with_boundaries(qcir);
    r.n_blocks = extracted.blocks.size();

    // k=1: per-block independent synthesis
    for (auto const& block : extracted.blocks) {
        auto const problem = phase_block_to_problem(block);
        auto const br      = benchmark_block(problem, cfg, max_rz);
        r.input_rz += br.rz;
        // When PhasePoly A* is skipped (block > max_rz), fall back to MST.
        bool const skipped = (br.cx_phasepoly == FAIL);
        r.cx_pp_k1  += skipped ? (br.cx_mst  == FAIL ? 0 : br.cx_mst)  : br.cx_phasepoly;
        r.rz_pp_k1  += skipped ? br.rz                                  : br.rz_phasepoly;
        r.cx_mst    += (br.cx_mst    == FAIL) ? 0 : br.cx_mst;
        r.cx_gstair += (br.cx_gstair == FAIL) ? 0 : br.cx_gstair;
        r.cx_gray   += (br.cx_gray   == FAIL) ? 0 : br.cx_gray;
        r.cx_naive  += (br.cx_naive  == FAIL) ? 0 : br.cx_naive;
        if (skipped) r.pp_has_skipped = true;
    }

    // k=2,3,5: group synthesis (SSA merge or joint A*)
    auto const gr2 = synthesize_grouped(extracted, 2, cfg);
    r.cx_pp_k2 = gr2.num_cx;  r.rz_pp_k2 = gr2.num_rz;
    auto const gr3 = synthesize_grouped(extracted, 3, cfg);
    r.cx_pp_k3 = gr3.num_cx;  r.rz_pp_k3 = gr3.num_rz;
    auto const gr5 = synthesize_grouped(extracted, 5, cfg);
    r.cx_pp_k5 = gr5.num_cx;  r.rz_pp_k5 = gr5.num_rz;

    // Full-circuit Todd+Naive comparison
    if (run_todd) {
        auto tableau_opt = to_tableau(qcir);
        if (tableau_opt) {
            optimize_phase_polynomial(*tableau_opt, ToddPhasePolynomialOptimizationStrategy{});
            auto todd_qcir = to_qcir(*tableau_opt,
                                     HOptSynthesisStrategy{},
                                     NaivePauliRotationsSynthesisStrategy{});
            if (todd_qcir) {
                r.cx_todd_naive = count_cx(*todd_qcir);
                r.rz_todd_naive = count_rz(*todd_qcir);
            }
        }
    }

    r.ok = true;
    return r;
}

// ---------------------------------------------------------------------------
// main
// ---------------------------------------------------------------------------

int main(int argc, char** argv) {
    spdlog::set_level(spdlog::level::warn);

    if (argc < 2) {
        fmt::println("Usage: {} [options] <file.qc> ...", argv[0]);
        fmt::println("       {} [options] --set {{table1|extended}}", argv[0]);
        fmt::println("Options:");
        fmt::println("  --set NAME          Run a predefined circuit set (table1: 19 circuits, extended: 28)");
        fmt::println("  --benchmark-dir D   Base directory for --set (default: benchmark/qc/optimized)");
        fmt::println("  --output FILE       Write results to FILE (default: results/benchmark_TIMESTAMP.txt)");
        fmt::println("  --max-rz N          Skip PhasePoly A* for blocks with > N Rz terms (default: 30)");
        fmt::println("  --max-queue N       A* open-set cap (default: 5000)");
        fmt::println("  --max-exp N         A* expansion cap (default: 100000)");
        fmt::println("  --no-todd           Skip full-circuit Todd comparison");
        fmt::println("  --joint-astar       Use joint A* on pairs (default: SSA merge, paper §3.3)");
        return 1;
    }

    PhasePolyConfig cfg;
    cfg.max_queue_size = 5000;
    cfg.max_expansions = 100000;
    cfg.max_solutions  = 5;

    size_t max_rz  = 30;
    bool   no_todd = false;
    fs::path bench_dir  = "benchmark/qc/optimized";
    fs::path output_path;  // empty = auto-generate

    std::vector<fs::path> circuit_paths;

    for (int i = 1; i < argc; ++i) {
        std::string arg{argv[i]};
        if (arg == "--max-rz" && i + 1 < argc) {
            max_rz = std::stoul(argv[++i]);
        } else if (arg == "--max-queue" && i + 1 < argc) {
            cfg.max_queue_size = std::stoul(argv[++i]);
        } else if (arg == "--max-exp" && i + 1 < argc) {
            cfg.max_expansions = std::stoul(argv[++i]);
        } else if (arg == "--no-todd") {
            no_todd = true;
        } else if (arg == "--joint-astar") {
            cfg.multi_block_strategy = MultiBlockStrategy::joint_astar;
        } else if (arg == "--benchmark-dir" && i + 1 < argc) {
            bench_dir = argv[++i];
        } else if (arg == "--output" && i + 1 < argc) {
            output_path = argv[++i];
        } else if (arg == "--set" && i + 1 < argc) {
            std::string set_name{argv[++i]};
            auto add_set = [&](std::vector<std::string> const& names) {
                for (auto const& n : names)
                    circuit_paths.push_back(bench_dir / (n + ".qc"));
            };
            if (set_name == "table1") {
                add_set(TABLE1_CIRCUITS);
            } else if (set_name == "extended") {
                add_set(TABLE1_CIRCUITS);
                add_set(EXTENDED_EXTRA_CIRCUITS);
            } else {
                fmt::println(stderr, "Unknown set '{}'. Use: table1, extended", set_name);
                return 1;
            }
        } else if (arg.rfind("--", 0) == 0) {
            fmt::println(stderr, "Unknown option: {}", arg);
            return 1;
        } else {
            circuit_paths.emplace_back(arg);
        }
    }

    if (circuit_paths.empty()) {
        fmt::println(stderr, "No circuit files specified.");
        return 1;
    }

    // ── Auto-generate output path if not specified ────────────────────────────
    if (output_path.empty()) {
        std::time_t now = std::time(nullptr);
        char ts[32];
        std::strftime(ts, sizeof(ts), "%Y%m%d_%H%M%S", std::localtime(&now));
        fs::create_directories("results");
        output_path = fs::path("results") / fmt::format("benchmark_{}.txt", ts);
    } else {
        fs::create_directories(output_path.parent_path().empty() ? fs::path(".") : output_path.parent_path());
    }
    std::ofstream outfile{output_path};
    if (!outfile) {
        fmt::println(stderr, "Warning: cannot open output file '{}'", output_path.string());
    }

    // Helper: print to stdout and file simultaneously.
    auto emit = [&](std::string const& line) {
        fmt::println("{}", line);
        if (outfile) outfile << line << '\n';
    };

    std::string const strategy_str =
        cfg.multi_block_strategy == MultiBlockStrategy::joint_astar ? "joint-astar" : "ssa-merge";
    emit(fmt::format("Config: max_rz={} max_queue={} max_exp={} todd={} strategy={}",
                     max_rz, cfg.max_queue_size, cfg.max_expansions,
                     no_todd ? "off" : "on", strategy_str));
    emit("");

    std::vector<CircuitResult> results;
    results.reserve(circuit_paths.size());
    for (auto const& p : circuit_paths) {
        fmt::print("  {:.<50}", p.stem().string());
        fflush(stdout);
        auto r = benchmark_circuit(p, cfg, max_rz, !no_todd);
        if (r.ok) {
            std::string tag = r.pp_has_skipped ? " [pp:partial]" : "";
            fmt::println(" done  (blocks={}, Rz={}{})", r.n_blocks, r.input_rz, tag);
        } else {
            fmt::println(" FAILED to read");
        }
        results.push_back(std::move(r));
    }
    fmt::println("");

    // Columns: circuit | q | blk | original | pp(1) | pp(2) | pp(3) | pp(5) | mst+P | gray+P | todd+naive
    // Each cell shows "CX/total" where total = CX + Rz.
    static constexpr int NW = 28;  // circuit name width
    static constexpr int CW = 10;  // data column width

    std::string const sub    = fmt::format("{:>{}}", "CX/tot", CW);
    std::string const hdr1   = fmt::format("{:<{}}  {:>3} {:>4}  {:>{}}  {:>{}}  {:>{}}  {:>{}}  {:>{}}  {:>{}}  {:>{}}  {:>{}}",
                 "circuit", NW, "q", "blk",
                 "pyzx", CW, "pp(1)", CW, "pp(2)", CW, "pp(3)", CW, "pp(5)", CW,
                 "mst+P", CW, "gray+P", CW, "todd+naive", CW);
    std::string const hdr2   = fmt::format("{:<{}}  {:>3} {:>4}  {}  {}  {}  {}  {}  {}  {}  {}",
                 "", NW, "", "",
                 sub, sub, sub, sub, sub, sub, sub, sub);
    int const LINE = static_cast<int>(hdr1.size());
    std::string const divider(LINE, '-');

    emit(hdr1);
    emit(hdr2);
    emit(divider);

    // accumulators: [0]=CX, [1]=total(CX+Rz)
    size_t tot_orig[2]{}, tot_pp_k1[2]{}, tot_pp_k2[2]{}, tot_pp_k3[2]{}, tot_pp_k5[2]{};
    size_t tot_mst[2]{}, tot_gray[2]{}, tot_todd[2]{};
    size_t n_todd = 0;

    // Format "CX/total" cell; append * if skipped.
    auto cell = [&](size_t cx, size_t rz, bool skipped = false) -> std::string {
        std::string s = fmt::format("{}/{}", cx, cx + rz);
        if (skipped) s += '*';
        return fmt::format("{:>{}}", s, CW);
    };
    auto cell_fail = [&]() -> std::string {
        return fmt::format("{:>{}}", "n/a", CW);
    };
    auto tot_cell = [&](size_t cx, size_t tot) -> std::string {
        return fmt::format("{:>{}}", fmt::format("{}/{}", cx, tot), CW);
    };

    for (auto const& r : results) {
        if (!r.ok) {
            emit(fmt::format("{:<{}}  (read failed)", r.name, NW));
            continue;
        }

        size_t const ref_rz  = r.input_rz;  // reference methods keep input Rz unchanged
        std::string todd_cell = (r.cx_todd_naive == FAIL) ? cell_fail()
                              : cell(r.cx_todd_naive, r.rz_todd_naive);

        emit(fmt::format("{:<{}}  {:>3} {:>4}  {}  {}  {}  {}  {}  {}  {}  {}",
                         r.name, NW, r.n_qubits, r.n_blocks,
                         cell(r.orig_cx,   r.orig_rz),
                         cell(r.cx_pp_k1,  r.rz_pp_k1, r.pp_has_skipped),
                         cell(r.cx_pp_k2,  r.rz_pp_k2),
                         cell(r.cx_pp_k3,  r.rz_pp_k3),
                         cell(r.cx_pp_k5,  r.rz_pp_k5),
                         cell(r.cx_mst,    ref_rz),
                         cell(r.cx_gray,   ref_rz),
                         todd_cell));

        tot_orig[0]  += r.orig_cx;   tot_orig[1]  += r.orig_cx  + r.orig_rz;
        tot_pp_k1[0] += r.cx_pp_k1; tot_pp_k1[1] += r.cx_pp_k1 + r.rz_pp_k1;
        tot_pp_k2[0] += r.cx_pp_k2; tot_pp_k2[1] += r.cx_pp_k2 + r.rz_pp_k2;
        tot_pp_k3[0] += r.cx_pp_k3; tot_pp_k3[1] += r.cx_pp_k3 + r.rz_pp_k3;
        tot_pp_k5[0] += r.cx_pp_k5; tot_pp_k5[1] += r.cx_pp_k5 + r.rz_pp_k5;
        tot_mst[0]   += r.cx_mst;   tot_mst[1]   += r.cx_mst   + ref_rz;
        tot_gray[0]  += r.cx_gray;  tot_gray[1]  += r.cx_gray  + ref_rz;
        if (r.cx_todd_naive != FAIL) {
            tot_todd[0] += r.cx_todd_naive;
            tot_todd[1] += r.cx_todd_naive + (r.rz_todd_naive == FAIL ? 0 : r.rz_todd_naive);
            ++n_todd;
        }
    }

    emit(divider);
    std::string todd_tot = (n_todd == 0) ? cell_fail()
                         : tot_cell(tot_todd[0], tot_todd[1]);
    emit(fmt::format("{:<{}}  {:>3} {:>4}  {}  {}  {}  {}  {}  {}  {}  {}",
                     "TOTAL", NW, "-", "-",
                     tot_cell(tot_orig[0],  tot_orig[1]),
                     tot_cell(tot_pp_k1[0], tot_pp_k1[1]),
                     tot_cell(tot_pp_k2[0], tot_pp_k2[1]),
                     tot_cell(tot_pp_k3[0], tot_pp_k3[1]),
                     tot_cell(tot_pp_k5[0], tot_pp_k5[1]),
                     tot_cell(tot_mst[0],   tot_mst[1]),
                     tot_cell(tot_gray[0],  tot_gray[1]),
                     todd_tot));

    emit("");
    emit("Notes:");
    emit("  Each cell: CX/total where total = CX + Rz (paper Table 1 format).");
    emit("  pyzx      = PyZX-pre-optimized input circuit (before PhasePoly synthesis).");
    emit(fmt::format("  pp(k)     = PhasePoly A* ({}) multi-block synthesis:", strategy_str));
    emit("              k=1: per-block independent A*");
    emit("              k>1: groups of k blocks SSA-merged into one problem, solved by single A*.");
    emit("              Use --joint-astar to switch to the experimental joint A* strategy.");
    emit("  mst+P / gray+P = block-level synthesis + PMH output-matrix pass (Rz = input Rz).");
    emit("  todd+naive = full-circuit Todd T-count opt + naive rotation synthesis.");
    if (n_todd > 0 && n_todd < results.size())
        emit(fmt::format("  todd+naive: {} / {} circuits succeeded.", n_todd, results.size()));
    emit(fmt::format("  * = PhasePoly A* skipped for blocks with > {} Rz terms; MST used.", max_rz));
    emit(fmt::format("Output saved to: {}", output_path.string()));

    return 0;
}
