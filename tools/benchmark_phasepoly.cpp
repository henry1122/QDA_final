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
    pp(k)       PhasePoly A*    warm-start group synthesis (k=group size)
                k=1: per-block independent; k>1: cross-block parity reuse
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
#include <filesystem>
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
// helpers
// ---------------------------------------------------------------------------

static size_t count_cx(QCir const& qcir) {
    size_t n = 0;
    for (auto* g : qcir.get_gates())
        if (g->get_operation() == CXGate()) ++n;
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
    size_t cx_pp = 0;
    if (m == 0) {
        cx_pp = pmh_o;
    } else if (m > max_rz) {
        cx_pp = FAIL;  // too large; marked as n/a
    } else {
        cx_pp = synthesize_phasepoly(problem, cfg).num_cx;
    }

    if (m == 0) {
        return {m, cx_pp, pmh_o, pmh_o, pmh_o, pmh_o};
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

    return {m, cx_pp, cx_mst, cx_gstair, cx_gray, cx_naive};
}

// ---------------------------------------------------------------------------
// per-circuit result
// ---------------------------------------------------------------------------

struct CircuitResult {
    std::string name;
    size_t n_qubits        = 0;
    size_t n_blocks        = 0;
    size_t total_rz        = 0;
    size_t cx_pp_k1        = 0;   // PhasePoly A* k=1 (per-block independent)
    size_t cx_pp_k2        = 0;   // PhasePoly A* k=2 warm-start groups
    size_t cx_pp_k3        = 0;   // PhasePoly A* k=3 warm-start groups
    size_t cx_pp_k5        = 0;   // PhasePoly A* k=5 warm-start groups
    size_t cx_mst          = 0;
    size_t cx_gstair       = 0;
    size_t cx_gray         = 0;
    size_t cx_naive        = 0;
    size_t cx_todd_naive   = FAIL;  // full-circuit Todd+Naive
    bool   pp_has_skipped  = false; // any block exceeded max_rz
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
    auto const extracted = extract_phase_blocks_with_boundaries(qcir);
    r.n_blocks = extracted.blocks.size();

    // k=1: per-block independent synthesis (same as before)
    for (auto const& block : extracted.blocks) {
        auto const problem = phase_block_to_problem(block);
        auto const br      = benchmark_block(problem, cfg, max_rz);
        r.total_rz    += br.rz;
        // When PhasePoly A* is skipped (block > max_rz), fall back to MST so
        // the pp total is still a fair whole-circuit number, not artificially low.
        size_t pp_this = (br.cx_phasepoly == FAIL) ? br.cx_mst : br.cx_phasepoly;
        r.cx_pp_k1    += (pp_this == FAIL) ? 0 : pp_this;
        r.cx_mst      += (br.cx_mst    == FAIL) ? 0 : br.cx_mst;
        r.cx_gstair   += (br.cx_gstair == FAIL) ? 0 : br.cx_gstair;
        r.cx_gray     += (br.cx_gray   == FAIL) ? 0 : br.cx_gray;
        r.cx_naive    += (br.cx_naive  == FAIL) ? 0 : br.cx_naive;
        if (br.cx_phasepoly == FAIL) r.pp_has_skipped = true;
    }

    // k=2,3,5: warm-start group synthesis
    r.cx_pp_k2 = synthesize_grouped(extracted, 2, cfg);
    r.cx_pp_k3 = synthesize_grouped(extracted, 3, cfg);
    r.cx_pp_k5 = synthesize_grouped(extracted, 5, cfg);

    // Full-circuit Todd+Naive comparison
    if (run_todd) {
        auto tableau_opt = to_tableau(qcir);
        if (tableau_opt) {
            optimize_phase_polynomial(*tableau_opt, ToddPhasePolynomialOptimizationStrategy{});
            auto todd_qcir = to_qcir(*tableau_opt,
                                     HOptSynthesisStrategy{},
                                     NaivePauliRotationsSynthesisStrategy{});
            if (todd_qcir) r.cx_todd_naive = count_cx(*todd_qcir);
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
        fmt::println("Options:");
        fmt::println("  --max-rz N    Skip PhasePoly A* for blocks with > N Rz terms (default: 30)");
        fmt::println("  --max-queue N A* open-set cap (default: 5000)");
        fmt::println("  --max-exp N   A* expansion cap (default: 100000)");
        fmt::println("  --no-todd     Skip full-circuit Todd comparison");
        return 1;
    }

    PhasePolyConfig cfg;
    cfg.max_queue_size = 5000;
    cfg.max_expansions = 100000;
    cfg.max_solutions  = 5;

    size_t max_rz  = 30;
    bool   no_todd = false;

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

    fmt::println("Config: max_rz={} max_queue={} max_exp={} todd={}",
                 max_rz, cfg.max_queue_size, cfg.max_expansions,
                 no_todd ? "off" : "on");
    fmt::println("");

    std::vector<CircuitResult> results;
    results.reserve(circuit_paths.size());
    for (auto const& p : circuit_paths) {
        fmt::print("  {:.<50}", p.stem().string());
        fflush(stdout);
        auto r = benchmark_circuit(p, cfg, max_rz, !no_todd);
        if (r.ok) {
            std::string tag = r.pp_has_skipped ? " [pp:partial]" : "";
            fmt::println(" done  (blocks={}, Rz={}{})", r.n_blocks, r.total_rz, tag);
        } else {
            fmt::println(" FAILED to read");
        }
        results.push_back(std::move(r));
    }

    fmt::println("");
    // ---- Header ----
    // Columns: circuit | q | blk | Rz | pp(1) | pp(2) | pp(3) | pp(5) | mst+P | gstair+P | gray+P | naive+P | todd+naive
    // k = group size for warm-start multi-block synthesis
    static constexpr int NW = 38;
    fmt::println("{:<{}} {:>3} {:>4} {:>5}  {:>7}  {:>6}  {:>6}  {:>6}  {:>7}  {:>8}  {:>7}  {:>8}  {:>10}",
                 "circuit", NW, "q", "blk", "Rz",
                 "pp(1)", "pp(2)", "pp(3)", "pp(5)",
                 "mst+P", "gstair+P", "gray+P", "naive+P", "todd+naive");
    int const LINE = NW + 3 + 4 + 5 + 7 + 6 + 6 + 6 + 7 + 8 + 7 + 8 + 10 + 13 * 2;
    fmt::println("{}", std::string(LINE, '-'));

    size_t tot_pp_k1 = 0, tot_pp_k2 = 0, tot_pp_k3 = 0, tot_pp_k5 = 0;
    size_t tot_mst = 0, tot_gstair = 0, tot_gray = 0, tot_naive = 0, tot_todd = 0;
    size_t n_todd = 0;

    for (auto const& r : results) {
        if (!r.ok) {
            fmt::println("{:<{}}  (read failed)", r.name, NW);
            continue;
        }

        std::string todd_str = (r.cx_todd_naive == FAIL) ? "       n/a" :
                               fmt::format("{:>10}", r.cx_todd_naive);
        std::string pp1_str  = r.pp_has_skipped ?
                               fmt::format("{:>7}*", r.cx_pp_k1) :
                               fmt::format("{:>7}", r.cx_pp_k1);

        fmt::println("{:<{}} {:>3} {:>4} {:>5}  {}  {:>6}  {:>6}  {:>6}  {:>7}  {:>8}  {:>7}  {:>8}  {}",
                     r.name, NW, r.n_qubits, r.n_blocks, r.total_rz,
                     pp1_str, r.cx_pp_k2, r.cx_pp_k3, r.cx_pp_k5,
                     r.cx_mst, r.cx_gstair, r.cx_gray, r.cx_naive,
                     todd_str);

        tot_pp_k1  += r.cx_pp_k1;
        tot_pp_k2  += r.cx_pp_k2;
        tot_pp_k3  += r.cx_pp_k3;
        tot_pp_k5  += r.cx_pp_k5;
        tot_mst    += r.cx_mst;
        tot_gstair += r.cx_gstair;
        tot_gray   += r.cx_gray;
        tot_naive  += r.cx_naive;
        if (r.cx_todd_naive != FAIL) { tot_todd += r.cx_todd_naive; ++n_todd; }
    }

    fmt::println("{}", std::string(LINE, '-'));
    std::string todd_total = (n_todd == 0) ? "       n/a" :
                             fmt::format("{:>10}", tot_todd);
    fmt::println("{:<{}} {:>3} {:>4} {:>5}  {:>7}  {:>6}  {:>6}  {:>6}  {:>7}  {:>8}  {:>7}  {:>8}  {}",
                 "TOTAL", NW, "-", "-", "-",
                 tot_pp_k1, tot_pp_k2, tot_pp_k3, tot_pp_k5,
                 tot_mst, tot_gstair, tot_gray, tot_naive,
                 todd_total);

    fmt::println("");
    fmt::println("Notes:");
    fmt::println("  pp(k) = PhasePoly A* with warm-start group size k (k=1: per-block independent)");
    fmt::println("  mst+P / gstair+P / gray+P / naive+P = block-level synthesis + PMH output-matrix pass");
    fmt::println("  todd+naive = full-circuit Todd T-count opt (Tableau pipeline) + naive rotation synthesis");
    if (n_todd > 0 && n_todd < results.size())
        fmt::println("  todd+naive: {} / {} circuits succeeded", n_todd, results.size());
    fmt::println("  * = PhasePoly A* skipped for blocks with > {} Rz terms; MST used as fallback", max_rz);

    return 0;
}
