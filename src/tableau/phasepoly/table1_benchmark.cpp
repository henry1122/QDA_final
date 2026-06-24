/****************************************************************************
  PackageName  [ tableau/phasepoly ]
  Synopsis     [ PhasePoly paper Table 1 benchmark harness ]
  Author       [ Design Verification Lab ]
  Copyright    [ Copyright(c) 2025 DVLab, GIEE, NTU, Taiwan ]
****************************************************************************/

#include "./table1_benchmark.hpp"

#include <iomanip>
#include <sstream>

#include "./multiblock.hpp"
#include "./table1_synthesizer.hpp"
#include "./todd_preprocess.hpp"
#include "qcir/qcir_io.hpp"

namespace qsyn::experimental::phasepoly {

namespace {

size_t count_rz_segments(std::vector<PhasePolySegment> const& segments) {
    size_t n = 0;
    for (auto const& seg : segments) {
        n += seg.block.num_rz();
    }
    return n;
}

std::optional<size_t> run_strategy(
    qcir::QCir const& circuit,
    Table1Strategy strategy,
    PhasePolyConfig const& config,
    Table1BenchmarkOptions const& options,
    size_t* phase_gate_out = nullptr) {
    size_t phase_gates = 0;
    if (auto out = synthesize_table1(circuit, strategy, config, &phase_gates, options.use_multiblock)) {
        if (phase_gate_out) *phase_gate_out = phase_gates;
        // Paper Table 1: CNOT count only (block synthesis + PMH for O; excludes Rz/H/T).
        return table1_cnot_count(*out);
    }
    return std::nullopt;
}

}  // namespace

std::vector<std::string> table1_circuit_names() {
    return {
        "tof_3_pyzx.qc", "tof_4_pyzx.qc", "tof_5_pyzx.qc", "tof_10_pyzx.qc",
        "barenco_tof_3_pyzx.qc", "barenco_tof_4_pyzx.qc", "barenco_tof_5_pyzx.qc", "barenco_tof_10_pyzx.qc",
        "grover_5_pyzx.qc",
        "ham15-low_pyzx.qc", "ham15-med_pyzx.qc", "ham15-high_pyzx.qc",
        "adder_8_pyzx.qc", "Adder8_pyzx.qc", "Adder16_pyzx.qc", "Adder32_pyzx.qc", "Adder64_pyzx.qc",
        "vbe_adder_3_pyzx.qc", "rc_adder_6_pyzx.qc", "qcla_adder_10_pyzx.qc",
        "mod_red_21_pyzx.qc", "mod_mult_55_pyzx.qc", "mod_adder_1024_pyzx.qc", "mod5_4_pyzx.qc",
        "hwb6_pyzx.qc", "hwb8_pyzx.qc",
    };
}

Table1Row benchmark_table1_circuit(
    std::filesystem::path const& qc_path,
    PhasePolyConfig const& config,
    Table1BenchmarkOptions const& options) {
    Table1Row row;
    row.circuit = qc_path.stem().string();

    auto const circuit = qcir::from_file(qc_path);
    if (!circuit) {
        row.error = "failed to read circuit";
        return row;
    }

    row.qubits = circuit->get_num_qubits();

    PhasePolyConfig run_config = config;
    run_config.apply_block_todd = options.use_todd;

    auto segments = segment_phase_poly_regions(*circuit);
    row.num_blocks = segments.size();
    row.num_rz     = run_config.apply_block_todd ? count_rz_after_todd(segments)
                                                  : count_rz_segments(segments);

    if (auto g = run_strategy(*circuit, Table1Strategy::phasepoly, run_config, options, &row.pp_phase)) {
        row.pp = *g;
    } else {
        row.error = "PhasePoly failed";
        return row;
    }

    if (auto g = run_strategy(*circuit, Table1Strategy::mst, run_config, options, &row.mst_phase)) {
        row.mst = *g;
    } else {
        row.error = "MST failed";
        return row;
    }

    if (auto g = run_strategy(*circuit, Table1Strategy::gray, run_config, options, &row.gray_phase)) {
        row.gray = *g;
    } else {
        row.error = "Gray failed";
        return row;
    }

    if (auto g = run_strategy(*circuit, Table1Strategy::gstair, run_config, options, &row.gstair_phase)) {
        row.gstair = *g;
    } else {
        row.error = "Gstair failed";
        return row;
    }

    if (auto g = run_strategy(*circuit, Table1Strategy::naive, run_config, options, &row.naive_phase)) {
        row.naive = *g;
    } else {
        row.error = "Naive failed";
        return row;
    }

    row.todd_naive = row.naive;
    row.ok         = true;
    return row;
}

Table1Summary benchmark_table1(
    std::vector<std::filesystem::path> const& circuits,
    PhasePolyConfig const& config,
    Table1BenchmarkOptions const& options) {
    Table1Summary summary;
    for (auto const& path : circuits) {
        auto row = benchmark_table1_circuit(path, config, options);
        if (row.ok) {
            summary.total_pp += row.pp;
            summary.total_mst += row.mst;
            summary.total_gstair += row.gstair;
            summary.total_gray += row.gray;
            summary.total_naive += row.naive;
            summary.total_todd_naive += row.todd_naive;
            summary.total_pp_phase += row.pp_phase;
            summary.total_mst_phase += row.mst_phase;
            summary.total_gstair_phase += row.gstair_phase;
            summary.total_gray_phase += row.gray_phase;
            summary.total_naive_phase += row.naive_phase;
        }
        summary.rows.push_back(std::move(row));
    }
    return summary;
}

void print_table1_markdown(Table1Summary const& summary, std::ostream& os) {
    os << "# gate metric: CNOT only (paper Table 1; block synthesis + PMH for O)\n\n";
    os << "| circuit | q | blk | Rz | pp | mst+P | gstair+P | gray+P | naive+P | pp/mst | pp/gray |\n";
    os << "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|\n";

    for (auto const& r : summary.rows) {
        if (!r.ok) {
            os << "| " << r.circuit << " | ERR (" << r.error << ") | | | | | | | | | |\n";
            continue;
        }
        auto const pct = [](size_t a, size_t b) {
            if (b == 0) return std::string("-");
            std::ostringstream ss;
            ss << std::fixed << std::setprecision(1) << (100.0 * static_cast<double>(a) / static_cast<double>(b)) << "%";
            return ss.str();
        };
        os << "| " << r.circuit << " | " << r.qubits << " | " << r.num_blocks << " | " << r.num_rz
           << " | " << r.pp << " | " << r.mst << " | " << r.gstair << " | " << r.gray << " | " << r.naive
           << " | " << pct(r.pp, r.mst) << " | " << pct(r.pp, r.gray) << " |\n";
    }

    os << "\n| **TOTAL** | | | | **" << summary.total_pp << "** | **" << summary.total_mst
       << "** | **" << summary.total_gstair << "** | **" << summary.total_gray << "** | **"
       << summary.total_naive << "** | ";
    if (summary.total_mst) os << (100 * summary.total_pp / summary.total_mst) << "%";
    else os << "-";
    os << " | ";
    if (summary.total_gray) os << (100 * summary.total_pp / summary.total_gray) << "%";
    else os << "-";
    os << " |\n";

    os << "\nPaper reference TOTAL: pp=3869, mst+P=4891 (79%), gray+P=5511 (70%)\n";
}

void print_table1_csv(Table1Summary const& summary, std::ostream& os) {
    os << "# gate metric: CNOT only (paper Table 1; block synthesis + PMH for O)\n";
    os << "circuit,q,blk,Rz,pp,mst,gstair,gray,naive,todd_naive,pp_mst_pct,pp_gray_pct,ok,error\n";
    for (auto const& r : summary.rows) {
        os << r.circuit << ',' << r.qubits << ',' << r.num_blocks << ',' << r.num_rz << ','
           << r.pp << ',' << r.mst << ',' << r.gstair << ',' << r.gray << ',' << r.naive << ','
           << r.todd_naive << ',';
        if (r.ok && r.mst) os << (100 * r.pp / r.mst);
        os << ',';
        if (r.ok && r.gray) os << (100 * r.pp / r.gray);
        os << ',' << (r.ok ? "1" : "0") << ',';
        if (!r.error.empty()) {
            std::string esc = r.error;
            for (auto& c : esc) {
                if (c == ',') c = ';';
            }
            os << esc;
        }
        os << '\n';
    }
}

void print_table1_phase_csv(Table1Summary const& summary, std::ostream& os) {
    os << "# phase-only: CX+Rz inside synthesized phase-polynomial blocks (excludes H/Toffoli/T boundaries)\n";
    os << "circuit,q,blk,Rz,pp,mst,gstair,gray,naive,todd_naive,pp_mst_pct,pp_gray_pct,ok,error\n";
    for (auto const& r : summary.rows) {
        os << r.circuit << ',' << r.qubits << ',' << r.num_blocks << ',' << r.num_rz << ','
           << r.pp_phase << ',' << r.mst_phase << ',' << r.gstair_phase << ',' << r.gray_phase << ','
           << r.naive_phase << ',' << r.naive_phase << ',';
        if (r.ok && r.mst_phase) os << (100 * r.pp_phase / r.mst_phase);
        os << ',';
        if (r.ok && r.gray_phase) os << (100 * r.pp_phase / r.gray_phase);
        os << ',' << (r.ok ? "1" : "0") << ',';
        if (!r.error.empty()) {
            std::string esc = r.error;
            for (auto& c : esc) {
                if (c == ',') c = ';';
            }
            os << esc;
        }
        os << '\n';
    }
}

}  // namespace qsyn::experimental::phasepoly
