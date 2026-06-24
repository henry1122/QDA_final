/****************************************************************************
  PackageName  [ tableau/phasepoly ]
  Synopsis     [ PhasePoly paper Table 1 benchmark harness ]
  Author       [ Design Verification Lab ]
  Copyright    [ Copyright(c) 2025 DVLab, GIEE, NTU, Taiwan ]
****************************************************************************/

#pragma once

#include <cstddef>
#include <filesystem>
#include <optional>
#include <string>
#include <vector>

#include "./config.hpp"

namespace qsyn::experimental::phasepoly {

/// @brief One row of paper Table 1 (CNOT counts after Todd + synthesis).
struct Table1Row {
    std::string circuit;
    size_t qubits       = 0;
    size_t num_blocks   = 0;
    size_t num_rz       = 0;

    size_t pp           = 0;
    size_t mst          = 0;
    size_t gstair       = 0;
    size_t gray         = 0;
    size_t naive        = 0;
    size_t todd_naive   = 0;

    /// CX+Rz counts inside phase-polynomial blocks only (excludes boundary gates).
    size_t pp_phase     = 0;
    size_t mst_phase    = 0;
    size_t gstair_phase = 0;
    size_t gray_phase   = 0;
    size_t naive_phase  = 0;

    bool ok = false;
    std::string error;
};

struct Table1Summary {
    std::vector<Table1Row> rows;
    size_t total_pp         = 0;
    size_t total_mst        = 0;
    size_t total_gstair     = 0;
    size_t total_gray       = 0;
    size_t total_naive      = 0;
    size_t total_todd_naive = 0;
    size_t total_pp_phase   = 0;
    size_t total_mst_phase  = 0;
    size_t total_gstair_phase = 0;
    size_t total_gray_phase = 0;
    size_t total_naive_phase = 0;
};

struct Table1BenchmarkOptions {
    /// Apply block-level Todd preprocessing before synthesis (paper pipeline).
    bool use_todd = true;
    /// Use multiblock SSA merging for PhasePoly (§3.3). Off = paper Table 1 single-block.
    bool use_multiblock = false;
};

/// @brief Paper Table 1 benchmark circuits (`benchmark/qc/optimized/*_pyzx.qc`).
std::vector<std::string> table1_circuit_names();

/**
 * @brief Run the full Table 1 pipeline on one circuit:
 *   read → tableau → Todd (+P) → synthesize {pp,mst,gstair,gray,naive,todd+naive}.
 */
Table1Row benchmark_table1_circuit(
    std::filesystem::path const& qc_path,
    PhasePolyConfig const& config = {},
    Table1BenchmarkOptions const& options = {});

/// @brief Run all Table 1 circuits and aggregate totals.
Table1Summary benchmark_table1(
    std::vector<std::filesystem::path> const& circuits,
    PhasePolyConfig const& config = {},
    Table1BenchmarkOptions const& options = {});

void print_table1_markdown(Table1Summary const& summary, std::ostream& os);
void print_table1_csv(Table1Summary const& summary, std::ostream& os);
/// @brief CSV with CX+Rz counts for phase-polynomial blocks only (aligned with paper P).
void print_table1_phase_csv(Table1Summary const& summary, std::ostream& os);

}  // namespace qsyn::experimental::phasepoly
