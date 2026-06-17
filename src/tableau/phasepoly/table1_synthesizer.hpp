/****************************************************************************
  PackageName  [ tableau/phasepoly ]
  Synopsis     [ Per-block synthesis for Table 1 (MST / Gray / PhasePoly / ...) ]
  Author       [ Design Verification Lab ]
  Copyright    [ Copyright(c) 2025 DVLab, GIEE, NTU, Taiwan ]
****************************************************************************/

#pragma once

#include <optional>
#include <string>

#include "./config.hpp"
#include "qcir/qcir.hpp"

namespace qsyn::experimental::phasepoly {

enum class Table1Strategy {
    phasepoly,
    mst,
    gray,
    gstair,
    naive,
};

/// @brief clifford + t-family on the full synthesized circuit (includes boundary gates).
size_t table1_gate_count(qcir::QCir const& circuit);

/// @brief CX + Rz/PZ gates inside synthesized phase-polynomial blocks only.
size_t table1_phase_cx_rz_count(qcir::QCir const& circuit);

/// @brief Synthesize a full circuit using the paper Table 1 pipeline on PhaseBlock IR.
/// When `phase_gate_out` is non-null, writes CX+Rz count for phase blocks only.
std::optional<qcir::QCir> synthesize_table1(
    qcir::QCir const& circuit,
    Table1Strategy strategy,
    PhasePolyConfig const& config = {},
    size_t* phase_gate_out        = nullptr);

}  // namespace qsyn::experimental::phasepoly
