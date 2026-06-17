/****************************************************************************
  PackageName  [ tableau/phasepoly ]
  Synopsis     [ Multi-block SSA merging for QCir-level optimization (§3.3) ]
  Author       [ Design Verification Lab ]
  Copyright    [ Copyright(c) 2025 DVLab, GIEE, NTU, Taiwan ]
****************************************************************************/

#pragma once

#include <cstddef>
#include <optional>
#include <utility>
#include <vector>

#include "./config.hpp"
#include "./phase_block.hpp"
#include "qcir/operation.hpp"
#include "qsyn/qsyn_type.hpp"

namespace qsyn::qcir {
class QCir;
}

namespace qsyn::experimental::phasepoly {

struct PhasePolyOptimizeStats;

struct BoundaryGate {
    qcir::Operation op;
    QubitIdList qubits;
};

/// @brief A phase-polynomial block plus the non-phase gates that follow it.
struct PhasePolySegment {
    PhaseBlock block;
    std::vector<BoundaryGate> trailing;
};

/**
 * @brief Split a circuit into maximal phase-polynomial blocks separated by
 *        non-phase gates (e.g. H, Toffoli).
 */
std::vector<PhasePolySegment> segment_phase_poly_regions(qcir::QCir const& circuit);

/**
 * @brief Merge `k` consecutive segments (blocks + SSA boundaries) into one
 *        lifted `PhaseBlock` when every interstitial gate is an H (SSA fork).
 *        Returns `nullopt` when the window is invalid or contains non-SSA gates.
 */
std::optional<PhaseBlock> merge_segments_ssa(
    std::vector<PhasePolySegment> const& segments,
    size_t start,
    size_t k);

/**
 * @brief Map synthesized SSA-row qubit indices back to logical qubit ids.
 */
std::vector<size_t> ssa_row_to_logical_map(
    std::vector<PhasePolySegment> const& segments,
    size_t start,
    size_t k);

/**
 * @brief Optimize a circuit with progressive multi-block grouping (k = 1,2,3,5).
 */
std::optional<qcir::QCir> optimize_qcir_phasepoly_multiblock(
    qcir::QCir const& circuit,
    PhasePolyConfig const& config           = {},
    bool replace_only_if_better             = true,
    PhasePolyOptimizeStats* stats           = nullptr,
    size_t* phase_gate_out                  = nullptr);

}  // namespace qsyn::experimental::phasepoly
