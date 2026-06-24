/****************************************************************************
  PackageName  [ tableau/phasepoly ]
  Synopsis     [ Multi-block SSA merging and warm-start group synthesis ]
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
#include "./phase_poly_problem.hpp"
#include "qcir/operation.hpp"
#include "qcir/qcir.hpp"
#include "qsyn/qsyn_type.hpp"

namespace qsyn::experimental::phasepoly {

struct PhasePolyOptimizeStats;

// ── SSA-merge API (used by qcir_optimizer and table1_synthesizer) ─────────────

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

// ── Warm-start group synthesis API (used by benchmark_phasepoly) ─────────────

/// Boundary between two consecutive phase blocks.
struct BlockBoundary {
    std::vector<size_t> h_qubits;  ///< qubits with H gate at this boundary
    bool has_non_h_gate = false;   ///< true if any non-H gate appears here (blocks SSA merge)
};

/// Result of extract_phase_blocks_with_boundaries.
struct ExtractedCircuit {
    size_t n_qubits = 0;
    std::vector<PhaseBlock> blocks;
    std::vector<BlockBoundary> boundaries;  ///< size == blocks.size() - 1
};

/// Result of synthesizing a group of k consecutive blocks.
struct GroupResult {
    size_t num_cx       = 0;
    size_t num_rz       = 0;
    bool used_fallback  = false;
};

/**
 * @brief Extract phase blocks and inter-block boundary information.
 *        Like extract_phase_blocks but also records which qubits have H gates
 *        at each boundary (needed for SSA warm-start parity reset).
 */
ExtractedCircuit extract_phase_blocks_with_boundaries(qcir::QCir const& circuit);

/**
 * @brief Synthesize k consecutive blocks.
 *        SSA-merge strategy (default): merges all k blocks into a single problem.
 *        Joint-A* strategy: synthesizes consecutive pairs jointly.
 *        Falls back to per-block independent synthesis when boundaries contain
 *        non-H gates.
 */
GroupResult synthesize_block_group(
    std::vector<PhaseBlock> const& blocks,
    std::vector<BlockBoundary> const& boundaries,
    PhasePolyConfig const& cfg);

/**
 * @brief Partition all blocks into non-overlapping groups of size `group_size`
 *        and synthesize each group. Returns aggregate CX and Rz counts.
 */
GroupResult synthesize_grouped(ExtractedCircuit const& extracted,
                               size_t group_size,
                               PhasePolyConfig const& cfg);

}  // namespace qsyn::experimental::phasepoly
