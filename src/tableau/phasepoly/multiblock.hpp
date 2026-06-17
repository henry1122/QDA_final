/****************************************************************************
  PackageName  [ tableau/phasepoly ]
  Synopsis     [ Multi-block warm-start group synthesis (paper §3.3) ]
  Author       [ Design Verification Lab ]
  Copyright    [ Copyright(c) 2025 DVLab, GIEE, NTU, Taiwan ]
****************************************************************************/
#pragma once
#include <cstddef>
#include <vector>
#include "./config.hpp"
#include "./phase_block.hpp"
#include "./phase_poly_problem.hpp"
#include "qcir/qcir.hpp"

namespace qsyn::experimental::phasepoly {

/// Boundary between two consecutive phase blocks.
struct BlockBoundary {
    std::vector<size_t> h_qubits;  ///< qubits with H gate at this boundary
    bool has_non_h_gate = false;   ///< true if any non-H gate appears here (blocks SSA merge)
};

/// Result of extract_phase_blocks_with_boundaries
struct ExtractedCircuit {
    size_t n_qubits = 0;
    std::vector<PhaseBlock> blocks;
    std::vector<BlockBoundary> boundaries;  ///< size == blocks.size() - 1
};

/// Result of synthesizing a group of k consecutive blocks
struct GroupResult {
    size_t num_cx  = 0;
    size_t num_rz  = 0;  ///< Rz gates in the synthesis output (may be < input Rz after SSA merging)
    bool used_fallback = false;
};

/**
 * @brief Like extract_phase_blocks but also captures which qubits have H gates
 *        at each inter-block boundary (needed for warm-start parity reset).
 */
ExtractedCircuit extract_phase_blocks_with_boundaries(qcir::QCir const& circuit);

/**
 * @brief Synthesize k consecutive blocks with warm-start parity reuse.
 *
 * Each block starts from the accumulated parity state M left by the previous
 * block's synthesis (paper §3.3). H-gated qubits are reset in M at each
 * boundary. Returns total CX across all k blocks.
 *
 * When k == 1, this is identical to independent per-block synthesis.
 */
GroupResult synthesize_block_group(
    std::vector<PhaseBlock> const& blocks,
    std::vector<BlockBoundary> const& boundaries,
    PhasePolyConfig const& cfg);

/**
 * @brief Partition all blocks into non-overlapping groups of size k and
 *        synthesize each group. Returns total CX and Rz across all groups.
 */
GroupResult synthesize_grouped(ExtractedCircuit const& extracted,
                               size_t group_size,
                               PhasePolyConfig const& cfg);

}  // namespace qsyn::experimental::phasepoly
