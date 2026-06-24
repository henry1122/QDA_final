/****************************************************************************
  PackageName  [ tableau/phasepoly ]
  Synopsis     [ Configuration for the PhasePoly memory-bounded A* search ]
  Author       [ Design Verification Lab ]
  Copyright    [ Copyright(c) 2025 DVLab, GIEE, NTU, Taiwan ]
****************************************************************************/

#pragma once

#include <cstddef>
#include <limits>
#include <vector>

#include "./gaussian.hpp"

namespace qsyn::experimental::phasepoly {

/// Multi-block synthesis strategy (used by `synthesize_block_group`).
enum class MultiBlockStrategy {
    /// SSA rename → merge all k blocks → single A* (paper §3.3, default).
    ssa_merge,
    /// Joint A* on consecutive pairs — experimental cross-block alternative.
    joint_astar,
};

/**
 * @brief Tunables for `synthesize_phasepoly` (paper §3.2).
 */
struct PhasePolyConfig {
    /// Cap on the open set; when exceeded, the lowest-priority states are dropped.
    size_t max_queue_size = 10000;
    /// Stop after collecting this many goal states; the cheapest is returned.
    size_t max_solutions = 5;
    /// Hard cap on state expansions (a deterministic stand-in for a timeout).
    size_t max_expansions = 500000;

    /// Strategy used to finish the output matrix `O -> I` for a candidate.
    LinearSynthesisMode finish_mode = LinearSynthesisMode::patel_markov_hayes;
    /// Strategy used for the per-state `h2` estimate (kept cheap for speed).
    LinearSynthesisMode h2_mode = LinearSynthesisMode::gauss_jordan;

    /// Adjacent phase-poly block group sizes to try (paper §3.3; Stage 5).
    /// The optimizer keeps the cheapest result across all sizes.
    std::vector<size_t> group_sizes = {1, 2, 3, 5};

    /// Strategy for multi-block synthesis (k > 1 groups).
    MultiBlockStrategy multi_block_strategy = MultiBlockStrategy::ssa_merge;

    // ── Search improvements ───────────────────────────────────────────────────
    /// Sort phase columns before hashing → deduplicate column-permuted states.
    bool canonical_state_key = true;
    /// Prune states where g + h2 ≥ best solution found so far (h2 is admissible).
    bool f_cutoff_prune = true;
    /// Keep only the top-K active pairs ranked by net benefit; SIZE_MAX = all pairs.
    size_t max_candidates = std::numeric_limits<size_t>::max();
    /// Scale max_expansions proportional to merged SSA block size (k > 1 groups).
    bool scale_budget_ssa = false;

    /// Apply Todd phase-polynomial optimization to each phase block before synthesis.
    bool apply_block_todd = true;
};

}  // namespace qsyn::experimental::phasepoly
