/****************************************************************************
  PackageName  [ tableau/phasepoly ]
  Synopsis     [ Configuration for the PhasePoly memory-bounded A* search ]
  Author       [ Design Verification Lab ]
  Copyright    [ Copyright(c) 2025 DVLab, GIEE, NTU, Taiwan ]
****************************************************************************/

#pragma once

#include <cstddef>

#include "./gaussian.hpp"

namespace qsyn::experimental::phasepoly {

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
};

}  // namespace qsyn::experimental::phasepoly
