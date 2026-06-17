/****************************************************************************
  PackageName  [ tableau/phasepoly ]
  Synopsis     [ Block-level Todd phase-polynomial preprocessing ]
  Author       [ Design Verification Lab ]
  Copyright    [ Copyright(c) 2025 DVLab, GIEE, NTU, Taiwan ]
****************************************************************************/

#pragma once

#include <optional>
#include <vector>

#include "./multiblock.hpp"
#include "./phase_block.hpp"
#include "./phase_poly_problem.hpp"

namespace qsyn::experimental::phasepoly {

/// @brief Todd-reduce a `(P, Theta, O)` problem; nullopt if no improvement.
std::optional<PhasePolyProblem> todd_optimize_problem(PhasePolyProblem const& problem);

/// @brief Todd-optimize every non-empty block inside `segments` (in place).
void todd_preprocess_segments(std::vector<PhasePolySegment>& segments);

/// @brief Count Rz terms after Todd preprocessing (for benchmark stats).
size_t count_rz_after_todd(std::vector<PhasePolySegment> const& segments);

}  // namespace qsyn::experimental::phasepoly
