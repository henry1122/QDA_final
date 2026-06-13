/****************************************************************************
  PackageName  [ tableau/phasepoly ]
  Synopsis     [ Memory-bounded A* search over the joint [P | O] matrix ]
  Author       [ Design Verification Lab ]
  Copyright    [ Copyright(c) 2025 DVLab, GIEE, NTU, Taiwan ]
****************************************************************************/

#pragma once

#include <cstddef>
#include <vector>

#include "./config.hpp"
#include "./phase_block.hpp"
#include "./phase_poly_problem.hpp"

namespace qsyn::experimental::phasepoly {

/**
 * @brief A synthesized phase-polynomial block plus search statistics.
 *
 * `gates` is the ordered, interleaved list of emitted CNOTs and `Rz`s. Emitting
 * them as a circuit (see `build_qcir`) realizes the same phase polynomial and
 * output basis as the input problem.
 */
struct SynthesisResult {
    size_t n_qubits = 0;
    std::vector<PhaseOp> gates;

    size_t num_cx = 0;
    size_t num_rz = 0;

    // statistics
    size_t expansions    = 0;  ///< number of states expanded
    size_t max_queue     = 0;  ///< peak open-set size
    size_t num_solutions = 0;  ///< goal states collected
    bool used_fallback   = false;

    size_t total_gates() const { return num_cx + num_rz; }
};

/**
 * @brief Co-optimize a phase-polynomial block with a memory-bounded A* search.
 *
 * Explores CNOT row operations on the joint `[P | O]` matrix (paper §3.1-3.2):
 * each CNOT reduces phase-parity columns (emitting an `Rz` when a column hits
 * weight 1) while transforming `O` toward the identity. The cost function is
 * `f = g + h1 + h2`; ties break on `(f, h1, h2, -g)`. Returns the cheapest of
 * up to `max_solutions` goal states; falls back to a guaranteed greedy
 * synthesis if the search budget is exhausted before any solution is found.
 */
SynthesisResult synthesize_phasepoly(PhasePolyProblem const& problem,
                                     PhasePolyConfig const& config = {});

}  // namespace qsyn::experimental::phasepoly
