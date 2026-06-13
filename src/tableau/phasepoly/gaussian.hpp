/****************************************************************************
  PackageName  [ tableau/phasepoly ]
  Synopsis     [ CNOT-only linear-reversible synthesis (output-matrix finish) ]
  Author       [ Design Verification Lab ]
  Copyright    [ Copyright(c) 2025 DVLab, GIEE, NTU, Taiwan ]
****************************************************************************/

#pragma once

#include <cstddef>
#include <vector>

#include "./parity_matrix.hpp"

namespace qsyn::experimental::phasepoly {

/**
 * @brief A CNOT gate emitted by linear-reversible synthesis / the search.
 *
 * Indices follow the locked PhasePoly convention: a matrix row operation
 * `apply_cnot(control, target)` is emitted as the physical `CNOT(control,
 * target)` with the same indices, in the same order.
 */
struct CnotOp {
    size_t control;
    size_t target;
    bool operator==(CnotOp const&) const = default;
};

/**
 * @brief Strategy for reducing the output matrix `O` to identity.
 *   - gauss_jordan:       convention-native column-by-column elimination with a
 *                         minimum-weight pivot seed. Simple and self-contained.
 *   - patel_markov_hayes: block elimination via `dvlab::BooleanMatrix`
 *                         (`gaussian_elimination_skip`), searching block sizes
 *                         for the fewest CNOTs. Default -- fewer gates.
 */
enum class LinearSynthesisMode {
    gauss_jordan,
    patel_markov_hayes,
};

/**
 * @brief Synthesize a CNOT sequence that reduces the square, invertible
 *        `matrix` to the identity (the PhasePoly "Gaussian finish", paper §3.2).
 *
 * Emitting the returned ops as physical `CNOT(control, target)` in order
 * reconstructs the linear map `matrix` encodes. Precondition: `matrix` is
 * square and invertible over GF(2) (always true for an output matrix derived
 * from a CNOT network).
 */
std::vector<CnotOp> synthesize_linear_reversible(
    ParityMatrix const& matrix,
    LinearSynthesisMode mode = LinearSynthesisMode::patel_markov_hayes);

/**
 * @brief The CNOT count `synthesize_linear_reversible` would emit -- the `h2(n)`
 *        output-parity cost estimate of the A* heuristic (paper §3.2.2).
 */
size_t linear_reversible_cnot_cost(
    ParityMatrix const& matrix,
    LinearSynthesisMode mode = LinearSynthesisMode::patel_markov_hayes);

}  // namespace qsyn::experimental::phasepoly
