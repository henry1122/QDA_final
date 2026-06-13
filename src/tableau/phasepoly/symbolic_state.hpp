/****************************************************************************
  PackageName  [ tableau/phasepoly ]
  Synopsis     [ Forward symbolic simulator for phase-polynomial parities ]
  Author       [ Design Verification Lab ]
  Copyright    [ Copyright(c) 2025 DVLab, GIEE, NTU, Taiwan ]
****************************************************************************/

#pragma once

#include <sul/dynamic_bitset.hpp>

#include <cstddef>
#include <vector>

#include "./parity_matrix.hpp"

namespace qsyn::experimental::phasepoly {

/**
 * @brief Tracks, for each qubit, the parity over the input variables that the
 *        qubit currently holds while scanning a CNOT network forward.
 *
 * This implements the *physical* circuit semantics and is used during block
 * extraction (paper §"Single-block extraction") and as the ground truth in the
 * convention unit tests:
 *   - Initially qubit `q` holds input variable `q`   (value[q] = e_q).
 *   - A CNOT with control `c` and target `t` performs `value[t] ^= value[c]`,
 *     i.e. the target qubit becomes the XOR of control and target.
 *   - An `Rz` on qubit `q` is associated with the parity `value[q]`.
 *
 * Contrast with `ParityMatrix::apply_cnot`, which XORs the *control* row: the
 * two are dual representations. The output basis matrix produced by
 * `to_output_matrix()` satisfies `O = T^T`, where `T` is the linear transform
 * realised by the scanned CNOT network; reducing `O` to identity with
 * `ParityMatrix::apply_cnot(c, t)` operations and emitting a physical
 * `CNOT(control = c, target = t)` for each, in the same order, reconstructs the
 * map. This round trip is asserted in the Stage-0 tests.
 */
class SymbolicState {
public:
    explicit SymbolicState(size_t n_qubits);

    size_t n_qubits() const { return _values.size(); }

    /// @brief The parity (over input variables) currently held by `qubit`.
    sul::dynamic_bitset<> const& parity_of(size_t qubit) const { return _values[qubit]; }

    /// @brief Physical CNOT: value[target] <- value[target] XOR value[control].
    void apply_cnot(size_t control, size_t target);

    /// @brief Build the output parity matrix `O`: column `q` is `parity_of(q)`.
    ParityMatrix to_output_matrix() const;

private:
    std::vector<sul::dynamic_bitset<>> _values;  // value[q] over input variables
};

/**
 * @brief Free-function alias matching the recommended helper naming in the
 *        PhasePoly spec. Applies the physical CNOT semantics
 *        `value[target] <- value[target] XOR value[control]`.
 */
inline void apply_cnot_to_state(SymbolicState& state, size_t control, size_t target) {
    state.apply_cnot(control, target);
}

}  // namespace qsyn::experimental::phasepoly
