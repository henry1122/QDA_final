/****************************************************************************
  PackageName  [ tableau/phasepoly ]
  Synopsis     [ Forward symbolic simulator for phase-polynomial parities ]
  Author       [ Design Verification Lab ]
  Copyright    [ Copyright(c) 2025 DVLab, GIEE, NTU, Taiwan ]
****************************************************************************/

#include "./symbolic_state.hpp"

#include <cassert>

namespace qsyn::experimental::phasepoly {

SymbolicState::SymbolicState(size_t n_qubits)
    : _values(n_qubits, sul::dynamic_bitset<>(n_qubits)) {
    for (size_t q = 0; q < n_qubits; ++q) {
        _values[q].set(q);  // qubit q initially holds input variable q
    }
}

void SymbolicState::apply_cnot(size_t control, size_t target) {
    assert(control < _values.size());
    assert(target < _values.size());
    if (control == target) return;
    _values[target] ^= _values[control];
}

ParityMatrix SymbolicState::to_output_matrix() const {
    size_t const n = _values.size();
    ParityMatrix output(n, n);
    for (size_t q = 0; q < n; ++q) {
        for (size_t r = 0; r < n; ++r) {
            if (_values[q].test(r)) output.set(r, q, true);
        }
    }
    return output;
}

}  // namespace qsyn::experimental::phasepoly
