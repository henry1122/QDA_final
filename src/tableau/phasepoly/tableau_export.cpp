/****************************************************************************
  PackageName  [ tableau/phasepoly ]
  Synopsis     [ Todd preprocessing and structure-preserving tableau export ]
  Author       [ Design Verification Lab ]
  Copyright    [ Copyright(c) 2025 DVLab, GIEE, NTU, Taiwan ]
****************************************************************************/

#include "./tableau_export.hpp"

#include <optional>
#include <vector>

#include "./symbolic_state.hpp"
#include "convert/qcir_to_tableau.hpp"
#include "convert/tableau_to_qcir.hpp"
#include "qcir/basic_gate_type.hpp"
#include "tableau/stabilizer_tableau.hpp"
#include "tableau/tableau_optimization.hpp"

namespace qsyn::experimental::phasepoly {

namespace {

bool clifford_is_cnot_only(StabilizerTableau const& clifford) {
    AGSynthesisStrategy strategy;
    for (auto const& op : extract_clifford_operators(clifford, strategy)) {
        if (op.first != CliffordOperatorType::cx) return false;
    }
    return true;
}

void replay_cnots_to_state(StabilizerTableau const& clifford, SymbolicState& state) {
    AGSynthesisStrategy strategy;
    for (auto const& op : extract_clifford_operators(clifford, strategy)) {
        if (op.first == CliffordOperatorType::cx && op.second[0] != op.second[1]) {
            state.apply_cnot(op.second[0], op.second[1]);
        }
    }
}

std::optional<size_t> find_parity_qubit(SymbolicState const& state, sul::dynamic_bitset<> const& parity) {
    for (size_t q = 0; q < state.n_qubits(); ++q) {
        if (state.parity_of(q) == parity) return q;
    }
    return std::nullopt;
}

bool append_diagonal_rotations(
    qcir::QCir& circuit,
    std::vector<PauliRotation> const& rotations,
    SymbolicState const& state) {
    for (auto const& rotation : rotations) {
        if (!rotation.is_diagonal()) return false;

        sul::dynamic_bitset<> parity(state.n_qubits());
        for (size_t q = 0; q < state.n_qubits(); ++q) {
            if (rotation.is_z(q)) parity.set(q);
        }
        if (parity.count() == 0) continue;

        auto const qubit = find_parity_qubit(state, parity);
        if (!qubit) return false;
        circuit.append(qcir::PZGate(rotation.phase()), {static_cast<QubitIdType>(*qubit)});
    }
    return true;
}

}  // namespace

std::optional<Tableau> todd_optimize_tableau(qcir::QCir const& circuit) {
    auto tableau = to_tableau(circuit);
    if (!tableau) return std::nullopt;
    optimize_phase_polynomial(*tableau, ToddPhasePolynomialOptimizationStrategy{});
    return tableau;
}

std::optional<qcir::QCir> tableau_to_phase_qcir(Tableau const& tableau) {
    if (tableau.is_empty()) return qcir::QCir{tableau.n_qubits()};

    qcir::QCir result{tableau.n_qubits()};
    AGSynthesisStrategy strategy;

    for (size_t i = 0; i < tableau.size(); ++i) {
        if (auto const* clifford = std::get_if<StabilizerTableau>(&tableau[i])) {
            if (auto clifford_circ = to_qcir(*clifford, strategy)) {
                result.compose(*clifford_circ);
            } else {
                return std::nullopt;
            }

            if (i + 1 < tableau.size()) {
                auto const* rotations = std::get_if<std::vector<PauliRotation>>(&tableau[i + 1]);
                if (!rotations || rotations->empty()) continue;

                if (clifford_is_cnot_only(*clifford)) {
                    SymbolicState state(tableau.n_qubits());
                    replay_cnots_to_state(*clifford, state);
                    if (!append_diagonal_rotations(result, *rotations, state)) {
                        return std::nullopt;
                    }
                } else {
                    NaivePauliRotationsSynthesisStrategy naive;
                    if (auto rot_circ = naive.synthesize(*rotations)) {
                        result.compose(*rot_circ);
                    } else {
                        return std::nullopt;
                    }
                }
                ++i;
            }
            continue;
        }

        if (auto const* rotations = std::get_if<std::vector<PauliRotation>>(&tableau[i])) {
            SymbolicState state(tableau.n_qubits());
            if (!append_diagonal_rotations(result, *rotations, state)) {
                NaivePauliRotationsSynthesisStrategy naive;
                if (auto rot_circ = naive.synthesize(*rotations)) {
                    result.compose(*rot_circ);
                } else {
                    return std::nullopt;
                }
            }
        }
    }

    return result;
}

std::optional<qcir::QCir> todd_preprocess_qcir(qcir::QCir const& circuit) {
    auto const tableau = todd_optimize_tableau(circuit);
    if (!tableau) return std::nullopt;
    return tableau_to_phase_qcir(*tableau);
}

}  // namespace qsyn::experimental::phasepoly
