/****************************************************************************
  PackageName  [ tableau/phasepoly ]
  Synopsis     [ Extract phase-polynomial blocks and their (P, Theta, O) ]
  Author       [ Design Verification Lab ]
  Copyright    [ Copyright(c) 2025 DVLab, GIEE, NTU, Taiwan ]
****************************************************************************/

#include "./extractor.hpp"

#include <cstddef>
#include <optional>
#include <string>
#include <unordered_map>
#include <utility>

#include "./parity_matrix.hpp"
#include "./symbolic_state.hpp"
#include "qcir/basic_gate_type.hpp"
#include "qcir/operation.hpp"
#include "qcir/qcir.hpp"
#include "qcir/qcir_gate.hpp"
#include "qsyn/qsyn_type.hpp"
#include "util/phase.hpp"

namespace qsyn::experimental::phasepoly {

namespace {

/// @brief A stable key for a parity bit-vector, used to merge identical columns.
std::string bits_key(sul::dynamic_bitset<> const& bits) {
    std::string key(bits.size(), '0');
    for (size_t i = 0; i < bits.size(); ++i) {
        if (bits.test(i)) key[i] = '1';
    }
    return key;
}

/// @brief If `op` is a single-control CNOT, return (control, target) qubits.
std::optional<std::pair<size_t, size_t>>
as_cx(qcir::Operation const& op, QubitIdList const& qubits) {
    auto const control_gate = op.get_underlying_if<qcir::ControlGate>();
    if (!control_gate || control_gate->get_num_ctrls() != 1) return std::nullopt;
    auto const target = control_gate->get_target_operation().get_underlying_if<qcir::PXGate>();
    if (!target || target->get_phase() != dvlab::Phase(1)) return std::nullopt;
    // CXGate pin order is (control, target).
    return std::pair{static_cast<size_t>(qubits[0]), static_cast<size_t>(qubits[1])};
}

/// @brief If `op` is a single-qubit Z-rotation (p / rz family), return its angle.
std::optional<dvlab::Phase> as_rz(qcir::Operation const& op) {
    if (auto const pz = op.get_underlying_if<qcir::PZGate>()) return pz->get_phase();
    if (auto const rz = op.get_underlying_if<qcir::RZGate>()) return rz->get_phase();
    return std::nullopt;
}

}  // namespace

PhasePolyProblem phase_block_to_problem(PhaseBlock const& block) {
    size_t const n = block.n_qubits();
    SymbolicState state(n);

    // Collect phase columns in first-appearance order, merging identical
    // parities by accumulating their angles.
    std::vector<sul::dynamic_bitset<>> columns;
    std::vector<dvlab::Phase> angles;
    std::unordered_map<std::string, size_t> column_of_parity;

    for (auto const& op : block.ops()) {
        if (op.is_cx()) {
            apply_cnot_to_state(state, op.control, op.target);
            continue;
        }
        auto const parity = state.parity_of(op.qubit());
        if (parity.count() == 0) continue;  // an all-zero parity carries no phase
        auto const key = bits_key(parity);
        if (auto const it = column_of_parity.find(key); it != column_of_parity.end()) {
            angles[it->second] += op.phase;
        } else {
            column_of_parity.emplace(key, columns.size());
            columns.push_back(parity);
            angles.push_back(op.phase);
        }
    }

    // Drop columns whose accumulated angle vanishes (mod 2*pi).
    std::vector<sul::dynamic_bitset<>> kept_columns;
    std::vector<dvlab::Phase> kept_angles;
    for (size_t j = 0; j < columns.size(); ++j) {
        if (angles[j] != dvlab::Phase()) {
            kept_columns.push_back(std::move(columns[j]));
            kept_angles.push_back(angles[j]);
        }
    }

    PhasePolyProblem problem;
    problem.n_qubits      = n;
    problem.phase_matrix  = ParityMatrix::from_columns(n, kept_columns);
    problem.phase_angles  = std::move(kept_angles);
    problem.output_matrix = state.to_output_matrix();
    return problem;
}

std::vector<PhaseBlock> extract_phase_blocks(qcir::QCir const& circuit) {
    size_t const n = circuit.get_num_qubits();
    std::vector<PhaseBlock> blocks;
    PhaseBlock current(n);

    auto flush = [&]() {
        if (!current.empty()) {
            blocks.push_back(std::move(current));
            current = PhaseBlock(n);
        }
    };

    for (auto const* gate : circuit.get_gates()) {
        auto const& op     = gate->get_operation();
        auto const qubits  = gate->get_qubits();
        if (auto const cx = as_cx(op, qubits)) {
            current.append_cx(cx->first, cx->second);
        } else if (auto const angle = as_rz(op)) {
            current.append_rz(static_cast<size_t>(qubits[0]), *angle);
        } else {
            flush();  // a non-phase gate (e.g. H) is a block boundary
        }
    }
    flush();

    return blocks;
}

}  // namespace qsyn::experimental::phasepoly
