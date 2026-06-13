/****************************************************************************
  PackageName  [ tableau/phasepoly ]
  Synopsis     [ Build a QCir from synthesized gates; symbolic verification ]
  Author       [ Design Verification Lab ]
  Copyright    [ Copyright(c) 2025 DVLab, GIEE, NTU, Taiwan ]
****************************************************************************/

#include "./synthesizer.hpp"

#include <algorithm>
#include <string>
#include <utility>
#include <vector>

#include "./extractor.hpp"
#include "./parity_matrix.hpp"
#include "qcir/basic_gate_type.hpp"
#include "qcir/qcir.hpp"
#include "qsyn/qsyn_type.hpp"

namespace qsyn::experimental::phasepoly {

namespace {

/// @brief The sorted (parity-bits, angle) multiset of a (P, Theta) pair.
std::vector<std::pair<std::string, std::string>>
phase_terms(ParityMatrix const& phase, std::vector<dvlab::Phase> const& angles) {
    std::vector<std::pair<std::string, std::string>> terms;
    terms.reserve(phase.n_cols());
    for (size_t c = 0; c < phase.n_cols(); ++c) {
        std::string bits(phase.n_rows(), '0');
        for (size_t r = 0; r < phase.n_rows(); ++r) {
            if (phase.get(r, c)) bits[r] = '1';
        }
        terms.emplace_back(std::move(bits), angles[c].get_print_string());
    }
    std::ranges::sort(terms);
    return terms;
}

}  // namespace

qcir::QCir build_qcir(size_t n_qubits, std::vector<PhaseOp> const& gates) {
    qcir::QCir circuit{n_qubits};
    for (auto const& op : gates) {
        if (op.is_cx()) {
            circuit.append(qcir::CXGate(), QubitIdList{op.control, op.target});
        } else {
            circuit.append(qcir::PZGate(op.phase), QubitIdList{op.qubit()});
        }
    }
    return circuit;
}

bool verify_synthesis(PhasePolyProblem const& problem, SynthesisResult const& result) {
    PhaseBlock block(problem.n_qubits);
    for (auto const& op : result.gates) {
        if (op.is_cx())
            block.append_cx(op.control, op.target);
        else
            block.append_rz(op.qubit(), op.phase);
    }
    auto const rebuilt = phase_block_to_problem(block);

    if (rebuilt.output_matrix != problem.output_matrix) return false;
    return phase_terms(rebuilt.phase_matrix, rebuilt.phase_angles) ==
           phase_terms(problem.phase_matrix, problem.phase_angles);
}

}  // namespace qsyn::experimental::phasepoly
