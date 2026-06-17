/****************************************************************************
  PackageName  [ tableau/phasepoly ]
  Synopsis     [ Per-block synthesis for Table 1 (MST / Gray / PhasePoly / ...) ]
  Author       [ Design Verification Lab ]
  Copyright    [ Copyright(c) 2025 DVLab, GIEE, NTU, Taiwan ]
****************************************************************************/

#include "./table1_synthesizer.hpp"

#include "./extractor.hpp"
#include "./gaussian.hpp"
#include "./multiblock.hpp"
#include "./qcir_optimizer.hpp"
#include "./search.hpp"
#include "./synthesizer.hpp"
#include "./todd_preprocess.hpp"
#include "convert/tableau_to_qcir.hpp"
#include "qcir/basic_gate_type.hpp"
#include "tableau/pauli_rotation.hpp"

namespace qsyn::experimental::phasepoly {

namespace {

std::vector<PauliRotation> problem_to_rotations(PhasePolyProblem const& problem) {
    std::vector<PauliRotation> rotations;
    for (size_t c = 0; c < problem.num_phase_terms(); ++c) {
        std::vector<Pauli> paulis(problem.n_qubits, Pauli::i);
        for (size_t r = 0; r < problem.n_qubits; ++r) {
            if (problem.phase_matrix.get(r, c)) paulis[r] = Pauli::z;
        }
        rotations.emplace_back(paulis.begin(), paulis.end(), problem.phase_angles[c]);
    }
    return rotations;
}

std::optional<qcir::QCir> finish_output_matrix(ParityMatrix const& output) {
    if (output.is_square_identity()) {
        return qcir::QCir{output.n_rows()};
    }
    auto const ops = synthesize_linear_reversible(output, LinearSynthesisMode::patel_markov_hayes);
    qcir::QCir circuit{output.n_rows()};
    for (auto const& op : ops) {
        if (op.control == op.target) continue;
        circuit.append(qcir::CXGate(), QubitIdList{op.control, op.target});
    }
    return circuit;
}

std::optional<qcir::QCir> synthesize_block_rotations(
    PhasePolyProblem const& problem,
    experimental::PauliRotationsSynthesisStrategy const& strategy) {
    auto const rotations = problem_to_rotations(problem);
    qcir::QCir result{problem.n_qubits};

    if (!rotations.empty()) {
        auto rot_circ = strategy.synthesize(rotations);
        if (!rot_circ) return std::nullopt;
        result.compose(*rot_circ);
    }

    if (auto o_circ = finish_output_matrix(problem.output_matrix)) {
        result.compose(*o_circ);
    }
    return result;
}

std::optional<qcir::QCir> synthesize_block_phasepoly(
    PhasePolyProblem const& problem,
    PhasePolyConfig const& config) {
    if (problem.is_trivial()) {
        return finish_output_matrix(problem.output_matrix);
    }
    auto const result = synthesize_phasepoly(problem, config);
    if (!verify_synthesis(problem, result)) return std::nullopt;
    return build_qcir(result.n_qubits, result.gates);
}

std::optional<qcir::QCir> synthesize_single_block(
    PhaseBlock const& block,
    Table1Strategy strategy,
    PhasePolyConfig const& config) {
    auto problem = phase_block_to_problem(block);
    if (config.apply_block_todd) {
        if (auto const opt = todd_optimize_problem(problem)) {
            problem = *opt;
        }
    }
    switch (strategy) {
        case Table1Strategy::phasepoly:
            return synthesize_block_phasepoly(problem, config);
        case Table1Strategy::mst:
            return synthesize_block_rotations(problem, experimental::MstSynthesisStrategy{});
        case Table1Strategy::gray:
            return synthesize_block_rotations(
                problem, experimental::GraySynthPauliRotationsSynthesisStrategy{});
        case Table1Strategy::gstair:
            return synthesize_block_rotations(
                problem,
                experimental::GraySynthPauliRotationsSynthesisStrategy{
                    experimental::GraySynthPauliRotationsSynthesisStrategy::Mode::staircase});
        case Table1Strategy::naive:
            return synthesize_block_rotations(problem, experimental::NaivePauliRotationsSynthesisStrategy{});
    }
    return std::nullopt;
}

std::optional<qcir::QCir> assemble_from_segments(
    std::vector<PhasePolySegment> segments,
    Table1Strategy strategy,
    PhasePolyConfig const& config,
    size_t* phase_gate_out = nullptr) {
    size_t const n = segments.empty() ? 0 : segments.front().block.n_qubits();
    qcir::QCir result{n};

    for (auto const& seg : segments) {
        if (!seg.block.empty()) {
            auto block_circ = synthesize_single_block(seg.block, strategy, config);
            if (!block_circ) return std::nullopt;
            if (phase_gate_out) *phase_gate_out += table1_phase_cx_rz_count(*block_circ);
            result.compose(*block_circ);
        }
        for (auto const& b : seg.trailing) {
            result.append(b.op, b.qubits);
        }
    }
    return result;
}

}  // namespace

size_t table1_phase_cx_rz_count(qcir::QCir const& circuit) {
    size_t n = 0;
    for (auto const* gate : circuit.get_gates()) {
        auto const& op = gate->get_operation();
        if (op.get_underlying_if<qcir::ControlGate>()) {
            ++n;
        } else if (op.get_underlying_if<qcir::PZGate>() || op.get_underlying_if<qcir::RZGate>()) {
            ++n;
        }
    }
    return n;
}

size_t table1_gate_count(qcir::QCir const& circuit) {
    auto const stat = qcir::get_gate_statistics(circuit);
    auto const clifford = stat.contains("clifford") ? stat.at("clifford") : 0;
    auto const t_family = stat.contains("t-family") ? stat.at("t-family") : 0;
    return clifford + t_family;
}

std::optional<qcir::QCir> synthesize_table1(
    qcir::QCir const& circuit,
    Table1Strategy strategy,
    PhasePolyConfig const& config,
    size_t* phase_gate_out) {
    if (strategy == Table1Strategy::phasepoly) {
        return optimize_qcir_phasepoly_multiblock(circuit, config, true, nullptr, phase_gate_out);
    }
    return assemble_from_segments(segment_phase_poly_regions(circuit), strategy, config, phase_gate_out);
}

}  // namespace qsyn::experimental::phasepoly
