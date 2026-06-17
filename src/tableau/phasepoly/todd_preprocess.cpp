/****************************************************************************
  PackageName  [ tableau/phasepoly ]
  Synopsis     [ Block-level Todd phase-polynomial preprocessing ]
  Author       [ Design Verification Lab ]
  Copyright    [ Copyright(c) 2025 DVLab, GIEE, NTU, Taiwan ]
****************************************************************************/

#include "./todd_preprocess.hpp"

#include "./extractor.hpp"
#include "tableau/stabilizer_tableau.hpp"
#include "tableau/tableau_optimization.hpp"

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

}  // namespace

std::optional<PhasePolyProblem> todd_optimize_problem(PhasePolyProblem const& problem) {
    auto const rotations = problem_to_rotations(problem);
    if (rotations.size() < 2) return std::nullopt;
    if (!is_phase_polynomial(rotations)) return std::nullopt;

    StabilizerTableau clifford(problem.n_qubits);
    auto const [opt_clifford, opt_rotations] =
        ToddPhasePolynomialOptimizationStrategy{}.optimize(clifford, rotations);

    if (opt_rotations.size() >= rotations.size()) return std::nullopt;

    if (auto const opt = tableau_block_to_problem(opt_clifford, opt_rotations, AGSynthesisStrategy{})) {
        return opt;
    }
    return rotations_to_problem(opt_rotations, problem.output_matrix);
}

void todd_preprocess_segments(std::vector<PhasePolySegment>& segments) {
    for (auto& segment : segments) {
        if (segment.block.empty()) continue;
        auto const problem = phase_block_to_problem(segment.block);
        if (auto const opt = todd_optimize_problem(problem)) {
            // Keep the block IR unchanged; Todd is applied at synthesis time.
            (void)opt;
        }
    }
}

size_t count_rz_after_todd(std::vector<PhasePolySegment> const& segments) {
    size_t n = 0;
    for (auto const& segment : segments) {
        if (segment.block.empty()) continue;
        auto const problem = phase_block_to_problem(segment.block);
        if (auto const opt = todd_optimize_problem(problem)) {
            n += opt->num_phase_terms();
        } else {
            n += segment.block.num_rz();
        }
    }
    return n;
}

}  // namespace qsyn::experimental::phasepoly
