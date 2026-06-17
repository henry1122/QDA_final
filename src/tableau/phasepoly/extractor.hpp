/****************************************************************************
  PackageName  [ tableau/phasepoly ]
  Synopsis     [ Extract phase-polynomial blocks and their (P, Theta, O) ]
  Author       [ Design Verification Lab ]
  Copyright    [ Copyright(c) 2025 DVLab, GIEE, NTU, Taiwan ]
****************************************************************************/

#pragma once

#include <vector>

#include "./phase_block.hpp"
#include "./phase_poly_problem.hpp"

namespace qsyn::experimental {
class PauliRotation;
class StabilizerTableau;
class StabilizerTableauSynthesisStrategy;
}  // namespace qsyn::experimental

namespace qsyn::qcir {
class QCir;
}

namespace qsyn::experimental::phasepoly {

/**
 * @brief Convert a phase-polynomial block into its `(P, Theta, O)` problem.
 *
 * Scans the block forward with the physical CNOT semantics (paper
 * §"Single-block extraction"): a CX updates the carried parities, and each `Rz`
 * contributes a phase column equal to the parity its qubit currently holds.
 * Identical parity columns are merged (angles summed) and zero-angle (mod 2*pi)
 * columns are dropped. The output matrix `O` is the final basis transformation.
 */
PhasePolyProblem phase_block_to_problem(PhaseBlock const& block);

/**
 * @brief Partition a circuit into maximal phase-polynomial blocks.
 *
 * Gates are scanned in topological order; CX and single-qubit Z-rotations
 * (`p`/`rz`-family) accumulate into the current block, while any other gate
 * (e.g. H) closes the current block and acts as a boundary. Blocks span all
 * `n` qubits of the circuit so qubit ids map directly to matrix rows.
 */
std::vector<PhaseBlock> extract_phase_blocks(qcir::QCir const& circuit);

/**
 * @brief Build `(P, Theta, O)` from a diagonal Pauli-rotation block.
 *
 * Phase columns are the Z-supports of each rotation (with angle normalization).
 * `output_matrix` defaults to identity when the block has no CNOT prefix.
 */
PhasePolyProblem rotations_to_problem(
    std::vector<PauliRotation> const& rotations,
    ParityMatrix const& output_matrix);

/**
 * @brief Derive the output basis `O` from a CNOT-only stabilizer prefix.
 *
 * Returns `nullopt` when the decomposed Clifford string contains non-CNOT gates.
 */
std::optional<ParityMatrix> output_matrix_from_stabilizer(
    StabilizerTableau const& clifford,
    StabilizerTableauSynthesisStrategy const& strategy);

/**
 * @brief Co-extract a tableau segment `[StabilizerTableau | PauliRotations]`.
 *
 * Returns `nullopt` when the stabilizer prefix is not CNOT-only or the
 * rotations are not diagonal.
 */
std::optional<PhasePolyProblem> tableau_block_to_problem(
    StabilizerTableau const& clifford,
    std::vector<PauliRotation> const& rotations,
    StabilizerTableauSynthesisStrategy const& strategy);

}  // namespace qsyn::experimental::phasepoly
