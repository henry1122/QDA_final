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

}  // namespace qsyn::experimental::phasepoly
