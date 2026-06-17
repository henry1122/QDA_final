/****************************************************************************
  PackageName  [ tableau/phasepoly ]
  Synopsis     [ Build a QCir from synthesized gates; symbolic verification ]
  Author       [ Design Verification Lab ]
  Copyright    [ Copyright(c) 2025 DVLab, GIEE, NTU, Taiwan ]
****************************************************************************/

#pragma once

#include <cstddef>
#include <vector>

#include "./phase_block.hpp"
#include "./phase_poly_problem.hpp"
#include "./search.hpp"

namespace qsyn::qcir {
class QCir;
}

namespace qsyn::experimental::phasepoly {

/// @brief Build a `QCir` (CX + p-gate) from an ordered list of phase-poly ops.
qcir::QCir build_qcir(size_t n_qubits, std::vector<PhaseOp> const& gates);

/**
 * @brief Symbolically verify a synthesis result against its problem.
 *
 * Re-extracts `(P', Theta', O')` from the synthesized gates and checks that the
 * phase-parity/angle multiset and the output basis match the original problem
 * -- exactly the "compare symbolic phase polynomial" criterion (plan §"Important
 * implementation details"). This is exact for phase-polynomial circuits and
 * works at any size (no unitary contraction needed).
 */
bool verify_synthesis(PhasePolyProblem const& problem, SynthesisResult const& result);

}  // namespace qsyn::experimental::phasepoly
