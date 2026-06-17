/****************************************************************************
  PackageName  [ tableau/phasepoly ]
  Synopsis     [ Todd preprocessing and structure-preserving tableau export ]
  Author       [ Design Verification Lab ]
  Copyright    [ Copyright(c) 2025 DVLab, GIEE, NTU, Taiwan ]
****************************************************************************/

#pragma once

#include <optional>

#include "qcir/qcir.hpp"
#include "tableau/tableau.hpp"

namespace qsyn::experimental::phasepoly {

/// @brief Convert a circuit to tableau, run Todd phase-polynomial optimization.
std::optional<Tableau> todd_optimize_tableau(qcir::QCir const& circuit);

/**
 * @brief Export a Todd-optimized tableau back to a QCir suitable for multiblock
 *        segmentation (CX + PZ/RZ + boundary Cliffords).
 *
 * Returns nullopt when a rotation block cannot be expressed in the current parity
 * basis (e.g. after Hadamard-bearing Clifford prefixes).
 */
std::optional<qcir::QCir> tableau_to_phase_qcir(Tableau const& tableau);

/// @brief Todd-optimize then export; nullopt if either step fails.
std::optional<qcir::QCir> todd_preprocess_qcir(qcir::QCir const& circuit);

}  // namespace qsyn::experimental::phasepoly
