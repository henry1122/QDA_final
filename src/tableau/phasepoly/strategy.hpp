/****************************************************************************
  PackageName  [ tableau/phasepoly ]
  Synopsis     [ PauliRotationsSynthesisStrategy adapter for the PhasePoly A* ]
  Author       [ Design Verification Lab ]
  Copyright    [ Copyright(c) 2025 DVLab, GIEE, NTU, Taiwan ]
****************************************************************************/

#pragma once

#include <optional>
#include <vector>

#include "convert/tableau_to_qcir.hpp"
#include "tableau/phasepoly/config.hpp"

namespace qsyn::experimental::phasepoly {

/**
 * @brief Wraps the PhasePoly co-optimizer as a PauliRotationsSynthesisStrategy.
 *
 * Accepts a list of diagonal PauliRotations (the phase-polynomial segment
 * between two Clifford layers in the Tableau IR), converts them to a
 * PhasePolyProblem with identity output matrix, runs the A* search, and
 * returns the synthesized QCir.
 *
 * Usage: pass to `to_qcir(Tableau, clifford_strategy, PhasePolySynthesisStrategy{})`.
 */
struct PhasePolySynthesisStrategy : public experimental::PauliRotationsSynthesisStrategy {
    explicit PhasePolySynthesisStrategy(PhasePolyConfig config = {})
        : _config(std::move(config)) {}

    std::optional<qcir::QCir> synthesize(
        std::vector<experimental::PauliRotation> const& rotations) const override;

private:
    PhasePolyConfig _config;
};

}  // namespace qsyn::experimental::phasepoly
