/****************************************************************************
  PackageName  [ tableau/phasepoly ]
  Synopsis     [ The (P, Theta, O) problem fed to the PhasePoly synthesizer ]
  Author       [ Design Verification Lab ]
  Copyright    [ Copyright(c) 2025 DVLab, GIEE, NTU, Taiwan ]
****************************************************************************/

#pragma once

#include <cstddef>
#include <string>
#include <vector>

#include "./parity_matrix.hpp"
#include "util/phase.hpp"

namespace qsyn::experimental::phasepoly {

/**
 * @brief The co-optimization problem extracted from a phase-polynomial block.
 *
 * It bundles the three components the joint search operates on (paper §3.1):
 *   - `phase_matrix`  (P): n x m, one column per surviving `Rz` parity term;
 *   - `phase_angles`  (Theta): the angle of each phase column (P and Theta are
 *                     kept in lockstep -- column `j` carries angle `Theta[j]`);
 *   - `output_matrix` (O): n x n, the basis transformation the synthesized
 *                     block must reproduce (goal: reduce it to identity).
 *
 * The phase columns have already been angle-normalized: identical parities are
 * merged (angles summed) and zero-angle (mod 2*pi) terms are dropped.
 */
struct PhasePolyProblem {
    size_t n_qubits = 0;
    ParityMatrix phase_matrix;               // P : n_qubits x num_phase_terms()
    std::vector<dvlab::Phase> phase_angles;  // Theta : length num_phase_terms()
    ParityMatrix output_matrix;              // O : n_qubits x n_qubits

    size_t num_phase_terms() const { return phase_angles.size(); }

    /// @brief True when there is nothing left to synthesize (no phase terms and
    ///        the output is already the identity basis).
    bool is_trivial() const {
        return phase_matrix.has_no_columns() && output_matrix.is_square_identity();
    }

    std::string to_string() const {
        std::string out = "PhasePolyProblem(" + std::to_string(n_qubits) + " qubits, " +
                          std::to_string(num_phase_terms()) + " phase terms)\n";
        out += "P:\n" + phase_matrix.to_string();
        out += "Theta: [";
        for (size_t j = 0; j < phase_angles.size(); ++j) {
            if (j) out += ", ";
            out += phase_angles[j].get_print_string();
        }
        out += "]\nO:\n" + output_matrix.to_string();
        return out;
    }
};

}  // namespace qsyn::experimental::phasepoly
