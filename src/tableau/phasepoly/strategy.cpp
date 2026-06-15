/****************************************************************************
  PackageName  [ tableau/phasepoly ]
  Synopsis     [ PauliRotationsSynthesisStrategy adapter for the PhasePoly A* ]
  Author       [ Design Verification Lab ]
  Copyright    [ Copyright(c) 2025 DVLab, GIEE, NTU, Taiwan ]
****************************************************************************/

#include "./strategy.hpp"

#include <spdlog/spdlog.h>

#include <map>
#include <string>

#include "./parity_matrix.hpp"
#include "./phase_poly_problem.hpp"
#include "./search.hpp"
#include "./synthesizer.hpp"

namespace qsyn::experimental::phasepoly {

/**
 * @brief Synthesize a phase-polynomial segment from the Tableau IR.
 *
 * Converts the list of diagonal PauliRotations to a PhasePolyProblem with
 * identity output matrix (the flanking StabilizerTableau layers handle the
 * Clifford basis changes), deduplicates identical parity columns, drops
 * zero-angle terms, then runs the A* search.
 */
std::optional<qcir::QCir> PhasePolySynthesisStrategy::synthesize(
    std::vector<experimental::PauliRotation> const& rotations) const {
    if (rotations.empty()) return qcir::QCir{0};

    size_t const n = rotations.front().n_qubits();
    if (n == 0) return qcir::QCir{0};

    if (!std::ranges::all_of(rotations, &experimental::PauliRotation::is_diagonal)) {
        spdlog::error("PhasePoly strategy only supports diagonal rotations");
        return std::nullopt;
    }

    // Deduplicate: identical parities accumulate their angles.
    // Key: string of '0'/'1' of length n encoding the Z-support bitmask.
    std::map<std::string, size_t> key_to_idx;
    std::vector<std::string>      parity_keys;
    std::vector<dvlab::Phase>     angles;

    for (auto const& rot : rotations) {
        std::string key(n, '0');
        for (size_t q = 0; q < n; ++q)
            if (rot.is_z(q)) key[q] = '1';

        auto [it, inserted] = key_to_idx.emplace(key, parity_keys.size());
        if (inserted) {
            parity_keys.push_back(key);
            angles.push_back(rot.phase());
        } else {
            angles[it->second] += rot.phase();
        }
    }

    // Drop zero-angle (mod 2π) terms.
    std::vector<std::string>  live_keys;
    std::vector<dvlab::Phase> live_angles;
    for (size_t j = 0; j < angles.size(); ++j) {
        if (angles[j] != dvlab::Phase(0)) {
            live_keys.push_back(parity_keys[j]);
            live_angles.push_back(angles[j]);
        }
    }

    size_t const m = live_keys.size();

    // Build PhasePolyProblem: identity output matrix because the Clifford
    // layers in the Tableau handle all basis transformations.
    PhasePolyProblem problem;
    problem.n_qubits      = n;
    problem.output_matrix = ParityMatrix::identity(n);
    problem.phase_matrix  = ParityMatrix(n, m);
    problem.phase_angles  = live_angles;

    for (size_t j = 0; j < m; ++j)
        for (size_t q = 0; q < n; ++q)
            problem.phase_matrix.set(q, j, live_keys[j][q] == '1');

    auto const result = synthesize_phasepoly(problem, _config);
    return build_qcir(result.n_qubits, result.gates);
}

}  // namespace qsyn::experimental::phasepoly
