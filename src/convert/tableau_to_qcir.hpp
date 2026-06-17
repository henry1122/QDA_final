/****************************************************************************
  PackageName  [ tableau ]
  Synopsis     [ Define pauli rotation class ]
  Author       [ Design Verification Lab ]
  Copyright    [ Copyright(c) 2023 DVLab, GIEE, NTU, Taiwan ]
****************************************************************************/

#pragma once

#include "qcir/qcir.hpp"
#include "tableau/pauli_rotation.hpp"
#include "tableau/phasepoly/config.hpp"
#include "tableau/tableau.hpp"

namespace qsyn {

namespace experimental {

struct PauliRotationsSynthesisStrategy {
public:
    virtual ~PauliRotationsSynthesisStrategy() = default;

    virtual std::optional<qcir::QCir> synthesize(std::vector<PauliRotation> const& rotations) const = 0;
};

struct NaivePauliRotationsSynthesisStrategy : public PauliRotationsSynthesisStrategy {
    std::optional<qcir::QCir> synthesize(std::vector<PauliRotation> const& rotations) const override;
};

struct TParPauliRotationsSynthesisStrategy : public PauliRotationsSynthesisStrategy {
    std::optional<qcir::QCir> synthesize(std::vector<PauliRotation> const& rotations) const override;
};

struct GraySynthPauliRotationsSynthesisStrategy : public PauliRotationsSynthesisStrategy {
    enum class Mode { star,
                      staircase };
    Mode mode;
    GraySynthPauliRotationsSynthesisStrategy(Mode mode = Mode::star) : mode(mode) {}
    std::optional<qcir::QCir> synthesize(std::vector<PauliRotation> const& rotations) const override;
};

/**
 * @brief Synthesize Pauli rotations using minimum spanning arborescence.
 * This method is based on the following paper by Vandaele et al.:
 * https://arxiv.org/abs/2104.00934
 *
 */
struct MstSynthesisStrategy : public PauliRotationsSynthesisStrategy {
    std::optional<qcir::QCir> synthesize(std::vector<PauliRotation> const& rotations) const override;
};

/**
 * @brief Synthesize Pauli rotations with PhasePoly co-optimization (paper §3).
 *
 * When paired with a preceding CNOT-only stabilizer block in `to_qcir(Tableau)`,
 * the joint `[P | O]` search is used; otherwise rotations are synthesized with
 * `O = I`.
 */
struct PhasePolySynthesisStrategy : public PauliRotationsSynthesisStrategy {
    phasepoly::PhasePolyConfig config;
    bool replace_only_if_better = true;

    std::optional<qcir::QCir> synthesize(std::vector<PauliRotation> const& rotations) const override;
};

std::optional<qcir::QCir> synthesize_cooptimized_block(
    StabilizerTableau const& clifford_prefix,
    std::vector<PauliRotation> const& rotations,
    PhasePolySynthesisStrategy const& strategy,
    StabilizerTableauSynthesisStrategy const& st_strategy);

std::optional<qcir::QCir> to_qcir(
    StabilizerTableau const& clifford,
    StabilizerTableauSynthesisStrategy const& strategy);
std::optional<qcir::QCir> to_qcir(
    std::vector<PauliRotation> const& rotations,
    PauliRotationsSynthesisStrategy const& strategy);
std::optional<qcir::QCir> to_qcir(
    Tableau const& tableau,
    StabilizerTableauSynthesisStrategy const& st_strategy,
    PauliRotationsSynthesisStrategy const& pr_strategy);

}  // namespace experimental

}  // namespace qsyn
