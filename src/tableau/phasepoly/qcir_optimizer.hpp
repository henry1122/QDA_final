/****************************************************************************
  PackageName  [ tableau/phasepoly ]
  Synopsis     [ QCir-level PhasePoly optimization (per-block resynthesis) ]
  Author       [ Design Verification Lab ]
  Copyright    [ Copyright(c) 2025 DVLab, GIEE, NTU, Taiwan ]
****************************************************************************/

#pragma once

#include <cstddef>
#include <optional>

#include "./config.hpp"

namespace qsyn::qcir {
class QCir;
}

namespace qsyn::experimental::phasepoly {

struct PhasePolyOptimizeStats {
    size_t num_blocks       = 0;
    size_t blocks_optimized = 0;
    size_t cx_before        = 0;
    size_t cx_after         = 0;
};

/**
 * @brief Resynthesize every phase-polynomial region (CX + Rz) in a circuit.
 *
 * Non-phase gates (e.g. H, Toffoli) are preserved in place. Each block is
 * optimized with `synthesize_phasepoly`; the original block is kept when
 * synthesis fails verification or does not reduce the CNOT count.
 */
std::optional<qcir::QCir> optimize_qcir_phasepoly(
    qcir::QCir const& circuit,
    PhasePolyConfig const& config           = {},
    bool replace_only_if_better             = true,
    PhasePolyOptimizeStats* stats           = nullptr);

}  // namespace qsyn::experimental::phasepoly
