/****************************************************************************
  PackageName  [ tableau/phasepoly ]
  Synopsis     [ QCir-level PhasePoly optimization (per-block resynthesis) ]
  Author       [ Design Verification Lab ]
  Copyright    [ Copyright(c) 2025 DVLab, GIEE, NTU, Taiwan ]
****************************************************************************/

#include "./qcir_optimizer.hpp"

#include "./multiblock.hpp"

namespace qsyn::experimental::phasepoly {

std::optional<qcir::QCir> optimize_qcir_phasepoly(
    qcir::QCir const& circuit,
    PhasePolyConfig const& config,
    bool replace_only_if_better,
    PhasePolyOptimizeStats* stats) {
    return optimize_qcir_phasepoly_multiblock(circuit, config, replace_only_if_better, stats);
}

}  // namespace qsyn::experimental::phasepoly
