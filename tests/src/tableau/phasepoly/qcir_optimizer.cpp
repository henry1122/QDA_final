/****************************************************************************
  PackageName  [ tests/tableau/phasepoly ]
  Synopsis     [ QCir-level PhasePoly optimization tests ]
  Author       [ Design Verification Lab ]
  Copyright    [ Copyright(c) 2025 DVLab, GIEE, NTU, Taiwan ]
****************************************************************************/

#include <catch2/catch_test_macros.hpp>

#include "qcir/basic_gate_type.hpp"
#include "qcir/qcir.hpp"
#include "qcir/qcir_equiv.hpp"
#include "tableau/phasepoly/qcir_optimizer.hpp"
#include "util/phase.hpp"

using namespace qsyn;
using namespace qsyn::qcir;
using namespace qsyn::experimental::phasepoly;
using dvlab::Phase;

namespace {

size_t count_cx(QCir const& circuit) {
    size_t n = 0;
    for (auto const* gate : circuit.get_gates()) {
        if (gate->get_operation().get_underlying_if<ControlGate>()) ++n;
    }
    return n;
}

}  // namespace

TEST_CASE("optimize_qcir_phasepoly reduces redundant CX in a block", "[phasepoly]") {
    QCir qcir{2};
    qcir.append(CXGate(), {0, 1});
    qcir.append(CXGate(), {0, 1});
    qcir.append(PZGate(Phase(1, 4)), {1});

    size_t const cx_before = count_cx(qcir);

    PhasePolyOptimizeStats stats{};
    auto optimized = optimize_qcir_phasepoly(qcir, {}, true, &stats);
    REQUIRE(optimized);
    CHECK(count_cx(*optimized) <= cx_before);
    CHECK(qsyn::qcir::is_equivalent(qcir, *optimized));
    CHECK(stats.blocks_optimized >= 1);
}

TEST_CASE("optimize_qcir_phasepoly preserves H boundary gates", "[phasepoly]") {
    QCir qcir{2};
    qcir.append(CXGate(), {0, 1});
    qcir.append(PZGate(Phase(1, 4)), {1});
    qcir.append(HGate(), {0});
    qcir.append(CXGate(), {0, 1});
    qcir.append(PZGate(Phase(1, 2)), {1});

    auto optimized = optimize_qcir_phasepoly(qcir);
    REQUIRE(optimized);
    CHECK(qsyn::qcir::is_equivalent(qcir, *optimized));
    CHECK(optimized->get_num_gates() >= 3);
}
