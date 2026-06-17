/****************************************************************************
  PackageName  [ tests/tableau/phasepoly ]
  Synopsis     [ Multi-block SSA merge and pyzx gate-alias tests ]
  Author       [ Design Verification Lab ]
  Copyright    [ Copyright(c) 2025 DVLab, GIEE, NTU, Taiwan ]
****************************************************************************/

#include <catch2/catch_test_macros.hpp>

#include "qcir/basic_gate_type.hpp"
#include "qcir/operation.hpp"
#include "qcir/qcir.hpp"
#include "qcir/qcir_equiv.hpp"
#include "qcir/qcir_io.hpp"
#include "tableau/phasepoly/config.hpp"
#include "tableau/phasepoly/multiblock.hpp"
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

TEST_CASE("str_to_operation accepts pyzx T* and S* aliases", "[phasepoly]") {
    auto const tdg = str_to_operation("T*");
    REQUIRE(tdg);
    auto const pz = tdg->get_underlying_if<PZGate>();
    REQUIRE(pz);
    CHECK(pz->get_phase() == Phase(-1, 4));

    auto const sdg = str_to_operation("S*");
    REQUIRE(sdg);
    auto const pz2 = sdg->get_underlying_if<PZGate>();
    REQUIRE(pz2);
    CHECK(pz2->get_phase() == Phase(-1, 2));
}

TEST_CASE("read pyzx qc file with T* gates", "[phasepoly]") {
    auto const qcir = from_qc("benchmark/qc/optimized/tof_3_pyzx.qc");
    REQUIRE(qcir);
    CHECK(qcir->get_num_qubits() == 5);
    CHECK(qcir->get_num_gates() > 0);
}

TEST_CASE("merge_segments_ssa joins blocks across H boundaries", "[phasepoly]") {
    QCir qcir{2};
    qcir.append(CXGate(), {0, 1});
    qcir.append(PZGate(Phase(1, 4)), {1});
    qcir.append(HGate(), {0});
    qcir.append(CXGate(), {0, 1});
    qcir.append(PZGate(Phase(1, 2)), {1});

    auto const segments = segment_phase_poly_regions(qcir);
    REQUIRE(segments.size() >= 2);

    auto const merged = merge_segments_ssa(segments, 0, 2);
    REQUIRE(merged);
    CHECK(merged->n_qubits() >= 2);
    CHECK(merged->num_ops() > segments.front().block.num_ops());
}

TEST_CASE("optimize_qcir_phasepoly multiblock preserves equivalence on H-split circuit", "[phasepoly]") {
    QCir qcir{2};
    qcir.append(CXGate(), {0, 1});
    qcir.append(PZGate(Phase(1, 4)), {1});
    qcir.append(HGate(), {0});
    qcir.append(CXGate(), {0, 1});
    qcir.append(CXGate(), {0, 1});
    qcir.append(PZGate(Phase(1, 2)), {1});

    PhasePolyConfig config;
    config.group_sizes = {1, 2};

    auto optimized = optimize_qcir_phasepoly(qcir, config);
    REQUIRE(optimized);
    CHECK(qsyn::qcir::is_equivalent(qcir, *optimized));
    CHECK(count_cx(*optimized) <= count_cx(qcir));
}
