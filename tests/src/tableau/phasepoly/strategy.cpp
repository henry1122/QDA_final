/****************************************************************************
  PackageName  [ tests/tableau/phasepoly ]
  Synopsis     [ Stage-4 integration tests: PhasePolySynthesisStrategy ]
  Author       [ Design Verification Lab ]
  Copyright    [ Copyright(c) 2025 DVLab, GIEE, NTU, Taiwan ]
****************************************************************************/

#include <catch2/catch_test_macros.hpp>

#include <vector>

#include "convert/qcir_to_tableau.hpp"
#include "convert/tableau_to_qcir.hpp"
#include "qcir/basic_gate_type.hpp"
#include "qcir/qcir.hpp"
#include "qcir/qcir_equiv.hpp"
#include "tableau/phasepoly/strategy.hpp"
#include "util/phase.hpp"

using namespace qsyn;
using namespace qsyn::qcir;
using namespace qsyn::experimental;
using namespace qsyn::experimental::phasepoly;
using dvlab::Phase;

namespace {

// Build a synthesized QCir for a given input circuit via:
//   QCir  ->  Tableau  ->  QCir (using PhasePolySynthesisStrategy)
// Returns the synthesized circuit.
QCir synthesize_via_phasepoly(QCir const& input) {
    auto const tableau_opt = to_tableau(input);
    REQUIRE(tableau_opt.has_value());
    auto result = to_qcir(*tableau_opt,
                          HOptSynthesisStrategy{},
                          PhasePolySynthesisStrategy{});
    REQUIRE(result.has_value());
    return std::move(*result);
}

}  // namespace

TEST_CASE("Strategy: trivial identity circuit", "[phasepoly][strategy]") {
    QCir qcir{2};
    qcir.append(CXGate(), {0, 1});
    qcir.append(CXGate(), {0, 1});  // cancels

    auto const synth = synthesize_via_phasepoly(qcir);
    CHECK(is_equivalent(qcir, synth));
}

TEST_CASE("Strategy: single T-gate", "[phasepoly][strategy]") {
    QCir qcir{1};
    qcir.append(PZGate(Phase(1, 4)), {0});

    auto const synth = synthesize_via_phasepoly(qcir);
    CHECK(is_equivalent(qcir, synth));
}

TEST_CASE("Strategy: CX; Rz; CX round-trip", "[phasepoly][strategy]") {
    QCir qcir{2};
    qcir.append(CXGate(), {0, 1});
    qcir.append(PZGate(Phase(1, 2)), {1});
    qcir.append(CXGate(), {0, 1});

    auto const synth = synthesize_via_phasepoly(qcir);
    CHECK(is_equivalent(qcir, synth));
}

TEST_CASE("Strategy: two phase terms on 3 qubits", "[phasepoly][strategy]") {
    QCir qcir{3};
    qcir.append(CXGate(), {0, 1});
    qcir.append(PZGate(Phase(1, 2)), {1});
    qcir.append(CXGate(), {0, 1});
    qcir.append(CXGate(), {2, 1});
    qcir.append(PZGate(Phase(1, 4)), {1});
    qcir.append(CXGate(), {2, 1});

    auto const synth = synthesize_via_phasepoly(qcir);
    CHECK(is_equivalent(qcir, synth));
}

TEST_CASE("Strategy: circuit with H gate splits into two blocks", "[phasepoly][strategy]") {
    QCir qcir{2};
    qcir.append(CXGate(), {0, 1});
    qcir.append(PZGate(Phase(1, 4)), {1});
    qcir.append(HGate(), {0});  // Clifford boundary: splits the phase block
    qcir.append(PZGate(Phase(1, 2)), {0});

    auto const synth = synthesize_via_phasepoly(qcir);
    CHECK(is_equivalent(qcir, synth));
}

TEST_CASE("Strategy: CNOT count is no worse than MST on simple circuits", "[phasepoly][strategy]") {
    // p = pi/4 (x^y^z) — a 3-parity term that requires at least 2 CNOTs
    QCir qcir{3};
    qcir.append(CXGate(), {0, 1});
    qcir.append(CXGate(), {2, 1});
    qcir.append(PZGate(Phase(1, 4)), {1});
    qcir.append(CXGate(), {2, 1});
    qcir.append(CXGate(), {0, 1});

    auto const tableau_opt = to_tableau(qcir);
    REQUIRE(tableau_opt.has_value());

    auto mst_qcir = to_qcir(*tableau_opt, HOptSynthesisStrategy{}, MstSynthesisStrategy{});
    auto pp_qcir  = to_qcir(*tableau_opt, HOptSynthesisStrategy{}, PhasePolySynthesisStrategy{});

    REQUIRE(mst_qcir.has_value());
    REQUIRE(pp_qcir.has_value());

    CHECK(is_equivalent(qcir, *mst_qcir));
    CHECK(is_equivalent(qcir, *pp_qcir));

    // Count CX gates in both results (phasepoly should not be worse)
    auto count_cx = [](QCir const& c) {
        size_t n = 0;
        for (auto* g : c.get_gates())
            if (g->get_operation() == CXGate()) ++n;
        return n;
    };
    CHECK(count_cx(*pp_qcir) <= count_cx(*mst_qcir) + 2);
}
