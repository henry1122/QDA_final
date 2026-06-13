/****************************************************************************
  PackageName  [ tests/tableau/phasepoly ]
  Synopsis     [ Stage-0 convention lock-in tests for the PhasePoly co-optimizer ]
  Author       [ Design Verification Lab ]
  Copyright    [ Copyright(c) 2025 DVLab, GIEE, NTU, Taiwan ]
****************************************************************************/

#include <catch2/catch_test_macros.hpp>

#include <random>
#include <string_view>
#include <utility>
#include <vector>

#include "tableau/pauli_rotation.hpp"
#include "tableau/phasepoly/parity_matrix.hpp"
#include "tableau/phasepoly/symbolic_state.hpp"
#include "util/boolean_matrix.hpp"
#include "util/phase.hpp"

using namespace qsyn::experimental;
using namespace qsyn::experimental::phasepoly;
using dvlab::Phase;

namespace {

// Reduce a square ParityMatrix to identity using only paper-convention CNOT row
// operations `apply_cnot(control, target)` (row[control] ^= row[target]) and
// record them in order. By the duality `O = T^T`, replaying each recorded op as
// a physical CNOT(control, target) in the same order reconstructs the original
// linear map -- this is the property exercised by the round-trip test below.
std::vector<std::pair<size_t, size_t>> reduce_to_identity(ParityMatrix m) {
    size_t const n = m.n_rows();
    std::vector<std::pair<size_t, size_t>> ops;
    for (size_t q = 0; q < n; ++q) {
        if (!m.get(q, q)) {
            // rows < q are already unit vectors e_r with a 0 in column q, so a
            // pivot for column q must live in some row r > q (m is invertible).
            for (size_t r = q + 1; r < n; ++r) {
                if (m.get(r, q)) {
                    m.apply_cnot(q, r);  // row[q] ^= row[r] : set the pivot
                    ops.emplace_back(q, r);
                    break;
                }
            }
        }
        for (size_t r = 0; r < n; ++r) {
            if (r != q && m.get(r, q)) {
                m.apply_cnot(r, q);  // row[r] ^= row[q] : clear the entry
                ops.emplace_back(r, q);
            }
        }
    }
    REQUIRE(m.is_square_identity());
    return ops;
}

}  // namespace

TEST_CASE("ParityMatrix::apply_cnot XORs the control row", "[phasepoly]") {
    auto m = ParityMatrix::from_rows({{1, 0, 0},
                                      {0, 1, 0},
                                      {0, 0, 1}});
    m.apply_cnot(/*control=*/1, /*target=*/0);  // row1 <- row1 XOR row0
    CHECK(m.get(1, 0));                          // row1 now has a 1 in column 0
    CHECK(m.get(1, 1));                          // row1 keeps its own 1
    CHECK(m.row(0) == ParityMatrix::from_rows({{1, 0, 0}}).row(0));  // row0 unchanged
}

TEST_CASE("PhasePoly reproduces the paper Eq. (4) row operation", "[phasepoly]") {
    // Joint [P | O] matrix from the paper (Fig. 6 example). Columns:
    //   0: x^y, 1: y^z  (phase parities)   2: x, 3: x^z, 4: x^y^z  (outputs)
    auto m = ParityMatrix::from_rows({{1, 0, 1, 1, 1},
                                      {1, 1, 0, 0, 1},
                                      {0, 1, 0, 1, 1}});
    // Paper "CNOT(2,1)" = control row 2, target row 1 (1-based) -> (1, 0) here.
    m.apply_cnot(/*control=*/1, /*target=*/0);

    auto const expected = ParityMatrix::from_rows({{1, 0, 1, 1, 1},
                                                   {0, 1, 1, 1, 0},
                                                   {0, 1, 0, 1, 1}});
    CHECK(m == expected);
}

TEST_CASE("PhasePoly row op matches PauliProduct::cx Z-action", "[phasepoly]") {
    // A diagonal (Z-only) Pauli rotation's support is exactly a parity column.
    auto rotation = PauliRotation(std::string_view{"ZZIZ"}, Phase(1, 4));

    // Build the same parity as a single-column ParityMatrix (row = qubit).
    auto column = ParityMatrix(4, 1);
    for (size_t q = 0; q < 4; ++q) {
        if (rotation.is_z(q)) column.set(q, 0, true);
    }

    size_t const control = 0, target = 2;
    rotation.cx(control, target);            // z[control] ^= z[target]
    column.apply_cnot(control, target);      // row[control] ^= row[target]

    for (size_t q = 0; q < 4; ++q) {
        CHECK(rotation.is_z(q) == column.get(q, 0));
    }
}

TEST_CASE("PhasePoly row op is the mirror of BooleanMatrix::row_operation", "[phasepoly]") {
    // dvlab::BooleanMatrix::row_operation(ctrl, targ) does row[targ] ^= row[ctrl].
    // Our apply_cnot(control, target) does row[control] ^= row[target]. Hence
    // apply_cnot(c, t) == row_operation(ctrl = t, targ = c).
    std::vector<std::vector<unsigned char>> data = {{1, 0, 1},
                                                    {1, 1, 0},
                                                    {0, 1, 1}};
    dvlab::BooleanMatrix bm;
    for (auto const& row : data) bm.push_row(dvlab::BooleanMatrix::Row(row));

    auto pm = ParityMatrix::from_rows({{1, 0, 1},
                                       {1, 1, 0},
                                       {0, 1, 1}});

    size_t const control = 2, target = 0;
    pm.apply_cnot(control, target);
    bm.row_operation(/*ctrl=*/target, /*targ=*/control);

    for (size_t r = 0; r < 3; ++r) {
        for (size_t c = 0; c < 3; ++c) {
            CHECK(static_cast<bool>(bm[r][c]) == pm.get(r, c));
        }
    }
}

TEST_CASE("Forward simulator reproduces the Fig. 1 phase-polynomial parities", "[phasepoly]") {
    // Fig. 1(a): inputs x, y, z on qubits 0, 1, 2.
    SymbolicState state(3);

    // CNOT(control = x, target = y): the Rz(pi/2) parity becomes x ^ y on qubit 1.
    apply_cnot_to_state(state, /*control=*/0, /*target=*/1);
    CHECK(state.parity_of(1) == ParityMatrix::from_columns(3, {[] {
                                    sul::dynamic_bitset<> b(3);
                                    b.set(0);
                                    b.set(1);
                                    return b;
                                }()})
                                    .column(0));

    // CNOT(control = y, target = z): qubit 2 now carries x ^ y ^ z.
    apply_cnot_to_state(state, /*control=*/1, /*target=*/2);
    auto const q2 = state.parity_of(2);
    CHECK(q2.test(0));  // x
    CHECK(q2.test(1));  // y
    CHECK(q2.test(2));  // z
    CHECK(q2.count() == 3);
}

TEST_CASE("PhasePoly column utilities behave as specified", "[phasepoly]") {
    auto m = ParityMatrix::from_rows({{1, 1, 0},
                                      {0, 1, 0},
                                      {0, 1, 0}});
    CHECK(m.column_weight(0) == 1);
    CHECK(m.column_weight(1) == 3);
    CHECK(m.column_weight(2) == 0);

    CHECK(m.single_one_row(0) == std::optional<size_t>{0});  // ready: emit Rz on qubit 0
    CHECK(m.single_one_row(1) == std::nullopt);              // weight 3
    CHECK(m.single_one_row(2) == std::nullopt);              // weight 0

    m.remove_column(0);
    CHECK(m.n_cols() == 2);
    CHECK(m.column_weight(0) == 3);  // old column 1 shifted into position 0
}

TEST_CASE("PhasePoly convention round-trip reconstructs the linear map", "[phasepoly]") {
    // For random CNOT networks, build the output matrix O by forward simulation,
    // reduce O to identity with apply_cnot row operations, then emit a physical
    // CNOT(control, target) per recorded op (same order) and forward-simulate
    // again. The reconstructed output matrix must equal O.
    std::mt19937 rng(20250611);

    for (size_t trial = 0; trial < 200; ++trial) {
        size_t const n = 3 + (rng() % 4);  // 3..6 qubits
        std::uniform_int_distribution<size_t> qubit(0, n - 1);

        SymbolicState forward(n);
        size_t const n_gates = 5 + (rng() % 20);
        for (size_t g = 0; g < n_gates; ++g) {
            size_t c = qubit(rng);
            size_t t = qubit(rng);
            while (t == c) t = qubit(rng);
            apply_cnot_to_state(forward, c, t);
        }

        auto const output = forward.to_output_matrix();
        auto const ops    = reduce_to_identity(output);

        SymbolicState replay(n);
        for (auto const& [control, target] : ops) {
            apply_cnot_to_state(replay, control, target);
        }
        CHECK(replay.to_output_matrix() == output);
    }
}
