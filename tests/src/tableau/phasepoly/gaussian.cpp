/****************************************************************************
  PackageName  [ tests/tableau/phasepoly ]
  Synopsis     [ Stage-2 tests: CNOT-only linear-reversible synthesis (O -> I) ]
  Author       [ Design Verification Lab ]
  Copyright    [ Copyright(c) 2025 DVLab, GIEE, NTU, Taiwan ]
****************************************************************************/

#include <catch2/catch_test_macros.hpp>

#include <random>
#include <vector>

#include "tableau/phasepoly/extractor.hpp"
#include "tableau/phasepoly/gaussian.hpp"
#include "tableau/phasepoly/parity_matrix.hpp"
#include "tableau/phasepoly/phase_block.hpp"
#include "tableau/phasepoly/symbolic_state.hpp"

using namespace qsyn::experimental::phasepoly;
using dvlab::Phase;

namespace {

constexpr LinearSynthesisMode kModes[] = {LinearSynthesisMode::gauss_jordan,
                                          LinearSynthesisMode::patel_markov_hayes};

/// @brief Apply the ops as matrix row operations to a copy of `matrix`.
ParityMatrix reduce_with(ParityMatrix matrix, std::vector<CnotOp> const& ops) {
    for (auto const& op : ops) matrix.apply_cnot(op.control, op.target);
    return matrix;
}

/// @brief Emit the ops as physical CNOTs from identity and read back the map.
ParityMatrix replay(size_t n, std::vector<CnotOp> const& ops) {
    SymbolicState state(n);
    for (auto const& op : ops) state.apply_cnot(op.control, op.target);
    return state.to_output_matrix();
}

/// @brief A random invertible matrix, built as the map of a random CNOT network.
ParityMatrix random_invertible(size_t n, std::mt19937& rng) {
    SymbolicState state(n);
    std::uniform_int_distribution<size_t> qubit(0, n - 1);
    size_t const n_gates = 4 + (rng() % (4 * n));
    for (size_t g = 0; g < n_gates; ++g) {
        size_t c = qubit(rng);
        size_t t = qubit(rng);
        while (t == c) t = qubit(rng);
        state.apply_cnot(c, t);
    }
    return state.to_output_matrix();
}

}  // namespace

TEST_CASE("Linear synthesis of the identity emits no CNOTs", "[phasepoly]") {
    for (auto mode : kModes) {
        auto const ops = synthesize_linear_reversible(ParityMatrix::identity(5), mode);
        CHECK(ops.empty());
        CHECK(linear_reversible_cnot_cost(ParityMatrix::identity(5), mode) == 0);
    }
}

TEST_CASE("Linear synthesis reproduces the paper's output basis", "[phasepoly]") {
    // O = |x, x^z, x^y^z>, the right block of Eq. (4).
    auto const target = ParityMatrix::from_rows({{1, 1, 1},
                                                 {0, 0, 1},
                                                 {0, 1, 1}});
    for (auto mode : kModes) {
        auto const ops = synthesize_linear_reversible(target, mode);
        CHECK(reduce_with(target, ops).is_square_identity());  // ops reduce O -> I
        CHECK(replay(3, ops) == target);                       // emitting them rebuilds O
    }
}

TEST_CASE("Linear synthesis reduces and rebuilds random invertible matrices", "[phasepoly]") {
    std::mt19937 rng(20250612);
    for (size_t trial = 0; trial < 150; ++trial) {
        size_t const n = 2 + (rng() % 6);  // 2..7 qubits
        auto const target = random_invertible(n, rng);

        for (auto mode : kModes) {
            auto const ops = synthesize_linear_reversible(target, mode);
            CHECK(reduce_with(target, ops).is_square_identity());
            CHECK(replay(n, ops) == target);
            CHECK(ops.size() <= n * n);  // reasonable: well under the naive bound
            CHECK(linear_reversible_cnot_cost(target, mode) == ops.size());
        }
    }
}

TEST_CASE("Linear synthesis finishes a Stage-1 extracted output matrix", "[phasepoly]") {
    // The 4-CNOT network from the Eq. (4) output-basis extraction test.
    PhaseBlock block(3);
    block.append_cx(0, 1);
    block.append_cx(1, 2);
    block.append_cx(2, 1);
    block.append_cx(0, 1);

    auto const problem = phase_block_to_problem(block);
    REQUIRE(problem.num_phase_terms() == 0);

    for (auto mode : kModes) {
        auto const ops = synthesize_linear_reversible(problem.output_matrix, mode);
        CHECK(reduce_with(problem.output_matrix, ops).is_square_identity());
        CHECK(replay(problem.n_qubits, ops) == problem.output_matrix);
    }
}
