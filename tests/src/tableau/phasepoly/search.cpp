/****************************************************************************
  PackageName  [ tests/tableau/phasepoly ]
  Synopsis     [ Stage-3 tests: memory-bounded A* phase-polynomial synthesis ]
  Author       [ Design Verification Lab ]
  Copyright    [ Copyright(c) 2025 DVLab, GIEE, NTU, Taiwan ]
****************************************************************************/

#include <catch2/catch_test_macros.hpp>

#include <random>

#include "qcir/qcir.hpp"
#include "qcir/qcir_equiv.hpp"
#include "tableau/phasepoly/config.hpp"
#include "tableau/phasepoly/extractor.hpp"
#include "tableau/phasepoly/phase_block.hpp"
#include "tableau/phasepoly/search.hpp"
#include "tableau/phasepoly/synthesizer.hpp"
#include "util/phase.hpp"

using namespace qsyn::experimental::phasepoly;
using dvlab::Phase;

namespace {

/// @brief Synthesize a block and assert it is correct: symbolic equivalence
///        always, full unitary equivalence for circuits small enough to verify.
SynthesisResult run(PhaseBlock const& block, PhasePolyConfig cfg = {}) {
    auto const problem = phase_block_to_problem(block);
    auto const result  = synthesize_phasepoly(problem, cfg);

    INFO("num_cx=" << result.num_cx << " num_rz=" << result.num_rz
                   << " fallback=" << result.used_fallback);
    CHECK(verify_synthesis(problem, result));
    CHECK(result.num_rz == problem.num_phase_terms());

    if (block.n_qubits() <= 7) {
        auto const original = build_qcir(block.n_qubits(), block.ops());
        auto const synth    = build_qcir(result.n_qubits, result.gates);
        CHECK(qsyn::qcir::is_equivalent(original, synth));
    }
    return result;
}

PhaseBlock random_block(size_t n, std::mt19937& rng) {
    PhaseBlock block(n);
    std::uniform_int_distribution<size_t> qubit(0, n - 1);
    size_t const n_gates = 6 + (rng() % 14);
    for (size_t g = 0; g < n_gates; ++g) {
        if (rng() % 3 == 0) {
            block.append_rz(qubit(rng), Phase(1, 1 << (1 + (rng() % 3))));  // pi/2, pi/4, pi/8
        } else {
            size_t c = qubit(rng);
            size_t t = qubit(rng);
            while (t == c) t = qubit(rng);
            block.append_cx(c, t);
        }
    }
    return block;
}

}  // namespace

TEST_CASE("Synthesis of a single Rz emits no CNOTs", "[phasepoly]") {
    PhaseBlock block(1);
    block.append_rz(0, Phase(1, 4));
    auto const result = run(block);
    CHECK(result.num_cx == 0);
    CHECK(result.num_rz == 1);
}

TEST_CASE("Synthesis of CX; Rz; CX is optimal (2 CNOTs)", "[phasepoly]") {
    PhaseBlock block(2);
    block.append_cx(0, 1);
    block.append_rz(1, Phase(1, 2));
    block.append_cx(0, 1);
    auto const result = run(block);
    CHECK(result.num_cx == 2);
    CHECK(result.num_rz == 1);
}

TEST_CASE("Synthesis removes redundant CNOTs", "[phasepoly]") {
    PhaseBlock block(2);
    block.append_cx(0, 1);
    block.append_cx(0, 1);  // cancels the first
    block.append_rz(0, Phase(1, 4));
    auto const result = run(block);
    CHECK(result.num_cx == 0);  // the parity is just q0; no CNOTs needed
    CHECK(result.num_rz == 1);
}

TEST_CASE("Synthesis merges identical-parity rotations", "[phasepoly]") {
    PhaseBlock block(2);
    block.append_cx(0, 1);
    block.append_rz(1, Phase(1, 4));
    block.append_rz(1, Phase(1, 4));
    block.append_cx(0, 1);
    auto const result = run(block);
    CHECK(result.num_rz == 1);   // merged into pi/2
    CHECK(result.num_cx == 2);
}

TEST_CASE("Synthesis cancels inverse rotations to the empty circuit", "[phasepoly]") {
    PhaseBlock block(2);
    block.append_cx(0, 1);
    block.append_rz(1, Phase(1, 4));
    block.append_rz(1, Phase(-1, 4));
    block.append_cx(0, 1);
    auto const result = run(block);
    CHECK(result.num_cx == 0);
    CHECK(result.num_rz == 0);
}

TEST_CASE("Synthesis of the paper's two-parity phase polynomial", "[phasepoly]") {
    PhaseBlock block(3);  // p = pi/2 (x^y) + pi/4 (y^z), output basis = identity
    block.append_cx(0, 1);
    block.append_rz(1, Phase(1, 2));
    block.append_cx(0, 1);
    block.append_cx(2, 1);
    block.append_rz(1, Phase(1, 4));
    block.append_cx(2, 1);
    auto const result = run(block);
    CHECK(result.num_rz == 2);
    CHECK(result.num_cx <= block.num_cx());  // no worse than the source
}

TEST_CASE("Synthesis preserves a non-identity output basis", "[phasepoly]") {
    PhaseBlock block(3);
    block.append_cx(0, 1);              // q1 = x ^ y
    block.append_rz(1, Phase(1, 4));    // phase on x ^ y
    block.append_cx(1, 2);              // q2 = x ^ y ^ z   (non-identity output)
    block.append_rz(2, Phase(1, 2));    // phase on x ^ y ^ z
    auto const result = run(block);
    CHECK(result.num_rz == 2);
}

TEST_CASE("Synthesis is correct on random phase-polynomial blocks", "[phasepoly]") {
    std::mt19937 rng(20250613);
    for (size_t trial = 0; trial < 20; ++trial) {
        size_t const n = 3 + (rng() % 4);  // 3..6 qubits
        run(random_block(n, rng));
    }
}

TEST_CASE("Greedy fallback still yields a correct circuit", "[phasepoly]") {
    PhaseBlock block(3);
    block.append_cx(0, 1);
    block.append_rz(1, Phase(1, 2));
    block.append_cx(2, 1);
    block.append_rz(1, Phase(1, 4));

    PhasePolyConfig cfg;
    cfg.max_expansions = 0;  // force the fallback path
    auto const result = run(block, cfg);
    CHECK(result.used_fallback);
}

TEST_CASE("A tiny memory bound still yields a correct circuit", "[phasepoly]") {
    std::mt19937 rng(20250614);
    PhasePolyConfig cfg;
    cfg.max_queue_size = 2;  // aggressively small open set
    cfg.max_solutions  = 1;
    for (size_t trial = 0; trial < 5; ++trial) {
        run(random_block(4, rng), cfg);
    }
}
