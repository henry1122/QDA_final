/****************************************************************************
  PackageName  [ tests/tableau/phasepoly ]
  Synopsis     [ Stage-5 tests: multi-block warm-start group synthesis ]
  Author       [ Design Verification Lab ]
  Copyright    [ Copyright(c) 2025 DVLab, GIEE, NTU, Taiwan ]
****************************************************************************/

#include <catch2/catch_test_macros.hpp>

#include <vector>

#include "qcir/basic_gate_type.hpp"
#include "qcir/qcir.hpp"
#include "tableau/phasepoly/config.hpp"
#include "tableau/phasepoly/extractor.hpp"
#include "tableau/phasepoly/multiblock.hpp"
#include "tableau/phasepoly/phase_block.hpp"
#include "tableau/phasepoly/search.hpp"
#include "util/phase.hpp"

using namespace qsyn::experimental::phasepoly;
using namespace qsyn::qcir;
using dvlab::Phase;

// ---------------------------------------------------------------------------
// Test 1: k=1 warm start equals single block synthesis
// ---------------------------------------------------------------------------
TEST_CASE("k=1 warm start equals single block synthesis", "[multiblock]") {
    // Build a 3-qubit block with one parity term {q0,q1}
    PhaseBlock block(3);
    block.append_cx(0, 1);
    block.append_rz(1, Phase(1, 4));  // pi/4 on q0^q1
    block.append_cx(0, 1);

    std::vector<PhaseBlock> blocks{block};
    std::vector<BlockBoundary> no_boundaries;

    PhasePolyConfig cfg;

    auto const gr = synthesize_block_group(blocks, no_boundaries, cfg);

    // Compare to independent synthesis
    auto const problem = phase_block_to_problem(block);
    auto const sr      = synthesize_phasepoly(problem, cfg);

    CHECK(gr.num_cx == sr.num_cx);
}

// ---------------------------------------------------------------------------
// Test 2: 2-block warm start reuses shared parity
// ---------------------------------------------------------------------------
TEST_CASE("2-block warm start reuses shared parity", "[multiblock]") {
    // Block A: parity {q1,q2}, angle pi/4
    // Output: non-identity (the CX is not undone)
    PhaseBlock blockA(3);
    blockA.append_cx(1, 2);
    blockA.append_rz(2, Phase(1, 4));  // pi/4 on q1^q2
    // Note: we don't restore, so output is non-identity

    // Boundary: H on qubit 0 (doesn't affect q1,q2 parity)
    BlockBoundary boundary;
    boundary.h_qubits = {0};

    // Block B: same parity {q1,q2}, angle pi/4
    // With warm start, M already has the CX(1,2) row op, so
    // the parity column should immediately be at weight 1.
    PhaseBlock blockB(3);
    blockB.append_cx(1, 2);
    blockB.append_rz(2, Phase(1, 4));  // pi/4 on q1^q2

    std::vector<PhaseBlock> blocks{blockA, blockB};
    std::vector<BlockBoundary> boundaries{boundary};

    PhasePolyConfig cfg;

    auto const gr_k2 = synthesize_block_group(blocks, boundaries, cfg);

    // Independent (k=1): each block synthesized separately
    size_t cx_k1 = 0;
    for (auto const& b : blocks) {
        auto const p  = phase_block_to_problem(b);
        auto const sr = synthesize_phasepoly(p, cfg);
        cx_k1 += sr.num_cx;
    }

    // Warm start should be no worse than independent
    INFO("k=2 CX=" << gr_k2.num_cx << " k=1 CX=" << cx_k1);
    CHECK(gr_k2.num_cx <= cx_k1);
}

// ---------------------------------------------------------------------------
// Test 3: H qubit reset — reuse only for unaffected qubits
// ---------------------------------------------------------------------------
TEST_CASE("H qubit reset: reuse only for unaffected qubits", "[multiblock]") {
    // Block A builds parity {q1,q2}
    PhaseBlock blockA(3);
    blockA.append_cx(1, 2);
    blockA.append_rz(2, Phase(1, 4));

    // H is on qubit 1 (resets q1's row in M — parity {q1,q2} can no longer reuse q1)
    BlockBoundary boundary;
    boundary.h_qubits = {1};

    // Block B: needs parity {q1,q2} — q1 row was reset, so full rebuild needed
    PhaseBlock blockB(3);
    blockB.append_cx(1, 2);
    blockB.append_rz(2, Phase(1, 4));

    std::vector<PhaseBlock> blocks{blockA, blockB};
    std::vector<BlockBoundary> boundaries{boundary};

    PhasePolyConfig cfg;

    auto const gr_k2 = synthesize_block_group(blocks, boundaries, cfg);

    // k=1 for comparison
    size_t cx_k1 = 0;
    for (auto const& b : blocks) {
        auto const p  = phase_block_to_problem(b);
        auto const sr = synthesize_phasepoly(p, cfg);
        cx_k1 += sr.num_cx;
    }

    // Even with qubit reset, k=2 should be no worse than k=1
    INFO("k=2 CX=" << gr_k2.num_cx << " k=1 CX=" << cx_k1);
    CHECK(gr_k2.num_cx <= cx_k1);
}

// ---------------------------------------------------------------------------
// Test 4: synthesize_grouped partitions correctly
// ---------------------------------------------------------------------------
TEST_CASE("synthesize_grouped partitions correctly", "[multiblock]") {
    // Build a 3-block circuit with 2 boundaries
    // blocks[0]: CX(0,1) + Rz(1, pi/4)
    // boundary[0]: H on q2
    // blocks[1]: CX(1,2) + Rz(2, pi/4)
    // boundary[1]: H on q0
    // blocks[2]: CX(0,2) + Rz(2, pi/4)

    ExtractedCircuit extracted;
    extracted.n_qubits = 3;

    PhaseBlock b0(3);
    b0.append_cx(0, 1);
    b0.append_rz(1, Phase(1, 4));
    extracted.blocks.push_back(b0);

    BlockBoundary bnd0;
    bnd0.h_qubits = {2};
    extracted.boundaries.push_back(bnd0);

    PhaseBlock b1(3);
    b1.append_cx(1, 2);
    b1.append_rz(2, Phase(1, 4));
    extracted.blocks.push_back(b1);

    BlockBoundary bnd1;
    bnd1.h_qubits = {0};
    extracted.boundaries.push_back(bnd1);

    PhaseBlock b2(3);
    b2.append_cx(0, 2);
    b2.append_rz(2, Phase(1, 4));
    extracted.blocks.push_back(b2);

    PhasePolyConfig cfg;

    // group_size=2: groups [0,1] and [2]
    size_t const cx_k2 = synthesize_grouped(extracted, 2, cfg);

    // Verify by manually computing
    // Group [0,1]
    std::vector<PhaseBlock> g01{b0, b1};
    std::vector<BlockBoundary> bnd01{bnd0};
    auto const gr01 = synthesize_block_group(g01, bnd01, cfg);

    // Group [2]
    std::vector<PhaseBlock> g2{b2};
    std::vector<BlockBoundary> bnd2{};
    auto const gr2 = synthesize_block_group(g2, bnd2, cfg);

    CHECK(cx_k2 == gr01.num_cx + gr2.num_cx);
}

// ---------------------------------------------------------------------------
// Test 5: extract_phase_blocks_with_boundaries captures H qubits
// ---------------------------------------------------------------------------
TEST_CASE("extract_phase_blocks_with_boundaries captures H qubits", "[multiblock]") {
    // Mirror the pattern from the existing extraction test (which passes):
    // CX(0,1) → PZ(1) → H(1) → CX(0,1) → PZ(0)
    // This circuit structure ensures H comes after Rz via qubit-1 wire dependency,
    // and the second CX syncs both wires so the second block has correct ordering.
    QCir qcir{2};
    qcir.append(CXGate(), {0, 1});
    qcir.append(PZGate(Phase(1, 4)), {1});  // Rz on qubit 1
    qcir.append(HGate(), {1});              // boundary: H on qubit 1 (same wire → proper ordering)
    qcir.append(CXGate(), {0, 1});
    qcir.append(PZGate(Phase(1, 2)), {0});  // second block Rz on qubit 0

    auto const extracted = extract_phase_blocks_with_boundaries(qcir);

    REQUIRE(extracted.blocks.size() == 2);
    REQUIRE(extracted.boundaries.size() == 1);

    CHECK(extracted.n_qubits == 2);
    // Block 0: CX + Rz
    CHECK(extracted.blocks[0].num_cx() == 1);
    CHECK(extracted.blocks[0].num_rz() == 1);
    // Block 1: CX + Rz
    CHECK(extracted.blocks[1].num_cx() == 1);
    CHECK(extracted.blocks[1].num_rz() == 1);

    // The boundary should record H on qubit 1
    REQUIRE(extracted.boundaries[0].h_qubits.size() == 1);
    CHECK(extracted.boundaries[0].h_qubits[0] == 1);
}

// ---------------------------------------------------------------------------
// Test 6: extract_phase_blocks_with_boundaries matches extract_phase_blocks
// ---------------------------------------------------------------------------
TEST_CASE("extract_phase_blocks_with_boundaries: blocks match extract_phase_blocks",
          "[multiblock]") {
    // Verify that extract_phase_blocks_with_boundaries produces the same blocks
    // as extract_phase_blocks (same gate contents, same count).
    // Use the same circuit as the existing extraction test that is known to work.
    QCir qcir{2};
    qcir.append(CXGate(), {0, 1});
    qcir.append(PZGate(Phase(1, 4)), {1});
    qcir.append(HGate(), {1});              // H on q1 (same wire as PZ — correct ordering)
    qcir.append(CXGate(), {0, 1});
    qcir.append(PZGate(Phase(1, 2)), {0});

    auto const blocks_orig   = extract_phase_blocks(qcir);
    auto const extracted     = extract_phase_blocks_with_boundaries(qcir);

    // Same number of blocks
    REQUIRE(extracted.blocks.size() == blocks_orig.size());
    CHECK(extracted.blocks.size() == 2);
    CHECK(extracted.boundaries.size() == 1);

    // Each block has same op count
    for (size_t i = 0; i < blocks_orig.size(); ++i) {
        CHECK(extracted.blocks[i].num_cx() == blocks_orig[i].num_cx());
        CHECK(extracted.blocks[i].num_rz() == blocks_orig[i].num_rz());
    }

    // Boundary records the H on qubit 1
    CHECK(extracted.boundaries[0].h_qubits.size() == 1);
    CHECK(extracted.boundaries[0].h_qubits[0] == 1);
}

// ---------------------------------------------------------------------------
// Test 7: empty circuit produces no blocks
// ---------------------------------------------------------------------------
TEST_CASE("extract_phase_blocks_with_boundaries: empty circuit", "[multiblock]") {
    QCir qcir{3};
    auto const extracted = extract_phase_blocks_with_boundaries(qcir);
    CHECK(extracted.blocks.empty());
    CHECK(extracted.boundaries.empty());
    CHECK(extracted.n_qubits == 3);
}

// ---------------------------------------------------------------------------
// Test 8: synthesize_grouped with group_size=1 matches block-level synthesis
// ---------------------------------------------------------------------------
TEST_CASE("synthesize_grouped k=1 matches per-block synthesis", "[multiblock]") {
    ExtractedCircuit extracted;
    extracted.n_qubits = 3;

    PhaseBlock b0(3);
    b0.append_cx(0, 1);
    b0.append_rz(1, Phase(1, 4));
    extracted.blocks.push_back(b0);

    BlockBoundary bnd;
    bnd.h_qubits = {2};
    extracted.boundaries.push_back(bnd);

    PhaseBlock b1(3);
    b1.append_cx(1, 2);
    b1.append_rz(2, Phase(1, 2));
    extracted.blocks.push_back(b1);

    PhasePolyConfig cfg;

    size_t const cx_grouped_k1 = synthesize_grouped(extracted, 1, cfg);

    // Compute per-block independently
    size_t cx_independent = 0;
    for (auto const& b : extracted.blocks) {
        auto const p  = phase_block_to_problem(b);
        auto const sr = synthesize_phasepoly(p, cfg);
        cx_independent += sr.num_cx;
    }

    CHECK(cx_grouped_k1 == cx_independent);
}
