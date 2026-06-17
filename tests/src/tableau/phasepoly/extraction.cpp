/****************************************************************************
  PackageName  [ tests/tableau/phasepoly ]
  Synopsis     [ Stage-1 extraction tests: PhaseBlock -> (P, Theta, O) ]
  Author       [ Design Verification Lab ]
  Copyright    [ Copyright(c) 2025 DVLab, GIEE, NTU, Taiwan ]
****************************************************************************/

#include <catch2/catch_test_macros.hpp>

#include <vector>

#include "qcir/basic_gate_type.hpp"
#include "qcir/qcir.hpp"
#include "tableau/pauli_rotation.hpp"
#include "tableau/phasepoly/extractor.hpp"
#include "tableau/phasepoly/parity_matrix.hpp"
#include "tableau/phasepoly/phase_block.hpp"
#include "tableau/phasepoly/phase_poly_problem.hpp"
#include "util/phase.hpp"

using namespace qsyn;
using namespace qsyn::experimental;
using namespace qsyn::qcir;
using namespace qsyn::experimental::phasepoly;
using dvlab::Phase;

TEST_CASE("Extract a single Rz", "[phasepoly]") {
    PhaseBlock block(1);
    block.append_rz(0, Phase(1, 4));

    auto const problem = phase_block_to_problem(block);
    CHECK(problem.n_qubits == 1);
    CHECK(problem.num_phase_terms() == 1);
    CHECK(problem.phase_matrix == ParityMatrix::from_rows({{1}}));
    CHECK(problem.phase_angles == std::vector<Phase>{Phase(1, 4)});
    CHECK(problem.output_matrix == ParityMatrix::identity(1));
}

TEST_CASE("Extract CX; Rz; CX exposes the XOR parity and restores the basis", "[phasepoly]") {
    // The Rz sees q0 ^ q1; the second CX restores the output basis to identity.
    PhaseBlock block(2);
    block.append_cx(0, 1);
    block.append_rz(1, Phase(1, 2));
    block.append_cx(0, 1);

    auto const problem = phase_block_to_problem(block);
    CHECK(problem.num_phase_terms() == 1);
    CHECK(problem.phase_matrix == ParityMatrix::from_rows({{1}, {1}}));  // column (1,1) = q0^q1
    CHECK(problem.phase_angles == std::vector<Phase>{Phase(1, 2)});
    CHECK(problem.output_matrix == ParityMatrix::identity(2));
}

TEST_CASE("Extract the paper's phase polynomial p = pi/2 (x^y) + pi/4 (y^z)", "[phasepoly]") {
    // Reproduces the left (phase) half of the paper's Eq. (4) joint matrix.
    PhaseBlock block(3);  // qubits x, y, z
    block.append_cx(0, 1);                 // q1 = x ^ y
    block.append_rz(1, Phase(1, 2));       // pi/2 . (x ^ y)
    block.append_cx(0, 1);                 // q1 = y      (restore)
    block.append_cx(2, 1);                 // q1 = y ^ z
    block.append_rz(1, Phase(1, 4));       // pi/4 . (y ^ z)
    block.append_cx(2, 1);                 // q1 = y      (restore)

    auto const problem = phase_block_to_problem(block);
    CHECK(problem.num_phase_terms() == 2);
    // Columns in first-appearance order: (x^y) = (1,1,0), then (y^z) = (0,1,1).
    CHECK(problem.phase_matrix == ParityMatrix::from_rows({{1, 0},
                                                           {1, 1},
                                                           {0, 1}}));
    CHECK(problem.phase_angles == std::vector<Phase>{Phase(1, 2), Phase(1, 4)});
    CHECK(problem.output_matrix == ParityMatrix::identity(3));
}

TEST_CASE("Extract the paper's output basis |x, x^z, x^y^z>", "[phasepoly]") {
    // Reproduces the right (output) half of the paper's Eq. (4) joint matrix.
    PhaseBlock block(3);
    block.append_cx(0, 1);  // q1 = x ^ y
    block.append_cx(1, 2);  // q2 = x ^ y ^ z
    block.append_cx(2, 1);  // q1 = z
    block.append_cx(0, 1);  // q1 = x ^ z

    auto const problem = phase_block_to_problem(block);
    CHECK(problem.num_phase_terms() == 0);
    CHECK(problem.phase_matrix.has_no_columns());
    // Output columns: x = (1,0,0), x^z = (1,0,1), x^y^z = (1,1,1).
    CHECK(problem.output_matrix == ParityMatrix::from_rows({{1, 1, 1},
                                                            {0, 0, 1},
                                                            {0, 1, 1}}));
}

TEST_CASE("Extraction merges identical parity columns", "[phasepoly]") {
    PhaseBlock block(2);
    block.append_cx(0, 1);
    block.append_rz(1, Phase(1, 4));  // pi/4 on q0^q1
    block.append_rz(1, Phase(1, 4));  // another pi/4 on the same parity
    block.append_cx(0, 1);

    auto const problem = phase_block_to_problem(block);
    CHECK(problem.num_phase_terms() == 1);                                   // merged into one
    CHECK(problem.phase_angles == std::vector<Phase>{Phase(1, 2)});          // pi/4 + pi/4
    CHECK(problem.phase_matrix == ParityMatrix::from_rows({{1}, {1}}));
}

TEST_CASE("Extraction cancels inverse-angle parity columns", "[phasepoly]") {
    PhaseBlock block(2);
    block.append_cx(0, 1);
    block.append_rz(1, Phase(1, 4));
    block.append_rz(1, Phase(-1, 4));  // cancels the previous rotation
    block.append_cx(0, 1);

    auto const problem = phase_block_to_problem(block);
    CHECK(problem.num_phase_terms() == 0);
    CHECK(problem.phase_matrix.has_no_columns());
    CHECK(problem.output_matrix == ParityMatrix::identity(2));
}

TEST_CASE("extract_phase_blocks splits at non-phase boundary gates", "[phasepoly]") {
    QCir qcir{2};
    qcir.append(CXGate(), {0, 1});
    qcir.append(PZGate(Phase(1, 4)), {1});
    qcir.append(HGate(), {1});  // boundary: ends the first block
    qcir.append(CXGate(), {0, 1});
    qcir.append(PZGate(Phase(1, 2)), {0});

    auto const blocks = extract_phase_blocks(qcir);
    REQUIRE(blocks.size() == 2);
    CHECK(blocks[0].num_cx() == 1);
    CHECK(blocks[0].num_rz() == 1);
    CHECK(blocks[1].num_cx() == 1);
    CHECK(blocks[1].num_rz() == 1);
}

TEST_CASE("extract_phase_blocks reads CX control/target and Rz correctly", "[phasepoly]") {
    QCir qcir{3};
    qcir.append(CXGate(), {0, 1});            // control 0, target 1
    qcir.append(CXGate(), {2, 1});            // control 2, target 1
    qcir.append(PZGate(Phase(1, 4)), {1});    // Rz on the exposed parity

    auto const blocks = extract_phase_blocks(qcir);
    REQUIRE(blocks.size() == 1);

    auto const& ops = blocks[0].ops();
    REQUIRE(ops.size() == 3);
    CHECK((ops[0].is_cx() && ops[0].control == 0 && ops[0].target == 1));
    CHECK((ops[1].is_cx() && ops[1].control == 2 && ops[1].target == 1));
    CHECK((ops[2].is_rz() && ops[2].qubit() == 1));

    // q1 carries (q0 ^ q1) ^ q2 = q0 ^ q1 ^ q2 when the Rz fires.
    auto const problem = phase_block_to_problem(blocks[0]);
    CHECK(problem.num_phase_terms() == 1);
    CHECK(problem.phase_matrix == ParityMatrix::from_rows({{1}, {1}, {1}}));
    CHECK(problem.phase_angles == std::vector<Phase>{Phase(1, 4)});
}

TEST_CASE("rotations_to_problem extracts diagonal Pauli rotations", "[phasepoly]") {
    std::vector<PauliRotation> rotations{
        PauliRotation({Pauli::z, Pauli::z, Pauli::i}, Phase(1, 2)),
        PauliRotation({Pauli::i, Pauli::z, Pauli::z}, Phase(1, 4)),
    };

    auto const problem = rotations_to_problem(rotations, ParityMatrix::identity(3));
    CHECK(problem.num_phase_terms() == 2);
    CHECK(problem.phase_matrix == ParityMatrix::from_rows({{1, 0},
                                                           {1, 1},
                                                           {0, 1}}));
    CHECK(problem.phase_angles == std::vector<Phase>{Phase(1, 2), Phase(1, 4)});
    CHECK(problem.output_matrix == ParityMatrix::identity(3));
}
