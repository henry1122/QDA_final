/****************************************************************************
  PackageName  [ tableau/phasepoly ]
  Synopsis     [ Phase-polynomial block IR (a sequence of CX / Rz operations) ]
  Author       [ Design Verification Lab ]
  Copyright    [ Copyright(c) 2025 DVLab, GIEE, NTU, Taiwan ]
****************************************************************************/

#pragma once

#include <cstddef>
#include <vector>

#include "util/phase.hpp"

namespace qsyn::experimental::phasepoly {

/**
 * @brief A single operation inside a phase-polynomial block: either a CNOT or
 *        an `Rz`/phase rotation. These are the only gates that may appear in a
 *        phase polynomial (paper §2).
 */
struct PhaseOp {
    enum class Kind { cx,
                      rz };

    Kind kind;
    size_t control;      ///< cx: control qubit; rz: the rotated qubit
    size_t target;       ///< cx: target qubit;  rz: equals `control`
    dvlab::Phase phase;  ///< rz: rotation angle; cx: unused

    static PhaseOp make_cx(size_t control, size_t target) {
        return PhaseOp{Kind::cx, control, target, dvlab::Phase()};
    }
    static PhaseOp make_rz(size_t qubit, dvlab::Phase const& phase) {
        return PhaseOp{Kind::rz, qubit, qubit, phase};
    }

    bool is_cx() const { return kind == Kind::cx; }
    bool is_rz() const { return kind == Kind::rz; }
    size_t qubit() const { return control; }  ///< readable accessor for `rz`
};

/**
 * @brief A maximal run of phase-polynomial gates (CX and Rz) over `n_qubits`
 *        qubits, as extracted from a circuit between non-phase boundary gates.
 *
 * This is the input IR to the PhasePoly synthesizer: `phase_block_to_problem`
 * turns it into the `[P | O]` matrices the search operates on.
 */
class PhaseBlock {
public:
    explicit PhaseBlock(size_t n_qubits) : _n_qubits{n_qubits} {}

    size_t n_qubits() const { return _n_qubits; }

    void append_cx(size_t control, size_t target) {
        _ops.push_back(PhaseOp::make_cx(control, target));
        ++_num_cx;
    }
    void append_rz(size_t qubit, dvlab::Phase const& phase) {
        _ops.push_back(PhaseOp::make_rz(qubit, phase));
        ++_num_rz;
    }

    std::vector<PhaseOp> const& ops() const { return _ops; }
    size_t num_ops() const { return _ops.size(); }
    size_t num_cx() const { return _num_cx; }
    size_t num_rz() const { return _num_rz; }
    bool empty() const { return _ops.empty(); }

    /// @brief Build a block from an explicit operation list.
    static PhaseBlock from_ops(size_t n_qubits, std::vector<PhaseOp> ops) {
        PhaseBlock block{n_qubits};
        for (auto const& op : ops) {
            if (op.is_cx())
                block.append_cx(op.control, op.target);
            else
                block.append_rz(op.qubit(), op.phase);
        }
        return block;
    }

    /// @brief Append another block's operations (same `n_qubits` required).
    void append_block(PhaseBlock const& other) {
        for (auto const& op : other.ops()) {
            if (op.is_cx())
                append_cx(op.control, op.target);
            else
                append_rz(op.qubit(), op.phase);
        }
    }

private:
    size_t _n_qubits;
    std::vector<PhaseOp> _ops;
    size_t _num_cx = 0;
    size_t _num_rz = 0;
};

}  // namespace qsyn::experimental::phasepoly
