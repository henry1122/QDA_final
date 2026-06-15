/****************************************************************************
  PackageName  [ tableau/phasepoly ]
  Synopsis     [ Multi-block warm-start group synthesis (paper §3.3) ]
  Author       [ Design Verification Lab ]
  Copyright    [ Copyright(c) 2025 DVLab, GIEE, NTU, Taiwan ]
****************************************************************************/

#include "./multiblock.hpp"

#include <cstddef>
#include <optional>
#include <utility>
#include <vector>

#include "./extractor.hpp"
#include "./parity_matrix.hpp"
#include "./search.hpp"
#include "qcir/basic_gate_type.hpp"
#include "qcir/operation.hpp"
#include "qcir/qcir.hpp"
#include "qcir/qcir_gate.hpp"
#include "qsyn/qsyn_type.hpp"

namespace qsyn::experimental::phasepoly {

namespace {

/// GF(2) matrix multiply: result[r][c] = XOR_k A[r][k] & B[k][c].
/// A is n×n square; B is n×m. Returns n×m result.
ParityMatrix gf2_matmul(ParityMatrix const& A, ParityMatrix const& B) {
    size_t const n = A.n_rows();
    size_t const m = B.n_cols();
    ParityMatrix result(n, m);
    for (size_t r = 0; r < n; ++r)
        for (size_t c = 0; c < m; ++c) {
            bool val = false;
            for (size_t k = 0; k < n; ++k) val ^= (A.get(r, k) & B.get(k, c));
            result.set(r, c, val);
        }
    return result;
}

bool is_h_gate(qcir::Operation const& op) {
    return op.get_underlying_if<qcir::HGate>().has_value();
}

std::optional<std::pair<size_t, size_t>>
as_cx(qcir::Operation const& op, QubitIdList const& qubits) {
    auto const cg = op.get_underlying_if<qcir::ControlGate>();
    if (!cg || cg->get_num_ctrls() != 1) return std::nullopt;
    auto const tgt = cg->get_target_operation().get_underlying_if<qcir::PXGate>();
    if (!tgt || tgt->get_phase() != dvlab::Phase(1)) return std::nullopt;
    return std::pair{static_cast<size_t>(qubits[0]), static_cast<size_t>(qubits[1])};
}

std::optional<dvlab::Phase> as_rz(qcir::Operation const& op) {
    if (auto const p = op.get_underlying_if<qcir::PZGate>()) return p->get_phase();
    if (auto const r = op.get_underlying_if<qcir::RZGate>()) return r->get_phase();
    return std::nullopt;
}

}  // namespace

// ─── Public API ───────────────────────────────────────────────────────────

ExtractedCircuit extract_phase_blocks_with_boundaries(qcir::QCir const& circuit) {
    size_t const n = circuit.get_num_qubits();
    ExtractedCircuit result;
    result.n_qubits = n;

    PhaseBlock current(n);
    std::vector<size_t> pending_h_qubits;
    bool has_flushed = false;

    auto flush = [&]() {
        if (!current.empty()) {
            if (has_flushed) {
                result.boundaries.push_back(BlockBoundary{pending_h_qubits});
                pending_h_qubits.clear();
            }
            result.blocks.push_back(std::move(current));
            current     = PhaseBlock(n);
            has_flushed = true;
        }
    };

    for (auto const* gate : circuit.get_gates()) {
        auto const& op    = gate->get_operation();
        auto const qubits = gate->get_qubits();
        if (auto const cx = as_cx(op, qubits)) {
            current.append_cx(cx->first, cx->second);
        } else if (auto const angle = as_rz(op)) {
            current.append_rz(static_cast<size_t>(qubits[0]), *angle);
        } else {
            flush();
            if (is_h_gate(op))
                pending_h_qubits.push_back(static_cast<size_t>(qubits[0]));
        }
    }
    flush();

    return result;
}

GroupResult synthesize_block_group(
    std::vector<PhaseBlock> const& blocks,
    std::vector<BlockBoundary> const& boundaries,
    PhasePolyConfig const& cfg) {
    if (blocks.empty()) return {};

    size_t const n = blocks[0].n_qubits();
    GroupResult result;

    // M tracks the accumulated A* row-operation matrix (paper §3.3).
    // Initially M = I.  After each block's synthesis, M is updated with the
    // emitted CNOT row operations.  At each H-boundary, H-gated qubit rows are
    // reset to e_j (the H gate restarts qubit j's parity in the new basis).
    //
    // Block i+1 then sees warm-started problem (M·P_{i+1}, M·O_{i+1}): parities
    // and output matrix expressed in the current accumulated basis.
    ParityMatrix M = ParityMatrix::identity(n);

    for (size_t i = 0; i < blocks.size(); ++i) {
        auto const base = phase_block_to_problem(blocks[i]);

        PhasePolyProblem warm;
        warm.n_qubits     = base.n_qubits;
        warm.phase_angles = base.phase_angles;
        warm.phase_matrix  = gf2_matmul(M, base.phase_matrix);
        warm.output_matrix = gf2_matmul(M, base.output_matrix);

        auto const sr = synthesize_phasepoly(warm, cfg);
        result.num_cx += sr.num_cx;
        if (sr.used_fallback) result.used_fallback = true;

        // Update M with all emitted CNOT row ops.
        for (auto const& gate : sr.gates) {
            if (gate.is_cx()) M.apply_cnot(gate.control, gate.target);
        }

        // H-boundary: for each H-gated qubit j, reset row j of M to e_j.
        // This models the H gate re-initialising qubit j's parity basis.
        if (i < boundaries.size()) {
            for (size_t j : boundaries[i].h_qubits) {
                for (size_t c = 0; c < n; ++c) M.set(j, c, c == j);
            }
        }
    }

    return result;
}

size_t synthesize_grouped(ExtractedCircuit const& extracted,
                          size_t group_size,
                          PhasePolyConfig const& cfg) {
    if (group_size == 0) group_size = 1;

    auto const& blocks     = extracted.blocks;
    auto const& boundaries = extracted.boundaries;
    size_t const nb        = blocks.size();
    size_t total_cx        = 0;

    for (size_t start = 0; start < nb; start += group_size) {
        size_t const end = std::min(start + group_size, nb);

        std::vector<PhaseBlock> group_blocks(blocks.begin() + start,
                                             blocks.begin() + end);
        std::vector<BlockBoundary> group_boundaries;
        if (end - start > 1) {
            size_t const bnd_start = start;
            size_t const bnd_end   = std::min(bnd_start + (end - start - 1),
                                              boundaries.size());
            group_boundaries.assign(boundaries.begin() + bnd_start,
                                    boundaries.begin() + bnd_end);
        }

        auto const gr = synthesize_block_group(group_blocks, group_boundaries, cfg);
        total_cx += gr.num_cx;
    }

    return total_cx;
}

}  // namespace qsyn::experimental::phasepoly
