/****************************************************************************
  PackageName  [ tableau/phasepoly ]
  Synopsis     [ Multi-block phase-poly synthesis: joint A* for two consecutive
                 blocks (Stage 6) — a single CNOT reduces parities from both
                 blocks simultaneously, achieving true cross-block reuse ]
  Author       [ Design Verification Lab ]
  Copyright    [ Copyright(c) 2025 DVLab, GIEE, NTU, Taiwan ]
****************************************************************************/

#include "./multiblock.hpp"

#include <cstddef>
#include <cstdint>
#include <limits>
#include <map>
#include <optional>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "./extractor.hpp"
#include "./gaussian.hpp"
#include "./parity_matrix.hpp"
#include "./search.hpp"
#include "qcir/basic_gate_type.hpp"
#include "qcir/operation.hpp"
#include "qcir/qcir.hpp"
#include "qcir/qcir_gate.hpp"
#include "qsyn/qsyn_type.hpp"

namespace qsyn::experimental::phasepoly {

namespace {

// ─── Gate-classification helpers (used by extract_phase_blocks_with_boundaries) ─

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

// ─── Stage-6 joint A* helpers ─────────────────────────────────────────────

/// Per-column weight sum — a lower bound on the CNOTs needed for a phase matrix.
size_t phase_cost_mat(ParityMatrix const& m) {
    size_t s = 0;
    for (size_t c = 0; c < m.n_cols(); ++c) s += m.column_weight(c);
    return s;
}

/// Pairs (i,j) where apply_cnot(i,j) reduces at least one column of `phase`.
std::vector<std::pair<size_t, size_t>> phase_active_pairs(ParityMatrix const& phase) {
    size_t const n = phase.n_rows();
    std::vector<std::pair<size_t, size_t>> pairs;
    for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j < n; ++j) {
            if (i == j) continue;
            for (size_t c = 0; c < phase.n_cols(); ++c)
                if (phase.get(i, c) && phase.get(j, c)) { pairs.emplace_back(i, j); break; }
        }
    return pairs;
}

/// Returns false if the square matrix `mat` has GF(2) rank < n_rows.
bool gf2_is_invertible(ParityMatrix mat) {
    size_t const n = mat.n_rows();
    size_t pivot_row = 0;
    for (size_t col = 0; col < n && pivot_row < n; ++col) {
        size_t found = n;
        for (size_t r = pivot_row; r < n; ++r)
            if (mat.get(r, col)) { found = r; break; }
        if (found == n) return false;
        if (found != pivot_row)
            for (size_t c = 0; c < n; ++c) {
                bool const a = mat.get(pivot_row, c), b = mat.get(found, c);
                mat.set(pivot_row, c, b);
                mat.set(found, c, a);
            }
        for (size_t r = 0; r < n; ++r)
            if (r != pivot_row && mat.get(r, col)) mat.apply_cnot(r, pivot_row);
        ++pivot_row;
    }
    return pivot_row == n;
}

/// Joint A* node for two consecutive phase-polynomial blocks.
///
/// Before h_done: every CNOT applied to the state acts on all four matrices
///   (phase0, output0, phase1, output1) simultaneously — the "joint benefit".
/// After  h_done: block 0 is complete (phase0 empty, output0 = I, H applied);
///   subsequent CNOTs act on phase1 and output1 only.
struct JointState {
    ParityMatrix phase0, phase1;    // remaining phase columns
    ParityMatrix output0, output1;  // output matrices
    ParityMatrix phase1_orig;       // original P1 (for H-boundary row restore)
    ParityMatrix output1_orig;      // original O1 (for H-boundary row restore)
    std::vector<dvlab::Phase> ang0, ang1;
    std::vector<PhaseOp> gates;
    size_t g_cost = 0;
    bool h_done   = false;
};

/// Eagerly emit every weight-1 phase column as an Rz gate and remove it.
/// `from_phase1 = false` → operate on phase0/ang0; `true` → phase1/ang1.
void emit_ready(JointState& s, bool from_phase1) {
    auto& ph  = from_phase1 ? s.phase1 : s.phase0;
    auto& ang = from_phase1 ? s.ang1   : s.ang0;
    bool changed = true;
    while (changed) {
        changed = false;
        for (size_t c = 0; c < ph.n_cols(); ++c) {
            if (auto const row = ph.single_one_row(c)) {
                s.gates.push_back(PhaseOp::make_rz(*row, ang[c]));
                ph.remove_column(c);
                ang.erase(ang.begin() + static_cast<std::ptrdiff_t>(c));
                changed = true;
                break;
            }
        }
    }
}

/// Deterministically restore output0 → I.
/// The same row operations also transform output1 and phase1 (cross-block
/// side-effect: output0-restore CNOTs can pre-reduce block 1 matrices for free).
void finish_output0_joint(JointState& s, LinearSynthesisMode mode) {
    for (auto const& op : synthesize_linear_reversible(s.output0, mode)) {
        s.output0.apply_cnot(op.control, op.target);
        s.output1.apply_cnot(op.control, op.target);
        s.phase1.apply_cnot(op.control, op.target);
        s.gates.push_back(PhaseOp::make_cx(op.control, op.target));
        ++s.g_cost;
    }
}

/// Deterministically restore output1 → I (called at goal; phase1 is empty).
void finish_output1_joint(JointState& s, LinearSynthesisMode mode) {
    for (auto const& op : synthesize_linear_reversible(s.output1, mode)) {
        s.output1.apply_cnot(op.control, op.target);
        s.gates.push_back(PhaseOp::make_cx(op.control, op.target));
        ++s.g_cost;
    }
}

/// Apply the H-boundary: for each H-gated qubit j, row j of phase1 and output1 is
/// restored from the ORIGINAL P1 and O1 values.  This models M_new[j] ← eⱼ in
/// accumulated-row-operation terms — qubit j starts fresh in the new H basis.
void apply_h_transition(JointState& s, std::vector<size_t> const& h_qubits) {
    for (size_t j : h_qubits) {
        for (size_t c = 0; c < s.phase1.n_cols(); ++c)
            s.phase1.set(j, c, s.phase1_orig.get(j, c));
        for (size_t c = 0; c < s.output1.n_cols(); ++c)
            s.output1.set(j, c, s.output1_orig.get(j, c));
    }
    s.h_done = true;
}

/// Canonical string key for the joint state (A* revisit detection).
std::string joint_key(JointState const& s) {
    std::string key;
    if (!s.h_done) {
        for (size_t c = 0; c < s.phase0.n_cols(); ++c) {
            for (size_t r = 0; r < s.phase0.n_rows(); ++r)
                key += s.phase0.get(r, c) ? '1' : '0';
            key += ':';
            key += s.ang0[c].get_print_string();
            key += ';';
        }
        key += '#';
        for (size_t r = 0; r < s.output0.n_rows(); ++r)
            for (size_t c = 0; c < s.output0.n_cols(); ++c)
                key += s.output0.get(r, c) ? '1' : '0';
        key += '|';
    }
    key += s.h_done ? '1' : '0';
    key += '|';
    for (size_t c = 0; c < s.phase1.n_cols(); ++c) {
        for (size_t r = 0; r < s.phase1.n_rows(); ++r)
            key += s.phase1.get(r, c) ? '1' : '0';
        key += ':';
        key += s.ang1[c].get_print_string();
        key += ';';
    }
    key += '#';
    for (size_t r = 0; r < s.output1.n_rows(); ++r)
        for (size_t c = 0; c < s.output1.n_cols(); ++c)
            key += s.output1.get(r, c) ? '1' : '0';
    return key;
}

/// Run joint A* for two consecutive blocks p0 and p1 with H-boundary h_qubits.
///
/// Feasibility pre-check: M_new = (O0_orig⁻¹ with H-qubit rows → eⱼ) must be
/// invertible.  Equivalent (cheaper) check: Q = identity with H-qubit rows
/// replaced by the corresponding rows of O0_orig must be invertible
/// (proof: M_new · O0_orig = Q, O0_orig invertible ⟹ M_new invertible ↔ Q invertible).
///
/// Falls back to independent per-block A* if infeasible or budget exhausted.
GroupResult synthesize_joint(PhasePolyProblem const& p0,
                             PhasePolyProblem const& p1,
                             std::vector<size_t> const& h_qubits,
                             PhasePolyConfig const& cfg) {
    using Priority = std::tuple<size_t, size_t, size_t, size_t, uint64_t>;

    // ── Feasibility check ──────────────────────────────────────────────────
    size_t const n = p0.n_qubits;
    {
        ParityMatrix Q = ParityMatrix::identity(n);
        for (size_t j : h_qubits)
            for (size_t c = 0; c < n; ++c)
                Q.set(j, c, p0.output_matrix.get(j, c));
        if (!gf2_is_invertible(Q)) {
            // H-boundary singularity: joint A* is infeasible; synthesize independently.
            GroupResult r;
            r.num_cx        = synthesize_phasepoly(p0, cfg).num_cx +
                              synthesize_phasepoly(p1, cfg).num_cx;
            r.used_fallback = true;
            return r;
        }
    }

    // ── Initial joint state ────────────────────────────────────────────────
    JointState root;
    root.phase0       = p0.phase_matrix;
    root.ang0         = p0.phase_angles;
    root.phase1       = p1.phase_matrix;
    root.ang1         = p1.phase_angles;
    root.output0      = p0.output_matrix;
    root.output1      = p1.output_matrix;
    root.phase1_orig  = p1.phase_matrix;
    root.output1_orig = p1.output_matrix;

    emit_ready(root, false);  // eager Rz from phase0

    // When phase0 becomes empty, deterministically finish output0 and cross the
    // H boundary.  Called on root and after every CNOT that empties phase0.
    auto finalize_block0 = [&](JointState& s) {
        if (!s.h_done && s.phase0.has_no_columns()) {
            finish_output0_joint(s, cfg.finish_mode);
            apply_h_transition(s, h_qubits);
            emit_ready(s, true);  // eager Rz from phase1 in the new H basis
        }
    };
    finalize_block0(root);

    // ── A* search ─────────────────────────────────────────────────────────
    uint64_t counter = 0;
    std::map<Priority, JointState> open;
    std::unordered_map<std::string, size_t> best_g;

    auto push = [&](JointState s) {
        auto const key = joint_key(s);
        if (auto const it = best_g.find(key);
            it != best_g.end() && it->second <= s.g_cost) return;
        best_g[key] = s.g_cost;
        // h1 = total remaining phase work (both blocks); guides the search toward
        // CNOTs that simultaneously reduce parities in both blocks.
        // h2 = output-matrix cost for the currently-active block.
        size_t const h1 = phase_cost_mat(s.phase0) + phase_cost_mat(s.phase1);
        size_t const h2 = s.h_done ? linear_reversible_cnot_cost(s.output1, cfg.h2_mode)
                                   : linear_reversible_cnot_cost(s.output0, cfg.h2_mode);
        Priority const prio{s.g_cost + h1 + h2, h1, h2,
                            std::numeric_limits<size_t>::max() - s.g_cost, counter++};
        open.emplace(prio, std::move(s));
    };

    push(root);

    std::vector<JointState> solutions;
    size_t expansions = 0;

    while (!open.empty() && expansions < cfg.max_expansions) {
        auto node    = open.extract(open.begin());
        JointState s = std::move(node.mapped());
        ++expansions;

        // Goal: h_done + phase1 empty.
        if (s.h_done && s.phase1.has_no_columns()) {
            if (!s.output1.is_square_identity()) finish_output1_joint(s, cfg.finish_mode);
            solutions.push_back(std::move(s));
            if (solutions.size() >= cfg.max_solutions) break;
            continue;
        }

        auto const& active = s.h_done ? s.phase1 : s.phase0;
        for (auto const& [i, j] : phase_active_pairs(active)) {
            JointState next = s;
            next.gates.push_back(PhaseOp::make_cx(i, j));
            ++next.g_cost;

            if (!next.h_done) {
                // Pre-H: CNOT acts on all four matrices — this is the joint benefit.
                next.phase0.apply_cnot(i, j);
                next.output0.apply_cnot(i, j);
                next.phase1.apply_cnot(i, j);
                next.output1.apply_cnot(i, j);
                emit_ready(next, false);
                finalize_block0(next);
            } else {
                // Post-H: CNOT acts on block-1 matrices only.
                next.phase1.apply_cnot(i, j);
                next.output1.apply_cnot(i, j);
                emit_ready(next, true);
            }

            push(std::move(next));
        }

        while (open.size() > cfg.max_queue_size) open.erase(std::prev(open.end()));
    }

    if (solutions.empty()) {
        // Budget exhausted without finding any solution: fall back to independent.
        GroupResult r;
        r.num_cx     = synthesize_phasepoly(p0, cfg).num_cx +
                       synthesize_phasepoly(p1, cfg).num_cx;
        r.used_fallback = true;
        return r;
    }

    JointState const* best = &solutions.front();
    for (auto const& sol : solutions)
        if (sol.g_cost < best->g_cost) best = &sol;

    return GroupResult{best->g_cost, false};
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

/// Synthesize a group of consecutive blocks using joint A* on consecutive pairs.
///
/// For a group of k blocks [B0, B1, B2, B3, ...]:
///   - Pairs (B0,B1), (B2,B3), ... are synthesized jointly: a single CNOT can
///     simultaneously reduce parity columns from both blocks, sharing the gate cost.
///   - A trailing singleton (if k is odd) is synthesized independently.
///
/// `synthesize_grouped` partitions the global block list into non-overlapping
/// groups of `group_size` and calls this function for each group, so:
///   group_size=1 → per-block independent A* (same as Stage 4)
///   group_size=2 → all pairs jointly synthesized
///   group_size=3 → triples: first pair jointly, third block independently
///   group_size=5 → quintuples: two pairs jointly, fifth block independently
GroupResult synthesize_block_group(
    std::vector<PhaseBlock> const& blocks,
    std::vector<BlockBoundary> const& boundaries,
    PhasePolyConfig const& cfg) {
    if (blocks.empty()) return {};

    GroupResult result;

    for (size_t i = 0; i + 1 < blocks.size(); i += 2) {
        auto const p0 = phase_block_to_problem(blocks[i]);
        auto const p1 = phase_block_to_problem(blocks[i + 1]);
        // boundaries[i] = boundary between blocks[i] and blocks[i+1]
        BlockBoundary const& bnd =
            (i < boundaries.size()) ? boundaries[i] : BlockBoundary{};
        auto const gr = synthesize_joint(p0, p1, bnd.h_qubits, cfg);
        result.num_cx += gr.num_cx;
        if (gr.used_fallback) result.used_fallback = true;
    }

    // Trailing singleton (odd group size).
    if (blocks.size() % 2 == 1) {
        auto const p  = phase_block_to_problem(blocks.back());
        auto const sr = synthesize_phasepoly(p, cfg);
        result.num_cx += sr.num_cx;
        if (sr.used_fallback) result.used_fallback = true;
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
