/****************************************************************************
  PackageName  [ tableau/phasepoly ]
  Synopsis     [ Multi-block SSA merging for QCir-level optimization (§3.3) ]
  Author       [ Design Verification Lab ]
  Copyright    [ Copyright(c) 2025 DVLab, GIEE, NTU, Taiwan ]
****************************************************************************/

#include "./multiblock.hpp"

#include "./qcir_optimizer.hpp"

#include <algorithm>
#include <numeric>
#include <optional>
#include <utility>
#include <vector>

#include "./extractor.hpp"
#include "./search.hpp"
#include "./synthesizer.hpp"
#include "./todd_preprocess.hpp"
#include "qcir/basic_gate_type.hpp"
#include "qcir/operation.hpp"
#include "qcir/qcir.hpp"
#include "qcir/qcir_gate.hpp"
#include "qsyn/qsyn_type.hpp"

#include <cstddef>
#include <cstdint>
#include <limits>
#include <map>
#include <string>
#include <unordered_map>

#include "./gaussian.hpp"
#include "./parity_matrix.hpp"

namespace qsyn::experimental::phasepoly {

namespace {

bool is_phase_op(qcir::Operation const& op, QubitIdList const& qubits, PhaseBlock& block) {
    if (auto const control_gate = op.get_underlying_if<qcir::ControlGate>()) {
        if (control_gate->get_num_ctrls() != 1) return false;
        auto const target = control_gate->get_target_operation().get_underlying_if<qcir::PXGate>();
        if (!target || target->get_phase() != dvlab::Phase(1)) return false;
        if (qubits[0] == qubits[1]) return false;
        block.append_cx(static_cast<size_t>(qubits[0]), static_cast<size_t>(qubits[1]));
        return true;
    }
    if (auto const pz = op.get_underlying_if<qcir::PZGate>()) {
        block.append_rz(static_cast<size_t>(qubits[0]), pz->get_phase());
        return true;
    }
    if (auto const rz = op.get_underlying_if<qcir::RZGate>()) {
        block.append_rz(static_cast<size_t>(qubits[0]), rz->get_phase());
        return true;
    }
    return false;
}

bool is_ssa_boundary(qcir::Operation const& op) {
    return op.get_underlying_if<qcir::HGate>().has_value();
}

size_t count_cx(qcir::QCir const& circuit) {
    size_t n = 0;
    for (auto const* gate : circuit.get_gates()) {
        if (gate->get_operation().get_underlying_if<qcir::ControlGate>()) ++n;
    }
    return n;
}

size_t count_cx_rz(qcir::QCir const& circuit) {
    size_t n = 0;
    for (auto const* gate : circuit.get_gates()) {
        auto const& op = gate->get_operation();
        if (op.get_underlying_if<qcir::ControlGate>()) {
            ++n;
        } else if (op.get_underlying_if<qcir::PZGate>() || op.get_underlying_if<qcir::RZGate>()) {
            ++n;
        }
    }
    return n;
}

class SSAContext {
public:
    explicit SSAContext(size_t n_logical)
        : _row_map(n_logical) {
        std::iota(_row_map.begin(), _row_map.end(), 0);
        _row_to_logical.resize(n_logical);
        std::iota(_row_to_logical.begin(), _row_to_logical.end(), 0);
    }

    size_t num_rows() const { return _row_to_logical.size(); }
    std::vector<size_t> const& row_to_logical() const { return _row_to_logical; }

    size_t remap(size_t logical_q) const { return _row_map.at(logical_q); }

    void h_gate(size_t logical_q) {
        _row_to_logical.push_back(logical_q);
        _row_map.at(logical_q) = _row_to_logical.size() - 1;
    }

private:
    std::vector<size_t> _row_map;
    std::vector<size_t> _row_to_logical;
};

void append_remapped(std::vector<PhaseOp>& ops, PhaseBlock const& block, SSAContext const& ctx) {
    for (auto const& op : block.ops()) {
        if (op.is_cx()) {
            ops.push_back(PhaseOp::make_cx(ctx.remap(op.control), ctx.remap(op.target)));
        } else {
            ops.push_back(PhaseOp::make_rz(ctx.remap(op.qubit()), op.phase));
        }
    }
}

qcir::QCir build_qcir_mapped(
    size_t n_logical,
    std::vector<PhaseOp> const& gates,
    std::vector<size_t> const& row_to_logical) {
    qcir::QCir circuit{n_logical};
    for (auto const& op : gates) {
        if (op.is_cx()) {
            auto const ctrl = row_to_logical.at(op.control);
            auto const targ = row_to_logical.at(op.target);
            if (ctrl == targ) continue;
            circuit.append(
                qcir::CXGate(),
                QubitIdList{
                    static_cast<QubitIdType>(ctrl),
                    static_cast<QubitIdType>(targ)});
        } else {
            circuit.append(
                qcir::PZGate(op.phase),
                QubitIdList{static_cast<QubitIdType>(row_to_logical.at(op.qubit()))});
        }
    }
    return circuit;
}

std::optional<qcir::QCir> synthesize_segment_group(
    std::vector<PhasePolySegment> const& segments,
    size_t start,
    size_t k,
    PhasePolyConfig const& config,
    bool replace_only_if_better) {
    auto merged_opt = merge_segments_ssa(segments, start, k);
    if (!merged_opt || merged_opt->empty()) return std::nullopt;

    PhaseBlock const& merged = *merged_opt;

    size_t cx_before = 0;
    for (size_t j = start; j < start + k; ++j) {
        cx_before += segments.at(j).block.num_cx();
    }

    auto const row_map = ssa_row_to_logical_map(segments, start, k);
    auto const fallback = [&]() {
        return build_qcir_mapped(segments.front().block.n_qubits(), merged.ops(), row_map);
    };

    auto problem = phase_block_to_problem(merged);
    if (config.apply_block_todd) {
        if (auto const opt = todd_optimize_problem(problem)) {
            problem = *opt;
        }
    }
    if (problem.is_trivial()) return fallback();

    auto const result = synthesize_phasepoly(problem, config);
    if (!verify_synthesis(problem, result)) return fallback();
    if (replace_only_if_better && result.num_cx >= cx_before) return fallback();

    return build_qcir_mapped(segments.front().block.n_qubits(), result.gates, row_map);
}

}  // namespace

std::vector<PhasePolySegment> segment_phase_poly_regions(qcir::QCir const& circuit) {
    size_t const n = circuit.get_num_qubits();
    std::vector<PhasePolySegment> segments;
    PhaseBlock current{n};

    auto flush_block = [&]() {
        if (current.empty()) return;
        segments.push_back(PhasePolySegment{std::move(current), {}});
        current = PhaseBlock{n};
    };

    for (auto const* gate : circuit.get_gates()) {
        if (is_phase_op(gate->get_operation(), gate->get_qubits(), current)) {
            continue;
        }
        flush_block();
        BoundaryGate boundary{gate->get_operation(), gate->get_qubits()};
        if (segments.empty()) {
            segments.push_back(PhasePolySegment{PhaseBlock{n}, {boundary}});
        } else {
            segments.back().trailing.push_back(boundary);
        }
    }
    flush_block();

    segments.erase(
        std::remove_if(
            segments.begin(),
            segments.end(),
            [](PhasePolySegment const& seg) { return seg.block.empty() && seg.trailing.empty(); }),
        segments.end());

    return segments;
}

std::optional<PhaseBlock> merge_segments_ssa(
    std::vector<PhasePolySegment> const& segments,
    size_t start,
    size_t k) {
    if (k == 0 || start + k > segments.size()) return std::nullopt;

    size_t const n_logical = segments.at(start).block.n_qubits();
    if (segments.at(start).block.empty()) return std::nullopt;

    SSAContext ctx{n_logical};
    std::vector<PhaseOp> ops;

    append_remapped(ops, segments.at(start).block, ctx);

    for (size_t j = start; j < start + k - 1; ++j) {
        if (segments.at(j).block.empty()) return std::nullopt;
        for (auto const& boundary : segments.at(j).trailing) {
            if (!is_ssa_boundary(boundary.op) || boundary.qubits.size() != 1) {
                return std::nullopt;
            }
            ctx.h_gate(static_cast<size_t>(boundary.qubits.front()));
        }
        if (segments.at(j + 1).block.empty()) return std::nullopt;
        append_remapped(ops, segments.at(j + 1).block, ctx);
    }

    return PhaseBlock::from_ops(ctx.num_rows(), std::move(ops));
}

std::vector<size_t> ssa_row_to_logical_map(
    std::vector<PhasePolySegment> const& segments,
    size_t start,
    size_t k) {
    size_t const n_logical = segments.at(start).block.n_qubits();
    SSAContext ctx{n_logical};
    for (size_t j = start; j < start + k - 1; ++j) {
        for (auto const& boundary : segments.at(j).trailing) {
            if (boundary.qubits.size() == 1) {
                ctx.h_gate(static_cast<size_t>(boundary.qubits.front()));
            }
        }
    }
    return ctx.row_to_logical();
}

std::optional<qcir::QCir> optimize_qcir_phasepoly_multiblock(
    qcir::QCir const& circuit,
    PhasePolyConfig const& config,
    bool replace_only_if_better,
    PhasePolyOptimizeStats* stats,
    size_t* phase_gate_out) {
    auto segments = segment_phase_poly_regions(circuit);
    if (segments.empty()) return circuit;

    size_t const n = circuit.get_num_qubits();
    qcir::QCir result{n};
    PhasePolyOptimizeStats local{};

    size_t i = 0;
    while (i < segments.size()) {
        if (segments.at(i).block.empty()) {
            for (auto const& boundary : segments.at(i).trailing) {
                result.append(boundary.op, boundary.qubits);
            }
            ++i;
            continue;
        }

        qcir::QCir chosen{n};
        size_t chosen_span = 1;
        size_t chosen_cx   = segments.at(i).block.num_cx();

        auto const try_group = [&](size_t k) {
            if (i + k > segments.size()) return;
            if (!merge_segments_ssa(segments, i, k)) return;
            auto synth = synthesize_segment_group(segments, i, k, config, replace_only_if_better);
            if (!synth) return;
            size_t const cx = count_cx(*synth);
            if (chosen.get_num_gates() == 0 || cx < chosen_cx) {
                chosen      = std::move(*synth);
                chosen_span = k;
                chosen_cx   = cx;
            }
        };

        for (size_t k : config.group_sizes) {
            try_group(k);
        }

        if (chosen.get_num_gates() == 0) {
            auto const row_map = ssa_row_to_logical_map(segments, i, 1);
            chosen             = build_qcir_mapped(n, segments.at(i).block.ops(), row_map);
            chosen_span        = 1;
            chosen_cx          = segments.at(i).block.num_cx();
        }

        local.num_blocks += chosen_span;
        size_t cx_orig = 0;
        for (size_t j = i; j < i + chosen_span; ++j) {
            cx_orig += segments.at(j).block.num_cx();
        }
        local.cx_before += cx_orig;
        local.cx_after += chosen_cx;
        if (chosen_cx < cx_orig) local.blocks_optimized += chosen_span;

        if (phase_gate_out) *phase_gate_out += count_cx_rz(chosen);
        result.compose(chosen);

        for (size_t j = i; j < i + chosen_span - 1; ++j) {
            for (auto const& boundary : segments.at(j).trailing) {
                result.append(boundary.op, boundary.qubits);
            }
        }
        for (auto const& boundary : segments.at(i + chosen_span - 1).trailing) {
            result.append(boundary.op, boundary.qubits);
        }

        i += chosen_span;
    }

    if (stats) *stats = local;
    return result;
}

// ── Warm-start group synthesis (extract + synthesize_block_group + synthesize_grouped) ──


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
            auto const sr0 = synthesize_phasepoly(p0, cfg);
            auto const sr1 = synthesize_phasepoly(p1, cfg);
            GroupResult r;
            r.num_cx        = sr0.num_cx + sr1.num_cx;
            r.num_rz        = sr0.num_rz + sr1.num_rz;
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
        auto const sr0 = synthesize_phasepoly(p0, cfg);
        auto const sr1 = synthesize_phasepoly(p1, cfg);
        GroupResult r;
        r.num_cx        = sr0.num_cx + sr1.num_cx;
        r.num_rz        = sr0.num_rz + sr1.num_rz;
        r.used_fallback = true;
        return r;
    }

    JointState const* best = &solutions.front();
    for (auto const& sol : solutions)
        if (sol.g_cost < best->g_cost) best = &sol;

    size_t rz_count = 0;
    for (auto const& op : best->gates)
        if (op.is_rz()) ++rz_count;

    return GroupResult{best->g_cost, rz_count, false};
}

// ─── SSA rename → merge → single A* (paper §3.3) ──────────────────────────

/// Merge k consecutive blocks (blocks[start..start+k-1]) into a single PhaseBlock
/// using SSA renaming at H-boundaries. Returns nullopt if any boundary contains
/// a non-H gate (unsafe to merge) or the index range is out of bounds.
std::optional<PhaseBlock> merge_blocks_ssa(
    std::vector<PhaseBlock> const& blocks,
    std::vector<BlockBoundary> const& boundaries,
    size_t start,
    size_t k) {
    if (k == 0 || start + k > blocks.size()) return std::nullopt;
    if (k == 1) return blocks[start];

    size_t const n_logical = blocks[start].n_qubits();
    SSAContext ctx{n_logical};
    std::vector<PhaseOp> ops;

    auto append_remapped = [&](PhaseBlock const& block) {
        for (auto const& op : block.ops()) {
            if (op.is_cx())
                ops.push_back(PhaseOp::make_cx(ctx.remap(op.control), ctx.remap(op.target)));
            else
                ops.push_back(PhaseOp::make_rz(ctx.remap(op.qubit()), op.phase));
        }
    };

    append_remapped(blocks[start]);

    for (size_t j = 0; j + 1 < k; ++j) {
        size_t const bnd_idx = start + j;
        if (bnd_idx >= boundaries.size()) return std::nullopt;
        BlockBoundary const& bnd = boundaries[bnd_idx];
        if (bnd.has_non_h_gate) return std::nullopt;

        for (size_t h_qubit : bnd.h_qubits)
            ctx.h_gate(h_qubit);

        append_remapped(blocks[start + j + 1]);
    }

    PhaseBlock merged{ctx.num_rows()};
    for (auto const& op : ops) {
        if (op.is_cx()) merged.append_cx(op.control, op.target);
        else            merged.append_rz(op.qubit(), op.phase);
    }
    return merged;
}

/// Synthesize a group of k consecutive blocks by merging them via SSA renaming
/// and running a single A* on the merged problem. Falls back to per-block
/// independent synthesis if SSA merge is infeasible (non-H boundary gate).
GroupResult synthesize_ssa_group(
    std::vector<PhaseBlock> const& blocks,
    std::vector<BlockBoundary> const& boundaries,
    PhasePolyConfig const& cfg) {
    if (blocks.empty()) return {};

    auto merged_opt = merge_blocks_ssa(blocks, boundaries, 0, blocks.size());
    if (merged_opt) {
        auto const p = phase_block_to_problem(*merged_opt);

        // Scale budget proportionally to the merged block's row count.
        // The merged problem has more rows (n_logical + H_extras) and more
        // phase columns (k blocks combined), so the same budget explores less.
        PhasePolyConfig scaled_cfg = cfg;
        if (cfg.scale_budget_ssa && blocks.size() > 1) {
            size_t const orig_rows   = blocks[0].n_qubits();
            size_t const merged_rows = merged_opt->n_qubits();
            size_t const k           = blocks.size();
            // Scale by max(k, merged/original) capped at 4× to avoid runaway.
            size_t const row_scale   = (merged_rows + orig_rows - 1) / orig_rows;
            size_t const scale       = std::min(std::max(k, row_scale), size_t{4});
            scaled_cfg.max_expansions = cfg.max_expansions * scale;
        }

        auto const sr = synthesize_phasepoly(p, scaled_cfg);
        return GroupResult{sr.num_cx, sr.num_rz, sr.used_fallback};
    }

    // Non-H boundary prevents SSA merge: fall back to independent per-block A*.
    GroupResult result;
    result.used_fallback = true;
    for (auto const& block : blocks) {
        auto const p  = phase_block_to_problem(block);
        auto const sr = synthesize_phasepoly(p, cfg);
        result.num_cx += sr.num_cx;
        result.num_rz += sr.num_rz;
        if (sr.used_fallback) result.used_fallback = true;
    }
    return result;
}

}  // namespace

// ─── Public API ───────────────────────────────────────────────────────────

ExtractedCircuit extract_phase_blocks_with_boundaries(qcir::QCir const& circuit) {
    size_t const n = circuit.get_num_qubits();
    ExtractedCircuit result;
    result.n_qubits = n;

    PhaseBlock current(n);
    std::vector<size_t> pending_h_qubits;
    bool pending_has_non_h = false;
    bool has_flushed       = false;

    auto flush = [&]() {
        if (!current.empty()) {
            if (has_flushed) {
                result.boundaries.push_back(
                    BlockBoundary{pending_h_qubits, pending_has_non_h});
                pending_h_qubits.clear();
                pending_has_non_h = false;
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
            else
                pending_has_non_h = true;
        }
    }
    flush();

    return result;
}

/// Synthesize a group of consecutive blocks using the strategy in `cfg`.
///
/// SSA merge (default, paper §3.3):
///   All k blocks are SSA-renamed and merged into a single problem, then solved
///   by one A* call. Falls back to per-block independent synthesis if any
///   boundary contains a non-H gate.
///
/// Joint A* (experimental):
///   Consecutive pairs (B0,B1), (B2,B3), ... are synthesized jointly; a single
///   CNOT can reduce parities from both blocks simultaneously. A trailing
///   singleton (odd k) is synthesized independently.
GroupResult synthesize_block_group(
    std::vector<PhaseBlock> const& blocks,
    std::vector<BlockBoundary> const& boundaries,
    PhasePolyConfig const& cfg) {
    if (blocks.empty()) return {};

    if (cfg.multi_block_strategy == MultiBlockStrategy::ssa_merge) {
        return synthesize_ssa_group(blocks, boundaries, cfg);
    }

    // ── Joint A* on consecutive pairs ─────────────────────────────────────
    GroupResult result;

    for (size_t i = 0; i + 1 < blocks.size(); i += 2) {
        auto const p0 = phase_block_to_problem(blocks[i]);
        auto const p1 = phase_block_to_problem(blocks[i + 1]);
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
        result.num_rz += sr.num_rz;
        if (sr.used_fallback) result.used_fallback = true;
    }

    return result;
}

GroupResult synthesize_grouped(ExtractedCircuit const& extracted,
                               size_t group_size,
                               PhasePolyConfig const& cfg) {
    if (group_size == 0) group_size = 1;

    auto const& blocks     = extracted.blocks;
    auto const& boundaries = extracted.boundaries;
    size_t const nb        = blocks.size();
    GroupResult total;

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
        total.num_cx += gr.num_cx;
        total.num_rz += gr.num_rz;
        if (gr.used_fallback) total.used_fallback = true;
    }

    return total;
}

}  // namespace qsyn::experimental::phasepoly
