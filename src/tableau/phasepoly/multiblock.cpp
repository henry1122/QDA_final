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

}  // namespace qsyn::experimental::phasepoly
