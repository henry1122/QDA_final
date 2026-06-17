/****************************************************************************
  PackageName  [ tableau/phasepoly ]
  Synopsis     [ Memory-bounded A* search over the joint [P | O] matrix ]
  Author       [ Design Verification Lab ]
  Copyright    [ Copyright(c) 2025 DVLab, GIEE, NTU, Taiwan ]
****************************************************************************/

#include "./search.hpp"

#include <cstddef>
#include <cstdint>
#include <limits>
#include <map>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

#include "./gaussian.hpp"
#include "./parity_matrix.hpp"

namespace qsyn::experimental::phasepoly {

namespace {

/// @brief A node in the search tree: the current matrices plus the gates so far.
struct SearchState {
    ParityMatrix phase;                // P
    std::vector<dvlab::Phase> angles;  // Theta (parallel to phase columns)
    ParityMatrix output;               // O
    std::vector<PhaseOp> gates;        // emitted CNOTs / Rz, in order
    size_t g_cost = 0;                 // number of CNOTs emitted
};

/// @brief h1: sum of Hamming weights of the remaining phase columns (paper §3.2.2).
size_t phase_cost(SearchState const& s) {
    size_t sum = 0;
    for (size_t c = 0; c < s.phase.n_cols(); ++c) sum += s.phase.column_weight(c);
    return sum;
}

/**
 * @brief Emit an `Rz` for every weight-1 phase column and remove it.
 *
 * When a column becomes `e_q`, qubit `q` holds exactly that column's parity, so
 * the rotation can be applied there (proven in the Stage-3 notes). Repeats until
 * no weight-1 column remains, since removing one may not expose others but the
 * scan is cheap.
 */
void remove_ready_phase_columns(SearchState& s) {
    bool changed = true;
    while (changed) {
        changed = false;
        for (size_t c = 0; c < s.phase.n_cols(); ++c) {
            if (auto const row = s.phase.single_one_row(c)) {
                s.gates.push_back(PhaseOp::make_rz(*row, s.angles[c]));
                s.phase.remove_column(c);
                s.angles.erase(s.angles.begin() + static_cast<std::ptrdiff_t>(c));
                changed = true;
                break;  // indices shifted; restart the scan
            }
        }
    }
}

/**
 * @brief Active row pairs (i, j): applying `apply_cnot(i, j)` reduces the
 *        Hamming weight of at least one phase column -- equivalently, rows i and
 *        j both carry a 1 in some column (paper §3.2.1).
 */
std::vector<std::pair<size_t, size_t>> active_row_pairs(ParityMatrix const& phase) {
    size_t const n = phase.n_rows();
    std::vector<std::pair<size_t, size_t>> pairs;
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            if (i == j) continue;
            for (size_t c = 0; c < phase.n_cols(); ++c) {
                if (phase.get(i, c) && phase.get(j, c)) {
                    pairs.emplace_back(i, j);
                    break;
                }
            }
        }
    }
    return pairs;
}

/// @brief A canonical key over (P with angles, O) for the visited cache.
std::string state_key(SearchState const& s) {
    std::string key;
    key.reserve(s.phase.n_rows() * s.phase.n_cols() + s.output.n_rows() * s.output.n_cols());
    for (size_t c = 0; c < s.phase.n_cols(); ++c) {
        for (size_t r = 0; r < s.phase.n_rows(); ++r) key += s.phase.get(r, c) ? '1' : '0';
        key += ':';
        key += s.angles[c].get_print_string();
        key += ';';
    }
    key += '#';
    for (size_t r = 0; r < s.output.n_rows(); ++r) {
        for (size_t c = 0; c < s.output.n_cols(); ++c) key += s.output.get(r, c) ? '1' : '0';
    }
    return key;
}

/// @brief Append the Gaussian-elimination finish that drives O to identity.
void finish_output(SearchState& s, LinearSynthesisMode mode) {
    for (auto const& op : synthesize_linear_reversible(s.output, mode)) {
        s.output.apply_cnot(op.control, op.target);
        s.gates.push_back(PhaseOp::make_cx(op.control, op.target));
        ++s.g_cost;
    }
}

/// @brief Guaranteed-correct greedy synthesis, used when the search budget runs out.
SearchState greedy_synthesize(PhasePolyProblem const& problem, LinearSynthesisMode finish_mode) {
    SearchState s{problem.phase_matrix, problem.phase_angles, problem.output_matrix, {}, 0};
    remove_ready_phase_columns(s);
    // Complete one column at a time by folding pairs of set rows together.
    while (!s.phase.has_no_columns()) {
        size_t const col = 0;
        while (s.phase.column_weight(col) > 1) {
            size_t first = s.phase.n_rows(), second = s.phase.n_rows();
            for (size_t r = 0; r < s.phase.n_rows(); ++r) {
                if (s.phase.get(r, col)) {
                    if (first == s.phase.n_rows())
                        first = r;
                    else {
                        second = r;
                        break;
                    }
                }
            }
            s.phase.apply_cnot(first, second);
            s.output.apply_cnot(first, second);
            s.gates.push_back(PhaseOp::make_cx(first, second));
            ++s.g_cost;
        }
        remove_ready_phase_columns(s);
    }
    finish_output(s, finish_mode);
    return s;
}

SynthesisResult to_result(SearchState const& s, PhasePolyProblem const& problem) {
    SynthesisResult result;
    result.n_qubits = problem.n_qubits;
    result.gates    = s.gates;
    result.num_cx   = s.g_cost;
    result.num_rz   = problem.num_phase_terms();
    return result;
}

}  // namespace

SynthesisResult greedy_synthesize_problem(PhasePolyProblem const& problem, PhasePolyConfig const& config) {
    return to_result(greedy_synthesize(problem, config.finish_mode), problem);
}

SynthesisResult synthesize_phasepoly(PhasePolyProblem const& problem, PhasePolyConfig const& config) {
    using Priority = std::tuple<size_t, size_t, size_t, size_t, uint64_t>;  // f, h1, h2, (max-g), counter

    SearchState root{problem.phase_matrix, problem.phase_angles, problem.output_matrix, {}, 0};
    remove_ready_phase_columns(root);

    std::map<Priority, SearchState> open;
    std::unordered_map<std::string, size_t> best_g;
    uint64_t counter = 0;

    auto const priority_of = [&](SearchState const& s) -> Priority {
        size_t const h1 = phase_cost(s);
        size_t const h2 = linear_reversible_cnot_cost(s.output, config.h2_mode);
        size_t const f  = s.g_cost + h1 + h2;
        return {f, h1, h2, std::numeric_limits<size_t>::max() - s.g_cost, counter++};
    };

    auto const push = [&](SearchState s) {
        auto const key = state_key(s);
        if (auto const it = best_g.find(key); it != best_g.end() && it->second <= s.g_cost) return;
        best_g[key] = s.g_cost;
        open.emplace(priority_of(s), std::move(s));
    };

    push(root);

    std::vector<SearchState> solutions;
    size_t expansions = 0;
    size_t peak_queue = open.size();

    while (!open.empty() && expansions < config.max_expansions) {
        auto node      = open.extract(open.begin());
        SearchState s  = std::move(node.mapped());
        ++expansions;

        if (s.phase.has_no_columns()) {
            if (!s.output.is_square_identity()) finish_output(s, config.finish_mode);
            solutions.push_back(std::move(s));
            if (solutions.size() >= config.max_solutions) break;
            continue;
        }

        for (auto const& [i, j] : active_row_pairs(s.phase)) {
            SearchState next = s;
            next.phase.apply_cnot(i, j);
            next.output.apply_cnot(i, j);
            next.gates.push_back(PhaseOp::make_cx(i, j));
            ++next.g_cost;
            remove_ready_phase_columns(next);
            push(std::move(next));
        }

        while (open.size() > config.max_queue_size) open.erase(std::prev(open.end()));
        peak_queue = std::max(peak_queue, open.size());
    }

    if (solutions.empty()) {
        auto result          = to_result(greedy_synthesize(problem, config.finish_mode), problem);
        result.expansions    = expansions;
        result.max_queue     = peak_queue;
        result.num_solutions = 0;
        result.used_fallback = true;
        return result;
    }

    SearchState const* best = &solutions.front();
    for (auto const& sol : solutions) {
        if (sol.g_cost < best->g_cost) best = &sol;
    }

    auto result          = to_result(*best, problem);
    result.expansions    = expansions;
    result.max_queue     = peak_queue;
    result.num_solutions = solutions.size();
    return result;
}

}  // namespace qsyn::experimental::phasepoly
