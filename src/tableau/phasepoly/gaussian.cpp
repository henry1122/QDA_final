/****************************************************************************
  PackageName  [ tableau/phasepoly ]
  Synopsis     [ CNOT-only linear-reversible synthesis (output-matrix finish) ]
  Author       [ Design Verification Lab ]
  Copyright    [ Copyright(c) 2025 DVLab, GIEE, NTU, Taiwan ]
****************************************************************************/

#include "./gaussian.hpp"

#include <cassert>
#include <cstddef>
#include <limits>
#include <optional>
#include <vector>

#include "util/boolean_matrix.hpp"

namespace qsyn::experimental::phasepoly {

namespace {

/**
 * @brief Convention-native Gauss-Jordan reduction of `m` to identity.
 *
 * Processes columns left to right, keeping each finished column a unit vector
 * on the diagonal (no row swaps, so the result is exactly the identity). When a
 * diagonal pivot is missing, it is seeded from the minimum-weight eligible row
 * below it to limit fill-in. Records each operation as the CNOT it emits:
 * `apply_cnot(c, t)` -> `CnotOp{control = c, target = t}`.
 */
std::vector<CnotOp> gauss_jordan_reduce(ParityMatrix m) {
    size_t const n = m.n_rows();
    std::vector<CnotOp> ops;

    for (size_t col = 0; col < n; ++col) {
        if (!m.get(col, col)) {
            // Seed the pivot from the sparsest row > col carrying a 1 here.
            // (Rows < col are already unit vectors with a 0 in this column, so
            //  using them would corrupt an earlier column.)
            std::optional<size_t> pivot;
            size_t best_weight = std::numeric_limits<size_t>::max();
            for (size_t r = col + 1; r < n; ++r) {
                if (m.get(r, col) && m.row(r).count() < best_weight) {
                    best_weight = m.row(r).count();
                    pivot       = r;
                }
            }
            assert(pivot.has_value() && "output matrix must be invertible");
            m.apply_cnot(col, *pivot);  // row[col] ^= row[pivot]
            ops.push_back({col, *pivot});
        }
        for (size_t r = 0; r < n; ++r) {
            if (r != col && m.get(r, col)) {
                m.apply_cnot(r, col);  // row[r] ^= row[col]
                ops.push_back({r, col});
            }
        }
    }
    return ops;
}

/// @brief Copy a ParityMatrix into a dvlab::BooleanMatrix.
dvlab::BooleanMatrix to_boolean_matrix(ParityMatrix const& matrix) {
    std::vector<dvlab::BooleanMatrix::Row> rows;
    rows.reserve(matrix.n_rows());
    for (size_t r = 0; r < matrix.n_rows(); ++r) {
        std::vector<unsigned char> bits(matrix.n_cols(), 0);
        for (size_t c = 0; c < matrix.n_cols(); ++c) {
            if (matrix.get(r, c)) bits[c] = 1;
        }
        rows.emplace_back(bits);
    }
    return dvlab::BooleanMatrix{std::move(rows)};
}

/**
 * @brief Patel-Markov-Hayes block elimination via dvlab::BooleanMatrix.
 *
 * `BooleanMatrix::row_operation(ctrl, targ)` performs `row[targ] ^= row[ctrl]`,
 * the mirror of our `apply_cnot`. So a recorded op (ctrl, targ) is emitted as
 * the physical CNOT(control = targ, target = ctrl), preserving order. We search
 * block sizes and keep the sequence with the fewest CNOTs.
 */
std::vector<CnotOp> patel_markov_hayes_reduce(ParityMatrix const& matrix) {
    size_t const n          = matrix.n_rows();
    auto const base         = to_boolean_matrix(matrix);
    std::vector<dvlab::BooleanMatrix::RowOperation> best;
    bool have_best = false;

    for (size_t block_size = 1; block_size <= n; ++block_size) {
        auto candidate = base;
        candidate.gaussian_elimination_skip(block_size, /*do_fully_reduced=*/true, /*track=*/true);
        auto const& row_ops = candidate.get_row_operations();
        if (!have_best || row_ops.size() < best.size()) {
            best      = row_ops;
            have_best = true;
        }
    }

    std::vector<CnotOp> ops;
    ops.reserve(best.size());
    for (auto const& [ctrl, targ] : best) {
        ops.push_back({/*control=*/targ, /*target=*/ctrl});
    }
    return ops;
}

}  // namespace

std::vector<CnotOp> synthesize_linear_reversible(ParityMatrix const& matrix, LinearSynthesisMode mode) {
    assert(matrix.n_rows() == matrix.n_cols() && "output matrix must be square");
    if (matrix.n_rows() == 0) return {};
    switch (mode) {
        case LinearSynthesisMode::gauss_jordan:
            return gauss_jordan_reduce(matrix);
        case LinearSynthesisMode::patel_markov_hayes:
            return patel_markov_hayes_reduce(matrix);
    }
    return {};  // unreachable
}

size_t linear_reversible_cnot_cost(ParityMatrix const& matrix, LinearSynthesisMode mode) {
    return synthesize_linear_reversible(matrix, mode).size();
}

}  // namespace qsyn::experimental::phasepoly
