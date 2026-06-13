/****************************************************************************
  PackageName  [ tableau/phasepoly ]
  Synopsis     [ GF(2) parity matrix used by the PhasePoly co-optimizer ]
  Author       [ Design Verification Lab ]
  Copyright    [ Copyright(c) 2025 DVLab, GIEE, NTU, Taiwan ]
****************************************************************************/

#include "./parity_matrix.hpp"

#include <cassert>

namespace qsyn::experimental::phasepoly {

ParityMatrix::ParityMatrix(size_t n_rows, size_t n_cols)
    : _n_cols{n_cols}, _rows(n_rows, sul::dynamic_bitset<>(n_cols)) {}

ParityMatrix ParityMatrix::identity(size_t n) {
    ParityMatrix m(n, n);
    for (size_t i = 0; i < n; ++i) {
        m.set(i, i, true);
    }
    return m;
}

ParityMatrix ParityMatrix::from_rows(std::vector<std::vector<bool>> const& rows) {
    size_t const n_rows = rows.size();
    size_t const n_cols = n_rows == 0 ? 0 : rows.front().size();
    ParityMatrix m(n_rows, n_cols);
    for (size_t r = 0; r < n_rows; ++r) {
        assert(rows[r].size() == n_cols);
        for (size_t c = 0; c < n_cols; ++c) {
            if (rows[r][c]) m.set(r, c, true);
        }
    }
    return m;
}

ParityMatrix ParityMatrix::from_columns(size_t n_rows, std::vector<sul::dynamic_bitset<>> const& columns) {
    ParityMatrix m(n_rows, columns.size());
    for (size_t c = 0; c < columns.size(); ++c) {
        assert(columns[c].size() == n_rows);
        for (size_t r = 0; r < n_rows; ++r) {
            if (columns[c].test(r)) m.set(r, c, true);
        }
    }
    return m;
}

void ParityMatrix::set(size_t row, size_t col, bool value) {
    assert(row < _rows.size());
    assert(col < _n_cols);
    _rows[row].set(col, value);
}

sul::dynamic_bitset<> ParityMatrix::column(size_t col) const {
    assert(col < _n_cols);
    sul::dynamic_bitset<> result(_rows.size());
    for (size_t r = 0; r < _rows.size(); ++r) {
        if (_rows[r].test(col)) result.set(r);
    }
    return result;
}

void ParityMatrix::apply_cnot(size_t control, size_t target) {
    assert(control < _rows.size());
    assert(target < _rows.size());
    if (control == target) return;
    // Paper Eq. (5): the control row absorbs the target row.
    _rows[control] ^= _rows[target];
}

size_t ParityMatrix::column_weight(size_t col) const {
    assert(col < _n_cols);
    size_t weight = 0;
    for (auto const& r : _rows) {
        if (r.test(col)) ++weight;
    }
    return weight;
}

std::optional<size_t> ParityMatrix::single_one_row(size_t col) const {
    assert(col < _n_cols);
    std::optional<size_t> found;
    for (size_t r = 0; r < _rows.size(); ++r) {
        if (_rows[r].test(col)) {
            if (found.has_value()) return std::nullopt;  // more than one 1
            found = r;
        }
    }
    return found;  // nullopt if the column is all-zero, else the unique row
}

bool ParityMatrix::is_square_identity() const {
    if (_rows.size() != _n_cols) return false;
    for (size_t r = 0; r < _rows.size(); ++r) {
        if (_rows[r].count() != 1 || !_rows[r].test(r)) return false;
    }
    return true;
}

void ParityMatrix::remove_column(size_t col) {
    assert(col < _n_cols);
    for (auto& r : _rows) {
        sul::dynamic_bitset<> shrunk(_n_cols - 1);
        for (size_t c = 0; c < _n_cols; ++c) {
            if (c == col) continue;
            if (r.test(c)) shrunk.set(c < col ? c : c - 1);
        }
        r = std::move(shrunk);
    }
    --_n_cols;
}

std::string ParityMatrix::to_string() const {
    std::string out;
    for (auto const& r : _rows) {
        for (size_t c = 0; c < _n_cols; ++c) {
            out += r.test(c) ? '1' : '0';
        }
        out += '\n';
    }
    return out;
}

}  // namespace qsyn::experimental::phasepoly
