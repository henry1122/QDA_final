/****************************************************************************
  PackageName  [ tableau/phasepoly ]
  Synopsis     [ GF(2) parity matrix used by the PhasePoly co-optimizer ]
  Author       [ Design Verification Lab ]
  Copyright    [ Copyright(c) 2025 DVLab, GIEE, NTU, Taiwan ]
****************************************************************************/

#pragma once

#include <sul/dynamic_bitset.hpp>

#include <cstddef>
#include <optional>
#include <string>
#include <vector>

namespace qsyn::experimental::phasepoly {

/**
 * @brief A boolean (GF(2)) matrix where each row corresponds to a qubit and
 *        each column corresponds to a parity term over the input variables.
 *
 * This is the data structure behind the PhasePoly joint matrix `[P | O]`
 * (paper §3.1): the phase-parity matrix `P` (one column per `Rz` term) and the
 * output-parity matrix `O` (one column per qubit output) are both instances of
 * `ParityMatrix`. A CNOT is modelled as a row operation applied identically to
 * both matrices.
 *
 * == Row-operation convention (paper Eq. 5) ==
 * A CNOT whose control is row `i` and whose target is row `j` updates the
 * matrix as
 *      row_i <- row_i XOR row_j      (the *control* row changes; row_j is kept)
 * implemented by `apply_cnot(i, j)`.
 *
 * Note this is the *mirror* of `dvlab::BooleanMatrix::row_operation(ctrl, targ)`
 * which performs `row_targ <- row_targ XOR row_ctrl`. It is, however, exactly
 * the Z-support action of `qsyn::experimental::PauliProduct::cx(control, target)`
 * (`z[control] ^= z[target]`). Both relationships are asserted in the Stage-0
 * convention unit tests.
 *
 * Rows are stored as `sul::dynamic_bitset<>` of width `n_cols()`, making a CNOT
 * (whole-row XOR) a word-parallel operation.
 */
class ParityMatrix {
public:
    ParityMatrix() = default;
    ParityMatrix(size_t n_rows, size_t n_cols);

    /// @brief Build the n x n identity matrix (qubit q's output is input q).
    static ParityMatrix identity(size_t n);

    /// @brief Build a matrix from explicit row-major bits (mainly for tests).
    static ParityMatrix from_rows(std::vector<std::vector<bool>> const& rows);

    /// @brief Build a matrix whose j-th column is `columns[j]` (over `n_rows`).
    static ParityMatrix from_columns(size_t n_rows, std::vector<sul::dynamic_bitset<>> const& columns);

    size_t n_rows() const { return _rows.size(); }
    size_t n_cols() const { return _n_cols; }

    /// @brief True when there are no columns left (all phase terms completed).
    bool has_no_columns() const { return _n_cols == 0; }

    bool get(size_t row, size_t col) const { return _rows[row].test(col); }
    void set(size_t row, size_t col, bool value = true);

    sul::dynamic_bitset<> const& row(size_t r) const { return _rows[r]; }

    /// @brief Extract column `col` as a bitset indexed by row.
    sul::dynamic_bitset<> column(size_t col) const;

    /**
     * @brief Apply a CNOT row operation following the paper convention:
     *        row[control] <- row[control] XOR row[target].
     */
    void apply_cnot(size_t control, size_t target);

    /// @brief Number of 1s in a column (its Hamming weight).
    size_t column_weight(size_t col) const;

    /// @brief If column `col` has Hamming weight 1, return the row of that 1.
    std::optional<size_t> single_one_row(size_t col) const;

    /// @brief True iff the matrix is square and equals the identity.
    bool is_square_identity() const;

    /// @brief Remove a (completed) column, shrinking the matrix by one column.
    void remove_column(size_t col);

    bool operator==(ParityMatrix const& rhs) const = default;

    std::string to_string() const;

private:
    size_t _n_cols = 0;
    std::vector<sul::dynamic_bitset<>> _rows;  // each of width `_n_cols`
};

/**
 * @brief Free-function alias matching the recommended helper naming in the
 *        PhasePoly spec. Performs the paper-convention row operation
 *        `row[control_row] <- row[control_row] XOR row[target_row]`.
 */
inline void row_xor(ParityMatrix& matrix, size_t control_row, size_t target_row) {
    matrix.apply_cnot(control_row, target_row);
}

}  // namespace qsyn::experimental::phasepoly
