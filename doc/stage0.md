# PhasePoly — Stage 0 Report: Convention Lock-in & GF(2) Primitives

**Goal of this stage (from `PHASEPOLY_PLAN.md`):** Build the GF(2) matrix primitives, a forward symbolic simulator, and — most importantly — **pin down the row-operation convention** so every later stage (extraction, A\* search, Gaussian finish, emission) shares one verified, unambiguous semantics. This is the single biggest source of silent bugs in a phase-polynomial synthesizer, so it gets locked behind unit tests before any search code exists.

---

## 1. What I implemented

All under a new, self-contained module `qsyn::experimental::phasepoly`. New files auto-compile via the existing `GLOB_RECURSE src/**/*.cpp` in `CMakeLists.txt` (no CMake edits needed).

| File | Contents |
|---|---|
| `src/tableau/phasepoly/parity_matrix.{hpp,cpp}` | `ParityMatrix` — the GF(2) matrix behind the joint `[P\|O]` representation, plus `row_xor` helper |
| `src/tableau/phasepoly/symbolic_state.{hpp,cpp}` | `SymbolicState` — forward (physical) CNOT simulator, plus `apply_cnot_to_state` helper |
| `tests/src/tableau/phasepoly/convention.cpp` | Catch2 tests that lock in the convention (7 test cases) |

### 1.1 `ParityMatrix`

A boolean matrix where **rows = qubits** and **columns = parity terms** over the input variables. It backs both the phase-parity matrix `P` (one column per `Rz`) and the output-parity matrix `O` (one column per qubit output). Rows are stored as `sul::dynamic_bitset<>` (same bitset library `PauliProduct` uses), so a CNOT is a word-parallel whole-row XOR.

Key API:
- `apply_cnot(control, target)` — the **paper-convention row operation** `row[control] ^= row[target]` (Eq. 5). `row_xor(m, c, t)` is a free-function alias matching the spec's recommended helper name.
- `column_weight(col)` — Hamming weight of a column.
- `single_one_row(col)` — if the column has weight 1, the row of that single `1` (i.e. the column is *ready* and an `Rz` can be emitted there); else `nullopt`.
- `is_square_identity()` — the output-matrix goal test.
- `remove_column(col)` — drop a completed phase column.
- Constructors/builders: `identity(n)`, `from_rows(...)`, `from_columns(...)`, plus `get/set/column/to_string` and `operator==`.

### 1.2 `SymbolicState`

Tracks, for each qubit, the parity over input variables it currently holds while scanning a CNOT network forward — i.e. the **physical** circuit semantics used during block extraction.
- Constructed at identity: qubit `q` holds input variable `q`.
- `apply_cnot(control, target)` — physical CNOT `value[target] ^= value[control]` (target qubit changes). `apply_cnot_to_state(s, c, t)` is the spec-named alias.
- `to_output_matrix()` — builds `O` with column `q` = the parity held by qubit `q`.

---

## 2. The convention, locked

This is the crux of Stage 0. The two representations are **dual**, and conflating them is the classic phase-poly bug:

| | XOR direction | "which thing changes" |
|---|---|---|
| **Matrix row op** `ParityMatrix::apply_cnot(c, t)` | `row[c] ^= row[t]` | the **control row** changes |
| **Physical / forward** `SymbolicState::apply_cnot(c, t)` | `value[t] ^= value[c]` | the **target qubit** changes |

Derived relationships (all asserted in tests):

1. **`apply_cnot(c, t)` ≡ `PauliProduct::cx(c, t)` on the Z-support** (`z[control] ^= z[target]`). This anchors our convention to qsyn's already-tested Clifford machinery.
2. **`apply_cnot(c, t)` is the mirror of `dvlab::BooleanMatrix::row_operation(ctrl, targ)`** (which does `row[targ] ^= row[ctrl]`). So if we ever reuse `BooleanMatrix`, we must call `row_operation(ctrl = t, targ = c)`. Documented and tested so the swap is never gotten wrong by accident.
3. **Duality `O = Tᵀ`** where `T` is the linear map realised by the scanned network. Consequence — the **emission rule**:
   > A matrix reduction operation `apply_cnot(c, t)` is emitted as a physical `CNOT(control = c, target = t)`, with the **same indices and the same order**.

   The "inverse" the spec warns about is *only* that the matrix XORs the control row while the physical gate XORs the target qubit value — the index labels and ordering are identical. I derived this by hand and then verified it on random circuits (test #7).

---

## 3. Tests (`tests/src/tableau/phasepoly/convention.cpp`)

Run with `make unit-test` (tag `[phasepoly]`). Seven cases:

1. **`apply_cnot` XORs the control row** — basic semantics; control row absorbs target, target row untouched.
2. **Paper Eq. (4) reproduction** — builds the exact 3×5 joint `[P|O]` matrix from the paper's Fig. 6 example, applies `CNOT(2,1)`, and checks the result equals the paper's printed matrix bit-for-bit. This directly validates the row-op direction against the source of truth.
3. **Matches `PauliProduct::cx` Z-action** — a diagonal `PauliRotation`'s Z-support and the equivalent single-column `ParityMatrix` evolve identically under `cx`/`apply_cnot`.
4. **Mirror of `BooleanMatrix::row_operation`** — confirms relationship #2 above with swapped arguments.
5. **Forward simulator on Fig. 1** — reproduces the known parities (`Rz` fires on `x⊕y` after the first CNOT; qubit 2 carries `x⊕y⊕z`), validating the extraction-direction semantics.
6. **Column utilities** — `column_weight`, `single_one_row` (ready-column detection), `remove_column` with column shifting.
7. **Convention round-trip (the decisive one)** — for 200 random CNOT networks (3–6 qubits, up to ~24 gates): forward-simulate to build `O`, reduce `O` to identity with `apply_cnot` row ops, replay each recorded op as a physical `CNOT(control, target)` in the same order, and assert the reconstructed output matrix equals `O`. This closes the full convention loop (extraction → matrix reduction → emission) in pure GF(2) algebra, independent of any circuit builder.

> Test #7 includes a small **test-only** Gauss-Jordan helper (`reduce_to_identity`) to exercise the round trip. The production-grade minimum-degree Gaussian elimination is Stage 2 — this helper is intentionally simple and lives only in the test file.

**Build & test status: ✅ PASSING.**
- `[phasepoly]` tests: **All tests passed (430 assertions in 7 test cases)** — includes the 200-trial random round-trip.
- Full suite (no regressions): **All tests passed (302,854 assertions in 45 test cases)**.
- One necessary CMake edit: added `sul::dynamic_bitset` to the `unit-test` target's `target_link_libraries_system` (the test target links `libqsyn.a` but didn't inherit the header-only bitset include path that the new headers need). Source files themselves still auto-glob.

---

## 4. Design choices & rationale

- **Reuse `sul::dynamic_bitset<>`** (qsyn's existing bitset) rather than `std::vector<bool>` or `dvlab::BooleanMatrix`: word-parallel row XOR, and it keeps us bit-compatible with `PauliProduct` for the Stage-1 extractor.
- **Did *not* reuse `dvlab::BooleanMatrix` as the core type.** Its `row_operation` uses the opposite convention; wrapping it would invite exactly the confusion Stage 0 exists to prevent. We instead document and test the mapping so `BooleanMatrix`'s Gaussian routines can still be leveraged in Stage 2 via the proven argument swap.
- **Separate `ParityMatrix` (synthesis) from `SymbolicState` (extraction/simulation).** They are duals with opposite XOR direction; keeping them distinct types makes it impossible to accidentally call the wrong one, and names the boundary between "physical scan" and "matrix search."
- **`namespace qsyn::experimental::phasepoly`** — nests under the existing `experimental` namespace that already houses `Tableau`/`PauliRotation`, so the synthesizer composes naturally in later stages.

---

## 5. How to build & run

```sh
make unit-test                          # builds libqsyn.a + qsyn-unit-test, runs all unit tests
make -C build -j$(nproc) unit-test      # faster parallel rebuild of just the test binary
./build/qsyn-unit-test "[phasepoly]"    # run only the Stage-0 PhasePoly tests
```

Verified output:

```
Filters: [phasepoly]
All tests passed (430 assertions in 7 test cases)
```

---

## 6. What Stage 0 deliberately does *not* do (next up)

- No phase-block extraction from a real `QCir`/`Tableau` yet → **Stage 1** (`extractor`, angle normalization).
- No production Gaussian elimination / `h2` estimate → **Stage 2** (`gaussian`).
- No A\* search, search state, or successor generation → **Stage 3** (`state`, `search`).
- No `QCir` emission or qsyn command wiring → **Stage 3/4**.

Stage 0 gives every one of those a verified, shared convention to build on.
