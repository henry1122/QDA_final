# PhasePoly — Stage 1 Report: Block Extraction → `(P, Θ, O)`

**Goal (from `PHASEPOLY_PLAN.md`):** Turn a phase-polynomial block into the three matrices the co-optimizer searches over — the phase-parity matrix `P`, the angle list `Θ`, and the output-basis matrix `O` — with angle normalization (merge identical parities, cancel zero-angle terms). **Acceptance:** round-trip hand-built blocks and match a by-hand derivation of the paper's Fig. 1 / Fig. 6 (Eq. 4) examples.

Builds directly on the Stage 0 primitives (`ParityMatrix`, `SymbolicState`) and their locked convention.

---

## 1. What I implemented

All in `qsyn::experimental::phasepoly` (auto-globbed into the build).

| File | Contents |
|---|---|
| `src/tableau/phasepoly/phase_block.hpp` | `PhaseOp` + `PhaseBlock` — the phase-poly block IR (a sequence of CX / Rz over `n` qubits) |
| `src/tableau/phasepoly/phase_poly_problem.hpp` | `PhasePolyProblem` — the `{n_qubits, P, Θ, O}` bundle fed to the synthesizer |
| `src/tableau/phasepoly/extractor.{hpp,cpp}` | `phase_block_to_problem(...)` and `extract_phase_blocks(QCir)` |
| `tests/src/tableau/phasepoly/extraction.cpp` | 8 Catch2 tests (tag `[phasepoly]`) |

### 1.1 `PhaseBlock` (the block IR)

A maximal run of phase-polynomial gates over `n` qubits. `PhaseOp` is a small tagged record: a `cx` carries `(control, target)`; an `rz` carries `(qubit, dvlab::Phase)`. The block offers `append_cx`, `append_rz`, `ops()`, and `num_cx`/`num_rz` counters. This is the clean boundary between "circuit scanning" and "matrix building" — tests construct blocks directly without needing a `QCir`.

### 1.2 `PhasePolyProblem`

Bundles `P` (`n × m`), `Θ` (length `m`), and `O` (`n × n`), kept in lockstep so column `j` of `P` carries angle `Θ[j]`. Helpers: `num_phase_terms()`, `is_trivial()` (no phase terms and `O = I`), and `to_string()` for debugging.

### 1.3 `phase_block_to_problem` (extraction core)

Implements the paper's single-block extraction using the Stage-0 **forward** simulator (`SymbolicState`, physical semantics `value[target] ^= value[control]`):

1. Scan ops left to right, applying each CX to the running parities.
2. On each `Rz(θ)` on qubit `q`, take the parity `value[q]` as a phase column with angle `θ`.
3. **Angle normalization** (paper §"Important details"):
   - merge identical parity bit-vectors, summing angles (keyed by the parity bit-string, first-appearance order preserved);
   - drop any column whose accumulated angle is `0 (mod 2π)` — handled exactly by `dvlab::Phase`, which is normalized mod 2π (so `π/4 + (−π/4) == Phase()`).
4. `O = SymbolicState::to_output_matrix()` at block end (column `q` = the final parity carried by qubit `q`).

### 1.4 `extract_phase_blocks(QCir)`

Scans the circuit in topological order (`qcir.get_gates()`), classifying each gate:
- **CX**: a `ControlGate` with exactly one control whose target is `X` (`PXGate(π)`); pins are `(control, target)`.
- **Rz**: a single-qubit Z-rotation — `PZGate` (`p`/`z`/`s`/`t`/…) or `RZGate` (`rz`).
- **anything else** (e.g. `H`, `CZ`, `CCX`): a **boundary** that closes the current block.

Blocks span all `n` qubits, so qubit ids map directly to matrix rows (consistent with `to_tableau`). Multi-controlled gates and `CZ` are treated as boundaries here — benchmarks are pre-decomposed to CX+Rz+H, and maximal cross-boundary merging is the Stage 5 (SSA multi-block) job.

---

## 2. Faithfulness to the paper (Eq. 4)

Two tests reproduce the two halves of the paper's joint `[P | O]` matrix from Fig. 6 / Eq. (4), derived independently by forward simulation:

- **Phase half** — block `CX(0,1); Rz(π/2); CX(0,1); CX(2,1); Rz(π/4); CX(2,1)` yields exactly
  `p(x,y,z) = π/2·(x⊕y) + π/4·(y⊕z)`, i.e. `P` columns `(1,1,0)` then `(0,1,1)` with `Θ = [π/2, π/4]`, and `O = I`.
- **Output half** — the 4-CNOT network `CX(0,1); CX(1,2); CX(2,1); CX(0,1)` produces the output basis `|x, x⊕z, x⊕y⊕z⟩`, i.e. `O` equal bit-for-bit to the right block of Eq. (4):
  ```
  1 1 1
  0 0 1
  0 1 1
  ```

---

## 3. Tests (`tests/src/tableau/phasepoly/extraction.cpp`)

8 cases, covering the plan's testing-plan micro-cases:

1. **Single `Rz`** — `P = [1]`, `Θ = [θ]`, `O = I₁`.
2. **`CX; Rz; CX`** — exposes parity `q0⊕q1` (`P` column `(1,1)`) and restores `O = I₂`.
3. **Paper phase polynomial** — the `π/2·(x⊕y) + π/4·(y⊕z)` reproduction above.
4. **Paper output basis** — the Eq. (4) `O` reproduction above (no phase terms).
5. **Identical-parity merge** — two `π/4` rotations on the same parity merge into one `π/2` column.
6. **Inverse-angle cancel** — `π/4` then `−π/4` drops the column entirely (`P` empty, `O = I`).
7. **`extract_phase_blocks` boundary split** — `H` between two CX/Rz runs yields 2 blocks.
8. **`extract_phase_blocks` gate decoding** — CX control/target and Rz qubit are read correctly, and the block round-trips through `phase_block_to_problem` to the expected `P, Θ`.

**Result:**
```
./build/qsyn-unit-test "[phasepoly]"
All tests passed (465 assertions in 15 test cases)   # 7 Stage-0 + 8 Stage-1
```
Full suite, no regressions: **All tests passed (344,904 assertions in 53 test cases)**.

---

## 4. Build-system note

Stage-1 unit tests construct a real `QCir`, whose gate headers (`operation.hpp`) transitively pull in `xtensor-blas`/zx includes. The `unit-test` target previously linked only `fmt`/`spdlog`/`GSL`/`sul`, so it couldn't compile qcir headers. I mirrored `libqsyn.a`'s full system-include + lapack/blas link set onto the `unit-test` target in `CMakeLists.txt`. This is a one-time, test-target-local change that also unblocks Stages 3–4, which will unit-test `QCir` emission and `qcir_equiv` verification.

> Reminder for future stages: after adding a new `src/**` or `tests/src/**` file, re-run `cmake` configure once (`make configure` or the `cmake -S . -B build …` line) so `GLOB_RECURSE` picks it up before `make -C build unit-test`.

---

## 5. Design choices

- **`PhaseBlock` as an explicit CX/Rz IR**, decoupled from `QCir`. The core algorithm (extraction, and later the search) is tested on hand-built blocks with zero circuit-framework dependency; `extract_phase_blocks` is the thin, separately-tested adapter from real circuits.
- **Forward (physical) extraction, dual to the synthesis convention.** Extraction uses `SymbolicState` (`value[target] ^= value[control]`); synthesis will use `ParityMatrix::apply_cnot` (`row[control] ^= row[target]`). Keeping them distinct types prevents mixing the two XOR directions — the exact trap Stage 0 pinned down.
- **Angle bookkeeping via `dvlab::Phase`.** Its mod-2π normalization gives correct merge/cancel for free, and reusing it keeps angles bit-compatible with the rest of qsyn (`PauliRotation`, gate phases) for Stage 3 emission.
- **Heavy gates are boundaries, not errors.** `CZ`/`CCX`/`H` simply end a block, so extraction is total over any circuit; exploiting them is deferred to multi-block (Stage 5).

---

## 6. What Stage 1 does *not* do yet

- No synthesis: `(P, Θ, O)` is produced but not yet reduced to a circuit → **Stage 2** (Gaussian-elim finish + `h2`) and **Stage 3** (A\* search, `state`, successors, `QCir` emission).
- No gate-type fidelity decision (PZ vs RZ global phase) on emission — recorded as a `dvlab::Phase`; the emit-side choice is a Stage 3/4 concern and unobservable to the unitary-equivalence checks used there.
- No cross-boundary merging / SSA renaming → **Stage 5**.

Stage 1 hands Stage 2/3 a normalized, verified `PhasePolyProblem` for every block of any input circuit.
