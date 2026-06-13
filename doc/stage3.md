# PhasePoly — Stage 3 Report: Memory-Bounded A\* Single-Block Synthesis

**Goal (from `PHASEPOLY_PLAN.md`):** the core of the project — the memory-bounded A\* search that co-optimizes the joint `[P | O]` matrix (paper §3.2): apply CNOT row operations that reduce phase-parity columns (emitting an `Rz` when a column hits weight 1) while driving the output matrix `O` to identity, finishing with the Stage-2 Gaussian elimination. **Acceptance (the project gate):** the micro-cases (single `Rz`; `CX;Rz;CX`; identical-parity merge; inverse-angle cancel; multi-parity / Toffoli-style) are correct and minimal, and equivalence is verified.

---

## 1. What I implemented

In `qsyn::experimental::phasepoly`:

| File | Contents |
|---|---|
| `src/tableau/phasepoly/config.hpp` | `PhasePolyConfig` (queue cap, max solutions, expansion cap, finish/h2 modes) |
| `src/tableau/phasepoly/search.{hpp,cpp}` | `SynthesisResult`, `synthesize_phasepoly(...)` — the A\* search + greedy fallback |
| `src/tableau/phasepoly/synthesizer.{hpp,cpp}` | `build_qcir(...)`, `verify_synthesis(...)` |
| `tests/src/tableau/phasepoly/search.cpp` | 10 Catch2 tests (tag `[phasepoly]`) |

`synthesize_phasepoly(problem, config)` returns a `SynthesisResult` — the ordered, interleaved list of emitted CNOTs and `Rz`s, plus counts and search statistics.

---

## 2. The search (paper §3.2)

Each search node is a `SearchState { P, Θ, O, gates, g_cost }`. The root is the extracted problem after emitting any immediately-ready rotations.

**Successor generation.** For every *active row pair* `(i, j)` — rows that share a 1 in some phase column, i.e. `apply_cnot(i, j)` reduces that column's Hamming weight (paper §3.2.1) — copy the state, apply `apply_cnot(i, j)` to **both** `P` and `O`, append the physical `CNOT(i, j)`, then emit `Rz`s for any newly weight-1 columns and remove them.

**Cost function.** `f = g + h1 + h2` with the lexicographic tie-break `(f, h1, h2, −g)` (the paper's exact ordering; deeper states win ties), realized as a `std::map` keyed by `(f, h1, h2, SIZE_MAX−g, counter)` so `begin()` is the best node and `prev(end())` the worst.
- `g` = CNOTs emitted so far;
- `h1` = sum of Hamming weights of the remaining phase columns;
- `h2` = `linear_reversible_cnot_cost(O)` from Stage 2.

**Memory bound.** After each expansion, `while (open.size() > max_queue_size) erase the worst node` — the space-bounded A\* of §3.2.2.

**Visited cache.** A canonical key over `(P with angles, O)` maps to the best `g` seen; a successor is pruned if an equal-or-better `g` already reached that exact state. The key includes angles, so two states are merged only when genuinely identical (never wrongly).

**Multiple solutions.** Goal states (or `P`-empty states finished via Gaussian elimination) are collected up to `max_solutions`; the cheapest by CNOT count is returned (paper §3.2.2).

**Finish step.** When `P` is empty but `O ≠ I`, `synthesize_linear_reversible(O, finish_mode)` (Stage 2, PMH) appends the CNOTs that complete the basis transformation.

**Per-state vs. final cost.** Following the Stage-2 note, `h2` uses the fast `gauss_jordan` mode per state, while the one-shot finish uses block-searched `patel_markov_hayes` for quality.

### Why emitting `Rz` on the weight-1 row is correct

The key lemma (derived for this stage): writing `M_k` for the product of search row-ops and `T_k` for the circuit's running transform, `M_k = (T_k⁻¹)ᵀ`. A phase column `P_k[:,j] = M_k · p_j` equals the unit vector `e_q` **iff** `p_j = (T_kᵀ) e_q =` the parity qubit `q` currently holds. So when a column collapses to a single 1 at row `q`, qubit `q` is carrying exactly that column's parity and the rotation belongs there. Likewise `O` reaching identity means the final linear map equals the target. This is what makes the interleaved CNOT/`Rz` emission realize the right unitary — confirmed end-to-end by the equivalence tests.

---

## 3. Guarantees

- **Always returns a correct circuit.** Active row pairs exist whenever a column has weight ≥ 2, so the search always reaches `P`-empty along some path. If the search *budget* (`max_expansions`/`max_queue_size`) is exhausted first, a **greedy fallback** (`greedy_synthesize`) folds each column's set rows together and finishes `O` — guaranteed to terminate and be correct. The tiny-memory and `max_expansions = 0` tests exercise this.
- **No row swaps.** Both the search and the finish reduce `O` to the *exact* identity, so the emitted circuit contains no SWAPs.

---

## 4. Verification

Two independent checks, both used in every test via the `run()` helper:

1. **Symbolic** (`verify_synthesis`, all sizes): re-extract `(P', Θ', O')` from the synthesized gates and require the phase parity/angle **multiset** and the output matrix to match the problem. This is exact for phase-polynomial circuits — the paper's "compare symbolic phase polynomial" criterion — and needs no unitary.
2. **Unitary** (`qsyn::qcir::is_equivalent`, ≤ 7 qubits): build a `QCir` from both the original block and the synthesized gates and compare via qsyn's tableau-then-tensor equivalence checker — an independent oracle that also validates the QCir emission convention end-to-end.

---

## 5. Tests (`tests/src/tableau/phasepoly/search.cpp`)

10 cases. Every case runs symbolic verification and (for ≤ 7 qubits) full unitary equivalence:

1. **Single `Rz`** → 0 CNOTs, 1 `Rz`.
2. **`CX; Rz; CX`** → **exactly 2 CNOTs** (proved optimal) and 1 `Rz`.
3. **Redundant CNOTs** (`CX; CX; Rz`) → **0 CNOTs** (the cancelling pair is removed).
4. **Identical-parity merge** → 1 `Rz`, 2 CNOTs.
5. **Inverse rotations** → the empty circuit (0 CNOTs, 0 `Rz`).
6. **Paper two-parity polynomial** `π/2(x⊕y) + π/4(y⊕z)` → 2 `Rz`, CNOTs ≤ source.
7. **Non-identity output basis** (`q2 = x⊕y⊕z`) preserved.
8. **20 random blocks**, `n = 3..6` — symbolic + unitary equivalence.
9. **Greedy fallback** (`max_expansions = 0`) still produces a verified circuit.
10. **Tiny memory bound** (`max_queue_size = 2`) still produces verified circuits.

**Result:**
```
./build/qsyn-unit-test "[phasepoly]"
All tests passed (1791 assertions in 29 test cases)   # 7 + 8 + 4 + 10
```
Full suite, no regressions: **All tests passed (417,045 assertions in 67 test cases)**.

### Measured CNOT reduction

Average source vs. synthesized CNOTs on random blocks (30 each; measured with a temporary probe, since removed):

| n | source CX | synth CX |
|---|---|---|
| 3 | 8 | 3 |
| 4 | 8 | 5 |
| 5 | 9 | 7 |
| 6 | 8 | 7 |

The co-optimization consistently removes CNOTs; the gap is largest for small/redundant blocks and narrows as random blocks on more qubits get sparser.

---

## 6. Design choices

- **`SearchState` kept internal to `search.cpp`.** The public surface is just `PhasePolyConfig`, `SynthesisResult`, and `synthesize_phasepoly` — the search internals (priority encoding, visited cache, fallback) don't leak.
- **`std::map` open set** rather than a heap, because the memory bound needs cheap access to *both* the best node (`begin`) and the worst (`prev(end)`); the priority tuple's `counter` field guarantees unique keys.
- **`PhaseOp` reused as the output gate record.** The synthesized circuit is a `std::vector<PhaseOp>` — the same IR Stage 1 produces — so it round-trips straight back through `phase_block_to_problem` for symbolic verification, and through `build_qcir` for unitary checks.
- **Two-tier verification.** Symbolic checks give exactness at any size and run in the hot path of tests; the unitary check adds an independent oracle for small cases. Together they make a convention slip impossible to miss.

---

## 7. Known limitations / notes for later stages

- **State copying.** Each successor deep-copies the parent (`gates`, two matrices). Fine for single blocks at the tested sizes; large blocks in Stage 5 may want a parent-pointer / incremental-undo representation.
- **`h2` cost.** `gauss_jordan` per state is O(n²); acceptable now. If profiling shows it dominates on big blocks, cache or incrementally update it.
- **Active-column-set refinement.** The paper narrows the active set during a reduction "phase"; I use *all* current columns (equivalent after each completion resets) and rely on the visited cache to prevent livelock. A finer active set is a possible future pruning optimization.
- **Gate-type fidelity.** `Rz` is emitted as a `p`-gate (`PZGate`); this matches the benchmarks' decomposition and passes `is_equivalent`. If a circuit mixes `p` and `rz` with differing global phases, Stage 4 can track the source gate type per term.

---

## 8. What Stage 3 does *not* do yet

- No qsyn pass wiring — `synthesize_phasepoly` is not yet exposed as a `PauliRotationsSynthesisStrategy` or a command, and full-circuit replace-if-better is not done → **Stage 4** (qsyn integration + benchmark harness over the Table 1 circuits).
- No cross-block / SSA merging → **Stage 5**; no hardware-aware successors → **Stage 6**.

The single-block co-optimizer — the paper's central contribution — is implemented, correct, and verified.
