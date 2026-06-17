# PhasePoly — Stage 2 Report: Gaussian-Elimination Finish + `h2`

**Goal (from `PHASEPOLY_PLAN.md`):** CNOT-only linear-reversible synthesis that reduces the output matrix `O` to identity, emitting the CNOTs (the "Gaussian finish" used when the phase matrix `P` is empty, paper §3.2). The same routine, in counting mode, yields the heuristic `h2(n)` — the estimated output-parity cost. **Acceptance:** for random invertible GF(2) `O`, the synthesized CNOTs reproduce `O` (verified by replaying the ops); cost is reasonable vs. naive.

---

## 1. What I implemented

In `qsyn::experimental::phasepoly`:

| File | Contents |
|---|---|
| `src/tableau/phasepoly/gaussian.{hpp,cpp}` | `CnotOp`, `LinearSynthesisMode`, `synthesize_linear_reversible(...)`, `linear_reversible_cnot_cost(...)` |
| `tests/src/tableau/phasepoly/gaussian.cpp` | 4 Catch2 tests (tag `[phasepoly]`) |

### 1.1 API

```cpp
struct CnotOp { size_t control; size_t target; };

enum class LinearSynthesisMode { gauss_jordan, patel_markov_hayes };

std::vector<CnotOp> synthesize_linear_reversible(
    ParityMatrix const& matrix,
    LinearSynthesisMode mode = LinearSynthesisMode::patel_markov_hayes);

size_t linear_reversible_cnot_cost(
    ParityMatrix const& matrix,
    LinearSynthesisMode mode = LinearSynthesisMode::patel_markov_hayes);
```

`synthesize_linear_reversible` returns the CNOT sequence that reduces a square, invertible `O` to identity. Per the locked convention, **emitting each `CnotOp` as the physical `CNOT(control, target)` in order reconstructs the linear map `O` encodes.** `linear_reversible_cnot_cost` is the same computation returning only the count — this is `h2(n)`.

### 1.2 Two strategies

- **`gauss_jordan` (convention-native reference).** Column-by-column elimination on a `ParityMatrix`, keeping each finished column a unit vector *on the diagonal* (no row swaps → result is exactly the identity, never a permuted one). When a diagonal pivot is missing it is seeded from the **minimum-weight eligible row below it** (a min-degree-style fill-reducing choice), then all other 1s in the column are cleared. Records `apply_cnot(c, t)` directly as `CnotOp{control = c, target = t}`. Simple, self-contained, trivially correct by construction.

- **`patel_markov_hayes` (default, fewer gates).** Reuses qsyn's tested block elimination `dvlab::BooleanMatrix::gaussian_elimination_skip(block_size, /*fully_reduced=*/true, /*track=*/true)` (the PMH algorithm, the same one the ZX extractor uses). It **searches block sizes `1..n` and keeps the sequence with the fewest CNOTs.**

---

## 2. The convention bridge (the subtle part)

`BooleanMatrix::row_operation(ctrl, targ)` performs `row[targ] ^= row[ctrl]` — the **mirror** of our `apply_cnot(c, t)` (`row[c] ^= row[t]`), exactly as documented in Stage 0. Therefore a recorded PMH op `(ctrl, targ)` is the native operation `apply_cnot(control = targ, target = ctrl)`, and is emitted as the physical **`CNOT(control = targ, target = ctrl)`**, preserving order.

I re-derived this from the Stage-0 lemma `O = Tᵀ`: if a sequence of row ops reduces `O` to `I`, emitting each as a physical CNOT with the *native* `(control, target)` labels, in the same order, rebuilds `O`. The translation is then a pure index swap on the recorded PMH pairs. Both strategies are checked against this in two independent ways (matrix-space reduction **and** physical replay), so any orientation error fails loudly.

---

## 3. Tests (`tests/src/tableau/phasepoly/gaussian.cpp`)

Each test runs **both** strategies and applies two complementary checks:
- *reduce check* — applying the ops as matrix row operations to a copy of `O` yields the identity;
- *replay check* — emitting the ops as physical CNOTs from the identity and reading back the map equals `O`.

1. **Identity** — emits zero CNOTs; cost is 0.
2. **Paper output basis** — `O = |x, x⊕z, x⊕y⊕z⟩` (the right block of Eq. 4) reduces to `I` and replays back to `O`.
3. **Random invertible matrices** — 150 trials, `n = 2..7`, both modes: reduce + replay both hold, `cost == ops.size()`, and `ops.size() ≤ n²` (comfortably under the naive bound).
4. **Stage-1 integration** — synthesizes the `O` produced by the Stage-1 extractor for the 4-CNOT Eq. (4) network and verifies reduce + replay.

**Result:**
```
./build/qsyn-unit-test "[phasepoly]"
All tests passed (1678 assertions in 19 test cases)   # 7 + 8 + 4
```
Full suite, no regressions: **All tests passed (775,959 assertions in 57 test cases)**.

---

## 4. "Reasonable vs. naive" — measured

Average CNOTs to synthesize random invertible `n×n` maps (50 matrices each), Gauss-Jordan vs. PMH:

| n | gauss_jordan | patel_markov_hayes | reduction |
|---|---|---|---|
| 8  | 21  | 18  | ~14% |
| 12 | 39  | 34  | ~13% |
| 16 | 63  | 54  | ~14% |
| 24 | 142 | 116 | ~18% |

PMH is consistently cheaper and the gap widens with size (the expected `O(n²/log n)` benefit), which is why it is the default. (Measured with a temporary instrumented test that was removed afterward to keep the suite clean.)

---

## 5. Design choices

- **Two strategies, cross-validated.** The native Gauss-Jordan is the trivially-correct reference and a zero-dependency fallback; PMH is the quality default. Running both through identical reduce+replay checks means each guards the other against convention bugs.
- **Reuse qsyn's PMH rather than re-rolling it.** It is already tested and used in the extractor; via the proven argument-swap we get its low gate counts for free. This is exactly the reuse the Stage-0 report anticipated.
- **Emit in matrix-reduction order.** The Stage-0 `O = Tᵀ` round trip established that the reduction order *is* the emission order — no reversal — so neither strategy needs to flip the sequence.
- **No row swaps.** Both strategies reduce to the *exact* identity without permuting rows, so the emitted circuit needs no SWAPs (which would cost 3 CNOTs each) — important for the logical-synthesis quality target.

---

## 6. Notes for Stage 3

- **`h2` performance.** `linear_reversible_cnot_cost` currently calls `synthesize_linear_reversible` and, for PMH, searches all block sizes. That is fine for the finish step but is called per-state by the A* heuristic. Stage 3 should evaluate `h2` with a single fixed block size (or the `gauss_jordan` mode) for speed, reserving the block-size search for the final emission of the chosen solution.
- **Finish step.** When the search reaches a state with `P` empty but `O ≠ I`, it calls `synthesize_linear_reversible(state.O)` and appends those CNOTs to complete a candidate solution (paper §3.2 / plan Stage 3).
- `CnotOp` is the shared gate record the search will accumulate alongside `Rz` emissions.

---

## 7. What Stage 2 does *not* do

- No search yet — this is only the finish/cost building block → **Stage 3** (`state`, A* `search`, successor generation, ready-column `Rz` emission, `QCir` assembly).
- No hardware/distance-aware costing (`h2'`) → **Stage 6**.
