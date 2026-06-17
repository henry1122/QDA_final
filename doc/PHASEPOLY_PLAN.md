# PhasePoly Reproduction Plan (qsyn / C++)

**Paper:** *PhasePoly: An Optimization Framework for Phase Polynomials in Quantum Circuits* (Chen et al., arXiv:2506.20624, 2025), in `doc /phase-poly.pdf`.

**Decision (confirmed):** Implement natively in **C++ inside qsyn**, reusing existing infrastructure. **Single-block co-optimization is the first priority** — it must beat MST/GraySynth on the Table 1 logical benchmarks before moving on.

---

## 1. What the paper actually does (condensed)

A phase-polynomial circuit (CNOT + `Rz` only) realizes
`U|x⟩ = e^{i·p(x)} |g(x)⟩`, where `p(x)` is a sum of `θ_j · (parity_j of x)` and `g(x)` is an affine/linear reversible (CNOT) basis map.

Three contributions, in order of the paper's emphasis:

1. **Single-block co-optimization (§3.1, the core).** Instead of synthesizing the phase-parity network and the output basis transform separately (prior work: Gray-Synth, Vandaele/MST), maintain a **joint matrix `[P | O]`**:
   - `P` (n×m): each column = one `Rz` parity term (over input vars), each row = a qubit. Angles `Θ` (length m).
   - `O` (n×n): each column = the symbolic output currently carried by a qubit; starts as the target basis transform, goal = identity.
   - A single CNOT applies the **same row operation to both** `P` and `O`. Paper's convention (§3.1, Eq. 5): for `CNOT` with control row `i`, target row `j`:
     `row_i ← row_i ⊕ row_j` (control row changes; target row unchanged).
     > ⚠ This is the *inverse-transpose* of the usual "target changes" model, because rows store *parity-of-inputs*. It must be unit-tested against a forward symbolic simulator and against qsyn's `PauliRotation::cx` (see Stage 0).
   - When a `P` column reaches Hamming weight 1 at row `q`: emit `Rz(θ)` on qubit `q`, drop that column (§3.1 rule 2).
   - Goal: `P` empty **and** `O = I`.

2. **Memory-bounded A\* search (§3.2).** Tree where nodes = circuit states, edges = CNOTs.
   - **Active row pairs**: a pair `(i,j)` is active iff `row_i ← row_i ⊕ row_j` reduces the Hamming weight of ≥1 column in the **active column set** (prevents livelock). Active set re-initializes to all remaining columns after any column completes.
   - Cost `f(n) = g(n) + h1(n) + h2(n)`:
     - `g` = CNOTs emitted so far.
     - `h1` = Σ Hamming weights of remaining `P` columns (phase-parity estimate).
     - `h2` = CNOTs to reduce `O` to `I` (run Gaussian elimination on a copy, count row ops).
   - Tie-break tuple `[f, h1, h2, −g]` (lexicographic; deeper paths win ties).
   - **Space-bounded**: cap priority queue at `max_queue_size`; drop worst nodes when exceeded.
   - **Multiple-solution search**: collect up to `k` goal states, return the cheapest by actual cost.
   - When `P` becomes empty mid-search, **finish `O` with (minimum-degree) Gaussian elimination** and record a complete solution.

3. **Multi-block optimization + SSA-style IR (§3.3).** General circuits interleave non-`Rz`/non-CNOT gates (e.g. H), splitting the circuit into phase-poly blocks. Rename qubits SSA-style at each non-`Rz`/non-X boundary (`q3 → q3' → q3''`). Renaming lets unchanged parities (e.g. `q0⊕q2`) be **reused across merged blocks**, cutting CNOTs. A merged block produces one larger parity matrix; dependency rules force earlier qubit-versions to be synthesized before later ones; only one matrix row is "active" per qubit-line at a time. **Group size `k`** controls how many adjacent blocks merge (try 1,2,3,5; larger isn't always better — Table 3).

4. **Hardware-aware extension (§3.4, utility).** State also tracks logical↔physical mapping. Successors: (A) executable active-row-pair CNOTs satisfying connectivity; (B) CNOTs that shrink inter-qubit distance in active columns even without weight reduction; (C) SWAPs only when no (A) exists and they reduce a column's min pairwise distance. Distance = shortest path on coupling graph, edge weight 1 into an active qubit, 3 through an inactive one (SWAP = 3 CNOTs). Heuristics become `h1'` (Eq. 8: Σ_columns Σ_{active pairs} Dist/(|active|−1)) and `h2'` (Gaussian-elim CNOTs costed by mapped distance). Sabre-style iterative initial-mapping refinement (forward/reverse passes).

**Targets to reproduce:** Table 1 (single-block vs MST/Rotation-Merging: avg 34.92% gate / 28.54% CNOT reduction), Table 3 (group-size sensitivity), Fig. 11/12 (physical).

---

## 2. How this maps onto qsyn (what already exists)

qsyn (C++20) already contains ~80% of the needed scaffolding. Key reusable pieces found in the tree:

| Paper concept | qsyn artifact | Location |
|---|---|---|
| `Rz`-on-parity term | `PauliRotation` (diagonal Z-product + `dvlab::Phase`); `is_z(i)` gives parity support, `phase()` gives `θ` | `src/tableau/pauli_rotation.hpp` |
| Phase polynomial block | `std::vector<PauliRotation>`; `is_phase_polynomial(...)` checks all-diagonal | `src/tableau/tableau_optimization.hpp` |
| Multi-block IR (Clifford ⟷ rotation partition) | `Tableau` = `vector<SubTableau>`, `SubTableau = variant<StabilizerTableau, vector<PauliRotation>>` | `src/tableau/tableau.hpp` |
| Circuit → blocks (extraction) | `to_tableau(QCir)` | `src/convert/qcir_to_tableau.{hpp,cpp}` |
| Output basis transform `O` (linear/CNOT part) | `StabilizerTableau` (stabilizer/destabilizer rows = Z/X images); `cx`, `prepend_cx`, `swap` | `src/tableau/stabilizer_tableau.hpp` |
| Blocks → circuit (synthesis) | `to_qcir(Tableau, ...)`, `to_qcir(vector<PauliRotation>, PauliRotationsSynthesisStrategy)` | `src/convert/tableau_to_qcir.{hpp,cpp}` |
| **Baselines to beat** | `MstSynthesisStrategy` (Vandaele), `GraySynthPauliRotationsSynthesisStrategy`, `NaivePauliRotationsSynthesisStrategy`, `TParPauliRotationsSynthesisStrategy` | `src/convert/tableau_to_qcir.hpp` |
| Existing phase-poly opt hook | `PhasePolynomialOptimizationStrategy` (`ToddPhasePolynomialOptimizationStrategy`), `optimize_phase_polynomial(...)`; command `tableau ... phasepoly <strategy>` | `src/tableau/tableau_optimization.{hpp,cpp}`, `src/cmd/tableau_cmd.cpp` |
| Coupling graph (HW-aware) | `Device` | `src/device/device.{hpp,cpp}` |
| Verification | `qcir_equiv` (`src/qcir/qcir_equiv.cpp`), `qcir_to_tensor` for unitary compare | `src/qcir/`, `src/convert/` |
| Benchmarks | `tof_*`, `barenco_tof_*`, `adder_*`, `grover_5`, `ham15-*`, `qcla_*`, `mod_*`, `vbe_adder_3`, `hwb6` (paper's Table 1 set) | `benchmark/qc/`, `benchmark/qasm/` |

**Integration verdict:** The new synthesizer slots in as **a new `PauliRotationsSynthesisStrategy` subclass** (`PhasePolySynthesisStrategy`), with the *output-matrix co-optimization* fed from the adjacent `StabilizerTableau`'s **linear part**. This is the single most important design point that differentiates us from the existing MST strategy, which only synthesizes the rotations and leaves the basis transform to a separate engine.

Because `to_qcir(vector<PauliRotation>, strategy)` doesn't receive the output basis, we add a **co-optimizing entry point** that takes `(StabilizerTableau output_basis, vector<PauliRotation> rotations)` → `QCir`, and call it from `to_qcir(Tableau, ...)` so the joint `[P|O]` optimization happens block-by-block during tableau→circuit synthesis.

---

## 3. Proposed C++ module layout

New files under `src/tableau/phasepoly/` (mirrors the user's logical decomposition, but C++ and reusing qsyn types):

```
src/tableau/phasepoly/
    parity_matrix.hpp/.cpp     # P, O over GF(2) (dvlab dynamic_bitset rows); row_xor, hamming, ready-column extraction
    state.hpp/.cpp             # SearchState: P, Θ, O, gates, rz_ops, g_cost, active_columns, (HW) mapping; hashing
    search.hpp/.cpp            # memory-bounded A*: priority queue, successors, active row pairs, multi-solution
    gaussian.hpp/.cpp          # CNOT-only linear-reversible synthesis (minimum-degree Gaussian elim) + h2 estimate
    extractor.hpp/.cpp         # PauliRotation block + StabilizerTableau -> PhasePolyProblem {P, Θ, O}
    multiblock.hpp/.cpp        # SSA renaming IR, block grouping (k), dependency ordering   [Stage 5]
    hardware.hpp/.cpp          # coupling-graph distance, HW successors (A/B/C), h1'/h2', Sabre init-map  [Stage 6]
    config.hpp                 # PhasePolyConfig (queue size, k, solutions, timeout, cost metric, verify, ...)
    synthesizer.hpp/.cpp       # PhasePolySynthesisStrategy + co-optimizing to_qcir entry point; PhasePolyProblem
```

Wire-in points (edit existing files):
- `src/convert/tableau_to_qcir.{hpp,cpp}`: register `PhasePolySynthesisStrategy`; add co-optimizing overload that consumes the adjacent `StabilizerTableau`.
- `src/cmd/tableau_cmd.cpp` and/or `src/cmd/conversion_cmd.cpp`: expose `phasepoly` as a synthesis strategy string + config flags (queue size, k, timeout, hardware-aware, coupling graph).
- `CMakeLists.txt`: add the new sources.
- `tests/tableau/` and `tests/optimization/`: dofiles + reference outputs.

---

## 4. Staged plan

Each stage has a concrete **deliverable** and **acceptance criteria**. Stages 0–4 are the priority (single-block). 5 and 6 are follow-ons.

### Stage 0 — Convention lock-in & GF(2) primitives  *(foundation, ~0.5–1 day)*
- Implement `parity_matrix` (rows as bitsets), `row_xor(M, i, j)`, `hamming_weight(col)`, `ready_columns`, equality/hash.
- Build a tiny **forward symbolic simulator** (track each qubit's input-parity through CNOTs) used only in tests.
- **Unit-test the row-op convention**: assert paper's `row_i ← row_i ⊕ row_j` (control=i) reproduces the symbolic-simulator state *and* matches `PauliRotation::cx` conjugation and physical `emit_cnot`. Lock the control/target ↔ row mapping here; everything downstream depends on it.
- **Accept:** convention tests pass; Eq. (4) example from the paper (`CNOT(2,1)` on the given 3×5 matrix) reproduced exactly.

### Stage 1 — Extraction: block → `PhasePolyProblem {P, Θ, O}`  *(~1 day)*
- `extractor`: from a `vector<PauliRotation>` block (the `Rz`s) build `P` (columns = Z-supports) and `Θ` (phases). From the adjacent linear part of the `StabilizerTableau` build `O`.
- **Angle normalization (paper §"Important details"):** merge identical parity columns (`θ ← θ_a + θ_b`), drop columns with `θ ≡ 0 (mod 2π)`. Reuse `dvlab::Phase` arithmetic + existing `merge_rotations`/`remove_identities`.
- **Accept:** round-trip a hand-built block; `P,Θ,O` match a by-hand derivation for the Fig. 1 / Fig. 6 examples.

### Stage 2 — Gaussian-elimination finish + `h2`  *(~1 day)*
- `gaussian`: CNOT-only linear-reversible synthesis reducing `O → I`, emitting CNOTs (start greedy/minimum-degree; the paper says "minimum degree-based"). Same routine, in counting mode, gives `h2(n)`.
- **Accept:** for random invertible GF(2) `O`, synthesized CNOTs reproduce `O` (verified by replaying row ops); cost is reasonable vs. naive.

### Stage 3 — Memory-bounded A\* single-block search  *(core, ~2–3 days)*
- `state` + `search`: priority queue keyed `[f, h1, h2, −g, counter]`; active-row-pair successor generation; eager ready-column `Rz` emission; visited cache keyed on hash(`P`,`O`,active set) storing best `g`; `max_queue_size` pruning; multi-solution (`k`); timeout/expansion caps; fallback to MST/Gaussian if no goal found.
- `synthesizer`: assemble emitted CNOTs + `Rz`s + Gaussian finish into a `QCir`; `PhasePolySynthesisStrategy` implementing `PauliRotationsSynthesisStrategy`; co-optimizing `to_qcir(StabilizerTableau, vector<PauliRotation>)`.
- **Accept (the gate of the whole project):** on the testing-plan micro-cases (single `Rz`; `CX;Rz;CX`; identical-parity merge; inverse-angle cancel; Toffoli decomposition) results are correct and minimal. `verify_equivalence` passes via `qcir_equiv`/tensor.

### Stage 4 — qsyn wiring, verification, single-block benchmark harness  *(~1–2 days)*
- Register `phasepoly` strategy in the tableau/conversion commands with config flags. `replace_only_if_better`: keep original block if synthesized cost ≥ original.
- `verifier`: small blocks → exact/numeric unitary compare; large → symbolic (same parity-angle multiset + same affine output). Gate behind `config.verify`.
- **Benchmark script** (dofile + Python report) over Table 1 circuits printing: original/optimized total gates, original/optimized CNOTs, runtime, max queue size, #expanded states.
- **Accept:** on `tof_3, barenco_tof_3, tof_4/5, vbe_adder_3, barenco_tof_*, grover_5, ham15-*`, PhasePoly **matches or beats MST** in CNOT and total gates, trending toward the paper's ~28%/~35% averages; all verified.

### Stage 5 — Multi-block SSA IR + group merging  *(~3–4 days)*
- `multiblock`: SSA renaming at non-`Rz`/non-`X` boundaries; merge `k` adjacent phase-poly blocks into one expanded parity matrix; dependency ordering (earlier qubit-versions first); one active row per qubit-line at a time.
- Progressive strategy: optimize per-block (k=1) → try k=2 with best-so-far as bound → k=3,5; keep best.
- **Accept:** reproduce Table 3 trends — `adder_8` improves with larger `k`, `barenco_10` jumps at k=3, `ham15_low` best at k=1. Cross-block parity reuse (e.g. `q0⊕q2`) reduces CNOTs vs. independent blocks.

### Stage 6 — Hardware-aware synthesis (optional)  *(~3–4 days)*
- `hardware`: extend state with logical↔physical maps; successor cases (A)/(B)/(C); weighted coupling-graph distance (1 active / 3 inactive); `h1'` (Eq. 8) and `h2'`; SWAP costed as 3 CNOTs. Sabre-style iterative initial-mapping refinement.
- Integrate with qsyn `Device`; optionally feed into the existing Duostra mapper for comparison.
- **Accept:** SWAP-free synthesis on the Fig. 9 example; linear-connectivity Toffoli template generated; CNOT-after-mapping beats `Original+mapping` on a representative subset (Fig. 11/12 directionally).

---

## 5. Key risks & how the plan addresses them

1. **Row-op convention errors** (most likely source of silent bugs). → Stage 0 dedicates a test harness reconciling paper rows ↔ `PauliRotation::cx` ↔ physical CNOT before any search code exists.
2. **Output-matrix `O` orientation** (stabilizer vs destabilizer rows; column = qubit-output vs input). → Stage 1 pins it by reproducing the paper's Eq. (4) and Fig. 6 examples exactly, and Stage 3's `verify_equivalence` catches regressions.
3. **Search blow-up / timeouts.** → memory-bounded queue + visited cache + multi-solution cap + timeout, all configurable; fallback to MST so we never produce worse-than-baseline output.
4. **Heuristic admissibility vs. speed.** Paper's `h1`/`h2` are not strictly admissible (it's bounded best-first, not optimal A\*). → match the paper exactly first; treat `h1 = Σ max(0, w−1)` only as a documented experimental variant.
5. **Verification cost on large circuits.** → tiered verifier (unitary for small, symbolic for large), matching the paper's Qiskit/MQT-QCEC approach.

## 6. Open items to confirm as we build
- Exact "minimum-degree Gaussian elimination" variant (paper cites linear-reversible synthesis refs [10,39,45]); start greedy, refine if CNOT counts lag Table 1.
- Whether to expose PhasePoly also as a `PhasePolynomialOptimizationStrategy` (tableau→tableau rewrite) in addition to the synthesis strategy — decide after Stage 4 based on how cleanly the co-opt entry point composes.

---

## 7. Immediate next actions (on approval)
1. Stage 0: create `src/tableau/phasepoly/parity_matrix.{hpp,cpp}` + convention unit tests; add to CMake.
2. Stage 1: `extractor` + angle normalization, validated on Fig. 1/Fig. 6.
3. Proceed Stage 2 → 3 → 4 to clear the single-block acceptance gate before touching multi-block.
