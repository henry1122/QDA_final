# Stage 5 — Multi-Block Warm-Start Group Synthesis

## Goal

Extend the Stage-4 single-block PhasePoly A\* to synthesize **groups of k
consecutive phase-polynomial blocks** jointly, so that row operations applied
during one block's synthesis can reduce parity columns of the next block
"for free" — without extra CNOT cost.

---

## Algorithm: warm-start group synthesis (paper §3.3)

### Setup

A quantum circuit between two Clifford layers decomposes into a sequence of
maximal phase-polynomial blocks, separated by boundary gates (primarily H gates).
`extract_phase_blocks_with_boundaries` recovers these blocks together with, for
each inter-block boundary, the set of H-gated qubits.

For a group of k consecutive blocks B₀, B₁, …, B_{k-1}:

### Accumulated parity matrix M

Track an n×n accumulated row-operation matrix **M** (starts as identity I_n).
M represents the net GF(2) linear transformation applied so far by all emitted
CNOT gates.

### Processing each block i

1. **Compute warm-started problem:**

   ```
   P_eff = M · P_i          (phase columns in current accumulated basis)
   O_eff = M · O_i          (output matrix in current accumulated basis)
   ```

   where M · X is GF(2) matrix multiplication: `result[r][c] = XOR_k M[r][k] & X[k][c]`.

2. **Synthesize** with the A\*: `synthesize_phasepoly(P_eff, O_eff, cfg)`.

3. **Update M** with all emitted CNOT row operations: for each CNOT(ctrl=c, tgt=t)
   in the synthesis result, apply `M[c] ^= M[t]`.

4. **H-boundary reset**: for each H-gated qubit j at the boundary after block i,
   reset row j of M to eⱼ (the standard basis vector).  This models the fact
   that qubit j starts fresh in the new basis after the Hadamard.

### Why the warm start is correct

After block i's synthesis, M = O_i⁻¹ (in the A\* row-operation convention;
see Stage-0 convention notes).  The warm problem (M·P_{i+1}, M·O_{i+1}) expresses
block i+1's synthesis task in the accumulated basis:

- When a column of M·P_{i+1} is already at weight 1 (= eₒ for some qubit q),
  qubit q already holds the target parity — the Rz can be emitted immediately
  without any CNOT.
- M·O_{i+1} is the remaining output transformation that must still be reduced
  to identity.

The A\* convention: `ParityMatrix::apply_cnot(c, t)` does `row[c] ^= row[t]`,
and emits physical CNOT(ctrl=c, tgt=t).  This is dual to the physical state
convention (see symbolic_state.hpp).  The relation is M^T = S, where S is the
SymbolicState parity-of-qubit matrix.

### H-boundary reset correctness

After H on qubit j, the physical qubit j enters a fresh basis.  In the M
representation, this is modelled by resetting M[j] ← eⱼ, which sets the
"accumulated parity row for qubit j" back to the standard unit vector.

---

## Concrete example (3 qubits, k=2)

Block A synthesizes parity q₀⊕q₁ with output O_A = [[1,0,0],[0,1,1],[0,0,1]].
The A\* emits CNOT(1,2) to achieve weight-1 parity, then CNOT(1,2) again to
reduce O_A to identity.  After both CNOTs, M accumulates:

```
M = I → after CNOT(1,2) → [[1,0,0],[0,1,1],[0,0,1]]
      → after CNOT(1,2) → I
```

The synthesis CNOTs cancel in M, leaving M = I.  Block B then sees its
original problem unchanged — no benefit, no regression.

For more complex output matrices (O_A ≠ I), M = O_A⁻¹ after synthesis.
Block B's warm problem becomes (O_A⁻¹·P_B, O_A⁻¹·O_B).  If the warm start
is to help, we need synth(O_A⁻¹·P_B, O_A⁻¹·O_B) < synth(P_B, O_B).

---

## Empirical results

Running the Table-1 benchmark (`--set table1`, `--no-todd`):

```
circuit                      pp(1)  pp(2)  pp(3)  pp(5)   mst+P  gray+P
tof_3_pyzx                      16     17     18     17      28      26
tof_4_pyzx                      25     27     31     31      45      43
tof_5_pyzx                      35     40     42     45      63      60
barenco_tof_3_pyzx              21     23     22     22      31      29
ham15-low_pyzx                 227    264    281    290     293     294
mod5_4_pyzx                     17     17     19     19      24      22
…
TOTAL                         3869   4468   4662   4788    5638    5569
```

**pp(k=1) is best across all circuits.**  pp(k>1) is always ≥ pp(k=1).
The warm start never improves over independent per-block synthesis.

---

## Diagnosis: why warm start doesn't help

### The zero-sum trade

After block i's complete synthesis, M = O_i⁻¹ (because A\* reduces M·O_i to
identity).  For block i+1:

```
M · O_{i+1} = O_i⁻¹ · O_{i+1}
```

Any phase CNOT "saved" (because M·P_{i+1} is easier than P_{i+1}) is exactly
offset by an output CNOT "added" (because M·O_{i+1} is harder than O_{i+1}).
The net saving requires:

```
Δ_phase  >  Δ_output
```

where Δ_phase = CNOTs saved on phase synthesis, Δ_output = extra CNOTs for
output synthesis.

### Why Δ_phase ≈ 0 for benchmark circuits

For warm-start to save phase CNOTs, M must pre-build parities that appear in
P_{i+1}.  Since M = O_i⁻¹ is determined entirely by block i's output matrix
(not its phase structure), the chance of M sharing structure with P_{i+1}
is essentially zero for typical circuits.

### Why Δ_output > 0

O_i⁻¹ · O_{i+1} is a "random" invertible matrix when O_i and O_{i+1} are
unrelated.  Reducing it to identity generally requires as many (or more) CNOTs
as reducing O_{i+1} alone.

### When warm start WOULD help

The warm start provides net savings when consecutive blocks share parity
structure — e.g., if many columns of P_{i+1} appear in the rows of O_i⁻¹
(so they're immediately at weight 1 in M·P_{i+1}) and O_i ≈ O_{i+1} (so
O_i⁻¹·O_{i+1} ≈ I, nearly trivial output).  This structure is absent in the
Toffoli, adder, and random-logic circuits in our benchmark.

---

## Root cause: greedy vs. joint optimization

The warm-start approach is **greedy and sequential**: block i's A\* makes CNOT
choices that optimize block i alone, leaving M = O_i⁻¹ for block i+1.  Block
i's A\* has no knowledge of block i+1 and therefore cannot choose CNOTs that
are jointly beneficial for both blocks.

The correct algorithm for cross-block parity sharing requires a **joint A\***
that simultaneously processes:

- `[P₀ | P₁]` — phase columns from both blocks (tagged by which block each belongs to)
- `O₀` — block-0 output matrix (must reach I before H)
- `O₁` — block-1 output matrix (must reach I at end)

In the joint A\*, a single CNOT(i,j) reduces columns from both P₀ and P₁
simultaneously — this is the true cross-block parity reuse.  P₁ columns can
only be emitted after the H transition (all P₀ done and O₀ = I).

The joint A\* has a larger state space (phase columns from both blocks) and
requires careful handling of the H transition (the H gate resets qubit j's
row in all matrices using the ORIGINAL P₁ values, not the accumulated
M·P₁ values).

---

## Files changed

| File | Change |
|---|---|
| `src/tableau/phasepoly/multiblock.hpp` | New: `BlockBoundary`, `ExtractedCircuit`, `GroupResult`, three public functions |
| `src/tableau/phasepoly/multiblock.cpp` | New: `extract_phase_blocks_with_boundaries`, `synthesize_block_group`, `synthesize_grouped` |
| `tests/src/tableau/phasepoly/multiblock.cpp` | New: 8 Catch2 test cases (26 assertions) |
| `tools/benchmark_phasepoly.cpp` | Added k=2,3,5 columns via `synthesize_grouped`; added Todd+Naive full-circuit column |
| `scripts/run_benchmarks.sh` | Added `--set table1/extended/all` flag for circuit selection |

---

## Test coverage

```
All tests passed (26 assertions in 8 test cases)
```

Tests cover: GF(2) matrix multiplication, single/multi-block H-boundary
extraction, warm-start synthesis for trivial and non-trivial output matrices,
group partitioning for k=1,2,3,5.

---

## Next step (Stage 6 / future work)

Implement the **joint A\*** for two consecutive blocks:

1. State: `(phase0, phase1, output0, output1, h_done)` where all matrices share
   the same accumulated row basis.
2. CNOT(i,j): applied to all four matrices simultaneously.
3. H transition: allowed only when `phase0` is empty and `output0 = I`.  At
   the transition, row j of `phase1` and `output1` is restored to the original
   P₁ and O₁ values (because M[j] ← eⱼ makes row j = original).
4. Goal: `phase0` empty, `h_done`, `phase1` empty, `output1 = I`.

Expected benefit: the A\* can apply a single CNOT that helps both P₀ and P₁
simultaneously — something the sequential warm-start cannot do.
