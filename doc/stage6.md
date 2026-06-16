# Stage 6 — Joint A* for Two Consecutive Blocks

## Goal

Replace the Stage-5 greedy warm-start with a **joint A\*** that processes two
consecutive phase-polynomial blocks simultaneously.  In the joint A\*, a single
CNOT gate can reduce parity columns from **both** blocks at once, achieving true
cross-block CNOT sharing — something the sequential warm-start cannot do.

---

## Why Stage 5 failed (recap)

After block i's synthesis, the accumulated row-operation matrix M equals
O_i⁻¹.  Block i+1 then sees the warm-started problem (O_i⁻¹·P_{i+1},
O_i⁻¹·O_{i+1}).  Since O_i⁻¹·O_{i+1} is a "random" harder matrix for
typical circuits, the output restoration cost increases more than the phase cost
decreases — a zero-sum trade that always regresses:

```
Stage-5 warm-start:  pp(k=1) = 3869,  pp(k=2) = 4468  (+15% WORSE)
```

---

## Algorithm: joint A* (paper §3.3 — joint formulation)

### State

```
(phase0, ang0, phase1, ang1, output0, output1, phase1_orig, output1_orig, h_done)
```

where:
- `phase0`, `output0`: block 0 phase columns and output matrix
- `phase1`, `output1`: block 1 phase columns and output matrix
- `phase1_orig`, `output1_orig`: originals of block 1 (needed for H-boundary row reset)
- `h_done`: whether the H-boundary transition has been applied

### H-boundary feasibility pre-check

At H transition, the accumulated row-operation matrix M equals O0_orig⁻¹.
Resetting H-gated qubit j's row to eⱼ gives M_new, which must be invertible
for block 1's synthesis to be valid.

**Key identity**: M_new · O0_orig = Q, where Q is the matrix obtained by
replacing H-qubit rows of the identity with the corresponding rows of O0_orig.
Since O0_orig is always invertible: `M_new invertible ↔ Q invertible`.

Q is cheap to build (no matrix inversion needed):

```
Q[j] = O0_orig[j,:]  for j in h_qubits
Q[k] = eₖ           for k not in h_qubits
```

If `rank(Q) < n`, the block pair is **infeasible** for joint A\* and falls back
to independent per-block synthesis.

### Operations

**Before h\_done** — CNOT(i,j) acts on all four matrices simultaneously:
```
phase0.apply_cnot(i, j)    output0.apply_cnot(i, j)
phase1.apply_cnot(i, j)    output1.apply_cnot(i, j)
```
This is the **joint benefit**: a CNOT chosen to reduce a phase0 column also
reduces the corresponding phase1 column, halving the cost for shared parities.

**H transition** (when phase0 empty → deterministic finish\_output0, then reset):
1. `finish_output0`: apply Gaussian-elimination CNOTs to drive output0 → I;
   the same CNOTs ALSO transform output1 and phase1 (further cross-block benefit).
2. `apply_h_transition`: for each H-gated qubit j:
   - `phase1[j,:] ← phase1_orig[j,:]`
   - `output1[j,:] ← output1_orig[j,:]`
3. Set `h_done = true`; eagerly emit weight-1 phase1 columns as Rz.

**After h\_done** — CNOT(i,j) acts on phase1, output1 only.

### Heuristic

```
h1 = phase_cost(phase0) + phase_cost(phase1)   [both blocks' remaining work]
h2 = linear_reversible_cnot_cost(output0)       [before h_done]
   | linear_reversible_cnot_cost(output1)       [after h_done]
f  = g_cost + h1 + h2
```

Including phase1 in h1 (even before h\_done) guides the search toward CNOTs
that simultaneously benefit both blocks.

### Fallback guarantee

After the joint A\* completes, the result is compared against **independent
per-block synthesis**.  If the joint solution is no better, the independent
counts are returned.  This guarantees:

```
pp(k > 1, joint) ≤ pp(k=1, independent)   for every circuit
```

### Group partitioning (how k affects pairing)

`synthesize_grouped(k)` partitions the global block list into non-overlapping
groups of size k.  Within each group, **consecutive pairs** are synthesized
jointly; a trailing singleton is synthesized independently:

| k | Group of 5: [B0,B1,B2,B3,B4] |
|---|---|
| 1 | B0, B1, B2, B3, B4 (all independent) |
| 2 | joint(B0,B1), joint(B2,B3), B4 |
| 3 | joint(B0,B1), B2 │ joint(B3,B4), (no trailing) |
| 5 | joint(B0,B1), joint(B2,B3), B4 |

k=2 maximises pairing density (every pair of consecutive blocks is jointly
synthesized) and gives the best results.

---

## Empirical results — Table 1 benchmark (--set table1, --no-todd)

```
circuit                                pp(1)  pp(2)  pp(3)  pp(5)   mst+P  gray+P
tof_3_pyzx                               16     16     16     16      28      26
tof_4_pyzx                               25     25     25     25      45      43
tof_5_pyzx                               35     35     34     35      63      60
tof_10_pyzx                              85     82     83     82     153     145
barenco_tof_3_pyzx                       21     21     21     21      31      29
barenco_tof_4_pyzx                       40     37     39     39      58      55
barenco_tof_5_pyzx                       60     58     60     59      82      81
barenco_tof_10_pyzx                     161    157    159    156     212     207
grover_5_pyzx                           262    256    259    261     342     340
ham15-low_pyzx                          227    226    227    227     293     294
ham15-med_pyzx                          440    435    437    438     602     602
ham15-high_pyzx                        1798   1771   1785   1780    2649    2612
mod5_4_pyzx                              17     17     17     17      24      22
mod_mult_55_pyzx                         45     45     45     45      74      71
mod_red_21_pyzx                          89     89     88     88     129     130
qcla_com_7_pyzx                         137    137    137    137     204     212
qcla_mod_7_pyzx                         301    299    300    300     484     480
rc_adder_6_pyzx                          65     65     65     65     101      97
vbe_adder_3_pyzx                         45     45     45     45      64      63
TOTAL                                  3869   3816   3842   3836    5638    5569
```

**pp(k=2) = 3816 is best** — 1.4% fewer CNOTs than pp(k=1) = 3869.
No circuit regresses (pp(k>1) ≤ pp(k=1) for every row, guaranteed by fallback).

### Comparison with Stage 5

| Method | Table-1 TOTAL CNOTs | vs pp(k=1) |
|--------|---------------------|------------|
| Stage 5: warm-start k=2 | 4468 | +15% WORSE |
| Stage 6: joint A* k=2   | 3816 |  −1.4% BETTER |

The joint A* eliminates the Stage-5 zero-sum trade by allowing a single CNOT
to reduce parities in both blocks simultaneously.

### Best per-circuit improvements (k=2 vs k=1)

| Circuit | pp(1) | pp(2) | Savings |
|---------|-------|-------|---------|
| barenco_tof_4_pyzx | 40 | 37 | −3 (7.5%) |
| grover_5_pyzx | 262 | 256 | −6 (2.3%) |
| ham15-high_pyzx | 1798 | 1771 | −27 (1.5%) |
| tof_10_pyzx | 85 | 82 | −3 (3.5%) |
| barenco_tof_10_pyzx | 161 | 157 | −4 (2.5%) |

---

## Extended benchmark (28 circuits, --set extended, --no-todd)

```
circuit                                pp(1)  pp(2)  pp(3)  pp(5)   mst+P  gray+P
Table-1 subtotal (19 circuits)         3869   3816   3842   3836    5638    5569
----
Adder8_pyzx                             157    156    157    156     225     230
adder_8_pyzx                            319    319    317    317     444     463
csla_mux_3_original_pyzx                 78     76     77     78     110     123
csum_mux_9_corrected_pyzx               140    138    139    139     187     189
hwb6_pyzx                               103    100    103    100     153     158
mod_adder_1024_pyzx                    1352   1346   1349   1347    2211    2145
nth_prime6_pyzx                         391    389    389    387     611     607
qcla_adder_10_pyzx                      195    194    193    194     290     281
qft_4_pyzx                               42     42     42     42      69      68
----
TOTAL (28 circuits)                    6646   6576   6608   6596    9938    9833
```

Joint A* k=2 achieves **6576 / 6646 = 1.1% fewer CNOTs** across 28 circuits.
pp(k=2) is ≤ pp(k=1) for every single circuit (fallback guarantee holds).

---

## Why k=2 beats k=3 and k=5

For k=3, each group of 3 blocks processes pair (B_{3i}, B_{3i+1}) jointly
and leaves B_{3i+2} as an independent singleton.  B_{3i+2} and B_{3i+3}
(start of the next group) are NOT paired jointly across the group boundary.
This reduces pairing opportunities compared to k=2, where every consecutive
pair is jointly synthesized.

---

## Files changed

| File | Change |
|---|---|
| `src/tableau/phasepoly/multiblock.cpp` | Replaced warm-start with joint A*: `synthesize_joint`, `JointState`, `finish_output0_joint`, `apply_h_transition`, `gf2_is_invertible`, updated `synthesize_block_group` |
| `src/tableau/phasepoly/multiblock.hpp` | No API change (same public functions) |
| `tests/src/tableau/phasepoly/multiblock.cpp` | Added Stage-6 tests (Tests 9–12, 4 cases / 5 assertions) |
| `tools/benchmark_phasepoly.cpp` | Updated notes to describe joint A* semantics |
| `doc/stage6.md` | This file |

---

## Test coverage

```
All tests passed (31 assertions in 12 test cases)
```

Tests 9–12 cover: shared-parity CNOT savings, no-regression guarantee,
infeasible H-boundary fallback, and global group-size comparison.

---

## Key correctness properties

1. **No regression**: joint solution is always compared to independent; the
   minimum is returned.  `pp(k>1) ≤ pp(k=1)` holds for every circuit.

2. **H-boundary validity**: the feasibility check ensures M_new is invertible
   before joint synthesis; infeasible pairs fall back gracefully.

3. **Cross-block benefit captured**: `finish_output0_joint` applies
   output0-restoring CNOTs to output1 and phase1 as well, capturing additional
   cross-block benefit even for deterministic output synthesis.

---

## Limitations and future work

1. **Budget sensitivity**: the joint A\* state is larger (includes phase1 in
   the key), so the search explores fewer unique phase0 configurations within
   the same budget.  Increasing `max_expansions` could yield larger improvements.

2. **k > 2 joint**: a true 3-block or 5-block joint A\* (processing all blocks
   simultaneously) would capture more cross-block structure but has an even
   larger state space.  Not implemented here.

3. **Feasibility rate**: for some block pairs the Q-matrix is singular, forcing
   fallback to independent synthesis.  This limits the fraction of pairs that
   benefit from joint synthesis.

4. **Phase-dominated gains**: the improvement is largest for circuits with
   consecutive blocks that share parity columns (barenco_tof family, grover).
   Circuits where parities are unrelated across blocks (ham15) see smaller gains.
