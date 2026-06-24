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

### Fallback (correctness only)

Infeasible pairs (H-boundary singularity) and budget-exhausted pairs fall back
to independent A*.  The joint A* result is otherwise reported as-is — it may
exceed pp(k=1) if the joint state space causes the search to miss the optimum
within the expansion budget.

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
synthesized).  Whether this gives the best results depends on budget sensitivity
(see empirical results below).

---

## Empirical results — Table 1 benchmark (--set table1, --no-todd)

Honest results: joint A* result returned as-is; no comparison against
independent synthesis.  Values are CNOT counts only.

```
circuit                                  pp(1)   pp(2)   pp(3)   pp(5)   mst+P  gray+P
tof_3_pyzx                                 16      16      16      16      28      26
tof_4_pyzx                                 25      26      29      26      45      43
tof_5_pyzx                                 35      40      39      38      63      60
tof_10_pyzx                                85      89      93      91     153     145
barenco_tof_3_pyzx                         21      22      21      22      31      29
barenco_tof_4_pyzx                         40      37      40      40      58      55
barenco_tof_5_pyzx                         60      60      62      59      82      81
barenco_tof_10_pyzx                       161     159     164     159     212     207
grover_5_pyzx                             262     272     268     270     342     340
ham15-low_pyzx                            227     255     255     252     293     294
ham15-med_pyzx                            440     525     474     490     602     602
ham15-high_pyzx                          1798    1991    1925    1954    2649    2612
mod5_4_pyzx                                17      17      17      17      24      22
mod_mult_55_pyzx                           45      53      52      51      74      71
mod_red_21_pyzx                            89      96      91      95     129     130
qcla_com_7_pyzx                           137     153     149     156     204     212
qcla_mod_7_pyzx                           301     364     338     347     484     480
rc_adder_6_pyzx                            65      79      77      78     101      97
vbe_adder_3_pyzx                           45      49      45      48      64      63
TOTAL                                    3869    4303    4155    4209    5638    5569
```

**Honest finding**: joint A* is **worse** for most circuits.

| Method | TOTAL CNOTs | vs pp(k=1) |
|--------|-------------|------------|
| pp(k=1) independent | 3869 | baseline |
| pp(k=3) | 4155 | +7.4% WORSE |
| pp(k=5) | 4209 | +8.8% WORSE |
| pp(k=2) | 4303 | +11.2% WORSE |
| Stage 5: warm-start k=2 | 4468 | +15.4% WORSE |

Joint A* eliminates the warm-start zero-sum trade (Stage 5) but introduces a
**budget sensitivity problem**: the joint state space (including phase1 in the
state key) is larger, causing A* to explore fewer unique configurations within
the same expansion budget.  Many joint pairs find a valid but suboptimal
solution — one that exceeds what independent A* would produce.

Notably, k=3 slightly outperforms k=2: since every group of 3 leaves the last
block as an independent singleton, fewer pairs are jointly processed, reducing
exposure to the joint A*'s suboptimality.

### Circuits that improve with joint A* (k=2)

| Circuit | pp(1) | pp(2) | Savings |
|---------|-------|-------|---------|
| barenco_tof_4_pyzx | 40 | 37 | −3 (7.5%) |
| barenco_tof_10_pyzx | 161 | 159 | −2 (1.2%) |

These circuits have consecutive blocks sharing parity columns — the joint
search genuinely benefits from a CNOT that cancels parities in both blocks.

### Circuits that regress most (k=2)

| Circuit | pp(1) | pp(2) | Regression |
|---------|-------|-------|------------|
| qcla_mod_7_pyzx | 301 | 364 | +63 (+20.9%) |
| ham15-med_pyzx | 440 | 525 | +85 (+19.3%) |
| ham15-high_pyzx | 1798 | 1991 | +193 (+10.7%) |

These circuits have many blocks with unrelated parities; the joint A* wastes
budget considering cross-block interactions that yield no benefit.

### Comparison with Stage 5

| Method | Table-1 TOTAL CNOTs | vs pp(k=1) |
|--------|---------------------|------------|
| Stage 5: warm-start k=2 | 4468 | +15.4% WORSE |
| Stage 6: joint A* k=3   | 4155 |  +7.4% WORSE |
| Stage 6: joint A* k=5   | 4209 |  +8.8% WORSE |
| Stage 6: joint A* k=2   | 4303 | +11.2% WORSE |

Joint A* is better than warm-start but still worse than independent synthesis
with the current budget (max_exp=100000).  Increasing the budget or pruning
the joint state space would be needed to achieve net gains.

---

## Extended benchmark (28 circuits, --set extended, --no-todd)

Honest results (no comparison against independent synthesis).

```
circuit                                  pp(1)   pp(2)   pp(3)   pp(5)   mst+P  gray+P
Table-1 subtotal (19 circuits)           3869    4303    4155    4209    5638    5569
----
Adder8_pyzx                              157     172     172     170     225     230
adder_8_pyzx                             319     350     334     338     444     463
csla_mux_3_original_pyzx                  78      86      82      87     110     123
csum_mux_9_corrected_pyzx                140     149     146     144     187     189
hwb6_pyzx                                103     103     113     107     153     158
mod_adder_1024_pyzx                     1352    1602    1504    1538    2211    2145
nth_prime6_pyzx                          391     438     421     426     611     607
qcla_adder_10_pyzx                       195     210     206     218     290     281
qft_4_pyzx                                42      48      44      46      69      68
----
TOTAL (28 circuits)                      6646    7461    7177    7283    9938    9833
```

| Method | TOTAL CNOTs | vs pp(k=1) |
|--------|-------------|------------|
| pp(k=1) | 6646 | baseline |
| pp(k=3) | 7177 | +8.0% WORSE |
| pp(k=5) | 7283 | +9.6% WORSE |
| pp(k=2) | 7461 | +12.3% WORSE |

Budget sensitivity dominates across all circuit families.  `hwb6_pyzx` is the
only circuit where pp(k=2) = pp(k=1) = 103 (flat).  All other circuits regress.

---

## Why k=3 outperforms k=2 (under tight budget)

For k=3, each group of 3 blocks processes pair (B_{3i}, B_{3i+1}) jointly
and leaves B_{3i+2} as an independent singleton.  Since the joint A* is
budget-sensitive and often produces suboptimal results, having 1/3 of all
blocks processed independently (as in k=3) reduces the damage.

k=2 pairs every consecutive pair jointly — maximising exposure to suboptimality.
k=3 leaves every third block independent — a natural hedge.
k=5 behaves like k=3 for groups of 5 (same trailing-singleton pattern).

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
All tests passed (33 assertions in 12 test cases)
```

Tests 9–12 cover: shared-parity CNOT savings, joint A* completion (honest,
no regression guarantee), infeasible H-boundary fallback, and global
group-size comparison.

---

## Key correctness properties

1. **H-boundary validity**: the feasibility check ensures M_new is invertible
   before joint synthesis; infeasible pairs fall back to independent synthesis
   gracefully.

2. **Budget-exhausted fallback**: when the joint A* exhausts max_exp without
   finding any complete solution, it falls back to independent synthesis.
   This preserves correctness but the joint result (when found within budget)
   is returned as-is — it may be worse than independent synthesis.

3. **Cross-block benefit captured**: `finish_output0_joint` applies
   output0-restoring CNOTs to output1 and phase1 as well, capturing additional
   cross-block benefit even for deterministic output synthesis.

---

## Limitations and future work

1. **Budget sensitivity (primary bottleneck)**: the joint A\* state space is
   ~n× larger than independent (includes phase1 in the key), so the same
   expansion budget (max_exp=100000) explores a much smaller fraction of the
   optimal path.  Empirically, this causes net regression in 17/19 Table-1
   circuits.  Increasing `max_expansions` substantially (e.g. 10×) could
   recover the gains, at proportional runtime cost.

2. **k > 2 joint**: a true 3-block or 5-block joint A\* (processing all blocks
   simultaneously) would capture more cross-block structure but has an even
   larger state space.  Not implemented here.

3. **Feasibility rate**: for some block pairs the Q-matrix is singular, forcing
   fallback to independent synthesis.  This limits the fraction of pairs that
   can benefit from joint synthesis at all.

4. **Phase-dominated gains**: improvement requires consecutive blocks to share
   parity columns (barenco_tof family).  Circuits with unrelated parities
   across blocks (ham15, qcla_mod_7) see the largest regressions.
