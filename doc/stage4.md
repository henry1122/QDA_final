# Stage 4 — CLI Integration & Benchmark Harness

## Goal

Wire the PhasePoly A\* synthesizer (Stages 0–3) into qsyn as a first-class
`PauliRotationsSynthesisStrategy` so it can be invoked from the command line and
compared directly against existing baselines on standard benchmarks.

---

## What was implemented

### 1. `PhasePolySynthesisStrategy` (`src/tableau/phasepoly/strategy.hpp/.cpp`)

A new class deriving from `qsyn::experimental::PauliRotationsSynthesisStrategy`
that adapts the A\* solver to the Tableau pipeline's interface.

**Input** (from `to_qcir` pipeline):
A `vector<PauliRotation>` — the pure phase-polynomial segment between two
`StabilizerTableau` Clifford layers in the Tableau IR.  All rotations are
diagonal (Z-only) by construction.

**Conversion to `PhasePolyProblem`**:
- **Output matrix = identity.**  The flanking Clifford layers handle all basis
  transformations; the phase-polynomial block only needs to realize the phase
  terms and return to the same basis.
- **Deduplication.**  Identical Z-support patterns are merged by accumulating
  their angles, then zero-angle (mod 2π) terms are dropped — same normalization
  as `phase_block_to_problem`.
- The parity column for rotation `j` is read directly from `PauliRotation::is_z(q)`.

**Then** calls `synthesize_phasepoly(problem, _config)` and wraps the result
with `build_qcir`.

### 2. CLI registration (`src/cmd/conversion_cmd.cpp`)

Added `"phasepoly"` to the `--rotation` strategy choices of
`convert tableau qcir`:

```
convert tableau qcir --rotation phasepoly
```

alongside the existing `naive`, `tpar`, `graysynth`, `gstair`, `mst`.

Full benchmark pipeline from the qsyn CLI:

```
qcir read benchmark/qc/optimized/tof_3_pyzx.qc
convert qcir tableau
convert tableau qcir --rotation phasepoly
qcir print --stat
quit -f
```

### 3. Benchmark script (`scripts/benchmark_phasepoly.py`)

A Python script that runs the full pipeline for every `.qc` file in
`benchmark/qc/optimized/` across all five strategies and prints a comparison
table of 2-qubit (CX) and T-family gate counts.

```bash
python scripts/benchmark_phasepoly.py \
    --qsyn ./build/qsyn \
    --bench ./benchmark/qc/optimized \
    --out results.csv
```

Optional arguments:
- `--strategies naive graysynth mst phasepoly` — select a subset
- `--out FILE` — write a CSV alongside the table

---

## Tests (`tests/src/tableau/phasepoly/strategy.cpp`, 6 tests)

| Test | What it checks |
|---|---|
| `trivial identity circuit` | CX;CX cancels → round-trip through Tableau is equivalent |
| `single T-gate` | 1-qubit single Rz passes through unchanged |
| `CX; Rz; CX round-trip` | 2-qubit phase term synthesized correctly |
| `two phase terms on 3 qubits` | full pi/2(x⊕y) + pi/4(y⊕z) example |
| `circuit with H gate splits into two blocks` | H boundary → two independent phase blocks each synthesized |
| `CNOT count ≤ MST + 2` | sanity-bound: our result is not drastically worse than MST on a simple 3-qubit circuit |

All 35 phasepoly tests (Stages 0–4) pass.

---

## How output matrix = identity is justified

`to_qcir(Tableau, st_strategy, pr_strategy)` processes each `SubTableau` in
order.  A `vector<PauliRotation>` sub-tableau is sandwiched between two
`StabilizerTableau` parts.  The Clifford layers emit basis-change gates
independently; the phase-polynomial layer must realize its phase terms and
leave the qubit basis unchanged.  Setting `output_matrix = I_n` enforces this
and the A\* search finds the minimum-CX circuit that satisfies both the phase
terms and the identity output constraint.

---

## Comparison baselines available

| Flag | Algorithm | Paper reference |
|---|---|---|
| `naive` | Independent Rz diagonalisation | — |
| `graysynth` | Gray-Synth star | Amy & Mosca 2016 |
| `gstair` | Gray-Synth staircase | Amy & Mosca 2016 |
| `mst` | Vandaele MST | Vandaele et al. 2021 — **Table 1 baseline** |
| `phasepoly` | PhasePoly A\* (this work) | Chen et al. 2025 |

The paper's primary baseline is **`mst`** (Table 1).  `graysynth` is the
secondary baseline.

---

## Next step (Stage 5)

Multi-block co-optimization: currently each `vector<PauliRotation>` segment in
the Tableau is synthesized independently.  The paper's §4 shows that merging
adjacent segments across Clifford boundaries can expose additional sharing
opportunities.
