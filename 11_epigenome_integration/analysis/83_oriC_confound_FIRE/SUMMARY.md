# C18: Is the methylome–FIRE / Exposed–FIRE signal an oriC-proximity confound? (2026-06-29)

oriC ≈ dnaA (SCO3879) midpoint = 4,271,763 bp. Seed 42.

## Result — the confound is LARGELY REAL
- Windows Spearman(meth_T1, FIRE) raw **ρ=0.319** (p=5.6e-42) → **partial ρ=0.026** controlling distance-to-oriC = collapses to ~0.
- FIRE vs distance-to-oriC: **ρ=−0.934** — FIRE is almost entirely an oriC-distance gradient.
- meth_T1 vs distance-to-oriC: ρ=−0.332.
- Exposed vs Shielded FIRE_tss overall: medE 1.27 vs medS 0.94, p=8.6e-6 — BUT **within oriC-distance quartiles, NOT significant in any (Q1 p=0.42, Q2 p=0.13, Q3 p=0.055, Q4 p=0.46)**.

## Interpretation (honest)
The methylome–FIRE correlation and the Exposed-higher-FIRE difference are **largely explained by shared proximity to oriC** (both methylation density and FIRE peak in the oriC-centred core). Controlling distance-to-oriC removes the methylome–FIRE partial correlation (0.32→0.03) and the Exposed/Shielded FIRE gap (vanishes within matched bins). NOTE: the manuscript's earlier control was only the coarse core/arm **binary** (left partial r=0.18); the continuous oriC-distance is a much stronger control and FIRE tracks it at ρ=−0.934.

## Implication for the paper
The "vegetative methylome independently tracks FIRE / could pre-screen integration neighbourhoods without Hi-C" sub-claim (Discussion 3D para, Supp Fig 8, the engineering proxy) is **not supported as an independent coupling** — it is largely the oriC/replication-compartment gradient that both signals share. This is CONSISTENT with the central permissive/spatial thesis (methylation, FIRE, and compartment all track one oriC-centred core architecture) but UNDERCUTS the standalone methylome→FIRE predictive/engineering claim.
