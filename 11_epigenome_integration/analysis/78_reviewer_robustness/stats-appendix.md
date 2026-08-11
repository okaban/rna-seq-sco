# Stats appendix — Phase-1 reviewer-robustness

Seed = 42 throughout. Software: Python 3.x, pandas 2.3.3, scipy 1.17.0, statsmodels 0.14.6, diptest 0.11.0.

## Units & samples
- Regulatory-gene table: 1,051 genes (62 Exposed / 989 Shielded by methylation map); **1,019 with a
  DESeq2 LFC_T2vsT1** (59 Exposed / 960 Shielded) — 32 dropped by DESeq2 independent filtering (low count).
  This matches the manuscript's "59 with a defined fold-change".
- Genome-wide dose-response: 7,496 genes with GCCGGC dose group + LFC + position.
- Per-read: T1 replicate 1 (`1-1_modonly.tsv`, 8.24 M mod rows), high-confidence prob ≥ 0.75, fail≠true.
- Pileup QC: T1 replicate 1, positions with coverage ≥ 10 (6.20 M cytosine positions).

## P1-1 cross-talk
- Statistic: per-offset 4mC co-call rate around high-conf 6mA = (#4mC hits at offset set)/(|offsets|×n_6mA).
  NEAR offsets {±3,±4}; FAR offsets {±20…±50}. Reported separately for AAGCCCG vs non-AAGCCCG 6mA.
- Motif specificity = AAGCCCG near-rate / non-AAGCCCG near-rate = 430×. Generic near/far = 1.09 (no bleed).
- Assumption/limitation: call_prob-level analysis cannot fully exclude signal cross-talk; the motif-specificity
  contrast (vs 2.1 M generic controls) is the strongest available evidence short of single-mod re-calling.
  OR=138,440 is from a saturated 2×2 (off-diagonal=2) and is not reported as an effect size.

## P1-2 m4C vs 5mC
- Global modified fraction = Σ Nmod / Σ Nvalid over cov≥10 positions. 5mC 0.0014 %, 4mC 0.0226 %.
- Single replicate; qualitative result stable across replicates.

## P1-3 / P1-3b spatial autocorrelation
- Null: circular shift (np.roll) of the position-ordered LFC vector, n=10,000; preserves 1-D spatial
  autocorrelation, breaks association. Two-sided empirical p = (#|null|≥|obs| + 1)/(N+1).
- Genome-wide dose-response: obs ρ=−0.087; block-perm p=0.219 (rho), 0.225 (KW H=61.03). Null |ρ| 95th=0.104.
- Regulatory weak bias (distance proxy): partial Spearman (residual-rank, control core/arm) ρ=+0.128;
  block-perm p=0.0021. Null |ρ| 95th=0.079. Raw ρ=+0.183, block-perm p=0.0032.
- Caveat: P1-3b uses log nearest-distance as the methylation proxy (canonical −0.09 uses ±2 kb occupancy);
  direction identical. Direct occupancy re-test is the optional confirmation.

## P1-4 equivalence (TOST)
- Cliff's δ via vectorised rank counts; 95% CI by stratified bootstrap (10,000 resamples, percentile).
- Welch TOST: two one-sided t-tests on (mean_E − mean_S) of |LFC|, equivalence margin ±0.5 log2, 90% CI.
- |LFC_T2vsT1|: δ=−0.010 CI[−0.164,+0.147]; TOST p=0.010 (equivalent). |LFC_T3vsT1|: δ=+0.091
  CI[−0.050,+0.232]; TOST p=0.016. DEG-proportion TOST (margin ±0.10): p=0.074 (borderline).
- Honest bound: medium/large differences excluded; small difference not excluded (δ CI just exceeds 0.147).

## P1-5 continuous / threshold
- Hartigan dip on log10(nearest_GCCGGC_T1+1): D=0.0078, p=0.981 (unimodal; reproduces D=0.007).
- Partial Spearman residual-rank, control core/arm: (dist, signed LFC) +0.128 p=4e-5; (dist, |LFC|) +0.008 p=0.81.
- Threshold sweep 100–500 bp (50 bp steps): Exposed-vs-Shielded |LFC| MWU p∈[0.51,0.96], Cliff δ∈[−0.10,+0.04].

## Blockers / not done here
- Single-mod Dorado re-calling (definitive cross-talk resolution) needs raw POD5 + GPU — not run; motif-specificity
  test substitutes. m4C orthogonal chemical validation out of scope (sequencing-only study).
- Occupancy-based spatial-permutation of the canonical −0.09 pending (distance proxy used).
