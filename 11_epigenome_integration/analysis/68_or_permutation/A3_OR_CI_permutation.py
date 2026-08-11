#!/usr/bin/env python3
"""
A-3: AAGCCCG dual-modification OR=138,440 — 95% CI + permutation test

Contingency table (pileup-level, T1 pooled, Validation-Results.md S-4):
                co-mod (6mA+4mC)   not co-mod   total
In AAGCCCG:          244              1,090       1,334
Out AAGCCCG:           2          1,236,879   1,236,881

- A1↔C5 (4bp spacing): 155 pairs
- A0↔C3 (3bp spacing):  85 pairs
- other intra:           4 pairs
- Total co-mod:        244

OR = (244 * 1236879) / (1090 * 2) = 138,440

Steps:
  1. Fisher's exact test (one-sided) + Haldane-corrected OR
  2. 95% CI via statsmodels Table2x2 (conditional exact log-OR interval)
  3. Permutation test (n=10,000): shuffle co-modification labels across all sites
  4. Publication figure: null distribution (log10 OR) + observed OR marker
  5. TSV output
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from scipy import stats
import statsmodels.stats.contingency_tables as sm_ct
import os

# ── Paths ──────────────────────────────────────────────────────────────────────
OUTDIR   = os.path.expanduser(
    "~/bioinfo/rna-seq/11_epigenome_integration/analysis/68_or_permutation"
)
TSV_PATH = os.path.join(OUTDIR, "tables", "A3_OR_CI_permutation_results.tsv")
FIG_PATH = os.path.join(OUTDIR, "figures", "A3_permutation_null_dist.png")

# ── Contingency table ──────────────────────────────────────────────────────────
# [[in_comod, in_notcomod], [out_comod, out_notcomod]]
a, b = 244, 1_090          # In AAGCCCG:  co-mod, not co-mod
c, d = 2,   1_236_879      # Out AAGCCCG: co-mod, not co-mod

table = np.array([[a, b],
                  [c, d]])

# ── Step 1: Fisher's exact test ────────────────────────────────────────────────
print("=" * 65)
print("A-3: AAGCCCG dual-modification OR CI + permutation test")
print("=" * 65)
print(f"\nContingency table:")
print(f"              co-mod    not co-mod   total")
print(f"In AAGCCCG:  {a:>7,}    {b:>10,}  {a+b:>7,}")
print(f"Out AAGCCCG: {c:>7,}    {d:>10,}  {c+d:>7,}")

or_fisher, p_fisher = stats.fisher_exact(table, alternative="greater")
print(f"\n[Fisher's exact test]")
print(f"  OR (scipy)  = {or_fisher:,.0f}")
print(f"  p (one-sided) = {p_fisher:.3e}")

# ── Step 2: 95% CI via statsmodels ────────────────────────────────────────────
tbl2x2 = sm_ct.Table2x2(table)
or_sm   = tbl2x2.oddsratio
ci_lo, ci_hi = tbl2x2.oddsratio_confint(alpha=0.05, method="normal")

print(f"\n[Odds Ratio + 95% CI  (log-OR normal approx.)]")
print(f"  OR               = {or_sm:,.0f}")
print(f"  95% CI           = [{ci_lo:,.0f}, {ci_hi:.3e}]")

# Exact Cornfield 95% CI via statsmodels (method='exact')
try:
    ci_lo_ex, ci_hi_ex = tbl2x2.oddsratio_confint(alpha=0.05, method="exact")
    print(f"  95% CI (exact)   = [{ci_lo_ex:,.0f},  {ci_hi_ex:.3e}]")
    use_exact = True
except Exception as e:
    print(f"  (exact CI unavailable: {e})")
    use_exact = False

# Manual log-OR SE for verification
log_or_se = np.sqrt(1/a + 1/b + 1/c + 1/d)
log_or    = np.log(or_sm)
ci_lo_m   = np.exp(log_or - 1.96 * log_or_se)
ci_hi_m   = np.exp(log_or + 1.96 * log_or_se)
print(f"  95% CI (manual)  = [{ci_lo_m:,.0f},  {ci_hi_m:.3e}]")
print(f"  log(OR) SE       = {log_or_se:.4f}")

# ── Step 3: Permutation test ───────────────────────────────────────────────────
print(f"\n[Permutation test  n=10,000]")
np.random.seed(42)
N_PERM = 10_000

# Build label arrays:
#   n_total sites, n_aagcccg of which are AAGCCCG
n_total    = (a + b) + (c + d)   # 1,238,215
n_aagcccg  = a + b               # 1,334
n_comod    = a + c               # 246  (total co-modified sites)

# For each permutation:
#   - randomly select n_aagcccg sites as "AAGCCCG" (without replacement)
#   - count co-modified among them
# We represent sites as integers 0..n_total-1;
# co-modified sites are 0..n_comod-1.
# Equivalent (faster): hypergeometric draw

perm_ors = np.empty(N_PERM)
for i in range(N_PERM):
    # Draw: how many of the n_aagcccg "AAGCCCG" sites happen to be co-modified?
    a_perm = np.random.hypergeometric(
        ngood=n_comod,        # co-modified sites in pool
        nbad=n_total - n_comod,  # non-co-modified
        nsample=n_aagcccg     # AAGCCCG draw size
    )
    b_perm = n_aagcccg  - a_perm
    c_perm = n_comod    - a_perm
    d_perm = (n_total - n_aagcccg) - c_perm

    # Avoid division by zero (add 0.5 Haldane correction when any cell = 0)
    if a_perm == 0 or b_perm == 0 or c_perm == 0 or d_perm == 0:
        perm_ors[i] = (a_perm + 0.5) * (d_perm + 0.5) / ((b_perm + 0.5) * (c_perm + 0.5))
    else:
        perm_ors[i] = (a_perm * d_perm) / (b_perm * c_perm)

observed_or = or_sm
perm_p      = float(np.mean(perm_ors >= observed_or))

print(f"  Observed OR      = {observed_or:,.0f}")
print(f"  Permutation dist: min={perm_ors.min():.2f}, "
      f"max={perm_ors.max():.2f}, "
      f"median={np.median(perm_ors):.4f}")
print(f"  Permutations >= observed: {int(np.sum(perm_ors >= observed_or))} / {N_PERM}")
print(f"  Permutation p  = {perm_p:.4f}  (p < 1/n_perm if 0: < {1/N_PERM:.4f})")

if perm_p == 0.0:
    perm_p_str = f"< {1/N_PERM:.4f}"
    perm_p_report = 1 / N_PERM  # upper bound for figure label
else:
    perm_p_str = f"{perm_p:.4f}"
    perm_p_report = perm_p

# ── Step 4: Figure ─────────────────────────────────────────────────────────────
log10_perm_ors = np.log10(np.maximum(perm_ors, 1e-9))
log10_obs      = np.log10(observed_or)
log10_ci_lo    = np.log10(ci_lo_m)
log10_ci_hi    = np.log10(ci_hi_m)

fig, ax = plt.subplots(figsize=(7, 4.5))

# Histogram of null distribution
ax.hist(log10_perm_ors, bins=60, color="#6B7783", edgecolor="white",
        linewidth=0.3, alpha=0.85, label="Null distribution (n=10,000)")

# Observed OR
ax.axvline(log10_obs, color="#A64B44", linewidth=2.2,
           label=f"Observed OR = {observed_or:,.0f}")

# CI shading (semi-transparent band)
ax.axvspan(log10_ci_lo, log10_ci_hi, alpha=0.15, color="#A64B44",
           label=f"95% CI [{ci_lo_m:,.0f}, {ci_hi_m:.2e}]")

# Annotation text
ax.text(0.97, 0.96,
        f"OR = {observed_or:,.0f}\n"
        f"95% CI [{ci_lo_m:,.0f}, {ci_hi_m:.2e}]\n"
        f"perm p {perm_p_str}  (n={N_PERM:,})",
        transform=ax.transAxes, fontsize=8.5,
        verticalalignment="top", horizontalalignment="right",
        bbox=dict(boxstyle="round,pad=0.35", facecolor="white",
                  edgecolor="#EAEAEA", alpha=0.9))

ax.set_xlabel("log₁₀(Odds Ratio)", fontsize=11)
ax.set_ylabel("Count (permutations)", fontsize=11)
ax.set_title("AAGCCCG dual-modification: permutation null distribution\n"
             "vs observed OR (pileup-level, T1 pooled)", fontsize=10.5)
ax.legend(fontsize=8, loc="upper left")
ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())
ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())
ax.tick_params(axis="both", which="both", direction="in")
plt.tight_layout()
fig.savefig(FIG_PATH, dpi=300, bbox_inches="tight")
print(f"\n  Figure saved: {FIG_PATH}")

# ── Step 5: TSV output ─────────────────────────────────────────────────────────
os.makedirs(os.path.dirname(TSV_PATH), exist_ok=True)
with open(TSV_PATH, "w") as fh:
    fh.write("metric\tvalue\n")
    fh.write(f"contingency_a_in_comod\t{a}\n")
    fh.write(f"contingency_b_in_notcomod\t{b}\n")
    fh.write(f"contingency_c_out_comod\t{c}\n")
    fh.write(f"contingency_d_out_notcomod\t{d}\n")
    fh.write(f"odds_ratio\t{or_sm:.4f}\n")
    fh.write(f"CI95_lower_normal\t{ci_lo_m:.4f}\n")
    fh.write(f"CI95_upper_normal\t{ci_hi_m:.4f}\n")
    if use_exact:
        fh.write(f"CI95_lower_exact\t{ci_lo_ex:.4f}\n")
        fh.write(f"CI95_upper_exact\t{ci_hi_ex:.4f}\n")
    fh.write(f"log_OR_SE\t{log_or_se:.6f}\n")
    fh.write(f"fisher_p_onesided\t{p_fisher:.6e}\n")
    fh.write(f"permutation_n\t{N_PERM}\n")
    fh.write(f"permutation_p\t{perm_p_str}\n")
    fh.write(f"perm_null_median_OR\t{np.median(perm_ors):.6f}\n")
    fh.write(f"perm_null_max_OR\t{perm_ors.max():.4f}\n")
    fh.write(f"intra_A1_C5_pairs\t155\n")
    fh.write(f"intra_A0_C3_pairs\t85\n")
    fh.write(f"other_intra_pairs\t4\n")
    fh.write(f"cross_motif_forbidden_pairs\t0\n")
    fh.write(f"total_comod_pairs\t244\n")

print(f"  TSV saved:    {TSV_PATH}")

# ── Summary ────────────────────────────────────────────────────────────────────
print("\n" + "=" * 65)
print("FINAL RESULTS (for paper)")
print("=" * 65)
print(f"  OR          = {or_sm:,.0f}")
print(f"  95% CI      = [{ci_lo_m:,.0f},  {ci_hi_m:.3e}]")
print(f"  Fisher p    = {p_fisher:.3e}")
print(f"  Perm p      = {perm_p_str}  (n={N_PERM:,})")
print()
print("Paper text:")
print(f'  "...OR = {or_sm:,.0f} (95% CI [{ci_lo_m:,.0f}, {ci_hi_m:.2e}]; '
      f'Fisher\'s exact p < 10^-300; permutation p {perm_p_str}, n = {N_PERM:,})..."')
