"""C5/C6/C7 — FIRE x GCCGGC methylation hotspots, characterisation, resolution.

Reviewer requests:
  C5  HOTSPOTS: discrete windows that are simultaneously high-methylation (T1
      GCCGGC) AND high-FIRE (top quantile of both). Ranked candidate BGC-
      integration neighbourhoods, with coordinates, FIRE, methylation density,
      nearest gene/product, Compartment-A membership, distance to known BGCs.
  C6  WHAT ARE they: gene density, functional categories, compartment membership,
      class enrichment of the high-FIRE/high-methylation set.
  C7  RESOLUTION: at what spatial resolution the methylome can nominate insertion
      sites. Window-size sweep, candidate-window counts, methylation-vs-FIRE
      concordance (precision/recall/overlap/Jaccard) at multiple thresholds.

Inputs (read-only):
  - 76_FIRE_methylation_crossref/tables/fire_methyl_bins.tsv  (pre-joined 5-kb bins)
  - 37_defense_island_GCCGGC/tables/GCCGGC_sites_by_timepoint.tsv  (raw T1 sites)
  - 05_annotation .../gene_annotation_basic.tsv
  - 47_BGC_methylation_geographic_test/tables/BGC_gene_geography.tsv

FIRE source: Deng et al. 2023 PNAS sd04 (M-phase ~ our T1). Genome-wide 5-kb FIRE.
Compartment A = 2.3-6.2 Mb; chromosomal core = 1.5-7.17 Mb.
All values chromosome NC_003888.3 only (FIRE/methylome are chromosome-scale).
"""
import numpy as np
import pandas as pd
from pathlib import Path

BASE = Path("/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis")
OUT = BASE / "78_reviewer_robustness"
TAB = OUT / "tables"
TAB.mkdir(parents=True, exist_ok=True)

CHROM = "NC_003888.3"
COMP_A_LO, COMP_A_HI = 2_300_000, 6_200_000   # Compartment A
CORE_LO, CORE_HI = 1_500_000, 7_170_000       # chromosomal core
BIN = 5000                                    # native FIRE resolution

# Known antiSMASH BGCs (min start, max end of member genes)
BGCS = [
    ("cda", 3_519_449, 3_602_320),
    ("act", 5_513_809, 5_535_091),
    ("red", 6_432_812, 6_464_206),
    ("cpk", 6_900_898, 6_948_414),
]


# ----------------------------------------------------------------------------
# Load pre-joined 5-kb FIRE x methylation table
# ----------------------------------------------------------------------------
bins = pd.read_csv(BASE / "76_FIRE_methylation_crossref/tables/fire_methyl_bins.tsv",
                   sep="\t")
bins = bins[bins["FIRE_M"].notna()].copy()
bins["center"] = bins["center"].astype(float)
bins["start_bp"] = (bins["bin"] * BIN).astype(int)
bins["end_bp"] = bins["start_bp"] + BIN
bins["in_compA"] = ((bins["center"] >= COMP_A_LO) &
                    (bins["center"] <= COMP_A_HI)).astype(int)
bins["in_core"] = ((bins["center"] >= CORE_LO) &
                   (bins["center"] <= CORE_HI)).astype(int)

# Gene annotation (chromosome only)
ann = pd.read_csv(BASE.parent.parent / "05_annotation/analysis/"
                  "05_annotation_260128_v1/tables/gene_annotation_basic.tsv",
                  sep="\t")
ann = ann[ann["contig"] == CHROM].copy()
ann["mid"] = (ann["start"] + ann["end"]) / 2.0
ann = ann.sort_values("mid").reset_index(drop=True)
gene_mid = ann["mid"].to_numpy()


def nearest_gene(center: float):
    """Return (gene_id, product, distance_bp) of gene whose midpoint is closest."""
    i = int(np.argmin(np.abs(gene_mid - center)))
    row = ann.iloc[i]
    return row["gene_id"], row["product"], abs(row["mid"] - center)


def genes_in_window(start: int, end: int) -> pd.DataFrame:
    """Genes overlapping [start, end)."""
    return ann[(ann["end"] > start) & (ann["start"] < end)]


def dist_to_nearest_bgc(center: float):
    """0 if inside a BGC, else bp to nearest BGC edge; plus nearest BGC name."""
    best_d, best_n = np.inf, None
    for name, s, e in BGCS:
        if s <= center <= e:
            return 0.0, name
        d = min(abs(center - s), abs(center - e))
        if d < best_d:
            best_d, best_n = d, name
    return best_d, best_n


# ----------------------------------------------------------------------------
# C5 — HOTSPOTS: top-quantile FIRE AND top-quantile methylation
# ----------------------------------------------------------------------------
# Methylation density per bin = T1 GCCGGC site count (m_n_T1). Use site count as
# the primary density metric; freq-sum kept as secondary.
QF = 0.80   # top-20% FIRE
QM = 0.80   # top-20% methylation (among methylated bins, to avoid zero-inflation)

fire_thr = bins["FIRE_M"].quantile(QF)
# methylation threshold computed over methylated bins only (m_n_T1>0)
methyl_pos = bins.loc[bins["m_n_T1"] > 0, "m_n_T1"]
methyl_thr = methyl_pos.quantile(QM)

bins["hi_fire"] = (bins["FIRE_M"] >= fire_thr).astype(int)
bins["hi_methyl"] = (bins["m_n_T1"] >= methyl_thr).astype(int)
bins["hotspot"] = ((bins["hi_fire"] == 1) & (bins["hi_methyl"] == 1)).astype(int)

hot = bins[bins["hotspot"] == 1].copy()

# annotate hotspots
rows = []
for _, r in hot.iterrows():
    c = r["center"]
    gid, prod, gdist = nearest_gene(c)
    bgc_d, bgc_n = dist_to_nearest_bgc(c)
    win_genes = genes_in_window(r["start_bp"], r["end_bp"])
    rows.append({
        "chrom": CHROM,
        "start_bp": int(r["start_bp"]),
        "end_bp": int(r["end_bp"]),
        "center_bp": int(c),
        "FIRE_M": round(float(r["FIRE_M"]), 4),
        "PC_M": round(float(r["PC_M"]), 4),
        "methyl_n_T1": int(r["m_n_T1"]),
        "methyl_freqsum_T1": round(float(r["m_fq_T1"]), 1),
        "n_genes_in_window": int(len(win_genes)),
        "nearest_gene": gid,
        "nearest_product": prod,
        "nearest_gene_dist_bp": int(gdist),
        "in_compartmentA": int(r["in_compA"]),
        "in_core": int(r["in_core"]),
        "dist_to_BGC_bp": int(bgc_d) if np.isfinite(bgc_d) else -1,
        "nearest_BGC": bgc_n,
    })
hot_df = pd.DataFrame(rows)
# rank: FIRE percentile + methylation percentile (combined evidence)
hot_df["FIRE_pct"] = hot_df["FIRE_M"].rank(pct=True)
hot_df["methyl_pct"] = hot_df["methyl_n_T1"].rank(pct=True)
hot_df["combined_score"] = (hot_df["FIRE_pct"] + hot_df["methyl_pct"]) / 2.0
hot_df = hot_df.sort_values("combined_score", ascending=False).reset_index(drop=True)
hot_df.insert(0, "rank", hot_df.index + 1)
hot_df.to_csv(TAB / "C5_hotspots_ranked.tsv", sep="\t", index=False)

# Merge contiguous hotspot bins (gap <= BIN) into discrete LOCI so that reviewers
# see distinct neighbourhoods rather than many adjacent 5-kb windows.
hs = hot_df.sort_values("start_bp").reset_index(drop=True)
hs["grp"] = (hs["start_bp"].diff() > BIN).cumsum()
loci = hs.groupby("grp").agg(
    start_bp=("start_bp", "min"), end_bp=("end_bp", "max"),
    n_bins=("rank", "size"), max_FIRE=("FIRE_M", "max"),
    sum_methyl_T1=("methyl_n_T1", "sum"),
    in_compartmentA=("in_compartmentA", "max"),
    nearest_BGC=("nearest_BGC", "first"),
    dist_to_BGC_bp=("dist_to_BGC_bp", "min")).reset_index(drop=True)
loci["width_kb"] = (loci["end_bp"] - loci["start_bp"]) / 1000.0
# nearest gene + product at locus centre
lc_rows = []
for _, r in loci.iterrows():
    gid, prod, _ = nearest_gene((r["start_bp"] + r["end_bp"]) / 2.0)
    lc_rows.append((gid, prod))
loci["nearest_gene"] = [x[0] for x in lc_rows]
loci["nearest_product"] = [x[1] for x in lc_rows]
loci = loci.sort_values(["sum_methyl_T1", "max_FIRE"],
                        ascending=False).reset_index(drop=True)
loci.insert(0, "locus_rank", loci.index + 1)
loci.to_csv(TAB / "C5_hotspot_loci_merged.tsv", sep="\t", index=False)


# ----------------------------------------------------------------------------
# C6 — CHARACTERISE the hotspot set vs background
# ----------------------------------------------------------------------------
def classify(product: str) -> str:
    """Coarse functional bucket from product string."""
    if not isinstance(product, str):
        return "other"
    p = product.lower()
    if "hypothetical" in p:
        return "hypothetical"
    if any(k in p for k in ["transcriptional regulator", "transcription factor",
                            "sigma factor", "anti-sigma", "two-component",
                            "response regulator", "histidine kinase", "regulatory"]):
        return "regulation/signaling"
    if any(k in p for k in ["transposase", "integrase", "recombinase", "insertion",
                            "mobile", "phage", "transposon"]):
        return "mobile/MGE"
    if any(k in p for k in ["polyketide", "nonribosomal", "nrps", "pks",
                            "synthase", "synthetase", "biosynthesis",
                            "secondary metab", "tailoring"]):
        return "secondary metabolism/biosynthesis"
    if any(k in p for k in ["abc transporter", "transporter", "permease",
                            "efflux", "mfs"]):
        return "transport"
    if any(k in p for k in ["restriction", "methyltransferase", "modification",
                            "crispr", "toxin-antitoxin", "antitoxin",
                            "defense", "abortive"]):
        return "defense/RM/modification"
    if any(k in p for k in ["oxidoreductase", "dehydrogenase", "hydrolase",
                            "transferase", "kinase", "reductase", "oxidase",
                            "isomerase", "lyase", "ligase", "esterase",
                            "peptidase", "protease"]):
        return "enzyme/metabolism"
    return "other"


def characterise(df_bins: pd.DataFrame, label: str) -> dict:
    span_bp = len(df_bins) * BIN
    genes = pd.concat([genes_in_window(int(r.start_bp), int(r.end_bp))
                       for r in df_bins.itertuples()]) if len(df_bins) else ann.iloc[0:0]
    genes = genes.drop_duplicates("gene_id")
    n_genes = len(genes)
    gene_density = n_genes / (span_bp / 1000.0) if span_bp else 0.0  # genes per kb
    cats = genes["product"].apply(classify).value_counts(normalize=True)
    return {
        "set": label,
        "n_windows": len(df_bins),
        "span_Mb": round(span_bp / 1e6, 3),
        "n_genes": n_genes,
        "genes_per_kb": round(gene_density, 4),
        "frac_compA": round(df_bins["in_compA"].mean(), 3) if len(df_bins) else 0,
        "frac_core": round(df_bins["in_core"].mean(), 3) if len(df_bins) else 0,
        "median_FIRE": round(df_bins["FIRE_M"].median(), 3) if len(df_bins) else 0,
        "median_methyl_n": round(df_bins["m_n_T1"].median(), 2) if len(df_bins) else 0,
        "cat_dist": cats.round(3).to_dict(),
    }


hot_bins = bins[bins["hotspot"] == 1]
bg_bins = bins[bins["hotspot"] == 0]
char_hot = characterise(hot_bins, "hotspot")
char_bg = characterise(bg_bins, "background")
char_all = characterise(bins, "genome_all")

# enrichment of functional categories: hotspot genes vs all chromosomal genes
hot_genes = pd.concat([genes_in_window(int(r.start_bp), int(r.end_bp))
                       for r in hot_bins.itertuples()]).drop_duplicates("gene_id")
hot_genes = hot_genes.copy()
hot_genes["cat"] = hot_genes["product"].apply(classify)
ann2 = ann.copy()
ann2["cat"] = ann2["product"].apply(classify)

from scipy.stats import fisher_exact
enr_rows = []
all_cats = sorted(set(ann2["cat"]))
n_hot = len(hot_genes)
n_all = len(ann2)
for cat in all_cats:
    a = int((hot_genes["cat"] == cat).sum())     # hotspot in cat
    b = n_hot - a                                # hotspot not cat
    c = int((ann2["cat"] == cat).sum()) - a      # bg in cat
    d = (n_all - n_hot) - c                      # bg not cat
    if a + c == 0:
        continue
    orr, p = fisher_exact([[a, b], [c, d]])
    enr_rows.append({
        "category": cat,
        "hotspot_genes": a,
        "hotspot_frac": round(a / n_hot, 3) if n_hot else 0,
        "genome_genes": int((ann2["cat"] == cat).sum()),
        "genome_frac": round(int((ann2["cat"] == cat).sum()) / n_all, 3),
        "odds_ratio": round(orr, 3),
        "fisher_p": p,
    })
enr_df = pd.DataFrame(enr_rows).sort_values("odds_ratio", ascending=False)
enr_df.to_csv(TAB / "C6_functional_enrichment.tsv", sep="\t", index=False)

# write characterisation summary
char_df = pd.DataFrame([
    {k: v for k, v in char_hot.items() if k != "cat_dist"},
    {k: v for k, v in char_bg.items() if k != "cat_dist"},
    {k: v for k, v in char_all.items() if k != "cat_dist"},
])
char_df.to_csv(TAB / "C6_set_characterisation.tsv", sep="\t", index=False)


# ----------------------------------------------------------------------------
# C7 — RESOLUTION sweep: window size vs candidate counts vs concordance
# ----------------------------------------------------------------------------
raw = pd.read_csv(BASE / "37_defense_island_GCCGGC/tables/"
                  "GCCGGC_sites_by_timepoint.tsv", sep="\t")
raw = raw[(raw["chrom"] == CHROM) & (raw["timepoint"] == "T1")].copy()
raw["position"] = pd.to_numeric(raw["position"], errors="coerce")
raw = raw.dropna(subset=["position"])
t1_pos = raw["position"].to_numpy()

# FIRE is only native at 5 kb. For window sizes that are multiples of 5 kb we can
# aggregate FIRE (mean of constituent 5-kb bins) and methylation (site counts);
# for the sub-5-kb case we report methylation-only resolution.
fire_by_bin = bins.set_index("bin")["FIRE_M"]


def sweep(win: int, qf: float = QF, qm: float = QM) -> dict:
    """Tile chromosome into `win`-bp windows; count candidates; FIRE concordance."""
    genome_end = 8_700_000
    edges = np.arange(0, genome_end + win, win)
    idx = np.arange(len(edges) - 1)
    # methylation per window
    m_count = np.histogram(t1_pos, bins=edges)[0]
    # FIRE per window: mean over overlapping 5-kb bins (win multiple of 5kb)
    fire_win = np.full(len(idx), np.nan)
    if win % BIN == 0:
        k = win // BIN
        fm = bins.set_index("bin")["FIRE_M"].reindex(
            range(0, int(genome_end // BIN) + 1)).to_numpy()
        # mean of k consecutive 5-kb bins
        for j in idx:
            seg = fm[j * k:(j + 1) * k]
            seg = seg[~np.isnan(seg)]
            if len(seg):
                fire_win[j] = seg.mean()
    valid = ~np.isnan(fire_win) if win % BIN == 0 else np.ones(len(idx), bool)

    # restrict concordance computation to windows with defined FIRE
    m_v = m_count[valid]
    f_v = fire_win[valid]
    methyl_pos_v = m_v[m_v > 0]
    if len(methyl_pos_v) == 0:
        return None
    m_thr = np.quantile(methyl_pos_v, qm)
    hi_m = m_v >= m_thr
    out = {"window_bp": win,
           "n_windows_total": int(valid.sum()),
           "n_methyl_windows": int((m_v > 0).sum()),
           "n_top_methyl_windows": int(hi_m.sum())}
    if win % BIN == 0 and not np.all(np.isnan(f_v)):
        f_thr = np.nanquantile(f_v, qf)
        hi_f = f_v >= f_thr
        tp = int(np.sum(hi_m & hi_f))
        n_hot_w = int(np.sum(hi_m & hi_f))
        prec = tp / hi_m.sum() if hi_m.sum() else 0      # of top-methyl, frac also top-FIRE
        rec = tp / hi_f.sum() if hi_f.sum() else 0       # of top-FIRE, frac also top-methyl
        union = int(np.sum(hi_m | hi_f))
        jac = tp / union if union else 0
        # hypergeometric enrichment of overlap
        from scipy.stats import hypergeom
        N = len(f_v)
        K = int(hi_f.sum())
        n = int(hi_m.sum())
        p_over = hypergeom.sf(tp - 1, N, K, n)
        out.update({"n_hotspot_windows": n_hot_w,
                    "precision_methyl->FIRE": round(prec, 3),
                    "recall_FIRE_captured": round(rec, 3),
                    "jaccard": round(jac, 3),
                    "overlap_p_hypergeom": p_over})
    else:
        out.update({"n_hotspot_windows": np.nan,
                    "precision_methyl->FIRE": np.nan,
                    "recall_FIRE_captured": np.nan,
                    "jaccard": np.nan,
                    "overlap_p_hypergeom": np.nan})
    return out


win_sizes = [2000, 5000, 10000, 20000, 50000, 100000]
sweep_rows = [r for r in (sweep(w) for w in win_sizes) if r is not None]
sweep_df = pd.DataFrame(sweep_rows)
sweep_df.to_csv(TAB / "C7_resolution_sweep.tsv", sep="\t", index=False)


# ----------------------------------------------------------------------------
# Console report
# ----------------------------------------------------------------------------
print("=" * 74)
print("C5 — HOTSPOTS (top-%d%% FIRE  AND  top-%d%% methylation among methylated)"
      % (int((1 - QF) * 100), int((1 - QM) * 100)))
print("=" * 74)
print(f"FIRE_M threshold (q{QF}) = {fire_thr:.3f}")
print(f"methyl_n_T1 threshold (q{QM} of methylated bins) = {methyl_thr:.1f} sites/5kb")
print(f"hi-FIRE bins        : {bins['hi_fire'].sum()}")
print(f"hi-methylation bins : {bins['hi_methyl'].sum()}")
print(f"HOTSPOT bins (both) : {bins['hotspot'].sum()}  "
      f"({bins['hotspot'].mean()*100:.1f}% of {len(bins)} FIRE bins)")
print(f"  in Compartment A  : {hot_bins['in_compA'].sum()}/{len(hot_bins)} "
      f"({hot_bins['in_compA'].mean()*100:.0f}%)")
print(f"  in chromosomal core: {hot_bins['in_core'].sum()}/{len(hot_bins)} "
      f"({hot_bins['in_core'].mean()*100:.0f}%)")
print(f"\nMerged into {len(loci)} discrete loci "
      f"({(loci['n_bins']>1).sum()} multi-bin, {(loci['n_bins']==1).sum()} singleton)")
print("Top 10 merged loci (by methylation load):")
lcols = ["locus_rank", "start_bp", "end_bp", "width_kb", "n_bins", "max_FIRE",
         "sum_methyl_T1", "in_compartmentA", "dist_to_BGC_bp", "nearest_BGC",
         "nearest_product"]
with pd.option_context("display.width", 220, "display.max_colwidth", 40):
    print(loci[lcols].head(10).to_string(index=False))
print("\nTop 10 hotspot 5-kb bins (by combined FIRE+methyl percentile):")
cols = ["rank", "start_bp", "end_bp", "FIRE_M", "methyl_n_T1",
        "n_genes_in_window", "nearest_gene", "in_compartmentA",
        "dist_to_BGC_bp", "nearest_BGC", "nearest_product"]
with pd.option_context("display.width", 220, "display.max_colwidth", 40):
    print(hot_df[cols].head(10).to_string(index=False))

print("\n" + "=" * 74)
print("C6 — CHARACTERISATION")
print("=" * 74)
print(char_df.to_string(index=False))
print("\nFunctional category mix (hotspot genes):")
for k, v in sorted(char_hot["cat_dist"].items(), key=lambda x: -x[1]):
    print(f"  {k:38s} {v*100:5.1f}%")
print("\nEnrichment vs genome (Fisher, OR>1 = enriched in hotspots):")
with pd.option_context("display.width", 200):
    print(enr_df.to_string(index=False))

print("\n" + "=" * 74)
print("C7 — RESOLUTION SWEEP")
print("=" * 74)
with pd.option_context("display.width", 220):
    print(sweep_df.to_string(index=False))

print("\n[saved]")
for f in ["C5_hotspots_ranked.tsv", "C6_set_characterisation.tsv",
          "C6_functional_enrichment.tsv", "C7_resolution_sweep.tsv"]:
    print("  tables/" + f)
