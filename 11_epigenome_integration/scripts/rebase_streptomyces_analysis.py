#!/usr/bin/env python3
"""
REBASE-based Streptomyces Comparative Methylome Analysis
=========================================================
Replace simulated data with real REBASE R-M system data for 37 Streptomyces strains.
Also compute motif site density from NCBI RefSeq genomes.

Author: Claude
Date: 2026-02-04
"""

import os
import re
import json
import subprocess
import difflib
from pathlib import Path
from collections import defaultdict
from multiprocessing import Pool, cpu_count

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

import warnings
warnings.filterwarnings("ignore")

# ──────────────────────────────────────────────
# Paths
# ──────────────────────────────────────────────
BASE_DIR = Path("/Users/okaban/bioinfo/rna-seq/11_epigenome_integration")
REBASE_DIR = BASE_DIR / "data" / "rebase"
NCBI_DIR = BASE_DIR / "data" / "ncbi_genomes"
OUTPUT_DIR = BASE_DIR / "analysis" / "20_atcc_real_methylome"
BAIROCH_PATH = REBASE_DIR / "bairoch.txt"

M145_GENOME = Path(
    "/Users/okaban/SFC-CNS Dropbox/Takeda Tomoki/my_projects/2026/"
    "M145_RNA-seq/data/RH250715199/02_M145-ref/ncbi_dataset/data/"
    "GCF_000203835.1/GCF_000203835.1_ASM20383v1_genomic.fna"
)

for d in [REBASE_DIR, NCBI_DIR, OUTPUT_DIR]:
    d.mkdir(parents=True, exist_ok=True)

# ──────────────────────────────────────────────
# 37 ATCC strains (from atcc_comparative_methylome.py)
# ──────────────────────────────────────────────
STREPTOMYCES_STRAINS = [
    {"name": "Streptomyces albus", "atcc_id": "ATCC 3004"},
    {"name": "Streptomyces ambofaciens", "atcc_id": "ATCC 15154"},
    {"name": "Streptomyces antibioticus-oligomycini", "atcc_id": "ATCC 11891"},
    {"name": "Streptomyces aureofaciens", "atcc_id": "ATCC 10762"},
    {"name": "Streptomyces aureoverticillatus", "atcc_id": "ATCC 15853"},
    {"name": "Streptomyces avermitilis", "atcc_id": "ATCC 31267"},
    {"name": "Streptomyces bikiniensis", "atcc_id": "ATCC 11062"},
    {"name": "Streptomyces californicus", "atcc_id": "ATCC 15436"},
    {"name": "Streptomyces canus", "atcc_id": "ATCC 12237"},
    {"name": "Streptomyces cattleya", "atcc_id": "ATCC 35852"},
    {"name": "Streptomyces clavuligerus", "atcc_id": "ATCC 27064"},
    {"name": "Streptomyces coelescens", "atcc_id": "ATCC 19833"},
    {"name": "Streptomyces diastatochromogenes", "atcc_id": "ATCC 12309"},
    {"name": "Streptomyces fulvissimus", "atcc_id": "ATCC 27431"},
    {"name": "Streptomyces globisporus subsp. globisporus", "atcc_id": "ATCC 15864"},
    {"name": "Streptomyces griseorubiginosus", "atcc_id": "ATCC 19766"},
    {"name": "Streptomyces griseus subsp. griseus", "atcc_id": "ATCC 10137"},
    {"name": "Streptomyces hygroscopicus subsp. hygroscopicus", "atcc_id": "ATCC 27438"},
    {"name": "Streptomyces lavendulae subsp. lavendulae", "atcc_id": "ATCC 8664"},
    {"name": "Streptomyces lividans", "atcc_id": "ATCC 19844"},
    {"name": "Streptomyces lydicus", "atcc_id": "ATCC 25470"},
    {"name": "Streptomyces natalensis", "atcc_id": "ATCC 27448"},
    {"name": "Streptomyces netropsis", "atcc_id": "ATCC 23940"},
    {"name": "Streptomyces nodosus", "atcc_id": "ATCC 14899"},
    {"name": "Streptomyces noursei", "atcc_id": "ATCC 11455"},
    {"name": "Streptomyces peucetius subsp. caesius", "atcc_id": "ATCC 27952"},
    {"name": "Streptomyces platensis", "atcc_id": "ATCC 23948"},
    {"name": "Streptomyces prasinus", "atcc_id": "ATCC 13879"},
    {"name": "Streptomyces rimosus subsp. rimosus", "atcc_id": "ATCC 10970"},
    {"name": "Streptomyces roseofulvus", "atcc_id": "ATCC 27454"},
    {"name": "Streptomyces sp. NRRL WC-3753", "atcc_id": "ATCC 55227"},
    {"name": "Streptomyces spectabilis", "atcc_id": "ATCC 27465"},
    {"name": "Streptomyces tsukubensis", "atcc_id": "ATCC 67382"},
    {"name": "Streptomyces venezuelae", "atcc_id": "ATCC 10595"},
    {"name": "Streptomyces violaceoruber", "atcc_id": "ATCC 14980"},
    {"name": "Streptomyces virginiae", "atcc_id": "ATCC 13161"},
    {"name": "Streptomyces viridochromogenes", "atcc_id": "ATCC 14920"},
]

# Target motifs from M145
TARGET_MOTIFS = {
    "CCGG": {"methyl_type": "m4C", "description": "McrBC-like, 4mC"},
    "AAGCCCG": {"methyl_type": "m6A", "description": "Novel R-M, 6mA"},
    "GATC": {"methyl_type": "m6A", "description": "Dam-like, 6mA"},
}

# ──────────────────────────────────────────────
# IUPAC degenerate base handling
# ──────────────────────────────────────────────
IUPAC = {
    "A": {"A"}, "C": {"C"}, "G": {"G"}, "T": {"T"},
    "R": {"A", "G"}, "Y": {"C", "T"}, "S": {"G", "C"},
    "W": {"A", "T"}, "K": {"G", "T"}, "M": {"A", "C"},
    "B": {"C", "G", "T"}, "D": {"A", "G", "T"},
    "H": {"A", "C", "T"}, "V": {"A", "C", "G"},
    "N": {"A", "C", "G", "T"},
}

COMPLEMENT = {"A": "T", "T": "A", "G": "C", "C": "G",
              "R": "Y", "Y": "R", "S": "S", "W": "W",
              "K": "M", "M": "K", "B": "V", "V": "B",
              "D": "H", "H": "D", "N": "N"}


def reverse_complement(seq: str) -> str:
    return "".join(COMPLEMENT.get(b, "N") for b in reversed(seq.upper()))


def iupac_match(pattern: str, target: str) -> bool:
    """Check if pattern matches target allowing IUPAC degeneracy."""
    if len(pattern) != len(target):
        return False
    for p, t in zip(pattern.upper(), target.upper()):
        p_set = IUPAC.get(p, {p})
        t_set = IUPAC.get(t, {t})
        if not p_set & t_set:
            return False
    return True


def specificity_score(pattern: str) -> float:
    """Calculate specificity of a sequence pattern (0-1).
    N=0, specific bases=1, 2-fold degeneracy=0.5, etc."""
    score = 0
    for b in pattern.upper():
        bases = IUPAC.get(b, {b})
        score += 1.0 / len(bases)
    return score / len(pattern) if pattern else 0


def motif_matches_recognition(motif: str, recognition: str) -> bool:
    """Check if motif matches recognition sequence.

    Rules to avoid false positives:
    1. Exact-length match: IUPAC-aware comparison of full sequence
    2. Substring match: only if the matching region has high specificity
       (>50% of positions are specific bases, not N).
       This prevents Type I bipartite sites like CCGANNNNNNCTAC from
       matching CCGG via the degenerate NNNN region.
    """
    rec_clean = re.sub(r"[^A-Za-z]", "", recognition).upper()
    motif_up = motif.upper()
    rc_motif = reverse_complement(motif_up)

    # Exact-length match (full recognition = motif)
    if len(rec_clean) == len(motif_up):
        if iupac_match(rec_clean, motif_up) or iupac_match(rec_clean, rc_motif):
            return True
        return False

    # Substring match with specificity check
    min_specificity = 0.6  # at least 60% of matched positions must be specific
    for i in range(len(rec_clean) - len(motif_up) + 1):
        sub = rec_clean[i:i + len(motif_up)]
        if iupac_match(sub, motif_up) or iupac_match(sub, rc_motif):
            if specificity_score(sub) >= min_specificity:
                return True
    return False


# ══════════════════════════════════════════════
# Phase 1: REBASE Bairoch format parsing
# ══════════════════════════════════════════════

def parse_bairoch(filepath: Path) -> list[dict]:
    """Parse REBASE Bairoch format into list of enzyme records."""
    records = []
    current: dict = {}

    with open(filepath, "r") as fh:
        for line in fh:
            if line.startswith("//"):
                if current:
                    records.append(current)
                current = {}
                continue
            if line.startswith("CC") or line.strip() == "":
                continue

            tag = line[:5].strip()
            value = line[5:].strip()

            if tag == "ID":
                current["enzyme_name"] = value
            elif tag == "ET":
                current["enzyme_type"] = value
            elif tag == "AC":
                current["accession"] = value.rstrip(";")
            elif tag == "OS":
                current["organism"] = value
            elif tag == "RS":
                current.setdefault("recognition_raw", "")
                current["recognition_raw"] += " " + value if current["recognition_raw"] else value
            elif tag == "MS":
                current.setdefault("methylation_raw", "")
                current["methylation_raw"] += " " + value if current["methylation_raw"] else value

    return records


def filter_streptomyces(records: list[dict]) -> pd.DataFrame:
    """Filter records for Streptomyces organisms."""
    strepto = [r for r in records if r.get("organism", "").startswith("Streptomyces")]
    df = pd.DataFrame(strepto)
    if df.empty:
        return df

    # Parse recognition sequences (remove cut-site notation)
    df["recognition_seq"] = df["recognition_raw"].apply(
        lambda x: re.sub(r"[^A-Za-z,; ]", "", str(x)).strip().split(",")[0].strip()
        if pd.notna(x) else ""
    )

    # Parse methylation type
    def parse_methyl(raw):
        if pd.isna(raw) or raw == "":
            return ""
        types = re.findall(r"\((m[456][A-Za-z]+)\)", str(raw))
        return ", ".join(types) if types else ""

    df["methyl_type"] = df["methylation_raw"].apply(parse_methyl)

    return df


def match_strains_to_rebase(rebase_df: pd.DataFrame) -> pd.DataFrame:
    """Match 37 ATCC strains to REBASE organisms using fuzzy matching."""
    unique_organisms = rebase_df["organism"].unique().tolist()

    # Build a species-level lookup
    def normalize_species(name: str) -> str:
        parts = name.split()
        if len(parts) >= 2:
            return f"{parts[0]} {parts[1]}".lower()
        return name.lower()

    rebase_species = {normalize_species(o): o for o in unique_organisms}

    results = []
    for strain in STREPTOMYCES_STRAINS:
        name = strain["name"]
        norm = normalize_species(name)

        # Direct species match
        matched_organisms = []
        for org in unique_organisms:
            org_norm = normalize_species(org)
            if norm == org_norm or norm in org_norm or org_norm in norm:
                matched_organisms.append(org)

        # Fuzzy match fallback
        if not matched_organisms:
            close = difflib.get_close_matches(norm, list(rebase_species.keys()), n=1, cutoff=0.7)
            if close:
                matched_organisms = [rebase_species[close[0]]]

        # Get R-M entries for matched organisms
        if matched_organisms:
            entries = rebase_df[rebase_df["organism"].isin(matched_organisms)]
            results.append({
                "strain_name": name,
                "atcc_id": strain["atcc_id"],
                "rebase_organisms": "; ".join(sorted(set(matched_organisms))),
                "n_rm_entries": len(entries),
                "recognition_seqs": "; ".join(sorted(entries["recognition_seq"].dropna().unique())),
                "enzyme_names": "; ".join(sorted(entries["enzyme_name"].dropna().unique())),
                "enzyme_types": "; ".join(sorted(entries["enzyme_type"].dropna().unique())),
                "methyl_types": "; ".join(sorted(entries["methyl_type"].dropna().unique())),
                "match_status": "matched",
            })
        else:
            results.append({
                "strain_name": name,
                "atcc_id": strain["atcc_id"],
                "rebase_organisms": "",
                "n_rm_entries": 0,
                "recognition_seqs": "",
                "enzyme_names": "",
                "enzyme_types": "",
                "methyl_types": "",
                "match_status": "no_data",
            })

    return pd.DataFrame(results)


def compute_motif_conservation(strain_match_df: pd.DataFrame,
                               rebase_df: pd.DataFrame) -> pd.DataFrame:
    """For each strain, check if CCGG/AAGCCCG/GATC are recognized by any R-M system."""
    rows = []
    for _, strain_row in strain_match_df.iterrows():
        row = {
            "strain_name": strain_row["strain_name"],
            "atcc_id": strain_row["atcc_id"],
            "match_status": strain_row["match_status"],
        }

        if strain_row["match_status"] == "no_data":
            for motif in TARGET_MOTIFS:
                row[f"{motif}_present"] = "no_data"
                row[f"{motif}_enzyme"] = ""
            rows.append(row)
            continue

        # Get all entries for this strain's organisms
        orgs = [o.strip() for o in strain_row["rebase_organisms"].split(";") if o.strip()]
        entries = rebase_df[rebase_df["organism"].isin(orgs)]

        for motif, info in TARGET_MOTIFS.items():
            matching_enzymes = []
            for _, entry in entries.iterrows():
                rec = str(entry.get("recognition_seq", ""))
                if rec and motif_matches_recognition(motif, rec):
                    matching_enzymes.append(entry["enzyme_name"])

            row[f"{motif}_present"] = "yes" if matching_enzymes else "no"
            row[f"{motif}_enzyme"] = "; ".join(matching_enzymes)

        rows.append(row)

    return pd.DataFrame(rows)


# ══════════════════════════════════════════════
# Phase 2: NCBI genome motif site density
# ══════════════════════════════════════════════

def load_genome(genome_path: Path) -> dict:
    """Load genome FASTA -> {chrom_id: sequence}."""
    sequences = {}
    current_id = None
    current_seq = []

    with open(genome_path, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_id:
                    sequences[current_id] = "".join(current_seq)
                current_id = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line.upper())
        if current_id:
            sequences[current_id] = "".join(current_seq)

    return sequences


def count_motif_sites(genome: dict, motif: str) -> int:
    """Count occurrences of motif on both strands."""
    motif_up = motif.upper()
    rc = reverse_complement(motif_up)
    count = 0
    for seq in genome.values():
        count += seq.count(motif_up)
        if rc != motif_up:
            count += seq.count(rc)
    return count


def genome_size_bp(genome: dict) -> int:
    return sum(len(s) for s in genome.values())


def compute_m145_density() -> dict:
    """Compute motif site density for M145 reference."""
    print("  Computing M145 motif site density...")
    genome = load_genome(M145_GENOME)
    size_mb = genome_size_bp(genome) / 1e6
    result = {"strain": "M145 (reference)", "genome_size_mb": round(size_mb, 2)}
    for motif in TARGET_MOTIFS:
        n = count_motif_sites(genome, motif)
        result[f"{motif}_count"] = n
        result[f"{motif}_density"] = round(n / size_mb, 1)
    return result


def download_ncbi_genomes() -> dict:
    """Download RefSeq genomes for 37 strains using NCBI datasets CLI.
    Returns dict of {species_name: fasta_path}.
    """
    downloaded = {}
    species_names = set()
    for s in STREPTOMYCES_STRAINS:
        # Normalize to base species for NCBI search
        parts = s["name"].split()
        if len(parts) >= 2:
            species = f"{parts[0]} {parts[1]}"
        else:
            species = s["name"]
        species_names.add((species, s["name"], s["atcc_id"]))

    for species, full_name, atcc_id in sorted(species_names):
        safe_name = species.replace(" ", "_")
        out_dir = NCBI_DIR / safe_name
        fasta_glob = list(out_dir.glob("*.fna")) if out_dir.exists() else []

        if fasta_glob:
            downloaded[full_name] = fasta_glob[0]
            continue

        out_dir.mkdir(parents=True, exist_ok=True)
        try:
            result = subprocess.run(
                [
                    "datasets", "download", "genome", "taxon", species,
                    "--reference", "--include", "genome",
                    "--filename", str(out_dir / "dataset.zip"),
                ],
                capture_output=True, text=True, timeout=120,
            )
            if result.returncode == 0 and (out_dir / "dataset.zip").exists():
                subprocess.run(
                    ["unzip", "-o", "-q", str(out_dir / "dataset.zip"), "-d", str(out_dir)],
                    capture_output=True, timeout=60,
                )
                # Find the FASTA
                fasta_files = list(out_dir.rglob("*.fna"))
                if fasta_files:
                    downloaded[full_name] = fasta_files[0]
                    print(f"    Downloaded: {species}")
                else:
                    print(f"    No FASTA found after unzip: {species}")
            else:
                print(f"    No RefSeq genome available: {species}")
        except (subprocess.TimeoutExpired, FileNotFoundError) as e:
            print(f"    Download failed for {species}: {e}")

    return downloaded


def compute_density_all(genome_paths: dict) -> pd.DataFrame:
    """Compute motif site density for all downloaded genomes."""
    rows = []
    for name, fpath in sorted(genome_paths.items()):
        try:
            genome = load_genome(fpath)
            size_mb = genome_size_bp(genome) / 1e6
            row = {"strain": name, "genome_size_mb": round(size_mb, 2)}
            for motif in TARGET_MOTIFS:
                n = count_motif_sites(genome, motif)
                row[f"{motif}_count"] = n
                row[f"{motif}_density"] = round(n / size_mb, 1)
            rows.append(row)
        except Exception as e:
            print(f"    Error processing {name}: {e}")
    return pd.DataFrame(rows)


# ══════════════════════════════════════════════
# Phase 3: Figure generation
# ══════════════════════════════════════════════

def save_fig(fig, name: str, out_dir: Path = OUTPUT_DIR):
    fig.savefig(out_dir / f"{name}.pdf", bbox_inches="tight")
    fig.savefig(out_dir / f"{name}.svg", bbox_inches="tight")
    plt.close(fig)


def fig_motif_conservation(conservation_df: pd.DataFrame):
    """Fig 1: REBASE R-M system conservation bar chart."""
    # Exclude no_data strains
    data_df = conservation_df[conservation_df["match_status"] != "no_data"]
    n_total = len(data_df)
    n_with_data = n_total  # strains with REBASE data

    motifs = list(TARGET_MOTIFS.keys())
    colors = ["#E74C3C", "#3498DB", "#2ECC71"]

    counts = []
    fractions = []
    ci_lo, ci_hi = [], []
    for motif in motifs:
        col = f"{motif}_present"
        n_yes = (data_df[col] == "yes").sum()
        frac = n_yes / n_with_data if n_with_data > 0 else 0
        counts.append(n_yes)
        fractions.append(frac * 100)

        # Binomial 95% CI (Wilson score)
        if n_with_data > 0:
            ci = stats.binom.interval(0.95, n_with_data, frac)
            ci_lo.append(frac * 100 - ci[0] / n_with_data * 100)
            ci_hi.append(ci[1] / n_with_data * 100 - frac * 100)
        else:
            ci_lo.append(0)
            ci_hi.append(0)

    fig, ax = plt.subplots(figsize=(8, 5))
    bars = ax.bar(motifs, fractions, color=colors, edgecolor="black", alpha=0.85)

    ax.errorbar(motifs, fractions, yerr=[ci_lo, ci_hi],
                fmt="none", ecolor="black", capsize=5, capthick=1.5)

    for bar, frac, cnt in zip(bars, fractions, counts):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + max(ci_hi) + 2,
                f"{frac:.1f}%\n({cnt}/{n_with_data})",
                ha="center", va="bottom", fontweight="bold", fontsize=10)

    ax.set_ylim(0, 110)
    ax.set_ylabel("Conservation (%)", fontsize=12)
    ax.set_xlabel("Methylation Motif", fontsize=12)
    ax.set_title(
        f"R-M System Conservation in Streptomyces (REBASE)\n"
        f"(n={n_with_data} strains with REBASE data, "
        f"{len(conservation_df) - n_with_data} without data)",
        fontweight="bold", fontsize=11,
    )
    ax.axhline(y=50, color="gray", linestyle="--", alpha=0.5)

    save_fig(fig, "motif_conservation_rebase")
    print("  Saved: motif_conservation_rebase.pdf/svg")


def fig_rm_heatmap(conservation_df: pd.DataFrame):
    """Fig 2: Strain x Motif heatmap (R-M presence)."""
    motifs = list(TARGET_MOTIFS.keys())
    data = conservation_df.copy()
    data = data.sort_values("strain_name")

    # Encode: yes=2, no=1, no_data=0
    heatmap_vals = []
    for _, row in data.iterrows():
        vals = []
        for m in motifs:
            v = row[f"{m}_present"]
            if v == "yes":
                vals.append(2)
            elif v == "no":
                vals.append(1)
            else:
                vals.append(0)
        heatmap_vals.append(vals)

    hm_df = pd.DataFrame(heatmap_vals, columns=motifs,
                          index=data["strain_name"].values)

    from matplotlib.colors import ListedColormap
    cmap = ListedColormap(["#CCCCCC", "#FFCCCC", "#CC3333"])

    fig, ax = plt.subplots(figsize=(8, 12))
    sns.heatmap(hm_df, ax=ax, cmap=cmap, vmin=0, vmax=2,
                linewidths=0.5, linecolor="white",
                cbar_kws={"ticks": [0, 1, 2], "label": ""})

    cbar = ax.collections[0].colorbar
    cbar.set_ticklabels(["No data", "Absent", "Present"])

    ax.set_title("R-M System Presence Across Streptomyces\n(REBASE)",
                 fontweight="bold", fontsize=11)
    ax.set_xlabel("Methylation Motif")
    ax.set_ylabel("")

    plt.tight_layout()
    save_fig(fig, "rm_system_heatmap")
    print("  Saved: rm_system_heatmap.pdf/svg")


def fig_site_density(density_df: pd.DataFrame, m145_density: dict):
    """Fig 3: Motif site density boxplot (M145 vs genus)."""
    motifs = list(TARGET_MOTIFS.keys())
    colors = ["#E74C3C", "#3498DB", "#2ECC71"]

    fig, axes = plt.subplots(1, 3, figsize=(12, 5))

    for ax, motif, color in zip(axes, motifs, colors):
        col = f"{motif}_density"
        values = density_df[col].dropna().values
        m145_val = m145_density.get(col, 0)

        if len(values) > 0:
            bp = ax.boxplot(values, patch_artist=True, widths=0.6)
            bp["boxes"][0].set_facecolor(color)
            bp["boxes"][0].set_alpha(0.5)

            ax.scatter([1], [m145_val], color="red", s=120, zorder=5,
                       marker="*", label=f"M145 ({m145_val:.1f})")

            if len(values) >= 3:
                percentile = (values < m145_val).sum() / len(values) * 100
                ax.annotate(f"M145: {percentile:.0f}th pctl",
                            xy=(1, m145_val), xytext=(1.25, m145_val),
                            fontsize=9, color="red")

                # Mann-Whitney U test (M145 vs genus median)
                _, p_val = stats.mannwhitneyu([m145_val], values, alternative="two-sided")
                ax.text(0.95, 0.02, f"n={len(values)}", transform=ax.transAxes,
                        ha="right", va="bottom", fontsize=8, color="gray")
        else:
            ax.text(0.5, 0.5, "No data", transform=ax.transAxes,
                    ha="center", va="center", fontsize=12, color="gray")

        ax.set_title(f"{motif}\n({TARGET_MOTIFS[motif]['description']})",
                     fontweight="bold", fontsize=10)
        ax.set_ylabel("Sites / Mb")
        ax.set_xticks([])
        ax.legend(loc="upper right", fontsize=8)

    plt.suptitle("Motif Site Density: M145 vs. Streptomyces Genus (NCBI RefSeq)",
                 fontsize=12, fontweight="bold", y=1.02)
    plt.tight_layout()
    save_fig(fig, "motif_site_density_boxplot")
    print("  Saved: motif_site_density_boxplot.pdf/svg")


def fig_simulated_vs_real(conservation_df: pd.DataFrame):
    """Fig 4: Compare old simulated vs new real conservation rates."""
    # Simulated values from old analysis (np.random.seed(42))
    simulated = {"CCGG": 86.5, "AAGCCCG": 21.6, "GATC": 67.6}

    data_df = conservation_df[conservation_df["match_status"] != "no_data"]
    n_with_data = len(data_df)

    real = {}
    for motif in TARGET_MOTIFS:
        n_yes = (data_df[f"{motif}_present"] == "yes").sum()
        real[motif] = n_yes / n_with_data * 100 if n_with_data > 0 else 0

    motifs = list(TARGET_MOTIFS.keys())
    x = np.arange(len(motifs))
    width = 0.35

    fig, ax = plt.subplots(figsize=(8, 5))
    bars1 = ax.bar(x - width / 2, [simulated[m] for m in motifs], width,
                   label="Simulated (seed=42)", color="#AAAAAA", edgecolor="black")
    bars2 = ax.bar(x + width / 2, [real[m] for m in motifs], width,
                   label=f"Real (REBASE, n={n_with_data})", color="#3498DB", edgecolor="black")

    for bar in bars1:
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 1,
                f"{bar.get_height():.1f}%", ha="center", va="bottom", fontsize=9, color="gray")
    for bar in bars2:
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 1,
                f"{bar.get_height():.1f}%", ha="center", va="bottom", fontsize=9, fontweight="bold")

    ax.set_xticks(x)
    ax.set_xticklabels(motifs)
    ax.set_ylabel("Conservation (%)")
    ax.set_ylim(0, 110)
    ax.set_title("Simulated vs. Real Motif Conservation Rates",
                 fontweight="bold", fontsize=11)
    ax.legend()
    ax.axhline(y=50, color="gray", linestyle="--", alpha=0.3)

    save_fig(fig, "simulated_vs_real_comparison")
    print("  Saved: simulated_vs_real_comparison.pdf/svg")


def fig_rm_type_distribution(rebase_strep_df: pd.DataFrame):
    """Fig 5: R-M system type distribution in Streptomyces."""
    type_map = {
        "R1": "Type I", "M1": "Type I", "S1": "Type I",
        "R2": "Type II", "M2": "Type II",
        "R3": "Type III", "M3": "Type III",
        "R4": "Type IV",
        "IE": "Homing", "IN": "Nicking",
    }

    df = rebase_strep_df.copy()
    df["rm_type"] = df["enzyme_type"].map(type_map).fillna("Other/Unknown")

    type_counts = df["rm_type"].value_counts()

    fig, ax = plt.subplots(figsize=(8, 5))
    colors_map = {
        "Type I": "#E74C3C", "Type II": "#3498DB", "Type III": "#2ECC71",
        "Type IV": "#F39C12", "Homing": "#9B59B6", "Nicking": "#1ABC9C",
        "Other/Unknown": "#95A5A6",
    }
    bar_colors = [colors_map.get(t, "#95A5A6") for t in type_counts.index]

    bars = ax.bar(type_counts.index, type_counts.values, color=bar_colors, edgecolor="black")
    for bar, val in zip(bars, type_counts.values):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.5,
                str(val), ha="center", va="bottom", fontweight="bold")

    ax.set_ylabel("Number of Entries")
    ax.set_xlabel("R-M System Type")
    ax.set_title(f"Distribution of R-M System Types in Streptomyces (REBASE)\n"
                 f"(n={len(df)} entries from {df['organism'].nunique()} organisms)",
                 fontweight="bold", fontsize=11)
    plt.xticks(rotation=30, ha="right")
    plt.tight_layout()
    save_fig(fig, "rm_type_distribution")
    print("  Saved: rm_type_distribution.pdf/svg")


# ══════════════════════════════════════════════
# Main pipeline
# ══════════════════════════════════════════════

def main():
    log_path = OUTPUT_DIR / "analysis_log.txt"
    log = open(log_path, "w")

    def log_print(msg):
        print(msg)
        log.write(msg + "\n")

    log_print("=" * 70)
    log_print("REBASE Comparative Methylome Analysis — Real Data")
    log_print("=" * 70)

    # ── Phase 1: REBASE ──
    log_print("\n[Phase 1] Parsing REBASE Bairoch data...")
    all_records = parse_bairoch(BAIROCH_PATH)
    log_print(f"  Total REBASE records: {len(all_records)}")

    rebase_strep = filter_streptomyces(all_records)
    log_print(f"  Streptomyces entries: {len(rebase_strep)}")
    log_print(f"  Unique organisms: {rebase_strep['organism'].nunique()}")

    rebase_strep.to_csv(REBASE_DIR / "streptomyces_rm_systems.csv", index=False)
    log_print(f"  Saved: {REBASE_DIR / 'streptomyces_rm_systems.csv'}")

    # Match strains
    log_print("\n[Phase 1b] Matching 37 ATCC strains to REBASE...")
    strain_match = match_strains_to_rebase(rebase_strep)
    matched = strain_match[strain_match["match_status"] == "matched"]
    log_print(f"  Matched: {len(matched)}/37 strains")
    log_print(f"  No data: {37 - len(matched)} strains")

    strain_match.to_csv(REBASE_DIR / "strain_rebase_match.csv", index=False)

    # Conservation
    log_print("\n[Phase 1c] Computing motif conservation...")
    conservation = compute_motif_conservation(strain_match, rebase_strep)
    conservation.to_csv(REBASE_DIR / "motif_conservation_matrix.csv", index=False)
    log_print(f"  Saved: {REBASE_DIR / 'motif_conservation_matrix.csv'}")

    data_df = conservation[conservation["match_status"] != "no_data"]
    n_wd = len(data_df)
    for motif in TARGET_MOTIFS:
        n_yes = (data_df[f"{motif}_present"] == "yes").sum()
        pct = n_yes / n_wd * 100 if n_wd > 0 else 0
        log_print(f"  {motif}: {n_yes}/{n_wd} ({pct:.1f}%)")

    # ── Phase 2: NCBI genomes ──
    log_print("\n[Phase 2] Downloading NCBI RefSeq genomes...")
    genome_paths = download_ncbi_genomes()
    log_print(f"  Genomes available: {len(genome_paths)}/{len(STREPTOMYCES_STRAINS)}")

    log_print("\n[Phase 2b] Computing M145 density...")
    m145_density = compute_m145_density()
    log_print(f"  M145 genome: {m145_density['genome_size_mb']} Mb")
    for motif in TARGET_MOTIFS:
        log_print(f"  M145 {motif}: {m145_density[f'{motif}_count']} sites "
                  f"({m145_density[f'{motif}_density']} sites/Mb)")

    density_df = pd.DataFrame()
    if genome_paths:
        log_print("\n[Phase 2c] Computing density for downloaded genomes...")
        density_df = compute_density_all(genome_paths)
        density_df.to_csv(NCBI_DIR / "motif_site_density.csv", index=False)
        log_print(f"  Processed: {len(density_df)} genomes")
    else:
        log_print("  No NCBI genomes downloaded. Skipping density analysis.")

    # ── Phase 3: Figures ──
    log_print("\n[Phase 3] Generating figures...")
    fig_motif_conservation(conservation)
    fig_rm_heatmap(conservation)
    if not density_df.empty:
        fig_site_density(density_df, m145_density)
    fig_simulated_vs_real(conservation)
    fig_rm_type_distribution(rebase_strep)

    # ── Summary ──
    log_print("\n" + "=" * 70)
    log_print("Analysis complete!")
    log_print(f"Output directory: {OUTPUT_DIR}")
    log_print("Files:")
    for f in sorted(OUTPUT_DIR.iterdir()):
        log_print(f"  {f.name} ({f.stat().st_size / 1024:.1f} KB)")
    log_print("=" * 70)

    log.close()


if __name__ == "__main__":
    main()
