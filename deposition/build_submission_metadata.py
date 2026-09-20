#!/usr/bin/env python3
"""Emit the filled SRA/BioSample submission tables.

Author-supplied facts (2026-09-10 Q&A) are the CONSTANTS below; everything
else is measured from the files. Re-run after editing a constant.
"""
import csv, glob, os, subprocess

METHYL = "/Users/okaban/bioinfo/methyl/260102_M145/analysis"
RNARAW = ("/Users/okaban/Library/CloudStorage/Dropbox-SFC-CNS/Takeda Tomoki/"
          "my_projects/2026/M145_RNA-seq/data/RH250715199/01_RawData")
NANO   = "/Users/okaban/Library/CloudStorage/Dropbox-SFC-CNS/Takeda Tomoki/nanopore"
PROV   = "deposition/nanopore_run_provenance.tsv"
OUT    = "deposition"

# --- author-supplied (Q&A 2026-09-10) --------------------------------------
ORGANISM   = "Streptomyces coelicolor A3(2)"
STRAIN     = "M145"
ISOLATION  = "R5 liquid culture"
COLLECTED  = "2025-12"          # time-course sampling 2025-12-04..07
GEO        = "Japan"
LAB_HOST   = "not applicable"
# --- vendor report RH250715199 ---------------------------------------------
ILLUMINA_MODEL = "Illumina NovaSeq X Plus"
RNA_DESIGN = ("Total RNA; rRNA depletion with Illumina Ribo-Zero Plus; "
              "strand-specific library with NEBNext Ultra II Directional RNA "
              "Library Prep Kit; 150 bp paired-end")
# --- MinKNOW sample sheets / final_summary ---------------------------------
ONT_MODEL  = "MinION"
ONT_DESIGN = ("Native genomic DNA, no PCR; SQK-RBK114-96 rapid barcoding; "
              "FLO-MIN114 (R10.4.1) flow cell, 400 bps; basecalled with "
              "modification calling (6mA, 4mC, 5mC), aligned to NC_003888.3 "
              "with minimap2 -ax map-ont -y (MM/ML tags retained)")

TP = {"1": "12", "2": "24", "3": "50"}
SAMPLES = ["1-1","1-2","1-3","2-1","2-3","2-4","3-2","3-3","3-4"]

def size(p):
    try: return os.path.getsize(p)
    except OSError: return 0

def dsize(p):
    return sum(size(os.path.join(dp,f)) for dp,_,fn in os.walk(p) for f in fn)

# ---------------- BioSample (9 rows, one per biological sample) ------------
bios = []
for s in SAMPLES:
    bios.append(dict(
        sample_name=f"M145_{s.replace('-','_')}",
        sample_title=f"S. coelicolor M145, {TP[s[0]]} h, biological replicate {s.split('-')[1]}",
        organism=ORGANISM, strain=STRAIN, isolate="not applicable",
        collection_date=COLLECTED, geo_loc_name=GEO, isolation_source=ISOLATION,
        host=LAB_HOST, lat_lon="not applicable",
        culture_collection="not applicable",
        description=(f"Time point T{s[0]} ({TP[s[0]]} h after inoculation). "
                     "Genomic DNA and total RNA were extracted from the same "
                     "harvested culture, so the nanopore methylome and the "
                     "Illumina transcriptome for this sample derive from one "
                     "biological specimen.")))

# ---------------- SRA runs (18: 9 ONT BAM + 9 Illumina fastq pairs) --------
runs = []
for s in SAMPLES:
    bam = f"{METHYL}/{s}_mapped.bam"
    runs.append(dict(
        sample_name=f"M145_{s.replace('-','_')}",
        library_ID=f"M145_{s.replace('-','_')}_ONT",
        title=f"Nanopore methylome of S. coelicolor M145, {TP[s[0]]} h",
        library_strategy="OTHER", library_source="GENOMIC",
        library_selection="RANDOM", library_layout="single",
        platform="OXFORD_NANOPORE", instrument_model=ONT_MODEL,
        design_description=ONT_DESIGN, filetype="bam",
        filename=os.path.basename(bam), filename2="",
        file_size_gb=round(size(bam)/1073741824, 3)))
for s in SAMPLES:
    u = s.replace('-','_')
    f1, f2 = f"{RNARAW}/M145_{u}_1.fastq.gz", f"{RNARAW}/M145_{u}_2.fastq.gz"
    runs.append(dict(
        sample_name=f"M145_{u}",
        library_ID=f"M145_{u}_RNA",
        title=f"Strand-specific RNA-seq of S. coelicolor M145, {TP[s[0]]} h",
        library_strategy="RNA-Seq", library_source="TRANSCRIPTOMIC",
        library_selection="Inverse rRNA", library_layout="paired",
        platform="ILLUMINA", instrument_model=ILLUMINA_MODEL,
        design_description=RNA_DESIGN, filetype="fastq",
        filename=os.path.basename(f1), filename2=os.path.basename(f2),
        file_size_gb=round((size(f1)+size(f2))/1073741824, 3)))

# ---------------- pod5 upload manifest -------------------------------------
runroot = {}
for d in os.listdir(NANO):
    p = os.path.join(NANO, d, "no_sample_id")
    if os.path.isdir(p):
        for r in os.listdir(p):
            runroot[r] = os.path.join(p, r)
pod5 = []
for r in csv.DictReader(open(PROV), delimiter="\t"):
    d = os.path.join(runroot.get(r["run"], ""), "pod5_pass", r["barcode"])
    pod5.append(dict(sample_name=f"M145_{r['sample'].replace('-','_')}",
                     run=r["run"], barcode=r["barcode"],
                     pct_of_sample_reads=r["pct_of_sampled"],
                     path=d, size_gb=round(dsize(d)/1073741824, 2) if os.path.isdir(d) else 0))

for name, rows in (("biosample_attributes.tsv", bios),
                   ("sra_metadata.tsv", runs),
                   ("pod5_upload_manifest.tsv", pod5)):
    with open(f"{OUT}/{name}", "w", newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)
    print(f"{name}: {len(rows)} rows")

print(f"\nBioSamples {len(bios)} | SRA runs {len(runs)} "
      f"(ONT {sum(1 for r in runs if r['platform']=='OXFORD_NANOPORE')}, "
      f"Illumina {sum(1 for r in runs if r['platform']=='ILLUMINA')})")
print(f"run file volume: {sum(r['file_size_gb'] for r in runs):.2f} GB")
print(f"pod5 volume:     {sum(p['size_gb'] for p in pod5):.2f} GB "
      f"across {len(pod5)} dirs")
blank = [f"{r['sample_name']}:{k}" for r in runs+bios for k,v in r.items() if v==""
         and k!="filename2"]
print("empty required fields:", blank or "none")
