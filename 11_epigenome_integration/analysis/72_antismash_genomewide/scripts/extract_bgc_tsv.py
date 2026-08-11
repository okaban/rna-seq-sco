#!/usr/bin/env python3
"""Extract BGC/region summary from antiSMASH 8 JSON to TSV."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


def extract_regions(json_path: Path) -> list[dict]:
    with json_path.open() as f:
        data = json.load(f)

    rows: list[dict] = []
    for record in data.get("records", []):
        record_id = record.get("id", "?")
        features = record.get("features", [])

        regions = [f for f in features if f.get("type") == "region"]
        cdss = [f for f in features if f.get("type") == "CDS"]

        for region in regions:
            quals = region.get("qualifiers", {})
            loc = region.get("location", "")
            start, end = _parse_location(loc)

            region_num_list = quals.get("region_number", ["?"])
            region_num = region_num_list[0] if region_num_list else "?"
            products = quals.get("product", [])
            categories = quals.get("category", [])
            candidate_ids = quals.get("candidate_cluster_numbers", [])
            contig_edge = quals.get("contig_edge", ["False"])[0]

            n_genes = sum(1 for c in cdss if _overlaps(c.get("location", ""), start, end))

            rows.append({
                "region": f"{record_id}_region{region_num}",
                "record": record_id,
                "region_number": region_num,
                "start": start,
                "end": end,
                "length_bp": end - start,
                "products": ",".join(products),
                "categories": ",".join(categories),
                "candidate_clusters": ",".join(candidate_ids),
                "contig_edge": contig_edge,
                "n_genes": n_genes,
            })
    return rows


def _parse_location(loc: str) -> tuple[int, int]:
    """Parse antiSMASH/biopython location strings like '[123:456](+)' or 'join{...}'."""
    import re
    nums = re.findall(r"\d+", loc)
    if len(nums) >= 2:
        return int(nums[0]), int(nums[-1])
    return 0, 0


def _overlaps(loc: str, start: int, end: int) -> bool:
    a, b = _parse_location(loc)
    return not (b < start or a > end)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("json_path", type=Path)
    ap.add_argument("--out-tsv", type=Path, required=True)
    args = ap.parse_args()

    rows = extract_regions(args.json_path)
    cols = [
        "region", "record", "region_number",
        "start", "end", "length_bp",
        "products", "categories", "candidate_clusters",
        "contig_edge", "n_genes",
    ]
    with args.out_tsv.open("w") as f:
        f.write("\t".join(cols) + "\n")
        for r in rows:
            f.write("\t".join(str(r[c]) for c in cols) + "\n")
    print(f"Wrote {len(rows)} regions -> {args.out_tsv}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
