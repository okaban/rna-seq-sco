#!/usr/bin/env python3
"""
Figure 3 — Promoter protection architecture (single panel).

Promotes the former Figure 2b (TSS metagene of GCCGGC 4mC density, all annotated
genes) to a standalone main figure, per the 2026-08-11 figure review. This is the
first evidence for the promoter protection zone and stands on its own.
"""
import importlib.util, sys
from pathlib import Path
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent

def _load(name, fname):
    spec = importlib.util.spec_from_file_location(name, HERE / fname)
    mod = importlib.util.module_from_spec(spec); sys.modules[name] = mod
    spec.loader.exec_module(mod); return mod

f2 = _load("mod_f2b", "02b_figure2_RM_redistribution.py")
f6 = _load("mod_f6b", "30_figure6_spatial_organizer.py")

def main():
    f6.apply_style()
    print("=== Figure 3: promoter protection metagene (single panel) ===")
    df_spatial = f2.load_spatial_profile()

    fig = plt.figure(figsize=(f6.mm_to_inch(88), f6.mm_to_inch(72)))
    ax = fig.add_axes([0.235, 0.165, 0.725, 0.745])
    f2.panel_b(ax, df_spatial)
    ax.set_title("")

    out = f6.FIG_DIR / "Figure3_protection_metagene"
    f6.assert_no_text_collisions(fig, "Figure3_protection")
    # single-column class: the default "full" left the figure at 88 mm, 2 mm over
    f6.save_figure(fig, out, formats=("pdf", "svg", "png"), width_class="single")
    print(f"  written {out}.png")

if __name__ == "__main__":
    main()
