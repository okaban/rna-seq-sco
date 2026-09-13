#!/usr/bin/env python3
"""
Merged Figure 2 — GCCGGC 4mC geographic redistribution.

Combines, per the 2026-08-11 figure review (7 -> 6 main figures):
  a  GCCGGC 4mC density tracks T1/T2/T3       (from 02b_figure2_RM_redistribution.panel_a)
  b  core->arm stacked bars                   (from 30_figure6_spatial_organizer.draw_panelA)
  c  T1 core-methylated -> T2 arm-redistributed schematic
                                              (from 30_figure6_spatial_organizer.draw_schematic)

Former Figure 6b (occupancy vs log2FC scatter) is dropped: it plots the same
reg[['occ2k','LFC_T2vsT1','region']] columns as main Figure 4b (code-verified).
"""
import importlib.util, sys
from pathlib import Path
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent

def _load(name, fname):
    spec = importlib.util.spec_from_file_location(name, HERE / fname)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod

f2 = _load("mod_f2", "02b_figure2_RM_redistribution.py")
f6 = _load("mod_f6", "30_figure6_spatial_organizer.py")

def main():
    f6.apply_style()
    print("=== Merged Figure 2: geography (Fig2a + Fig6a + schematic) ===")

    _, df_sites = f2.load_geographic_redistribution()
    g, rowsA = f6._load_panelA()
    print("  panel b: " + ", ".join(
        f"{t} n={rowsA[t]['n']} core={rowsA[t]['pct_core']:.0f}%" for t in ("T1", "T2", "T3")))

    fig = plt.figure(figsize=(f6.mm_to_inch(174), f6.mm_to_inch(112)))
    ax_a = fig.add_axes([0.085, 0.630, 0.885, 0.300])
    # 2026-09-13 (FIG-08): panel b lifted 0.02 so its Core/Arms legend fits
    # BELOW the two-line x-tick labels instead of on top of them.
    ax_b = fig.add_axes([0.085, 0.175, 0.335, 0.310])
    # 2026-09-13 (FIG-08): panel c's axes top used to sit ABOVE panel b's top
    # (0.505 vs 0.485), so its panel letter landed on panel a's x-axis label
    # ("Chromosome position (Mb)"). Align the tops of b and c exactly.
    ax_s = fig.add_axes([0.525, 0.115, 0.465, 0.370])

    f2.panel_a(ax_a, df_sites)
    f6.draw_panelA(ax_b, rowsA)
    # 2026-09-13 (FIG-08): draw_panelA anchors its Core/Arms legend below the
    # axis at y=-0.30, where it overlapped the two-line x-tick labels
    # "(24 h)"/"(50 h)". Re-anchor it clear below the tick labels.
    ax_b.legend(frameon=False, fontsize=6, loc='upper center',
                bbox_to_anchor=(0.5, -0.30), ncol=2, handlelength=1.1,
                columnspacing=1.2, borderpad=0.2)
    f6.draw_schematic(ax_s)
    # draw_schematic reserves top whitespace for its in-plot mark legend; the
    # merged layout is tighter, so crop the unused band below the chromosomes.
    ax_s.set_ylim(2.6, 10.0)

    # panel titles collide with the neighbouring panel letters in this 3-panel
    # layout; the information is carried by the axis labels and the legend.
    ax_a.set_title("")
    ax_b.set_title("")

    for ax, L in ((ax_a, "a"), (ax_b, "b"), (ax_s, "c")):
        f6._utils.add_panel_label(ax, L)

    out = f6.FIG_DIR / "Figure2_merged_geography"
    f6.save_figure(fig, out, formats=("pdf", "svg", "png"))
    print(f"  written {out}.png")
    return out

if __name__ == "__main__":
    main()
