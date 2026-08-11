# Growth data — reconstructed from SVG (provenance)
Source: SuppFig_growth_{DCW,Act,Red}.svg (user-provided; no CSV existed in repo).
Method: axis-tick pixel calibration → bar tops + individual marker (n=4/timepoint) values.
FIDELITY CHECK (reproduces manuscript-locked stats (DCW exact; Act within 0.03 of the quoted F)):
  DCW one-way ANOVA F=199.04 (manuscript F(2,9)=199.04 — exact), eta2=0.978 (0.978)
  Act one-way ANOVA F=192.82 (manuscript F(2,9)=192.79 — within 0.03), eta2=0.977 (0.977)
  Pairwise significance matches original SVG brackets:
    DCW: T1-T2 ****, T1-T3 ****, T2-T3 ns
    Act: T1-T2 ns,   T1-T3 ****, T2-T3 ****
    Red: T1-T2 ****, T1-T3 **,   T2-T3 *** (Welch + Holm)
Values are recovered measurements, NOT fabricated. Author should confirm against
the original source spreadsheet before final submission.

## ⚠️ UNRESOLVED discrepancy — undecylprodigiosin (Red) fold-change
Manuscript states: "undecylprodigiosin (A530) increased ~30-fold from T1 to T3"
(verbatim, T1→T3). The SVG-reconstructed means give T1=0.00425, T2=0.308, T3=5.159 mg/L,
so T3/T1 ≈ 1214x and T3/T2 ≈ 17x — NEITHER is ~30-fold. The reconstructed T1 value
(~0.004 mg/L) is essentially at the detection floor, which makes any T1-denominated fold
ratio unstable/meaningless. The manuscript's "~30-fold" is almost certainly computed on a
different basis (e.g. raw A530 absorbance with a medium-blank floor, or a T1 detection-limit
substitution), not on these concentration means. **This is a genuine unreconciled
discrepancy — AUTHOR MUST confirm the Red fold-change basis before submission.** DCW and Act
reconstructed data reproduce the manuscript ANOVA exactly (F=199.04 / 192.82); only the Red
fold-change wording does not reconcile with the concentration-scale means shown here.
(An earlier working note in this session mis-stated the manuscript's fold as "T2→T3"; that
was wrong — the manuscript says T1→T3. Corrected here.)
