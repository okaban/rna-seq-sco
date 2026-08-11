# Deng et al. 2023 integration — comprehensive angle gap analysis (2026-06-29)

Deng 2023 (PNAS) local data = sd01–sd05 + our derived fire_methyl_bins (has FIRE_M, FIRE_L, PC_M, PC_L; meth T1/T2).

## Deng data dimensions available locally
- D1 FIRE per window, **M-phase AND L-phase** (sd04)
- D2 Compartment score PC, **M and L** (PC_M/PC_L in bins)
- D3 **CID (interaction-domain) boundaries + boundary-located genes, M and L** (sd02)  ← never used
- D4 TF motifs enriched at CID boundaries, M and L (sd01)  ← never used
- D5 BGC information (sd03)
- D6 HCR integration loci (Table S1)

## What we HAVE analyzed (all use STATIC M-phase only)
1. methylome × FIRE_M genome-wide (ρ=0.32) — now shown compartment-level only (oriC/within-core null)
2. Exposed vs Shielded FIRE_M
3. HCR-M1–10 integration concordance (fragile)
4. co-high FIRE∩methylation hotspots
5. compartment-A occupancy via PC_M (62%→4% T1→T2; OR 1.96/3.24)

## UNEXPLORED angles (gap)
| # | Angle | Data | Feasible now? | Value | Why |
|---|---|---|---|---|---|
| **G1** | **CID-boundary alignment**: do GCCGGC/AAGCCCG sites, the protection zone, or Exposed promoters coincide-with / avoid Deng's CID boundaries? | sd02 (local) | **YES** | **HIGH** | Paper itself names this as "future work" but data is in hand. Tests epigenome×3D at the *boundary/architecture* level — a structural finding that can REPLACE the dead FIRE-density proxy. |
| **G2** | **Dynamic co-change**: Δmethylation(T1→T2) vs Δcompartment(PC_M→PC_L) and ΔFIRE(M→L), per region | bins M+L (local) | **YES** | **HIGH** | Central claim "methylation relocation mirrors the 3D refolding" is currently only static overlap. Correlating the two *changes* quantifies it; partly answers the static-alignment caveat. |
| G3 | Do our methylation motifs (GCCGGC/AAGCCCG) appear among CID-boundary-enriched TF motifs? | sd01 | YES | MED | Direct motif-level link to boundaries |
| G4 | 3-way mediation: does methylation explain expression beyond FIRE, or is the weak meth–expr link mediated by FIRE? | local | YES | MED | Sharpens permissive thesis |
| G5 | compartment-switching genes (A↔B, M→L) ∩ Exposed-62 / demethylated | PC_M,PC_L per gene | YES | MED | Ties demethylation to refolding mechanistically |
| G6 | within-arm methylome–FIRE (ρ=0.264) driver | local | YES | LOW | Explain the residual arm correlation |
| G7 | BGC 3D positioning × methylation × expression | sd03 + local | YES | LOW | BGC methylation already shown compositional |

## Recommendation
Run **G1 (CID boundaries)** and **G2 (dynamic co-change)** now — both local, both HIGH value. G1 can convert the collapsed FIRE-proxy into a stronger, novel epigenome–architecture finding (boundary alignment); G2 quantitatively backs the paper's central "mirrors refolding" claim. G3–G5 as second tier.
