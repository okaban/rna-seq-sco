# Revision Log R1 — peer review response

**作業日:** 2026-05-08
**対象査読:** `reports/260508_peer_review_report.md`
**改訂対象ファイル:**
- `15_paper_figures/manuscript/full_manuscript.md`
- `15_paper_figures/manuscript/01_results.md`
- `15_paper_figures/manuscript/03_discussion.md`

---

## 概要

査読報告の Major Concerns (MC-1〜MC-6) と主な Minor Concerns を、優先順位 P1〜P7 として原稿 3 ファイルに反映した。捏造防止の方針として、knockout など本研究で実施していない実験データは追加せず、Limitations の明示と語彙置換のみで対応している。Supplementary_statistical_tests.tsv の登録値（特に T21 の TCS 置換検定 p = 0.83）と整合する形で記述を改めた。

---

## P1 — 因果言語 → 共在・相関フレームへの全面改訂

**MC-1 対応。** 全 3 ファイルの Abstract / Introduction / Results / Discussion で因果動詞・名詞を中立化した。ただし temporal decoupling (rho=0.136, p=0.299) を示す段落は「因果否定の根拠」として明示的に保持し、加えて Abstract と Discussion にもこの数値を補強的に追加して内部矛盾を解消した。

主要な置換（網羅ではない）:

| 改訂前 | 改訂後 | 箇所 |
| --- | --- | --- |
| four-layer "Gatekeeper" architecture | four-layer methylation-exclusion architecture | Abstract / Intro / Results / Discussion |
| channels methylation effects through a specific minority of regulators | spatially co-localizes methylation accessibility with a specific minority of regulators | Abstract |
| constitute a distributed vegetative-to-developmental switch | structurally co-localized with a distributed transcriptional transition pattern at the vegetative-to-developmental boundary | Abstract |
| DNA methylation does not control transcription directly | DNA methylation does not directly control transcription | Abstract |
| sole conduit through which methylation interacts with the transcriptional program | primary regulatory subset that structurally co-localizes with methylation-accessible DNA | Results / 01_results / Discussion |
| true, structurally mediated pathway through 57 specific regulatory genes | structurally mediated co-occurrence pattern between methylation accessibility and 57 specific regulatory genes | Introduction |
| Methylation's regulatory influence is thus channeled through a narrow bottleneck of 5.9% | Methylation accessibility is therefore associated with a narrow ~5.4% subset | Discussion |
| The activation bloc collectively drives the onset of morphogenesis | The activation bloc collectively tracks with the onset of morphogenesis | Results / 01_results |
| operates independently of all computationally characterized transcriptional networks | is not captured by any of the computationally characterized transcriptional networks examined here | Results / 01_results |
| coordinated, simultaneous switch | coordinated, simultaneous transition | Results / 01_results |
| Both blocs activate simultaneously at the vegetative-to-developmental transition (heading) | Both blocs show coordinated activation at the vegetative-to-developmental transition | Results / 01_results |
| developmental "gate" opens once at the vegetative-to-developmental boundary | Expression changes are concentrated at the vegetative-to-developmental boundary | Results / 01_results |
| an activation bloc (~35 genes) driving developmental programs and a repression bloc (~26 genes) shutting down vegetative programs | an activation bloc (~35 genes) whose upregulation accompanies developmental programs and a repression bloc (~26 genes) whose downregulation accompanies the shutdown of vegetative programs | Discussion summary / Abstract |
| constituting a distributed developmental switch | constituting a structurally co-localized regulatory pattern rather than a methylation-driven switch | Discussion summary |
| The entire apparent genome-wide methylation-expression correlation is a statistical artifact | The apparent genome-wide methylation-expression correlations are statistical artifacts | Discussion (Reframing) |
| The Gatekeeper model is shaped... / The Gatekeeper model was shaped... | The four-layer methylation-exclusion model is/was shaped... | Results / 01_results |

加えて Discussion の "Reframing" 段落に、AAGCCCG 系では geographic confound が小さいものの完全には排除できない旨を 1 文追加。

---

## P2 — TCS 7/7 asymmetry の再フレーミング

**MC-2 対応。** Results と Discussion 両方で機構的解釈を圧縮し、p = 0.83 (Supplementary_statistical_tests.tsv T21, permutation, label shuffle, n=10,000) を明記した上で、descriptive observation only として記述を残した。

- Results 旧文: "The probability of this perfect asymmetry occurring by chance is less than 10^-3 (binomial test). ... This asymmetric arrangement is functionally significant ..."
- Results 新文: "We emphasize that this is a descriptive observation, not a statistically supported pattern: a permutation test ... yielded p = 0.83 (Supplementary_statistical_tests.tsv, T21) ... The pattern is therefore preliminary and underpowered (n = 7 pairs); we report it for completeness without inferring that methylation accessibility selectively gates the sensory versus response arm of these systems."
- Discussion 旧段: 「7/7 perfect asymmetry が methylation system が TCS の sensing 機能を gating する可能性を示唆する」
- Discussion 新段: 「permutation test で p = 0.83。chance と区別不能。7 pairs では mechanistic interpretation を支持しない。larger TCS datasets での再検討の価値はあるが、現時点では selective TCS gating の根拠とはしない」

---

## P3 — Protection zone 33% / 67% の根拠提示と定性化

**MC-3 対応。** AUC 分解では 51%/49% となり 33%/67% を導出できないとの査読指摘に対し、(a) 数値を「partial component」「larger non-sequence component」と定性化し、(b) 2 種類の量化指標（motif depletion 0.74–0.76 ≈ 約 1/3 vs AUC-based decomposition 51%/49%）を本文中に明示し、(c) 単一 apportionment の主張を撤回した。

- Results: "These sequence contributions account for approximately 33% ..." → 「partial component」と再定義し、両指標の数値と意味の差異（DNA レベルの motif 欠落 vs 分類器の above-chance 精度の分解）を明示。
- Results 続き: "The remaining ~67% of the protection cannot be explained by DNA sequence" → "The non-sequence component of the protection cannot be explained by DNA sequence alone"
- Discussion (NAP 段落): "The two-tier nature of the protection — ~33% evolutionary counter-selection, ~67% protein occupancy" → 「partial evolutionary counter-selection of motifs combined with extensive promoter protein occupancy ... we therefore avoid asserting a single percentage apportionment」
- Abstract / Introduction / Discussion summary でも "(~33%) / (~67%)" の括弧表記を削除し、定性記述に統一。

---

## P4 — BGC 内部比較への geographic stratification 留意点

**MC-5 対応。** Results "Within the BGC compartment" 段落末に caveat を 1 文 (実際には 2 文に拡張) 追加。

新規追加文:
> "We caution that the within-BGC comparison was not stratified for chromosomal sub-position: the four major BGCs are all core-located while minor BGCs span both core and arms, so the same Simpson's-paradox geometry that nullified the genome-wide BGC enrichment may also contribute here. The within-BGC 4mC signature should therefore be interpreted as a tentative, descriptive pattern requiring sub-region-stratified validation."

加えて p = 0.062 が「above the conventional 0.05 threshold」であることを文中で明示。

(注) 01_results.md は BGC 段落を含んでいないため、修正は full_manuscript.md のみ。

---

## P5 — TSS proxy と感度解析への言及

**MC-6 対応。** Methods の "Transcription start site assignment" 段落と Discussion の Limitations 段落の両方に注記を追加。

Methods 追加文:
> "Because annotation start coordinates were used as a proxy for the true TSS for the remaining ~75% of regulatory genes, the nearest-methylation-distance metric and the 293 bp boundary depend on TSS-position accuracy; bacterial 5′ UTR lengths can vary from 0 to several hundred bp, which can shift the inferred boundary. To assess robustness, we performed a Pfam-domain-based replication using the experimental-TSS subset only (Supplementary Note), which yielded an AUC of 0.923 and an optimal boundary of 318 bp — quantitatively consistent with the keyword-based analysis on the full set (AUC = 0.917, 293 bp). A genome-wide experimental-TSS analysis is identified as a key future sensitivity test."

Discussion Limitations にも同等の "Fifth, ..." 項目として 1 段落追加。

(注) 01_results.md は Methods 節を含まないため、修正は full_manuscript.md と 03_discussion.md のみ。

---

## P6 — 遺伝学的介入実験なしの限界を Limitations に明記

**MC-4 対応。** Discussion の Limitations 段落の冒頭を全面改訂し、「No genetic intervention experiments — methyltransferase knockouts, TF knockouts, ChIP-seq, or DAP-seq — were performed in this study」と明記。さらに「the four-layer methylation-exclusion architecture should therefore be regarded as an observational organizational model rather than a causally established mechanism」を追加。

最後に "epiphenomenon の可能性は除外できない" を明示:
> "Until such experiments are available, the possibility that the observed protection zone and the developmental transition are co-occurring epiphenomena driven by an independent factor (e.g., NAP occupancy dynamics or chromosome 3D reorganization) cannot be excluded."

加えて、最終段落の「The central insight ... methylation can organize the regulatory genome through structural partitioning rather than direct transcriptional control」を「The central observation ... methylation can be spatially organized relative to the regulatory genome through structural partitioning, even in the absence of direct transcriptional control」と再フレーム。

---

## P7 — 数値・術語・プレースホルダー修正

| 項目 | 改訂前 | 改訂後 | 該当ファイル |
| --- | --- | --- | --- |
| 5.9% → 5.4% (Introduction) | "the remaining 5.9% of regulatory genes" | "the remaining 5.4% of regulatory genes" | full_manuscript.md |
| 5.9% → 5.4% (Discussion Reframing) | "narrow bottleneck of 5.9% of regulatory genes" | "narrow ~5.4% subset of regulatory genes" | full_manuscript.md, 03_discussion.md |
| 61 of 57 → 56 of 57 | "four modules covering 61 of 57 genes" | "four modules covering 56 of 57 genes" | full_manuscript.md, 01_results.md |
| Double comma | "genome-wide,, we asked" | "genome-wide, we asked" | full_manuscript.md (mn-7) |
| Gatekeeper 術語 | "Gatekeeper architecture" / "Gatekeeper model" (>10 instances) | "four-layer methylation-exclusion architecture" / "four-layer methylation-exclusion model" | 全 3 ファイル (mn-2) |

[CITE] プレースホルダーは末尾の Citation Mapping 表により全て参考文献番号にマップ済み (full_manuscript.md L313–L331)。本改訂では追加削減せず、最終整形時の機械的置換に委ねる。

---

## 確認事項 (今回の改訂では未対応)

以下の Minor Concerns は時間的・実験的制約から本ラウンドでは反映を見送った。次回以降の改訂で対応する候補:

- **mn-3:** SC_RS17645 を "orphan N6-adenine methyltransferase (HsdM-type fold; no cognate HsdR/HsdS)" と Methods で明記する書き換え。現状は「HsdM-type methyltransferase」のままで、Type I R-M に紐づく印象を残している。
- **mn-4:** STREME 解析結果との整合性 (unassigned 6mA sites に高信頼 motif なし) を 6mA セクションに明記する追記。現状の "suggesting the existence of at least one additional uncharacterized MTase" を「STREME による de novo motif discovery では新規 motif 同定に至らず」と補強する余地。
- **mn-5:** "Both blocs engage at the same time" の equivalence testing (TOST) 追加。現行版では "true simultaneity cannot be distinguished from a rapid cascade" と限界を明示する形で部分対応。
- **mn-6:** [CITE: ...] プレースホルダーの最終的な番号置換。Citation mapping 表は完備済みだが、本文中の置換は最終整形時に実施する想定。
- **mn-8:** Methods "Regulatory gene annotation" の重複記述 (keyword-matching) の整理。現状そのまま。
- **mn-9:** KEGG Quorum Sensing enrichment が S. coelicolor では non-AHL 制御を反映する旨の文脈説明。現状そのまま。
- **mn-10:** 3 timepoints 設計が「同時応答」主張に与える影響の議論拡張。Limitations 第 4 項目にて部分的に言及済み。
- **Pisciotta et al. 2023:** 本文中 (Discussion) では引用しているが References セクションに登録なし。最終整形時に追加が必要。

---

## 実施チェック

- [x] full_manuscript.md: P1–P7 全て反映
- [x] 01_results.md: 該当する P1, P2, P3, P7 を反映
- [x] 03_discussion.md: 該当する P1, P2, P3, P5, P6, P7 を反映
- [x] grep audit: "Gatekeeper", "5.9%", "61 of 57", "sole conduit", "channels methylation effects", "channeled through", "drives the onset", "operates independently of all" 全て 0 ヒット
- [x] Supplementary_statistical_tests.tsv の値 (T21: p = 0.83) と整合
- [x] knockout など捏造データなし

---

## 次のアクション候補

1. mn-3 (SC_RS17645 orphan MTase 表記) の Methods 修正
2. mn-4 (STREME negative result) の Results 補強
3. References セクションに Pisciotta et al. 2023 を追加
4. 査読レポート [MC-1]〜[MC-6] のうち [MC-4] (実験的検証) について、reviewer rebuttal letter のドラフトでより詳しい "future experiments" を提示
