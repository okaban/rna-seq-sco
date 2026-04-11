# T3 追加シーケンシングデータによる再解析レポート

**プロジェクト:** *Streptomyces coelicolor* A3(2) M145 エピゲノム-トランスクリプトーム統合解析
**日付:** 2026-02-24
**解析ディレクトリ:** `/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/`

---

## 1. 背景と目的

### 1.1 問題点

初回解析（2026年1月）において、T3タイムポイントのサンプル3-2のNanoporeシーケンシングカバレッジが14.3xと極めて低く、T3のメチル化サイト検出に深刻な影響を与えていた。特に4mCのT3検出数（227サイト）はT1（1,366）・T2（1,701）と比較して著しく少なく、時系列解析の信頼性が限定的であった。

### 1.2 対応

2026年2月に以下の追加シーケンシングを実施し、データを統合した。

| サンプル | 追加シーケンシングラン | 追加リード数 |
|---------|---------------------|------------|
| 3-2 (barcode52) | 260216_takeda + 260217_takeda | ~450,000 |
| 3-4 (barcode53) | 260216_takeda | ~38,000 |

### 1.3 解析パイプライン変更

旧データでは複数の解析手法を試行したが、今回の再解析では以下に統一した。

| パラメータ | 旧（v1b weighted） | 新（v2 current） |
|-----------|-------------------|-----------------|
| コンセンサス手法 | カバレッジ加重平均 | 単純平均（unweighted） |
| MIN_REPS | 2 | 3 |
| MIN_COVERAGE | 5 | 5 |
| MIN_TOTAL_COVERAGE | 30 | 30 |
| MIN_MOD_FREQ | 50% | 50% |
| T3レプリケート | 3（3-2旧, 3-3, 3-4旧） | 3（3-2新, 3-3, 3-4新） |

MIN_REPS=3（3/3レプリケート一致を要求）とすることで、カバレッジ加重の必要性を排除し、より保守的で再現性の高い基準とした。

---

## 2. シーケンシングデータの改善

### 2.1 アライメント統計

| サンプル | タイムポイント | 全リード数 | マップリード | アライメント率 | カバレッジ |
|---------|-------------|-----------|------------|-------------|----------|
| 1-1 | T1 | 152,076 | 133,171 | 87.6% | 87.7x |
| 1-2 | T1 | 154,159 | 141,771 | 92.0% | 82.0x |
| 1-3 | T1 | 218,526 | 190,806 | 87.3% | 124.9x |
| 2-1 | T2 | 365,498 | 336,297 | 92.0% | 99.2x |
| 2-3 | T2 | 434,470 | 382,580 | 88.1% | 89.9x |
| 2-4 | T2 | 370,101 | 340,948 | 92.1% | 86.1x |
| **3-2** | **T3** | **709,346** | **603,783** | **85.1%** | **57.6x** |
| 3-3 | T3 | 360,965 | 309,496 | 85.7% | 79.6x |
| **3-4** | **T3** | **222,525** | **193,108** | **86.8%** | **75.4x** |

### 2.2 サンプル3-2の改善

| 指標 | 旧データ | 新データ（マージ後） | 改善率 |
|------|---------|-------------------|--------|
| 総リード数 | ~259,000 | 709,346 | **+174%** |
| マップリード | ~216,000 | 603,783 | **+179%** |
| 推定カバレッジ | ~14.3x | 57.6x | **+303%** |
| modkit処理リード | ~140-200k | 606,803 | **~3-4倍** |

---

## 3. メチル化サイト検出への影響

### 3.1 High-Confidence Sites の変化

| メチル化型 | TP | v1b (旧データ, weighted, MIN_REPS=2) | v2 (新データ, unweighted, MIN_REPS=3) | 変化 |
|-----------|-----|--------------------------------------|---------------------------------------|------|
| 6mA | T1 | 1,889 | 1,934 | +45 (+2.4%) |
| 6mA | T2 | 2,102 | 2,120 | +18 (+0.9%) |
| 6mA | T3 | 2,257 | 2,295 | +38 (+1.7%) |
| 4mC | T1 | 1,995 | 1,987 | -8 (-0.4%) |
| 4mC | T2 | 2,458 | 2,446 | -12 (-0.5%) |
| **4mC** | **T3** | **1,077** | **1,647** | **+570 (+52.9%)** |
| **合計** | | **11,778** | **12,429** | **+651 (+5.5%)** |

最大の変化は**4mC T3の+570サイト**（+52.9%）である。追加シーケンシングによりサンプル3-2のカバレッジが大幅に改善され、3/3レプリケート一致の条件を満たすサイトが増加した。

### 3.2 ユニーク位置数（Census）

| メチル化型 | ユニーク位置数 |
|-----------|-------------|
| 4mC | 2,693 |
| 6mA | 3,248 |
| **合計** | **5,941** |

### 3.3 モチーフ帰属

**4mC（2,693位置）:**

| モチーフ | サイト数 | 割合 |
|---------|---------|------|
| TGGCCGGC | 1,717 | 63.8% |
| AAGCCCG | 814 | 30.2% |
| CCGG (その他) | 146 | 5.4% |
| GGCCGG | 16 | 0.6% |

**6mA（3,248位置）:**

| モチーフ | サイト数 | 割合 |
|---------|---------|------|
| 未帰属 | 2,130 | 65.6% |
| AAGCCCG | 362 | 11.1% |
| CCGKCA | 150 | 4.6% |
| CCGG | 126 | 3.9% |
| GCCG | 108 | 3.3% |
| CCGC | 106 | 3.3% |
| CCSGG | 102 | 3.1% |
| GAACCGG | 62 | 1.9% |
| CGGCAACC | 59 | 1.8% |
| GATC | 39 | 1.2% |

---

## 4. メチル化-発現相関解析

### 4.1 ゲノムワイド Spearman 相関

| メチル化型 | 比較 | n | Spearman r | p値 | 有意 |
|-----------|------|---|-----------|-----|------|
| **4mC** | **T2vsT1** | **626** | **0.125** | **0.0017** | **Yes** |
| 6mA | T2vsT1 | 621 | -0.038 | 0.348 | No |
| 6mA | T3vsT1 | 662 | -0.093 | 0.016 | Yes |
| 6mA | T3vsT2 | 667 | -0.092 | 0.018 | Yes |
| 4mC | T3vsT1 | 542 | -0.003 | 0.937 | No |
| 4mC | T3vsT2 | 600 | -0.030 | 0.468 | No |

### 4.2 旧データとの相関比較

| メチル化型 | 比較 | v0 r (n) | v1a r (n) | **v2 r (n)** | 変化 |
|-----------|------|----------|----------|-------------|------|
| 4mC | T2vsT1 | 0.172 (474) | 0.172 (474) | **0.125 (626)** | r低下、n増加、依然有意 |
| 6mA | T2vsT1 | 0.169 (361) | 0.169 (361) | **-0.038 (621)** | **有意性消失** |
| 6mA | T3vsT1 | -0.050 (330) | 0.064 (449) | **-0.093 (662)** | 負の相関が出現 |

**重要な変化:** 6mA T2vsT1の正の相関（r=0.169, p=0.001）は、サンプルサイズ増加後に**完全に消失**した（r=-0.038, p=0.348）。一方、4mC T2vsT1の正の相関はサンプルサイズ倍増後も維持（r=0.125, p=0.002）。これは4mC相関のロバスト性を支持する。

### 4.3 Spurious Correlation Validation（4mC T2vsT1）

| 検証テスト | 結果 | 解釈 |
|-----------|------|------|
| 観測 r | 0.1248 | 正の相関 |
| Permutation p | 0.0020 | 有意 |
| Partial r（交絡因子制御後） | 0.1388（111.2%保持） | ロバスト |
| Bootstrap 95% CI | [0.046, 0.194] | ゼロを含まない |
| Z vs random | 3.19 | 強いシグナル |

**結論:** 4mC T2vsT1相関は permutation test、partial correlation、bootstrap CI の全てで検証に耐え、**spurious（偽の相関）ではない**。

### 4.4 TSS-anchored 相関解析（67テスト、FDR補正後）

全67テスト中、FDR補正後に有意なのは**1件のみ**:

- **6mA, T3vsT1, gene body**: p_adj=0.049 (Spearman r=-0.081, n=936)

距離帯別解析、プロモーター vs gene body、SARP zone（-60 to -20 from TSS）いずれにおいても、FDR補正後に系統的な相関パターンは認められない。

---

## 5. Differentially Methylated Genes (DMGs)

### 5.1 DMG数

| 比較 | 6mA hyper | 6mA hypo | 4mC hyper | 4mC hypo | 合計 |
|------|-----------|----------|-----------|----------|------|
| T2vsT1 | 159 | 144 | 188 | 96 | 566 |
| T3vsT1 | 208 | 153 | 79 | 211 | 619 |
| T3vsT2 | 198 | 148 | 56 | 275 | 648 |

### 5.2 DEG-DMG Overlap

| 比較 | DEGs | DMGs | Overlap | Odds ratio | p値 | 有意 |
|------|------|------|---------|-----------|-----|------|
| T2vsT1 | 3,848 | 566 | 267 | 1.00 | 1.0 | No |
| **T3vsT1** | **4,841** | **619** | **407** | **1.35** | **5.6e-4** | **Yes** |
| T3vsT2 | 3,507 | 648 | 301 | 1.16 | 0.069 | No |

T3vsT1のみDEG-DMGの有意な重複あり（OR=1.35）。T2vsT1は完全にランダム期待値と一致（OR=1.00）。

### 5.3 DMG機能エンリッチメント

| 比較 | KEGG sig | GO sig | COG sig |
|------|----------|--------|---------|
| T2vsT1 | 0 | 0 | 0 |
| T3vsT1 | 0 | 1 | 0 |
| T3vsT2 | 0 | 1 | 0 |

**結論:** DMGsはKEGG/GO/COGいずれの機能カテゴリにおいてもほぼエンリッチメントを示さない。メチル化変動は遺伝子機能に対して**非選択的（non-selective）**である。

### 5.4 Coordinated Genes（メチル化・発現同時変動遺伝子）の機能エンリッチメント

| 比較 | グループ | 遺伝子数 | KEGG sig | GO sig | COG sig |
|------|---------|---------|----------|--------|---------|
| T2vsT1 | 全coordinated | 339 | 0 | 0 | 1 |
| T2vsT1 | positive | 114 | 0 | 0 | 0 |
| T2vsT1 | negative | 98 | 0 | 0 | 1 |
| T3vsT1 | 全coordinated | 273 | 0 | 0 | 0 |
| T3vsT2 | 全coordinated | 181 | 0 | 0 | 0 |

Coordinated genesにおいても機能的偏りは認められない。

---

## 6. Expanded Motif Search

### 6.1 Final Motif Judgment

| モチーフ | メチル化型 | サイト数 | スコア | 分類 | 発見方法 |
|---------|-----------|---------|--------|------|---------|
| **AAGCCCG** | **4mC+6mA (dual)** | **973+360** | **10/10** | **HIGH-CONFIDENCE NOVEL** | MEME-1 (6mA), STREME (4mC) |
| **CCGKCA** | 6mA | 153 | 7/10 | CANDIDATE | MEME-2 (6mA) |
| CCGG/GGCCGG/TGGCCGGC | 4mC | 2,034 | 4/10 | ESTABLISHED | MEME-1 (4mC) |
| GAACCGG | 6mA | 61 | 4/10 | ESTABLISHED/INSUFFICIENT | STREME-4 |
| GATC | 6mA | 37 | 4/10 | ESTABLISHED/INSUFFICIENT | Dam-type |
| CGGCAACC | 6mA | 57 | 3/10 | ESTABLISHED/INSUFFICIENT | STREME-3 |

### 6.2 AAGCCCGモチーフの特記事項

- **Dual modification**: 4mC（973サイト、36.3%）+ 6mA（360サイト、11.2%）の二重メチル化
- **O/E ratio**: 0.65（ゲノム中で回避される配列 → selection pressure）
- **REBASE**: 0/82マッチ（既知モチーフに該当なし → 新規）
- **候補MTase**: SC_RS17645（N-6 DNA methylase）
- **Coordinated genesプロモーター enrichment**: OR=12.89, p=1.1e-107

### 6.3 モチーフ時系列動態

| メチル化型 | モチーフ | T1 sites | T2 sites | T3 sites | 変動パターン |
|-----------|---------|----------|----------|----------|------------|
| 4mC | CCGG | 1,400 (70.5%) | 452 (68.1%) | 27 (64.3%) | 安定 |
| 4mC | AAGCCCG | 698 (35.1%) | 257 (38.7%) | 21 (50.0%) | 微増傾向 |
| 6mA | AAGCCCG | 260 (13.4%) | 64 (8.9%) | 38 (6.4%) | 低下傾向 |
| 6mA | GATC | 19 (1.0%) | 10 (1.4%) | 10 (1.7%) | 微増傾向 |

---

## 7. TF Binding Site メチル化解析

### 7.1 BS-メチル化空間的重複

| 解析窓 | BS数 | メチル化BSの数 | 割合 | 4mC fold | 6mA fold | 全体p値 |
|--------|------|-------------|------|---------|---------|---------|
| 直接重複 | 782 | 63 | 8.1% | 0.66 | 0.66 | 6.8e-5 |
| ±50bp | 782 | 100 | 12.8% | 0.76 | 0.73 | 3.6e-4 |
| ±200bp | 782 | 210 | 26.9% | 0.83 | 0.80 | 2.2e-4 |

**結論:** TF結合部位におけるメチル化は有意に**枯渇**（depleted）している。これはメチル化がTF結合を阻害する、あるいはTF結合がメチル化を排除する制約の存在を示唆する。

### 7.2 BS近傍メチル化の時系列変動

| パターン | サイト数 | 割合 |
|---------|---------|------|
| Lost_T2T3（T1で検出→T2/T3で消失） | 94 | 75.2% |
| Gained（T2/T3で新規出現） | 31 | 24.8% |

BS近傍のメチル化は大部分がT1特異的で、分化進行に伴い消失する。

### 7.3 BSメチル化 vs ターゲット発現相関

- Mann-Whitney: p=0.083 (T2vsT1), p=0.645 (T3vsT1)
- **有意な相関は認められない** → BSメチル化は直接的な発現制御の証拠を欠く

### 7.4 TFファミリー別メチル化率

| TFファミリー | BS数 | メチル化BS | 率 |
|-------------|------|----------|-----|
| LuxR | 13 | 4 | 30.8% |
| MarR | 7 | 2 | 28.6% |
| MerR | 5 | 1 | 20.0% |
| TCS response regulator | 80 | 14 | 17.5% |
| Sigma | 153 | 19 | 12.4% |
| GntR | 23 | 2 | 8.7% |
| TetR | 18 | 1 | 5.6% |

---

## 8. BGC・GRNの詳細解析

### 8.1 Act/Red BGCのメチル化状態

| BGC | 遺伝子数 | メチル化遺伝子 | メチル化サイト/遺伝子 | 発現変動 |
|-----|---------|-------------|-------------------|---------|
| Act | 20 | 3 | 各1サイト | T3で劇的に上昇（actII-ORF2 +7.7, actVA-ORF1 +10.9） |
| Red | 11 | 3 | 各1サイト | T2で劇的に上昇（redD +5.0, redL +4.2） |

**結論:** Act・Red両BGCは事実上**メチル化フリー**であり、二次代謝制御における直接的なエピジェネティック制御の証拠はない。

### 8.2 主要GRN TFのメチル化状態

37の文献既知TFのうち、プロモーターにメチル化を持つのは**5 TFのみ**:

| TF | メチル化状態 | 発現変動 (T2vsT1) | 備考 |
|----|------------|-------------------|------|
| bldN | Lost in T2 (2→1→2サイト) | -- | 形態分化の鍵TF |
| afsR | Stable (1サイト全TP) | -0.26 (NS) | グローバルレギュレーター |
| afsS | Stable (1サイト T1/T2) | -- | afsR下流 |
| redZ | Stable (2サイト T1/T2) | -2.25 (有意) | **唯一の協調変動TF** |

redZはプロモーターの6mAメチル化と発現低下の協調変動を示す唯一のGRN TFであるが、redDは逆にlog2FC=+4.77と大幅に上昇しており、absA2の脱抑制がredZ抑制を上回る。

### 8.3 -10 box メチル化

| 検定 | 観測値 | 期待値 | p値 | 解釈 |
|------|--------|--------|-----|------|
| 長さベースnull | 32遺伝子 | 58.1 | <0.001 | **有意に枯渇** |
| モチーフベースnull | 32遺伝子 | 32.1 | 0.536 | 期待値通り |
| SARP zone Fisher | 0 | 0.42 | 0.655 | NS |

-10 box内のメチル化は長さ期待値の55%しかなく有意に枯渇しているが、-10 box配列自体のモチーフ構成を考慮するとランダム期待値と一致する。これはメチル化モチーフ（CCGG等）と-10 boxモチーフの配列的非一致によるものと解釈できる。

---

## 9. MTase候補の発現動態

| MTase | 推定モチーフ | メチル化型 | 発現傾向 | 信頼度 |
|-------|-----------|-----------|---------|--------|
| SC_RS17645 | AAGCCCG | 6mA (+4mC dual) | DOWN_sustained (245→52→145) | HIGH |
| SC_RS19770 | CCGG/TGGCCGGC | 4mC | DOWN_then_UP (11→5→62) | MEDIUM |
| SC_RS36410 | CCGG/TGGCCGGC | 4mC | UP_sustained (1→4→65) | MEDIUM |
| SC_RS28835 | CCGKCA? | 6mA (BREX-2) | DOWN_sustained (604→395→382) | LOW-MEDIUM |
| SC_RS35335 | Unknown | 6mA (BREX-2) | STABLE (1666→1989→980) | LOW-MEDIUM |

SC_RS17645（AAGCCCG MTase候補）の発現低下はAAGCCCG部位の6mA頻度低下と時間的に一致する。

---

## 10. 主要結論のまとめ

### 10.1 新データによる結論の変化

| 知見 | 旧解析での結論 | 新解析での結論 | 変化 |
|------|--------------|--------------|------|
| 4mC T2vsT1 正の相関 | r=0.172, p=2e-4 | r=0.125, p=0.002 | **維持（ロバスト）** |
| 6mA T2vsT1 正の相関 | r=0.169, p=0.001 | r=-0.038, p=0.35 | **消失（非再現）** |
| T3サイト検出 | 4mC T3=227-1,077 | 4mC T3=1,647 | **大幅改善** |
| DMG機能非選択性 | エンリッチメントなし | エンリッチメントなし | **維持** |
| AAGCCCGの新規性 | Score 10/10 | Score 10/10 | **維持** |

### 10.2 確立された結論

1. **4mCプロモーターメチル化は遺伝子発現と正に相関する**（T2vsT1、permutation/bootstrap検証済み）
2. **メチル化変動は遺伝子機能に対して非選択的**（KEGG/GO/COGエンリッチメントなし）
3. **AAGCCCGは高信頼度の新規dual-modificationモチーフ**（4mC+6mA、REBASE未登録、O/E=0.65）
4. **TF結合部位はメチル化から保護されている**（0.66倍、p<0.001）
5. **BGC（Act/Red）はメチル化フリー** → 二次代謝の直接的エピジェネティック制御なし
6. **FDR補正後、TSS-anchored相関はほぼ全てNS**（67テスト中1テストのみ有意）
7. **-10 box内メチル化は枯渇**しているが、これはモチーフの配列的不一致による

### 10.3 新たに得られた知見

1. **6mA T2vsT1相関は非再現的** → 旧解析でのサンプルサイズ不足による偽陽性の可能性
2. **T3vsT1で有意なDEG-DMG重複**（OR=1.35, p=5.6e-4）→ 長期培養（T3）ではメチル化変動と発現変動が弱く連動
3. **BS近傍メチル化の75%がT1特異的**（Lost_T2T3）→ 初期培養での一過的メチル化
4. **4mC T3の大幅な検出増**（+570サイト）により、T3でのCCGG/AAGCCCG 4mCメチル化が過小評価されていたことが明らかに

---

## 11. 再解析で実行されたスクリプト一覧

計26スクリプトを再実行した。

| # | スクリプト | 出力ディレクトリ | 主要出力 |
|---|----------|----------------|---------|
| 1 | 01_expanded_motif_discovery.py | 23_expanded_motif_search | motif_temporal_dynamics.csv |
| 2 | 03_mtase_motif_prediction.py | 23_expanded_motif_search | orphan_mtase_inventory.csv |
| 3 | 04_conservation_analysis.py | 23_expanded_motif_search | conservation results |
| 4 | 05_integration_final_judgment.py | 23_expanded_motif_search | integration_final_scoring.csv |
| 5 | 06_rebase_refseq_conservation_figure.py | 23_expanded_motif_search | figures |
| 6 | 01_TF_BS_methylation_analysis.py | 13_TF_binding-site | T1-T5 tables, F1-F4 figures |
| 7 | 01_DMG_functional_enrichment.py | 14_DMG_enrichment | DMG_summary.tsv |
| 8 | 02_DMG_selectivity_deep_analysis.py | 14_DMG_selectivity | enrichment tables |
| 9 | T2vsT1_coordinated_genes.py | 01_integration | T2vsT1_coordinated_genes.csv (808 genes) |
| 10 | t3_coordinated_analysis.py | 16_t3_coordinated | T3vsT1 (283), T3vsT2 (184) |
| 11 | DEG_DMG_overlap_analysis.py | 03_overlap_analysis | overlap statistics + Venn |
| 12 | tss_analyses.py | 18_tss_analyses | 67 FDR-corrected tests |
| 13 | spurious_correlation_validation.py | 09_spurious_validation | validation report |
| 14 | null_model_simulation.py | 18_tss_analyses | -10 box null model |
| 15 | temporal_motif_expression_analysis.py | 18_tss_analyses | gained/lost motif bias |
| 16 | motif_tss_spatial_analysis.py | 18_tss_analyses | spatial analysis |
| 17 | mtase_methylation_dynamics.py | 23_expanded_motif_search | MTase dynamics |
| 18 | sarp_methylation_spatial.py | 18_tss_analyses | SARP zone analysis |
| 19 | aagcccg_promoter_analysis.py | 14_aagcccg_promoter_analysis | OR=12.89 enrichment |
| 20 | act_bgc_epigenetic_analysis.py | 15_act_bgc_epigenetic | BGC methylation |
| 21 | tf_methylation_enrichment.py | 18_tss_analyses | TF enrichment test |
| 22 | tf_methylation_target_triad.py | 12_grn_tf_methylation | redZ triad |
| 23 | tf_promoter_methylation_analysis.py | 12_grn_tf_methylation | 37 TF scan |
| 24 | 02_DMG_selectivity_deep_analysis.py | 14_DMG_selectivity | 再実行（coordinated更新後） |
| 25 | 03_mtase_motif_prediction.py | 23_expanded_motif_search | 再実行（Phase 1更新後） |
| 26 | 05_integration_final_judgment.py | 23_expanded_motif_search | 再実行（Phase 1更新後） |

---

## 12. 制限事項

1. **サンプル3-2のカバレッジ（57.6x）**は他のサンプル（80-125x）よりまだ低い。T3タイムポイントの4mC検出は引き続き保守的な推定値である可能性がある。
2. **タイムポイント数（n=3）**が限られているため、MTase発現-メチル化頻度の時系列相関は統計的検出力が不足している。
3. **STREME解析**（de novo motif discovery の一部）は外部ツール依存であり、今回の再解析ではSTREMEの再実行は行っていない。ただしSTREME結果に影響する入力データ（census files）は更新済み。
4. 相関係数の大きさ（r=0.125）は小さく、4mCメチル化が発現に与える影響は限定的である。生物学的意義の解釈には注意が必要。

---

## 13. 出力ファイル

本レポート: `11_epigenome_integration/analysis/01_integration/REANALYSIS_REPORT_260224.md`

主要出力テーブル:
- `01_integration/high_confidence_sites_weighted.csv` — 12,429サイト
- `01_integration/integrated_methyl_expression_weighted.csv` — 8,083遺伝子
- `01_integration/correlation_analysis_weighted.csv` — 6比較のSpearman相関
- `01_integration/T2vsT1_coordinated_genes.csv` — 808遺伝子ペア
- `03_overlap_analysis/DEG_DMG_overlap_statistics.csv` — DEG-DMG重複統計
- `09_spurious_validation/SPURIOUS_VALIDATION_REPORT.md` — 検証レポート
- `23_expanded_motif_search/integration_final_scoring.csv` — モチーフ最終判定
- `13_TF_binding-site/.../tables/T1-T5` — TF BS methylation 5テーブル
- `14_DMG_enrichment/.../tables/DMG_summary.tsv` — DMG機能エンリッチメント
