# Streptomyces属全種メチル化モチーフ保存率・ゲノム密度解析レポート

**作成日**: 2026-02-05
**プロジェクト**: *Streptomyces coelicolor* A3(2) M145 メチローム解析
**対象**: マニュスクリプト執筆者・共著者
**解析ディレクトリ**: `11_epigenome_integration/analysis/21_genuswide_motif_conservation/`
**スクリプト**: `11_epigenome_integration/scripts/streptomyces_motif_conservation_genuswide.py`

---

## Key Insights（主要なインサイト一覧）

| # | Insight | 対応Panel | 根拠となる数値 |
|---|---------|-----------|----------------|
| 1 | AAGCCCG認識R-M系はREBASE登録82種中未報告（配列自体は833種に普遍的に存在） | Panel A, B | 0/82 (0.0%, 95% CI: 0.0-4.5%); 833種中央値166 sites/Mb |
| 2 | CCGGとGCGGが属内最高保存率（各~23%）だが属共通ではない | Panel A | CCGG: 18/82 (22.0%), GCGG: 19/82 (23.2%) |
| 3 | GATCは属全体でO/E ~2.1と強い過剰表現を示す | Panel D | 中央値O/E=2.12, M145 O/E=1.99 |
| 4 | CGACNNNCTCCは最高のO/E比（~2.75）を示す | Panel D | 中央値O/E=2.75, M145 O/E=2.86 |
| 5 | AAGCCCGはゲノム上で過少表現（O/E=0.70）→ 回避シグナル | Panel B, D | M145 O/E=0.65, 属中央値O/E=0.70, パーセンタイル24% |
| 6 | R-M系は種特異的に進化し、属レベルでの共有は限定的 | Panel C | 最高保存率23.2%、大半の種は1-2モチーフのみ一致 |

---

## 1. 解析の目的

本研究（CCGG, AAGCCCG, GATC）および先行研究（Pisciotta 2023: GGCCGG, GCCCG; Fang 2022: GCGG, CGACNNNCTCC）で同定された7種のメチル化モチーフについて、2つの独立した視点から属内分布を評価する:

1. **R-M系保存率（REBASE）**: 各モチーフを認識するR-M系がStreptomyces属内でどの程度保存されているか？（酵素・機能レベル）
2. **ゲノムモチーフ密度（RefSeq）**: 各モチーフ配列がゲノム上にどの程度存在するか？GC含量から予想される密度との乖離は？（配列レベル）

---

## 2. データソース

| データソース | バージョン | Streptomyces種数 | 用途 |
|------------|-----------|----------------|------|
| **REBASE** | v602 (2026-01-28) | 224エントリ → **82種**（二名法集約） | R-M系保存率 (Panel A, C) |
| **NCBI RefSeq** | 2026-02 | **833種**（代表ゲノム） | モチーフ密度 (Panel B, D) |

### RefSeqゲノム選択基準
- NCBI Datasets API v2で全Streptomyces属complete genome assemblies取得
- 種レベルで重複排除（reference/representative genome優先 → N50最大）
- ゲノムサイズ: 3.9-15.0 Mb（中央値8.4 Mb）
- GC含量: 67.6-74.8%（中央値71.7%）

---

## 3. 解析対象モチーフ（7種）

| モチーフ | 修飾 | 出典 | 同定種 | 検出法 |
|---------|------|------|--------|--------|
| **CCGG** | 4mC | 本研究 | *S. coelicolor* M145 | Nanopore |
| **AAGCCCG** | 6mA | 本研究 | *S. coelicolor* M145 | Nanopore |
| **GATC** | 6mA | 本研究 | *S. coelicolor* M145 | Nanopore |
| **GGCCGG** | 5mC | Pisciotta 2023 | *S. coelicolor* M145 | BS-seq |
| **GCCCG** | 5mC | Pisciotta 2023 | *S. coelicolor* M145 | BS-seq |
| **GCGG** | 4mC | Fang 2022 | *S. roseosporus* L30 | SMRT |
| **CGACNNNCTCC** | 6mA | Fang 2022 | *S. roseosporus* L30 | SMRT |

**除外**: CG/CHG/CHH（汎用パターン、R-M系特異性なし。全ゲノムで飽和検出されるため比較に不適）

---

## 4. R-M系保存率（REBASE v602、82種）

### 【Figure】fig_motif_conservation_composite.pdf

**ファイル**: `analysis/21_genuswide_motif_conservation/fig_motif_conservation_composite.pdf`

### Panel A: R-M系認識配列の保存率

| モチーフ | 修飾 | 出典 | 一致種/82 | 保存率 | 95% CI | 代表的酵素 |
|---------|------|------|----------|--------|--------|-----------|
| **CCGG** | 4mC | 本研究 | 18/82 | **22.0%** | 14.4-32.1% | SauLPI, SfiI, SgrAI |
| **AAGCCCG** | 6mA | 本研究 | **0/82** | **0.0%** | 0.0-4.5% | — |
| **GATC** | 6mA | 本研究 | 11/82 | 13.4% | 7.7-22.4% | SalAI, SgrAII |
| **GGCCGG** | 5mC | Pisciotta 2023 | 2/82 | 2.4% | 0.7-8.5% | SfiI |
| **GCCCG** | 5mC | Pisciotta 2023 | 4/82 | 4.9% | 1.9-11.9% | SfiI, SrfI |
| **GCGG** | 4mC | Fang 2022 | 19/82 | **23.2%** | 15.4-33.4% | SacII, SgrBI |
| **CGACNNNCTCC** | 6mA | Fang 2022 | 1/82 | 1.2% | 0.2-6.6% | M.SfiL30I |

### Panel C: 種×モチーフ ヒートマップ

REBASE登録82種のうち、少なくとも1モチーフに一致するR-M系を持つ種を抽出し、階層的クラスタリング（Hamming距離、average法）で並び替え。

- 大半の種は1-2モチーフのみでR-M系を保有
- CCGG/GCGGの同時保有パターンが複数種で見られる
- AAGCCCGは全種で空白（R-M系なし）

---

## 5. ゲノムモチーフサイト密度（RefSeq 833種）

### ゲノムモチーフサイト密度とO/E比（7モチーフ × 833 RefSeq種）

**定義**: 密度 = モチーフ出現回数（順鎖+逆相補鎖）/ ゲノムサイズ (Mb)。O/E = 観測カウント / GC含量ベース期待カウント。

| モチーフ | 修飾 | 存在種/833 | 属中央値密度 (sites/Mb) [IQR] | M145密度 | M145 pctl | 属中央値O/E [IQR] | M145 O/E | 選択圧 |
|---------|------|-----------|----------------------------|---------|----------|----------------|----------|--------|
| **CCGG** | 4mC | 833/833 | 17,041 [15,684-18,208] | 17,441 | **59%** | 1.03 [0.98-1.07] | 1.03 | 中立 |
| **AAGCCCG** | 6mA | 833/833 | 166 [154-179] | 154 | **24%** | **0.70 [0.65-0.76]** | **0.65** | **過少（回避）** |
| **GATC** | 6mA | 833/833 | 5,498 [5,036-5,760] | 5,023 | **23%** | **2.12 [2.01-2.18]** | **1.99** | **過剰（維持）** |
| GGCCGG | 5mC | 833/833 | 4,230 [3,729-4,719] | 4,256 | **51%** | 1.00 [0.93-1.07] | 0.97 | 中立 |
| GCCCG | 5mC | 833/833 | 10,267 [9,713-10,865] | 10,422 | **56%** | 0.87 [0.85-0.89] | 0.86 | 軽度過少 |
| GCGG | 4mC | 833/833 | 32,013 [30,451-33,460] | 32,375 | **56%** | 0.97 [0.95-0.99] | 0.96 | 中立 |
| CGACNNNCTCC | 6mA | 833/833 | 233 [216-249] | 244 | **68%** | **2.75 [2.55-2.93]** | **2.86** | **過剰（強い正の選択）** |

→ 全7モチーフの配列が833種全てのゲノムに存在する。AAGCCCGは属全体でO/E=0.70（過少表現）、GATCとCGACNNNCTCCは強い過剰表現を示す。

### Panel D: ゲノムO/E比（GC含量ベース期待値との比較）

**O/E比の解釈**:
- O/E > 1: ゲノム上で**過剰表現**（GC含量から予想されるより多い）
- O/E < 1: ゲノム上で**過少表現**（回避シグナルの可能性）
- O/E = 1: ランダム期待通り

**主要な知見**:

1. **GATC（O/E=2.12）**: 属全体で強い過剰表現。ATリッチモチーフがGCリッチゲノム（~72%）中で期待以上に保持されており、機能的選択圧（Dam-like活性?）の存在を示唆
2. **CGACNNNCTCC（O/E=2.75）**: 最高のO/E比。長く複雑なモチーフが期待を大幅に上回る→ 強い正の選択圧
3. **AAGCCCG（O/E=0.70）**: 属全体で過少表現。R-M系の未同定と合わせ、このモチーフが回避されている可能性を示唆。M145でのパーセンタイルも24%と低い
4. **CCGG（O/E=1.03）、GGCCGG（O/E=1.00）**: ほぼランダム期待通り
5. **GCGG（O/E=0.97）、GCCCG（O/E=0.87）**: 軽度の過少表現

---

## 6. AAGCCCGの新規性に関する統合評価

| 評価軸 | 結果 | 解釈 |
|--------|------|------|
| R-M系保存率（REBASE） | **0/82種 (0.0%)** | REBASE登録種でこのモチーフ認識R-M系は未特性決定（属の9.8%のカバレッジ） |
| ゲノム密度パーセンタイル（M145） | **24%**（属内下位） | モチーフ配列自体は稀ではないが低め |
| O/E比 | **0.65**（M145）/ **0.70**（属中央値） | 属全体でゲノム回避傾向 |
| 他の0%保存率モチーフとの比較 | CGACNNNCTCC=1.2%で最低 | AAGCCCGのみ完全にゼロ |

→ **AAGCCCGは（1）配列は833種全てに存在するがR-M系がREBASE登録種で未報告、（2）ゲノム上で過少表現（O/E=0.70）を示す。SC_RS17645（M145）がこのモチーフ認識R-M系の初の候補MTase**

---

## 7. Limitation

- **REBASEカバレッジ**: 82/833種（9.8%）のみ実験的R-M系データあり。残り91%の種は未特性決定であり、保存率は下限推定
- **RefSeqゲノム品質**: Complete genomeのみ選択したが、アセンブリ品質のばらつきあり
- **種レベルの解像度**: 同種内の株間変異はカバーされていない（代表1ゲノム/種）
- **O/E計算の前提**: GC含量のみに基づく単純モデル。ジヌクレオチド頻度の偏りは考慮せず
- **密度 ≠ メチル化**: ゲノム上のモチーフ配列の存在は、実際にメチル化されていることを意味しない

---

## 8. 出力ファイル一覧

### Figure（PDF + SVG + PNG）
| ファイル | 内容 |
|---------|------|
| `fig_motif_conservation_composite.pdf/.svg/.png` | 4パネル統合図 |
| `panel_a_conservation.pdf` | Panel A: R-M保存率棒グラフ |
| `panel_b_density.pdf` | Panel B: モチーフ密度box+stripプロット |
| `panel_c_heatmap.pdf` | Panel C: 種×モチーフヒートマップ |

### データ
| ファイル | 内容 |
|---------|------|
| `rebase_motif_conservation_matrix.csv` | 7モチーフのREBASE保存率（一致酵素リスト付き） |
| `rebase_species_motif_matrix.csv` | 種×モチーフ二値マトリクス（Panel C用） |
| `motif_site_density_genuswide.csv` | 833種×7モチーフの密度・O/E（Panel B, D用） |
| `species_representative_accessions.tsv` | 833種の代表ゲノムアクセッション |
| `analysis_log.txt` | 解析サマリーログ |

---

## 9. 方法の詳細

### R-M系保存率の算出
1. REBASE Bairochフォーマットの`RS`フィールドから認識配列を抽出
2. 7モチーフそれぞれについて認識配列との一致を判定:
   - IUPAC曖昧塩基対応
   - 逆相補鎖も検索
   - 部分一致はspecificity score >= 0.6の場合のみ採用
3. 各種について認識配列が1つでも一致すればカウント
4. Wilson二項95%信頼区間を算出

### ゲノムモチーフ密度の算出
1. 833種の代表ゲノムFASTA配列を読み込み
2. 各モチーフについて順鎖+逆相補鎖でのカウント:
   - 固定モチーフ: `str.count()`
   - 縮退モチーフ（CGACNNNCTCC）: `re.findall()` with IUPAC→regex変換
3. 密度 = カウント / ゲノムサイズ(Mb)
4. GC含量ベース期待値 = P(motif_forward) + P(motif_rc) × ゲノム長
5. O/E = 観測カウント / 期待値
6. 並列処理: `multiprocessing.Pool`（8 workers）

---

## 参考文献

- Roberts, R. J., et al. (2023) REBASE. *Nucleic Acids Res.*, 51, D629-D635.
- Pisciotta, A., et al. (2023) *Scientific Reports* 13, 7038.
- Fang, G., et al. (2022) *Nature Biotechnology* 30, 1232-1239.

---

*最終更新: 2026-02-05*
