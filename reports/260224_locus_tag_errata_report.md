# Errata: TF locus_tag 系統的マッピングエラーの発見と修正

**作成日**: 2026-02-24
**最終更新**: 2026-02-24
**プロジェクト**: *Streptomyces coelicolor* A3(2) M145 トランスクリプトーム-エピゲノム統合解析
**関連ディレクトリ**: 全解析ステップ（05〜13）

---

## Key Insights（主要なインサイト一覧）

| # | Insight | 対応Figure | 根拠となる数値 |
|---|---------|-----------|---------------|
| 1 | literature_tf_master.csv の37 TF中36個のlocus_tagが誤っていた | — | 36/37 (97.3%) |
| 2 | 真のredZ (SC_RS31650/SCO5881) はT2で上昇発現し、メチル化なし | `redZ_paradox_multipanel.pdf` | LFC=+0.91, padj=3.1e-5 |
| 3 | AAGCCCGカスケードモデルのredZ経路は無効（誤った遺伝子に基づく） | — | SC_RS27300≠redZ |
| 4 | ゲノムワイド統計（4mC相関、DEG-DMG重複、モチーフ解析等）は影響なし | — | — |
| 5 | redDが55倍上昇（LFC=+5.79）でRed BGCの主要ドライバー | `balance_diagram.pdf` | LFC=+5.79 |

---

## 1. エラーの概要

### 1.1 発見の経緯

探索ループ（LOOP_PROMPT_M145.md）フェーズ3のH3仮説検証（redZパラドックスの定量的解決）において、DESeq2結果からredZ (SCO5881) の発現を直接参照したところ、包括的レポート（v2）の記載と矛盾する値が得られた。

- **v2レポートの記載**: redZ (SC_RS27300) = LFC -2.25（下方制御）
- **GFFからの正しいマッピング**: SCO5881 = **SC_RS31650**（LFC +0.91、上方制御）
- **SC_RS27300の正体**: SCO5027（winged helix DNA-binding protein、redZではない）

### 1.2 エラーの性質

`literature_tf_master.csv` において、SCO番号からSC_RS locus_tagへの変換に**系統的な誤り**が存在した。37遺伝子中36遺伝子で誤ったlocus_tagが割り当てられていた（afsR/SC_RS24295のみ正解）。

### 1.3 エラーの原因推定

SCO番号→SC_RS locus_tag変換において、GFFファイルの `old_locus_tag` フィールドを参照すべきところ、何らかの算術的オフセットまたは異なるアノテーションバージョンに基づく変換が行われた可能性が高い。

---

## 2. マッピング修正表

| 遺伝子名 | 誤locus_tag | 誤が指す遺伝子 | 正locus_tag | 正SCO |
|----------|-------------|--------------|-------------|-------|
| **redZ** | SC_RS27300 | SCO5027 | **SC_RS31650** | SCO5881 |
| **redD** | SC_RS27225 | SCO5012 | **SC_RS31630** | SCO5877 |
| **actII-ORF4** | SC_RS27570 | SCO5082 | **SC_RS27585** | SCO5085 |
| **absA1** | SC_RS02965 | SCO0203 | **SC_RS18240** | SCO3225 |
| **absA2** | SC_RS02970 | SCO0204 | **SC_RS18245** | SCO3226 |
| **afsS** | SC_RS22980 | SCO4158 | **SC_RS24290** | SCO4425 |
| **afsK** | SC_RS24290 | SCO4425 | **SC_RS24280** | SCO4423 |
| **cdaR** | SC_RS17785 | SCO3132 | **SC_RS18200** | SCO3217 |
| **cpkO/kasO** | SC_RS31605 | SCO6270 | **SC_RS33650** | SCO6280 |
| bldA | SC_RS26640 | SCO4896 | SC_RS28375 | SCO5239 |
| bldB | SC_RS25655 | SCO4700 | SC_RS30830 | SCO5723 |
| bldC | SC_RS36355 | SCO7427 | SC_RS22590 | SCO4091 |
| bldD | SC_RS25165 | SCO4602 | SC_RS09420 | SCO1489 |
| bldG | SC_RS25070 | SCO4583 | SC_RS19880 | SCO3549 |
| bldH/adpA | SC_RS14000 | SCO2378 | SC_RS16045 | SCO2792 |
| bldM | SC_RS18130 | SCO3097 | SC_RS26000 | SCO4768 |
| bldN | SC_RS26120 | SCO4793 | SC_RS18720 | SCO3323 |
| afsR | SC_RS24295 | **SCO4426（正解）** | SC_RS24295 | SCO4426 |
| crp | SC_RS10725 | SCO1673 | SC_RS19990 | SCO3571 |
| dasR | SC_RS19305 | SCO3436 | SC_RS28335 | SCO5231 |
| absB | SC_RS27825 | SCO5132 | SC_RS30075 | SCO5572 |
| afsQ1 | SC_RS33490 | SCO6472 | SC_RS26695 | SCO4907 |
| afsQ2 | SC_RS33495 | SCO6473 | SC_RS26690 | SCO4906 |
| nsdA | SC_RS35450 | SCO7156 | SC_RS30130 | SCO5582 |
| nsdB | SC_RS37265 | SCO7578 | SC_RS11445 | SCO1883 |
| wblA | SC_RS04915 | SCO0600 | SC_RS20030 | SCO3579 |
| papR2 | SC_RS18200 | SCO3217 | SC_RS35040 | SCO6569 |
| hrdB | SC_RS25850 | SCO4738 | SC_RS31340 | SCO5820 |
| hrdA | SC_RS23200 | SCO4205 | SC_RS19450 | SCO3465 |
| hrdC | SC_RS25845 | SCO4737 | SC_RS30330 | SCO5621 |
| hrdD | SC_RS18125 | SCO3096 | SC_RS18125 | SCO3202 |
| sigE | SC_RS24925 | SCO4558 | SC_RS18895 | SCO3356 |
| sigH | SC_RS25290 | SCO4627 | SC_RS28395 | SCO5243 |
| sigR | SC_RS26330 | SCO4835 | SC_RS28260 | SCO5216 |
| sigF | SC_RS04925 | SCO0600 | SC_RS22295 | SCO4035 |
| sigB | SC_RS26785 | SCO4926 | SC_RS04925 | SCO0600 |
| sigU | SC_RS04965 | SCO0608 | SC_RS16875 | SCO2954 |

---

## 3. 影響を受ける解析と受けない解析

### 3.1 影響を受けない解析（ゲノムワイド統計）

以下の解析は遺伝子IDの名前→locus_tag変換に依存せず、GFFから直接取得したSC_RS locus_tagを使用しているため、**影響なし**:

| 解析 | 理由 |
|------|------|
| DESeq2差次的発現解析（全DEG統計） | featureCounts→DESeq2はGFF locus_tag直接使用 |
| 4mC T2vsT1 Spearman相関 (r=0.125) | ゲノムワイドメチル化-発現相関 |
| 擬似相関検証（並べ替え、Bootstrap） | 上記の検証 |
| メチル化サイト検出（12,429サイト） | modkit pileupはゲノム座標ベース |
| モチーフ発見（AAGCCCG, CCGG等） | 配列モチーフベース |
| DEG-DMG重複検定（OR=1.35 for T3vsT1） | SC_RS locus_tagのセットベース |
| DMG機能エンリッチメント（enrichmentなし） | 遺伝子セットベース |
| FDR多重検定補正（67テスト中1件有意） | ゲノムワイド統計 |
| BGC発現ヒートマップ（act, red） | BGC遺伝子はゲノム座標で定義 |

### 3.2 影響を受ける解析（TF名→locus_tag変換に依存）

| 解析 | 影響の程度 | 具体的な影響 |
|------|----------|------------|
| **GRN-TFメチル化統合解析** | **致命的** | 37 TF全体のメチル化-発現協調判定が誤り |
| **AAGCCCGカスケードモデル** | **致命的** | redZ経由の経路は無効 |
| **TF BSメチル化変動解析** | **重大** | TF BSの位置が全て異なる遺伝子を指す |
| **レギュレーターネットワーク** | **重大** | 個別TFの発現値参照が誤り |
| **SARP統合解析** | **重大** | SARP TFのlocus_tag参照が誤り |
| **候補TF優先順位付け** | **重大** | 候補遺伝子リストが誤り |
| **BGC制御因子ランドスケープ** | **重大** | 27制御因子の大部分が誤った遺伝子 |

### 3.3 要修正のKey Insights

**v2包括的レポートの10 Key Insightsのうち:**

| # | Insight | 影響 | 修正要否 |
|---|---------|------|---------|
| 1 | AAGCCCG二重修飾モチーフ | 影響なし | 不要 |
| 2 | SC_RS17645 = Type I HsdM | 影響なし | 不要 |
| 3 | 4mC T2vsT1相関のrobustness | 影響なし | 不要 |
| 4 | 6mA T2vsT1相関の消失 | 影響なし | 不要 |
| 5 | FDR補正後67テスト中1件のみ有意 | 影響なし | 不要 |
| 6 | TF BS methylation depletion | **要再検証** | BS位置が異なる |
| 7 | DMG機能非選択性 | 影響なし | 不要 |
| 8 | **redZ = 唯一の協調変動TF** | **取り消し** | redZではなくSCO5027 |
| 9 | CCGKCA = BREX-2候補 | 影響なし | 不要 |
| 10 | **ゲートキーパーモデル** | **Layer 3要修正** | カスケード経路が無効 |

---

## 4. 修正対応

### 4.1 実施済み

| 対応 | ファイル | 日時 |
|------|--------|------|
| 修正済みliterature_tf_master.csv作成 | `12_grn_tf_methylation/literature_tf_master.csv` | 2026-02-24 |
| 旧ファイルアーカイブ | `literature_tf_master_ARCHIVED_260224.csv` | 2026-02-24 |
| 影響データファイルのアーカイブ | `*_ARCHIVED_260224.*` | 2026-02-24 |
| H3検証解析（真のredZ発現確認） | `26_redZ_paradox/` | 2026-02-24 |

### 4.2 今後必要な対応

| 優先度 | 対応 | 概要 |
|--------|------|------|
| **P1** | GRN-TFメチル化解析の再実行 | 修正済みmaster CSVで12_grn_tf_methylation全体を再解析 |
| **P1** | TF BS解析の再実行 | 正しいTF locus_tagでBS位置を再取得・メチル化重複を再計算 |
| **P1** | AAGCCCGカスケードモデルの再構築 | redZ経路を除去、SC_RS27300 (SCO5027) の正しい位置づけを評価 |
| **P2** | レギュレーターネットワークの再解析 | 07_regulator_network, 08_candidate_TF全体 |
| **P2** | 包括的レポートv3の作成 | 影響箇所の修正版 |
| **P3** | SARP統合解析の再確認 | 09_SARP_integration |

---

## 5. 発見の文脈：探索ループによる品質管理

このエラーは、LOOP_PROMPT_M145.md に基づく仮説駆動型探索ループの第1回実行（フェーズ3: H3仮説検証）において、DESeq2結果を直接参照することで発見された。

**教訓**:
- 遺伝子名→locus_tag変換は必ずGFFファイルの `old_locus_tag` フィールドで検証すべき
- 独立した経路（別のデータソース）から同じ結果を再現する「クロスバリデーション」が有効
- 探索ループによる仮説検証は、既存結果の品質管理としても機能する

---

## Figure解説

### redZ_paradox_multipanel.pdf

**ファイル**: `11_epigenome_integration/analysis/26_redZ_paradox/figures/redZ_paradox_multipanel.pdf`

**説明**: Red BGC制御因子の発現動態を示す3パネル図

**読み取り方**:
- Panel A: 各レギュレーターのlog2FC（T2vsT1）棒グラフ
- Panel B: Red BGC構造遺伝子の平均発現とレギュレーターバランスの対比
- Panel C: 全Red BGC遺伝子のタイムポイント別ヒートマップ

**注目点**:
- 真のredZ (SC_RS31650) はLFC=+0.91で上昇
- redD (SC_RS31630) がLFC=+5.79で圧倒的な活性化ドライバー

**サポートするInsight**: #2, #5

### balance_diagram.pdf

**ファイル**: `11_epigenome_integration/analysis/26_redZ_paradox/figures/balance_diagram.pdf`

**説明**: Activator/Repressor バランスの定量的図示

**サポートするInsight**: #5

---

*最終更新: 2026-02-24*
