# 10_SARP_motif_scan_M145_prompt.md

# 10_SARP_MOTIF_SCAN（M145, SARP 結合モチーフの in silico スキャン）

あなたは、*Streptomyces coelicolor* A3(2) M145 の RNA-seq 解析後段として、
SARP ファミリー（ActII-ORF4 など）の既知結合モチーフに基づき、
4 主要 BGC 内プロモーターおよび step 08 候補 TF 遺伝子のプロモーターにおいて
SARP 結合候補サイトを in silico で網羅的にスキャンし、
「どの SARP がどの BGC / TF を直接制御し得るか」の候補マップを作るバイオインフォマティクスエージェントです。

---

## 0. 前提・入力

以下のステップが完了している前提で動いてください:

- `COMMON_PROMPT_M145.md`
- steps 01–08（QC〜候補TF優先順位付け）
- `09_SARP_integration_M145_prompt.md`（v1）
- `09_SARP_integration_M145_prompt_v2.md`（v2 — actII-orf4 を機能的 SARP として統合済み）

### 対象ゲノム

- *Streptomyces coelicolor* A3(2) strain M145
- リファレンス: NCBI RefSeq GCF_000203835.1
- ゲノム配列・注釈は既存パスから利用（取得不要）

### 対象 BGC（4 主要クラスター）

| BGC | 産物 | CSR SARP（v2定義） |
|-----|------|-------------------|
| act | actinorhodin | actII-orf4（SC_RS27570 / SCO5082, functional_only）+ SCO5085（SC_RS27585, structural_only） |
| red | undecylprodigiosin | redD（SC_RS31630 / SCO5877, both）; redZ（SC_RS31650 / SCO5881, response regulator — SARP ではないが red CSR） |
| cda | CDA (calcium-dependent antibiotic) | SCO3217（SC_RS18200, both; 文献名 cdaR, GFF に gene_name なし） |
| cpk | coelimycin | SCO6288（SC_RS33690, both; 文献名 cpkN, GFF に gene_name なし） |

注意:
- **cpkO/kasO**（SC_RS33660 / SCO6282）は GFF で「SDR family oxidoreductase」とアノテーションされており、v1 hmmsearch でも SARP ドメインは未検出。SARP としては扱わない。
- **redZ**（SC_RS31650 / SCO5881）は GFF で「response regulator transcription factor」とアノテーションされ、SARP ではなく response regulator。ただし red クラスターの重要な CSR であるため、プロモータースキャン対象には含める。

### SARP リスト（09_v2 より）

モチーフスキャンの対象となる SARP は `SARP_list_M145_v2.tsv`（9遺伝子）を使用:

| gene_id | old_locus_tag | Subtype | BGC | Role | Category |
|---------|---------------|---------|-----|------|----------|
| SC_RS27570 | SCO5082 | small | act | CSR | functional_only |
| SC_RS27585 | SCO5085 | small | act | CSR | structural_only |
| SC_RS31630 | SCO5877 | small | red | CSR | both |
| SC_RS18200 | SCO3217 | medium | cda | CSR | both |
| SC_RS33690 | SCO6288 | small | cpk | CSR | both |
| SC_RS24295 | SCO4426 (afsR) | large | none | global | structural_only |
| SC_RS06435 | SCO0898 | small_long | none | global | structural_only |
| SC_RS22750 | SCO4116 | small_long | none | global | structural_only |
| SC_RS13320 | SCO2259 | small_long | none | global | structural_only |

※ AfsR（SCO4426）は **large SARP**（DBD+BTAD+NB-ARC+TPR）であり、LAL 型ではない。M145 には SARP-LAL 型は存在しない（v1 で確認済み）。

### Step 08 候補 TF（9遺伝子）

`TF_candidates_experimental_priority.tsv` より:

| # | gene_id | old_locus_tag | product |
|---|---------|---------------|---------|
| 1 | SC_RS27570 | SCO5082 | TetR family TR ActII (= actII-orf4) |
| 2 | SC_RS37190 | SCO6993 | LuxR family TR AbsR2 |
| 3 | SC_RS36895 | SCO6937 | LuxR C-terminal-related TR |
| 4 | SC_RS31650 | SCO5881 | response regulator TF (= redZ in literature) |
| 5 | SC_RS21215 | SCO3818 | response regulator |
| 6 | SC_RS29775 | SCO5518 | PucR family TR |
| 7 | SC_RS28860 | SCO5337 | XRE family TR |
| 8 | SC_RS26000 | SCO4768 | response regulator TF |
| 9 | SC_RS16870 | SCO2953 | anti-sigma U factor RsuA |

---

## 1. パスと変数定義

```bash
MOTIF_ROOT="/Users/okaban/bioinfo/rna-seq/10_SARP_motif_scan"
SARP_ROOT="/Users/okaban/bioinfo/rna-seq/09_SARP_integration"
ANNOT_ROOT="/Users/okaban/bioinfo/rna-seq/05_annotation"
TF_ROOT="/Users/okaban/bioinfo/rna-seq/08_candidate_TF_prioritization"

ANNOT_RUN_DIR="${ANNOT_ROOT}/analysis/05_annotation_260128_v1"
TF_RUN_DIR="${TF_ROOT}/analysis/08_candidate_TF_260128_v1"
SARP_V2_RUN_DIR="${SARP_ROOT}/analysis/09_SARP_integration_260128_v2"

# リファレンス
REF_DIR="/Users/okaban/SFC-CNS Dropbox/Takeda Tomoki/my_projects/2026/M145_RNA-seq/data/RH250715199/02_M145-ref/ncbi_dataset/data/GCF_000203835.1"
GENOME_FNA="${REF_DIR}/GCF_000203835.1_ASM20383v1_genomic.fna"
GENOME_GFF="${REF_DIR}/genomic.gff"

# 入力テーブル
GENE_MASTER_FILE="${ANNOT_RUN_DIR}/tables/gene_master_with_BGC_regulators.tsv"
BGC_DEF_FILE="${ANNOT_RUN_DIR}/tables/BGC_definition_manual.tsv"
SARP_V2_LIST="${SARP_V2_RUN_DIR}/tables/SARP_list_M145_v2.tsv"
SARP_V2_ACTIVITY="${SARP_V2_RUN_DIR}/tables/SARP_BGC_activity_summary_v2.tsv"
TF_CANDIDATES="${TF_RUN_DIR}/tables/TF_candidates_experimental_priority.tsv"

# 出力
MOTIF_OUT_ROOT="${MOTIF_ROOT}/analysis"
MOTIF_SCRIPT_DIR="${MOTIF_ROOT}/scripts"
MOTIF_RUN_DATE=$(date +%y%m%d)
MOTIF_RUN_ID="10_SARP_motif_scan_${MOTIF_RUN_DATE}_v1"
MOTIF_RUN_DIR="${MOTIF_OUT_ROOT}/${MOTIF_RUN_ID}"
```

ディレクトリ:

```text
${MOTIF_RUN_DIR}/
├── tables/
├── figures/
├── logs/
├── pwm/
├── promoters/
└── SARP_motif_scan_report_M145.md
```

---

## 2. 環境・ログセットアップ

```bash
conda activate rna-seq

# MEME Suite（FIMO）がインストールされていない場合:
# conda install -y -c bioconda meme
# ※ FIMO のバージョンを確認: fimo --version

mkdir -p "${MOTIF_RUN_DIR}/tables" \
         "${MOTIF_RUN_DIR}/figures" \
         "${MOTIF_RUN_DIR}/logs" \
         "${MOTIF_RUN_DIR}/pwm" \
         "${MOTIF_RUN_DIR}/promoters" \
         "${MOTIF_SCRIPT_DIR}"

LOG_FILE="${MOTIF_RUN_DIR}/logs/SARP_motif_scan_pipeline.log"
exec > >(tee -a "${LOG_FILE}") 2>&1

echo "=== 10_SARP_motif_scan started at $(date) ==="
echo "Using SARP v2 list: ${SARP_V2_LIST}"
echo "Using TF candidates: ${TF_CANDIDATES}"
```

---

## 3. 参考モチーフ情報の整理

SARP ファミリーの既知 Motif 情報を統合して PWM を構築する。

### 3.1 ActII-ORF4（act クラスター）

- DNase I footprint 解析により、*actVI-ORF1–ORFA* 間のプロモーター −35 領域に、5'-TCGAG-3' を含む保護領域が特定されている（Arias et al., 1999）。
- ActII-ORF4 は、複数の act プロモーター（例: *actVI-ORF1* など）の −35 に結合して転写を活性化する。
- 既報のコンセンサスはクラシカルな 11 bp SARP 配列（1 回転分、heptamer + 4 bp spacer）に対応するとされる。

→ 手順: Arias 1999（ActII-ORF4）で報告されている結合領域のシーケンスを抽出し、PWM を構築。

### 3.2 一般 SARP モチーフ（heptameric direct repeats）

- SARP-LAL 型 PimR では、CGGCAAG の 7-mer 直列リピート（4 bp spacer）を 2–3 コピー持つオペレーターが主要結合配列である。
- 他の SARP-LAL（SanG, PolR など）でも同様の heptamer DR が報告されており、「7-mer + 4 nt spacer（計 11 bp）」がクラシカルな SARP 結合単位と総括されている。
- 総説でも、SARP は標的プロモーター −35 近傍に、2–3 つの heptamer DR を認識して結合するのが一般像とされる。

→ 手順:
- PimR 型の heptamer（例: CGGCAAG）を seed とし、報告されている他 SARP（SanG, PolR 等）のオペレーター配列から「general SARP heptamer-PWM」を作成。
- 「厳しめ（高スコア）」と「ゆるめ（やや保存性低）」の 2 種類の PWM を用意して FIMO を二段階でかける。

### 3.3 SCO3217（cdaR）/ SCO6288（cpkN）

- SCO3217 は CDA クラスター内の複数プロモーターに結合する SARP で、クラスター内の複数 TSS が SCO3217/cdaR 依存であることが示されているが、明確な短いコンセンサスは確立されていない。
- SCO6288（cpkN）は CPK クラスターの SARP だが、結合部位は同定されておらず、「SARP であること」「cpk 遺伝子発現に必須であること」のみ明らかである。

→ 手順: これらについては個別の PWM を作るのではなく、「一般 SARP heptamer PWM」で cda/cpk クラスター内プロモーターをスキャンし、候補オペレーターを列挙する。

---

## 4. 解析ステップ

### ステップ 1: ゲノム・注釈の準備と座標整理

1. 既存の M145 ゲノム配列（`${GENOME_FNA}`）と注釈（`${GENOME_GFF}`）を使用（取得不要）。
2. `${BGC_DEF_FILE}` と `${GENE_MASTER_FILE}` から、4 BGC（act/red/cda/cpk）の座標と、v2 SARP 9遺伝子 + step 08 候補 TF 9遺伝子の座標を整理。

出力:
- `tables/BGC_and_TF_coordinates.tsv`（gene_id, old_locus_tag, gene_name, bgc_name, start, end, strand, gene_type [BGC_gene / SARP / TF_candidate]）

### ステップ 2: プロモーター領域の抽出

対象:
- 各 BGC 内の全遺伝子の開始位置の上流 300–500 bp
- v2 SARP 9遺伝子自身の上流 300–500 bp
- step 08 候補 TF 9遺伝子の上流 300–500 bp

手順:
1. `BGC_and_TF_coordinates.tsv` を元に、+strand なら start − 500〜start − 1、−strand なら end + 1〜end + 500 をプロモーター候補領域として抽出（ゲノム端の処理は適宜）。
2. 各領域を FASTA に書き出す（1 gene 1 sequence）。

出力:
- `promoters/promoters_BGCs_500bp.fasta`
- `promoters/promoters_SARPs_500bp.fasta`
- `promoters/promoters_TFs_500bp.fasta`

### ステップ 3: SARP PWM の構築

1. ActII-ORF4 結合領域（複数プロモーターの footprint データ）を ActII-ORF4 PWM の元データとする（Arias et al., 1999）。
2. PimR / SanG / PolR などの heptamer DR を整理し、「general SARP heptamer PWM」を作る。
3. 必要なら heptamer + spacer + heptamer の「複合 PWM」も定義（ただし初回は heptamer 単独モチーフ中心でもよい）。

出力:
- `pwm/PWM_ActII_ORF4.meme`
- `pwm/PWM_SARP_heptamer_strict.meme`
- `pwm/PWM_SARP_heptamer_relaxed.meme`

### ステップ 4: FIMO によるモチーフスキャン

1. BGC プロモーター FASTA を対象に FIMO 実行。
   - 入力 PWM: ActII-ORF4 + SARP heptamer（strict/relaxed）
   - 出力: モチーフ名, seq_id, start, end, strand, score, p-value, matched sequence を含む TSV
2. SARP プロモーター FASTA に対しても同様に実行。
3. 候補 TF プロモーター FASTA に対しても同様に実行。
4. p-value や q-value で閾値（例: q < 0.05）を設定し、「高信頼ヒット」にフラグ付け。

出力:
- `tables/fimo_BGC_promoters.tsv`
- `tables/fimo_SARP_promoters.tsv`
- `tables/fimo_TF_promoters.tsv`

---

## 5. 集約と解釈の方針

### 5.1 BGC 内プロモーターでの SARP モチーフ分布

- 各 BGC について:
  - 高スコア SARP ヒット数 / プロモーター数
  - プロモーター位置に対するモチーフ位置（TSS 推定位置からの距離）分布
- 期待:
  - act クラスター内では ActII-ORF4 PWM による結合サイトが、既報の *actVI-ORF1* など周辺で強く出るはず
  - cda/cpk では、未同定の SARP 結合サイト候補がいくつか浮上する可能性

→ 出力: `tables/BGC_SARP_motif_summary.tsv`（BGC, gene_id, old_locus_tag, motif_type, motif_count, top_score, mean_distance_from_start）

### 5.2 候補 TF プロモーターでの SARP モチーフ

- step 08 候補 TF のうち、「プロモーターに SARP ヒットが多いもの」を抽出
- これにより:
  - SARP → TF → BGC という階層制御の可能性を示唆
  - 逆に SARP ヒットがほとんどない TF は、SARP 直結とは考えにくい

→ 出力: `tables/TF_promoters_SARP_hits.tsv`（gene_id, old_locus_tag, product, motif_type, hit_count, best_score）

### 5.3 SARP ↔ BGC/TF 関係のネットワーク草案

次のような三者関係を表にする:

| SARP (gene_id) | 標的候補 (gene_id) | 種別 | SARP ヒット数 | best_score | コメント |
|-----------------|-------------------|------|-------------|-----------|---------|
| SC_RS27570 (actII-orf4) | actVI-ORF1 | BGC gene | n | x.xx | 既知ターゲット |
| SC_RS27570 (actII-orf4) | SC_RS21215 (SCO3818) | TF | n | x.xx | プロモーターに SARP motif |
| SC_RS18200 (SCO3217/cdaR) | cdaA etc. | BGC gene | n | x.xx | CDA 内の主要プロモーター |
| SC_RS33690 (SCO6288/cpkN) | cpk genes | BGC gene | n | x.xx | cpk 内の新規候補サイト |

（具体的な標的遺伝子名は解析結果に依存）

→ 出力: `tables/SARP_target_network.tsv`

---

## 6. レポートに含めるべき内容

最終レポート `SARP_motif_scan_report_M145.md` では少なくとも以下を盛り込む:

1. **既知知見との整合性**
   - ActII-ORF4 の既知結合部位（*actVI-ORF1* プロモーターなど；Arias et al., 1999）が、PWM/FIMO でも高スコアで再現されるか
   - CDA/CPK における SCO3217/SCO6288 の役割との整合性

2. **新規候補**
   - 4 BGC 内で、従来明確に報告されていない SARP モチーフが見つかるプロモーター
   - 候補 TF プロモーターにおける SARP モチーフの有無（「SARP の上にさらに SARP/TF が乗っている」可能性）

3. **v2 二層分類との統合**
   - structural SARP（ドメイン検出あり）vs functional_only（actII-orf4）のモチーフパターンの違い
   - SARP_category ごとのヒット分布の違い

4. **今後の wet 実験への接続**
   - EMSA/DNase I footprint を行うべきプロモーター候補
   - プロモーター置換・点変異実験（heptamer DR を壊す）に進める候補部位

---

## 7. 成果物チェックリスト

### ファイル・テーブル類

- [ ] `${MOTIF_RUN_DIR}/logs/SARP_motif_scan_pipeline.log`
- [ ] `${MOTIF_RUN_DIR}/tables/BGC_and_TF_coordinates.tsv`
  - [ ] act/red/cda/cpk クラスタの全遺伝子座標
  - [ ] v2 SARP 9遺伝子の座標
  - [ ] step 08 候補 TF 9遺伝子の座標
- [ ] `${MOTIF_RUN_DIR}/promoters/promoters_BGCs_500bp.fasta`
- [ ] `${MOTIF_RUN_DIR}/promoters/promoters_SARPs_500bp.fasta`
- [ ] `${MOTIF_RUN_DIR}/promoters/promoters_TFs_500bp.fasta`
- [ ] `${MOTIF_RUN_DIR}/pwm/PWM_ActII_ORF4.meme`
- [ ] `${MOTIF_RUN_DIR}/pwm/PWM_SARP_heptamer_strict.meme`
- [ ] `${MOTIF_RUN_DIR}/pwm/PWM_SARP_heptamer_relaxed.meme`
- [ ] `${MOTIF_RUN_DIR}/tables/fimo_BGC_promoters.tsv`
- [ ] `${MOTIF_RUN_DIR}/tables/fimo_SARP_promoters.tsv`
- [ ] `${MOTIF_RUN_DIR}/tables/fimo_TF_promoters.tsv`
- [ ] `${MOTIF_RUN_DIR}/tables/BGC_SARP_motif_summary.tsv`
- [ ] `${MOTIF_RUN_DIR}/tables/TF_promoters_SARP_hits.tsv`
- [ ] `${MOTIF_RUN_DIR}/tables/SARP_target_network.tsv`

### 図・可視化（任意だが推奨）

- [ ] BGC ごとの線形スキーマ図（各 gene 配置 + プロモーター位置 + SARP motif ヒット位置）
- [ ] ヒートマップ（横軸: BGC/TF、縦軸: motif 種類、値: ヒット数 or −log10(q)）

### レポート

- [ ] `${MOTIF_RUN_DIR}/SARP_motif_scan_report_M145.md`
  - [ ] ActII-ORF4 の既知結合部位が PWM で再現されることの確認
  - [ ] BGC 内の新規 SARP モチーフ候補リスト
  - [ ] SARP 直結が疑われる候補 TF のリスト
  - [ ] 追試（EMSA / footprint / プロモーター変異）候補の優先順位リスト
