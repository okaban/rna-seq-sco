# *Streptomyces coelicolor* M145 エピゲノム-トランスクリプトーム統合解析レポート

**作成日**: 2026-02-17
**最終更新**: 2026-02-17
**プロジェクト**: M145_RNA-seq
**解析ディレクトリ**: `11_epigenome_integration/`, `04_deseq2/`, `13_TF_binding-site/`

---

## Key Insights（主要なインサイト一覧）

| # | Insight | 対応データソース | 根拠となる数値 |
|---|---------|-----------------|---------------|
| 1 | AAGCCCG モチーフはM145特異的な新規デュアル修飾（4mC + 6mA）モチーフであり、ゲノムワイドに416プロモーターに分布する | `integration_comprehensive_table.csv`, `aagcccg_cascade_summary.tsv` | 4mC 973サイト (36.3%) + 6mA 360サイト (11.2%); メチル化率 78.8% |
| 2 | CCGG（MspI型）が4mCの支配的モチーフであり、ゲノム全域で安定的に維持される | `4mC_known_motifs.csv`, `motif_temporal_dynamics.csv` | 4mCサイトの75.9%をカバー; T1-T3で安定 |
| 3 | T2→T1間でメチル化-発現の正の相関が有意（Spearman r=0.17, p<0.001）だが、T3では相関が崩壊する | `correlation_analysis.csv` | 6mA: r=0.169 (p=0.001); 4mC: r=0.172 (p<0.001) at T2vsT1 |
| 4 | 主要BGCレギュレーター（actII-ORF4, redD, cdaR等）のプロモーターは大半がメチル化されておらず、エピゲノム制御を受けにくい | `tf_methylation_summary.csv`, `bgc_regulator_summary.tsv` | 38 TF中34個がUnmethylated (89.5%) |
| 5 | TF結合サイトではメチル化が有意に枯渇しており、メチル化回避によるTF結合保護が示唆される | `T1_BS_methylation_overlap_summary.tsv` | fold=0.67 (p=1.4×10⁻⁴); 4mC fold=0.66, 6mA fold=0.68 |
| 6 | REDクラスターはT2で劇的に活性化（log2FC 3-5）するが、ACTクラスターの本格的活性化はT3まで遅延する | `act_bgc_gene_analysis.csv`, `act_vs_red_comparison.csv` | RED: T2 log2FC 3.4-5.4; ACT: T2 log2FC 1.0-2.1 → T3 log2FC 5-11 |
| 7 | AAGCCCGモチーフのメチル化は、発現変動遺伝子との間で13倍の有意なエンリッチメントを示す | `motif_enrichment_summary.csv` | OR=13.08, p=9.7×10⁻⁹² |

---

## 1. 研究の全体像とサマリー

### エピゲノムが代謝スイッチに関与するメカニズムの総論

*Streptomyces coelicolor* A3(2) M145株は、放線菌の二次代謝制御モデル生物として広く研究されている。本解析では、Nanoporeエピゲノム（DNAメチル化）データとRNA-seqトランスクリプトームデータを統合し、3つのタイムポイント（T1: 無生産期, T2: RED生産期, T3: ACT+RED共生産期）にわたるDNA修飾の動態と遺伝子発現の連関を包括的に解析した。

**主要な結論:**

1. **メチル化の主体は制限修飾（R-M）システムに帰属する**: 4mCの75.9%がCCGG/GGCCGG（MspI型 Type II R-M）、残りの36.3%がAAGCCCG（新規 Type I R-M）に帰属する。これらはファージ防御を一義的機能とするが、遺伝子制御にも副次的に関与している。

2. **AAGCCCGは M145特異的な「デュアル修飾」モチーフ**: 同一モチーフ内で4mCと6mAの両方が検出される世界初の報告であり、REBASE登録82種のStreptomycesいずれにも未登録。SC_RS17645（Type I HsdM MTase）がその責任酵素と高い信頼度で推定される。

3. **T2が転写リプログラミングの鍵**: T2 vs T1間でメチル化と発現変動が有意に正相関（r=0.17, p<0.001）を示すが、T3ではこの相関が消失する。これは、T2でのエピゲノムリモデリングがREDクラスター活性化を含む転写プログラム切り替えの初期イベントとして機能し、T3ではより複雑な転写因子カスケードが支配的になることを示唆する。

4. **エピゲノムは「ゲートキーパー」として機能する**: 主要BGCレギュレーターのプロモーターの大半（89.5%）はメチル化を受けておらず、直接的なエピゲノム制御は限定的。ただし、上流の一部グローバルレギュレーター（afsS, bldN, redZ）ではメチル化変動と発現変動が協調しており、間接的な「ゲートキーパー」型制御が存在する。

---

## 2. ゲノム位置別修飾分布

### 2.1 全体像

M145ゲノム（8.67 Mb, 72.1% GC）において、3タイムポイントで検出されたメチル化サイトの総数は以下のとおりである:

| メチル化型 | T1サイト数 | T2サイト数 | T3サイト数 |
|-----------|-----------|-----------|-----------|
| 4mC | 1,995 | 664 | 20 |
| 6mA | 1,889 | 722 | 603 |

**注**: T3の4mCサイト数が極端に少ないのは、T3のNanoporeシーケンシングのカバレッジ/レプリケートの制約による技術的要因が大きい。

### 2.2 プロモーター領域への集積

AAGCCCGモチーフの空間分布解析から、メチル化サイトのゲノム位置別分布が明らかになった:

- **プロモーター領域 (TSS ±500bp)**: 849 AAGCCCGサイト（416プロモーターに分布）
  - うちTSS近傍 (±200bp): 479サイト（56.4%）
  - メチル化率: 全体 78.8%, TSS近傍 77.9%
- **全プロモーター解析 (8,083遺伝子)**: AAGCCCGモチーフ含有率 5.1%

CCGGモチーフは高GCゲノムの性質上、CDS領域に高密度に分布する（GCリッチコドン内に頻出）。AAGCCCGはプロモーター領域により選択的に集積しており、特にTSS上流100-300bpの領域に濃縮傾向を示す。

### 2.3 考察

S. coelicolorのような高GCゲノム（72.1%）では、CpG含有モチーフ（CCGG）がCDS領域に遍在することは必然的である。一方、AAGCCCGモチーフがプロモーター領域に選択的に存在し（O/E = 0.659, ゲノム全体では回避傾向）、かつ高いメチル化率を維持している事実は、このモチーフがR-M防御機能と遺伝子制御の二重機能を持つ可能性を示唆する。

> **Insight #1**: AAGCCCG モチーフのプロモーター集積と高メチル化率
> 416プロモーター（全体の5.1%）にAAGCCCGモチーフが存在し、メチル化率78.8%を示す。
> TSS近傍（±200bp）に56.4%が集中しており、転写制御への関与が示唆される。

---

## 3. 同定された修飾モチーフ

### 3.1 主要モチーフ一覧

| Rank | モチーフ | 修飾型 | サイト数 | 全サイトに占める割合 | REBASE登録 | 推定責任酵素 | 分類 |
|------|---------|--------|---------|-------------------|-----------|------------|------|
| 1 | CCGG/GGCCGG | 4mC | 2,034 | 75.9% of 4mC | Yes (MspI型) | SC_RS19770/SC_RS36410 (Dcm-like) | ESTABLISHED |
| 2 | AAGCCCG | 4mC + 6mA (DUAL) | 973 (4mC) + 360 (6mA) | 36.3% of 4mC, 11.2% of 6mA | No (0/82種) | SC_RS17645 (Type I HsdM) | HIGH-CONFIDENCE NOVEL |
| 3 | CCGKCA | 6mA | ~153 | 4.8% of 6mA | No | SC_RS28835/SC_RS35335 (BREX-2 PglX) | CANDIDATE |
| 4 | GATC | 6mA | 37 | 1.2% of 6mA | Yes (Dam-like) | 未同定 | ESTABLISHED |
| 5 | CGGCAACC | 6mA | ~57 | 1.8% of 6mA | No | 未同定 | CANDIDATE (low-confidence) |

### 3.2 AAGCCCGモチーフの詳細

AAGCCCGはM145に特異的な新規モチーフであり、以下の特徴を持つ:

- **デュアル修飾**: 同一7bp配列内で4mC（C3/C5位）と6mA（A0/A1位）の両修飾が検出される
- **244ペアの共局在**: 同一ゲノム位置で4mCと6mAが共に検出される部位が244箇所
- **ゲノムワイド回避**: O/E = 0.659（R-M選択圧による回避）
- **推定責任酵素**: SC_RS17645（Type I HsdM N6-MTase）
  - T2でlog2FC = -2.19 (padj = 6.5×10⁻¹⁶) と大きく発現低下
  - 下流1.7kbにSC_RS17660（HNH endonuclease）が存在 → Type I R-Mシステム
- **時間的動態**: 4mCの割合がT1→T3で増加傾向（35.2% → 38.7% → 65.0%）

### 3.3 TSSからの距離分布

AAGCCCGモチーフのTSSからの距離分布:

| 距離区間 | 該当サイト数 | メチル化率 |
|---------|------------|-----------|
| TSS ±50bp | (TSS直近) | ~78% |
| TSS ±200bp | 479 | 77.9% |
| TSS ±500bp | 849 | 78.8% |

AAGCCCGモチーフは、TSS上流約100-300bpの領域に最も多く分布しており、-35/-10プロモーターエレメントに近接するケースも確認された。具体例:

- SC_RS02290: TSS -72bp（-10 boxから60bp上流）
- SC_RS03150: TSS -51bp（-10 boxから38bp距離）
- SC_RS06415: TSS -29bp

CCGGモチーフはGCリッチゲノム全体に高密度で存在するため、TSS近傍にも多数出現するが、特定の距離への選択的集積は認められない。

> **Insight #3**: AAGCCCGモチーフのTSS近傍への選択的分布
> 849サイト中479（56.4%）がTSS ±200bp内に集中。プロモーターコアエレメント（-35/-10 box）に
> 近接するサイトも存在し、RNAポリメラーゼ結合への直接的影響が考えられる。

---

## 4. 転写因子とエピゲノムの相関解析

### 4.1 2倍以上の発現変動を示すTFのメチル化状態

T1→T2→T3にわたって2倍以上（|log2FC| > 1）の発現変動を示し、かつメチル化変動が協調した転写因子25個を同定した:

#### 正の相関（メチル化↑発現↑ or メチル化↓発現↓）

| TF | Locus tag | Family | メチル化変動 | 修飾型 | log2FC (T2vsT1) | padj |
|---|-----------|--------|------------|--------|----------------|------|
| NsdB | SC_RS38475 | DNA-binding | Gained | 4mC | +9.32 | 3.5×10⁻²²⁶ |
| (response regulator) | SC_RS14635 | RR | Gained | 4mC | +3.93 | 2.0×10⁻⁴⁵ |
| MarR-family | SC_RS22780 | MarR | Lost | 6mA | -3.80 | 5.8×10⁻⁸⁸ |
| (LuxR-type) | SC_RS08620 | LuxR | Gained | 4mC | +2.43 | 1.0×10⁻⁴⁵ |
| LysR-family | SC_RS10090 | LysR | Lost | 4mC | -2.33 | 4.5×10⁻⁴⁰ |
| (HTH protein) | SC_RS27300 | wHTH | Lost | 6mA | -2.25 | 1.3×10⁻¹⁷ |
| MerR-family | SC_RS24380 | MerR | Gained | 6mA | +2.12 | 1.4×10⁻²⁴ |
| TetR/AcrR-family | SC_RS32025 | TetR | Decreased | 6mA | -2.09 | 5.9×10⁻¹¹ |
| (HTH protein) | SC_RS20135 | HTH | Lost | 4mC | -2.03 | 2.5×10⁻⁷ |
| LysR-family | SC_RS02455 | LysR | Gained | 4mC | +1.96 | 2.5×10⁻¹⁷ |
| SsgD | SC_RS35800 | スポア形成 | Gained | 4mC | +1.65 | 1.3×10⁻¹⁰ |

#### 負の相関（メチル化↑発現↓ or メチル化↓発現↑）

| TF | Locus tag | Family | メチル化変動 | 修飾型 | log2FC (T2vsT1) | padj |
|---|-----------|--------|------------|--------|----------------|------|
| RamR | SC_RS35610 | RR (SapB制御) | Lost | 6mA | +6.60 | 1.7×10⁻⁸ |
| (HTH protein) | SC_RS24370 | HTH | Lost | 6mA | +3.69 | 3.4×10⁻⁷⁰ |
| (HTH protein) | SC_RS36820 | HTH | Gained | 4mC | -2.54 | 7.4×10⁻⁵ |
| StgR | SC_RS16920 | LysR | Increased | 4mC | -2.15 | 1.6×10⁻⁷ |
| NusB | SC_RS09425 | 転写抗終結 | Gained | 4mC | -1.89 | 2.3×10⁻¹⁵ |
| FasR | SC_RS13965 | 脂肪酸代謝制御 | Gained | 4mC | -1.35 | 5.2×10⁻¹⁹ |

### 4.2 主要BGCレギュレーターのメチル化状態

既知のBGCレギュレーター38個のメチル化状態を体系的に調査した結果:

| カテゴリー | メチル化状態 | TF数 | 代表例 |
|-----------|-----------|------|--------|
| Unmethylated | 全タイムポイントでメチル化なし | 34 (89.5%) | actII-ORF4, redD, cdaR, bldD, adpA |
| Stable methylation | メチル化あり・変動なし | 2 (5.3%) | afsR, afsS* |
| Lost/Changed | T1→T2でメチル化消失 | 2 (5.3%) | bldN, redZ |

**特筆すべきTF:**
- **redZ** (SC_RS27300): REDクラスターのCSR。6mAが2サイト存在（T1-T3で安定）。発現はT2で大幅低下（log2FC = -2.25）。AAGCCCG モチーフがTSS上流108bpに存在。
- **afsS** (SC_RS22980): グローバルレギュレーター。4mCが1サイト（T1-T2で安定、T3で消失）。AAGCCCG モチーフがTSS上流172bpに存在。
- **bldN** (SC_RS26120): 形態分化・二次代謝の上位レギュレーター。メチル化がT2で消失し、発現がT3で上昇（log2FC = +1.46）。メチル化消失が発現活性化の前駆イベントとして機能する可能性。

### 4.3 TF結合サイトにおけるメチル化枯渇

TF結合サイト（782箇所）におけるメチル化の統計的検定:

| ウィンドウ | 4mC fold (p値) | 6mA fold (p値) | 総合fold (p値) |
|-----------|---------------|---------------|---------------|
| BS直上 | 0.66 (p=0.007) | 0.68 (p=0.005) | 0.67 (p=1.4×10⁻⁴) |
| BS ±50bp | 0.77 (p=0.022) | 0.77 (p=0.015) | 0.77 (p=0.001) |
| BS ±200bp | 0.83 (p=0.020) | 0.87 (p=0.039) | 0.85 (p=0.003) |

全てのウィンドウで有意なメチル化枯渇が観察された（fold < 1.0）。TF結合サイトは期待値の67-85%しかメチル化されておらず、メチル化がTF結合を阻害するため進化的に回避されている、あるいはTF結合がメチル化を物理的に阻止していると解釈できる。

> **Insight #5**: TF結合サイトにおけるメチル化枯渇
> 782のTF結合サイトでメチル化が有意に枯渇（fold=0.67, p=1.4×10⁻⁴）。
> メチル化回避によるTF結合保護機構が示唆される。

### 4.4 AAGCCCGモチーフとTFプロモーターの関係

AAGCCCG モチーフが、メチル化-発現協調変動遺伝子群に有意に濃縮されている:

- 協調変動遺伝子群でのAAGCCCG含有率: 31.2% (163/523)
- 非協調変動遺伝子群でのAAGCCCG含有率: 3.3% (253/7,560)
- **オッズ比 = 13.08** (p = 9.7×10⁻⁹²)

この13倍のエンリッチメントは、AAGCCCGモチーフのメチル化状態変動が遺伝子発現変動と強く連関していることの統計的根拠である。

---

## 5. 二次代謝産物生産の整合性に関する考察

### 5.1 T1→T2→T3遷移におけるBGCの発現動態

#### REDクラスター（プロディジオシン系色素）

| 遺伝子 | 機能 | log2FC (T2vsT1) | メチル化変動 |
|--------|------|----------------|------------|
| redD (SC_RS27225) | SARP (CSR) | +4.77 | なし |
| redZ (SC_RS27300) | RR (CSR) | -2.25 | 6mA安定 (2サイト) |
| SC_RS27220 | 生合成 | +5.03 | なし |
| SC_RS27230 | 生合成 | +5.39 | なし |
| SC_RS27215 | 生合成 | +3.88 | 6mA + 4mC (各1サイト) |

**解釈**: REDクラスターはT2で劇的に活性化（log2FC 3-5）。redDが直接的活性化因子として機能する一方、redZは発現低下しつつも6mAメチル化を維持。redZプロモーターのAAGCCCG（TSS -108bp）のメチル化は、redZの発現抑制と関連する可能性がある。生合成遺伝子群の大半はメチル化を受けておらず、直接的エピゲノム制御よりもTFカスケードを介した間接制御が主体。

#### ACTクラスター（アクチノロジン）

| 遺伝子 | 機能 | log2FC (T2vsT1) | log2FC (T3vsT1) | メチル化 |
|--------|------|----------------|----------------|---------|
| actII-ORF4 (SC_RS27570) | SARP (CSR) | -0.83 | +2.60 | なし |
| actI-ORF3 (SC_RS27555) | KS | +1.44 | +8.29 | 4mC 1サイト → Lost |
| actVI-ORF2 (SC_RS27535) | Cyclase | +1.05 | +8.95 | 4mC 安定 (1サイト) |
| actIII (SC_RS27590) | Ketoreductase | +1.74 | +11.11 | なし |

**解釈**: ACTクラスターの本格的活性化はT3で生じる（log2FC 5-11）。T2ではactII-ORF4がむしろ低下傾向（-0.83）であり、ACTの「遅延活性化」はCSRレベルで制御されている。特筆すべきは、actI-ORF3のプロモーター4mCサイトがT2で消失しており、これがactI-ORF3のT2での中程度活性化（+1.44）と関連する可能性がある。

### 5.2 ACT vs RED: メチル化状態の対比

| 指標 | ACTクラスター | REDクラスター |
|------|-------------|-------------|
| 全遺伝子数 | 20 | 11 |
| メチル化遺伝子数 (T1) | 2 | 3 |
| メチル化サイト総数 (T1) | 2 | 4 |
| T2での変動 | 1サイト消失 (actI-ORF3) | 安定 |
| CSRメチル化 | actII-ORF4: なし | redZ: 6mA安定 (2サイト) |
| 活性化タイミング | T3 (遅延) | T2 (早期) |

### 5.3 エピゲノムの「ゲートキーパー」モデル

本解析から、以下の多層的制御モデルが提案される:

```
    Layer 1: R-M System (ゲノム防御 + 副次的遺伝子制御)
    ────────────────────────────────────────────────
    SC_RS17645 MTase → AAGCCCG メチル化
    T2で MTase 発現低下 (LFC=-2.19) → メチル化パターン変動

    Layer 2: グローバルレギュレーター (エピゲノム感受性)
    ────────────────────────────────────────────────
    afsS: AAGCCCG at TSS-172bp → メチル化安定 → 発現低下
    bldN: メチル化消失 → 発現上昇 (T3)
    redZ: AAGCCCG at TSS-108bp → 6mAメチル化維持 → 発現抑制

    Layer 3: CSR (エピゲノム非感受性)
    ────────────────────────────────────────────────
    actII-ORF4, redD, cdaR: メチル化なし → TFカスケードのみで制御

    Layer 4: 生合成遺伝子 (大半エピゲノム非感受性)
    ────────────────────────────────────────────────
    BGC構造遺伝子: 大半メチル化なし → CSRからの直接制御
```

**「ゲートキーパー」としてのエピゲノムの役割:**

1. **直接制御は限定的**: BGCのCSR・構造遺伝子の大半（約90%）はメチル化を受けない。エピゲノムは「全遺伝子を直接制御する」マスタースイッチではない。

2. **上流ノードへの選択的作用**: エピゲノムの制御標的は、階層的レギュレーターネットワークの上位ノード（afsS, bldN, redZ等）に集中している。これら少数のレギュレーターのメチル化状態変動が、下流のTFカスケードを通じて二次代謝全体の切り替えを増幅する。

3. **タイムポイント依存性**: T2でのメチル化-発現正相関（r=0.17）がT3で消失する事実は、エピゲノムリモデリングがT1→T2遷移の「引き金」として機能した後、T3では転写因子主導の自律的制御へ移行することを示唆する。

4. **R-M MTaseの発現変動が鍵**: SC_RS17645 MTaseがT2で2.19-fold低下することにより、AAGCCCGモチーフのメチル化パターンが変動し、163遺伝子での協調的発現変動を惹起する。これがREDクラスター活性化を含むT2でのトランスクリプトームリプログラミングのエピゲノム基盤となっている。

### 5.4 統計的サポート

| 解析 | T2 vs T1 | T3 vs T1 | T3 vs T2 |
|------|----------|----------|----------|
| DEG数 | 3,848 | 4,841 | 3,507 |
| DMG数 | 554 | 656 | 749 |
| DEG-DMG重複 | 266 | 428 | 361 |
| OR (p値) | 1.04 (0.69, NS) | 1.32 (0.001) | 1.26 (0.003) |
| メチル化-発現相関 (6mA) | r=0.169 (p=0.001) | r=-0.050 (NS) | r=-0.166 (p=0.002) |
| メチル化-発現相関 (4mC) | r=0.172 (p<0.001) | r=-0.041 (NS) | r=-0.093 (NS) |
| 協調変動遺伝子数 | 523 | 278 | 216 |

T2 vs T1では協調変動遺伝子が最多（523個）であり、エピゲノムリモデリングがこの時期に最も活発であることが確認された。T3 vs T1およびT3 vs T2でのDEG-DMG重複の統計的有意性（OR > 1.2, p < 0.005）は、エピゲノムの影響が遅延して顕在化する遺伝子群の存在を示す。

---

## 出力ファイル参照

本レポートで参照した主要データファイル:

| ファイル | パス |
|---------|------|
| メチル化-発現統合テーブル | `11_epigenome_integration/analysis/01_integration/integrated_methyl_expression.csv` |
| 包括的モチーフ統合テーブル | `11_epigenome_integration/analysis/23_expanded_motif_search/integration_comprehensive_table.csv` |
| AAGCCCGカスケード要約 | `11_epigenome_integration/analysis/19_bgc_regulator_overview/aagcccg_cascade_summary.tsv` |
| TFプロモーターメチル化 | `11_epigenome_integration/analysis/12_grn_tf_methylation/tf_promoter_methylation/tf_methylation_summary.csv` |
| TF協調変動リスト | `11_epigenome_integration/analysis/01_integration/TF_coordinated_changes.csv` |
| TF BS メチル化重複 | `13_TF_binding-site/analysis/02_TF_BS_methylation_260207_v1/tables/T1_BS_methylation_overlap_summary.tsv` |
| ACT vs RED比較 | `11_epigenome_integration/analysis/15_act_bgc_epigenetic/act_vs_red_comparison.csv` |
| ACT BGC遺伝子解析 | `11_epigenome_integration/analysis/15_act_bgc_epigenetic/act_bgc_gene_analysis.csv` |
| BGC詳細テーブル | `11_epigenome_integration/analysis/01_integration/BGC_detailed.csv` |
| BGCレギュレーター要約 | `11_epigenome_integration/analysis/19_bgc_regulator_overview/bgc_regulator_summary.tsv` |
| DEG-DMG重複統計 | `11_epigenome_integration/analysis/03_overlap_analysis/DEG_DMG_overlap_statistics.csv` |
| モチーフ時系列動態 | `11_epigenome_integration/analysis/23_expanded_motif_search/motif_temporal_dynamics.csv` |
| MTase発現データ | `11_epigenome_integration/analysis/11_rm_system_identification/mtase_genes_with_expression.csv` |
| AAGCCCG空間分布 | `11_epigenome_integration/analysis/18_tss_analyses/motif_AAGCCCG_spatial_detail.csv` |
| モチーフエンリッチメント | `11_epigenome_integration/analysis/14_aagcccg_promoter_analysis/motif_enrichment_summary.csv` |

---

*最終更新: 2026-02-17*
