# 20遺伝子 配列・構造類似性解析サマリー

- 日付: 2026-05-06
- プロジェクト: S. coelicolor A3(2) M145 RNA-seq + Epigenome 統合
- 解析ディレクトリ: `analysis/57_temporal_dynamics_exposed_TF/`
- 性質: 既存結果ファイルの集約レポート（新規解析なし）

---

## 1. 結論

- n = **20 / 20 遺伝子** がカテゴリーAリストから抽出済み（欠損なし）
- **10 / 20** が文献由来 6 リファレンスのいずれかと **Pfam ファミリーレベルで一致**
- **3 / 20** が BLASTP E < 1e-3 のヒットを保有（pident ≥ 30% を満たすのは **1 / 20** のみ）
- 結論: 「20 遺伝子は、他生物でメチル化結合制御因子として記載されている TF ファミリーを **Pfam レベルで** カバーする」。**オーソロジー主張は不可**（Streptomyces TF は配列レベルで分岐が大きく、強い BLAST 一致なし）

---

## 2. 入力 20 遺伝子の TF ファミリー構成（`tf_family_input` 集計）

| TF ファミリー | 件数 |
|---|---|
| HTH (other) | 5 |
| TetR | 4 |
| Other regulatory | 4 |
| Sensor kinase | 2 |
| Response regulator | 2 |
| MarR | 1 |
| LysR | 1 |
| Sigma factor | 1 |
| **合計** | **20** |

選択基準: TSS ±500 bp メチル化変化 AND \|log2FC\| ≥ 1 AND 同方向（`tables/exposed_TF_notable_list.tsv`）

---

## 3. 文献リファレンス 6 件との Pfam-overlap 集計

`notable_genes_pfam_overlap_matrix.tsv` および `notable_genes_similarity_summary.txt` より:

| 文献ファミリー（リファレンス） | Pfam 一致数 |
|---|---|
| TetR/AcrR (Casadesús & Low 2006) | 4 |
| TCS RR GerE/Receiver (S. coelicolor RamR analogue) | 3 |
| MarR (E. coli MarR; methylation-context regulator) | 1 |
| LysR/OxyR (E. coli OxyR; agn43 Dam-switch) | 1 |
| Sigma-70 ECF (M. tuberculosis SigB; m6A-near-σ-box) | 1 |
| Lrp/AsnC (E. coli Lrp; pap GATC switch) | 0 |

リファレンス 6 件の出典は `data/reference_proteins.tsv`（MarR/OxyR/Lrp E. coli K-12、HP1021 H. pylori 26695、SigB M. tuberculosis H37Rv、RamR S. coelicolor）

---

## 4. 個別ヒット（Pfam ファミリー一致したもの）

| locus_tag | old_locus_tag | tf_family_input | top_pfam | 文献ファミリー一致 |
|---|---|---|---|---|
| SC_RS22780 | SCO4122 | MarR | MarR (PF01047, E=5.2e-09) | MarR |
| SC_RS02455 | SCO0089 | LysR | LysR_substrate (PF03466, E=1.0e-35) | LysR/OxyR |
| SC_RS38785 | SCO7314 | Sigma factor | Sigma70_r2 (PF04542, E=2.8e-18) | Sigma-70 ECF |
| SC_RS17795 | SCO3134 | Response regulator | GerE (PF00196, E=1.3e-17) | TCS RR GerE |
| SC_RS35610 | SCO6685 (ramR) | Response regulator | GerE (PF00196, E=8.5e-14) | TCS RR GerE |
| SC_RS06310 | SCO0877 | HTH (other) | AAA_16 (PF13191, E=6.4e-18) | TCS RR GerE |
| SC_RS32025 | SCO5956 | TetR | TetR_N (PF00440, E=3.8e-09) | TetR/AcrR |
| SC_RS13145 | SCO2223 | TetR | TetR_C_46 (PF21943, E=4.1e-36) | TetR/AcrR |
| SC_RS25370 | SCO4639 | TetR | TetR_N (PF00440, E=2.9e-07) | TetR/AcrR |
| SC_RS33745 | SCO6299 | TetR | TetR_C_6 (PF13977, E=2.0e-25) | TetR/AcrR |

唯一 BLASTP で pident ≥ 30% を満たすのは SC_RS02455 (SCO0089, LysR) → OxyR_E._coli pident=32.0%, len=297, E=4.4e-29。次点が SC_RS22780 (SCO4122, MarR) → MarR_E._coli pident=22.5%, len=129, E=3.1e-06。

---

## 5. HTH_24 / Sigma70_r4_2 を除外した場合の分類

両ドメインは多数のファミリーで共有されるため、ファミリー判別性が弱い。`all_pfam_domains` 列での内訳:

- **HTH_24 を持つもの**: SC_RS22780 (SCO4122, MarR) **1 件のみ**
  - HTH_24 を除外しても MarR ファミリー判定は `MarR (PF01047)` で維持される
- **Sigma70_r4_2 を持つもの**: 4 件
  - SC_RS17795 (SCO3134) — `GerE, Response_reg` 残存 → TCS RR 判定維持
  - SC_RS35610 (SCO6685, ramR) — `GerE, HTH_23, HTH_28, HTH_IclR, Sigma70_r4` 残存 → TCS RR 判定維持
  - SC_RS38785 (SCO7314) — `Sigma70_r2, Sigma70_r3, Sigma70_r4` 残存 → Sigma-70 ECF 判定維持
  - SC_RS06310 (SCO0877) — `AAA_16, GerE, HTH_20, HTH_23` 残存 → TCS RR 判定維持

**結論: HTH_24 と Sigma70_r4_2 を除外しても、表 4 の文献ファミリー一致 10/20 の判定は不変。**判別はそれぞれの top_pfam（MarR, LysR_substrate, Sigma70_r2, GerE, TetR_N, TetR_C_46, TetR_C_6 など）に依拠している。

---

## 6. 主要な留保（既存サマリーから転記）

- 全 BLASTP top hit は弱い一致（pident < 35% が大半、E > 1e-3 が大半）。**オーソロジー主張は不可**
- Pfam ファミリー共有は **「メカニズムが同一」を意味しない**。Dam/Dcm 様メチル化スイッチが S. coelicolor で実際に作動するかは未検証（メチル部位変異・ChIP-seq・EMSA が必要）
- リファレンスセットの caveats:
  - UniProt Q9XAP2 (一部 DB で "RamR" と記載) は実体が **RamC (Radical_SAM + CofH/MqnC)**。SCO6685 RR の正対照には不適。GerE+Response_reg Pfam ペアを TCS-RR の判定基準として代替使用
  - UniProt O25617 (HP1021) は E ≤ 1e-3 で強い Pfam ヒットなし。HP1021 ファミリー判定は文献依拠

---

## 7. 結果ファイルの場所

### 一次出力（`analysis/57_temporal_dynamics_exposed_TF/`）

| ファイル | 内容 |
|---|---|
| `data/categoryA_genes.tsv` | 20 遺伝子の入力リスト + tf_family_input + log2FC + Δmethyl |
| `data/categoryA_proteins.faa` | 20 遺伝子のタンパク質配列（GenBank `ref.gbk` 由来） |
| `data/categoryA_pfam.tblout` | hmmscan vs Pfam-A.hmm 結果（E ≤ 1e-3） |
| `data/categoryA_pfam.domtblout` | hmmscan ドメイン詳細出力 |
| `data/categoryA_vs_refs.blast.tsv` | BLASTP vs 文献 6 リファレンス |
| `data/reference_proteins.faa` / `.tsv` | UniProt 由来文献リファレンス 6 件 |
| `results/notable_genes_similarity_analysis.tsv` | 統合表（20 行 × Pfam/BLAST/literature_family_match） |
| `results/notable_genes_pfam_overlap_matrix.tsv` | 20 × 6 Pfam-family overlap マトリクス |
| `results/notable_genes_similarity_summary.txt` | 既存サマリー（method, headline, caveat） |
| `results/notable_genes_for_literature_check.txt` | 文献検討用クエリリスト |

### 主要図

- `analysis/02_publication_figures/Figure3_20genes_dynamics.{pdf,png}`
- `analysis/02_publication_figures/Figure3_20genes_dynamics_updown.{pdf,png}`

### 関連スクリプト

- `analysis/52_shielded_exposed_boundary/scripts/H29c_pfam_regulatory_gene_selection.py` (Pfam-based 制御因子選定の上流スクリプト)

---

## 8. 解析メソッド（既存サマリーから転記）

- 入力: 20 遺伝子（カテゴリーA = `tables/exposed_TF_notable_list.tsv` の "category contains A"）
- タンパク質配列: `/Users/okaban/bioinfo/methyl/260102_M145/data/ref.gbk` (NCBI RefSeq GCF_000203835.1, NC_003888.3) から locus_tag/old_locus_tag で抽出
- Pfam: hmmscan vs Pfam-A.hmm（E ≤ 1e-3）
- 配列類似性: blastp vs UniProt 6 リファレンス（緩い E-value 10、tabular 出力）
- ファミリー判定ルール: クエリの Pfam accession 集合がリファレンスのファミリーグループ（`FAMILY_GROUPS` in script）と重なる場合に "family analogue" と判定

---

## 9. 次の検討（任意、今回は実施せず）

- HTH_24 / Sigma70_r4_2 を陽に除外したサブセットでの再計上は本レポート 5 節で実施済み（結論不変）
- 配列が著しく分岐するため、必要なら **structural similarity (Foldseek vs AlphaFold DB)** で再評価
- "メチル化スイッチが作動するか" の検証は別途 wet 実験 / 過去 ChIP-seq との照合が必要
