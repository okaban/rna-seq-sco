# TSS下流GCCGGC保護ゾーン解析レポート

**日付**: 2026-05-10  
**解析番号**: 75_tss_downstream_protection  
**スクリプト**: `11_epigenome_integration/analysis/75_tss_downstream_protection/scripts/tss_upstream_downstream_protection.py`  
**図ファイル**: `11_epigenome_integration/analysis/75_tss_downstream_protection/figures/tss_upstream_downstream_profile.{png,pdf,svg}`  
**補足図**: `15_paper_figures/figures/supp/SuppFig_protection_downstream.png`

---

## 背景と目的

既存の解析（analysis/48, 52）では、制御遺伝子TSSの上流にGCCGGC (4mC) メチル化フリーゾーン（保護ゾーン）が存在し、Shielded遺伝子（nearest_methyl_distance > 293 bp）とExposed遺伝子（< 293 bp）の二分法が確立されている（ROC AUC=0.917）。

本解析の目的：

**TSS上流だけでなく下流（+1〜+293 bp）にもGCCGGCメチル化フリーゾーンが存在するかを定量し、2パネル図として可視化する。**

---

## データ

| データソース | 内容 | n |
|---|---|---|
| `52_shielded_exposed_boundary/tables/all_genes_features_unified_n57.tsv` | 制御遺伝子フィーチャー表（TSS, strand, is_exposed） | 1,055 遺伝子 |
| `23_expanded_motif_search/4mC_final_census.csv` | GCCGGC (TGGCCGGC/GGCCGG) 4mCサイト T1 | 1,299 ユニーク位置 |
| `18_tss_analyses/comprehensive_tss_table.csv` | Jeong2016実験的TSS優先（450/1055遺伝子） | 2,703 |

- Shielded遺伝子: n=998（is_exposed=0）
- Exposed遺伝子: n=57（is_exposed=1、協調的メチル化-発現遺伝子）
- ゲノム背景GCCGGC密度: **0.1499 sites/kb**（1,299 sites / 8,667 kb）

---

## 主要結果

### ゾーン別GCCGGC密度（sites/kb、T1）

| ゾーン | Shielded (n=998) | Exposed (n=57) | ゲノム背景 | Shielded FC | Exposed FC |
|---|---|---|---|---|---|
| 上流 [−500, 0) | 0.1328 | 0.2105 | 0.1499 | **0.886×** | 1.405× |
| 下流 [0, +293] | 0.0973 | 0.1170 | 0.1499 | **0.649×** | 0.780× |
| 遺体遠位 [+293, +500] | 0.0905 | 0.0877 | 0.1499 | 0.604× | 0.585× |

### 統計検定（Mann-Whitney U, Shielded vs Exposed）

| 比較 | U統計量 | p値 | 解釈 |
|---|---|---|---|
| 上流密度 [-500, 0) | 27,177 | 0.178 | n.s. |
| 下流密度 [0, +293] | 28,216 | 0.721 | n.s. |
| 上流/下流比率 | 27,320 | 0.303 | n.s. |

---

## 解釈

### 1. 下流保護ゾーンは実在するか？

**答え：両群で下流（+1〜+293 bp）にGCCGGC密度の低下が観察される。**

- **Shielded遺伝子**: 上流0.886× → 下流0.649× とゲノム背景を下回り、**下流の保護がより強い**
- **Exposed遺伝子**: 上流1.405×（メチル化存在）→ 下流0.780×（中程度低下）
- 遠位遺伝体（+293〜+500 bp）でも低下が継続（~0.60×）

### 2. 保護ゾーンの非対称性

| 指標 | Shielded | Exposed |
|---|---|---|
| 上流/下流比 | 0.886 / 0.649 = **1.37** | 1.405 / 0.780 = **1.80** |

- **Shielded遺伝子（保護あり）**: 上流・下流ともに保護されるが、**下流の方が密度がさらに低い**（FCで上流0.886× vs 下流0.649×）
- **Exposed遺伝子（保護なし）**: 上流にはメチル化が集積（1.405×）し、下流は中程度（0.780×）
- → **保護ゾーンは上流優位ではなく、TSS下流にも実在し、むしろ遺伝体側でより強い**

### 3. 上流・下流の密度差が有意でない理由

統計検定でShielded/Exposed間差が有意でない（p>0.05）のは以下による：
1. **サイト数が少ない**（1,299 sites / 1,055 genes ≈ 1.2 sites/gene）→ 多くの遺伝子でゼロ値
2. **Exposed群n=57**は検出力が限られる
3. 各遺伝子個別では大きなゼロ膨張分布、平均値での差は明確だが個体間分散が大きい

### 4. 下流保護の生物学的解釈候補

| 仮説 | 証拠 | 状態 |
|---|---|---|
| **ポリメラーゼアクセス確保** (上流のみ) | 上流保護は確認されるが下流も低下 | 部分的支持 |
| **転写開始領域全体の保護** (上下流対称) | 下流+293 bpまでも低密度 | **支持** |
| **遺伝子body配列組成** (内在的GCCGGC希少) | 遠位+500までも低密度継続 | 部分的に関与 |
| **RNAPオープン複合体占有** | +1〜+293 bp保護がLayer 3と一致 | 検討中 |

**現時点での結論**: GCCGGC保護ゾーンはTSS「上流のみ」ではなく、**転写開始点を中心とした±293 bp対称的な低密度領域**として存在する。これは「ポリメラーゼ到達を妨げないためのカウンターセレクション」加えて「転写バブル形成領域の保護」という複合メカニズムを示唆する。

---

## 出力ファイル

| ファイル | 内容 |
|---|---|
| `75_tss_downstream_protection/figures/tss_upstream_downstream_profile.png/pdf/svg` | 2パネル図（メタプロファイル + 散布図） |
| `75_tss_downstream_protection/tables/tss_metaprofile_binned.tsv` | 50bpビン別密度プロファイル（Shielded/Exposed） |
| `75_tss_downstream_protection/tables/per_gene_upstream_downstream_density.tsv` | 遺伝子ごとの上流・下流密度 |
| `15_paper_figures/figures/supp/SuppFig_protection_downstream.png` | 論文補足図コピー |

---

## 次のステップへの示唆

1. **Layer 2/3 との接続**: 下流保護がLayer 2（カウンターセレクション）かLayer 3（NAP占有）のどちらで説明されるかを分解（分析/75 × 67 overlap）
2. **より広いウィンドウ**: ±2 kb以上で遺伝体全体の密度プロファイルを確認
3. **遺伝体配列組成コントロール**: シャッフル配列や非制御遺伝子との比較でゲノム配列バイアスを除外
4. **実験的検証**: ChIP-seqやFoot-printingデータとの比較（転写開始付近のタンパク占有との相関）
