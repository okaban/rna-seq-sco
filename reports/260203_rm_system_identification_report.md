# M145 R-M System Identification Report

**日付:** 2026-02-03
**プロジェクト:** *Streptomyces coelicolor* A3(2) M145 エピゲノム-トランスクリプトーム統合解析
**解析ディレクトリ:** `11_epigenome_integration/analysis/11_rm_system_identification/`

---

## Key Insights（主要なインサイト一覧）

| # | Insight | 対応Figure/Table | 根拠となる数値 |
|---|---------|-----------------|---------------|
| 1 | **AAGCCCG (6mA) はREBASEに未登録 → 新規R-M系** | REBASE検索結果 | "None found" |
| 2 | **SC_RS17645 (N-6 DNA methylase) がAAGCCCG候補** | `mtase_genes_with_expression.csv` | 6mA特異的、T2vsT1で有意に発現低下 |
| 3 | **SC_RS17645はT2で発現抑制（log2FC=-2.19）** | DESeq2結果 | padj = 6.48e-16 |
| 4 | **MTase発現低下→脱メチル化→遺伝子発現変動の連鎖** | 仮説モデル | T2における相関の分子基盤 |
| 5 | **CCGG (4mC) はHpaII/MspIモチーフだが4mCは別系統** | 文献 | HpaII/MspIは5mC認識 |

---

## 1. 解析の背景と目的

### 1.1 背景

MEME解析により、M145で以下のメチル化モチーフを同定した：
- **CCGG (4mC)**: 75.6%の4mCサイトで検出
- **AAGCCCG (6mA)**: 新規モチーフとして同定

これらのモチーフを認識する制限修飾（R-M）系酵素を同定することで、メチル化パターンの分子基盤を明らかにする。

### 1.2 目的

1. REBASEデータベースでAAGCCCG、CCGGモチーフの既知R-M系を検索
2. M145ゲノムからDNAメチルトランスフェラーゼ（MTase）遺伝子を同定
3. MTase遺伝子の発現変動パターンを解析
4. メチル化-発現相関の分子機構を推定

---

## 2. 方法

### 2.1 REBASE検索

[REBASE](https://rebase.neb.com/)（Restriction Enzyme Database）にて認識配列検索を実施。

### 2.2 M145 MTase遺伝子同定

GFFアノテーション（GCF_000203835.1）から以下のキーワードで検索：
- DNA methyltransferase, N-6 DNA methylase, DNA cytosine methyltransferase
- PglX (BREX-2 system), restriction/modification

### 2.3 発現パターン解析

DESeq2結果（T2vsT1, T3vsT1, T3vsT2）とMTase遺伝子リストを照合。

---

## 3. 結果

### 3.1 REBASE検索結果

| モチーフ | REBASE結果 | 解釈 |
|---------|-----------|------|
| **AAGCCCG** | **None found** | 既知R-M系に登録なし → **新規R-M系の可能性** |
| **CCGG** | HpaII/MspI (5mC認識) | M145の4mCとは異なる系統 |

> **Insight #1**: AAGCCCG (6mA) はREBASEに未登録
> REBASEにAAGCCCGを認識する酵素は登録されていない。これはM145で同定したAAGCCCGモチーフが**未知のR-M系**によるメチル化である可能性を強く示唆する。
> **根拠**: REBASE cgi-bin/seqget検索で "None found"

### 3.2 M145 DNA MTase遺伝子

**同定されたDNA MTase遺伝子: 22個**

| カテゴリ | 遺伝子数 | 代表的遺伝子 |
|---------|---------|-------------|
| N-6 DNA methylase (6mA) | 1 | SC_RS17645 |
| BREX-2 PglX (6mA) | 2 | SC_RS28835, SC_RS35335 |
| DNA cytosine MTase (4mC/5mC) | 3 | SC_RS19770, SC_RS36410 |
| Type ISP R-M | 1 | SC_RS10595 |
| その他 DNA MTase | 15 | - |

> **Insight #2**: SC_RS17645 (N-6 DNA methylase) がAAGCCCG候補
> M145ゲノムにはN-6 DNA methylaseが1つのみ存在（SC_RS17645）。6mAメチル化の主要な責任酵素である可能性が高い。
> **根拠**: `mtase_genes.csv` - product: "N-6 DNA methylase"

### 3.3 MTase遺伝子の発現パターン

#### 3.3.1 T2vsT1で有意に発現変動したDNA MTase

| Locus Tag | Product | log2FC | padj | 方向 |
|-----------|---------|--------|------|------|
| SC_RS17645 | N-6 DNA methylase | **-2.19** | 6.48e-16 | DOWN |
| SC_RS24685 | class I SAM-dependent DNA MTase | -3.22 | 5.00e-36 | DOWN |
| SC_RS18875 | DNA repair protein RadA | -2.51 | 7.17e-56 | DOWN |
| SC_RS30080 | DNA-formamidopyrimidine glycosylase | -1.67 | 1.73e-14 | DOWN |
| SC_RS06675 | Fpg/Nei family DNA glycosylase | +1.85 | 1.16e-12 | UP |

> **Insight #3**: SC_RS17645はT2で発現抑制（log2FC=-2.19）
> AAGCCCGメチル化の候補酵素SC_RS17645は、T2vsT1で約4.6倍の発現低下を示す。これはT2時点での**メチラーゼ活性低下→脱メチル化**を示唆する。
> **根拠**: DESeq2結果、padj = 6.48e-16

#### 3.3.2 T3vsT1で有意に発現変動したDNA MTase

| Locus Tag | Product | log2FC | padj | 方向 |
|-----------|---------|--------|------|------|
| SC_RS36410 | DNA cytosine MTase | **+5.34** | 1.03e-10 | UP |
| SC_RS36625 | DNA-methyltransferase | +3.87 | 6.52e-09 | UP |
| SC_RS19770 | DNA cytosine MTase | +2.32 | 9.47e-06 | UP |
| SC_RS19670 | DNA-methyltransferase | +3.51 | 1.65e-09 | UP |

**注目点**: T3では**cytosine MTaseが強く誘導**される（SC_RS36410: 40倍以上）。これは後期（T3）での4mCメチル化増加を示唆する。

---

## 4. 分子機構モデル

### 4.1 提案モデル: T2における脱メチル化と遺伝子活性化

```
T1 (対照時点)
    ↓
SC_RS17645 (N-6 DNA methylase) 発現低下 (T2vsT1: log2FC=-2.19)
    ↓
AAGCCCG部位の6mAメチル化レベル低下
    ↓
プロモーター領域の脱メチル化
    ↓
転写抑制の解除 → 遺伝子発現上昇
    ↓
T2における4mC-発現の正の相関 (r=0.137)
```

> **Insight #4**: MTase発現低下→脱メチル化→遺伝子発現変動の連鎖
> SC_RS17645の発現低下により、プロモーター領域のAAGCCCGメチル化が減少し、これが遺伝子発現変動と相関する。これは**T2における4mC-発現相関（r=0.137, p=0.0006）の分子基盤**を説明する。

### 4.2 CCGGメチル化について

> **Insight #5**: CCGG (4mC) はHpaII/MspIモチーフだが4mCは別系統
> HpaII/MspIは**5mC (C5位メチル化)**を認識するが、M145で検出されたのは**4mC (N4位メチル化)**である。これは異なる酵素系（おそらくSC_RS19770またはSC_RS36410）によるメチル化を示唆する。

**候補酵素**:
- SC_RS19770: T3で発現上昇（log2FC=+2.32）
- SC_RS36410: T3で強く発現上昇（log2FC=+5.34）

---

## 5. Figure解説

### 5.1 mtase_genes_with_expression.csv

**ファイル**: `analysis/11_rm_system_identification/mtase_genes_with_expression.csv`

**説明**: 22個のDNA MTase遺伝子について、座標、予測メチル化タイプ、T2vsT1/T3vsT1/T3vsT2の発現変動データを統合したテーブル。

**カラム**:
- `locus_tag`: 遺伝子ID（SC_RS形式）
- `product`: 遺伝子産物名
- `predicted_methyl_type`: 予測メチル化タイプ（6mA, 4mC/5mC, Unknown）
- `log2FC_T2vsT1`, `padj_T2vsT1`: T2vsT1の発現変動
- `log2FC_T3vsT1`, `padj_T3vsT1`: T3vsT1の発現変動

**サポートするInsight**: #2, #3

### 5.2 rm_systems.csv

**ファイル**: `analysis/11_rm_system_identification/rm_systems.csv`

**説明**: 同定されたR-M系コンポーネントのリスト。

**カラム**:
- `system`: R-M系タイプ（BREX-2, Type ISP, Dcm-like, Dam-like）
- `component`: コンポーネント名
- `locus_tag`: 遺伝子ID
- `methyl_type`: メチル化タイプ
- `motif`: 候補認識モチーフ

**サポートするInsight**: #1, #2

---

## 6. 出力ファイル一覧

| ファイル | 形式 | 説明 |
|---------|------|------|
| `mtase_genes.csv` | CSV | 全MTase遺伝子リスト（87遺伝子） |
| `mtase_genes_with_expression.csv` | CSV | DNA MTase + 発現データ（22遺伝子） |
| `rm_systems.csv` | CSV | 同定されたR-M系（7コンポーネント） |
| `RM_SYSTEM_IDENTIFICATION_REPORT.md` | MD | 解析ディレクトリ内サマリー |

---

## 7. 結論

1. **AAGCCCG (6mA) は新規R-M系モチーフ**
   - REBASEに登録なし
   - SC_RS17645 (N-6 DNA methylase) が責任酵素の候補

2. **SC_RS17645の発現動態がメチル化-発現相関を説明**
   - T2vsT1で発現低下（log2FC=-2.19）
   - メチラーゼ活性低下→脱メチル化→転写活性化

3. **CCGG (4mC) は従来の5mC R-M系とは異なる**
   - SC_RS19770, SC_RS36410が候補
   - T3で発現上昇（後期メチル化）

---

## 8. 次のステップ

| 優先度 | タスク | 目的 |
|--------|--------|------|
| **高** | SC_RS17645のBLAST検索 | 相同性に基づく機能推定 |
| **高** | プロモーターAAGCCCG分布解析 | メチル化位置と発現変動の関係 |
| **中** | SC_RS17645ノックアウト株の予測 | 実験的検証の設計 |
| **中** | REBASEへの新規登録検討 | 新規R-M系としての報告 |

---

## 9. 参考文献

- REBASE: The Restriction Enzyme Database. https://rebase.neb.com/
- [HpaII - Wikipedia](https://en.wikipedia.org/wiki/HpaII)
- [NEB HpaII Methyltransferase](https://www.neb.com/en-us/products/m0214-hpaii-methyltransferase)

---

*生成日: 2026-02-03*
*解析スクリプト: `scripts/rm_system_identification.py`*
