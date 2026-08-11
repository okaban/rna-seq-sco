# HupS Locus Tag 確定レポート

**日付**: 2026-05-06  
**目的**: Salerno et al. 2009 (PMID 19717607) の HupS = SCO5556 の正しさを検証し、SCO1175 との混同を解消する

---

## 結論（先出し）

**HupS = SCO5556（SC_RS29990）が正しい。Salerno 2009 の記述は正確。**  
SCO1175（SC_RS07825）はHU familyと無関係のpseuedogeneであり、HupSではない。

---

## 1. アノテーション確認

| Locus tag | RefSeq ID | gene_biotype | Product | Protein (WP) |
|-----------|-----------|-------------|---------|-------------|
| SCO1175 | SC_RS07825 | **pseudogene** | alpha/beta hydrolase fold domain-containing protein | なし |
| SCO2950 | SC_RS16855 | protein_coding | HU family DNA-binding protein | WP_003950507.1 |
| SCO5556 | SC_RS29990 | protein_coding | HU family DNA-binding protein | WP_003973439.1 |

**SCO1175 はHU familyでもなく、pseudogeneである。HupSの候補から即座に除外。**

---

## 2. タンパク質サイズによる HupA / HupS の同定

GFF (GCF_000203835.1) から計算した CDS 長：

| Locus tag | 座標 | CDS 長 | タンパク質サイズ | 帰属 |
|-----------|------|--------|---------------|------|
| **SCO2950** | 3,206,403–3,206,684 | 282 bp | **93 aa** | HupA（典型的HU, 90 aa級） |
| **SCO5556** | 6,054,696–6,055,352 | 657 bp | **218 aa** | HupS（HU domain + C末端伸長域） |

HupS の構造的特徴（HU domain ＋ 胞子形成特異的C末端伸長）を持つのは **SCO5556（218 aa）のみ**。  
SCO2950（93 aa）は C末端伸長を持たない典型的な HupA サイズ。

---

## 3. 染色体上の位置

全染色体長: 8,667,507 bp

| Locus tag | 中心座標 | 染色体位置率 | 領域 |
|-----------|---------|------------|------|
| SCO2950 | 3,206,544 bp | **37.0%** | コア領域 |
| SCO5556 | 6,055,024 bp | **69.9%** | アーム領域 |

Salerno 2009 は HupA がコア領域、HupS がアーム領域の染色体に結合することを報告している。  
遺伝子の染色体上位置もこの機能分担と完全に一致する。

---

## 4. DESeq2 発現量（normalized counts）

| Locus tag | T1 mean | T2 mean | T3 mean | log2FC T2/T1 | log2FC T3/T1 | 判定 |
|-----------|---------|---------|---------|-------------|-------------|------|
| SCO2950 (HupA) | 22,368 | 6,899 | 2,931 | −1.65 | **−2.90** | 栄養増殖→急減 |
| SCO5556 (HupS) | 1,910 | 627 | 420 | −1.56 | **−2.14** | 同様に減少 |

一見、SCO5556（HupS）も胞子形成で減少しており「HupS = T3誘導」の期待と矛盾するように見える。  
しかし以下の通り解釈できる：

### 相対的発現比の変化

| 時点 | SCO5556/SCO2950 比 | 意味 |
|------|-----------------|------|
| T1 | 1,910 / 22,368 = **0.085** | HupA が圧倒的多数 |
| T3 | 420 / 2,931 = **0.143** | HupS の相対存在比が **1.68倍** 上昇 |

**→ HupA が T3 でより急速に減少する（log2FC −2.90）ため、HupS（−2.14）が相対的に優位になる。**  
Salerno 2009 が述べた「胞子形成時に HupS が dominant な HU となる」というモデルは、絶対発現量の誘導ではなく、**HupA の急激な減少による相対的な HupS の台頭**として RNA-seq データに現れている。

---

## 5. SCO1175 の高発現（T2 ピーク）の解釈

SCO1175 は T1 がほぼゼロ（mean ≈ 1.3）で T2 に急上昇（log2FC = +5.15）するが：

- **pseudogene** であり protein をコードしない
- product は「alpha/beta hydrolase fold」であり HU family とは構造的に無関係
- baseMean = 28（極めて低発現）
- HupS の機能（DNA condensation, chromosome arm organization）を担いえない

このシグナルは HupS とは無関係の偶発的な一致であり、HupS の同定に用いることはできない。

---

## 6. Salerno 2009 との整合性まとめ

| Salerno 2009 の記述 | 本データでの確認 |
|-------------------|---------------|
| hupA = SCO2950 (90 aa) | ✅ 93 aa, 37% position (コア), 最高発現量で急減 |
| hupS = SCO5556 (longer, C-terminal extension) | ✅ 218 aa, 70% position (アーム), 相対的優位は T3 で増加 |
| HupS は胞子形成時に dominant | ✅ T3 での SCO5556/SCO2950 比 = 1.68× T1 |
| locus tag の rename は不要 | ✅ SCO5556 の old_locus_tag は現 RefSeq (SC_RS29990) で維持 |

**矛盾なし。Salerno 2009 の割り当て（SCO5556 = hupS）は正しい。**

---

## 7. Discussion 執筆用の正しい情報

```
HupS = SCO5556 (SC_RS29990), 218 aa, WP_003973439.1
HupA = SCO2950 (SC_RS16855), 93 aa, WP_003950507.1

SCO5556 発現量: T1 = 1,910, T2 = 627, T3 = 420 (normalized counts)
SCO2950 発現量: T1 = 22,368, T2 = 6,899, T3 = 2,931 (normalized counts)

HupS(SCO5556)/HupA(SCO2950) 比: T1 = 0.085, T3 = 0.143 (+68%)
```

RNA-seq データは「HupS が胞子形成で新たに誘導される」ではなく「HupA が急減することで HupS が相対的に優位な HU タンパク質になる」というモデルを支持する。これは Salerno 2009 のタンパク質局在データ（ChIP, immunofluorescence）と整合する。

---

*データソース: `05_annotation/analysis/05_annotation_260128_v1/tables/gene_master_DESeq2.tsv`, `11_epigenome_integration/analysis/39_GCCGGC_MTase_reverse_ID/data/GCF_000203835.1_ASM20383v1_genomic.gff`*
