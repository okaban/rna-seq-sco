# 15_paper_figures

## 概要

論文用の最終図表・原稿を管理するディレクトリ。全解析ステップの結果を統合し、出版用のメインフィギュア（Figure 1-5）、補足図表、および原稿テキストを生成する。

## ディレクトリ構造

```
15_paper_figures/
├── figures/
│   ├── main/                              # メインフィギュア（PDF/SVG）
│   │   ├── new_Figure1_methylation_landscape.*
│   │   ├── new_Figure2_gatekeeper_overview.*
│   │   ├── new_Figure3_exposed_TFs.*
│   │   ├── new_Figure4_negative_results.*
│   │   └── new_Figure5_switch_model.*
│   └── supplementary/                     # 補足図（FigS1-S8）
├── tables/
│   └── supplementary/                     # 補足テーブル
├── scripts/
│   ├── 00_shared_utils.py                 # 共通ユーティリティ
│   ├── 10-16_new_figure*.py               # 図表生成スクリプト群
│   └── 15-16_new_supplementary_*.py       # 補足図表スクリプト
├── manuscript/
│   ├── full_manuscript.md                 # 統合原稿
│   ├── 00_outline.md                      # 構成案
│   ├── 01_results.md ~ 06_references.md   # セクション別原稿
│   └── 07_style_guide.md                  # スタイルガイド
```

## 主要な出力

- `figures/main/` - 出版用メインフィギュア5点（PDF + SVG）
- `figures/supplementary/` - 補足図（FigS1-S16、PDF + SVG）
- `tables/supplementary/` - 補足テーブル
- `manuscript/full_manuscript.md` - 論文原稿（統合版）

## 2026-04-11 更新内容

今セッションで追加された補足図（FigS14–16）および原稿更新：

### 新規補足図
| 図 | スクリプト | 内容 |
|----|---------|------|
| FigS14 | `24_figS14_motif_reliability.py` | モチーフ内メチル化位置（AAGCCCG pos1=63%、GCCGGC回文）、モチーフ充足率（AAGCCCG=31.6%/GCCGGC=7.96%）、6mA Unassigned率（T1=72.2%） |
| FigS15 | `25_figS15_window_methylation.py` | ゲノムワイド50kbウィンドウ密度（4mCはT2ピーク→T3急減、6mAは単調増加、コア/アーム比） |
| FigS16 | `26_figS16_timepoint_TSS.py` | タイムポイント別TSS保護ゾーン（4mC制御遺伝子T3=1.046崩壊、6mAはT2最深保護=0.732） |

### 原稿更新箇所（01_results.md / full_manuscript.md）
- Section 1: モチーフ位置・占有率・Unassigned段落（A-1/A-2/A-3）を追加
- Section 1: 50kbウィンドウ密度動態段落（B-1）を追加
- Section 2（保護ゾーン）: タイムポイント別保護ゾーン動態段落（C-2）を追加
- Section 4（62 exposed TF）: GO/KEGG機能エンリッチメント段落（E-1）を追加
