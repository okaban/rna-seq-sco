# レビューカード — 原稿整合監査（数値・用語）
_stage `s2_consistency` · 2026-07-21T00:58:51_

## 主要数値
- **files_scanned**: 12
- **violations_critical**: 0
- **violations_informational**: 17
- **critical_files**: []
- **informational_files**: ['Key-Claims.md', 'Reviewer-Response-Draft.md', 'Story-Draft.md']
- **canonical_presence**: {'Manuscript_EN_consolidated_260615.md': {'has_62': True, 'has_989': True, 'has_1051': True, 'has_293bp': True}, 'Manuscript_EN_NAR_numbered.md': {'has_62': True, 'has_989': True, 'has_1051': True, 'has_293bp': True}}
- **violations_file**: out/consistency_violations.csv

## 要確認フラグ
- ✅ 投稿に入る章（build sources + consolidated）に撤回値/語なし。
- ℹ️ 参考: ビルド非対象の作業ドラフトに 17 件（['Key-Claims.md', 'Reviewer-Response-Draft.md', 'Story-Draft.md']）— アーカイブ候補。

## 次のアクション
違反行を確認・修正提案は s5 で。図/解析起因なら[PC]へ。

---
## 追加検証（agent, 2026-07-16）: 編集済みファイルの直接スキャン
自動監査は build source として `Manuscript_EN_consolidated_260615.md` と
`Manuscript_EN_NAR_numbered.md` を対象とするが、今回編集したのは `_inline` 版と2つの凡例
ファイルなので、これらを直接スキャンした。撤回語/値の全ヒットを文脈確認した結果、**すべて
「撤去済み」「主張しない」「WITHDRAWN」「retracted … を置換」という否定・撤回参照**であり、
生きた主張として復活したものは皆無。正準値(62/989/1,051/293)は全ファイルに存在。
- Manuscript_EN_NAR_numbered_inline.md: AUC×8 = 全て否定文（"no AUC is claimed/reported/removed"）
- Manuscript_JP_NAR_numbered_inline.md: 同上
- Figure-Legends.md: "Gatekeeper switch"/"four-layer" ×1 = 「撤回した…図を置換」の説明文; AUC×5 全否定
- Supplementary-Materials-List.md: "57 Exposed" ×1 = 取り消し線付きWITHDRAWN凡例; AUC×9 全否定
結論: **投稿パッケージに撤回内容のdrift無し（CRITICAL 0件）。**

## 凡例ファイルの更新（C7表記統一 + C10/C12/C13 反映）
- Figure-Legends.md: Fig2凡例を3パネル(A/B/C)へ更新; パネル(b)/(c)の相補性ノート(C13)追加;
  293 bpが二値分類ではなく連続分布上の操作閾値である旨(C12)明記; 表記 m4C→4mC 統一(全9箇所);
  Fig2(b)の n を 2,646 に訂正。
- Supplementary-Materials-List.md: 表記 m4C→4mC 統一(全9箇所);
  新規 **Supplementary Figure 19**（旧Fig2D＝陰性対照）の凡例を追加。
