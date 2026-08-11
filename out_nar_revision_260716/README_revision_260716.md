# NAR revision — 16-comment triage + HITL review (2026-07-16) — FINAL
# S. coelicolor M145 epigenome × transcriptome / 中心テーゼ = permissive spatial organiser

全16コメント対応完了。**本文edit（EN+JP）は著者承認のうえ適用済み**（in-place, git-reversible;
コミットは著者側）。図はスクリプトを編集して非破壊で再生成し、原稿スロットへ反映済み。10の
著者決定点のうち9が解決済み。残り1（決定点2＝Intro削除済みの「~2,200 bpゾーン」主張を完全に落とすか、
制御遺伝子パネルへ配置し直して帰属するか）は著者のサブ判断待ち（09参照）。全成果物は artifact として保存済み。

## 読む順序（推奨）
1. **11_cover_letter_draft.md** — エディタ宛カバーレター草案（6章＋照会2点）
2. **09_response_skeleton.md** — 16コメント point-by-point 対応（FINAL・全決定点RESOLVED）
3. **10_applied_body_edits.md** — EN+JPに適用した本文edit全6件の正確な差分＋判断根拠
4. **02_strategy_memo.md** / **08_strategist_review.md** — 戦略診断と専門家レビュー
5. 個別: 07(中心メッセージ) / 06(図2) / 05(図修正ログ) / 04(C2/C8/C15/C16) / 03(サーベイ)

## 成果物一覧
### 戦略・対応文書
| ファイル | 内容 | version_id |
|---|---|---|
| 00_comment_ledger.md | 16コメント正規化台帳 | d8e5f840-eeab-4c13-9ef8-172fe07e91ae |
| 01_grounding_contract.md | 全エージェント共通の接地契約 | 24d82d33-d753-4aad-a4cf-6ed4c47c0b7c |
| 02_strategy_memo.md | 中心メッセージ再定義＋トリアージ | c7b59117-7eb1-4d20-a54e-70549a2eeb7f |
| 03_prokaryote_epigenetics_survey.md | 原核エピジェネ先行研究サーベイ（41文献検証） | 971d6935-5f40-4b50-9323-3a8329f65746 |
| 04_analysis_text_proposals.md | C2/C8/C15/C16 検証＋diff | 8ca33b09-c6ac-4d9a-b65b-2dac19713c25 |
| 05_figure_fixes_log.md | 図修正の証拠ログ（RESOLVED判定反映済） | 4ef599a1-f52d-4671-bdf6-8ee64ba3d713 |
| 06_figure2_restructure_proposal.md | C9/C10/C12/C13/C14 図2再構成 | 785e7347-fa45-46cf-8266-c6a043f383f9 |
| 07_central_message_reframe.md | C1中心メッセージ再フレーミング（ROUND2） | d3ac420f-197e-4690-b920-7eb609b9f03d |
| 08_strategist_review.md | 戦略専門家レビュー | 3c58e534-b08b-4b3a-baf2-4e04671fe31f |
| 09_response_skeleton.md | 査読対応（9/10決定点RESOLVED・点2は著者サブ判断） | b449afd6-ac3c-46be-b2d9-a6150bc3f117 |
| 10_applied_body_edits.md | 適用済み本文edit6件（EN+JP）の差分＋判断 | 22f0780a-a136-402a-8901-66e54752f5bb |
| 11_cover_letter_draft.md | カバーレター草案 | 97cc149a-216e-4879-9080-f190ff5b2f25 |
| 12_commit_prep.sh | 2リポジトリ分のcommit準備（著者が実行） | cfabf8ef-7204-4e56-b780-8694a357d491 |

### 最終図（原稿スロットへ反映済み）
| 図 | コメント | version_id |
|---|---|---|
| Figure 1 | C4/C5/C6＋パネル番号統一 | 16eeaa4d-08ad-4bde-aa43-a289921ef0f8 |
| Figure 2（3パネルA/B/C） | C14/C10/凡例余白化＋参照線クリップ＋パネル番号統一 | 37d286f6-9c4f-41ff-a749-a24faa07351d |
| Figure 3 | パネル番号統一 | 614c3975-1d1d-46df-8b65-39c4b5c553bb |
| Figure 4 | タイトル衝突解消＋パネル番号統一 | 8b9f8c43-72a7-460c-b078-6d556a904b28 |
| Figure 5 | タイトル衝突解消＋パネル番号統一 | 3a5cc494-96b8-411a-b19d-021451442c32 |
| Figure 6 | パネル番号統一 | 683fda18-5a7b-4ba9-8267-728393d754b7 |
| Figure 8 | パネル番号統一 | f3e480a5-4547-43b3-a81c-ed2e577b94d0 |
| SuppFigure 10 | C7（6mA=441） | f269637c-7cdf-461f-9a9a-dc7e3013919b |
| SuppFigure 15（成長統合） | C3 | 4ec67870-84fd-4ca2-8496-4e2080c69c5c |
| SuppFigure_neg_control（旧Fig2D） | C10 | b924d439-fb5f-4f9a-958d-f23d02051c7c |
| SuppFig_6mA_companion | C9 | fb16c7a0-9352-412b-a132-a76c1d371847 |
| SuppFig10_caption.md（移設したレジェンド） | C7 | f4541873-d559-465c-99e9-68e386cd700a |

### 成長データ（元CSV不在→SVGから実測値復元、原稿統計を再現）
| ファイル | version_id |
|---|---|
| growth_reconstructed_from_svg.csv | be7f59fb-93c4-40b0-8e64-63e5277ed45c |
| make_growth_merged.py | 5c16a7f2-67b9-41b3-b012-f569449e9a49 |
| PROVENANCE.md | acd58fbc-98df-4988-b608-405364dfed9d |

## 適用済みの変更
### 本文（EN: Manuscript_EN_NAR_numbered_inline.md / JP: Manuscript_JP_NAR_numbered_inline.md）
各6件、in-place・git-reversible: C8見出し改名 / C14 ±300 bp訂正 / C1 Intro階層化 /
C1 Abstract再フレーム（関係動詞）/ Abstract共動態(ρ=0.43)降格 / Red fold floor-aware文言。
※注意: EN原稿の `_inline` 版はgit未追跡（tracked版は `_numbered.md`）。commit時に要判断（12参照）。

### 図スクリプト（granted repo, git-reversible）
- `15_paper_figures/scripts/01_figure1_landscape.py` — C4/C5/C6
- `15_paper_figures/scripts/02b_figure2_RM_redistribution.py` — C14帯±300bp＋Fig2Dを標準補足へ分離（主図3パネル化）
- `11_epigenome_integration/analysis/77_reviewer_figures/fig_modpos_and_redistribution.py` — C7（4mC統一・6mA 441・グレー整理）

## 専門家プロファイル（作成済み）
- **PAPER_STRATEGY_STRATEGIST**（論文戦略ストラテジスト）: フル権限。ハンドリングエディタ視点で
  中心メッセージ・図ナラティブをトリアージ。canon固定・過剰主張抑制を内蔵。08レビュー実施。

### パネル番号の統一（全主図）
共有ヘルパー `00_shared_utils.add_panel_label` を「外側マージン固定・オフセット定数配置」に
書き換え、Fig1/2/3/6/8（および F4_F5 のローカル版）を全てこの1ヘルパー経由に統一。
太字・小文字 a/b/c（凡例の (a)(b)(c) と一致）・全パネル同サイズ・y軸ラベル幅に依存しない同一視覚位置。
旧状態: 図ごとに大文字/小文字・サイズ・per-panel x/y がバラバラだった。

## 作成したスキル（公開済み）
- **figure-comment-triage** — 査読コメントの正規化→接地→disposition→データ検証→階層化ワークフロー
- **svg-figure-data-recovery** — SVGから数値データを復元し既知統計量で検証するワークフロー
- **figure-legibility-qc** — 凡例×データ重なり・凡例を貫通する参照線・目盛り衝突・パネル番号不統一の
  修正/予防。「レンダリングを必ず目視する」工程を必須化。ヘルパー: `qc_add_panel_label` /
  `qc_headroom` / `qc_clip_span` / `qc_kilo_ticks`。

## 残タスク（著者/手動）
- git commit/push（12_commit_prep.sh; 認証は著者側。2リポジトリ: bioinfo/rna-seq と obsidian vault）
- EN `_inline` 原稿のgit追跡方針の判断（_numbered.md の再生成 or _inline を作業正本とする）
- Fig 2C/2B のキャプション補足（C12/C13; Figure-Legends へ; non-blocking）
- 照会2点: Mutilka同定（C16）/ 4mC直交検証の位置づけ
- 既存実験タスク: 4mC直交検証、k-mer対照2nd replicate、データ寄託(PRJNA891940)

## Google Drive アップロード（2026-07-16, cleanフォルダ 1yIayPiOqHZzSDRNrIntZxXTNeT0y8ud4）
著者チェック用に260716版をアップ済み（v8版も同フォルダに併存。混同回避のためファイル名にversion付与）。
| ファイル | Drive ID |
|---|---|
| Manuscript_EN_NAR_v260716.pdf（39p, 図修正反映） | 1saUT5T5A2R5JLAHfvNr25Yaa5t7rG4RZ |
| Manuscript_JP_NAR_v260716.pdf（38p） | 1_o68tkgBE15_l1hoWVP3RF6oC5x1nXuS |
| 解説レポート_epi-trans_260716_v9.pdf（22p） | 1_QuSMd65dCdTIOUVdZMlYll2L1_nXiqZ |
| 解説レポート_epi-trans_260716_v9.md | 1Zn3HpVrkKcH6AkHYsyqk1-TqG7tS8TAc |

日本語解説v9の主更新: v8→中心結論+2補助観察に階層化 / ρ=0.43をAbstractから降格 / 図2を3パネル化+±300bp訂正（決定点2は著者サブ判断未決を明記）/ 成長3パネル統合図に差替 / 全図パネル番号統一を反映。
アーティファクト: EN=08bdea57-ed0c-4cdb-ae8b-28ea408fe972 / JP=ea0a55e2-8008-4565-afcf-debeea62125d / 解説pdf=0b58640b-71a8-441d-8dd8-04b2e35b13f0 / 解説md=80daa220-0b2f-46a5-a609-2665d7eeeafb。

## 重要な原則（全作業で遵守）
- canon値は固定、数値は成果物/原稿行から取得（記憶からハードコードしない）
- 撤回内容(0.908 AUC・「57 Exposed」・two-blocs・eigengene・four-layer Gatekeeper)は再導入せず
- null結果は「反証テストが通った」証拠（本文維持、図パネルのみ補足送り）
- 本文editは著者承認のうえ適用（差分と根拠は10に全記録）
