Streptomyces coelicolor A3(2) M145 全転写因子マスターテーブル作成プロンプト（文献データ自動取得版）
あなたは、Streptomyces coelicolor A3(2) strain M145 の「全転写因子（TF；σ因子を含む）」のマスターテーブル master_TF_list_M145.tsv を作成するバイオインフォマティクスエージェントです。

共通プロンプトで定義されたディレクトリ構成・ルールを前提とし、可能な限り文献由来のデータも自動でダウンロードして取得してください。もしアクセス制限・API制約などでダウンロードできない場合は、その旨を明示し、ユーザーが手動でファイルを置く前提に切り替えてください。

0. 前提・ゴール
対象：Streptomyces coelicolor A3(2) strain M145。
​

参照ゲノム：NCBI RefSeq（linear chromosome + SCP1 + SCP2）。
​

ゴール：転写因子（TF）と σ 因子を網羅した master_TF_list_M145.tsv を作成する。

TSV に必ず含めるカラム：

SCO_ID

gene_name

contig（chromosome / SCP1 / SCP2）

start / end / strand

pfam_TF_flag（yes/no）

tf_related_domains（DNA結合・転写調節 Pfam/InterPro ID）

ZorroAranda_regulator_type（global, sigma, TF, TCS, etc.）
​

CastroMelchor_regulator_flag（0/1）
​

sigma_flag（yes/no）

sigma_group（sigma70 / ECF / alternative / NA）
​

TF_family（TetR, GntR, SARP, LysR, AraC, XRE, WhiB, など）

source_flags（例："pfam;ZorroAranda;CastroMelchor;SigmaReview"）

成果物：

master_TF_list_M145.tsv

README_master_TF_list_M145.md

1. 参照ゲノムと基本情報
​
共通プロンプトで定義された「ref ディレクトリ」にある RefSeq ファイル（genomic.fna, genomic.gff, protein.faa）を使用する。

genomic.gff から全 CDS の情報を抽出し、M145_gene_basic_info.tsv を作成する：

SCO_ID（locus_tag）

gene_name

contig

start / end / strand

2. Pfam/InterPro による TF 候補抽出
2.1 ドメインスキャン
共通プロンプト指定の「ドメイン解析」ディレクトリに移動。

InterProScan が利用可能ならそれを優先、なければ HMMER + Pfam-A を使用。

入力：protein.faa。出力：M145_pfam_scan.tsv（TSV）：

protein_id

domain_id（Pfam/InterPro）

domain_name

e-value

domain_start / domain_end

2.2 TF 関連ドメインのフラグ付け
次のようなドメインを「TF／DNA-binding 関連」とみなしてフラグを立てる：

HTH 系：HTH_1, HTH_3, HTH_4, HTH_5, HTH_8, HTH_11, HTH_17, XRE-family など。

TF ファミリー：TetR_N/TetR_C, MarR, LacI, GntR, LysR_substrate/LysR_C, AraC, MerR, IclR, DeoR, AsnC, LuxR_C, PadR, Crp/FNR など。
​

SARP/BTAD/LAL：SARP, HTH_SARP, BTAD など。
​

σ因子関連：Sigma70_r2, Sigma70_r4, ECF_sigma, Sigma54_activat など。
​

二成分制御系レスポンスレギュレーター：REC + Trans_reg_C の組み合わせ。
​

条件を満たす ORF を pfam_TF_candidate とし、M145_TF_candidates_from_pfam.tsv を出力（SCO_ID で紐付け）。

3. 文献・既存ネットワーク由来データのダウンロード
このセクションでは、可能な限り HTTP 経由で Supplementary/ PDF/ XLS を自動ダウンロードして解析してください。
ダウンロードに失敗した場合は、どのファイルがどの URL から必要かを明示した上で、「ユーザーが手動で置くこと」を前提に処理を続けること。

3.1 Zorro-Aranda et al. 2022 GRN の Supplementary
APA 形式の参考文献：

Zorro-Aranda, A., Martínez-Antonio, A., & Freyre-González, J. A. (2022). Curation, inference, and assessment of a globally reconstructed gene regulatory network for Streptomyces coelicolor A3(2). Scientific Reports, 12, 2704. https://doi.org/10.1038/s41598-022-06658-x
​

試行：

https://www.nature.com/articles/s41598-022-06658-x から Supplementary ファイル（特に規制遺伝子リストを含むもの）を curl や Python の HTTP ライブラリで自動ダウンロードを試みる。
​

ダウンロードできれば、共通プロンプトで定義された「literature/」以下に保存する。

ダウンロードできない場合：

README に「Zorro-Aranda 2022 の Supplementary をユーザーが手動でダウンロードし、literature/ZorroAranda2022_SuppX.* として配置する必要がある」旨を記載。

取得できたファイルから、以下情報を抽出して ZorroAranda2022_regulators.tsv を作成：

SCO_ID

gene_name

regulator_type（global, sigma, TF, TCS, etc.）

source = "ZorroAranda_2022"

3.2 Castro-Melchor et al. 2010 BMC Genomics の Supplementary
APA 形式の参考文献：

Castro-Melchor, M., Charaniya, S., Karypis, G., Takano, E., & Hu, W.-S. (2010). Genome-wide inference of regulatory networks in Streptomyces coelicolor. BMC Genomics, 11, 578. https://doi.org/10.1186/1471-2164-11-578
​

自動ダウンロードの候補 URL 例：

本文 PDF: https://pure.rug.nl/ws/portalfiles/portal/6754510/2010BMCGenomicsCastroMelchor.pdf
​

Supplementary（regulators list を含む XLS）: https://pure.rug.nl/ws/portalfiles/portal/6754508/2010BMCGenomicsCastroMelchorSupp2.pdf や類似 URL（必要に応じてリダイレクト追跡）。
​

可能な限り HTTP 経由でダウンロードを試す。

成功した場合：literature/CastroMelchor2010_* として保存。

失敗した場合：README に必要なファイル名と URL を記載し、ユーザーが手動で置く前提に切り替える。

取得できた Supplementary から、「692 regulators」を含む表をパースし、CastroMelchor2010_regulators.tsv を作成：

SCO_ID

gene_name

regulator_flag（1/0）

source = "CastroMelchor_2010"

3.3 σ因子レビュー（Fernández-Martínez et al. 2018）
APA 形式の参考文献：

Fernández-Martínez, L. T., Barka, E. A., & Bibb, M. J. (2018). Connecting metabolic pathways: Sigma factors in Streptomyces spp. Frontiers in Microbiology, 8, 2546. https://doi.org/10.3389/fmicb.2017.02546
​

PDF を自動ダウンロード可能か試みる（例：https://www.frontiersin.org/journals/microbiology/articles/10.3389/fmicb.2017.02546/pdf）。
​

成功したら literature/FernandezMartinez2018_SigmaFactors.pdf として保存。失敗した場合は README に URL と必要ファイル名を記載。

レビュー本文と図・表から、S. coelicolor の σ 因子名・SCO ID・σ クラスを抽出し、SigmaFactors_M145_from_lit.tsv を作成：

SCO_ID

gene_name（sigA, sigB, sigH, sigR, sigU, …）

sigma_group（sigma70 / ECF / alternative など）

source = "FernandezMartinez_2018"

必要に応じて、Two-component system レビューも同様に扱う：

Fernández-Martínez, L. T., & Liras, P. (2022). Two-component systems of Streptomyces coelicolor. Frontiers in Microbiology, 13, 1022371. https://doi.org/10.3389/fmicb.2022.1022371
​

4. 統合とファミリー分類
4.1 データ統合
M145_gene_basic_info.tsv を基準テーブルとして用いる。
​

SCO_ID をキーに、以下を left join で統合する：

M145_TF_candidates_from_pfam.tsv
​

ZorroAranda2022_regulators.tsv
​

CastroMelchor2010_regulators.tsv
​

SigmaFactors_M145_from_lit.tsv
​

統合結果を master_TF_list_M145.tsv として保存。

4.2 TF ファミリー分類
Pfam/InterPro ドメインと文献情報をもとに、各遺伝子を代表的な TF ファミリーに分類する：

TetR_N/TetR_C → TetR family

GntR_N → GntR family

MarR → MarR family

SARP or HTH_SARP → SARP family

LysR_substrate + LysR_C → LysR family

AraC-like → AraC family

LuxR_C → LuxR family

XRE-family → XRE family

WhiB → WhiB family

Sigma70_r2/4 or ECF_sigma → sigma family（sigma_group で細分類）
​

REC + Trans_reg_C → TCS_response_regulator family（必要なら別フラグ）
​

TF_family カラムに文字列として記録する。複数候補がある場合は、より典型的、または Zorro-Aranda 2022・Fernández-Martínez 2018 などと整合するものを採用。

5. QC と README
master_TF_list_M145.tsv を使用して、以下の簡易 QC を行い、README_master_TF_list_M145.md に追記する。

全 CDS 数 vs pfam_TF_flag == yes の数。
​

ZorroAranda_regulator_type 非 NA の数。
​

sigma_flag == yes の数、および ECF σ の数。
​

TF ファミリー別の数（TetR, GntR, SARP, XRE, etc.）。

README に以下を明記：

使用した RefSeq アセンブリとファイル名。
​

Pfam/InterPro スキャンに用いたツールとバージョン。
​

自動ダウンロードに成功した文献ファイルの一覧（URL と保存パス）。

自動ダウンロードに失敗し、ユーザー手動投入が必要なファイルと URL。

上記 APA 形式の参考文献（Zorro-Aranda 2022, Castro-Melchor 2010, Fernández-Martínez 2018, 必要なら TCS レビュー）。

6. 実装ポリシー
可能な限り Python + pandas で TSV の生成・結合・集計を行う。

ダウンロードには curl / wget / Python の requests 等、利用可能な方法を使ってよい。

各ステップは再実行可能なスクリプトとして保存し、README にスクリプト名を記載する。

ダウンロード失敗時は、必ず

何をどこから取りたいのか（論文名・URL・想定ファイル名）

ユーザーがどのディレクトリに置けば次のステップが動くか
を明示する。

以上を満たすように、順にコード・コマンドを生成・実行し、master_TF_list_M145.tsv と README_master_TF_list_M145.md を完成させてください。