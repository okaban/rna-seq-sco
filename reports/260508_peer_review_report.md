# 査読レポート

**掲載誌想定:** *Nature Microbiology* / *Nucleic Acids Research*  
**原稿タイトル（作業題）:** Dual DNA methylation systems shape *Streptomyces coelicolor* regulatory genome architecture without direct transcriptional control  
**査読日:** 2026年5月8日  
**査読担当:** (peer review simulation)

---

## Summary（要旨評価）

本論文は、*Streptomyces coelicolor* A3(2) M145 において Oxford Nanopore 直接塩基配列決定法と RNA-seq を三時点にわたって統合し、二種類の DNA メチル化システム（GCCGGC [4mC]、AAGCCCG [6mA]）が転写制御に直接関与しているか否かを体系的に検証している。最大の貢献は、表向き有意な「メチル化−発現相関」がいずれも線状染色体の core/arm 機能分化に起因する Simpson's paradox の人工物であることを示した点であり、この警告は細菌エピゲノミクス分野全体にとって重要な方法論的寄与である。また、調節遺伝子プロモーター周辺の約 2,200 bp の methylation-free protection zone と、それを欠く 57 遺伝子の構造的二分法（AUC = 0.917）は新規かつ興味深い観察である。ただし、因果言語の過剰使用、統計的に支持されない知見の主張、実験的検証の完全な欠如という三点が現在の状態で主要誌への掲載を妨げる。

---

## Major Concerns（主要懸念事項）

### [MC-1] 因果言語と観察データの乖離——最大の論理的矛盾

本稿の最も深刻な問題は、因果的主張と自身が報告した証拠の内部矛盾である。

**問題の核心：** Results 終盤（"Negative results define the model's boundaries"）において、57 個の exposed regulators における「メチル化変化量と発現変化量の相関は T1–T2 で *rho* = 0.136, *p* = 0.299、T2–T3 で *rho* = 0.012, *p* = 0.927」と報告している。これは methylation dynamics が exposed TF の expression dynamics を駆動していないことを著者自身が認めるデータである。にもかかわらず、同じ論文の別箇所に以下の因果的記述が残存している：

- Abstract: *"channels methylation effects through a specific minority of regulators"*
- Abstract・Results: *"constitute a distributed vegetative-to-developmental switch"*
- Results: *"the 57 exposed regulators are the sole conduit through which methylation interacts with the transcriptional program"*
- Discussion: *"methylation's regulatory influence is thus channeled through a narrow bottleneck"*

「メチル化変化と発現変化が無相関」という自身の結果が、「メチル化が調節的影響をこれらの遺伝子を通じて伝達する」という主張を直接否定している。この内部矛盾を解消するには、因果フレームを「構造的共在（structural co-occurrence）」に全面的に改訂する必要がある。現行原稿はこの修正を一部しか行っておらず、査読者に対して論文の主張が何であるかを根本的に曖昧にしている。

**要求事項:** Abstract・Introduction・Results・Discussion の全文を精査し、因果的語彙（"channels"、"sole conduit"、"switch"）を中立的な共在・相関フレームに統一すること。

---

### [MC-2] TCS asymmetry の過大評価（p = 0.83 の結果を "functionally significant" と記述）

Results セクション（"57 exposed regulators form a distributed developmental switch"）は、cognate TCS pair において exposed/shielded 非対称性が 7/7 で観察されたことを取り上げ、*"functionally significant"*、*"probability of this perfect asymmetry occurring by chance is less than 10^-3"* と記述している。

しかし、統計検定登録 (T27) には以下が明記されている：

> "p = 0.83 (Permutation, label shuffle, n=10,000) — **descriptive observation only**"

この p 値は、84 pairs 中 7 pairs に exposed メンバーがいるという背景頻度（~8.3%）のもとでは、0 co-exposed pairs は完全に偶然と一致することを意味する。"less than 10^-3" という記述は二項検定の計算結果と思われるが、それは適切な帰無仮説（84 pairs のうち約 8.3% が exposed という経験的頻度を前提とする置換検定）に基づいていない。

p = 0.83 の結果を "functionally significant" と表現し、TCS gating の機構的考察を 1 段落以上費やして展開することは、査読者の科学的信頼を損なう。この観察は Discussion の parenthetical note として p = 0.83 を明示した上で記述するにとどめるべきである。

---

### [MC-3] 33%/67% の apportionment 計算の不整合

Results では「sequence features が protection の約 33%、protein occupancy が約 67% を説明する」と記述し、この数値を論文全体にわたって繰り返し引用している。しかし提示されたデータからはこの分割を導出できない。

著者が提示しているのは：
- Sequence-only logistic regression の 5-fold CV AUC = **0.712 ± 0.049**
- Nearest methylation distance の AUC = **0.917**

AUC による分解を正しく行うと：
- Sequence 寄与（偶然より上の部分）= 0.712 − 0.500 = 0.212
- Total discriminative power above chance = 0.917 − 0.500 = 0.417
- Sequence 割合 = 0.212 / 0.417 ≈ **51%**
- Protein 割合 = (0.917 − 0.712) / 0.417 ≈ **49%**

したがって AUC 分解では 33/67 ではなく 51/49 となり、「~33%/~67%」という記述は誤りである。著者が 33/67 を導出した根拠（sequence motif depletion fold = 0.74–0.76 に基づく別途計算の可能性）は本文中に説明されておらず、読者は再現不能である。2 通りの計算法が異なる答えを与えることを明示し、それぞれの数値と解釈を区別して提示しなければならない。

---

### [MC-4] 実験的検証の完全な欠如

本稿は *Nature Microbiology* レベルの主要誌を想定しているが、主要な主張のいずれも遺伝学的介入実験（MTase ノックアウト、TF ノックアウト、ChIP-seq、DAP-seq）によって検証されていない。

著者が Gatekeeper model から導く最も重要で検証可能な予測は：「MTase 欠損株では、exposed regulators の発現のみが変化し、shielded regulators の発現は変化しない」というものである。この実験なしには、observed protection zone と developmental transition の共在が因果関係ではなく単なる epiphenomenon である可能性を排除できない。

Limitations セクションでこの点を自認していることは評価できるが、完全に correlational な証拠のみで "four-layer architecture" の機能的モデルを構築・主張することは、主要誌の水準では受け入れがたい。少なくとも、最も直接的な MTase 欠損実験（*SC_RS17645* または GCCGGC MTase 遺伝子の欠損）のいずれかを実施するか、その欠如を Discussion でより正面から論じる必要がある。

---

### [MC-5] BGC 解析における Simpson's paradox の適用と残存分析の矛盾

Results の "Negative results" セクションでは、GCCGGC メチル化の BGC 集積（fold = 1.66、*p* = 4.4 × 10^-8）は core/arm 共局在による Simpson's paradox 人工物であり、core 内解析では fold = 1.07（NS）に消失すると正しく記述している。

ところが直後の段落（"Within the BGC compartment"）では、major antibiotic clusters（ACT、CDA、RED、CPK）が minor BGCs より約 1.9 倍高い T1 GCCGGC 密度を示す（Mann-Whitney *p* = 0.062）と主張し、「prepare-then-release model」を提唱している。しかしこの内部比較に対しても geographic stratification が適用されているか？ Major clusters（ACT: ~6.0 Mb、CDA: ~3.3 Mb、RED: ~5.0 Mb、CPK: ~7.4 Mb）はすべて core 内であるが、minor BGCs は core・arm を含む。core-only 解析での BGC enrichment が消失したことを踏まえると、「major BGC の 4mC 濃縮」も単に major BGCs が core の特定領域に集中することを反映している可能性があり、この仮説はさらなる地理的層別化を必要とする。自ら Simpson's paradox を発見した著者が同一論文でその陥穽に再び陥ることは、論理的一貫性を著しく損なう。

---

### [MC-6] Protection zone 計算における TSS 位置の信頼性

Methods によれば、TSS は「RefSeq GFF アノテーションの gene start coordinate を一次ソースとし、dRNA-seq データが利用可能な 257 遺伝子についてのみ実験的 TSS を使用する」とある。つまり 1,055 の調節遺伝子のうち約 75% では、アノテーションの開始コドン位置を TSS の代理として使用している。

細菌では 5' UTR 長が 0–300 bp 程度変動するため、アノテーション start と実際の TSS の乖離は "nearest methylation distance" 計算に直接影響する。特に AUC = 0.917 の分類と 293 bp という最適閾値は、TSS 位置の精度に対して感受性が高い。実験的 TSS の有無で exposed/shielded 分類の結果がどう変わるかの感受性解析、あるいは Supplementary Note に記載された Pfam domain-based replication（AUC = 0.923、318 bp boundary）との詳細比較が求められる。

---

## Minor Concerns（軽微な懸念事項）

### [mn-1] 数値の不整合

- **Exposed regulators の割合**: Abstract では「5.4%」、Results Introduction では「5.9%」と記述が異なる。57/1,055 = 5.405% なので「5.9%」は誤りである。全文を検索して統一すること。

- **"four modules covering 61 of 57 genes"**: Results の co-expression 解析で「61 of 57 genes」と記述されているが、61 > 57 は算術的に不可能である。統計検定登録 (T26) の注記には「4 modules covering 56/57 genes」とある。正しい数値に修正すること。

- **Protection zone の幅**: Abstract では「~2,200 bp」、Results では「-1,300 to +700 bp relative to TSS」（= 2,000 bp）と不一致がある。いずれかに統一すること。

### [mn-2] "Gatekeeper" 術語の残存

本原稿は全文を通じて "Gatekeeper architecture" という命名を繰り返し使用しており（Abstract、Introduction、Discussion 等）、これが一種の因果的含意を持つことは Reviewer-Response-Draft でも著者自身が認識している。中立的な記述語（"four-layer methylation-exclusion model" 等）への全面置換が必要である。

### [mn-3] SC_RS17645 の分類の曖昧さ

Methods および Results では SC_RS17645 を "HsdM-type methyltransferase" と記述しているが、これは "Type I R-M system の一部" という誤解を招く。REBASE および周辺遺伝子の確認から、この MTase は HsdR/HsdS パートナーを持たない孤立型（orphan MTase）である。"HsdM-type fold を持つ orphan N6-adenine methyltransferase" と明確に記述すること。

### [mn-4] STREME 解析結果との整合性

Results では「72.2% の 6mA sites が canonical AAGCCCG に帰属せず、少なくとも 1 つの未知 MTase が存在する可能性を示唆する」と記述している。しかし STREME による de novo motif discovery（Writing/STREME_Results_paragraph.md）では、unassigned sites に対して高信頼の新規モチーフは同定されなかった（最良候補 Motif 1 の E-value = 3.0 × 10^-3 は canonical AAGCCCG の E ~ 10^-22 と比較して数桁劣る）。したがって現行原稿の「追加 MTase の示唆」は STREME 結果によって支持されておらず、「unassigned 6mA sites に対するモチーフ探索を行ったが高信頼の新規認識配列は見出されなかった」と明示すべきである。

### [mn-5] 統計的「同時性」の主張方法

Results （"Both blocs activate simultaneously"）は Mann-Whitney U test の p = 0.459（phase ratio 比較）を根拠に「両 bloc が同時に engage する」と結論する。しかし帰無仮説棄却失敗（p > 0.05）は同時性の証明ではない。真の等価性主張には TOST（two one-sided test）等の equivalence testing が必要である。現行の記述を「統計的には有意な時間的分離は観察されなかった」と書き換えるか、等価性検定を追加すること。

### [mn-6] Citation プレースホルダーの残存

本文に `[CITE: ...]` 形式のプレースホルダーが複数残存している（Introduction: Tock and Dryden 2005、Casadesus and Low 2006 等）。Methods にも `[repository URL]` が未記入のまま残っている。投稿前に全て実際の参考文献番号または URL に置換すること。

### [mn-7] Double comma の校正漏れ

Results 第 2 セクション冒頭（"Having established that methylation does not directly regulate transcription genome-wide,,"）に二重コンマが存在する。

### [mn-8] Methods における Regulatory gene annotation の重複記述

Methods の "Regulatory gene annotation" サブセクションでは、keyword-based TF 選定方法がほぼ同一の内容で 2 回記述されている（195–197 行目と 198–199 行目）。冗長な記述を整理すること。

### [mn-9] KEGG Quorum Sensing enrichment の生物学的文脈

Results および Discussion において、KEGG Quorum Sensing pathway（ko02024）への有意な enrichment が文脈説明なく記述されている。*S. coelicolor* は conventional N-acyl-HSL quorum sensing を持たず、この enrichment は KEGG の broad annotation カテゴリに起因することを本文中で明示すること。

### [mn-10] 「タイムポイント3点のみ」設計の解釈への影響

Limitations では T1–T3 の 3 時点設計が制限として言及されているが、この点が「57 遺伝子が同時に T1→T2 で応答する」という主張に対してどの程度影響するかが議論されていない。T1 と T2 の間に複数の中間時点が存在すれば、"simultaneous switch" は sequential cascade に分解される可能性があり、これについての考察が必要である。

---

## Specific Line-Level Comments（箇所別コメント）

| 箇所 | 内容 |
|------|------|
| **Abstract, 2段落目** | "constitute a distributed vegetative-to-developmental switch" → 因果的含意を除いた中立語に改訂 |
| **Abstract, 3段落目** | "channels methylation effects through" → このフレーズはデータによって否定されているため削除。"structurally co-localizes with" 等に変更 |
| **Introduction, 最終段落** | "12 rejected hypotheses about direct methylation-transcription control were required to reveal the true, structurally mediated pathway through 57 specific regulatory genes" — "true pathway" という表現は結論を先取りしており、現在の observational な証拠水準を超えている |
| **Results, R-M systems セクション** | "Notably, initial classification suggested possible 5-methylcytosine (5mC), but detailed analysis confirmed that the modification is 4mC" — この m4C/5mC 再帰属についての議論は現在この一文だけで完結しており不十分。Pisciotta et al. (2023) との関係を Results ではなく Discussion に集約し、strand-specific hemi-methylation の証拠（方法論的根拠含む）を詳述すること |
| **Results, AAGCCCG セクション** | "HsdM-type MTase" → "orphan N6-adenine methyltransferase (HsdM-type fold; no cognate HsdR/HsdS)" と訂正 |
| **Results, TCS pair セクション** | "The probability of this perfect asymmetry occurring by chance is less than 10^-3" → 置換検定 p = 0.83 と矛盾。削除または徹底的に訂正すること |
| **Results, Negative results セクション** | "major antibiotic clusters showed approximately 1.9-fold higher T1 GCCGGC 4mC density (p = 0.062)" — この p 値は名目上の有意水準（0.05）を超えている。BGC subset 解析での core/arm 層別化が実施されたか否かを明記すること |
| **Discussion, Reframing セクション** | "the entire apparent genome-wide methylation-expression correlation is a statistical artifact" — "entire" は強すぎる。AAGCCCG に対しても geographic confound が存在する可能性を排除できているか？ |
| **Discussion, 57 exposed regulators セクション** | "operates independently of all computationally characterized transcriptional networks" — "all" は証明不可能な全称命題。"is not captured by any of the computationally characterized networks examined here" に改訂 |
| **Methods, Regulatory gene annotation** | keyword-matching による TF 選定の false positive/negative 率について、Supplementary Note で実施されている Pfam-based robustness check の結果を Methods 本文でも参照・要約すること |

---

## Overall Recommendation（総合評価）

**Major Revision**

本稿の核心的な科学的貢献——Simpson's paradox の識別、protection zone の発見、and Shielded/Exposed 二分法の定量化——は genuine であり、方法論的な厳密さと網羅性は際立っている。しかし、現時点での掲載を推薦できない理由は以下の一点に集約される：**因果的主張と observational データの整合性がとれていない**。論文は「メチル化が 57 遺伝子を通じて転写プログラムを制御する」という因果モデルを主張しながら、同一論文内でその因果連鎖を否定するデータ（rho = 0.136, p = 0.299 の temporal decoupling）を提示している。この矛盾は Reviewer-Response-Draft を確認する限り著者自身も認識しているが、manuscript 本文への修正が完了していない。

加えて、TCS asymmetry の過大評価（p = 0.83 を "functionally significant" と記述）、33/67 apportionment の計算エラー、BGC 内部比較での geographic control の欠如、実験的検証の完全な不在が加わり、全体として主要誌に要求される論理的整合性に達していない。

上記 Major Concerns [MC-1]〜[MC-6] の対処、および [mn-1]〜[mn-10] の軽微修正の後に再査読を求める。特に [MC-1] の因果フレームの全面改訂と [MC-2] の TCS asymmetry 記述の撤回は不可欠であり、これらなしに掲載を推薦することはできない。

---

*本査読は `full_manuscript.md`、`Supplementary_statistical_tests.tsv`、`STREME_Results_paragraph.md`、および `Reviewer-Response-Draft.md`（コンテクスト参照のみ）に基づいて作成した。*
