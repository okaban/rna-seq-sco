# SRA 提出ガイド（2026-09-10 作成）

すべてのメタデータは埋まっている。**アカウント操作だけが残っている**。
所要は入力 30 分 ＋ アップロード（86.7 GB、回線次第で数時間〜1 日）。

## 0. 先に確認する 1 件 — pod5 を SRA が受け付けるか

SRA の公式フォーマットガイドは Nanopore ネイティブ形式として **fast5 しか挙げていない**
（pod5 の記載なし、2026-09 時点）。72 GB を上げてから弾かれるのを避けるため、
**先にヘルプデスクへ問い合わせる**。下書き: `sra_helpdesk_email.txt`（送信先 sra@ncbi.nlm.nih.gov）。

返答までの間、BAM と FASTQ（14.3 GB）の提出は先に進めてよい。pod5 は同じ
BioProject に後から run を追加できる。ダメと言われた場合の代替は
**Zenodo に pod5 をアーカイブ**（無料枠 50 GB/レコード。72 GB は 2 レコードに分割か、
上限引き上げをリクエスト）。

## 1. 提出の構造

```
BioProject（1 件）
 └ BioSample × 9（同一培養から DNA/RNA を分取したので 9 個で両アッセイを共有）
    ├ SRA run: ONT BAM（9）
    ├ SRA run: Illumina FASTQ ペア（9）
    └ SRA run: pod5（9、ヘルプデスク OK なら）
```

## 2. 手順

1. https://submit.ncbi.nlm.nih.gov/ にログイン（NCBI アカウント。ORCID 連携可）
2. **BioProject** → New submission
   - Project type: Raw sequence reads
   - Target: `Streptomyces coelicolor A3(2)`、Sample scope: Multispecies ではなく **Monoisolate**
   - Title / Description は `.zenodo.json` の title/description を流用してよい
   - Release date: **論文 acceptance まで hold**（"Release on specified date" か "upon publication"）
3. **BioSample** → New submission → Package: **Microbe; version 1.0**
   - `biosample_attributes.tsv` をそのままアップロード（9 行、必須項目に空欄なし）
4. **SRA** → New submission → 上記 BioProject / BioSample を紐付け
   - `sra_metadata.tsv` をアップロード（18 行）
   - pod5 を出す場合は同シートに 9 行追加（`pod5_upload_manifest.tsv` を参照。
     ディレクトリごと tar にまとめる場合はファイル名を `<sample>_pod5.tar` にし、
     ヘルプデスクの指示に従って個別ファイル名の列挙要否を決める）
5. **ファイル転送** — Aspera（推奨、`ascp`）または FTP。SRA 画面に出る
   一時ディレクトリのパスと鍵を使う。転送対象:
   - ONT BAM 9 件: `ont_bam_md5.tsv` のパス列（6.18 GB）
   - Illumina FASTQ 18 件: `raw_checksum_report.tsv` のファイル名。実体は
     `~/Library/CloudStorage/Dropbox-SFC-CNS/Takeda Tomoki/my_projects/2026/M145_RNA-seq/data/RH250715199/01_RawData/`（8.13 GB）
   - pod5: `pod5_upload_manifest.tsv` の 17 ディレクトリ（72.42 GB）
6. 転送後、SRA 側の md5 と手元（`ont_bam_md5.tsv`、`raw_checksum_report.tsv`）を照合
7. 査読者に見せる場合: 提出完了後 SRA から **reviewer link** を発行できる

## 3. 注意点

- **organism は必ず *Streptomyces coelicolor* A3(2)**。受託先報告書の生物種欄
  （*S. albidoflavus* J1074）は転記ミスで、アラインメント率 98.4–98.7% が根拠
- ONT run の `library_strategy` は OTHER / `library_selection` RANDOM / `library_layout` single
- BAM の Reference assembly 欄は `GCF_000203835.1`（NCBI assembly なので fasta 添付不要）
- replicate 番号は元のまま（T1: 1,2,3 / T2: 1,3,4 / T3: 2,3,4）。振り直さない
- アクセッションが出たら:
  `python deposition/apply_manuscript_updates.py --accessions PRJNAxxxxxx --doi 10.5281/zenodo.xxxxxxx --write`
  → `make_bootstrap.py` 再実行 → PASS 確認。
  `.zenodo.json` の related_identifiers に BioProject も追加してリリースを更新
