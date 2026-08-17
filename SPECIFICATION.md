# CPM Paired Fastq IgBLAST AIRR tsv Ver3.0 仕様書

実装version: **3.0.0**
仕様日: **2026-08-17**

## 1. 目的

CPM社のpaired FASTQを非マージのR1/R2としてIgBLAST AIRR TSVへ変換し、
同一の注釈結果から次の2単位を別々に出力する。

1. RG版と同じ定義のread-pair support
2. BCR clonotype内で数えた完全一致raw UMI family support

Ver3.0ではIgBLAST前にUMIや全長read配列を用いてpairをcollapseしない。
全保持pairを1回だけIgBLASTへ渡し、BCR clonotypeを決定した後にUMIを数える。

標準portable folder名:

```text
CPM Paired Fastq IgBLAST AIRR tsv Ver3.0
```

GitHub:

```text
YoheiFunakoshi/Cpm-paired-fastq-igblast-airr-tsv-Ver3.0
```

## 2. 解析単位と用語

| 用語 | 定義 |
|---|---|
| raw FASTQ pair | 入力R1/R2の対応する1組 |
| retained pair | FASTQ構造、pair ID、選択readのQCを通り、queryへ書かれたpair |
| integrated pair | IgBLAST後にR1/R2のcallを統合した1行 |
| BCR clonotype | 正規化V候補集合、J候補集合、canonical junction AAの組 |
| raw UMI | CPM R2から抽出して補正せず記録した12塩基文字列 |
| UMI family | 1つのBCR clonotype内で観測されたdistinctな有効raw UMI |
| UMI missing pair | 有効12-mer UMIを得られなかったが保持されたpair |

read pair、UMI family、UMI missing pairは異なる観測単位である。

## 3. 標準処理順

```text
R1/R2 FASTQ
  -> pairを検証
  -> raw R2からUMIを抽出しquery名へ記録（collapseしない）
  -> R1/R2を別queryとしてIgBLAST（1回）
  -> R1/R2 AIRR rowをpair単位で統合
  -> RG版と同じBCR clonotypeを決定
  -> read-pair countsを作成
  -> 各clonotype内で完全一致raw UMIを数え、UMI countsを作成
  -> 全出力成功後だけ正式出力へ公開
  -> completion manifestを最後に公開
```

UMI抽出は物理的にR2を読むためIgBLAST前に行うが、UMI familyのgroupingと
countingはBCR annotation後に行う。

## 4. 入力

必須:

- R1 FASTQ（plainまたはgzip）
- R2 FASTQ（plainまたはgzip）
- `igblastn`
- V/D/J germline DB prefix
- IgBLAST auxiliary data

FASTQは4行1record、sequence/quality同長でなければならない。`.gz`の大文字小文字を
問わずgzipとして読む。R1/R2のrecord数、正規化read ID、取得可能なmate metadataを
検証する。既定ではR1/R2 ID不一致を失敗にする。

## 5. UMI抽出

raw R2先頭を次の構造として扱う。

```text
TATCAACGCAGAGTGGCCAT + NNNN + T + NNNN + T + NNNN + TCTT + insert
```

- anchor: `TATCAACGCAGAGTGGCCAT`
- 既定anchor許容差: Hamming distance 2
- UMI: 3個の4塩基blockを連結した12塩基
- separatorおよび後続`TCTT`: 既定の抽出合否には使用しない

抽出値はquery名へraw値のまま保存する。familyとして有効なのは
`A/C/G/T`だけからなる12-merである。anchor不成立、短すぎるR2、`N`など曖昧塩基を
含む値、長さ不正はUMI missingとして扱う。pair自体は削除しない。

1塩基違いを含む近傍UMI補正は行わない。サンプル全体でUMIを統合しない。

## 6. query作成

既定:

- read selection: `both`
- R1 orientation: `forward`
- R2 orientation: `reverse-complement`
- R1/R2 trim: 0
- minimum length: 0
- maximum N rate: 1.0
- query名: `{read_id}|{read}|UMI={umi}`

R1/R2はマージせず、別々のFASTA recordとして出力する。UMIはraw R2から抽出し、
同じpairのR1/R2 query名へ記録する。

custom query名は下流parserと往復可能でなければならない。`{read_id}`、`{read}`、
CPM解析では`{umi}`を必須とし、format conversion/specifier、予約済みR1/R2 component、
複数の`UMI=` componentなど、曖昧に解釈されるtemplateを拒否する。

## 7. IgBLAST

既定:

- output format: AIRR rearrangement TSV (`-outfmt 19`)
- organism: `human`
- domain system: `imgt`
- sequence type: `Ig`
- threads: 4
- batch size: 10,000 query records

pipeline管理対象の`-query`、`-out`、`-outfmt`、DB、aux、thread等を
`extra args`で上書きすることを拒否する。

非空queryに対し、AIRR data row数がquery record数と一致しない場合は失敗とする。
header-only、batchの部分欠落、重複rowを完了解析として公開しない。

## 8. R1/R2統合

統合規則はRG paired-FASTQ版と同じである。

### 8.1 junction AA

- R1/R2が同値: その値、status=`match`
- 片側のみ値あり: その側、status=`r1_only`または`r2_only`
- 両側が異なり片側だけcanonical: canonical側
- 両側ともcanonicalまたは両側とも非canonical: R1優先、status=`conflict`

canonical junction AA:

- 値あり
- stop `*`なし
- `C`開始
- `W`または`F`終了
- 5～40 amino acids

### 8.2 他のcall

- V call: R2優先、なければR1
- J/D/productive/junction nucleotide: junction AAで採用したreadを優先
- 優先側に値がなければ反対側へfallback

元のR1/R2 call、採用元、判断理由、競合状態を`integrated.tsv`へ保存する。

## 9. BCR clonotype

key:

```text
(unique_v_gene_set, unique_j_gene_set, final_junction_aa)
```

- V/Jはallele suffixを除去する。
- comma区切り候補を重複除去し、文字列昇順で正規化する。
- 候補集合は完全一致で比較し、部分的な共通geneだけでは同じkeyにしない。
- D callは監査用に保持するがkeyに含めない。
- locusはkeyに含めない。
- `final_productive`は基本countsの必須条件にしない。

countsへ含めるにはV集合、J集合、canonical final junction AAが必要である。
条件を満たさないpairも`integrated.tsv`に残し、`include_in_counts=false`と
`exclude_reason`を記録する。

このkeyはBCR全塩基配列の完全一致ではない。同義置換やkey外のSHM差は同じ
clonotypeに入ることがある。

## 10. read-pair counts

`integrated_counts`はRG版と同じ10列で、同じclonotypeに入ったintegrated pairを
1 pairずつ数える。

```text
unique_v_gene_set
unique_j_gene_set
final_junction_aa
read_pair_count
match_count
conflict_count
r1_only_count
r2_only_count
productive_true_count
canonical_junction_aa_count
```

`read_pair_count`はraw FASTQ総pair数ではない。FASTQ/query QC、IgBLAST、統合、
counts条件を通過した保持pair数であり、PCR duplicateを含み得る。

`final_productive_counts`は`final_productive=true`のintegrated pairだけを同じkeyで
再集計した二次確認表である。

## 11. UMI counts

有効raw UMIはclonotypeごとに独立したsetへ入れる。

```text
umi_family_count = count(distinct valid raw UMI within this clonotype)
```

同じUMIが別clonotypeに現れた場合、各clonotypeで1 familyとして数える。
サンプル全体で共有UMIを解決・統合・削除せず、共有UMIフラグも必要としない。

出力列:

```text
unique_v_gene_set
unique_j_gene_set
final_junction_aa
umi_family_count
read_pair_count
umi_known_read_pair_count
umi_missing_read_pair_count
inclusive_support_count
inclusive_support_percent
match_count
conflict_count
r1_only_count
r2_only_count
productive_true_count
canonical_junction_aa_count
```

関係式:

```text
read_pair_count
  = umi_known_read_pair_count + umi_missing_read_pair_count

inclusive_support_count
  = umi_family_count + umi_missing_read_pair_count
```

`inclusive_support_count`はUMI familyと補正不能read pairを加えたhybridな包含支持量で、
厳密な分子数ではない。必ず構成列とともに解釈する。

`inclusive_support_percent`は表内の`inclusive_support_count`合計に対する割合を
小数点以下6桁でTSVへ書く。XLSXでは0～1の数値percent cellとして保存する。

既定sort:

1. `inclusive_support_count`降順
2. `umi_family_count`降順
3. `read_pair_count`降順
4. V、J、junction AA昇順

`final_productive_umi_counts`はproductive対象行からUMI set、missing数、割合を
独立に再計算する。

## 12. ユーザー例

Pattern 1:

```text
A, A, A, B, missing, missing
```

結果:

```text
read_pair_count = 6
umi_family_count = 2
umi_missing_read_pair_count = 2
inclusive_support_count = 4
```

Pattern 2:

```text
C, C, D, E, A, A, missing
```

結果:

```text
read_pair_count = 7
umi_family_count = 4
umi_missing_read_pair_count = 1
inclusive_support_count = 5
```

UMI AはPattern 1とPattern 2のそれぞれで1 familyである。どちらも捨てない。

## 13. 出力

基準出力を`sample.airr.tsv`とした場合:

```text
sample.airr.tsv
sample.R1.airr.tsv
sample.R2.airr.tsv
sample.integrated.tsv
sample.integrated_counts.tsv
sample.integrated_counts.xlsx
sample.final_productive_counts.tsv
sample.final_productive_counts.xlsx
sample.umi_counts.tsv
sample.umi_counts.xlsx
sample.final_productive_umi_counts.tsv
sample.final_productive_umi_counts.xlsx
sample.queries.fasta                 # 保存指定時のみ
sample.run.json
```

Ver3.0は`_umiSeq5` / `_umiNoCollapse`などmode suffixを新規出力名に付けない。
1回のIgBLASTから1組のpair-level/summary出力を作る。

確認順:

1. `integrated_counts.xlsx`: RG互換read-pair比較
2. `umi_counts.xlsx`: CPMのUMI familyとmissing支持
3. `integrated.tsv`: pair単位の採用根拠・除外理由
4. `final_productive_*`: productive限定感度確認

## 14. manifest

`sample.run.json`はcompletion markerであり、全データ出力の後に最後に公開する。

必須識別情報:

- `manifest_schema_version = 2`
- `software_version = 3.0.0`
- `counting_semantics = cpm_v3_read_pair_and_exact_raw_umi_per_clonotype_v1`
- 開始時入力path/size/mtime
- query/QC/UMI抽出/IgBLAST設定
- prepare/pair-summary統計
- 全正式出力pathとsize

入力FASTQ内容やhashはmanifestへ保存しない。manifest不在は未完了として扱う。

## 15. 出力安全性

- R1/R2同一pathを拒否する。
- 入力、query、全出力、IgBLAST executable/DB/auxのpath衝突を拒否する。
- Windowsのcase-insensitive alias、relative/absolute alias、legacy suffix aliasも保護する。
- 既存結果は既定で上書きしない。
- GUIの上書きは明示確認、CLIは`--overwrite`必須。
- query、AIRR、派生TSV/XLSX、manifestをscratchへ作成する。
- 全処理成功後だけ正式出力へ順次replaceし、失敗時は旧一式へrollbackする。
- completion manifestを最後にcommitする。
- 出力ごとのlockで同一結果への同時実行を拒否する。
- 入力FASTQのresolved path、size、mtime、file identityを開始時に記録し、解析中と
  publish直前に不変を確認する。
- IgBLAST resourceはsource identity/role/content fingerprint別のimmutable bundleへ
  SHA-256検証付きでstageする。
- scratchは成功・失敗後に削除する。

## 16. GUI/portable

GUIは必須入力、R1/R2同一、出力衝突、IgBLAST、V/D/J DB component、aux、数値範囲を
worker開始前に一覧検証する。解析中はRunを無効化し、ウィンドウ終了を拒否する。

portable標準起動:

```text
Open CPM Paired Fastq IgBLAST AIRR tsv Ver3.0.cmd
```

portable配布にはPython、IgBLAST、refdataを含めてよいが、GitHubへbinary/DBや
研究FASTQ/実結果をcommitしない。

## 17. Ver2.0との非互換

Ver2.0の`umiSeq5`は、IgBLAST前に補正UMIとR1/R2全長Hamming距離を用いて
sequence clusterを残した。Ver3.0は全pairを注釈し、clonotype内exact raw UMIを
数える。UMI family数、read-pair数、順位、ファイル名、manifest semanticsは
数値互換ではない。

旧結果と統合・比較する場合は、元FASTQをVer3.0で再解析するか、Ver2/Ver3を
別methodとして明示する。

## 18. 方法上の限界

Ver3.0はannotation-first, clonotype-aware UMI countingであり、UMI consensus caller
ではない。

- annotation errorがV/J/junction AAを変えると1真分子を複数clonotypeへ分け得る。
- 同一UMI・同一clonotypeの真のUMI collisionは観測情報だけでは分離できない。
- UMI missing pairは分子補正できない。
- exact raw UMIはUMI sequencing errorを補正しない。
- read-pair表はPCR duplicateを含む。

これらをpair-level audit tableと分離列で可視化する。将来consensusや近傍UMI補正を
追加する場合は、Ver3.0既定値を黙って変えず、別名・別versionの方式として検証する。

## 19. 受入試験

最低限、次を自動試験する。

- 全pairについてprepare/IgBLASTが1回だけ実行される。
- RG互換countsが10列である。
- 同一clonotype内の同一UMI反復が1 familyになる。
- 1塩基違いUMIが別familyになる。
- 同じUMIが別clonotypeでそれぞれ数えられる。
- UMI missing/曖昧UMIのpairが削除されず別計数される。
- counts外BCR rowが`integrated.tsv`に残る。
- productive限定表が対象行から独立再計算される。
- XLSX percentが数値cellである。
- AIRR row部分欠落を公開しない。
- input/output/resource衝突を拒否する。
- publish途中失敗・Ctrl+C・rollback異常で偽completion manifestを残さない。
- 入力が解析中に置換された場合、正式出力を公開しない。
