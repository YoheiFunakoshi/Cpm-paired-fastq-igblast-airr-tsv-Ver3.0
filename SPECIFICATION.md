# CPM Paired Fastq IgBLAST AIRR tsv Ver3.0 仕様書

実装version: **3.0.1**
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
| inclusive UMI表 | UMI familyとUMI missing pairを分離列で残す包含集計 |
| exact UMI family表 | 有効UMI familyだけをstrict family countと割合分母にする集計 |

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

### 8.3 Productive

R1/R2の`productive`はIgBLAST AIRR出力を保持し、Ver3.0で独自に再判定しない。
値は`true`、`false`、空欄を区別し、空欄を`false`として扱わない。AIRR TSVでは
通常`T`、`F`、空欄として表現される。

`final_productive`はR1/R2のAND判定でもOR判定でもない。8.1でjunction AAを採用したreadを
`preferred_read`とし、次の規則で決める。

1. R1/R2の両側にproductive値があり同値なら、その値とsource=`both`を採用する。
2. `preferred_read`に値があれば、その値を採用する。
3. `preferred_read=both`で値が異なる場合はR1を優先する。
4. 優先側が空欄ならR2、次にR1の順で値のある側へfallbackする。
5. 両側とも空欄なら`final_productive`も空欄とする。

したがって、片側が`true`で反対側が`false`または空欄でも、`true`側が採用されれば
`final_productive=true`になり得る。R1/R2は同一分子の部分配列であり、両側で
productive判定を得ることをBCRの生物学的productive条件とはしない。

IgBLASTのproductiveは、V(D)J配列がin-frame junction、内部stop codonなし、
内部V frame shiftなし等からタンパク質をコードできると予測されたことを表す。
実際のタンパク質発現を直接証明する値ではない。`final_junction_aa`は保存CとW/Fを
含むAIRR JUNCTIONのアミノ酸配列であり、BCR全タンパク質配列ではない。
`complete_vdj=true`はproductive表の追加条件にしない。

参照: [AIRR Rearrangement Schema](https://docs.airr-community.org/en/stable/datarep/rearrangements.html)、
[NCBI IgBLAST productive update](https://blast.ncbi.nlm.nih.gov/doc/blast-news/2021-BLAST-News.html)

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
再集計したproductive限定の派生集計表である。

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

UMI missingは処理段階と分母が異なる次の値を混同しない。

```text
全入力pairのUMI未割当率（prepare段階）
  = umi_missing_pairs / total_pairs

集計表内のpair単位UMI missing率
  = sum(umi_missing_read_pair_count) / sum(read_pair_count)

inclusive support内のUMI missing割合
  = sum(umi_missing_read_pair_count) / sum(inclusive_support_count)
```

検証runのproductive対象では、`read_pair_count=413,106`、
`umi_known_read_pair_count=348,928`、`umi_missing_read_pair_count=64,178`で、
pair単位のUMI missing率は15.535%であった。有効UMI pairをexact UMIで数えると
`umi_family_count=74,968`、`inclusive_support_count=139,146`となるため、
inclusive support内のUMI missing割合は46.123%となった。これはUMI抽出失敗率ではなく、
有効UMI側だけがclonotype内exact UMIで約4.65 pair/familyへ縮約され、PCR重複の影響を
抑えることで生じる分母効果である。

同runの全入力pairでは96,030 / 644,824 = 14.892%がUMI missingであった。
raw FASTQを別途監査した結果、全R2は301塩基で、anchor通過後に非ACGT UMIとなった
pairは0であり、全missing pairは
期待anchorとの距離が既定許容2を超えたことによる。うち89,481 / 96,030 = 93.18%は
6塩基以上異なった。したがって「UMI塩基を読めなかった」ではなく、「期待anchorを
確認できずUMIを安全に割り当てなかった」と記述する。単純な1～2塩基の読み誤りだけで
大半を説明せず、library構造、位置ずれ、off-targetなども候補とするが原因は断定しない。

この実測値は1 runのQC例であり、全protocol共通の正常範囲を定義しない。同一protocolの
検体間でpair単位のmissing率を比較し、大きな差がある場合はlibrary構造、anchor設定、
read向き、品質を再確認する。missing率が異なる検体間では、hybridな
`inclusive_support_count`の比較にbiasが入り得る。

`inclusive_support_percent`は表内の`inclusive_support_count`合計に対する割合を
小数点以下6桁でTSVへ書く。XLSXでは0～1の数値percent cellとして保存する。

既定sort:

1. `inclusive_support_count`降順
2. `umi_family_count`降順
3. `read_pair_count`降順
4. V、J、junction AA昇順

`final_productive_umi_counts`はproductive対象行からUMI set、missing数、割合を
独立に再計算する。

### 11.1 Exact UMI family counts

包含UMI表を変更せず、同じ集計からstrictなexact UMI family表を追加する。

収載条件:

```text
umi_family_count > 0
```

出力列:

```text
unique_v_gene_set
unique_j_gene_set
final_junction_aa
umi_family_count
umi_family_percent
```

割合:

```text
umi_family_percent
  = umi_family_count / sum(umi_family_count in the same exact-UMI-family table)
```

UMI missing pairは`umi_family_count`にも割合分母にも加えない。UMI missingだけで
支持されたclonotypeはexact UMI family表には収載しない。ただし元pairを
`integrated.tsv`や包含UMI表から削除せず、既存の`umi_counts`と
`final_productive_umi_counts`も変更しない。

`final_productive_exact_umi_family_counts`は、productive対象pairからUMI setを
独立に再計算した後、`umi_family_count > 0`条件を適用する。通常exact表の単純な
行filterではない。

既定sort:

1. `umi_family_count`降順
2. V、J、junction AA昇順

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

Pattern 3:

```text
missing
```

結果:

```text
read_pair_count = 1
umi_family_count = 0
umi_missing_read_pair_count = 1
inclusive_support_count = 1
```

Pattern 1～3を互いに異なるBCR clonotypeとする。UMI AはPattern 1とPattern 2の
それぞれで1 familyであり、clonotypeを越えて統合・相殺・共有UMI処理をしない。

全pairが`final_productive=true`の場合、包含UMI表とstrict UMI表は次の値になる。

包含UMI表:

| BCR clonotype | `umi_family_count` | `umi_missing_read_pair_count` | `inclusive_support_count` | `inclusive_support_percent` |
|---|---:|---:|---:|---:|
| Pattern 1 | 2 | 2 | 4 | 40.000000% |
| Pattern 2 | 4 | 1 | 5 | 50.000000% |
| Pattern 3 | 0 | 1 | 1 | 10.000000% |

strict UMI表:

| BCR clonotype | `umi_family_count` | `umi_family_percent` |
|---|---:|---:|
| Pattern 1 | 2 | 33.333333% |
| Pattern 2 | 4 | 66.666667% |

包含UMI表の分母は`4 + 5 + 1 = 10`、strict UMI表の分母は`2 + 4 = 6`である。
Pattern 3はUMI missingだけで支持されるため、包含UMI表では1支持として残り、
strict UMI表では`umi_family_count = 0`のため収載しない。この例ではstrict表の
BCR clonotype行数は3から2へ1つ減る。元pairは`integrated.tsv`および包含UMI表から
削除しない。

この例では全pairがproductiveなので、productive限定なしの対応するUMI表にも同じ値が
出る。productive/nonproductiveが混在する場合は、productive 3表を
`final_productive=true`のpairだけから独立に再集計する。

count可能なV/J/canonical junction AA keyを作れないpairは、除外理由付きで
`integrated.tsv`に保持するがcounts Excelへ収載しない。

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
sample.exact_umi_family_counts.tsv
sample.exact_umi_family_counts.xlsx
sample.final_productive_exact_umi_family_counts.tsv
sample.final_productive_exact_umi_family_counts.xlsx
sample.queries.fasta                 # 保存指定時のみ
sample.run.json
```

Ver3.0は`_umiSeq5` / `_umiNoCollapse`などmode suffixを新規出力名に付けない。
1回のIgBLASTから1組のpair-level/summary出力を作る。

技術QC・監査の確認順（主解析の優先順位ではない）:

1. `integrated_counts.xlsx`: RG互換read-pair比較
2. `umi_counts.xlsx`: CPMのUMI familyとmissing支持
3. `exact_umi_family_counts.xlsx`: valid exact UMI familyだけのstrict支持と割合
4. `integrated.tsv`: pair単位の採用根拠・除外理由
5. `final_productive_*`: productiveな再構成候補を対象とする派生集計

6つのExcelは同じintegrated pairから作る。

| Excel | 列数 | 収載する主な支持量 | 収載条件 |
|---|---:|---|---|
| `integrated_counts.xlsx` | 10 | read pair | V集合、J集合、canonical final junction AAあり |
| `final_productive_counts.xlsx` | 10 | read pair | 上記かつ`final_productive=true` |
| `umi_counts.xlsx` | 15 | exact raw UMI family、UMI missing pair、read pair | 基本countsと同じ |
| `final_productive_umi_counts.xlsx` | 15 | exact raw UMI family、UMI missing pair、read pair | `final_productive=true` |
| `exact_umi_family_counts.xlsx` | 5 | exact raw UMI family、同表family割合 | productive限定なしで`umi_family_count > 0` |
| `final_productive_exact_umi_family_counts.xlsx` | 5 | productive対象のexact raw UMI family、同表family割合 | productive対象で`umi_family_count > 0` |

本解析でproductiveなレパトアを主目的とし、UMI missingのclonotypeも捨てない場合、
主解析表と大小・割合の指標を次のように定める。

| 目的 | 使用するExcel | 大小 | 割合 |
|---|---|---|---|
| CPM productive主解析（捨てなし） | `final_productive_umi_counts.xlsx` | 包含hybrid支持量として`inclusive_support_count` | `inclusive_support_percent` |
| productive strict感度解析 | `final_productive_exact_umi_family_counts.xlsx` | `umi_family_count` | `umi_family_percent` |
| RGとのproductive read-pair比較 | `final_productive_counts.xlsx` | `read_pair_count` | 必要時に同表合計から別途計算 |

`inclusive_support_count`はexact UMI familyとUMI missing pairを加えたhybrid支持量であり、
厳密な元分子数ではない。試料間比較ではUMI missing率を併記し、strict感度解析で
結論の安定性を確認する。QASAS等で包含表の一致clonotype数と割合を主解析とする場合、
両方を包含表から次のように算出する。

```text
一致clonotype数 = 一致したclonotype行の数
一致inclusive support割合
  = sum(inclusive_support_count of matched clonotype rows)
    / sum(inclusive_support_count of all clonotype rows) * 100
```

一致clonotype数は種類数、一致inclusive support割合は包含支持量の割合であり、
同じ意味の数値ではない。strict表の`umi_family_percent`は別名の感度解析として
併記できるが、包含表のclonotype数とstrict表の割合を同じ単位の1結果として混ぜない。

`final_productive_umi_counts`は通常UMI表の行を単純削除するのではなく、productive対象pair
だけからUMI set、missing数、inclusive support、percentを独立に再計算する。
UMI family全体のconsensusをproductive判定する処理ではない。同じclonotype・同じUMIに
productive/nonproductive pairが混在する場合、productive pairが1件以上あれば、そのUMIは
productive subset内で1 familyとして数える。

exact UMI family表では、UMI missing pairはstrict family count・割合・分母へ入らない。
`exact_umi_family_counts`のproductive限定なしは`final_productive=false`限定ではなく、
通常counts対象のtrue/false/空欄を含む。
包含表のclonotype数とexact表の`umi_family_percent`をQASAS等で併記する場合は、
前者をUMI missing-onlyも残す検出指標、後者を有効exact UMI familyだけの支持割合と
明記し、1つの共通単位として混ぜない。strictな大小関係には`umi_family_count`または
`umi_family_percent`を使う。

productive限定を最終評価に使う場合も、通常表との減少を次の同単位で確認する。

- clonotype行数
- `read_pair_count`合計
- `umi_family_count`合計
- `umi_missing_read_pair_count`合計
- `inclusive_support_count`合計

clonotype行数の減少率と支持量の減少率は同義ではない。低支持singleton clonotypeが
除外されると、clonotype行数の方が大きく減り得る。

## 14. manifest

`sample.run.json`はcompletion markerであり、全データ出力の後に最後に公開する。

必須識別情報:

- `manifest_schema_version = 2`
- `software_version = 3.0.1`
- `counting_semantics = cpm_v3_read_pair_and_exact_raw_umi_per_clonotype_v1`
- `settings.umi.exact_umi_family_views`にrow filter、counting unit、割合分母を記録
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

`exact_umi_family_counts`は注釈後にvalid exact UMI familyだけを5列へ投影した表であり、
Ver2.0の`umiStrict`や全長read Hamming距離collapseを復活させるものではない。

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
- `productive=true`はIgBLASTによる配列上の予測であり、実タンパク質発現の証明ではない。
- R1/R2は部分配列であるため、片側のproductiveが空欄でも元分子がnonproductiveとは
  限らない。両側trueを必須にする解析は、標準表とは別の追加QCとして扱う。

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
- exact UMI family表が5列で、`umi_family_count=0`行を含まない。
- `umi_family_percent`の分母が各exact UMI family表の`umi_family_count`合計である。
- productive exact UMI family表がproductive対象pairから独立再計算される。
- counts外BCR rowが`integrated.tsv`に残る。
- productive限定表が対象行から独立再計算される。
- XLSX percentが数値cellである。
- AIRR row部分欠落を公開しない。
- input/output/resource衝突を拒否する。
- publish途中失敗・Ctrl+C・rollback異常で偽completion manifestを残さない。
- 入力が解析中に置換された場合、正式出力を公開しない。
