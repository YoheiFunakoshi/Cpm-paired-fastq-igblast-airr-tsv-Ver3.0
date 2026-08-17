# CPM Paired Fastq IgBLAST AIRR tsv Ver3.0

CPM社のpaired FASTQを、非マージのR1/R2としてIgBLAST AIRR TSVへ変換し、
RG版と同じread-pair集計と、BCR clonotype内のUMI family集計を1回の解析から
作成するWindows向けツールです。

実装version: **3.0.0**

## Ver3.0の基本方針

Ver3.0では、UMIを理由にIgBLAST前のread pairをcollapseしません。

```text
R1/R2 FASTQ
  -> R2からUMIを抽出して各pairへ記録
  -> 全保持pairをIgBLAST（1回）
  -> RG版と同じ規則でR1/R2を統合
  -> BCR clonotypeを決定
  -> そのclonotype内で完全一致raw UMIを数える
```

このため、同じ解析から次の2単位を明確に分けて確認できます。

- **read-pair support**: RG版と同じ形式。PCR重複を含む保持pair数。
- **UMI-family support**: 同一clonotype内のdistinctな完全一致raw UMI数。

Ver2.0の`umiSeq5`で使用していたR1/R2全長Hamming距離によるcollapseと、
サンプル全体での近傍UMI統合は、Ver3.0の標準処理では行いません。
Ver2.0とVer3.0のUMI値は数値互換ではありません。

詳しい定義は[仕様書](SPECIFICATION.md)と
[解析方式の説明](docs/CPM_Ver3_method.md)を参照してください。

## BCR clonotypeの定義

RG版と同じ3項目をkeyにします。

```text
unique_v_gene_set
+ unique_j_gene_set
+ final_junction_aa
```

- V/Jのアリル番号を外し、候補セットを正規化します。
- 候補セットは完全一致で判定します。
- D callは追跡用に保存しますが、keyには使いません。
- `final_junction_aa`は、stopなし、C開始、W/F終了、5～40 aaをcanonicalとします。
- R1/R2でjunction AAが異なる場合は、片側だけcanonicalならcanonical側、
  両側canonicalならR1を採用します。

これはBCR全塩基配列の完全一致ではなく、clonotypeの定義です。同義塩基置換や
key外のSHM差は同じclonotypeに入ることがあります。

## Productiveの定義とR1/R2の扱い

IgBLAST/AIRRの`productive`は、そのqueryのV(D)J配列が、in-frame junction、
内部stop codonなし、V領域の内部frame shiftなしなどの条件から、タンパク質を
コードできると**予測**されたことを表します。実際にタンパク質として発現した
ことを直接証明する値ではありません。

定義の背景は[AIRR Rearrangement Schema](https://docs.airr-community.org/en/stable/datarep/rearrangements.html)と
[NCBI IgBLAST productive update](https://blast.ncbi.nlm.nih.gov/doc/blast-news/2021-BLAST-News.html)を参照してください。

このソフトで集計keyに使う`final_junction_aa`は、保存CとW/Fを含むAIRRの
JUNCTIONアミノ酸配列であり、BCR全タンパク質配列ではありません。また、
`complete_vdj=true`はproductive表の追加条件にしていません。

Ver3.0の`final_productive`は、R1とR2のAND判定でもOR判定でもありません。

1. `final_junction_aa`を決めたreadを優先する。
2. そのreadの`productive`値を`final_productive`へ採用する。
3. 優先readに値がなければ、反対側の値へfallbackする。
4. R1/R2のjunction AAが同じでproductive値だけが異なる場合はR1を優先する。

したがって、R1=`true`、R2=`false`または空欄でも、R1が採用側なら
`final_productive=true`になり得ます。空欄は`false`ではなく、IgBLASTがそのreadに
値を出力しなかったことを表します。判定に必要な範囲の不足などが原因として
考えられますが、空欄だけから原因は断定しません。AIRR TSVの実値は通常`T`、`F`、
空欄であり、この説明の`true`/`false`はそれぞれの論理値を表します。

CPM R2の301塩基には先頭38塩基のanchor/UMI/区切り配列が含まれるため、BCR insertは
最大約263塩基です。R2はV callを得てもJ/junctionまで届かず、productiveが空欄に
なることがあります。このため、R1/R2両方のproductiveを主解析の必須条件には
しません。両方trueは必要に応じて`integrated.tsv`から確認する追加QCです。

## UMIの数え方

- 同一clonotype内で同じ有効raw UMIが複数回あっても1 familyです。
- 同じUMI文字列が別clonotypeにあれば、それぞれで1 familyです。
- サンプル全体で同じUMIを1つにまとめません。
- 1塩基違いのUMIを自動統合しません。
- familyとして有効なのは`A/C/G/T`だけからなる12-merです。
- UMIを抽出できないpairや、`N`など曖昧塩基を含むpairも捨てません。
- UMI集計によって`integrated.tsv`からread pairを削除しません。

UMIがないpairは元分子数へ補正できないため、列を分けます。

```text
inclusive_support_count
= umi_family_count + umi_missing_read_pair_count
```

`inclusive_support_count`は、UMI familyと補正不能read pairを合わせた包含的な
支持量です。厳密な分子数とは呼びません。

## 主な出力

例として出力基準名を`sample.airr.tsv`とした場合:

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
sample.queries.fasta             # 保存を指定した場合のみ
sample.run.json                  # 完了manifest
```

### 4つのExcelの内容

4つは別々のIgBLAST解析ではありません。同じR1/R2統合結果を、集計する支持量と
productive限定の有無で分けた2×2の集計表です。

| Excel | 収載する主な支持量 | 対象 | 主な用途 |
|---|---|---|---|
| `integrated_counts.xlsx` | read pair | V/J/canonical junction AAを持つ全採用pair | RG互換のread-pair比較、全体QC |
| `final_productive_counts.xlsx` | read pair | 上記のうち`final_productive=true` | productive限定のRG比較 |
| `umi_counts.xlsx` | exact raw UMI familyとUMI missing pair | 全採用pair | CPMのUMI支持とPCR重複の確認 |
| `final_productive_umi_counts.xlsx` | exact raw UMI familyとUMI missing pair | `final_productive=true`のpair | productive pairから観測されたCPM UMI支持の評価 |

どのExcelでも1行は次の3項目で定義した1 BCR clonotypeです。

```text
unique_v_gene_set
+ unique_j_gene_set
+ final_junction_aa (canonical)
```

1行の単位はclonotypeですが、`read_pair_count`、`umi_family_count`、
`umi_missing_read_pair_count`は互いに異なる支持単位です。

read-pair表は次の10列です。

| 列 | 意味 |
|---|---|
| `unique_v_gene_set` / `unique_j_gene_set` | allele suffixを除去し正規化したV/J候補集合 |
| `final_junction_aa (canonical)` | 採用されたcanonical junction AA |
| `read_pair_count` | そのclonotypeを支持する保持pair数。PCR重複を含み、分子数ではない |
| `match_count` | R1/R2のjunction AAが一致したpair数 |
| `conflict_count` | R1/R2のjunction AAが異なり、統合規則で片側を採用したpair数 |
| `r1_only_count` / `r2_only_count` | junction AAを片側だけから得たpair数 |
| `productive_true_count` | `final_productive=true`だったpair数 |
| `canonical_junction_aa_count` | canonical junction AAを持つpair数 |

UMI表は上記のstatus/count列に、次の列を加えた15列です。

| 列 | 意味 |
|---|---|
| `umi_family_count` | clonotype内のdistinctな完全一致raw 12-mer UMI数 |
| `read_pair_count` | counts条件を満たした保持pair数 |
| `umi_known_read_pair_count` | 有効UMIを持つ保持pair数。PCR重複を含む |
| `umi_missing_read_pair_count` | 有効UMIを得られず、UMI familyへの分子補正はできないが削除せず保持したpair数 |
| `inclusive_support_count` | `umi_family_count + umi_missing_read_pair_count` |
| `inclusive_support_percent` | 同じ表のinclusive support合計に占める割合 |

`final_productive_umi_counts`では、productive対象pairだけからUMI set、missing数、割合を
独立に再計算します。UMI表の既定順は`inclusive_support_count`降順です。
UMI familyだけで順位を見たい場合は`umi_family_count`で並べ替えます。

これはUMI family全体のconsensusをproductive判定した表ではありません。同じclonotype・
同じUMIに`T`と`F`のpairが混在しても、`T`のpairが1つ以上あれば、そのUMIは
productive対象内で1 familyとして数えられます。

目的別の確認先:

- RGと同じread-pair単位: `integrated_counts.xlsx`
- CPMのexact UMI family: `umi_counts.xlsx`
- productiveなread-pairだけ: `final_productive_counts.xlsx`
- productiveを主目的とするCPM解析の主解析候補: `final_productive_umi_counts.xlsx`

この場合も、通常表は除外されたpairやclonotypeを確認する監査・感度参照として残します。

productive限定による減少率はデータごとに異なります。必ず通常表とproductive表の
clonotype数、`read_pair_count`合計、`umi_family_count`合計、
`inclusive_support_count`合計を同じ単位同士で比較します。稀なsingleton clonotypeが
除かれると、支持量の減少率よりclonotype種類数の減少率が大きくなることがあります。

### 追跡用表

`integrated.tsv`には、R1/R2の元call、採用元、競合状態、UMI、counts除外理由を
保存します。有効なBCR patternを決められないpairもこの表から削除しません。

## RG版との比較

推奨する比較順:

1. RG `integrated_counts` vs CPM Ver3.0 `integrated_counts`
   - 同じclonotype定義・read-pair単位で比較
2. CPM Ver3.0 `integrated_counts` vs `umi_counts`
   - CPM内でPCR重複の影響を確認
3. RG vs CPM Ver3.0 `umi_counts`
   - CDR3/V/Jや順位は比較可能。ただしRGはread pair、CPMはUMI familyなので、
     絶対数を同一単位として扱わない

ライブラリ作製、PCR、シーケンス、trim条件の違いも残るため、同じ検体でも両社の
値が完全一致することを期待する解析ではありません。

## CPM R2 UMI抽出

既定ではR2先頭を次の構造として扱います。

```text
TATCAACGCAGAGTGGCCAT + NNNN + T + NNNN + T + NNNN + TCTT + insert
```

3個の4塩基blockを連結した12塩基をraw UMIとして記録します。anchorは既定で
2 mismatchまで許容します。UMIを認識できない場合は`NA`として保持します。
`N`を含む12-merもfamily数には入れず、UMI missing pairとして保持します。

## ポータブル作業フォルダ

標準配置:

```text
CPM Paired Fastq IgBLAST AIRR tsv Ver3.0
├─ Open CPM Paired Fastq IgBLAST AIRR tsv Ver3.0.cmd
├─ Launch CPM Paired Fastq IgBLAST AIRR tsv Ver3.0.ps1
├─ PORTABLE_USAGE.txt
├─ app
├─ python
├─ tools
│  └─ igblast-1.21.0
├─ refdata
│  └─ IgBlast_refdata_edit_imgt
└─ Results of CPM Paired Fastq IgBLAST AIRR tsv Ver3.0
```

通常は`Open CPM Paired Fastq IgBLAST AIRR tsv Ver3.0.cmd`をダブルクリックします。
移動可能な標準起動は`.cmd`です。`.lnk`は絶対パスを持つことがあるため標準配布には
使用しません。

## 開発環境での起動

Python 3.11以上を使用します。

Ver2.0とVer3.0は同じPython module名を使用するため、同じ仮想環境へ同時に
installしないでください。Ver3.0専用のvirtual environment、または同梱portable
Pythonを使用します。

```powershell
python -m pip install -e .
cpm-paired-fastq-igblast-airr-tsv-v3 gui
```

インストールせずに確認する場合:

```powershell
$env:PYTHONPATH = "src"
python -m airr_igblast_paired --help
python -m airr_igblast_paired gui
```

## 安全性

- R1/R2、IgBLAST executable/DB/aux、query、全出力のパス衝突を事前拒否
- scratchへ作成し、全処理成功後だけ正式出力を公開
- 既存結果は既定で上書きしない
- 上書き時も全出力をrollback可能なtransactionとして扱う
- 入力FASTQのpath、size、mtime、file identityを処理中に再確認
- 同じ出力先への同時実行をlockで拒否
- manifestを最後に公開し、manifest不在を未完了判定に使用
- IgBLAST resourceは内容fingerprint別の不変cacheへ配置

## テスト

```powershell
$env:PYTHONPATH = "src"
python -m unittest discover -s tests -v
```

## GitHubに含めないもの

- 研究用FASTQ
- 実解析AIRR TSV、Excel、manifest
- 研究資料や検体情報
- IgBLAST/IMGT参照DB本体
- portable Python/IgBLAST binary

GitHubにはソース、テスト、仕様書、launcher templateのみを置きます。

## 方法上の位置付け

Ver3.0は**annotation-first, clonotype-aware UMI counting**です。UMI consensus
配列を作るツールではありません。一般的なfeature内UMI countingや、alignment後に
UMIを利用してassembleするBCR workflowと整合する考え方ですが、consensus-firstの
BCR解析とは目的が異なります。

- [UMI-tools: gene/feature内でのUMI grouping/counting](https://umi-tools.readthedocs.io/en/latest/Single_cell_tutorial.html)
- [MiXCR: targeted BCR UMI workflow](https://mixcr.com/mixcr/guides/generic-umi-bcr/)
- [pRESTO: consensus-first BCR workflowの例](https://presto.readthedocs.io/en/stable/workflows/Stern2014_Workflow.html)
