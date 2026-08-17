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

### RG互換read-pair表

最初に次を確認します。

```text
sample.integrated_counts.xlsx
```

列とclonotype keyはRG版の`integrated_counts`と同じです。RGとCPM Ver3.0を
read-pair単位で比較する場合に使います。両方ともPCR重複を含むため、分子数とは
解釈しません。

### UMI表

次に確認します。

```text
sample.umi_counts.xlsx
```

主な列:

- `umi_family_count`: clonotype内のdistinctな完全一致raw UMI数
- `read_pair_count`: counts条件を満たした保持pair数
- `umi_known_read_pair_count`: UMIを抽出できた保持pair数
- `umi_missing_read_pair_count`: UMIを抽出できなかった保持pair数
- `inclusive_support_count`: family数とUMI missing pair数の和
- `inclusive_support_percent`: 同じ表のinclusive support合計に対する割合

### Productive限定表

`final_productive_counts`と`final_productive_umi_counts`は、採用側の
`final_productive=true`だけを同じkeyで再集計した二次確認表です。R1/R2を
マージせず部分配列としてIgBLASTへ渡すため、基本解析ではproductiveを除外条件に
しません。

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
