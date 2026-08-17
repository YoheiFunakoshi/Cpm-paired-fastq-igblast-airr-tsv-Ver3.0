# Changelog

## Unreleased

- clonotypeごとにUMIを独立集計する具体例を追加。同じUMIが別clonotypeで干渉しないこと、
  UMI missing-only clonotypeが包含UMI表には残りstrict表からは除外されること、割合分母を明記。
- 6つのExcelそれぞれについて、大小・割合を見る列を追記。
- productiveかつUMI missingを捨てない解析では`final_productive_umi_counts.xlsx`を主解析、
  `inclusive_support_count` / `inclusive_support_percent`を主指標とし、productive strict表を
  感度解析、productive read-pair表をRG比較に使う方針を明記。

## 3.0.1 - 2026-08-17

- `productive`と`final_productive`の定義、R1/R2の採用・fallback規則、空欄と
  `false`の違い、実タンパク質発現との解釈境界を文書化。
- 既存4つのExcelをread-pair/UMIとproductive限定の2×2として整理し、全列の意味、
  目的別の使い分け、productive限定前後を同じ観測単位で比較する方法を追記。
- 既存の包含UMI表を変更せず、valid exact UMI familyだけをcount・割合分母にする
  `exact_umi_family_counts`と`final_productive_exact_umi_family_counts`をTSV/XLSXで追加。
- exact UMI family表は`umi_family_count > 0`のclonotypeだけを5列で出力し、
  `umi_family_percent`を各表のfamily合計から独立計算する。

## 3.0.0 - 2026-08-17

- IgBLAST前のUMI/全長配列collapseを廃止し、全保持read pairを1回だけ注釈する方式へ変更。
- RG版と同じR1/R2統合規則およびclonotype keyで、read-pair countsを作成。
- 注釈後、各clonotype内の完全一致raw UMIを独立に数えるUMI countsを追加。
- 同じUMI文字列が異なるclonotypeに現れた場合は、それぞれで1 familyとして保持。
- UMI missing pairを削除せず、family数とは別のread-pair支持として出力。
- UMI missingのpair単位率と、exact UMI縮約後のinclusive support内割合を分離して説明し、
  検証runの15.535%と46.123%が異なる分母によることを文書化。
- `integrated_counts` / `final_productive_counts`をRG互換のread-pair表とし、
  `umi_counts` / `final_productive_umi_counts`を別表として追加。
- Ver2.0の`umiSeq5`、近傍UMI補正、複数modeの同時IgBLASTを標準処理と公開UIから削除。
- 1解析・1出力セット・1completion manifestへ簡素化。
- Ver2.0で追加した入力保護、transactional publish、出力lock、参照DBの
  content-addressed staging、manifestを継承。
- パッケージ名、GUI、launcher、仕様書をVer3.0へ更新。

### 互換性

Ver3.0のUMI family数はVer2.0のcollapse後sequence-cluster数と定義が異なり、
数値互換ではありません。旧結果と列名だけで混在させず、manifestの
`software_version`と`counting_semantics`を確認してください。
