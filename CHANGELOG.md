# Changelog

## 3.0.0 - 2026-08-17

- IgBLAST前のUMI/全長配列collapseを廃止し、全保持read pairを1回だけ注釈する方式へ変更。
- RG版と同じR1/R2統合規則およびclonotype keyで、read-pair countsを作成。
- 注釈後、各clonotype内の完全一致raw UMIを独立に数えるUMI countsを追加。
- 同じUMI文字列が異なるclonotypeに現れた場合は、それぞれで1 familyとして保持。
- UMI missing pairを削除せず、family数とは別のread-pair支持として出力。
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
