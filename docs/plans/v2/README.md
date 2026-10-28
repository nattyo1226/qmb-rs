# v2（0.2.0）の範囲

作成日: 2026-10-29
更新日: 2026-10-29

0.2.0 に向けた機能追加の範囲を記録する。今後の合意に応じて、このディレクトリへ詳細を追記する。

## 確定範囲

### 先行する構造整理

以下の3機能の実装に先立ち、各 crate 内のモジュール構造を整理する。具体的な分割・配置と公開範囲は順次決定する。

現在のモジュール構成と責任範囲は[モジュール構造整理：現状調査](module-structure.md)に記録する。

構造整理の[決定事項と構成案](module-structure-plan.md)に、primitive の次元能力を `Dimension` へ改名して `Identity` とともに `ops` へ移す方針、恒等元の扱い、および core 全体の配置方針を記録する。有限線型結合・共通列は `expr::parts`、評価処理はトップレベルの `evaluate`、基底処理は `basis` にまとめる。

geometry は `boundary`・`geometry`・`relation`・`connection` に責任を分ける。幾何の種類は `geometry`、関係の種類は `relation` 配下で拡張する。

export は公開データ型・出力引数・変換手順・backend を分離する。spin-half は `utility` を `models` に改名し、既存の primitive・基底の分離を維持する。

normalize は `norm`・`expectation`・`error` に分け、既存の公開 API を crate root から提供する。

### 追加機能

以下の3機能を v2 の実装対象とする。

| 機能 | 内容 |
| --- | --- |
| Graded tensor product | Grading と parity に基づく tensor product・合成則への対応 |
| Fermion | Fermion の局所代数、反交換関係・符号規則への対応 |
| Spin-one | Spin-1 の局所代数の提供 |

各機能の詳細仕様、公開 API、crate 構成、および3機能間の実装順序は今後決定する。

## 関連文書

- [今後の課題](../../future-work.md)：初期版から延期した課題の一覧
- [要件](../../requirements.md)：既存の数学的契約と検証条件
- [複合式の設計](../../design/composite-expression.md)：ordinary / graded の合成則
- [Hilbert 空間の設計](../../design/hilbert-space.md)：因子構造と tensor product
