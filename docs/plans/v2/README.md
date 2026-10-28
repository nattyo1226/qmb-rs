# v2（0.2.0）の範囲

作成日: 2026-10-29
更新日: 2026-10-29

0.2.0 に向けた機能追加の範囲を記録する。今後の合意に応じて、このディレクトリへ詳細を追記する。

## 確定範囲

以下の3機能を v2 の実装対象とする。

| 機能 | 内容 |
| --- | --- |
| Graded tensor product | Grading と parity に基づく tensor product・合成則への対応 |
| Fermion | Fermion の局所代数、反交換関係・符号規則への対応 |
| Spin-one | Spin-1 の局所代数の提供 |

各機能の詳細仕様、公開 API、crate 構成、および3機能間の実装順序は今後決定する。

## 詳細検討

- [Ordinary / graded tensor product の実装修正案](tensor-products.md)：現行実装との差分、情報の保持、評価・基底・export の変更、移行方針と受入条件。

## 関連文書

- [アーキテクチャ概要](../../architecture.md#22-crate-内のモジュール構造)：実装済みのモジュール構成
- [今後の課題](../../future-work.md)：初期版から延期した課題の一覧
- [要件](../../requirements.md)：既存の数学的契約と検証条件
- [複合式の設計](../../design/composite-expression.md)：ordinary / graded の合成則
- [Hilbert 空間の設計](../../design/hilbert-space.md)：因子構造と tensor product
