# v2（0.2.0）の範囲

作成日: 2026-10-29
更新日: 2026-10-29

0.2.0 に向けた機能追加の範囲を記録する。今後の合意に応じて、このディレクトリへ詳細を追記する。

## 確定範囲

以下の3機能を v2 の実装対象とする。

| 機能 | 内容 |
| --- | --- |
| Graded tensor product | 共通機構を実装済み。[共通評価と移行](../../design/composite-expression.md#12-ordinary--graded-の共通評価と移行)を参照 |
| Fermion | Fermion の局所代数、反交換関係・符号規則への対応 |
| Spin-one | Spin-1 の局所代数の提供 |

Fermion / spin-one の詳細仕様・公開 API・crate 構成と実装順序は、以下の実装案で検討する。Graded tensor product の実装済み契約は設計文書に記録する。

## 詳細検討

- [Fermion の局所 primitive 実装案](fermion-primitives.md)：単一 mode の局所規則・parity・基底・係数要件と、spinful への接続範囲。
- [Spin-one の局所 primitive 実装案](spin-one-primitives.md)：三状態・全九次元演算子空間の表現、局所規則・完全基底・√2 の係数能力。

## 関連文書

- [アーキテクチャ概要](../../architecture.md#22-crate-内のモジュール構造)：実装済みのモジュール構成
- [今後の課題](../../future-work.md)：初期版から延期した課題の一覧
- [要件](../../requirements.md)：既存の数学的契約と検証条件
- [複合式の設計](../../design/composite-expression.md)：ordinary / graded の合成則
- [Hilbert 空間の設計](../../design/hilbert-space.md)：因子構造と tensor product
