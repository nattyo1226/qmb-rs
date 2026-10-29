# v2（0.2.0）の範囲

作成日: 2026-10-29
更新日: 2026-10-30

0.2.0 に向けた機能追加の範囲を記録する。今後の合意に応じて、このディレクトリへ詳細を追記する。

## 確定範囲

以下の3機能を v2 の実装対象とする。

| 機能 | 内容 |
| --- | --- |
| Graded tensor product | 共通機構を実装済み。[共通評価と移行](../../design/composite-expression.md#12-ordinary--graded-の共通評価と移行)を参照 |
| Fermion | `qmb-fermion` に実装済み。[局所代数・基底・検証](../../design/local-structure.md#9-fermion-の局所代数)と[標準／symmetric Hubbard](../../design/geometry.md#6-fermion-の-hubbard-モデル)を参照 |
| Spin-one | `qmb-spin-one` に実装済み。[局所代数・基底・検証](../../design/local-structure.md#10-spin-one-の局所代数)と[XYZ・Heisenberg・BLBQ](../../design/geometry.md#7-spinone-のモデル構築)を参照 |

3機能とも実装済みであり、仕様・検証記録は上表の設計文書に集約する。この文書は v2 の対象範囲と参照先を案内する。

## 関連文書

- [アーキテクチャ概要](../../architecture.md#22-crate-内のモジュール構造)：実装済みのモジュール構成
- [今後の課題](../../future-work.md)：初期版から延期した課題の一覧
- [要件](../../requirements.md)：既存の数学的契約と検証条件
- [複合式の設計](../../design/composite-expression.md)：ordinary / graded の合成則
- [Hilbert 空間の設計](../../design/hilbert-space.md)：因子構造と tensor product
