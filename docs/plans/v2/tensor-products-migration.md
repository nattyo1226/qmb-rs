# 0.2.0: ordinary / graded tensor product の移行

2026-10-29。共通の tensor 機構を実装した。Fermion / spin-one の製品実装、異種 primitive、実質的な ordinary / graded 混在、permutation、regrouping、partial trace は今回の対象外である。

## 空間と状態の構築

`TensorKind::{Ordinary, Graded}` を構築時に選ぶ。既存の `tensor` / `tensor_product` は ordinary の短縮形であり、子や primitive から種類を推測しない。

| 対象 | 種類指定 | ordinary | graded |
| --- | --- | --- | --- |
| HilbertSpace / State | `tensor_with_kind(kind, children)` | `tensor(children)` | `tensor_g(children)` |
| HilbertSpace / State の二項積 | `tensor_product_with_kind(kind, &rhs)` | `tensor_product(&rhs)` | `tensor_product_g(&rhs)` |
| fallible な State | `try_tensor_with_kind(kind, states)` | `try_tensor(states)` | `try_tensor_g(states)` |
| QuantumSystem | `state_tensor_with_kind(kind, states)` | `state_tensor(states)` | `state_tensor_g(states)` |
| fallible な QuantumSystem | `try_state_tensor_with_kind(kind, states)` | `try_state_tensor(states)` | `try_state_tensor_g(states)` |
| 配置を伴う QuantumSystem | `state_tensor_at_with_kind(kind, assignments)` | `state_tensor_at(assignments)` | `state_tensor_g_at(assignments)` |

`kind()` は保持した Tensor 節点の種類、`effective_kind()` は空・単因子を同一視した後の種類を返す。Factor の `kind()` は `None`。空・単因子の有効な種類はないが、複数葉を持つ唯一の非空の子の種類は wrapper に継承する。木の grouping は保持する。

互換性は葉順と有効な種類、複数葉の部分系は相対順序と有効な種類を照合する。`SpaceError::MixedTensorKinds` を追加したため、`SpaceError` の網羅的 match は更新が必要である。`ExpressionError::Space`、`BasisError::Expression`、`ExportError` への既存の伝播経路を使う。Fallible 集約は入力の最初のエラーで消費を止め、全入力の成功後に空間を検証する。

## 独自代数の trait 移行

複合評価の共通 API は ordinary 空間でも必要な grading 能力を要求する。暗黙に even とする blanket impl はない。

```rust,ignore
impl HomogeneousParity for MyState {
    fn parity(&self) -> Parity { Parity::Even }
}
impl HomogeneousParity for MyOperator {
    fn parity(&self) -> Parity { Parity::Even }
}
impl ParityOperator for MyOperator {
    fn parity_operator() -> Self { Self::identity() }
}
```

上は自明な grading の代数の場合である。非自明な grading では各 primitive を homogeneous とし、混合 parity は式の和で表す。P は自己随伴な involution であり、状態の even / odd 成分にそれぞれ +1 / -1 で作用する。随伴は parity を保存し、局所積・作用は parity を modulo 2 で加える。法則の保証は代数実装の責任となる。Spin-half は上の自明な grading を明示実装する。

| 操作 | 従来の能力に追加するもの |
| --- | --- |
| 複合状態の evaluate / 内積 / 期待値 / 展開 / 座標、対応する export / normalize | O の `HomogeneousParity + ParityOperator` |
| 複合演算子の trace / Hilbert–Schmidt 内積 | O の `HomogeneousParity + ParityOperator` |
| 複合演算子の展開 / 座標 | O の `HomogeneousParity`、C の `Neg<Output=C>` |

構築・配置・束縛・遅延した作用や随伴には評価能力を追加していない。局所だけの API の bound も維持する。作用は局所 `Act` に P を渡すので、O の `Product` や交差符号のためだけの C の `One` / `Neg` は不要。S の parity を直接問い合わせる必要もない。

## 評価・基底・数値接続

`evaluate`、`trace`、`expand_in`、`coordinates_in` と各明示空間版、`Export`、normalize の公開メソッドを共通入口として使う。明示空間は旧束縛を上書きする。作用を構築した時点の選択済み空間は Act 節点に保存する。

graded で因子 j の odd primitive を作用させるときは、その Act の scope 内の先行因子へ P を作用させる。先に単因子上で構築した作用を外側の tensor に組み込んでも、scope は拡大しない。状態内積に追加の符号はなく、full trace は通常の非正規化 trace である。

演算子基底は因子順の配置積 `E_0(b_0)…E_n(b_n)`。その intrinsic 座標では odd 成分同士の交差だけに負号を加える。trace 用の通常表現には P string を挿入し、別の private な走査型で扱う。たとえば第二因子の消滅演算子の intrinsic 座標は `I ⊗g a`、状態基底上の行列は `P ⊗ a` になる。

密疎 export は基底状態への作用と状態座標から列を作る既存経路を共用し、backend に符号処理を追加しない。状態基底は右端が最速、metadata は選択空間の kind と因子順を保持する。

実行例は [tensor_products.rs](../../../examples/src/bin/tensor_products.rs)。非自明な grading は core のテスト専用二状態代数で検証し、製品用 fermion API は追加しない。
