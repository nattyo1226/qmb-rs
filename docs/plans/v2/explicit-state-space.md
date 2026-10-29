# 所属空間を明示する状態のテンソル構築

作成日: 2026-10-30
更新日: 2026-10-30
状態: 設計方針を採用済み、未実装

## 1. 採用する方針

Hilbert 空間を先に構築し、状態のテンソル構築には所属先の空間を明示的に渡す。テンソル積の種類の情報源を `HilbertSpace` に集約し、状態構築側の `tensor_g`・`tensor_product_g`・種類指定版を廃止する。状態式の `Tensor` 節点からも `TensorKind` を除く。

`State` は所属空間を必ず持つ。演算はその空間、または当該演算で選択した対象空間から ordinary / graded を取得する。遅延作用は構築時に選択した空間を保持し、外側への合成で作用範囲を変更しない。

この文書は [関数名整理の N2・L2・L4](api-function-naming.md) に関連する設計変更を記録する。状態・QuantumSystem の tensor 構築については、同文書の旧来の `_g` / `_with` の改名案より本方針を優先する。Hilbert 空間自体の種類指定入口は引き続き必要であり、空間側の `tensor_g` を廃止する決定ではない。

ここまでが採用した設計方針である。以下の `State::tensor_in`・`F::product_state` は具体化のための API 案とし、正式な関数名・シグネチャは実装前に確定する。実装時期は未定とする。

## 2. 現在の実装と変更理由

現在の `State::tensor_with_kind` は、入力状態から新しい Hilbert 空間を構築し、式木の Tensor 節点にも同じ種類を保存する。[構築処理](../../../crates/qmb-core/src/expr/composite.rs#L340)と[節点の定義](../../../crates/qmb-core/src/expr/composite.rs#L104)を参照。このため、状態を合成する入口にも ordinary / graded の選択が現れる。

一方、演算子の作用評価は、[選択した空間の effective_kind](../../../crates/qmb-core/src/evaluate/composite.rs#L187) から種類を取得している。状態の Tensor 節点の種類は再構築で引き継がれるが、現在の作用の符号判定の情報源ではない。空間の選択を状態構築より前に移すことで、状態構築 API と式木に種類を重複して持たせる必要がなくなる。

演算子は既に空間を新設せず、束縛先や作用先の空間を使って解釈される。状態側も既存の所属空間を使って構築・解釈する流れに揃える。ただし未所属の State や、構築済みの状態を任意の空間へ読み替える機構は追加しない。

## 3. 情報を保持する場所

| 対象 | 保持するもの・役割 |
| --- | --- |
| `HilbertSpace` | 因子構造・因子順序・ordinary / graded の種類。有効な種類の計算と互換性の契約は維持する |
| `State` | 所属する `HilbertSpace` の共有ハンドルと状態式。空間を必須とし、状態全体の even / odd フラグは追加しない |
| `StateNode::Tensor` | 順序付きの子の列。種類のフィールドを除き、構築・再構築で独立した種類を指定しない |
| `StateAction` | 作用の構築時に選択した Hilbert 空間の共有ハンドル。局所作用の範囲と因子順序を固定する |
| `Operator` | 現行の任意の束縛を維持する。評価時の空間選択と明示空間の優先順位も維持する |

「空間に情報を集約する」とは、種類を表す独立したフラグを式木へ複製しないことを指す。状態や作用節点が必要な Hilbert 空間を共有することは認める。primitive の parity と空間の合成規則は別の情報であり、混合 parity の状態の和も引き続き表現する。

## 4. API と使用イメージ

以下は変更後の使用イメージであり、現時点でコンパイルできるコードではない。空間構築と局所式の構築には現行 API を使用し、状態合成の入口だけを提案している。

```rust
use num_complex::Complex;
use qmb_core::{HilbertSpace, QuantumSystem, State};
use qmb_fermion::{Fermion, FermionOperator as O, FermionState as S};

type F = Fermion<Complex<f64>>;

let a = HilbertSpace::factor("a");
let b = HilbertSpace::factor("b");

// 種類・順序・構造を空間構築で確定する。
let h = HilbertSpace::tensor_g([a.clone(), b.clone()])?;

// 提案する基礎 API：配置済みの状態を、指定した空間上で合成する。
let left = F::state(S::Occupied).at(&a)?;
let right = F::state(S::Occupied).at(&b)?;
let psi = State::tensor_in(&h, [left, right])?;

// 提案する補助 API：葉の順序に沿って primitive を配置して合成する。
let same = F::product_state(&h, [S::Occupied, S::Occupied])?;

// 現行と同じ作用 API。graded の規則は所属空間 h から取得する。
let remove_b = F::operator(O::Annihilate).at(&b)?;
let result = psi.act(&remove_b)?;
// result は -|Occupied, Empty> を表す。
```

ordinary 空間では空間構築を `HilbertSpace::tensor(...)` に変更するだけで、状態合成の呼び出し方は同じになる。

基礎 API は配置済みの式を受け取るため、局所的な重ね合わせや、複数因子にまたがるもつれた部分系も合成できる。

```rust
let left = F::state(S::Occupied).at(&a)?;
let right = (F::state(S::Empty) + F::state(S::Occupied)).at(&b)?;
let psi = State::tensor_in(&h, [left, right])?;
```

補助 API の `product_state` は primitive の列専用とし、`h.leaves()` の順に一つずつ配置する案を推奨する。基礎 API と異なる入力を同じ名前の overload で受け取らない。QuantumSystem が primitive と係数型を選ぶ現在の役割を維持し、HilbertSpace に具体的代数の知識を持たせない。

## 5. 構築時の契約

- 指定された空間を結果の所属空間として共有する。入力から別の空間を作り直さず、指定空間の grouping を保持する。
- 基礎 API では、入力状態の葉の Index 列を連結したものが、指定空間の葉の Index 列と順序も含めて完全に一致することを検証する。重複・不足・余分な因子を拒否する。既存の空間互換性と同様に、入力と対象の grouping 自体の一致は要求しない。
- 複数因子を持つ入力部分系の有効な種類は、指定空間と互換でなければならない。ordinary の部分系を graded として再解釈しない。単因子・単位空間の種類の扱いは現行の互換性に従う。
- 補助 API では primitive の数と葉の数を完全一致させる。暗黙の補完・切り捨て・並べ替えを行わない。入力順序が葉の順序への割り当てを意味し、交換操作や交換符号を自動で追加しない。
- 空入力は指定空間が零因子の単位空間である場合に受理し、単位状態を返す。非零因子の空間では入力不足として拒否する。単位状態にも指定空間のハンドルを保持する。零状態の構築は既存の `State::zero(&space)` を使う。
- 単項入力も空間との互換性を検証し、結果には指定空間を保持する。既存の部分木の共有を利用し、単項を理由に指定空間を無視しない。
- 構築に primitive の評価や係数演算を要求しない。Index の生成・Clone・ソートを追加せず、検証は Index の同一性と順序に基づく。

fallible な列集約も設ける場合は、例えば `State::try_tensor_in(&space, states)` とし、最初の入力 `ExpressionError` で消費を止め、全入力成功後に空間を検証する既存契約を維持する。正確な入口の一覧、エラー variant と診断情報は実装前に確定する。

## 6. 遅延作用の対象空間を保持する

既存の [StateAction の構築](../../../crates/qmb-core/src/expr/composite.rs#L380) は、選択した対象空間を作用節点に保存する。この仕組みを維持する。

```rust
// 提案 API による比較。s_a / s_b はそれぞれ a / b 上の状態。
let local_first = State::tensor_in(&h, [s_a.clone(), s_b.act(&remove_b)?])?;
let global_first = State::tensor_in(&h, [s_a, s_b])?.act(&remove_b)?;
```

`s_a` と `remove_b` が odd なら、全空間で作用を構築した後者だけに、先行する a の parity に由来する負号が現れる。前者の作用は b の空間内で既に構築されているため、外側の合成で範囲を h へ拡大しない。各作用節点の空間から種類を取得し、常に最外側の空間だけで評価する方式にはしない。

この区別は評価前後の合成、基底展開、部分木の共有・再構築でも維持する。

## 7. 移行範囲と検証

State の `tensor` / `tensor_product` と各 `try_`・`_g`・`_with_kind` 版、および QuantumSystem の `state_tensor*` の入口を、明示空間を受け取る構築へ移行する。単なる改名ではなく入力契約の変更として扱い、種類指定を状態から空間の構築へ移す。旧入口から空間を自動生成して恒久的に併設する方針は採らない。

HilbertSpace のテンソル構築・種類指定、Operator の束縛と作用時の空間選択、局所式の構築・配置は維持する。空間側の関数名の短縮は命名文書で引き続き検討する。基底列挙・展開で State を再構築する箇所も既存空間を引き継ぎ、Tensor 節点に種類を復元しない。

実装時は次を検証する。

- ordinary / graded の同一の構築入口と、graded 作用の交差符号。
- 指定空間の共有・grouping・Index 順序、重複・不足・過剰・順序違いと種類の不整合。
- 単位空間・空入力・単項・零状態、および単因子・単位 wrapper の互換性。
- もつれた部分系の合成と、混合 parity の状態の和。
- 局所作用と全空間での作用の差、評価前後の合成と展開・再構築での意味の保存。
- fallible 入力の消費停止とエラー優先順位、構築に不要な数値能力・primitive / Index の Clone を要求しないこと。
- 既存の内積・期待値・基底処理・normalize・export・モデル構築・examples の移行。

この文書の追加時点では実装・実行テストは行っていない。既存の設計文書は現行実装の記録として維持し、実装の移行時に更新する。

## 関連文書

- [関数名・公開 API 名の整理候補](api-function-naming.md)
- [複合式の設計](../../design/composite-expression.md)
- [Hilbert 空間の設計](../../design/hilbert-space.md)
