# Fermion の局所 primitive 実装案

作成日: 2026-10-29
更新日: 2026-10-29
状態: 検討案。記載する crate・型・追加 API は未実装。

## 1. 対象と既存実装との関係

単一 fermion mode の二次元占有空間、その上の全演算子代数、および局所基底を実装する案を記録する。局所の積・作用・内積・随伴・通常 trace・基底展開は、現在の `qmb-core` の能力だけで実装・検証できる。因子間の符号、空間の種類、Jordan–Wigner string は [tensor product の計画](tensor-products.md)の担当とし、本書の局所規則に埋め込まない。

spinful fermion は [要件 2.4 節](../../requirements.md)に従い二つの基本 mode の局所複合系とする。単一 mode の primitive を四状態の別 primitive へ置き換えることは、本案の前提にしない。単一 mode の完成と、graded 合成を使った spinful / 多 mode 系の完成を分けて判定する。

参照した既存契約と実装は以下のとおり。

| 参照 | 本案に適用する内容 |
| --- | --- |
| [局所構造の設計](../../design/local-structure.md) | S・O が意味と局所規則を所有。Identity は primitive の値、Dimension は span の次元 |
| [primitive / 式の設計](../../design/primitive-expression.md)・[ops.rs](../../../crates/qmb-core/src/ops.rs) | `Product<C>`・`Act<C>`・`Adjoint<C>` の関連型と `LinearCombination`、評価能力を構築へ要求しない |
| [基底の設計](../../design/basis.md)・[basis.rs](../../../crates/qmb-core/src/basis.rs) | 完全な基底、expand / index / elements の整合性、通常の式への展開 |
| [spin-half の operator.rs](../../../crates/qmb-spin-half/src/operator.rs)・[state.rs](../../../crates/qmb-spin-half/src/state.rs)・[basis.rs](../../../crates/qmb-spin-half/src/basis.rs) | 数値精度を持たない enum、個別能力の実装、基底 marker の構成 |
| [spin-half の lib.rs](../../../crates/qmb-spin-half/src/lib.rs)・[system.rs](../../../crates/qmb-core/src/system.rs) | `QuantumSystem` による係数型の選択、共通の局所式を返す入口 |
| [要件](../../requirements.md) 2.4、5.5 節 | 単一 mode の CAR、homogeneous な primitive、P の規約、spinful の位相 |

調査時点の workspace に fermion crate と core の parity trait は存在しない。spin-half は通常のスピン行列を提供するが、別因子の反交換関係を持つ fermion と同一の型にはしない。

## 2. 提案する crate と公開型

`qmb-fermion` を追加し、`no_std + alloc`、`unsafe` 禁止、crate root からの主要型・基底の再公開という既存の具体的代数の構成に合わせる。

```rust
enum FermionState { Empty, Occupied }

enum FermionOperator {
    I,
    Annihilate,
    Create,
    Number,
    Parity,
}
```

全 variant は有効な局所値なので、整数ラベルや占有数を受ける検証付き constructor は不要である。`Clone, Copy, Debug, PartialEq, Eq, Hash` を実装する。Zero は primitive に追加せず、`LinearCombination::zero()` と共通式の zero を使う。

`Fermion<C>` を `PhantomData<fn() -> C>` の型選択 marker とし、`QuantumSystem` の State / Operator / Coeff を上記型と C に関連付ける。marker や primitive の保持・`state` / `operator` に評価・数値能力の bound を課さない。独自の式木・evaluator・演算子構文は作らない。

想定ファイルは `src/{lib,state,operator,basis}.rs`、必要なら小さな定数生成をまとめる `coefficient.rs`、`tests/{contracts,minimal_scalar,local_expressions}.rs` とする。本番依存は既存 workspace の `qmb-core` と `num-traits`、`num-complex` は両精度の検証用 dev-dependency とする。Geometry・Index 発行・数値 backend に依存しない。

## 3. 規約・次元・parity

状態の順序を `(|0>, |1>) = (Empty, Occupied)` とし、内積を `⟨m|n⟩ = δ_mn` とする。`a = Annihilate`、`a† = Create`、`n = Number`、`P = Parity` と書く。

```math
a=\begin{pmatrix}0&1\\0&0\end{pmatrix},\quad
a^\dagger=\begin{pmatrix}0&0\\1&0\end{pmatrix},\quad
n=\begin{pmatrix}0&0\\0&1\end{pmatrix},\quad
P=\begin{pmatrix}1&0\\0&-1\end{pmatrix}=I-2n.
```

これらの行列で規約を固定し、局所の match 規則とは独立したテストの参照値にする。製品の局所規則を密行列の乗算として実装する指定ではない。

| 対象 | 値 |
| --- | --- |
| `FermionState::dim()` | 2 |
| `FermionOperator::dim()` | 4。5個の primitive は線型従属 |
| `Identity::identity()` | `FermionOperator::I`。係数型に依存しない |
| 状態の parity | Empty は even、Occupied は odd |
| 演算子の parity | I・Number・Parity は even、Annihilate・Create は odd |

P を primitive に持つ理由は、局所 parity の生成・作用を係数の除算なしで利用するためである。`P = I - 2n` を保持しても全演算子空間は四次元のままで、Dimension や基底要素数を5にしない。

## 4. 局所能力の実装

### 4.1. Product の全規則

表の行を左 operand、列を右 operand とし、右端が先に作用する。結果は同じ O の `LinearCombination<O,C>` とする。

| 左 × 右 | I | a | a† | n | P |
| --- | --- | --- | --- | --- | --- |
| I | I | a | a† | n | P |
| a | a | 0 | I − n | a | −a |
| a† | a† | n | 0 | 0 | a† |
| n | n | 0 | a† | n | −n |
| P | P | a | −a† | −n | I |

`Product<C, Term=FermionOperator>` は25組の直接 match とし、identity の積では元の自然な primitive を返す。例えば `a a†` は `I − n` の二項、nilpotent な積は空の線型結合とする。通常の Product に graded 交換則を追加しない。ここでの局所反交換関係そのものが具体的代数の規則である。

表から `a²=(a†)²=0`、`aa†+a†a=I`、`n=a†a`、`n²=n`、`P²=I`、`Pa=-aP`、`Pa†=-a†P` が得られる。P と a の掛け順を取り違えない。

### 4.2. Act・内積

`Act<C>` は `Operator=FermionOperator`、`Term=FermionState` とする。

| 演算子 | Empty への作用 | Occupied への作用 |
| --- | --- | --- |
| I | Empty | Occupied |
| a | 0 | Empty |
| a† | Occupied | 0 |
| n | 0 | Occupied |
| P | Empty | −Occupied |

状態ごとの parity が明確なので、P の作用と even / odd の定義が一致する。混合 parity の状態は共通の局所状態式の和で表す。

### 4.3. Adjoint・Trace

`Adjoint<C, Term=FermionOperator>` は a と a† を交換し、I・n・P を保存する。primitive の係数は1とし、式の係数の共役は core が担当する。

| 演算子 | I | a | a† | n | P |
| --- | --- | --- | --- | --- | --- |
| 通常の非正規化 trace | 2 | 0 | 0 | 1 | 0 |

`Tr(I)=2` とし、正規化 trace や supertrace に変更しない。局所 trace の実装に状態基底・行列を要求しない。

## 5. 基底と展開

状態基底 `OccupationBasis` は `[Empty, Occupied]`、既定の演算子基底 `OccupationOperatorBasis` は `[I, Annihilate, Create, Number]` の順とする。型を分けて `Basis` の対象を明確にする。

演算子基底の四要素は独立で全空間を張る。parity は homogeneous だが、I と n の Hilbert–Schmidt 内積は0ではなく、正規直交基底とは扱わない。

`OccupationOperatorBasis::expand` は基底 primitive をそのまま返し、P だけを `I − 2n` に展開する。`index(P)` は `None`、`elements()` は基底の四要素だけを一度ずつ返す。線型従属な全 primitive を列挙しない。

P を使う表現のために、第二の完全基底 `ParityBasis = [I, Parity, Annihilate, Create]` も提供する案とする。この基底では `n = (I − P)/2`、`index(Number)=None` である。両基底間の直接展開・再展開と、局所式の作用・trace の一致を検証する。基底の順序を演算子の行列成分と同一視しない。

## 6. 係数型と能力の最小化

既存 `LinearCombination::primitive` が `C: One` を要求することと、`One`・`Zero` 自体が持つ算術上の制約を前提とする。次は具体的 impl の想定 bound であり、式の保持へは伝播させない。

| 能力 | 想定する C の能力 |
| --- | --- |
| Identity、Dimension、parity 問い合わせ、QuantumSystem | 数値制約なし |
| InnerProduct、Trace | `Zero + One`。2は `one + one` で生成 |
| Adjoint、状態基底 | `One` |
| Product、Act | `One + Neg<Output=C>`。零は空の線型結合で返す |
| OccupationOperatorBasis | `One + Add<Output=C> + Neg<Output=C>`。2n を生成 |
| ParityBasis | 上記に `Div<Output=C>`。1/2 を表せることも数学的契約 |

局所規則には Float・平方根・ImaginaryUnit・C の Clone / Copy を要求しない。したがって既定基底と局所演算は符号付き整数や独自の厳密係数でも表せる。一方、ParityBasis に整数の切り捨て除算を使ってよいという意味ではない。共通式の評価・集約・export には、それぞれ core が要求する Clone・共役・算術能力が別途必要である。

精度は C が決め、`Complex<f32>` と `Complex<f64>` を暗黙に混ぜない。定数のために f64 経由のキャストやパースを追加しない。

## 7. Tensor 実装から独立した接続方針

局所 primitive には係数に依存しない `is_odd() -> bool` を用意できる。core の parity trait が確定した段階で、この値に対応する trait を具体的代数側で実装する。係数型に依存しない P の取得能力には `FermionOperator::Parity` 自身を返す。仮称の core trait を本書の局所実装の前提にせず、同じ数学的契約を重複定義しない。

多 mode の順序・交換符号・先行区間の P は core の合成機構に委ねる。`FermionOperator` に Index・Site・mode 番号・parity string・global ordering を保持しない。既存 ordinary tensor へ複数の a を配置するだけでは fermion の異なる mode の CAR は得られないため、その構成を fermion 多体系の利用例として案内しない。

spinful の内部順序は `(up,down)`、二重占有の位相は `|up down>=f_up† f_down†|0>` とする。合成側との接続試験では、`f_up|up down>=|down>`、`f_down|up down>=-|up>`、総 parity `diag(1,-1,-1,1)` を確認する。表示順 `[empty, up, down, double]` は、二 mode の mixed-radix 順 `[empty, down, up, double]` と異なるので明示的に対応付ける。この表示順の違いを単一 mode primitive の符号へ埋め込まない。

## 8. 実装順序と受入条件

1. 新 crate、enum、Dimension・Identity、型選択 marker と再公開を追加する。
2. Product・Act・InnerProduct・Adjoint・Trace と二種類の演算子基底を実装する。
3. core の局所式を通した評価・展開・再展開・座標取得、既存 export / normalize の単一 mode 利用を検証する。
4. tensor 側の API が確定してから parity trait を接続し、多 mode / spinful の統合検証を別途行う。

局所実装の完成条件は以下とする。これらは実装後の検証計画であり、現時点の Rust テスト結果ではない。

- 25組の Product、10組の Act、5種の Adjoint / Trace、4組の内積を、3節の独立な2×2行列・二成分ベクトルと照合する。
- primitive の全三重積125組の結合則、全 Product / 状態の50組の作用整合性、単位元、随伴の反準同型性・二回の随伴、trace の巡回性を検証する。構造の一致ではなく値を比較する。
- 局所 CAR、n の射影性、P の自己随伴性・平方・作用、各非零返却項の parity を確認する。P を省いた表と P を含めた表を別物として扱わない。
- 状態基底2要素・演算子基底4要素の列挙 / 番号 / 完全性と、全 primitive の展開・再展開を確認する。`P=I−2n` を経由した評価も同じ値になる。
- `Complex<f32>` / `Complex<f64>` の局所式の作用・期待値・通常 trace・正規化・座標・2×2行列 export を照合する。Vec と任意 feature の ndarray / sprs を同じ基底順で比較する。
- Float・ImaginaryUnit・除算のない係数で既定の局所規則を利用できること、Clone / Copy のない係数で primitive の規則を利用できること、評価能力のない型でも式を構築できることを compile / capability テストで確認する。
- spin-half と独立した型として動作し、core に fermion 固有分岐を追加しない。多 mode の CAR・spinful の接続は局所完成条件と別に記録する。

実装時の検証入口は `cargo test -p qmb-fermion`、workspace 登録後の `cargo test --workspace --all-features`、`cargo check -p qmb-fermion --no-default-features`、`cargo fmt --all -- --check` とする。ホスト上の no-default-features check と専用 no_std target の検証は区別し、未実施事項を明記する。

## 9. 設計判断と残る詳細

推奨案は、a・a†・n・P・I を自然な primitive として保持する。四つの行列単位だけにすると規則は単純になるが、現行 `Identity` に必要な I は別途 primitive の値として追加する必要がある。物理的な名前と P の直接保持を優先し、五つの線型従属な primitive を採る。

型・variant・基底の最終命名、第二基底の提供時期、core の parity trait の具体名は実装前に固定する。局所代数の規約・表・次元は tensor API の決定を待たず確定・検証できる。Hubbard / hopping モデル、spinful の構築 utility、粒子数 sector は本書の primitive 実装と別に計画する。

## 10. 本ドキュメント作成時の確認

現行の設計文書・core 契約・spin-half 実装を照合した。補助的な Python 標準ライブラリの行列計算では、本書の表を読み取って25組の積と10組の作用を参照行列と照合し、局所 CAR、P の平方・交換関係、trace の値も確認した。相対ファイルリンク・コードフェンス・差分の空白を確認している。

変更は計画文書のみであり、Rust の fermion 実装・crate 追加・ビルド・テストは行っていない。8節の受入条件は実装後に検証する。
