# LocalStructure の役割・依存と採用方針

更新日: 2026-10-27

## 1. 位置づけと結論

[整理メモの3.2節](refactoring-notes.md#32-localstructure-の必要性)に対する調査と設計判断である。コミット `6fc01b6` の実装を調査し、2026-10-26 の議論で以下の方針を採用した。第2〜3節は移行前の実装、第4節以降は採用した移行方針を記録する。製品コード・テスト・利用例の移行を完了した。

**LocalStructure を廃止し、式に primitive 型 S・O を直接指定する。追加の定義ラベル・定義 ID は設けない。** primitive の意味・規則は S・O の型と各能力の実装で定める。保持と評価の分離を維持し、因子構造を扱う空間・Geometry／Chain から局所定義への依存を外す。

その後、同日に Chain が Site ごとの Index 列と境界条件を直接保持し、Geometry 構造体と空間保持を廃止する方針を採用した。本書の Geometry／Chain の空間保持と検証結果は LocalStructure 廃止時点の実装記録である。Chain の新方針への移行は2026-10-27に実施し、[Chain の設計](design/geometry.md)と[実装計画第9節](implementation-plan.md#9-chain-の-index-列保持への移行)に記録する。

### 1.1. 採用事項

| 責任 | 採用した配置 |
| --- | --- |
| primitive 型の関連付け | 演算子式に O、状態式に S・O を直接指定する。作用の型対応は評価時の Act 等の能力で要求する |
| 次元 | primitive 型の要素全体が張る空間の次元を提供する。現在の完全な状態・演算子空間では S が d、O が d²。恒等元と次元の配置は先に暫定採用した方針を引き継ぐ |
| 恒等元 | 演算子 primitive O に含め、係数型に依存せず O の値として生成する能力を提供する |
| 定義の区別 | S・O の型と各能力の実装が担う。同じ S・O に別の意味を与えるための追加ラベルは設けない |
| 因子の構造 | 空間が Index・順序・grouping を扱う。空間比較・bind のために L を保持しない |

数学的な契約は[局所構造の設計](design/local-structure.md#7-localstructure-廃止の採用方針2026-10-26)、実装の移行範囲は[実装計画](implementation-plan.md#8-localstructure-の廃止)にも記録する。State／Operator の一般の写像としての整理、空間の Grading を含む最小構成、Chain の内部表現の詳細は後続の議論とする。

## 2. 移行前の機能（調査記録）

以下は移行前の定義である。旧 module は削除し、現在の能力は [primitive.rs](../crates/qmb-core/src/primitive.rs) で提供する。

```rust
pub trait LocalStructure {
    type State;
    type Operator;
    const DIMENSION: core::num::NonZeroUsize;
}

pub trait Identity<C>: LocalStructure {
    fn identity() -> LinearCombination<Self::Operator, C>;
}
```

| 機能 | 現在の提供方法 | 保証の範囲 |
| --- | --- | --- |
| 状態・演算子 primitive 型の関連付け | `State`・`Operator` の関連型 | 式節点と、評価で要求する入出力の型を揃える |
| 局所状態空間の固定次元 | `DIMENSION: NonZeroUsize` | 正の定数を提供する。primitive の個数や規則との数学的整合性は自動検査しない |
| 固定定義・固定パラメータの区別 | 式・空間が型引数 `L` を共有する | 異なる `L` の接続を型で拒否する。同じ primitive 型・同じ次元でも区別できる |
| 恒等演算子の具体的な表現 | 別 trait `Identity<C>` | 演算子 primitive の有限線型結合を返す。LocalStructure 自体の必須能力ではない |

`L` の実行時の値は式・空間・Geometry に保持しない。空間は `PhantomData<fn() -> L>` を持ち、式は `L` から得た primitive と式節点を保持する。定義用の値を `Arc<L>` で共有・照合する構成ではない。`L: Clone + Send + Sync + 'static` も要求していない。借用を含む primitive と、Send／Sync でないメタデータ型を扱う既存テストがある。

### 2.1. 既に別の場所が提供している機能

| 機能 | 現在の担当 |
| --- | --- |
| 積・作用・内積・随伴・trace | primitive 側の [Product・Act・InnerProduct・Adjoint・Trace](../crates/qmb-core/src/ops.rs) |
| primitive の有効な値の検査 | 具体的な primitive の constructor と実装者の契約 |
| 基底の展開・順序番号・列挙 | [Basis<P, C>](../crates/qmb-core/src/basis.rs) と具体的な基底型 |
| 線型結合・式の保持・共有・遅延構築 | [LinearCombination](../crates/qmb-core/src/linear_combination.rs)・[局所式](../crates/qmb-core/src/expr/local.rs)・[複合式](../crates/qmb-core/src/expr/composite.rs) |
| 係数能力・数値精度 | 係数型 `C` と操作ごとの bound |
| 因子の Index・順序・grouping・重複検査 | [HilbertSpace](../crates/qmb-core/src/space.rs) |
| Index の発行 | [qmb-basic](../crates/qmb-basic/src/lib.rs) と新規発行時の `I: Default` |
| Site の対応と近接関係 | [Geometry・Chain](../crates/qmb-geometry/src/lib.rs) |
| norm・正規化・正規化期待値 | [qmb-normalize](../crates/qmb-normalize/src/lib.rs) |

数学的な整合性は、[局所構造の設計](design/local-structure.md)に契約として記述され、具体的代数のテストで検証する。LocalStructure の実装だけで結合則・正定値性・`Tr(I) = d` が証明されるわけではない。Grading・parity・実行時の定義パラメータの照合・数値 export も、この trait が現在提供している機能ではない。

## 3. 依存箇所の棚卸し

### 3.1. 関連型を実際に使う箇所

| 箇所 | 必要としている機能 | 分離した場合の影響 |
| --- | --- | --- |
| [local.rs](../crates/qmb-core/src/expr/local.rs) の StateNode・OperatorNode | 葉の `L::State`・`L::Operator`。StateNode::Act は同じ定義の演算子式も保持する | S・O を直接型引数にしても保持できる。次元・Identity は構築に不要 |
| [evaluate.rs](../crates/qmb-core/src/expr/evaluate.rs) | 局所規則の入力と返却要素の型の対応 | `Act<C, Operator = L::Operator, Term = L::State>` 等の制約を評価メソッドで維持する |
| [composite_eval.rs](../crates/qmb-core/src/expr/composite_eval.rs) | 配置した局所式への評価の委譲 | 同じ関連付けを使う。次元の値は使わない |
| [inner_product.rs](../crates/qmb-core/src/expr/inner_product.rs) | 評価済みの局所状態の内積 | primitive に必要なのは `L::State: InnerProduct<C>`。この縮約自体に Identity や次元は不要 |
| [qmb-normalize](../crates/qmb-normalize/src/lib.rs) | core の内積・期待値に必要な型対応 | 次元・Identity は不要。分割後も正規化の実数・係数能力だけを追加する |

LocalOperator 単体の構文には状態 primitive が現れないが、現在は共通の `L` を通じて状態型も関連付けられている。LocalState は遅延した Act 節点を保持するため、状態 primitive 型だけでは式全体の型を決められない。関連付けの廃止は単に型引数を一つ削る変更ではない。

### 3.2. DIMENSION を直接参照する箇所

製品コードで値を直接読む箇所は次の3箇所である。定数の定義・テスト・doc test はこの数に含めない。

| 箇所 | 処理 | 必要な情報 |
| --- | --- | --- |
| [space.rs の HilbertSpace::dimension](../crates/qmb-core/src/space.rs) | 葉ごとに `L::DIMENSION` を掛け、overflow を検査する | 局所次元。primitive 型・規則は不要 |
| [basis.rs の state_basis_len](../crates/qmb-core/src/basis.rs) | `d` を返す | 局所次元だけ |
| [basis.rs の operator_basis_len](../crates/qmb-core/src/basis.rs) | `d²` を checked arithmetic で返す | 局所次元だけ |

Basis 自体は LocalStructure に依存しない。基底の要素数を得る共通関数が依存している。基底を列挙して数える方式に変える必要はなく、primitive の種類数を次元とみなすこともできない。

### 3.3. 恒等演算子の依存は二種類ある

`LocalOperator::identity()` と `Operator::identity()` は、[局所式](../crates/qmb-core/src/expr/local.rs)・[複合式](../crates/qmb-core/src/expr/composite.rs)の記号的な節点を作るだけである。`L: Identity<C>` も係数の One も不要である。状態への作用では恒等節点をそのまま通過し、積の評価でも単位元として処理する。

一方、[LocalOperator::trace](../crates/qmb-core/src/expr/evaluate.rs) は恒等節点の trace を得るときに `L::identity()` を呼び、その線型結合に primitive の Trace を適用する。[複合 trace](../crates/qmb-core/src/expr/composite_eval.rs) は、演算子のない spectator 因子にもこの局所恒等演算子の trace を使う。Hilbert–Schmidt 内積もこの経路に依存する。総次元を `usize` に整数化して掛ける処理ではない。

具体例では、[SpinHalf](../crates/qmb-spin-half/src/operator.rs) の Identity は単一の I を返す。[三次元代数の fixture](../crates/qmb-core/tests/support/three_level.rs) は対角の matrix unit 三つの和を返し、その基底展開からも Identity を利用する。現在の API は単一の primitive 値を必須にしていない。1.1節の採用案では、この API を O の値の生成へ変更し、選択基底への展開は Basis が担う。

局所 trace は係数の One や整数変換を要求しない。`Tr(I) = d` という数学的等式だけを理由に、次元定数から係数を生成する実装へ置き換えると、`One`・`From<usize>` 等の追加条件や変換方針が必要になる。既存の最小限の型制約を保つには、次元取得と恒等元の表現を別の能力として扱う必要がある。

### 3.4. 空間側と式側に重複する型制約

[HilbertSpace](../crates/qmb-core/src/space.rs) は、次元取得以外では State・Operator・DIMENSION を使わない。Index の保持・列挙・重複検査・順序・grouping・共有・破棄は因子構造の処理である。現実装では `is_compatible_with` と `is_subsystem_of` の引数が同じ `Self` なので、空間側にも同じ L という制約がある。ただし、これは空間側に残す必要がある制約ではない。

[複合式](../crates/qmb-core/src/expr/composite.rs)は `State<I, L, C>`・`Operator<I, L, C>` として L を持ち、その内部の局所式は `L::State`・`L::Operator` を保持する。State::tensor は `IntoIterator<Item = Self>`、加算は同じ Self、作用は `Operator<I, L, C>` を受け取る。**空間の型を HilbertSpace<I> に変更しても、式側のこれらのシグネチャを維持すれば異なる L・primitive 型の式は合成できない。**

[配置と束縛](../crates/qmb-core/src/expr/composite.rs)では、place は Factor かを調べ、bind の validate は Index の存在を調べる。空間が因子構造だけを表すなら、ここでは式の所属先・作用先の構造を確認すればよい。L のない空間を spin-1/2 と spin-1 の式がそれぞれ利用できても、その二つの式を混在させる操作が新たに可能になるわけではない。

[Index の契約](design/index.md)にある因子の同一性と局所定義の区別は、因子構造と型の付いた式の役割分担で維持できる。因子構造だけを tensor する操作には、そもそも局所定義の種類という情報がない。この操作と、異種の代数を含む式の tensor product の対応は別に考える。

### 3.5. 外側へ伝播している依存

| 箇所 | 実際の利用 | 判定 |
| --- | --- | --- |
| [Geometry・Chain](../crates/qmb-geometry/src/lib.rs) | `HilbertSpace<I, L>` を Site 列・全空間として保持し、Index を列挙する | 関連型・次元・Identity を一切参照しない。幾何の処理としては因子構造だけで足りる |
| [spin-half の utility](../crates/qmb-spin-half/src/utility.rs) | Chain の既存 Factor へ spin 式を配置し、全空間へ bind する | 局所定義は構築する式の型で指定できる。配置先の空間へ同じ定義型を付け直す必要はない |
| [qmb-basic](../crates/qmb-basic/src/lib.rs) | `HilbertSpace<Id128, L>` の型別名 | Index の発行自体は既に独立。型別名の追従が必要 |
| [qmb](../crates/qmb/src/lib.rs) | 型・trait の再公開と利用例 | 分割後の trait と、モデル構築・空間参照 API の再公開を追従させる |

## 4. 採用した責任分担と API への影響

### 4.1. 式は S・O を直接指定する

実装した型引数は以下のとおりである。

```rust
LocalOperator<O, C>
LocalState<S, O, C>
Operator<I, O, C>
State<I, S, O, C>
```

演算子式には状態を葉とする節点がないため、S は不要である。状態式は遅延した作用節点を保持するため、S と O の両方を指定する。OperatorPrimitive・StatePrimitive という全能力をまとめた必須 trait は導入しない。

保持・primitive の取り込み・共有・和・スカラー倍・遅延した積／作用の構築には演算能力を要求しない。状態の評価時には、現在の制約を次のように置き換えられる。

```rust
S: Act<C, Operator = O, Term = S>,
O: Adjoint<C, Term = O>,
C: Clone + Conjugate<Output = C>,
```

内積にはさらに `S: InnerProduct<C>` と集計用の係数能力、演算子の積評価には `O: Product<C, Term = O>` を要求する。作用評価に Product を要求しない分離も維持する。

状態式の tensor・和・内積は同じ S・O・C、演算子式の和・積は同じ O・C の入力に限定する。作用は `State<I, S, O, C>` と `Operator<I, O, C>` を組み合わせる。同じ O が異なる S に作用することは、各 S の Act 実装で対応できる。O が構築時の値から推論できない場合には、型注釈・型別名や後続の利用箇所から指定する。

### 4.2. 恒等元は O の値として提供する

係数型から独立した生成能力として、以下を core と facade から公開する。

```rust
pub trait Identity: Sized {
    fn identity() -> Self;
}
```

恒等元を含めても primitive 全体の線型独立性は要求しない。選択基底にその値が含まれない場合は、Basis が基底の primitive の線型結合へ展開する。O が提供する恒等元と Trace から、局所 identity や spectator 因子の trace を取得する。係数の One や整数変換を、局所 trace の追加条件にしない。

旧 `Identity<C>: LocalStructure` は削除し、O の Identity へ移行した。外部の O に利用者が新しい trait 実装を追加できない場合は、提供元の実装か、利用者が所有する wrapper 等で接続する。

式木の記号的な Identity 節点は維持した。恒等 primitive の生成能力なしで遅延構築でき、`Iψ = ψ` と積の単位元を直接処理する既存経路を保つためである。具体的な恒等元は trace で必要になった時に `O::identity()` から得る。

### 4.3. 次元は primitive 型全体の span の次元とする

次元は個々の primitive 値や variant の個数ではなく、その型の要素全体が張る線型空間の次元とする。現在の完全な状態空間では S が d、全演算子空間では O が d² を提供する。spin-half では状態空間が2、演算子空間が4であり、演算子 primitive の6種類という個数とは異なる。

状態基底の要素数は S、演算子基底の要素数は O の span の次元から取得する。`Basis<S, C>`・`Basis<O, C>` の展開・順序番号・列挙は維持する。次元取得にだけ能力を要求し、保持・遅延構築へ波及させない。整数化時の overflow、特に d² が usize に収まらない場合のエラーを維持する。API は `SpanDimension::span_dimension() -> Result<NonZeroUsize, DimensionOverflow>` とした。

O の span の次元 d² と、その作用先の Hilbert 空間の次元 d は区別する。Operator は S を持たないため、演算子単独で作用先の次元が必要になる API は、`space.dimension_with(local_dimension)` に正の局所次元を明示できるようにした。`space.dimension::<S>()` は状態空間、`space.dimension::<O>()` は演算子空間の次元を返す。総次元の整数化を通常の作用・内積・trace の前提にしない。

### 4.4. 追加の定義ラベルを設けない

primitive の意味や規則が違う場合は S・O の型と能力の実装で表す。同じ S・O と同じ能力の実装を、別の L を付けて区別する機構は残さない。全 primitive に定義用の D、独立した定義 ID、実行時の定義オブジェクトを追加することも要求しない。

現在の Product・Act・InnerProduct・Adjoint・Trace は S・O・C で実装が決まり、L を変更するだけで規則を切り替える構成ではない。L 固有だった次元と Identity も primitive 側へ移すため、型の組を指定するために LocalStructure を残す必要はない。同じ primitive 型を使う異なる L の混在だけを拒否する既存の doc test は、新しい型契約に合わせて変更する。

選択基底・表示名・Index・共有方法の違いは定義ラベルの導入理由にしない。異なる primitive 型や数値精度の暗黙の混在・自動変換は導入しない。実行時の表現パラメータや型消去による異種代数の合成は、今回の採用に含めない。

### 4.5. 空間・Geometry・Chain は因子の構造を担う

空間から L を外し、Index・因子順序・grouping・積の構造を扱う。局所定義の識別を代替する型付き空間や、新しい必須の FactorLayout 層を導入する必要はない。式側の型制約が異種の式の混在を拒否する。

place・bind・空間比較は Factor の形、Index の存在・重複、因子順序を確認する。zero・identity・unit でも式の S・O は型として残り、primitive の値がないことを理由に空間へ L を戻す必要はない。

Geometry／Chain から L への依存も外す。モデル utility は既存の Factor に spin の局所式を配置し、全空間へ bind する。Index の再発行や共有構造の全複製はしない。Site の分割・全空間の因子順序・空 Site の契約を維持する。Grading を含む空間の最小構成と Chain の内部表現の詳細は、整理メモの3.4・3.5で詰める。

## 5. 採用理由と不採用にした構成

| 構成 | 判断 |
| --- | --- |
| LocalStructure を縮小して型の組と定義ラベルを残す | 不採用。S・O を直接指定でき、追加の定義ラベルも不要としたため |
| 空間側に L や新しい定義ラベルを残す | 不採用。式の型による混在拒否と責任が重複する |
| Act の関連型から O を導き、保持型に Act を要求する | 不採用。評価に必要な能力が、保持・遅延構築へ波及するため |
| 規則を実行時の Definition／provider に集約する | 今回は導入しない。固定した primitive の規則を個別 trait で提供する構成を維持する |
| S・O を直接指定し、次元・恒等元を primitive 側へ移す | 採用。必要な型対応と演算能力を利用箇所で要求できる |

## 6. 実施した移行と検証対象

1. 次元と恒等元の能力 API を具体化し、SpinHalf・三次元代数・最小係数の fixture を移行した。恒等 primitive の選択基底への展開、span の次元と列挙数、`Tr(I) = d` を確認する。
2. 局所／複合式・節点・評価器・演算子構文・反復的な破棄・qmb-normalize を S・O の直接指定へ移した。同じ primitive 型に別の定義ラベルを付ける検証は廃止し、異なる S・O・C の式を不適切に合成できないことを検証した。
3. 空間・Geometry／Chain・basic の型別名・spin のモデル utility・facade・例を移行した。新しい次元取得 API と、空間比較・配置・bind が因子構造を正しく扱うことを確認する。State／Operator と空間・Chain の後続設計と調整する。
4. 対象テスト、workspace の全 feature 有効／default 無効のテスト、format、Clippy、ドキュメント生成、core 単独の no_std + alloc と依存 feature を確認し、移行完了を記録した。

保持・構築の最小限の型制約、Clone のない primitive と借用値、両精度・独自係数・非正規直交基底、深い木の反復的な評価・破棄、不変な構造の共有を維持する。総次元の overflow は次元取得時に扱い、Index 付きの次元1の Factor と空の tensor の区別も維持する。

## 7. 移行前の調査時の検証

移行前の製品コードと既存テストの依存を調査した。Rust 1.98.0 で次の既存テストを実行し、成功した。これは新しい設計の実装検証ではない。

- `cargo test -p qmb-core --no-default-features -j 2 --quiet`
- `cargo test -p qmb-geometry -j 2 --quiet`

## 8. 実装した API と移行結果（2026-10-26）

全6 crate の製品コード・テスト・doc test・利用例を移行した。旧 LocalStructure と SpinHalf の定義用 marker、両 local_structure module は削除し、互換の定義ラベルや能力をまとめた置換 trait は導入していない。

| 移行前 | 移行後 |
| --- | --- |
| `LocalState<L, C>` | `LocalState<S, O, C>` |
| `LocalOperator<L, C>` | `LocalOperator<O, C>` |
| `State<I, L, C>`・`Operator<I, L, C>` | `State<I, S, O, C>`・`Operator<I, O, C>` |
| `HilbertSpace<I, L>`・`Geometry<I, L>`・`Chain<I, L>` | いずれも Index 型 I のみ |
| `qmb_basic::HilbertSpace<L>` | `qmb_basic::HilbertSpace`（Id128 を使う型別名） |
| `Identity<C>` が恒等元の線型結合を返す | O の `Identity::identity() -> Self` |
| `L::DIMENSION` | S・O の `SpanDimension::span_dimension()`。S は d、O は d² |
| `state_basis_len<L>()` | `state_basis_len<S>() -> Result<usize, DimensionOverflow>` |
| `operator_basis_len<L>()` | `operator_basis_len<O>() -> Result<usize, DimensionOverflow>`。再度二乗しない |
| 引数なしの `space.dimension()` | `space.dimension::<P>()` または `space.dimension_with(NonZeroUsize)` |

標準の局所状態は `LocalState::<SpinHalfState, SpinHalfOperator, Complex<f64>>::primitive(SpinHalfState::Up)` のように指定する。core と qmb root は Identity・SpanDimension を再公開する。既定の正規化 feature と旧 libm alias は維持し、外部依存は追加していない。

Rust 1.98.0 で次を実行し、すべて成功した。

- `cargo test --workspace --all-features -j 2`：doc test を含む108テスト。
- `cargo test --workspace --no-default-features -j 2`：doc test を含む107テスト。
- `cargo fmt --all -- --check`
- `cargo clippy --workspace --all-targets --all-features -j 2 -- -D warnings`
- `RUSTDOCFLAGS='-D warnings' cargo doc --workspace --all-features --no-deps -j 2`
- `cargo check -p qmb-core --no-default-features -j 2`
- `cargo run -p qmb --example heisenberg -j 2` と `--no-default-features` の両構成。期待値は `1.25+0i`、既定構成の norm は1。

core 単独の `cargo tree -p qmb-core --no-default-features -e normal,features` で、数値依存の std・libm と normalize が有効にならないことを確認した。core は引き続き no_std + alloc である。crate 内には LocalStructure の定義・利用・再公開が残っていない。

[primitive の契約テスト](../crates/qmb-core/tests/primitive_contracts.rs)は、異なる primitive 型の式が同じ因子構造を共有して配置・bind できること、数値能力のない係数での保持・遅延構築、d と d² の区別、次元取得にだけ現れる overflow、空の tensor の次元1を検証する。局所／複合式の compile-fail は S・O・C の不適切な混在と Send 制約を検証する。既存の行列照合・非正規直交基底・独自係数・零ノルム・深い式の128 KiBスタック上の評価／破棄・割り当て・Chain の各契約も成功した。

MSRV 1.85 と専用 no_std target の実ビルド、性能ベンチマークは未実施である。State／Operator の一般の写像としての整理、Grading、Chain の内部表現の追加の見直しは、この移行の完了とは別の後続課題とする。
