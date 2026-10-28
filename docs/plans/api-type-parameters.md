# 型引数と係数型の指定を集約する API 案

作成日: 2026-10-27
更新日: 2026-10-28
状態: 2026-10-28 実装済み。以下は採用理由と API の設計記録であり、「提案」の表記を含むコードも、任意の補助候補を除いて製品 API に反映した。

実装先は `qmb_core::system`（crate root から再公開）と `qmb_spin_half::SpinHalf<C>`。公開モデル入口は関連関数へ統一し、旧 utility 自由関数は撤去した。`UtilityError` は spin-half の crate root から公開する。純粋な Chain・with/on の接続・leaves と併せて移行し、正式な契約は[要件](../requirements.md#11-型選択の構築-api)と[Chain の設計](../design/geometry.md)に反映した。

`coefficient`・zero/identity 補助・`MatrixArgs::for_state`・`find_factor` は今回追加しない。当初は `.scaled(C)` と constructor を利用していた。後続の実装では `.scaled(impl Into<C>)` に変更し、MatrixArgs の状態 primitive 型と出力先コンテナは引き続き明示する。式型全体の型引数を束ね直す変更は行わず、戻り値には投影別名を使う。

後続の改善として、`scaled` の引数を `impl Into<C>` に広げ、失敗し得る項を集約する `try_sum` を追加する方針を採用した。実装と検証項目は[式構築 API の実装計画](expression-construction.md)に分離した。2026-10-28に scaled と try_sum、および QuantumSystem の初回範囲の構築・集約メソッドを実装した。

構築・集約の関連関数を QuantumSystem に揃える方針と、配置・積・和・tensor・零・単位元の追加候補は[QuantumSystem の構築 API](quantum-system-construction.md)に記録する。本書冒頭の「補助を追加しない」は当初の実装範囲を指し、後続候補とその追加順は新しい計画で扱う。

## 1. 推奨する方向

状態・演算子 primitive と係数型を関連付ける `QuantumSystem` trait を用意し、primitive を受け取って既存の式木を返す既定メソッドを提供する案を推奨する。spin-1/2 では、`SpinHalf<C>` がその trait を実装する。

```rust
// 実装済みの API。
use qmb_core::QuantumSystem;
use qmb_spin_half::{SpinHalfState, SpinHalfOperator};

type Spin = qmb_spin_half::SpinHalf<num_complex::Complex<f64>>;

let up = Spin::state(SpinHalfState::Up);
let sx = Spin::operator(SpinHalfOperator::Sx);
let result = up.act(&sx);
```

係数型の選択は `Spin` の定義に集約する。状態・演算子の primitive 型は spin-half 側で指定する。各構築結果の型はこの時点で決まり、後から期待値を計算することや戻り値を注釈することを型推論の前提にしない。

`LocalStateOf<Q>`・`StateOf<Q, I>` などの汎用の型別名で、戻り値の型も選択済みの `Q` から書けるようにする。式木・評価器を変更せずに導入できる範囲から着手し、型引数そのものを core 全体で組み替える案は、残る負担を確認してから判断する。

[使いやすさの論点一覧](api-usability.md)のうち、型引数と係数型の指定に対象を絞った検討である。関連型名は `State`・`Operator`・`Coeff` とし、これらは状態 primitive・演算子 primitive・係数の型を表す。構築の入口は `state(primitive)`・`operator(primitive)` とし、`up()`・`sx()` など個別 primitive 用メソッドは今回追加しない。基本 API とモデルの公開名・シグネチャは製品実装へ反映し、任意の補助候補は未採用のままとする。

## 2. 移行前 API で負担が生じる理由

型指定には、役割の違う情報が混在している。

| 情報 | 現行の指定箇所 | 集約の方針 |
| --- | --- | --- |
| 状態・演算子 primitive の組 `S, O` | `LocalState<S, O, C>`、`State<I, S, O, C>` など | `QuantumSystem` の関連型でまとめる |
| 式の係数型 `C` | 状態と演算子の構築、モデル utility の型引数・推論 | `QuantumSystem` の関連型として一度選び、各戻り値へ引き継ぐ |
| Index 型 `I` | `HilbertSpace<I>`、`State`、`Operator` | 配置先や Chain から推論する |
| モデルの実結合定数の型 `R` | `heisenberg`・`xyz` の数値引数 | `C` への変換元として引数から推論する |
| 行列が作用する状態型 `S` | `MatrixArgs<S, B>` | 係数型とは別の論点として、状態からの推論を補助する案を検討する |
| 数値出力先 | `Exported<Vec<_>, _>` などの戻り値注釈 | データ型の選択として明示する。係数の精度選択は引き継ぐ |

現在も、後段で状態と演算子を組み合わせれば `C` を推論できる場合がある。[Heisenberg の利用例](../../examples/src/bin/heisenberg.rs)はその例である。しかし、個別の状態・演算子を独立に構築する際の型指定は残り、コードの一部を切り出した際に後段の型推論へ依存することもある。

状態型から `Act<C>::Operator` を求めれば `O` を省けるが、保持・遅延構築に `Act<C>` の実装を要求することになる。現在は、演算能力を持たない型でも式を保持・構築できるため、この方法だけで core の型引数を削る案は勧めない。[現在の保持契約のテスト](../../crates/qmb-core/tests/primitive_contracts.rs)を維持する。

## 3. 候補の比較

| 案 | 型指定への効果 | 主な代償・限界 | 評価 |
| --- | --- | --- | --- |
| builtin の型別名だけを追加する | `S, O` の繰り返しを省ける | 状態と演算子で `C` を別々に指定する入口が残る | 推奨案の補助として使う |
| `QuantumSystem` と投影した型別名を使う | `S, O, C` を一度選び、構築と戻り値型の両方に引き継げる | primitive の組を表す概念の導入と、その責任を定義する必要がある | 第一候補 |
| `SpinHalf<C>` 固有の関連関数だけを追加する | 代数と `C` を一度選び、状態・演算子・モデルに引き継げる | 型の関連付けが具体的代数の中に留まり、汎用コードから参照できない | 共通 trait を設ける前の案 |
| `let spin = SpinHalf::<C>::new()` のメソッドを使う | 関連関数の案と同じ型指定削減ができる | 保持する情報がない場合も構築用の値を管理する | 値を渡して構築処理を抽象化したい場合の候補 |
| 代数・係数型を束ねた型を core の式型引数にする | `State<I, F>` など、明示的な式型も短くできる | 公開型・trait 制約・演算実装・export まで広く変更する | 明示的な型の負担が残る場合に再検討 |
| 係数型を評価時まで決めない | 初期の primitive 構築で `C` を書かずに済む | 構築時に保持する係数の型、精度、演算の整合性を再設計する必要がある | 今回の改善には採用しない |

構築用の値は、将来、実行時に変わる設定を保持する API には適する可能性がある。現在の対象では代数と係数型を型で選ぶだけなので、値の生成が不要な関連関数を第一候補とする。将来の設定や独自代数への適用を、この入口に先回りして追加しない。

core の型引数を束ねる案では、束ねた型自体を式の同一性に含めるか、関連型を通じて既存の式へ展開するかも決める必要がある。前者では同じ `S, O, C` を指定した別の束ね方が合成できなくなる可能性があり、後者では関連型から元の束ねた型を推論できるかが論点になる。構築時の負担を減らすために、この判断まで必須にする必要はない。

## 4. 推奨案の利用イメージ

### 4.1. 型の関連付けと共通の構築メソッド

```rust
// 提案 API。State・Operator は primitive 型を表す。
pub trait QuantumSystem: Sized {
    type State;
    type Operator;
    type Coeff;

    fn state(primitive: Self::State) -> LocalStateOf<Self> {
        qmb_core::LocalState::primitive(primitive)
    }

    fn operator(primitive: Self::Operator) -> LocalOperatorOf<Self> {
        qmb_core::LocalOperator::primitive(primitive)
    }
}

pub type LocalStateOf<Q> = qmb_core::LocalState<
    <Q as QuantumSystem>::State,
    <Q as QuantumSystem>::Operator,
    <Q as QuantumSystem>::Coeff,
>;

pub type LocalOperatorOf<Q> = qmb_core::LocalOperator<
    <Q as QuantumSystem>::Operator,
    <Q as QuantumSystem>::Coeff,
>;
```

`SpinHalf<C>` の実装では、次のように関連型を指定する。

```rust
impl<C> QuantumSystem for SpinHalf<C> {
    type State = SpinHalfState;
    type Operator = SpinHalfOperator;
    type Coeff = C;
}
```

trait の既定メソッドを使い、`Spin::state(SpinHalfState::Up)`・`Spin::operator(SpinHalfOperator::Sx)` と構築する。これらの呼び出しでは trait の import が必要になる。

trait は関連型の選択を表し、値として状態・演算子・係数を保持するわけではない。`Act`・`Product`・数値演算などの supertrait は課さず、実際の評価メソッドで対応する能力を要求する。`Q: QuantumSystem` だけで任意の評価ができることは保証しない。

### 4.2. 局所式と 2 因子の系

```rust
// 提案 API。配置・合成・評価は現在の core API を使用する。
use num_complex::Complex;
use qmb_core::{HilbertSpace, QuantumSystem};
use qmb_spin_half::{SpinHalf, SpinHalfState as Ket, SpinHalfOperator as Op};

type Spin = SpinHalf<Complex<f64>>;

let a = HilbertSpace::factor("a");
let b = HilbertSpace::factor("b");
let space = a.tensor_product(&b)?;

let psi = Spin::state(Ket::Up).at(&a)?
    .tensor_product(&Spin::state(Ket::Down).at(&b)?)?;
let zz = Spin::operator(Op::Sz).at(&a)?
    .product(&Spin::operator(Op::Sz).at(&b)?)?
    .bind(&space)?;

let value = psi.expectation(&zz)?; // -0.25 + 0i
```

すべての primitive を `state(primitive)`・`operator(primitive)` の同じ入口から構築する。zero・identity の補助を追加する場合も選択済みの `C` を持つ式を返す案だが、今回の基本の入口とは別の検討事項とする。

### 4.3. モデル構築の関連関数（合意済みの方針）

現在 `qmb_spin_half::utility` にある `xyz`・`heisenberg`・`neel` は、新しく実装する `SpinHalf<C>` の関連関数とする。いずれも `self` を受け取らず、係数型は `SpinHalf<C>` の `C` から決める。これらのモデル固有の関数を `QuantumSystem` trait へ追加しない。

Chain と Site 順の局所空間列 `&[HilbertSpace<I>]` を受け取る。各要素は葉に限定せず、`leaves()` 相当から配置先の葉を取得できるものとする。モデル関数への全空間 `&space` の引数は不要とする。呼び出し方は次のとおりとする。

| 関連関数 | 引数 | 成功時の戻り値 |
| --- | --- | --- |
| `Spin::xyz` | `&chain, &factors, jx, jy, jz` | `OperatorOf<Spin, I>` |
| `Spin::heisenberg` | `&chain, &factors, coupling` | `OperatorOf<Spin, I>` |
| `Spin::neel` | `&chain, &factors, first` | `StateOf<Spin, I>` |

すべて失敗を表す `Result` を返す。`factors.len() == chain.site_count()` と、全要素の葉にわたる Index の重複がないことを検証する。一要素を一 Site の局所空間とし、列順を Site 順とする。複数葉・空・単項 Tensor も入力として受け入れ、配置先は中の葉から取得する。葉数による制限は設けない。ペアがない 0・1 Site でも検証を省略しない。

XYZ・Heisenberg は `nn_pairs_with(&factors)` を使い、最近接 Site 間の葉の全組み合わせへ結合を与える。Site 内の葉同士には結合を自動追加しない。Neel の向きは Site ごとに交互に決め、その Site の全葉に同じ向きを割り当てる。空 Site は位置を保持し、その Site から項・状態を生成せず、前後を直接接続し直さない。詳細は [空間・幾何の計画第3.2節](space-geometry-connection.md#32-モデルは全空間ではなく-factor-列を受け取る)に記録する。

全空間は各局所空間のハンドルを共有し、列順に tensor してモデル側で構築する。内部の grouping はそのまま保持する。XYZ・Heisenberg はその空間に束縛した演算子を返し、Neel 状態も Site 順に葉をたどった因子順序で構築する。Index の Clone・再発行は不要で、現在と同じ `I: Eq` を要求する。0 Site は単位空間に束縛した zero Hamiltonian と単位状態を返す。

Chain と HilbertSpace の責任分担は [接続・保持構造の検討](space-geometry-connection.md)に記録する。Chain はサイト数と境界条件だけを保持し、外部の Factor 列から近接ペアを列挙する。モデルは元の全空間を受け取らないため、その root・grouping や列順と異なる全空間順序を保存する契約は持たない。既存の状態や行列を暗黙に並べ替える機能は追加しない。

```rust
// 提案 API。QuantumSystem / SpinHalf と、純粋な Chain の案を組み合わせる。
use num_complex::Complex;
use qmb_basic::Id128;
use qmb_core::HilbertSpace;
use qmb_geometry::{Boundary, Chain};
use qmb_normalize::NormalizedExpectation;
use qmb_spin_half::{SpinHalf, SpinHalfState};

type Spin = SpinHalf<Complex<f64>>;

let chain = Chain::new(6, Boundary::Open);
let factors: Vec<_> = (0..6)
    .map(|_| HilbertSpace::factor(Id128::default()))
    .collect();
let h = Spin::heisenberg(&chain, &factors, 1.0)?;
let psi = Spin::neel(&chain, &factors, SpinHalfState::Up)?;
let energy = psi.normalized_expectation(&h)?;
```

`C` は構築用の型から、`I` は Factor 列から決まる。結合定数は現在と同じ実数型 `R` の引数として受け取り、型制約を各関数へ限定する。

- `xyz` は `C: Clone + From<R> + Neg<Output = C>` を要求する。
- `heisenberg` は同じ `C` の能力に加え `R: Clone` を要求し、等方的な XYZ として構築する。
- `neel` は係数の数値演算や `Clone` を要求しない。

`impl<C> SpinHalf<C>` 全体にこれらの能力を課さない。既存の utility に実装を委譲するか、別の実装用モジュールへ移すかは公開の呼び出し方から独立に決める。

単精度へ切り替える際は `Spin` の定義を変更する。構築用の型に既定の `f64` を設定せず、意図しない精度選択を避ける。明示的に `f64` とした変数を `Complex<f32>` へ暗黙に変換する API は追加しない。

### 4.4. 明示的な戻り値の型

複合式も関連型を投影して、既存 core 型を指す別名を提供する。

```rust
pub type StateOf<Q, I> = qmb_core::State<
    I,
    <Q as QuantumSystem>::State,
    <Q as QuantumSystem>::Operator,
    <Q as QuantumSystem>::Coeff,
>;

pub type OperatorOf<Q, I> = qmb_core::Operator<
    I,
    <Q as QuantumSystem>::Operator,
    <Q as QuantumSystem>::Coeff,
>;
```

例えば、複数の関数で同じ精度を使うアプリケーションでは、次のように選択をモジュールに集約できる。

```rust
// 提案 API を使用する側の型別名。
type Spin = qmb_spin_half::SpinHalf<num_complex::Complex<f64>>;
type State<I> = qmb_core::StateOf<Spin, I>;
type Operator<I> = qmb_core::OperatorOf<Spin, I>;
```

具体的な精度と primitive の組の選択を `Spin` の定義に集約し、戻り値型では `Spin` と Index 型だけを指定できる。汎用コードでも `Q: QuantumSystem` と `StateOf<Q, I>` を使える。ただし関連型の投影は逆向きの型推論を保証しないため、既存の式から元の `Q` を推論できることは前提にしない。入口の `Q` は明示的に選ぶ。

### 4.5. スカラー値と数値出力

スカラー値の型は `.scaled(...)` などの引数から推論できる。実数から `C` への変換を明示したい場合には、`Spin::coefficient(value)` を `C: From<R>` の範囲で提供する案がある。この補助は任意とし、数値変換の規則を追加しない。

数値出力は、現在の `Exported<Vec<_>, _>`・`Exported<Array2<_>, _>` などの型指定で係数型を引き継げる。出力先コンテナの選択は残す。

行列 export の `MatrixArgs::<SpinHalfState, _>` は別の型指定である。状態が既にある場面には、`MatrixArgs::for_state(&psi, &basis)` のように状態型を推論する補助を `qmb-export` 側に追加する候補がある。この場合、状態は型推論に使い、作用空間は現在と同じ束縛・明示指定に従う、と契約を明示する。状態から自動的に作用空間も選ぶかは別に判断する。

この補助のために spin-half から export への依存を増やす必要はない。局所状態にも同じ補助を提供するか、状態が存在しない行列出力をどう案内するかは、export API の検討事項として残す。

## 5. 実装する場合の境界と契約

構築用の型は、`PhantomData<fn() -> C>` などで係数型だけを指定する。インスタンスの生成を必要とせず、各関数が既存の constructor・utility を呼ぶ構成にできる。

| 項目 | 推奨案での扱い |
| --- | --- |
| 公開先 | `QuantumSystem` と投影した型別名を `qmb-core`、`SpinHalf<C>` とその実装を `qmb-spin-half` で提供する案。具体的な module 配置は実装時に決める |
| 戻り値 | 既存の core の式型。新たな式ラッパーは保持しない |
| 構築時の制約 | `state(primitive)`・`operator(primitive)` に数値演算や `Clone` を要求しない |
| 評価時の制約 | 現在の式メソッドが必要とする能力を維持する |
| モデル構築 | `xyz`・`heisenberg`・`neel` を `SpinHalf<C>` の関連関数とし、Chain と Site 順の局所空間列を受け取る。配置先は各要素の葉から取得し、その操作に必要な能力だけを要求する |
| 空間と Index | Factor 列の `HilbertSpace<I>` を共有して空間を構築する。構築用の型や Chain に登録・保持しない |
| 同一性 | 構築用の型や値の識別子を式へ追加しない。同じ `S, O, C` なら従来どおり合成できる |
| 数値精度 | `C` は構築時に確定する。異なる精度の式の暗黙変換は追加しない |
| 独自代数 | `QuantumSystem` の関連型を指定すれば共通 constructor を使える。既存 core API の直接利用も維持する |
| 依存関係 | 本案のための新たな production dependency は不要。spin-half の `no_std` を維持する |

本案は、過去の LocalStructure が担っていた primitive の組の関連付けを、係数型とともに明示的な trait として再導入する側面がある。その責任を「関連型の選択と共通の式構築」として定義し直す。局所規則・次元・恒等元の数学的定義は、現在と同じ primitive と個別能力が担う。

`QuantumSystem` という名前は全 Hilbert 空間やモデルの所有も連想させるため、正式な命名では責任の範囲を確認する。現案の trait には空間・幾何・Index 発行を含めない。

## 6. 実装前の試作による確認と限界

Rust 1.98.0 を使い、リポジトリ外の検討用パッケージから現行 crate を参照して、関連関数と値を使う構築 API の双方を試作した。製品ソース・依存設定は変更していない。初回の試作は個別 primitive 用メソッドを含む旧案で、`cargo test --offline --locked --jobs 1` に検証用の target directory 指定を加え、9 件のテストと 2 件の compile-fail doctest が成功した。

確認した範囲は次のとおりである。

- `Complex<f32>`・`Complex<f64>` の両方で、一度選んだ `C` から局所式と Heisenberg・Neel を構築できる。
- 両精度で、接尾辞のない実数リテラル `1.0` を結合定数として使える。任意のユーザ定義係数でも同じ推論が必ず成立することまでは保証しない。
- 2 因子の期待値が `-0.25`、6 Site の Neel 状態の Heisenberg エネルギーが `1.25` になる。
- 密・疎行列 export に選択済みの係数型を引き継げる。状態から行列引数の型を推論する補助関数も試作できる。
- 数値能力を持たない係数型、非 `Clone` の primitive でも保持・遅延構築できる。
- 作成元の値の生存期間に依存せず、既存の式を返せる。
- 単精度の状態と倍精度の演算子の作用は、関連関数・値の両案でコンパイル時に拒否される。

試作ライブラリは `#![no_std]` とし、ホスト上でコンパイルした。MSRV 1.85 と専用の `no_std` target では検証していない。公開 API 全体の実装、利用者による評価、性能測定は未実施である。

その後、`QuantumSystem` の既定 constructor と投影した 4 種類の型別名を追加試作した。数値・代数能力のない型での遅延構築、選択した係数型での既存演算、同じ関連型を持つ異なる trait 実装間での式の合成を、追加の対象テスト 3 件で確認した。

## 7. 実装へ進む場合の範囲

1. `QuantumSystem` と投影した型別名、`SpinHalf<C>` の実装を追加し、primitive の組と係数型を一箇所で選ぶ利用例を用意する。
2. primitive の構築は `state(primitive)`・`operator(primitive)` を使い、数値能力の制約を各操作へ限定する。zero・identity の補助は別途範囲を決める。
3. `xyz`・`heisenberg`・`neel` を `SpinHalf<C>` の関連関数として提供し、Chain と Site 順の局所空間列を受け取る。第4.3節の構築・重複検証・複数葉と空 Site・列順の契約に合わせて移行する。
4. 基底展開・数値出力まで一連の利用コードを確認し、型の再指定がどこに残るか記録する。
5. 構築・精度不一致・既存 core 型との接続に回帰検証を付ける。
6. 戻り値型を書く負担や独自代数での繰り返しがなお大きければ、core の型引数を束ねる案を再検討する。

既存の式 constructor は維持した。モデル構築の公開入口は `SpinHalf<C>` の関連関数へ統一し、旧 utility の自由関数は撤去した。要件・設計文書へ正式な契約を反映済みである。

## 8. 製品実装の検証

回帰テストは core の `tests/systems.rs`、geometry の `tests/geometry.rs`・`tests/space.rs`、spin-half の `tests/systems.rs`、examples の `tests/chain.rs` に置いた。非 Clone primitive・数値能力のない係数の構築、同じ関連型を持つ別の系との合成、既存 core 型との接続、f32/f64、深い入れ子、複数葉・空 Site、借用の同一性、Index 重複拒否を確認する。精度不一致は SpinHalf の compile-fail doctest で拒否する。

基底展開の bin は一つの Spin 別名から状態・演算子を構築し、Heisenberg の bin はその係数選択を密疎 export まで引き継ぐ。MatrixArgs の状態 primitive 型、出力コンテナ、スカラー値は引き続き指定する。後段の期待値や戻り値注釈なしでも入口で式の係数型が確定する。

製品実装は Rust 1.98.0 の devShell で、`cargo test --workspace --all-features --offline --locked`、`cargo test --workspace --no-default-features --offline --locked`、全 target・全 feature の Clippy（`-D warnings`）、format、ドキュメント生成、および core・geometry・spin-half の `--no-default-features` check を確認した。実行例では基底座標 `[5, -1]` と6 Site のエネルギー `1.25`、密疎 export の一致を確認した。MSRV 1.85 と専用 no_std target の検証は今回も未実施であり、ホストでの no_std crate の確認にとどまる。
