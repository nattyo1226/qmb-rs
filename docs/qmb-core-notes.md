# qmb-core 共通式・評価の実装検討メモ

更新日: 2026-10-22

## 1. 位置付けと提供範囲

本書は、[実装計画](implementation-plan.md)の Step 2「局所式と評価」に向けた一時的な作業用記録である。恒久的な仕様書として維持することは想定せず、既存の恒久ドキュメントからはリンクしない。Step 2 の局所式・評価と API 移行は実装済みであり、後続段階の機能と区別する。

局所式の数学的契約は、[局所 primitive・式木の数学的能力](design/primitive-expression.md)と[要件](requirements.md)に従う。共通の代数処理は `qmb-core` に集約し、局所式・評価に加え、後続の Index・空間・複合式・基底処理も core に置く。具体的代数・Geometry・数値 export はそれぞれの crate に分離する。

crate の責任、公開する式型、共通の構文パーツ、`Arc` による共有、作用の `act` と評価の `evaluate` の役割は確定事項とする。`LinearCombination<P, T>` は薄い有限線型結合型として維持し、演算 trait は `T` を型引数、返却要素を関連型 `Term` とする。局所規則・局所基底展開の結果は、式木側で `map` して取り込む。これらの API 整備と接続を Step 2 で行う。以下は Step 2 で実装した採用済みの設計であり、後続段階の設計事項は第6節にまとめる。

局所状態・演算子の構築、加算、係数乗算、zero・identity、演算子積、作用、内積、随伴、局所 trace を実装する。これらから作れる commutator・anticommutator、局所期待値、Hilbert–Schmidt 内積も含める。

norm・正規化・bra は、実数値の扱いと API を別途検討してから追加する。既存の局所基底について、primitive 一つの展開は線型結合を返す契約を維持し、局所式への取り込みを Step 2 に含める。Index・空間・配置は Step 3、式全体の基底展開・複合基底・座標取得は Step 5 に置く。これは初期版全体からそれらの機能を除外するものではない。

`qmb-core` の `no_std + alloc` を維持する。`Arc`・`Vec` など既存の標準機構で実装し、Step 2 では新しい production 依存や feature は追加しない。

core 内を局所契約、共通の構文パーツ、式・構築構文、空間、評価、基底処理に分け、主要な公開型を crate root から再公開する。Step 2 では局所式と評価に必要な部分を追加し、具体的代数は検証に利用する。通常依存へ具体的代数・Geometry・数値 backend を持ち込まない。

source の自動保持、式の構造比較、`Debug`・`Display` を含む式表示は後続へ延期する。現在の式と source を別々に保持する二重のルートや、それらの参照・比較 API は導入しない。操作は入力を変更せず、利用者は必要なら評価前の式を変数に保持できる。不変実体の共有は維持する。

### 1.1. 実装・検証状況（Step 2 完了）

2026-10-22 に Step 2 の実装と検証を完了した。

1. `LocalStructure` と `Identity<T>` を分離し、演算 trait を `T` と関連型 `Term`・`Act::Operator` の契約に移行した。`StateBasis`・`OperatorBasis` は定義・再公開・実装宣言・コード内の利用を削除した。
2. 公開の共通パーツ、薄い線型結合と `map`、`Term` の型別名化、局所式の constructor と所有値・借用値の演算子構文を実装した。spin-half と非公開テスト用の三次元代数・基底も移行した。
3. `Arc` による部分式の共有、基底を指定しない状態評価・作用・内積・期待値・随伴・積・trace・導出操作を実装した。評価・破棄は明示的な作業スタックで処理する。状態への作用評価は `Product` 能力を要求しない。
4. spin と三次元代数の式を単精度・倍精度で独立した行列と照合し、非正規直交基底、型制約、共有、項順序・重複・零係数、空・単一項の無確保、2万段の式を 128 KiB のスタックで評価・破棄するケースを検証した。
5. 既存の direnv 環境（Rust 1.98.0）で workspace の全46テスト（doc test 含む）、format check、Clippy、警告をエラーとするドキュメント生成の成功を確認した。通常依存の `std`・`libm` を有効にせず、新しい依存や feature を追加していない。MSRV（Rust 1.85）での実行と `std` のない target 向けビルドは未確認である。

`Act::act` は状態を受け手として局所規則を計算し、線型結合を返す。状態式の inherent method `act` は遅延した作用を構築し、`evaluate` が必要な規則を呼んで通常の状態式を返す。両者の役割を区別する。

## 2. 式の型・内部表現・共有

### 2.1. 公開する式型

`LocalState<L, T>` と `LocalOperator<L, T>` の二型を採用する。`L: LocalStructure` から primitive 型を取得し、局所構造と係数精度の整合性を型で保証する。`LocalStructure` 自体は係数型を持たず、状態・演算子の関連型と次元だけを定める。演算能力は必要なメソッドに付け、式の保持・構築に数値制約を波及させない。

式ハンドルが `Arc<StateNode<L, T>>` または `Arc<OperatorNode<L, T>>` を保持し、enum が節点の種類を表す。和・スカラー倍・積・作用・随伴の各 variant の内容は `ScaleExpr`・`SumExpr`・`ProductExpr`・`ActExpr`・`AdjointExpr` の共通パーツが保持する。これらは公開 module `qmb_core::expr::parts` に置き、演算能力の `Product`・`Act` trait と名前を区別する。式ハンドルは crate root から再公開する。節点 enum と共有の内部表現は非公開とする。

将来の複合状態・演算子とは公開型を分ける。共通式木を使う契約は、局所式と複合式の公開型を一つに統一することまでは要求しない。

### 2.2. Arc による所有権と不変共有

`Arc` は不変な部分式を複数箇所から共有するために用いる。例えば `&a + &a` の二つの operand は同じ `a` の実体を参照し、`a.scaled(c)` は新しい係数乗算の節点から既存の `a` を参照する。新しい式を構築しても元の式は変更せず、木全体や primitive の deep copy を行わない。

primitive は節点内に所有でき、節点の共有によって primitive も共有する。primitive ごとに別の `Arc` を置くことは必須としない。式の `Clone` は手動実装などで `Arc` の共有参照だけを複製し、primitive・`L`・`T` の `Clone` を要求しない。評価中は借用を優先し、不要な `Arc::clone` を避ける。

共有は同じ式ハンドルの再利用によって成立する。別々に構築した同形の式を自動検出して共有する機構は導入しない。概念上は式木であり、共有した実体は DAG になり得る。これは source の自動保持や評価結果のキャッシュとは別の仕組みである。

primitive に一律の `Clone`・`Copy`・`Send`・`Sync`・`'static` を要求しない。借用を含む primitive も、その型の lifetime の範囲で利用できるようにする。係数は節点内に保持し、数値演算の制約は必要な操作に付ける。

`Send + Sync` は保持する型に応じて成立させる。非 `Send` な primitive も格納でき、その場合は式のスレッド間共有ができない。将来のスレッド間共有への接続を重視して `Arc` を採用するが、並列 evaluator の実装は後続へ延期する。

共有の内部表現は非公開とし、`Rc`／`Arc` 切替 feature は設けない。`Arc` は参照の複製・破棄に atomic 操作を用いるため、`Rc` より参照カウントのコストがある。借用による読み取りのたびに参照カウントを操作する必要はない。

長い式の破棄は、最後の所有者になった節点を明示的な作業スタックで反復的に解放する。共有先が残る節点は、その共有を維持する。式の深さに比例した destructor の関数再帰に依存しない。この保証は core が管理する式節点を対象とし、ユーザ定義の要素自身の内部構造には及ばない。

`alloc::sync::Arc` は `no_std + alloc` で利用できるが、ポインタ幅の atomic 操作を持つ target が必要となる。これらの制約は既存環境の Rust 1.98.0 の `alloc` 実装で確認した。性能比較は未実施である。

### 2.3. 共通の構文パーツと薄い線型結合型

式木に共通する構造を、要素型に依存しないパーツとして切り出す。構造の保持・構築・線型部分の走査を共用し、積・作用などの数学的な評価規則は各式の評価実装が担当する。

| パーツ | 保持するもの |
| --- | --- |
| `ScaleExpr<E, T>` | `Complex<T>` の係数と operand `E` |
| `SumExpr<E>` | 空・単一・複数の順序付き要素 |
| `ProductExpr<E>` | 左右二つの operand。積の順序を保持する |
| `ActExpr<S, O>` | 作用先の状態 `S` と演算子 `O` |
| `AdjointExpr<E>` | 随伴を取る operand `E` |

共通パーツは保持内容を公開する単純なデータ型とし、式木に依存せず利用できるようにする。

```rust
pub struct ScaleExpr<E, T> {
    pub coefficient: Complex<T>,
    pub operand: E,
}

pub enum SumExpr<E> {
    Empty,
    One(E),
    Many(Vec<E>),
}

pub struct ProductExpr<E> {
    pub lhs: E,
    pub rhs: E,
}

pub struct ActExpr<S, O> {
    pub state: S,
    pub operator: O,
}

pub struct AdjointExpr<E> {
    pub operand: E,
}
```

`LinearCombination<P, T>` は、要素 `P` と `Complex<T>` の係数からなる有限項列だけを表す薄い型として維持する。式木・局所構造・評価規則に依存させない。共通パーツを使った保持内容は次の形とし、既存の公開 API も Step 2 で更新した。

```rust
pub struct LinearCombination<P, T> {
    terms: SumExpr<ScaleExpr<P, T>>,
}
```

項順序・重複項・零係数を保持し、空の項列は零を表す。空・単一要素でコンテナ自体の heap allocation が不要という性質を維持する。要素や式節点の確保は別に数える。

既存の項型 `Term<P, T>` は次の型別名に統合し、crate root の再公開を維持する。要素フィールド名は `primitive` から `operand` へ揃え、core・spin-half・テスト・doc test の利用箇所を Step 2 で更新する。

```rust
pub type Term<P, T> = ScaleExpr<P, T>;
```

この項型は係数と要素の組であり、演算 trait の関連型 `Term` は係数を含まない要素型を指す。両者を区別する。`terms()`・`into_terms()`・`from_terms()`、`FromIterator`・所有値と借用値の `IntoIterator` を維持し、走査の実装は `SumExpr` 側へ集約する。既存 iterator の順序・両端走査・正確な長さ・終了後の振る舞いも維持する。

要素だけを変換する所有権付きの `map` を追加する。

```rust
impl<P, T> LinearCombination<P, T> {
    pub fn map<Q>(self, f: impl FnMut(P) -> Q) -> LinearCombination<Q, T>;
}
```

`map` は各要素を項順に一度ずつ変換し、係数・項数・順序・重複・零係数を保持する。要素と係数を移動するため、この操作自体に `Clone` や数値演算の制約は不要である。

`SumExpr` は係数を付けない要素列、`ScaleExpr` は係数と対象の組とする。局所式ではそれぞれを `Sum`・`Scale` variant の内容に用い、線型構造を明示する。primitive の線型結合を不透明な葉に格納する `PrimitiveSum` や、別の `LinearExpr` 層は設けない。式とは別の汎用出力 trait も追加しない。

### 2.4. 局所式の節点

以下が採用する構成である。コードでは `L: LocalStructure` の制約を省略する。`parts` は `qmb_core::expr::parts` を指す。内部のハンドル保持方法は、反復的な破棄を可能にする範囲で調整できる。

```rust
pub struct LocalState<L, T>(Arc<StateNode<L, T>>);
pub struct LocalOperator<L, T>(Arc<OperatorNode<L, T>>);

enum StateNode<L, T> {
    Primitive(L::State),
    Scale(parts::ScaleExpr<LocalState<L, T>, T>),
    Sum(parts::SumExpr<LocalState<L, T>>),
    Act(parts::ActExpr<LocalState<L, T>, LocalOperator<L, T>>),
}

enum OperatorNode<L, T> {
    Primitive(L::Operator),
    Identity,
    Scale(parts::ScaleExpr<LocalOperator<L, T>, T>),
    Sum(parts::SumExpr<LocalOperator<L, T>>),
    Product(parts::ProductExpr<LocalOperator<L, T>>),
    Adjoint(parts::AdjointExpr<LocalOperator<L, T>>),
}
```

zero は空の `Sum`、演算子の恒等元は `Identity` で表す。内積・trace は scalar を求める操作とし、今回は scalar 用の式木を導入しない。

### 2.5. 演算規則の返却型と式木への接続

演算 trait の出力は `LinearCombination<Self::Term, T>` に限定する。係数精度 `T` は trait の型引数、返却要素は関連型 `Term` とし、汎用の `Output` 関連型や返却要素に `Primitive` という名前は用いない。局所規則では `Term` を対応する primitive 型とし、複合式の演算規則でも適切な要素型を指定して同じ契約を利用できる。演算 trait 自体は `LocalState`・`LocalOperator`・`LocalStructure` に依存しない。

| 局所規則 | Step 2 での返却形式 |
| --- | --- |
| 状態 primitive の Act | `LinearCombination<L::State, T>` |
| 演算子 primitive の Product・随伴 | `LinearCombination<L::Operator, T>` |
| 局所構造の恒等演算子 | `LinearCombination<L::Operator, T>` |
| 状態 primitive の基底展開 | `LinearCombination<L::State, T>` |
| 演算子 primitive の基底展開 | `LinearCombination<L::Operator, T>` |
| 内積・trace | 従来どおり `Complex<T>` |

実装した trait の契約は以下のとおりとする。

```rust
pub trait Act<T> {
    type Operator;
    type Term;

    fn act(&self, operator: &Self::Operator) -> LinearCombination<Self::Term, T>;
}

pub trait Product<T> {
    type Term;

    fn product(&self, rhs: &Self) -> LinearCombination<Self::Term, T>;
}

pub trait Adjoint<T> {
    type Term;

    fn adjoint(&self) -> LinearCombination<Self::Term, T>;
}

pub trait LocalStructure {
    type State;
    type Operator;

    const DIMENSION: NonZeroUsize;
}

pub trait Identity<T>: LocalStructure {
    fn identity() -> LinearCombination<Self::Operator, T>;
}
```

`impl<T>` 側で代数・演算ごとの必要な数値制約を付ける。Primitive に係数型を持たせたり、係数型付きの演算用ラッパーを導入したりする必要はない。`LocalStructure` の関連型には演算能力の制約を置かず、評価メソッドなどの利用箇所で能力と入出力の対応を要求する。例えば状態の作用規則は `L::State: Act<T, Operator = L::Operator, Term = L::State>` とする。恒等演算子の規則は `L: Identity<T>` として必要な箇所だけで要求する。

局所規則と恒等演算子は、同じ局所構造の primitive の有限線型結合を返す。未評価の作用・積・随伴を返して計算を先送りしない。指定基底への展開や同類項の結合は要求しない。

式木への取り込みでは、線型結合の `map` で要素を対応する式へ変換し、各項を `Scale`、項列を `Sum` として組み立てる。例えば `Sx` の `Up` への作用は、まず係数 `1/2` と `Down` の一項を返し、式木側が `Scale(1/2, Primitive(Down))` にする。変換は共通 constructor に集約し、不透明な線型結合の葉は追加しない。

```rust
let terms = primitive_state.act(&primitive_operator);
let result = LocalState::from_linear_combination(
    terms.map(LocalState::primitive),
);
```

基底展開も線型結合を返し、各要素は選択した基底の primitive に限る。入力は primitive 一つとし、式全体への展開は core の共通処理が担当する。Step 2 では `Basis<P, T>` を維持し、`StateBasis`・`OperatorBasis` は削除する。利用側で `L: LocalStructure` と、状態基底には `B: Basis<L::State, T>`、演算子基底には `B: Basis<L::Operator, T>` を直接要求する。基底情報や特別な展開済み状態は結果に付随させない。

```rust
// Basis<L::State, T> の展開部分
fn expand(&self, primitive: &L::State) -> LinearCombination<L::State, T>;

// Basis<L::Operator, T> の展開部分
fn expand(&self, primitive: &L::Operator) -> LinearCombination<L::Operator, T>;
```

core の marker trait 定義・再公開、spin-half とテスト用局所代数の実装宣言、テスト・doc test の型制約を同じ Step 2 で更新し、互換用の marker trait や別名は残さない。`Basis` の `index() -> Option<usize>` と、所有要素を遅延列挙する `elements()` は維持する。次元だけを使う `state_basis_len<L>()`・`operator_basis_len<L>()` からは不要になる型引数 `T` を外し、後者の二乗時の overflow 検証を維持する。基底の完全性・線型独立性と、要素数が状態では `d`、演算子では `d²` となる数学的契約・検証条件も維持する。

spin-half の `SzBasis`・`CartesianBasis`・`LadderBasis` は線型結合を返す形式を維持し、core の非公開テスト用局所代数・非正規直交基底とともに局所式への取り込みを検証する。演算 trait の関連型への移行と線型結合・局所式の接続は Step 2 に含める。

線型結合を返す段階では式節点の確保は不要であり、式木への取り込み時に `Arc` による確保が生じる。`SumExpr` のコンテナに関する空・単一要素の割り当て特性を、式全体の保証には拡張しない。`map` のコンテナ再確保を省けるかは要素型に依存するため、無確保を保証しない。性能比較・割り当ての最適化は別途検証する。

### 2.6. 複合式への接続

初期版では、一つの系に含まれる全基本因子を同じ局所構造型 `L` に揃える。Step 3 でも型付き局所式を共有して保持し、同じ `L` の局所式と Index を組にして複合式の葉とする。複合式の葉から必要な操作を局所 evaluator へ静的に委譲する。

異種局所構造の混在・合成と、そのための型消去 API・保持機構は後続へ延期し、Step 2・3 の設計対象に含めない。独自代数は `L` の別実装として別の系で単独利用できる。同種因子の入れ子や一つの Site への複数自由度の配置は維持する。

将来の異種構造への対応や、型消去後の `Send + Sync` の扱いは、その機能を導入する段階で設計する。初期版のために `'static` 等の追加制約を先に課さない。

複合式も共通の構文パーツを再利用し、配置・束縛・状態の tensor product などの固有の構造を組み合わせる。所属空間・束縛情報の保持と整合性判定は複合式側で担当する。具体的な節点構成は Step 3 で設計する。

## 3. 公開 API と構築規約

### 3.1. メソッドと演算子構文

名前付き constructor・メソッドを基本 API とし、演算子構文は同じ処理への入口とする。primitive の取り込みは局所規則を実行せず、`Primitive` 節点を構築する。局所規則・基底展開は第2.5節のとおり線型結合を返し、要素を `map` した後に共通 constructor で式木へ取り込む。

primitive の取り込みは所有値を受け取り、式同士の操作は `&self` と共有借用を基本とする。変更されない operand は `Arc` で共有する。名前付き API は inherent method とする。演算 trait は線型結合を返す規則を表し、式の `act`・`product`・`adjoint` メソッドは遅延構築を担当する。

以下の `P` は、状態式では `L::State`、演算子式では `L::Operator` を表す。型制約は第4.3節に従い、必要な操作ごとに付ける。

| 対象 | API・引数・返却値 |
| --- | --- |
| 両式型の primitive 取り込み | `primitive(p: P) -> Self` |
| 両式型の線型結合の取り込み | `from_linear_combination(terms: LinearCombination<Self, T>) -> Self`。要素の式への変換は呼び出し前に `map` で行う |
| 両式型の零・和 | `zero() -> Self`、`sum(iter: impl IntoIterator<Item = Self>) -> Self` |
| 両式型の係数乗算 | `scaled(&self, coefficient: Complex<T>) -> Self` |
| 演算子の恒等元・積・随伴 | `identity() -> Self`、`product(&self, rhs: &Self) -> Self`、`product_all(iter: impl IntoIterator<Item = Self>) -> Self`、`adjoint(&self) -> Self` |
| 状態式への遅延作用 | `act(&self, operator: &LocalOperator<L, T>) -> Self` |
| 状態式の評価 | `evaluate(&self) -> Self` |
| 状態式の内積 | `inner_product(&self, rhs: &Self) -> Complex<T>` |
| 演算子式の局所 trace | `trace(&self) -> Complex<T>` |
| 演算子式の導出操作 | `commutator(&self, rhs: &Self) -> Self`、`anticommutator(&self, rhs: &Self) -> Self`、`hilbert_schmidt_inner_product(&self, rhs: &Self) -> Complex<T>` |
| 状態式の局所期待値 | `expectation(&self, operator: &LocalOperator<L, T>) -> Complex<T>` |

`+`・`-`・単項 `-`、演算子積、演算子と状態の積、複素係数による左右の乗算を提供する。式同士の二項演算は所有値・所有値、所有値・借用値、借用値・所有値、借用値・借用値の4通りに対応する。内部は所有値を受け取る構築処理へ集約し、所有ハンドルを移動、借用ハンドルは共有参照だけを複製する。名前付きメソッドは共有借用を基本とする。primitive・係数の `Clone` は構築の条件にしない。実数係数はまず `Complex<T>` への明示変換を基本とし、暗黙の精度変換は導入しない。

### 3.2. constructor の意味

内部表現によらず、次の意味を満たす。

| 入力 | 意味 |
| --- | --- |
| primitive | 係数 `1` の生成を必要としない `Primitive` 節点 |
| 線型結合 | 各項の係数と式を `Scale` にし、項順を保って `Sum` にする。空列は zero |
| 空の和 | 空の `Sum` による zero |
| 空の演算子積 | `Identity` |
| 一要素の和・積 | その要素を返す |
| 二つの式の加算 | operand のハンドルを保持する `Sum` |
| 係数乗算 | 係数と operand を保持する `Scale` |
| 二つの演算子の積 | 順序付きの `Product` |
| 複数要素の和・積 | 入力順を保つ和・積。積は二項の `Product` を組み合わせる |

構築時に既存の和を平坦化したり、係数の乗算・比較・項の結合を行ったりしない。係数 `0`・`1` もそのまま保持する。既存の大きな式はハンドルを共有し、primitive・係数の `Clone` を構築の一律の条件にしない。

## 4. 評価

### 4.1. 遅延構築と作用の評価

`product`・`act`・`adjoint` は遅延式を構築し、`inner_product`・`trace` は必要な処理を実行して `Complex<T>` を返す。`LocalState::act(&self, operator)` は `Act(ActExpr { state, operator })` 節点を含む状態式を構築し、その式に対する `LocalState::evaluate(&self) -> Self` が作用を評価する。

以下は利用形の例である。同じ API の実行例を spin-half の doc test でも検証している。

```rust
let expr = state.act(&operator); // 遅延した作用を構築する
let result = expr.evaluate();    // 構築済みの状態式を評価する

let result = state.act(&b).act(&a).evaluate(); // B、A の順に作用させ、A(Bψ) を得る

let expr = &operator * &state; // 演算子構文も同じ遅延構築を行う
let result = expr.evaluate();
```

`act` と `evaluate` はともに共有借用を受け取り、元の式を変更せず通常の `LocalState<L, T>` を返す。加算・作用・内積を続けて適用でき、`act` のメソッドチェーンを構築するだけでは局所規則を呼ばない。局所規則の `LinearCombination` は式木へ取り込み、`evaluate` 自体の返り値は状態式とする。汎用の `simplify` や局所積だけを整理する公開操作は追加しない。

`evaluate()` は状態式の和・係数乗算・入れ子の作用をたどり、必要な primitive の `Act::act` に委譲し、作用後の状態式を得る。例えば `Act(Sx, Up)` は `Down` の係数 `1/2` の式になる。基底を入力に取らず、選択基底への展開・座標化や密・疎行列の生成も行わない。入力を変更せず、すでに作用が評価されている部分は共有して再利用する。内積・期待値は、利用者が事前に `evaluate()` を呼ぶことを必須とせず、必要な部分を内部で評価する。

### 4.2. 評価経路と結果の項

状態の評価、演算子の作用、演算子項の走査を目的別に分ける。線型構文の分解と係数処理は内部で共用し、すべての評価を演算子の完全展開へ集約しない。走査には明示的な作業スタックを使い、ユーザ式の深さに比例した関数再帰に依存しない。

`(AB)ψ` は `B`、`A` の順に作用させ、不要な primitive 間の `Product` を呼ばずに処理する。一方、`trace(AB)` などで必要な局所積は `Product` に委譲する。随伴は構築時に分配せず、評価時に係数の共役と積順序の反転を処理する。

`evaluate()` の結果では、同じ primitive の項を自動結合せず、零係数の項も自動削除しない。明示した空の `Sum` は zero として処理する。和は左から右、演算子積の作用は右から左に処理し、局所規則が返す項の順序を保存する。評価に primitive の比較を必須とせず、数学的に零を表す式と明示した zero を区別する。同類項の結合は Step 5 の基底展開・座標取得で行う。

評価中の一時表現には共通パーツを必要に応じて利用し、既存の葉や恒等作用など入力をそのまま利用できる経路では共有実体を再利用する。局所規則が返す線型結合は `map` と共通 constructor で式木へ取り込む。係数を変更するために primitive 本体を複製せず、数値係数の複製は各演算で必要な箇所に限定する。特定の一時表現をすべての評価経路に共通する完全展開形として強制しない。

### 4.3. 数値制約とエラー

式の保持・与えられた係数や operand を保存する遅延構築・`map` には数値制約を付けない。係数計算・共役・複製に必要な数値制約は、それを実行するメソッドにだけ付ける。係数 `1` を生成する線型結合の constructor や、係数 `−1` を生成する式の単項負号などには、その生成に必要な最小限の制約を付ける。既存依存の数値演算が要求する制約も必要な箇所に限定し、全操作へ一律に `Clone + Num` を要求する方針にはしない。primitive の規則は現在の演算ごとの最小限の制約を維持する。全能力をまとめた係数 trait は導入せず、`Float`・`Copy` も一律に要求しない。

状態への作用評価は primitive の `Act` と必要な随伴規則などで処理し、演算子の `Product` 能力を要求しない。局所積が必要な trace などのメソッド側で `Product` を要求する。`Identity<T>` も恒等演算子の規則を実際に取得する操作だけで要求し、`Identity` 節点の構築自体には要求しない。

局所構造と精度は型で一致し、primitive の規則も infallible なので、基本的な局所評価には包括的な `Result` を付けない。NaN・infinity を独自に拒否せず、内部の項除去にも近似許容誤差を使わない。

実数値を返す norm 等のエラー設計は、その機能を追加する段階で決める。

## 5. 完了条件と検証

spin と独立した三次元局所代数の双方で、それぞれの代数だけを使う別の系を構成し、積・作用・内積・随伴・trace を小規模行列と照合する。両代数を一つの系に混在させない。単精度・倍精度、虚部のある係数、線型従属な primitive、zero・identity・空和・空積を含める。

規則の呼び出し回数と primitive の複製回数を観測し、構築時の未評価、式の共有、入力の不変性、`(AB)ψ` で不要な `Product` が呼ばれないことを確認する。精度・局所構造の混在は compile-fail テストで検証する。

局所規則・恒等演算子・基底展開が対応する要素の有限線型結合を返し、式木への取り込みで `Primitive`・`Scale`・`Sum` に構成されることを確認する。基底展開に現れる primitive が選択基底の要素であること、基底順序番号・遅延列挙と整合することを検証する。`map` が係数・項数・順序・重複・零係数を保持し、要素・係数の `Clone` を要求しないことも確認する。空・単一要素のコンテナ割り当て特性を維持する。これらの検証のために式の構造比較を公開することは要求しない。

core・spin-half・テスト・doc test の演算 trait を関連型 `Term` と `Act::Operator` に揃え、局所規則・恒等演算子・既存の局所基底の線型結合を共通 constructor で式木へ取り込めることを Step 2 の完了条件に含める。現行の数値規約・演算ごとの最小限の型制約・基底順序・非正規直交基底の数学的な結果は維持する。

`LocalStructure` と `Identity<T>` の分離、`Term<P, T>` の型別名化と `operand` への変更、`StateBasis`・`OperatorBasis` の削除と `Basis` の直接利用、次元取得 API の不要な `T` の除去も Step 2 の移行対象とする。core・spin-half・テスト・doc test のコードから両 marker trait の定義・再公開・実装宣言・利用を除去する。数値能力を持たない係数型でも式を保持・構築できることと、`Product` を実装しない演算子でも状態への作用評価が成立することを型検証する。演算子構文の所有値・借用値の4通りは数学的な結果と共有の両方を検証する。

`&a + &a` と `a.scaled(c)` が元の部分式を共有すること、式の `Clone` が primitive・`L`・`T` の `Clone` を要求しないことを確認する。別々に構築した同形の式を自動的に共有することは完成条件に含めない。

`ψ.act(&b).act(&a)` と `ψ.act(&a.product(&b))` の数学的な一致を確認し、`act` のチェーンが遅延構築、`evaluate` が評価であることと、左作用の順序を検証する。

評価結果を再合成・再評価でき、入力の式を保持した利用者が引き続きそれを計算に利用できることを確認する。基底を与えずに評価できること、constructor の空・一要素・複数要素の意味、係数 `0`・`1` の保持、評価結果の重複項・零係数の保持、評価済み部分の再利用も検証する。source の自動保持・構造比較・式表示は受入条件に含めない。

長い積・作用のチェーンを評価・破棄し、式の深さに比例した関数再帰を避けられることを検証する。破棄時には共有先が残る部分式を引き続き利用できることも確認する。Step 2 の性能上の合格条件は、部分式の共有、不要な primitive 複製の抑制、空・単一項のコンテナ無確保とする。`map` 全般の無確保や特定の実行時間は保証せず、速度比較・追加の割り当て最適化は後続で検証する。

既存の direnv 環境で対象テスト、workspace の test・format check・Clippy・警告をエラーとするドキュメント生成を実行し、成功を確認した。詳細は第1.1節に記録した。

共通パーツを使う再帰型、数値・`Clone` 制約を追加しない構築と共有、借用を含む primitive の保持、関連型・`map`・公開 API の移行を実際のコードとテストで確認した。性能の合格条件は上記の共有・複製・割り当て特性とし、速度比較・追加の最適化は未実施である。

## 6. 確定事項と後続の設計

共通パーツの役割と公開範囲、局所式の節点構成、薄い `LinearCombination` と項型の統合、演算 trait の `T` と関連型 `Term`、`map` と共通 constructor による式木への取り込み、および `Arc` の共有目的は確定事項とする。局所構造と演算能力の分離、基底の marker trait の削除と `Basis` の直接利用、メソッドごとの型制約、所有値・借用値の構築経路、反復的な評価・破棄と性能上の合格条件も採用済みとし、core・spin-half の API 整備は Step 2 で完了した。

非公開の作業スタック・補助関数・module と各操作の型制約は実装に反映した。今後の最適化でも、この方針と検証条件を維持する。

後続段階では次を扱う。

- Step 3 で同じ局所構造型 `L` の局所式を配置・合成する具体的な型と API。
- Step 5 の式全体の基底展開・複合基底・座標取得。局所基底の線型結合を局所式へ取り込む接続は Step 2 に含める。
- norm・正規化・bra の詳細。

非公開の module 分割などは、確定方針とテスト条件の範囲内で実装時に判断する。異種局所構造の混在・合成と型消去、source の自動保持・式の構造比較・表示は初期版の対象外とする。
