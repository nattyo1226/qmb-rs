# qmb-algebra 実装検討メモ

更新日: 2026-10-21

## 1. 位置付けと提供範囲

本書は、[実装計画](implementation-plan.md)の Step 2「局所式と評価」に向けた一時的な作業用記録である。恒久的な仕様書として維持することは想定せず、既存の恒久ドキュメントからはリンクしない。ここに記す API や機能が実装済みであることは意味しない。

局所式の数学的契約は、[局所 primitive・式木の数学的能力](design/primitive-expression.md)と[要件](requirements.md)に従う。Step 1 の `qmb-core`・`qmb-spin-half` に対して第1.1節の作用 API の変更を先に完了し、その公開契約を使って共通式木と評価を `qmb-algebra` に実装する。

本書の設計方針・型名・操作名は確定事項とする。実装時に具体化する内部の詳細と、後続段階の設計事項は第7節にまとめる。

局所状態・演算子の構築、加算、係数乗算、zero・identity、演算子積、作用、内積、随伴、局所 trace を実装する。これらから作れる commutator・anticommutator、局所期待値、Hilbert–Schmidt 内積も含める。

norm・正規化・bra は、実数値の扱いと API を別途検討してから追加する。Index・空間・配置は Step 3、基底展開・座標取得は Step 5 に置く。これは初期版全体からそれらの機能を除外するものではない。

`qmb-algebra` も `no_std + alloc` とし、通常依存はまず `qmb-core` のみで開始する。`Arc`・`Vec` など既存の標準機構で実装し、新しい production 依存や feature は追加しない。

内部を局所式、構築構文、評価、比較・表示に分け、主要な公開型を crate root から再公開する。`qmb-spin-half` は開発用依存に置き、通常依存へ具体的代数を持ち込まない。

### 1.1. 実装前の前提作業：状態 primitive の Act の整備

`qmb-algebra` の実装に着手する前に、`qmb-core` の作用 API を以下の状態 primitive 側の `Act<O, T>` の契約に揃える。この整備は今後行う前提作業であり、本書の記録時点で完了したことを意味しない。

`Act` は「この状態 primitive に、指定した演算子 primitive を作用させられる」能力を表す。

```rust
pub trait Act<O, T>: Sized {
    /// 数学的には operator * self を計算する。
    fn act(&self, operator: &O) -> LinearCombination<Self, T>;
}
```

`LocalStructure<T>` の関連型の制約は次の形に変更する。`DIMENSION` と `identity()` の契約は維持する。

```rust
type State: InnerProduct<T> + Act<Self::Operator, T>;
type Operator: Product<T> + Adjoint<T> + Trace<T>;
```

前提作業には、以下を含める。

1. core の trait 定義・再公開・関連型の制約を変更する。
2. spin の作用規則を `Act<SpinHalfOperator, T> for SpinHalfState` として実装し、局所規則の数値制約と数学的結果を維持する。
3. 三次元局所代数などのテスト実装、呼び出し箇所、doc test、作用の実装先を説明する既存ドキュメントを新しい契約に揃える。
4. 局所規則と型制約のテスト、workspace の test・format check・Clippy・ドキュメント生成を既存の direnv 環境で確認してから、`qmb-algebra` の実装を開始する。

Act は状態をメソッドの受け手とし、数学的には演算子の左作用を計算する。primitive の `Act::act` は局所規則を計算して有限線型結合を返す。一方、状態式の `act` は遅延した作用を構築し、`evaluate` が必要な primitive の `Act::act` を呼ぶ。両者の役割を区別する。

## 2. 式の型・内部表現・共有

### 2.1. 公開する式型

`LocalState<L, T>` と `LocalOperator<L, T>` の二型を採用する。`L: LocalStructure<T>` から primitive 型を取得し、局所構造と係数精度の整合性を型で保証する。

将来の複合状態・演算子とは公開型を分ける。共通式木を使う契約は、局所式と複合式の公開型を一つに統一することまでは要求しない。

### 2.2. Arc による所有権と不変共有

式節点・primitive・source は `Arc` を基本として、不変な実体を共有する。式の `Clone` は共有参照の複製だけにし、primitive や `T` の `Clone` を要求しない。評価中は借用を優先し、不要な `Arc::clone` を避ける。

primitive に一律の `Clone`・`Copy`・`Send`・`Sync`・`'static` を要求しない。借用を含む primitive も、その型の lifetime の範囲で利用できるようにする。係数は節点内に保持し、数値演算の制約は必要な操作に付ける。

`Send + Sync` は保持する型に応じて成立させる。非 `Send` な primitive も格納でき、その場合は式のスレッド間共有ができない。将来のスレッド間共有への接続を重視して `Arc` を採用するが、並列 evaluator の実装は後続へ延期する。

共有の内部表現は非公開とし、`Rc`／`Arc` 切替 feature は設けない。`Arc` は参照の複製・破棄に atomic 操作を用いるため、`Rc` より参照カウントのコストがある。借用による読み取りのたびに参照カウントを操作する必要はない。

`alloc::sync::Arc` は `no_std + alloc` で利用できるが、ポインタ幅の atomic 操作を持つ target が必要となる。これらの制約は既存環境の Rust 1.98.0 の `alloc` 実装で確認した。性能比較は未実施である。

### 2.3. 式木の節点

| 対象 | 節点 |
| --- | --- |
| 状態式 | `Primitive / Zero / Add / Scale / Apply` |
| 演算子式 | `Primitive / Zero / Identity / Add / Scale / Product / Adjoint` |

内積・trace は scalar を求める操作とし、今回は scalar 用の式木を導入しない。

`Add`・`Product` は二分節点とし、既存の長い項列を構築のたびにコピーせず、括り方も保持する。任意個の和・積は iterator を受ける constructor で提供し、空和は zero、演算子の空積は identity とする。

### 2.4. 複合式への接続

初期版では、一つの系に含まれる全基本因子を同じ局所構造型 `L` に揃える。Step 3 でも型付き局所式を共有して保持し、同じ `L` の局所式と Index を組にして複合式の葉とする。複合式の葉から必要な操作を局所 evaluator へ静的に委譲する。

異種局所構造の混在・合成と、そのための型消去 API・保持機構は後続へ延期し、Step 2・3 の設計対象に含めない。独自代数は `L` の別実装として別の系で単独利用できる。同種因子の入れ子や一つの Site への複数自由度の配置は維持する。

将来の異種構造への対応や、型消去後の `Send + Sync` の扱いは、その機能を導入する段階で設計する。初期版のために `'static` 等の追加制約を先に課さない。

## 3. 公開 API と構築規約

### 3.1. メソッドと演算子構文

名前付き constructor・メソッドを基本 API とし、演算子構文は同じ処理への入口とする。primitive・`LinearCombination` の取り込みは局所規則を実行せず、通常の式木を構築する。

primitive・有限線型結合の取り込みは所有値を受け取り、式同士の操作は `&self` と共有借用を基本とする。変更されない operand は `Arc` で共有する。名前付き API は inherent method とし、primitive 用で返却形式が異なる core の `Product`・`Act` 等を局所式へそのまま実装しない。

以下の `P` は、状態式では `L::State`、演算子式では `L::Operator` を表す。型制約は第4.3節に従い、必要な操作ごとに付ける。

| 対象 | API・引数・返却値 |
| --- | --- |
| 両式型の primitive 取り込み | `primitive(p: P) -> Self` |
| 両式型の有限線型結合の取り込み | `from_linear_combination(terms: LinearCombination<P, T>) -> Self` |
| 両式型の零・和 | `zero() -> Self`、`sum(iter: impl IntoIterator<Item = Self>) -> Self` |
| 両式型の係数乗算 | `scaled(&self, coefficient: Complex<T>) -> Self` |
| 演算子の恒等元・積・随伴 | `identity() -> Self`、`product(&self, rhs: &Self) -> Self`、`product_all(iter: impl IntoIterator<Item = Self>) -> Self`、`adjoint(&self) -> Self` |
| 状態式への遅延作用 | `act(&self, operator: &LocalOperator<L, T>) -> Self` |
| 状態式の評価 | `evaluate(&self) -> Self` |
| 状態式の内積 | `inner_product(&self, rhs: &Self) -> Complex<T>` |
| 演算子式の局所 trace | `trace(&self) -> Complex<T>` |
| 演算子式の導出操作 | `commutator(&self, rhs: &Self) -> Self`、`anticommutator(&self, rhs: &Self) -> Self`、`hilbert_schmidt_inner_product(&self, rhs: &Self) -> Complex<T>` |
| 状態式の局所期待値 | `expectation(&self, operator: &LocalOperator<L, T>) -> Complex<T>` |

`+`・`-`・単項 `-`、演算子積、演算子と状態の積、複素係数による左右の乗算を提供する。式同士は所有値・借用値の組み合わせに対応させる。実数係数はまず `Complex<T>` への明示変換を基本とし、暗黙の精度変換は導入しない。

### 3.2. constructor が作る構造

次の規約に統一する。

| 入力 | 構築する式 |
| --- | --- |
| 空の和 | `Zero` |
| 空の演算子積 | `Identity` |
| 一要素の和・積 | その要素をそのまま返す |
| 複数要素の和・積 | 入力順の左結合で構築する |
| 有限線型結合の各項 | `Scale(coefficient, Primitive(p))` |
| 空の有限線型結合 | `Zero` |

有限線型結合の項列も入力順の左結合で和にする。係数が `1` や `0` でも `Scale` を残し、構築に係数比較・数値演算を要求しない。直接 `primitive(p)` で作った式には `Scale(1, …)` を追加しない。両者は数学的には同じでも構造比較では異なる式となる。

## 4. 評価

### 4.1. 遅延構築と作用の評価

`product`・`act`・`adjoint` は遅延式を構築し、`inner_product`・`trace` は必要な処理を実行して `Complex<T>` を返す。`LocalState::act(&self, operator)` は `Apply(operator, state)` 節点を含む状態式を構築し、その式に対する `LocalState::evaluate(&self) -> Self` が作用を評価する。`Apply` は内部節点名であり、公開メソッドに `apply` は設けない。

以下は利用形の例であり、実装・検証済みのコードではない。

```rust
let expr = state.act(&operator); // 遅延した作用を構築する
let result = expr.evaluate();    // 構築済みの状態式を評価する

let result = state.act(&b).act(&a).evaluate(); // B、A の順に作用させ、A(Bψ) を得る

let expr = &operator * &state; // 演算子構文も同じ遅延構築を行う
let result = expr.evaluate();
```

`act` と `evaluate` はともに共有借用を受け取り、元の式を変更せず通常の `LocalState<L, T>` を返す。加算・作用・内積を続けて適用でき、`act` のメソッドチェーンを構築するだけでは局所規則を呼ばない。内部の有限線型結合を専用の公開結果型にはしない。汎用の `simplify` や局所積だけを整理する公開操作は追加しない。

`evaluate()` は状態式の和・係数乗算・入れ子の作用をたどり、必要な primitive の `Act::act` に委譲し、作用後の状態式を得る。例えば `Apply(Sx, Up)` は `Scale(1/2, Down)` になる。入力と source は変更せず、すでに作用が評価されている部分は共有して再利用する。内積・期待値は、利用者が事前に `evaluate()` を呼ぶことを必須とせず、必要な部分を内部で評価する。

### 4.2. 評価経路と結果の項

状態項の走査、演算子の作用、演算子項の走査を目的別に分ける。線型構文の分解と係数処理は内部で共用し、すべての評価を演算子の完全展開へ集約しない。

`(AB)ψ` は `B`、`A` の順に作用させ、不要な primitive 間の `Product` を呼ばずに処理する。一方、`trace(AB)` などで必要な局所積は `Product` に委譲する。随伴は構築時に分配せず、評価時に係数の共役と積順序の反転を処理する。

`evaluate()` の結果では、同じ primitive の項を自動結合せず、零係数の項も自動削除しない。明示的な `Zero` と、局所規則が返す空の有限線型結合は zero として処理する。和は左から右、演算子積の作用は右から左に処理し、局所規則が返す項の順序を保存する。評価に primitive の比較を必須とせず、数学的な零と構文上の `Zero` を区別する。同類項の結合は Step 5 の基底展開・座標取得で行う。

一時的な項をまとめる必要がある箇所では、非公開の内部表現として `LinearCombination<Arc<P>, T>` を用いる。既存の葉は `Arc` を共有し、局所規則が返す所有値の primitive は新しい `Arc` に移す。係数を変更しても primitive 本体を複製せず、恒等作用など入力をそのまま利用できる経路では既存の共有実体を再利用する。数値係数の複製は各演算で必要な箇所に限定する。この内部形式を、すべての評価経路に共通する完全展開形として強制しない。

### 4.3. 数値制約とエラー

型の保持と構築には最小限の制約を付ける。評価では `T: Clone + Num`、共役を使う操作には追加で `Neg<Output = T>` を要求する方針とし、実際の制約は各操作の必要条件を確認して付ける。`Float`・`Copy` は一律に要求しない。

局所構造と精度は型で一致し、primitive の規則も infallible なので、基本的な局所評価には包括的な `Result` を付けない。NaN・infinity を独自に拒否せず、内部の項除去にも近似許容誤差を使わない。

実数値を返す norm 等のエラー設計は、その機能を追加する段階で決める。

## 5. source・構造比較・表示

### 5.1. 現在の式と source の引き継ぎ

各式に「現在の式」と「元の source」の二つの共有ルートを保持する。現在の式は以後の計算に使い、source は評価による置き換えを行う前の式を表す。通常の構築では両者は同じ実体を共有し、`evaluate()` は現在の式だけを結果へ置き換えて source を保持する。

評価結果をさらに合成するときは、現在の式同士、source 同士からそれぞれのルートを構築する。別の操作履歴リストは持たず、式木の共有だけで実現する。次の `current`・`source` は説明用の名前であり、公開フィールドを意味しない。

| 操作 | 現在の式 | source |
| --- | --- | --- |
| `φ = ψ.act(&A).evaluate()` | 評価した作用結果 | `Apply(A.source, ψ.source)` |
| `χ = φ + η` | `Add(φ.current, η.current)` | `Add(φ.source, η.source)` |
| `ξ = χ.act(&B).evaluate()` | 評価した作用結果 | `Apply(B.source, χ.source)` |

計算は現在の式から行い、保持した source を毎回評価し直さない。入力は変更せず、source の保持は評価結果のキャッシュとは分ける。

### 5.2. 参照・比較・表示の API

| API | 意味 |
| --- | --- |
| `source(&self) -> Self` | 保持した source を通常の式として取得する |
| `source_eq(&self, rhs: &Self) -> bool` | 両者の source の構造を比較する |
| `structural_eq(&self, rhs: &Self) -> bool` | 両者の現在の式の構造を比較する |

`source()` は共有ルートを複製し、返す式の現在の式と source をともに元の source にする。繰り返し呼んでも意味は変わらず、評価は実行しない。比較には対象となる primitive と係数の `PartialEq` を必要な範囲で要求し、`L: PartialEq` は要求しない。Step 2 では式型への `PartialEq` 実装を設けず、比較の意味をメソッド名で明示する。

構造比較では構文・係数・primitive を比較し、構築時の平坦化、項順序変更、zero 除去は行わない。NaN を含む比較を壊さないよう、共有ポインタの一致だけで真を返さない。

表示はまず `Debug` を提供し、現在の式を評価せず表示する。元の source は `expr.source()` の `Debug` 表示で確認する。表示対象の primitive と係数に必要な `Debug` 制約を付け、`L: Debug` は要求しない。`Display` は後続へ回す。

## 6. 完了条件と検証

spin と独立した三次元局所代数の双方で、それぞれの代数だけを使う別の系を構成し、積・作用・内積・随伴・trace を小規模行列と照合する。両代数を一つの系に混在させない。単精度・倍精度、虚部のある係数、線型従属な primitive、zero・identity・空和・空積を含める。

規則の呼び出し回数と primitive の複製回数を観測し、構築時の未評価、式の共有、source の保持、`(AB)ψ` で不要な `Product` が呼ばれないことを確認する。精度・局所構造の混在は compile-fail テストで検証する。

`ψ.act(&b).act(&a)` と `ψ.act(&a.product(&b))` の数学的な一致を確認し、`act` のチェーンが遅延構築、`evaluate` が評価であることと、左作用の順序を検証する。

評価結果の再合成と再評価でも、現在の式と source が第5.1節のとおりに引き継がれること、`source()` の反復、二種類の構造比較、非評価の `Debug` 表示を確認する。constructor の空・一要素・複数要素の構造、係数 `0`・`1` の保持、評価結果の重複項・零係数の保持、評価済み部分の再利用も検証する。

既存の direnv 環境で対象テスト、workspace の test・format check・Clippy・ドキュメント生成を確認する。これは今後の検証方針であり、Step 2 のチェックを実行済みとするものではない。

## 7. 実装時の具体化と後続の設計事項

本書の設計方針・API・構築規約は再選択の対象とせず、以下の未記述の詳細を具体化する。

- 数値演算・比較・表示ごとの正確な trait 制約と、演算子構文の所有値・借用値の実装。
- 二つの共有ルートを保持する非公開フィールド、内部走査・一時領域の具体的な実装、非公開の補助型・関数名。
- Step 3 で同じ局所構造型 `L` の局所式を配置・合成する具体的な型と API。

非公開の module 分割や走査処理の詳細は、確定方針とテスト条件の範囲内で実装時に判断する。異種局所構造の混在・合成と型消去は初期版の対象外とする。norm・正規化・bra の詳細は、それらを追加する段階で設計する。
