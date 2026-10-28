# scaled の変換対応と try_sum の実装計画

作成日・更新日: 2026-10-28
状態: 初回範囲を2026-10-28に実装。scaled の impl Into<C> 対応、4式型の try_sum、State::try_tensor、Operator::try_product_all、および QuantumSystem の構築・集約入口を追加した。検証結果は末尾に記録する。

後続の合意として、利用者向けの構築・集約の入口を QuantumSystem に揃える。メソッド候補と優先順位は[QuantumSystem の構築 API](quantum-system-construction.md)に記録する。本書は scaled と try_sum の基礎処理・検証を扱い、関連関数の転送入口も合わせて実装する。初回の範囲は scaled の変更と構築 API 計画の優先順位1〜3までとすることで2026-10-28に合意した。

## 1. 目的と対象

[QuantumSystem の導入](api-type-parameters.md)により、primitive と係数型を一度選んで式を構築できるようになった。[two_spins の利用例](../../examples/src/bin/two_spins.rs)では、さらに次の記述を減らしたい。

- 実数のスカラー倍にも `Complex::new(value, 0.0)` を書く必要がある。
- 配置・積が `Result` を返すため、項を生成して和にする際に、利用側で `Vec` や `collect::<Result<...>>()?` を挟む必要がある。

既存の式型・評価器・精度選択・空間の契約を維持し、共通の構築 API を改善する。配置込みの積、構築済みの式の tensor・積と try 版も初回の範囲に含め、その契約と命名の比較・推奨案は別の構築 API 計画に記録する。零・単位元、単独配置、QuantumSystem の局所式集約用の転送メソッドは後続候補とする。個別 primitive 用メソッド、DSL、構築用 builder は対象に含めない。macro と Mul の入力拡張は保留とする。`scaled_real` や `QuantumSystem::coefficient` も追加せず、変換を `scaled` の入口へ集約する。

## 2. scaled の引数を impl Into<C> にする

### 2.1. 対象とシグネチャ

次の全公開 `scaled` を対象とする。局所式・複合式だけでなく、同名 API の一貫性のため `LinearCombination` も揃える。

| 型 | 実装したシグネチャ | 係数型に必要な能力 |
| --- | --- | --- |
| `LocalState<S,O,C>` | `scaled(&self, coefficient: impl Into<C>) -> Self` | 変換以外は不要 |
| `LocalOperator<O,C>` | 同上 | 同上 |
| `State<I,S,O,C>` | 同上 | 同上 |
| `Operator<I,O,C>` | 同上 | 同上 |
| `LinearCombination<P,C>` | `scaled(self, coefficient: impl Into<C>) -> Self` | 現行と同じ `C: Clone + Mul<Output=C>` |

入力を入口で一度だけ `.into()` し、以後は既存の C を受け取る処理へ渡す。式の内部 `into_scaled`・節点・保持する係数型は変更しない。`LinearCombination` は変換後の C を使って従来どおり各項の係数を乗算する。元の入力型に Clone を要求しない。

```rust
// 実装した形。式ハンドル4型に同じ方針を適用する。
pub fn scaled(&self, coefficient: impl Into<C>) -> Self {
    self.clone().into_scaled(coefficient.into())
}
```

既存の C を直接渡す経路は恒等変換によって受け入れる。`From<R> for C` がある場合には、その blanket な `Into<C>` 実装を利用できる。利用者が直接実装した `Into<C>` も受け入れる。新たなキャスト・変換 trait・production dependency は導入しない。

### 2.2. 数値精度と型推論

変換先は式の C であり、入力によって式の精度を変更しない。標準 `Complex<f32>` に対する f32、`Complex<f64>` に対する f64 と、同じ C の入力を検証する。接尾辞なしの `0.5` も両精度で確認する。一方、明示的な f64 を `Complex<f32>` に渡す場合や、単精度・倍精度の複素係数を混ぜる場合は、対応する既存の Into がない限りコンパイル時に拒否する。整数を含む任意の数値から変換できるとは約束しない。

`scaled(C)` によって C 自体を推論していたコードは、引数が `impl Into<C>` になると推論できなくなる可能性がある。これは単に受け入れる値が増えるだけの互換変更とは扱わない。

2026-10-28 に、利用者と型推論の変化を許容して `scaled(impl Into<C>)` で進める方針を合意した。従来の型推論との完全互換は実装の条件にせず、必要な構築入口・入力値の型注釈や、余分な `.into()` の除去で対応する。QuantumSystem で primitive と係数型を先に選ぶ使い方を基本とする。既存の C を直接渡す経路、選択した精度、数学的な結果、最小能力の契約は維持する。

以下は移行箇所を把握して調整するために確認し、推論差があること自体を理由に API の変更を差し戻す条件にはしない。

- QuantumSystem から C が確定している状態・演算子の、実数リテラル・複素係数によるスカラー倍。
- 直接の constructor で C を省略し、scaled の引数や後段の評価・戻り値注釈から推論する既存コード。
- `Complex::new(...)` の成分型、`.into()`・`Default::default()` などの変換先を、旧 scaled の引数型に依存していたコード。
- generic な C、メソッドの関数値としての利用、および core・spin-half・normalize・export の既存呼び出し。

注釈が必要になった場合は、構築入口の QuantumSystem 選択や既存式型の注釈で補い、移行方法を記録する。推論のために既定精度や新たな数値能力の bound を導入しない。構築時の最小能力と no_std + alloc を維持する。

`*` のスカラー演算、`scalar`・`term` など別の constructor の引数は今回広げない。従来どおり `scaled` を任意の係数入力に対する基本の入口とする。

### 2.3. 型推論の変化と移行例

以下は scaled の引数変更に対する移行案内である。`Spin = SpinHalf<Complex<f64>>`、Op は SpinHalfOperator とする。

1. scaled の引数から式の C を決めていたコードは、構築入口で C を指定する。

   ```rust
   // 移行前は C = Complex<f64> と推論できるが、変更後は推論できない。
   let op = LocalOperator::primitive(Op::Sz)
       .scaled(Complex::new(0.5_f64, 0.0));

   // 移行後の案。QuantumSystem で先に C を選択する。
   let op = Spin::operator(Op::Sz).scaled(0.5);
   // 直接 constructor を使う場合は型引数や戻り値注釈でも指定できる。
   let op = LocalOperator::<Op, Complex<f64>>::primitive(Op::Sz)
       .scaled(Complex::new(0.5, 0.0));
   ```

2. 引数側の `.into()` は変換先の期待型を失う場合がある。既に C が選択済みでも発生するため、不要な事前変換を外して入力を直接渡す。

   ```rust
   // 移行前は成功。変更後は内側の into() の変換先が曖昧になる。
   Spin::operator(Op::Sz).scaled(0.5_f64.into());
   // 移行後の案。
   Spin::operator(Op::Sz).scaled(0.5_f64);
   ```

3. `Default::default()` のように引数の期待型で生成型を決めていた場合は、その入力型を明示する。

   ```rust
   // 移行前は C の default。変更後は入力型を推論できない。
   Spin::operator(Op::Sz).scaled(Default::default());
   // 移行後の案。
   Spin::operator(Op::Sz).scaled(Complex::<f64>::default());
   ```

これらは Rust 1.98.0 で、製品コードを変更せず、実際の LocalOperator に `scaled_into(&self, value: impl Into<C>)` という補助メソッドを実装して確認した。上記3件は移行前の scaled ではコンパイルでき、補助メソッドでは推論エラーになり、記載した移行方法ではコンパイルできた。

また、SpinHalf<Complex<f32>>／SpinHalf<Complex<f64>> で C を選択済みなら、`.scaled_into(0.5)` と `.scaled_into(Complex::new(0.5, 0.0))` は両精度とも追加の型注釈なしでコンパイルできた。この確認は引数型の変更を模した小規模な検証であり、製品実装や workspace 全体の回帰検証を完了したことは意味しない。

## 3. try_sum を追加する

2026-10-28 に、try 版のエラー契約を「ExpressionError 固定・最初の入力 Err で停止・全入力成功後に合成検証」とする方針を利用者と合意した。本節は未決の比較案ではなく、実装・回帰検証の基準とする。製品 API へ実装済みである。

### 3.1. 対象とシグネチャ

既存の `sum` は維持し、式ハンドル4型に関連関数を追加する。入力は所有式の `Result` 列、エラー型は既存 constructor と同じ `ExpressionError` に固定する。入力・合成のエラーを新たなラッパーで包まず、そのまま返す。任意のエラー E や `E: From<ExpressionError>` への一般化は今回行わない。独自エラーを入力へ渡したい場合は、ExpressionError として表現できるものを呼び出し側で map_err するか、独自エラーを呼び出し側で処理する。

try という接頭辞は入力要素が Result であることを表す。通常の sum も複合式では合成検証により Result を返すため、try は「このメソッドだけが失敗し得る」という意味ではない。

```rust
// 実装済みの API。各型の inherent impl に置く。
// LocalState / LocalOperator: primitive・係数の能力は要求しない。
pub fn try_sum(
    operands: impl IntoIterator<Item = Result<Self, ExpressionError>>,
) -> Result<Self, ExpressionError>;

// State<I,S,O,C>: 現行 sum と同じ I: Eq を要求する。
pub fn try_sum(
    space: &HilbertSpace<I>,
    states: impl IntoIterator<Item = Result<Self, ExpressionError>>,
) -> Result<Self, ExpressionError>;

// Operator<I,O,C>: 現行 sum と同じ I: Eq を要求する。
pub fn try_sum(
    operators: impl IntoIterator<Item = Result<Self, ExpressionError>>,
) -> Result<Self, ExpressionError>;
```

`State` の明示空間引数は省略しない。特に空の入力列でも、どの空間に属する zero state かを保つ。`OperatorOf<Q,I>` などは既存型への別名なので、別名自体に独立した実装は追加しない。利用者向けには QuantumSystem の `state_sum`・`try_state_sum`・`operator_sum`・`try_operator_sum` を既定メソッドとして用意し、既存の sum と本節の try_sum に委譲する。局所式用の転送入口は構築 API 計画の後続候補とする。`LinearCombination` の構築には既存の Result の collect を利用し、try_sum は追加しない。

### 3.2. 入力の集約とエラーの順序

処理は次の順序に固定する。

1. 入力を順番に取り出し、最初の Err で停止して、そのエラーを変更せずに返す。その後の入力へ next を呼ばない。
2. 全要素が Ok の場合だけ、集約した式を既存の sum の構築処理へ渡す。局所式の構築結果は Ok で包む。
3. 複合式の空間・束縛の検証は既存の sum と共通化し、その ExpressionError をそのまま返す。

このため、途中の Ok に不整合があっても、その後に入力の Err があれば入力の Err を優先する。全入力が Ok なら、検証エラーが返る場合にも入力生成は最後まで行われる。生成側で発生した副作用を巻き戻す契約は持たない。部分的な和は返さず、primitive の評価・基底展開・数値的な加算は行わない。

例えば A・B を互換でない空間に束縛された演算子、U を未束縛演算子とすると、次の結果になる。

| 入力列 | 結果 |
| --- | --- |
| `[Ok(A), Err(NotAFactor), Ok(B)]` | NotAFactor。3要素目は Iterator から取り出さない |
| `[Ok(A), Ok(B)]` | 合成時の IncompatibleSpaces |
| `[Ok(A), Ok(B), Err(NotAFactor)]` | 入力の NotAFactor を優先。先に空間不整合が潜んでいても、入力エラーがあるため合成検証は行わない |
| `[Ok(A), Ok(B), Ok(U)]` | 成功し、未束縛の和になる。既存 sum の全入力に基づく束縛規則を維持 |

これは、入力エラーと合成エラーを合わせて「入力位置が最も早いエラー」を返す契約ではない。State の和や tensor、演算子積の try 版にも同じ二段階の優先順位を適用する。

既存の sum の契約を維持する。

- 項順、重複、零係数、共有構造を保ち、簡約しない。
- 空・単項入力の扱いは各型の既存 sum と同じとする。局所式は zero、State は指定空間の zero、Operator は未束縛の zero となる。
- State は各項の空間と明示空間の順序付き互換性を検証する。
- Operator は全項が束縛済みのときに互換性を検証する。未束縛項が一つでもあれば和は未束縛となる。

特に Operator を二項 add の fold で構築しない。互換でない二つの束縛済み演算子の後に未束縛演算子が現れる場合、既存 sum は全体を未束縛として受理するため、逐次 add では同じ結果にならない。

集約には既存の `SumExpr<Self>` を利用する方針とし、入力を二度走査したり、式・primitive・係数・Index の Clone を要求したりしない。必要なら既存 sum の検証・構築処理を、集約済み SumExpr を受け取る内部ヘルパーへ切り出して共有する。これにより、公開 sum への委譲のためだけに再収集する処理や、不要な一時 Vec を避ける。空・単項を inline に保持する既存の性質を維持する。

### 3.3. 短絡で止められる範囲

停止の保証は、最初の Err を取り出した後に、入力 Iterator の next を呼ばないことである。`items.into_iter().map(build_operator)` のような遅延生成なら、その後の要素に対する build_operator も呼ばれない。

一方、配列・Vec の各要素を関数呼び出しで先に生成してから try_sum に渡した場合、その生成処理は呼び出し前に実行済みである。try_sum は既に実行された処理を止めたり、副作用を巻き戻したりしない。全入力が Ok の場合は、後続の合成検証が失敗するとしても入力を最後まで生成する。

## 4. 利用例の移行

以下は two_spins の a・b・space・up_down・down_up を使う実装済み API の記述である。

```rust
let singlet = (&up_down - &down_up)?
    .scaled(1.0 / 2.0_f64.sqrt());
let triplet = (&up_down + &down_up)?
    .scaled(1.0 / 2.0_f64.sqrt());

let h = Spin::try_operator_sum(
    [Op::Sx, Op::Sy, Op::Sz].into_iter().map(|axis| {
        Spin::operator_product_at([(&a, axis), (&b, axis)])
    }),
)?.bind(&space)?;
```

利用側の項収集用 Vec をなくし、実数係数を直接渡す。演算子の配置と全作用空間への bind は引き続き明示する。既存の減算構文も使い、API の追加前から短く書ける箇所と、新 API による改善を区別する。

同じ実数のスカラー倍を basis_expansion の例にも適用する。既存の複素係数を使う例は残し、虚部を持つ入力も従来どおり使えることを案内する。モデル utility・正規化などの内部呼び出しは、型推論の確認に必要な範囲で調整し、数学的な契約を変更しない。

## 5. 実装手順

1. 全公開 scaled の呼び出しと既存の型推論・係数能力テストを洗い出す。2.2 節の推論例を小さく確認し、必要な注釈・互換性への影響を特定する。
2. local.rs・composite.rs の4型と linear_combination.rs の scaled を変更する。内部の C 専用処理を維持し、最小能力・変換回数を確認する。
3. 4式型の try_sum を追加する。既存 SumExpr の Result 集約と sum の構築処理を使い、再収集を避けるために必要な処理を内部で共通化する。QuantumSystem の state_sum・try_state_sum・operator_sum・try_operator_sum も追加し、構築例の入口を揃える。
4. 下記の回帰検証を追加し、既存コードで必要な注釈を最小限に補う。
5. two_spins と basis_expansion、examples の README、rustdoc を更新する。製品 API への反映時に primitive-expression・composite-expression・要件の関連箇所へ契約を記録する。
6. 検証を実施して結果を本書へ記録し、未実装という状態を更新する。

## 6. 検証と完了条件

### scaled

- 5対象型で、C そのもの・実数・ユーザ定義の Into 入力を受け取り、既存の C を渡した場合と同じ数学的結果になる。
- f32/f64 の実数・接尾辞なしリテラル・複素係数を検証する。許可しない精度変換は、入力先 C を明示した compile-fail doctest で拒否を確認する。
- 変換回数が一度であることを、複数項・空の LinearCombination も含めて確認する。非 Clone の変換元を受け入れる。
- 式の保持・遅延構築には Clone・数値演算を要求せず、LinearCombination のみ既存の乗算用 bound を維持する。
- 既存の推論例・generic 呼び出しを確認し、破壊的な推論差があれば注釈の移行例とともに記録する。

### try_sum

- 4式型で空・単項・複数項の成功結果が既存 sum と一致する。項順・重複・零係数を保持する。
- 先頭・途中・末尾の Err を確認し、Iterator の観測で最初の Err より後を消費しないことを検証する。
- 遅延 map では Err より後の項生成を行わず、先に配列・Vec を生成した場合は既に行った処理を巻き戻さないことを確認する。
- 入力 Err と空間不整合が共存する場合の入力 Err 優先を確認する。
- State の明示空間・不整合・空の和、Operator の全束縛・不整合・未束縛混在を確認する。不整合な束縛済み2項の後に未束縛項が続くケースも含める。
- 非 Clone の primitive・係数・Index を保持でき、構築時に評価用 trait を要求しないことを確認する。
- QuantumSystem の和の既定メソッドが独自の系でも使え、対応する core constructor と同じ空間・束縛・エラー・数学的結果になることを確認する。

### 横断確認

リポジトリの devShell（必要なら `direnv exec .`）で次を実施する。下記は実装時の検証項目であり、実施結果は末尾で区別して記録する。

```sh
cargo fmt --all -- --check
cargo test --workspace --all-features --offline --locked
cargo test --workspace --no-default-features --offline --locked
cargo clippy --workspace --all-targets --all-features --offline --locked -- -D warnings
cargo doc --workspace --all-features --no-deps --offline --locked
cargo check -p qmb-core -p qmb-geometry -p qmb-spin-half --no-default-features --offline --locked
cargo run -p examples --bin two_spins --offline --locked
cargo run -p examples --bin two_spins --no-default-features --offline --locked
cargo run -p examples --bin basis_expansion --offline --locked
```

two_spins の一重項・三重項のエネルギー −0.75・0.25、単位ノルム、固有状態の残差と内積、および basis_expansion の座標を維持する。既存の正規化・密疎 export の回帰テストも通し、検証環境と未実施のチェックを区別して記録する。新規の production dependency や式ラッパーを追加せず、上記の契約・検証・利用例と移行案内が揃った時点を完了とする。

## 7. 実装結果と検証（2026-10-28）

初回範囲を実装した。QuantumSystem の配置込み入口は `state_tensor_at`・`operator_product_at` とし、構築済みの式は `state_sum`・`operator_sum`・`state_tensor`・`operator_product` と各 `try_` 版から集約する。4式型の `try_sum`、`State::try_tensor`、`Operator::try_product_all` を基礎 API として提供し、同じ検証・構築処理を共有する。後続候補の zero／identity・単独配置・局所式転送メソッド、macro、Mul の拡張は追加していない。

新規・追加の回帰検証は core の `tests/fallible_construction.rs`・`tests/systems.rs`・`tests/coefficient_contracts.rs`・`tests/allocation.rs`、内部の局所／複合式テスト、spin-half の `tests/systems.rs`、core の compile-fail doctest に置いた。全10個の QuantumSystem 追加メソッド、入力エラー優先と短絡、全入力に基づく束縛、最小能力、借用期間と共有、順序・重複・零係数、両精度、変換一回、単項の一時割り当て不要を確認した。

型推論の調整は core と spin-half の generic な代数則テストで必要になった。`LinearCombination::primitive(...).scaled(...)` の係数型と、一部の generic な係数生成関数の型引数を明示した。新たな数値制約や既定精度は導入していない。

Rust 1.98.0 の devShell で次を確認した（Cargo は offline・locked）。

- workspace の全 feature、および default feature 無効のテスト・doctest が成功。
- 全 target・全 feature の Clippy（`-D warnings`）、format、警告をエラー扱いにしたドキュメント生成が成功。
- core・geometry・spin-half の default feature 無効の単独 check が成功。
- two_spins を default／default feature 無効で実行し、エネルギー −0.75・0.25、単位ノルム、ゼロの残差と内積を確認。
- basis_expansion を実行し、状態座標 1:5・2:−1 と ladder 基底の一成分を確認。
- 既存の独立行列比較・正規化・密疎 export の両精度の回帰テストが成功。

MSRV 1.85 と専用 no_std target は今回も未検証であり、no_std crate の確認はホストでの単独 check とソースの契約にとどまる。依存 manifest と lockfile は変更していない。
