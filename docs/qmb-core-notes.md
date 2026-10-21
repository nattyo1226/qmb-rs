# qmb-core 実装検討メモ

更新日: 2026-10-21

本書は、実装計画の step 1 に向けた一時的な検討メモである。議論した選択肢と現在の合意を区別し、実装に反映するまでの作業用記録とする。恒久的な仕様書として維持することは想定せず、他のドキュメントからはリンクしない。ここに記した型・API が実装済みであることは意味しない。

## 1. core の内部構造・API の検討対象

`qmb-core` は抽象契約と最小限の共通型を担当し、共通の式木・評価は `qmb-algebra`、具体的な局所代数は専用 crate に置く。この境界を前提に、以下を具体化する。

| 対象 | 検討状況・決定すべき事項 |
| --- | --- |
| 1. module と公開範囲 | 決定済み。`no_std + alloc`、5 module、主要型と数値 crate の再公開。1.1・1.2 節 |
| 2. 数値契約 | 決定済み。`Complex<T>` と操作単位の制約。実数値の評価ポリシーは上位層へ引き継ぐ。2・3 節 |
| 3. 局所構造 | 決定済み。`LocalStructure<T>` の関連型、型ごとの正の固定次元、実行時パラメータなし。4 節 |
| 4. 局所構造の一致 | 決定済み。同じ Rust 型で表し、`Definition` は導入しない。別実装間の照合は今回の対象外。4・8 節 |
| 5. primitive | 決定済み。精度から独立した型、有効性は具体型の構築時に保証。一律の `Clone` などは要求しない。4・7 節 |
| 6. 数学的能力 trait | 決定済み。内積・積・作用・随伴・trace を分離し、恒等演算子は局所構造から取得。4 節 |
| 7. 有限線型結合 | `LinearCombination<P, T>` の表現、構築・列挙・加算・係数乗算、trait 制約、構造比較は合意済み。5 節に記録する |
| 8. 局所基底 | 決定済み。展開・順序番号・遅延列挙、対象別 marker trait、次元からの要素数導出。6 節 |
| 9. エラーと検証 | 決定済み。core は次元計算の桁あふれを扱い、`thiserror` を使う。不正値の構築は具体 crate の責務。7 節 |
| 10. 異種構造への接続 | 決定済み。core は静的 dispatch。共有・型消去・実行時照合は上位層。7・8 節 |
| 11. 契約の検証 | 方針決定済み。最小制約、独立した局所代数、法則、MSRV・`no_std` を実装時に確認。8 節 |

有限線型結合に未評価の積・作用・束縛節点や evaluator は持たせない。空の項列と、primitive の線型従属性によって数学的に零になる線型結合は区別する。Index・Hilbert 空間・完全な式木・複合基底・数値 export の実装は core の対象に含めない。

今回の core の実装に先立って利用者に選択を求める設計上の未決事項は、現時点ではない。実装による検証事項と、上位層・将来の機能に残す事項は 8 節に分けて記録する。

### 1.1. module と公開範囲

crate は `#![no_std]` と `extern crate alloc;` を使い、複数項の保持に `alloc::vec::Vec` を利用する。`std` を既定で要求せず、初期実装では独自 feature や prelude は設けない。実行ファイル側で必要になる allocator の提供は利用側の責務とする。

| module | 主な公開対象 |
| --- | --- |
| `linear_combination` | `LinearCombination`・`Term`・`IntoTerms` |
| `ops` | `InnerProduct`・`Product`・`ActOn`・`Adjoint`・`Trace` |
| `local_structure` | `LocalStructure` |
| `basis` | `Basis`・`StateBasis`・`OperatorBasis`・基底要素数の補助関数 |
| `error` | `DimensionOverflow` |

主要な公開型・trait・関数は crate root でも再公開する。`Complex` は `pub use num_complex::Complex;`、数値 trait は `pub use num_traits;` で参照できるようにする。内部表現や補助実装は公開しない。

### 1.2. workspace と依存関係

依存バージョンと既定 feature は workspace ルートの `Cargo.toml` で管理する。初期設定は以下とする。このコード例は今後の実装用であり、今回 Cargo ファイルは変更しない。

```toml
[workspace.package]
rust-version = "1.85"

[workspace.dependencies]
num-complex = { version = "0.4.6", default-features = false }
num-traits = { version = "0.2.19", default-features = false }
thiserror = { version = "2.0.21", default-features = false }
```

既存の workspace の edition・version を引き継ぎ、個別 crate では次のように利用する。

```toml
[package]
name = "qmb-core"
version.workspace = true
edition.workspace = true
rust-version.workspace = true

[dependencies]
num-complex.workspace = true
num-traits.workspace = true
thiserror.workspace = true
```

初期実装では `libm` を有効にしない。上位 crate が追加 feature を必要とする場合もバージョンは workspace から継承する。core 単独の依存閉包が `std` を要求しないことを検証する。

workspace での継承方法は [Cargo の依存関係の共有](https://doc.rust-lang.org/cargo/reference/workspaces.html#the-dependencies-table)、feature の設定は [num-complex](https://docs.rs/crate/num-complex/0.4.6/features)・[num-traits](https://docs.rs/crate/num-traits/0.2.19/features)・[thiserror](https://docs.rs/crate/thiserror/2.0.21/features) を参照する。Rust 1.85 は既存の edition 2024 に合わせる。[Rust 2024](https://doc.rust-lang.org/edition-guide/rust-2024/index.html)

## 2. 数値契約の現在の合意

### 2.1. 係数を Complex<T> に統一する

係数は当面 `num_complex::Complex<T>` に統一する。実数型 `T` は汎用化し、単精度・倍精度の複素数を標準の対応・検証対象とする。独自の複素数表現を直接受け入れる抽象化は、現在の構成では導入しない。

| 必要なもの | 利用する表現・API |
| --- | --- |
| 対応する実数型 | `Complex<T>` の型パラメータ `T` |
| 実部・虚部 | `re`・`im` |
| 実数から複素数への埋め込み | `Complex::new(value, T::zero())` |
| 共役・複素数演算 | `Complex` の既存 API・演算実装 |
| 絶対値の二乗 | `T: Clone + Num` のもとで `norm_sqr()` を利用可能 |
| 絶対値・平方根 | 必要な制約と feature を満たす上位層で `norm()`・`sqrt()` などを検討 |

実数型との対応付け用に検討した `HasReal`、実数の埋め込み用に検討した `FromReal` は採用しない。複素数型を差し替える公開 trait を設計する代わりに、`Complex` の API を使って型の対応関係を処理する。

同じ演算の係数には同じ実数型 `T` を要求し、単精度・倍精度の混在と暗黙の精度昇格は認めない。精度変更は明示的に行う。値ごとに精度を持つ型を対象にする場合の検証方法は、将来の拡張事項とする。

### 2.2. 能力は操作ごとに要求する

実数型には、利用する操作と `Complex` の実装に必要な能力を要求する。例えば既存の複素数の乗算には `T: Clone + Num`、共役にはさらに `Neg<Output = T>` が必要となる。全操作に一律の `Float`・`Copy`・数学関数の実装を要求する方針にはしない。

`num_traits::Num` は、以下を要求する。

| 要求 | 内容 |
| --- | --- |
| `PartialEq` | 等値比較 |
| `Zero` | `zero()`・`is_zero()` |
| `One` | `one()` |
| `NumOps` | 加減乗除と剰余。既定では入力・出力とも自身の型 |
| `FromStrRadixErr`・`from_str_radix()` | 基数を指定した文字列からの変換とエラー型 |

`Num` 自体には `Clone`・`Copy`・`Neg`・大小比較・平方根・絶対値・数値型間の変換は含まれない。整数型も実装しており、実数としての除算などの数学的性質を保証する契約ではない。型制約は操作を記述するための必要条件であり、選んだ数値型と局所規則が必要な数学的性質を満たすことは実装者の契約とする。core はこれを新たな数値 marker trait や実行時検査で証明しようとはしない。

剰余と文字列変換は qmb の基本的な代数処理に直接必要な能力ではないが、`Complex<T>` の既存演算を使うための制約として `Num` を受け入れる。

`std`・`libm` を無効にした初期設定では、`Complex` の `norm()`・`sqrt()` は利用できない。core はこれらを呼ばず、必要な基本演算だけを利用する。`norm_sqr()` は数学関数の feature を必要としない。[num-complex の実装](https://docs.rs/num-complex/0.4.6/src/num_complex/lib.rs.html)

## 3. 上位層へ引き継ぐ実数値の評価

複素絶対値とその二乗は、本質的に実数値を返す演算である。近似比較でも、複素絶対値・許容誤差・大小比較を実数型で扱う。これらを、一般の複素数の虚部を捨てて実数と解釈する操作とは区別する。

### 3.1. norm と正規化

状態の norm の数学的な定義は `sqrt(⟨ψ,ψ⟩)` とする。非零状態では自己内積は正の実数だが、汎用の内積が返す複素数値には数値誤差による虚部が残り得る。

norm などで結果を実数として扱う際には、虚部を検査する方針とする。`re` の取得だけでは検査にならない。検査を通過した値の実部を用いる補助操作を検討するが、具体的な API・許容条件・失敗時の返却形式は未決である。

実数 norm を用いる正規化の候補は、自己内積の計算、虚部の検査、実部の非負性の確認、実数の平方根、零ノルムの検出、実数 norm による状態の除算という流れになる。複素係数が必要な箇所では、実数の逆数を `Complex::new` で埋め込むこともできる。評価手順と必要な trait 制約はまだ確定していない。

非独立 primitive・非正規直交基底を許すため、一般の状態の norm を、係数の絶対値の二乗の和だけから計算してはならない。交差項を含む内積を扱い、通常の norm 評価に指定基底への全展開を強制しない。

複素平方根によって、自己内積から複素係数のまま正規化する案も検討した。ただし採用は決定しておらず、正規化のために複素平方根の能力を必須にする合意はない。

### 3.2. qmb-algebra で決める数値ポリシー

以下は norm・正規化などの評価を設計する際に決める。今回の core に検査 API や数値ポリシーは追加しない。

- 虚部の検査を厳密比較・許容誤差のどちらで行うか、許容誤差の指定方法と判定式。`|Im z| <= atol + rtol |Re z|` は検討案であり、確定していない。
- 負の実部と数値誤差の扱い。小さな負値を自動的に零へ補正する方針は決めていない。
- 検査時の NaN・infinity の扱い。通常の算術を数値型に委譲する既存の契約との関係も明確にする。
- 検査失敗・零ノルムなどのエラー型と、検査 API の配置・公開範囲。
- 平方根、絶対値、比較、定数生成、明示的な精度変換に必要な具体的な制約。

検査・比較の許容誤差を、通常の式評価での項の削除に流用しない。一般の内積・trace・期待値など、複素数値を返す演算を一律に実数化する方針ではない。

## 4. 局所構造と primitive の関連付け

### 4.1. 局所構造側の関連型を使う

`LocalStructure<T>` が、状態・演算子 primitive を関連型で指定する。局所構造型と実数型 `T` が決まれば、その組に対応する primitive 型と必要な演算能力を参照できる構成とする。primitive 側から局所構造を逆向きに関連付ける trait は、いったん導入しない。

採用する局所構造の定義は以下とする。

```rust
pub trait LocalStructure<T> {
    type State: InnerProduct<T>;
    type Operator: Product<T> + ActOn<Self::State, T> + Adjoint<T> + Trace<T>;

    const DIMENSION: core::num::NonZeroUsize;

    fn identity() -> LinearCombination<Self::Operator, T>;
}
```

局所構造が対応する primitive と必要な能力の組を指定し、個々の規則は primitive 側で実装する。Product・ActOn は演算子 primitive 側の能力とする既存の責任分担を維持する。

### 4.2. 構造の値と数値精度を分離する

局所構造のインスタンスに実行時の定義パラメータを持たせない。次元や代数の規約は局所構造の型ごとに決まり、異なる構造は別の型で定義する。例えば `SpinHalf` と `SpinOne` を区別する。局所規則に構造のインスタンスを渡す構成は要求しない。primitive 自身が状態のラベルなどの値を持つことは、この判断とは別である。

局所構造・primitive の具体型は数値型 `T` を持たず、primitive に実装する演算 trait が `T` を持つ。例えば同じ `SpinHalf`・`SpinState`・`SpinOperator` を用いて、単精度・倍精度それぞれの規則を実装できる。`LocalStructure<T>` はその精度で必要な規則が実装されていることを要求する。

局所構造が完全な式木に触れるために `T` を必要とするわけではない。局所規則は係数 `Complex<T>` を使った有限線型結合を返し、完全な式木への取り込みは `qmb-algebra` が担当する。

`DIMENSION` は状態空間の有限な正の次元とし、`NonZeroUsize` で零を除外する。同じ局所構造型では、`T` を変えても次元と `State`・`Operator` の具体型を変えない。`LocalStructure<T>` の形式だけではこの一致は強制できないため、実装者の契約として記録し、対応する精度間のテストで確認する。

定義識別用の `Definition` 関連型・marker・独自 ID は導入しない。同じ構造の利用には同じ Rust 型を使う。別の Rust 実装型同士を同じ数学的構造として照合する仕組みは今回の対象外とする。

### 4.3. primitive に実装する能力

```rust
pub trait InnerProduct<T> {
    fn inner_product(&self, rhs: &Self) -> Complex<T>;
}

pub trait Product<T>: Sized {
    fn product(&self, rhs: &Self) -> LinearCombination<Self, T>;
}

pub trait ActOn<S, T> {
    fn act_on(&self, state: &S) -> LinearCombination<S, T>;
}

pub trait Adjoint<T>: Sized {
    fn adjoint(&self) -> LinearCombination<Self, T>;
}

pub trait Trace<T> {
    fn trace(&self) -> Complex<T>;
}
```

入力は借用し、結果は所有値で返す。`a.product(&b)` は演算子積 `ab`、すなわち状態への作用では右の `b` を先に適用する順序とする。内積は左引数について共役線型、右引数について線型とする。`identity()` は恒等演算子の有限線型結合を返し、単一 primitive で表せることを要求しない。

局所規則は有効な primitive に対して常に結果を返し、`Result` を返さない。不正なラベルなどは具体型の constructor で拒否し、有効な値だけを規則に渡せるようにする。全 primitive の列挙や実行時の代数法則の検証を core の API に要求しない。法則を破る実装は契約違反として扱う。

### 4.4. 既存の設計ドキュメントとの範囲の違い

既存の設計には、実行時の定義パラメータや、異なる実装にまたがる数学的な定義の一致も含まれる。今回の step 1 では、ここに記録した型ごとの固定構造と静的な関連付けに範囲を絞る。既存ドキュメントは今回変更せず、それらの広い構想を core の必須 API として持ち込まない。

## 5. 局所規則が返す有限線型結合

### 5.1. 浅い式木として保持する

`qmb-core` の返却型を `LinearCombination<P, T>` とし、状態・演算子の双方で同じ容器を使う。`P` は対象の primitive、係数は `Complex<T>` とする。

```rust
use alloc::vec::Vec;
use num_complex::Complex;

pub struct LinearCombination<P, T> {
    repr: Repr<P, T>,
}

pub struct Term<P, T> {
    pub coefficient: Complex<T>,
    pub primitive: P,
}

enum Repr<P, T> {
    Zero,
    Term(Term<P, T>),
    Sum(Vec<Term<P, T>>),
}
```

`Repr` は非公開とし、項数 0・1・2 以上に応じて `Zero`・`Term`・`Sum` を選ぶ。単一項の保持にヒープ確保は必要ない。再帰的な和や scalar multiplication の節点は導入しない。

構築時には項順序・重複項・零係数を保持し、同類項の結合や零係数の削除は行わない。`Zero` と、係数が零の一項は異なる構造として保持する。

### 5.2. 構築・列挙・演算 API

以下は `LinearCombination<P, T>` の API とする。数値型の制約は型そのものではなく、それを使う操作に付ける。

| API | 契約・制約 |
| --- | --- |
| `zero() -> Self` | 明示的な零。追加制約なし |
| `primitive(p: P) -> Self` | 係数 `1` の一項。`T: Zero + One`。係数は `Complex::new(T::one(), T::zero())` で作る |
| `term(coefficient: Complex<T>, primitive: P) -> Self` | 係数付き一項。追加制約なし |
| `from_terms(terms: impl IntoIterator<Item = Term<P, T>>) -> Self` | 入力順で構築し、項数に応じて表現を選ぶ。追加制約なし |
| `terms(&self) -> core::slice::Iter<'_, Term<P, T>>` | 項を借用して列挙する。追加制約なし |
| `into_terms(self) -> IntoTerms<P, T>` | 項を所有権ごと取り出す。追加制約なし |
| `len(&self) -> usize` | 保持する項数。零係数も一項として数える |
| `is_empty(&self) -> bool` | 項列が空かを返す。数学的な零判定ではない |
| `scaled(self, coefficient: Complex<T>) -> Self` | 各項の係数を乗算する。`T: Clone + Num`。`P: Clone` は要求しない |

`Zero`・単一項は空または一要素の slice として借用できる。所有列挙用の `IntoTerms<P, T>` は公開 iterator 型とし、その内部表現は非公開にする。零・一項の所有列挙でヒープ確保せず、複数項では保持している `Vec` の所有権を使う。

両 iterator は入力順で列挙し、`DoubleEndedIterator`・`ExactSizeIterator`・`FusedIterator` に対応する。

### 5.3. 標準 trait と構造比較

| trait | 契約・制約 |
| --- | --- |
| `Clone` | `P: Clone, T: Clone` の場合に提供する |
| `Debug` | `P: Debug, T: Debug` の場合に提供する |
| `PartialEq` | `P: PartialEq, T: PartialEq` の場合に提供する |
| `Default` | 追加制約なしで `zero()` を返す |
| `Add<Self, Output = Self>` | 所有値同士の項を左・右の順に連結する。係数加算・primitive 比較・clone は要求しない |
| `FromIterator<Term<P, T>>` | `from_terms` と同じ構築処理を使う |
| `IntoIterator` | 所有値では `Term<P, T>`、共有参照では `&Term<P, T>` を入力順で返す |

`Term` にも、同じ条件で `Clone`・`Debug`・`PartialEq` を提供する。`Eq`・`Hash` や借用値同士の加算は、初期の必須 API に含めない。

等値比較は、項数・項順序・各係数・各 primitive の厳密な構造比較とする。係数比較は採用する数値型の `PartialEq` に従い、NaN を含む係数は自身との比較でも等しくならない場合がある。数学的な零判定の `is_zero()` は提供しない。

| 比較する表現 | 構造比較 |
| --- | --- |
| 同じ順序で、係数・primitive がそれぞれ等しい項を持つ表現 | 等しい |
| `a·p + b·q` と `b·q + a·p` | 項列が異なる場合は等しくない |
| `p + p` と `2·p` | 二項と一項なので等しくない |
| `Zero` と `0·p` | 空と一項なので等しくない |

### 5.4. 未評価式との境界と検証対象

`ActOn(X, ket(0))` のような未評価の作用は、Index を付与する前の局所式でも、付与後の複合式でも `qmb-algebra` の式木で保持する。積・作用・随伴などの式を構築した時点で局所規則を呼ぶことは要求しない。

必要な評価時に局所規則を呼び、結果の `LinearCombination` を共通式木へ取り込む。元の source は保持する。この返却型に、未評価の積・作用・随伴・束縛節点や evaluator は追加しない。

実装時には、0・1・複数項の構築、項順序・重複項・零係数の保持、両 iterator の動作、所有操作が primitive の clone を要求しないこと、および構造比較を検証する。虚部検査・norm の数値ポリシーは上位層へ引き継ぐ。

## 6. 局所基底

### 6.1. 展開・番号・列挙

```rust
pub trait Basis<P, T> {
    fn expand(&self, primitive: &P) -> LinearCombination<P, T>;
    fn index(&self, primitive: &P) -> Option<usize>;
    fn elements(&self) -> impl Iterator<Item = P> + '_;
}

pub trait StateBasis<L, T>: Basis<L::State, T>
where
    L: LocalStructure<T>,
{
}

pub trait OperatorBasis<L, T>: Basis<L::Operator, T>
where
    L: LocalStructure<T>,
{
}
```

`Basis` は primitive を、同じ primitive 型で表した基底要素の有限線型結合へ展開する。`StateBasis`・`OperatorBasis` は対象となる局所構造と状態・演算子の区別を付ける marker trait とする。primitive 自体が基底であることや、基底が正規直交であることは要求しない。

`index()` は基底に含まれる要素の零始まりの番号を返す。有効な primitive でも、その基底の要素でなければ `None` となる。`elements()` は各基底要素を番号順に一度だけ、所有値として返す。事前の `Vec` 構築や `Box` は要求せず、遅延列挙を可能にする。iterator に `ExactSizeIterator`・`DoubleEndedIterator` は必須としない。

有効な基底は対象空間の基底を成し、`expand()` は有効な任意の対象 primitive に対して定義される。展開による数学的な値の保存、基底要素と番号・列挙順の一致、必要な要素数を実装者の契約とする。設定可能な基底がある場合、不正な設定は具体的な基底型の構築時に拒否する。汎用の基底検証や失敗する `expand()` は追加しない。

### 6.2. 要素数と型推論

公開する補助関数は以下とする。

```rust
pub fn state_basis_len<L, T>() -> usize
where
    L: LocalStructure<T>;

pub fn operator_basis_len<L, T>() -> Result<usize, DimensionOverflow>
where
    L: LocalStructure<T>;
```

状態基底の要素数は `L::DIMENSION.get()`、演算子基底の要素数はその二乗とする。二乗は `checked_mul` で計算し、`usize` の範囲を超えたら `DimensionOverflow` を返す。基底ごとに重複した `len()` や番号から要素を取得する API は初期実装に設けない。

`index()`・`elements()` の入出力や定数の参照だけでは `T` を推論できない場合がある。その場合は `<B as Basis<P, T>>::index(&basis, &primitive)`、`<B as Basis<P, T>>::elements(&basis)`、`<L as LocalStructure<T>>::DIMENSION` のように明示する。このための別の精度非依存 trait は導入しない。

## 7. エラー・所有権・上位層との接続

### 7.1. エラーの責務

core で共通に返すエラーは、演算子空間の次元計算の桁あふれとする。エラー定義には `thiserror` を利用し、包括的な `CoreError` は設けない。

```rust
#[derive(Debug, Clone, Copy, PartialEq, Eq, thiserror::Error)]
#[error("operator-space dimension overflows usize")]
pub struct DimensionOverflow;
```

`thiserror` の `std` feature を無効にし、`core::error::Error` に対応させる。不正な primitive・基底設定の構築エラーは具体 crate、式の評価・実行時の互換性・数値検査のエラーは `qmb-algebra` で定義する。有効な入力に対して数学的法則を破る局所規則は、利用者の入力エラーとして扱わない。

### 7.2. 制約を局所化する

primitive や数値型に、一律の `Clone`・`Copy`・`Debug`・`PartialEq`・`Send`・`Sync`・`'static` は要求しない。各操作・標準 trait の実装が必要とする箇所だけに制約を付ける。例えば構造比較は対応する `PartialEq` を必要とするが、線型結合の構築・所有加算・所有列挙では primitive の比較や clone は不要である。

core は静的 dispatch を前提とし、`dyn` で利用できることを全 trait の要件としない。`Rc`・`Arc`・`Box<dyn ...>` を保持方式として要求せず、異種構造を含む式木の型消去・共有・実行時照合は `qmb-algebra` が必要な制約を付けて実装する。

## 8. 未決事項の再点検と実装時の検証

### 8.1. 今回の core で実装確認する事項

ここまでの API と責務分担で実装を開始できる。以下は実装して確かめる事項であり、新たな利用者の選択を待つ項目ではない。コンパイルで矛盾や不要な制約が見つかった場合は、実装とこのメモを更新する。

- `LocalStructure`・能力 trait・基底 trait の関連型と制約を、公開 API だけを使う独立したテスト用局所代数で確認する。`T` によらない primitive 型と次元も確認する。
- `Copy` ではない `Clone + Num` の数値型、および `Clone` を実装しない primitive で、操作ごとの最小制約を確認する。
- 有限線型結合の 0・1・複数項、項順序・重複・零係数の保存、iterator の前後列挙・長さ・終端、構造比較を確認する。
- 三次元の状態と行列単位 `E_ij` などを使い、積・作用・内積・随伴・trace・複数項からなる恒等演算子の法則を確認する。非独立 primitive と非正規直交基底の展開も含める。
- 異なる精度の線型結合の混在や、対応しない状態への作用が型エラーになることを compile-fail の doc test で確認する。
- 基底の展開・番号・列挙・要素数の整合、次元の二乗の桁あふれ、エラーの表示と `core::error::Error` への適合を確認する。

実装後に行う主なチェックは次のとおりとする。

```sh
cargo test -p qmb-core
cargo fmt --all -- --check
cargo clippy -p qmb-core --all-targets -- -D warnings
cargo +1.85.0 check -p qmb-core --lib
cargo check -p qmb-core --lib --target thumbv7em-none-eabi
```

MSRV 用 toolchain と対象 target が利用可能な環境で実施する。最後のチェックでは `std` のない target 向けにライブラリ単独でビルドし、依存閉包も検証する。最終実行ファイルの allocator はこのチェックの対象ではない。[対象 target](https://doc.rust-lang.org/rustc/platform-support/thumbv7em-none-eabi.html)

ホスト側のテストでは `std` を利用できるが、テストのためだけに通常依存の `std` feature を有効にしない。数値例の確認にはホスト側の実数演算や feature 不要の `norm_sqr()` を使い、core が `Complex::norm()`・`sqrt()` を呼ばずに成立することを保つ。

### 8.2. qmb-algebra の設計に残す事項

- 3 節の虚部検査・許容誤差・負値・特殊値・norm・正規化と、それらに必要な数値能力・エラー API。
- 未評価式、局所規則の適用と有限線型結合の取り込み、source の保持。
- 異種構造の型消去、共有、必要な lifetime・スレッド間共有の制約、実行時の互換性検査。

### 8.3. 今回の実装範囲外とする拡張

- 別の Rust 実装型にまたがる同一定義の照合、永続化用の安定した識別子。
- 値ごとに精度を持つ数値型の整合性検査。
- graded な代数に固有の契約や、並列評価の追加契約。

これらのために `Definition` や汎用の識別・数値検査 API を先行して core に追加しない。今回の再点検では、上位層や将来の要求を除いた core の設計に新たな未決事項は見つかっていない。ただし型制約・MSRV・`no_std` の成立は、今後の実装と上記チェックによる確認を必要とする。

## 9. 既存ライブラリから確認したこと

- `num-complex::ComplexFloat` は `type Real`・共役・絶対値などを提供するが、sealed で独自型に直接実装できず、`Copy` や多数の数学関数も要求する。独自複素数型のための最小契約には適さなかった。[ComplexFloat](https://docs.rs/num-complex/latest/num_complex/trait.ComplexFloat.html)
- `simba::ComplexField` は `RealField` を関連付け、実部取得・絶対値・実数による拡大縮小などを提供するが、数学関数・表示・スレッド間共有なども含む広い契約である。[ComplexField](https://docs.rs/simba/latest/simba/scalar/trait.ComplexField.html)
- `cauchy::Scalar` も実数型・複素数型を関連付けるが、`Copy`・シリアライズ・数学関数などを要求する。[Scalar](https://docs.rs/cauchy/latest/cauchy/trait.Scalar.html)
- 確認した `nalgebra` の実装では、norm の二乗は自己内積の実部を取得して求め、norm は実数の平方根、正規化は実数 norm による除算を使う。`simba` の `Complex<N>` 向け `real()` は `self.re` を返すだけで、虚部を検査しない。[norm の実装](https://docs.rs/nalgebra/latest/src/nalgebra/base/norm.rs.html#238-262)、[正規化の実装](https://docs.rs/nalgebra/latest/src/nalgebra/base/norm.rs.html#371-380)、[実部取得の実装](https://docs.rs/simba/latest/src/simba/scalar/complex.rs.html#1035-1038)
- `ndarray-linalg` は成分の絶対値の二乗を `A::Real` として加算して平方根を取り、正規化では実数 norm を `from_real` で係数型に埋め込む。[実装](https://docs.rs/crate/ndarray-linalg/latest/source/src/norm.rs)
- `faer` の norm の二乗も `T::Real` を返し、一般の計算経路では成分ごとの `abs2()` を実数として加算する。[実装](https://docs.rs/faer/latest/src/faer/linalg/reductions/norm_l2_sqr.rs.html)

型が実数型の戻り値を保証することと、その実装が数学的な契約を満たすことは別である。既存実装は設計の参考とし、qmb の虚部検査の要否や非正規直交な状態への評価手順をそのまま決めるものではない。

数値型の制約を実装するときの参照先は、[Complex](https://docs.rs/num-complex/latest/num_complex/struct.Complex.html)、[Num](https://docs.rs/num-traits/latest/num_traits/trait.Num.html)、[NumOps](https://docs.rs/num-traits/latest/num_traits/trait.NumOps.html) とする。今回の更新では Rust コードや Cargo の依存関係は変更しない。
