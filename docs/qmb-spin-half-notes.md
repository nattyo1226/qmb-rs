# qmb-spin-half 実装検討メモ

更新日: 2026-10-21

本書は、実装計画の step 1 における `qmb-spin-half` の実装までの一時メモである。合意した API・局所規則・実装順序・検証方針を記録する。恒久的な仕様書として維持することは想定せず、他のドキュメントからはリンクしない。ここに記した API が実装済みであることは意味しない。

## 1. crate の責任と構成

`qmb-spin-half` は、既存の `qmb-core` の公開契約を実装する具体的な局所代数とする。局所構造、状態・演算子 primitive、局所規則、および三種類の局所基底を提供する。

crate は `#![no_std]` とする。局所規則は core の `LinearCombination` を返し、複数項に必要な確保はこの容器に委ねる。通常依存は `qmb-core` のみとし、`Complex`・数値 trait は core の再公開を利用する。`std`・`libm` を要求する feature は追加しない。

将来の実装時に `crates/qmb-spin-half` を workspace member に加え、workspace の依存定義へ次を追加する。

```toml
[workspace.dependencies]
qmb-core = { path = "crates/qmb-core" }
```

個別 crate は workspace の version・edition・rust-version を継承し、依存は `qmb-core.workspace = true` とする。

| module | 主な公開型 |
| --- | --- |
| `local_structure` | `SpinHalf` |
| `state` | `SpinHalfState` |
| `operator` | `SpinHalfOperator` |
| `basis` | `SzBasis`・`CartesianBasis`・`LadderBasis` |

上記の型は crate root でも再公開する。係数生成などの補助処理は非公開とする。

共通式木、未評価節点、式全体の評価・簡約、Index・配置・束縛、複合基底、norm・正規化、数値 export は上位層の責務とする。`qmb-algebra` への依存を追加しない。保留中の `+`・`-`・`*` による演算子構文の追加は、この実装計画に含めない。

## 2. 公開型と core への接続

| 公開型 | 定義・役割 |
| --- | --- |
| `SpinHalf` | 実行時のパラメータを持たない unit struct |
| `SpinHalfState` | `Up`・`Down` の enum |
| `SpinHalfOperator` | `I`・`Sx`・`Sy`・`Sz`・`SPlus`・`SMinus` の enum |
| `SzBasis` | 状態基底を表す unit struct |
| `CartesianBasis` | Cartesian 演算子基底を表す unit struct |
| `LadderBasis` | 昇降演算子基底を表す unit struct |

各型は `Clone`・`Copy`・`Debug`・`PartialEq`・`Eq`・`Hash` に対応する。unit struct には `Default` も提供する。状態・演算子 enum に既定値や実行時パラメータは設けず、公開 variant から有効な値だけを構築できるようにする。固有のエラー型は不要とする。

状態・演算子の具体型は数値型 `T` を持たず、規則の実装が `T` を持つ。`LocalStructure<T> for SpinHalf` は次を関連付ける。

- `State = SpinHalfState`
- `Operator = SpinHalfOperator`
- `DIMENSION = NonZeroUsize::new(2).unwrap()`
- `identity()` は係数 `1` の `SpinHalfOperator::I`

状態には `InnerProduct<T>`、演算子には `Product<T>`・`ActOn<SpinHalfState, T>`・`Adjoint<T>`・`Trace<T>` を実装する。入力・出力のシグネチャは core の trait に従う。

係数型を入力から推論できない呼び出しでは、返却型を指定するか、`<SpinHalfOperator as Product<f64>>::product(...)` のように明示する。精度を primitive に埋め込んだり、暗黙に既定精度へ変換したりしない。

## 3. 数値契約

係数は `Complex<T>` とする。型の保持自体に数値制約は加えず、各能力の実装に必要なものだけを要求する。

| 実装 | `T` の制約 |
| --- | --- |
| `InnerProduct`・`Adjoint`・`Trace`・`identity()`・`Basis<SpinHalfState, T> for SzBasis` | `Zero + One` |
| `Basis<SpinHalfOperator, T> for CartesianBasis` | `Zero + One + Neg<Output = T>` |
| `Product`・`ActOn`・`Basis<SpinHalfOperator, T> for LadderBasis` | `Zero + One + Div<Output = T> + Neg<Output = T>` |
| `LocalStructure<T> for SpinHalf` | 全能力を成立させる `Zero + One + Div<Output = T> + Neg<Output = T>` |

`Zero`・`One` は `qmb_core::num_traits`、`Div`・`Neg` は `core::ops` を使う。`StateBasis<SpinHalf, T>`・`OperatorBasis<SpinHalf, T>` の marker 実装には、対応する `LocalStructure<T>` を成立させる制約も付ける。

`Zero` は加算、`One` は乗算の能力も要求する。`2` は `one() + one()`、`1/2` は `one() / (one() + one())` で生成する。`1/4` も新たに生成した `one()` の和を分母にして構成する。`i`・負号は実部・虚部と `Neg` で表現し、係数の生成に `Float`・`Copy`・`Clone`・`Num`・文字列変換・数値キャストを一律には要求しない。

これらは必要な操作の型制約であり、数値型の数学的性質を証明するものではない。`T` は実数 scalar として `1/2`・`1/4` を適切に表せる必要がある。整数除算による代用は想定しない。単精度・倍精度を標準の対応・検証対象とし、それ以外の型も契約を満たせば利用できる構成とする。

## 4. 数学的規約と局所規則

### 4.1. 規約

`hbar = 1`、`Sα = σα/2`、`SPlus = Sx + i Sy`、`SMinus = Sx - i Sy` とする。状態基底の順序は `Up, Down`、それぞれの `Sz` 固有値は `+1/2, -1/2` とする。

局所積は左引数と右引数の順序を保存する。`a.product(&b)` は `ab` であり、状態へ適用する場合は右の `b` が先に作用する。

規則は有効な primitive に対して常に結果を返す。数学的に零と分かっている局所結果は `LinearCombination::zero()` とする。非零結果では、規則から直接必要な項を構築する。返却後の容器に同類項の結合や零判定を要求しない。

### 4.2. 局所積

局所積は、最も自然な primitive の線型結合で返す。例えば `I·SPlus` は `SPlus`、`SPlus·SMinus` は `I/2 + Sz` とし、特定の基底への全展開は別操作にする。一般的な最短表現の探索や最適化機構は導入しない。

以下の `X,Y,Z,P,M` はそれぞれ `Sx,Sy,Sz,SPlus,SMinus` を表す。行が左引数、列が右引数である。全36組をこの表の直接規則として実装する。

| × | I | X | Y | Z | P | M |
| --- | --- | --- | --- | --- | --- | --- |
| I | I | X | Y | Z | P | M |
| X | X | I/4 | iZ/2 | −iY/2 | I/4−Z/2 | I/4+Z/2 |
| Y | Y | −iZ/2 | I/4 | iX/2 | iI/4−iZ/2 | −iI/4−iZ/2 |
| Z | Z | iY/2 | −iX/2 | I/4 | P/2 | −M/2 |
| P | P | I/4+Z/2 | iI/4+iZ/2 | −P/2 | 0 | I/2+Z |
| M | M | I/4−Z/2 | −iI/4+iZ/2 | M/2 | I/2−Z | 0 |

二項の結果は表の記載順で返す。積の実装で指定基底への展開や行列乗算を必須にせず、行列表現は独立したテストの参照計算に使う。

### 4.3. 状態への作用

| 演算子 | Up への作用 | Down への作用 |
| --- | --- | --- |
| `I` | `Up` | `Down` |
| `Sx` | `Down/2` | `Up/2` |
| `Sy` | `i Down/2` | `−i Up/2` |
| `Sz` | `Up/2` | `−Down/2` |
| `SPlus` | `0` | `Up` |
| `SMinus` | `Down` | `0` |

消滅する作用は空、その他は係数付き一項の線型結合で返す。未評価の作用を表す状態 primitive は追加しない。

### 4.4. 内積・随伴・trace

内積は左引数に共役線型、右引数に線型である。primitive の内積は次の正規直交規則とする。

| 左引数 | 右引数 | 内積 |
| --- | --- | --- |
| `Up` | `Up` | `1` |
| `Up` | `Down` | `0` |
| `Down` | `Up` | `0` |
| `Down` | `Down` | `1` |

| 演算子 | 随伴 | trace |
| --- | --- | --- |
| `I` | `I` | `2` |
| `Sx` | `Sx` | `0` |
| `Sy` | `Sy` | `0` |
| `Sz` | `Sz` | `0` |
| `SPlus` | `SMinus` | `0` |
| `SMinus` | `SPlus` | `0` |

随伴は係数 `1` の一項、内積・trace は虚部が零の `Complex<T>` を返す。trace を次元で割らない。

## 5. 基底

### 5.1. 型・順序・展開

| 基底型 | core の実装 | 順序付き要素 |
| --- | --- | --- |
| `SzBasis` | `Basis<SpinHalfState, T>`・`StateBasis<SpinHalf, T>` | `Up, Down` |
| `CartesianBasis` | `Basis<SpinHalfOperator, T>`・`OperatorBasis<SpinHalf, T>` | `I, Sx, Sy, Sz` |
| `LadderBasis` | `Basis<SpinHalfOperator, T>`・`OperatorBasis<SpinHalf, T>` | `I, SPlus, SMinus, Sz` |

状態基底では両状態をそれぞれ係数 `1` の自身へ展開する。演算子基底では次の規則とし、二項は記載順で返す。

| primitive | Cartesian 基底への展開 | Ladder 基底への展開 |
| --- | --- | --- |
| `I` | `I` | `I` |
| `Sx` | `Sx` | `SPlus/2 + SMinus/2` |
| `Sy` | `Sy` | `−i SPlus/2 + i SMinus/2` |
| `Sz` | `Sz` | `Sz` |
| `SPlus` | `Sx + i Sy` | `SPlus` |
| `SMinus` | `Sx − i Sy` | `SMinus` |

基底要素自身は係数 `1` の一項として返す。どの展開結果も選択した基底の要素だけを含む。

### 5.2. 順序番号と列挙

`index()` は上表の順序に従う零始まりの番号を返す。Cartesian 基底の `SPlus`・`SMinus`、Ladder 基底の `Sx`・`Sy` は `None` とする。`SzBasis` の対象 primitive はすべて基底要素である。

`elements()` は固定配列の所有 iterator を返す。`Vec` や `Box` による列挙用の確保は行わない。数値型を推論できない箇所は `<CartesianBasis as Basis<SpinHalfOperator, f64>>::elements(&basis)` のように明示する。

要素数は core の `state_basis_len::<SpinHalf, T>()`・`operator_basis_len::<SpinHalf, T>()` からそれぞれ `2`・`4` と導出する。基底ごとの重複した要素数 API は追加しない。

## 6. 実装順序と完了条件

| 順序 | 実装対象 | 完了条件 |
| --- | --- | --- |
| 1 | crate の骨組み・公開型・係数生成 | workspace 継承と `no_std` が成立し、数値精度を持たない primitive を定義できる |
| 2 | 内積・随伴・trace・恒等元 | 各規則の全入力を確認できる |
| 3 | 状態への作用 | 全12組が既定の行列表現と一致する |
| 4 | 局所積 | 全36組が返却表と行列乗算の両方に一致する |
| 5 | 三種類の基底 | 展開・順序番号・列挙・要素数が整合する |
| 6 | 公開 API の結合テストとドキュメント | core の契約だけで全能力を利用でき、単精度・倍精度と最小数値制約を確認できる |

今回の計画では、実装開始前に追加の選択を必要とする未決事項はない。実装で型制約の矛盾などが見つかった場合は、実装とこのメモを更新する。

## 7. 検証計画

### 7.1. 数値参照と代数法則

- `f32`・`f64` の両方で、要件に定める `Up, Down` の二次元ベクトルと六種類の `2×2` 行列を独立に定義し、局所規則から得た結果と照合する。浮動小数点の比較では精度ごとの `atol`・`rtol` を明示する。
- 積36組、作用12組、内積4組、随伴6種、trace6種、および恒等元を網羅する。数学的な照合に加え、積の表で選んだ primitive・項数・項順序も構造として確認する。
- 全演算子の三つ組について積の結合則、全演算子二つと状態について `(AB)ψ = A(Bψ)`、左右の恒等元、恒等作用を確認する。線型結合への拡張はテスト内の補助処理で行う。
- 随伴の二重適用、積の順序反転、内積と随伴の整合、trace の巡回性を確認する。
- `[Sx,Sy] = i Sz` と循環する交換関係、`[SPlus,SMinus] = 2 Sz`、`[Sz,SPlus] = SPlus`、`[Sz,SMinus] = −SMinus` を確認する。
- `Sx² = Sy² = Sz² = I/4`、`SPlus² = SMinus² = 0`、`SPlus = Sx+iSy`・`SMinus = Sx−iSy` の線型関係の保存を確認する。

### 7.2. 基底・型制約・公開 API

- 各基底の列挙順、番号、要素数、自己展開、非基底 primitive の `None`、展開先の要素の所属を確認する。
- 全 primitive の基底展開を行列と照合し、両演算子基底を往復した結果が直接展開と同じ座標になることを確認する。`SPlus + SMinus = 2 Sx` もテスト内で座標を結合して確認する。
- 数値型に `Clone`・`Copy`・`Num` を実装しないテスト用の型を用意し、指定した演算能力だけで局所規則と基底が利用できることを確認する。内部で必要な分数を表現できる型を使い、新たな通常依存は追加しない。
- `SpinHalf` の primitive 型と次元が `f32`・`f64` で共通であることを確認する。
- integration test から公開 API だけを利用し、直接規則を呼ぶ例と、型推論に必要な返却型指定・完全修飾呼び出しを doc test に含める。

### 7.3. 実行するチェック

既存の direnv 環境を使い、次を実行する。

```sh
direnv exec . cargo test -p qmb-spin-half
direnv exec . cargo test --workspace
direnv exec . cargo fmt --all -- --check
direnv exec . cargo clippy --workspace --all-targets -- -D warnings
direnv exec . env RUSTDOCFLAGS='-D warnings' cargo doc -p qmb-spin-half --no-deps
direnv exec . cargo tree -p qmb-spin-half -e features
```

依存 feature に通常依存の `std`・`libm` が紛れ込んでいないことも確認する。テスト側の `std` 利用のために、通常依存へ `std` feature を追加しない。

検証には direnv が提供する既存の Rust 環境を使い、使用したバージョンを記録する。追加 toolchain の導入は前提としない。今回のメモ作成では Rust コードや Cargo ファイルを変更せず、この crate のテストを実行済みとはしない。
