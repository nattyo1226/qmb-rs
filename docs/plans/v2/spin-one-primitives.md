# Spin-one の局所 primitive 実装案

作成日: 2026-10-29
更新日: 2026-10-29
状態: primitive・基底・列挙順、局所積の返却表現、std / Float、係数の除算と match による局所規則の実装方針を決定。記載する crate・型・追加 API は未実装。

## 1. 対象と参照した実装

spin-1 の三状態、全九次元演算子空間、局所規則と完全な基底を提供する案を記録する。局所規則は core の既存契約で実装できる。Tensor 合成には[実装済みの共通評価機構](../../design/composite-expression.md#12-ordinary--graded-の共通評価と移行)を使う。

| 参照 | 本案に適用する内容 |
| --- | --- |
| [局所構造の設計](../../design/local-structure.md)・[primitive / 式の設計](../../design/primitive-expression.md) | S・O の直接指定、同じ型の有限線型結合への閉性、Identity・Dimension、個別能力の bound |
| [ops.rs](../../../crates/qmb-core/src/ops.rs)・[basis.rs](../../../crates/qmb-core/src/basis.rs) | Act の関連型、Product / Adjoint の返却型、完全基底の三能力 |
| [spin-half の operator.rs](../../../crates/qmb-spin-half/src/operator.rs)・[state.rs](../../../crates/qmb-spin-half/src/state.rs) | `hbar=1`、Cartesian / ladder の関係、数値精度を持たない primitive |
| [spin-half の basis.rs](../../../crates/qmb-spin-half/src/basis.rs)・[coefficient.rs](../../../crates/qmb-spin-half/src/coefficient.rs)・[lib.rs](../../../crates/qmb-spin-half/src/lib.rs) | 局所基底 marker、定数の生成、`QuantumSystem` の入口 |
| [テスト用三準位代数](../../../crates/qmb-core/tests/support/three_level.rs) | 行列単位と Identity の共存、delta による積・作用、基底番号、独立した行列参照 |
| [要件](../../requirements.md) 2.5、6 節・[基底の設計](../../design/basis.md) | 差し替え検証と公開 builtin の区別、全演算子基底、座標と行列の区別 |

既存の三準位テスト代数は spin-1 の公開実装ではない。特に spin 演算子・√2 の係数・spin の交換則を提供していないので、その存在をもって今回の機能を実装済みとは扱わない。

## 2. 規約と閉性の問題

`hbar=1`、状態の順序を `(|+1>, |0>, |-1>)` とする。正規直交内積と以下の ladder 規約を採る。

```math
S_z|m\rangle=m|m\rangle,\qquad
S_\pm|m\rangle=\sqrt{2-m(m\pm1)}|m\pm1\rangle,
\qquad S_\pm=S_x\pm iS_y.
```

端の状態への ladder 作用は zero とする。`r=√2` と書くと、独立な参照行列は次である。

```math
S_z=\begin{pmatrix}1&0&0\\0&0&0\\0&0&-1\end{pmatrix},\quad
S_+=r\begin{pmatrix}0&1&0\\0&0&1\\0&0&0\end{pmatrix},\quad S_-=S_+^\dagger,
```

```math
S_x=\frac{r}{2}\begin{pmatrix}0&1&0\\1&0&1\\0&1&0\end{pmatrix},\quad
S_y=\frac{r}{2}\begin{pmatrix}0&-i&0\\i&0&-i\\0&i&0\end{pmatrix}.
```

spin-half の6 variant `I,Sx,Sy,Sz,SPlus,SMinus` だけを移植しても、その span は高々4次元である。spin-1 では、例えば `Sz²=diag(1,0,1)` はその span に入らず、`SPlus²` も非零となる。全9次元の演算子空間と Product の閉性を満たすために、追加の primitive が必要である。

## 3. 採用する primitive と基底の方針

`qmb-spin-one` を `std` の利用を許可し、`unsafe` 禁止の具体的代数 crate として追加する。no_std 対応は初期要件に含めない。core の既存契約・数値型の汎用性は維持する。

```rust
enum SpinOneState { Plus, Zero, Minus }

enum SpinOneOperator {
    I, Sx, Sy, Sz, SPlus, SMinus,
    Dxy, Q0, Qxy, Qyz, Qzx,
    Q1Plus, Q1Minus, Q2Plus, Q2Minus,
}
```

公開 variant 名は上の列挙で確定する。二乗成分の差 `Sx²−Sy²` は `Dxy` とする。各演算子の定義は以下のとおりである。`SpinOneState::Zero` は磁気量子数0の正規化状態であり、zero vector ではない。Zero operator は共通の空の線型結合で表す。

### 3.1. 演算子 primitive

恒等演算子1個、スピン5個、四重極9個の計15値を用意する。反交換子を `{A,B}=AB+BA` と書く。

| Primitive | 定義・意味 |
| --- | --- |
| `I` | 恒等演算子 |
| `Sx`, `Sy`, `Sz` | スピンの Cartesian 3成分 |
| `SPlus`, `SMinus` | `S±=Sx±iSy`。磁気量子数を一つ変える昇降演算子 |
| `Dxy` | `Sx²−Sy²` |
| `Q0` | `3Sz²−2I` |
| `Qxy`, `Qyz`, `Qzx` | それぞれ `{Sx,Sy}`, `{Sy,Sz}`, `{Sz,Sx}` |
| `Q1Plus`, `Q1Minus` | `Q1±={Sz,S±}=Qzx±iQyz` |
| `Q2Plus`, `Q2Minus` | `Q2±=S±²=Dxy±iQxy` |

四重極の Cartesian 5値と昇降型4値は線型従属であり、四重極の span は5次元である。これらは[双極子・四重極による spin-one の構成（式3–7）](https://arxiv.org/pdf/1305.2560)と[球面テンソルによる多重極展開（式29–30）](https://arxiv.org/pdf/1409.1106)に対応する。今回は上表の係数・位相を使い、正規直交な球面テンソルの規格化は付けない。

状態・演算子は `Clone, Copy, Debug, PartialEq, Eq, Hash` を実装する。構造的な等値性は variant を区別する。`SpinOneState::dim()` は3、`SpinOneOperator::dim()` は9とし、15という値の個数と次元を混同しない。`Identity::identity()` は I を返し、数値能力を要求しない。

### 3.2. 初期提供する基底

状態基底 `SzBasis` は `[Plus,Zero,Minus]` とする。演算子基底は `CartesianBasis` と `LadderBasis` の2種類を提供する。spin-half と名称を揃え、各基底の列挙順は次の9要素の順で確定する。

| 基底 | 列挙順 | 用途 |
| --- | --- | --- |
| `CartesianBasis` | `[I,Sx,Sy,Sz,Dxy,Q0,Qxy,Qyz,Qzx]` | 全要素が Hermite。磁化・単イオン異方性・スピンネマティック秩序の成分を直接扱う |
| `LadderBasis` | `[I,SPlus,SMinus,Sz,Q0,Q1Plus,Q1Minus,Q2Plus,Q2Minus]` | 球面テンソルに対応する昇降型。磁気量子数の変化と選択則を明示する |

`LadderBasis` の各要素 B は `[Sz,B]=qB` を満たす。q は I・Sz・Q0 では0、S±・Q1± では±1、Q2± では±2である。回転に対して I は rank 0、スピンは rank 1、四重極は rank 2の成分をなす。ただし上の要素には標準的な球面テンソルの規格化・位相係数を付けていない。

どちらも全演算子空間に完全かつ線型独立で、Hilbert–Schmidt 内積 `Tr(A†B)` に関して直交する。正規化は行わず、CartesianBasis のノルム二乗は I が3、Q0 が6、他が2、LadderBasis では I が3、Sz が2、Q0 が6、他が4となる。基底を正規直交化するための √3・√6 の係数能力は追加しない。

両基底の和集合が公開 primitive の15値となる。`elements()` は選んだ基底の9値だけを列挙し、`index` はその順序番号を返す。基底に含まれない primitive には `None` を返す。この番号は export 行列の格納順序の指定ではない。

`SpinOne<C>` を `QuantumSystem` に接続し、共通の `LocalState` / `LocalOperator` を構築する。primitive や marker に C の bound を一律に要求しない。

### 3.3. 表現の選択理由

| 案 | 判断 |
| --- | --- |
| 6種の spin primitive のみ | 積に閉じず、全演算子空間を張らないため不採用 |
| I と9個の E のみ。Sx 等は式 constructor | 最小の局所規則になるが、既存 spin-half と異なり spin 演算子自体を primitive として保持できない |
| 6種の自然な primitive と9個の E | 旧推奨案。行列単位は規則の導出・検証用に参照し、初期の公開 primitive・基底には採用しない |
| I・spin・Cartesian 四重極の独立な9要素のみ | 完全基底になる。昇降型の要素も primitive として列挙するため、両基底の和集合へ拡張する |
| I・spin 5値・四重極9値の計15値 | 採用。Cartesian / ladder 両基底の要素を直接保持でき、物理的な入口と積の閉性を両立する |

Gell-Mann 基底、公開の MatrixUnitBasis、正規直交化した球面テンソル基底は後続候補とする。今回の2基底は選択を利用者が明示し、通常の局所評価には基底を要求しない。

## 4. 全 primitive を覆う局所規則の実装方針

### 4.1. 参照用の行列単位分解

基本順序を `+,0,−`、`E_ab=|a><b|` とし、以下を局所規則の導出・検証用の分解として記録する。実行時の局所規則は4.2節の match で直接記述する。

```math
\begin{aligned}
I &= E_{++}+E_{00}+E_{--}, & S_z &= E_{++}-E_{--},\\
S_+ &= r(E_{+0}+E_{0-}), & S_- &= r(E_{0+}+E_{-0}),\\
S_x &= \tfrac r2(E_{+0}+E_{0+}+E_{0-}+E_{-0}),\\
S_y &= \tfrac r2(-iE_{+0}+iE_{0+}-iE_{0-}+iE_{-0}).
\end{aligned}
```

四重極の分解は以下となる。

```math
\begin{aligned}
Q_0&=E_{++}-2E_{00}+E_{--},\\
D_{xy}&=E_{+-}+E_{-+}, & Q_{xy}&=-iE_{+-}+iE_{-+},\\
Q_{yz}&=\tfrac r2(-iE_{+0}+iE_{0+}+iE_{0-}-iE_{-0}),\\
Q_{zx}&=\tfrac r2(E_{+0}+E_{0+}-E_{0-}-E_{-0}),\\
Q_{1+}&=r(E_{+0}-E_{0-}), & Q_{1-}&=r(E_{0+}-E_{-0}),\\
Q_{2+}&=2E_{+-}, & Q_{2-}&=2E_{-+}.
\end{aligned}
```

全15値の参照行列をこの分解と照合できる。実行時には密3×3行列や backend を生成しない。

### 4.2. Product

`Product<C, Term=SpinOneOperator>` の規則の導出には、次の行列単位の積を参照できる。

```math
E_{ab}E_{cd}=\delta_{bc}E_{ad}.
```

局所積の返却表現は、入力 primitive が共通して含まれる基底の中で閉じる方針に決定する。次の表の順に判定する。

| 入力の組 | 返却する primitive |
| --- | --- |
| I との積 | 相手をそのまま一項で返す |
| 両方が共通要素 I・Sz・Q0 | 共通要素だけで返す。返却順は `[I,Sz,Q0]` |
| 両方が CartesianBasis の要素 | CartesianBasis の要素で返す |
| 両方が LadderBasis の要素 | LadderBasis の要素で返す |
| Cartesian 専用要素と Ladder 専用要素の混在 | 固定の規則として CartesianBasis の要素で返す |

共通要素は対角演算子の完全基底であり、積も共通部分に閉じる。例えば `Sz²=(2I+Q0)/3`、`Sz Q0=Q0 Sz=Sz`、`Q0²=2I−Q0` となる。Cartesian の例は `Sx Sy=(Qxy+iSz)/2`、ladder の例は `SPlus²=Q2Plus` である。

入力 primitive の値だけからこの規則を判定できるため、式に選択基底を保持させたり、Product に基底引数を追加したりしない。各基底内の線型結合同士の積も線型性により同じ基底の中に閉じる。混在時や異なる積の括り方で、数学的に同じ結果の式構造が一致することは要求しない。

実装は既存の spin-half / fermion と同様に `match (*self, *rhs)` で積の規則を直接記述する。各分岐で、上の返却規則に従った `LinearCombination<SpinOneOperator,C>` を生成する。同じ結果になる組は or pattern 等でまとめられるが、数学的に zero になる組も含めて全225組を覆う。返却順は対象基底の列挙順、項数は最大9個とする。

規則に現れる1/2・1/3・1/6等は C の One・加算・除算で直接生成し、虚係数は ImaginaryUnit・乗算・符号反転で生成する。必要に応じて小さな private helper を共有する。独自の有理数型や記号係数の演算機構、実行時の行列単位への分解・逆変換は導入しない。数学的に zero の項は規則に含めず、zero operator は空の線型結合を返す。浮動小数点の相殺結果を調べて項を削除する処理は行わない。

確認すべき例は以下である。

```math
S_z^2=\tfrac23 I+\tfrac13 Q_0,\quad
S_+^2=Q_{2+},\quad S_-^2=Q_{2-},\quad S_+^3=S_-^3=0,
```

```math
S_+S_-=\tfrac43 I-\tfrac13 Q_0+S_z,\qquad
S_-S_+=\tfrac43 I-\tfrac13 Q_0-S_z.
```

spin-half の `Sα²=I/4`、`SPlus²=0` という分岐を流用しない。

### 4.3. Act と内積

標準の作用実装は `Act<Complex<R>>`、`R: num_traits::Float` とし、`Operator=SpinOneOperator`、`Term=SpinOneState` とする。作用規則も演算子と状態の組に対する match で直接記述し、名前付き spin に対して次の表と一致させる。I は状態を保存する。行列単位による導出では `E_ab|c>=δ_bc|a>` を参照できる。

| 演算子 | Plus | Zero | Minus |
| --- | --- | --- | --- |
| Sz | Plus | 0 | −Minus |
| SPlus | 0 | r Plus | r Zero |
| SMinus | r Zero | r Minus | 0 |
| Sx | (r/2) Zero | (r/2)(Plus + Minus) | (r/2) Zero |
| Sy | i(r/2) Zero | −i(r/2) Plus + i(r/2) Minus | −i(r/2) Zero |

表中の `0` は zero vector、`Zero` は m=0 の状態である。InnerProduct は3状態の Kronecker delta とし、式の半双線型性と係数の共役は core に委ねる。

### 4.4. Adjoint と通常 Trace

Adjoint は `SPlus†=SMinus`、`Q1Plus†=Q1Minus`、`Q2Plus†=Q2Minus` とし、I・Sx・Sy・Sz・Cartesian 四重極5値は自身とする。今回の昇降型の定義では随伴に追加の符号を付けない。`Adjoint<C, Term=SpinOneOperator>` は primitive を保存し、余分な基底展開をしない。

Trace は `Tr(I)=3`、他の全14値は0を直接返す。正規化や supertrace を採用しない。

## 5. √2 の生成と係数能力

現行 [coefficient.rs](../../../crates/qmb-core/src/coefficient.rs) には ImaginaryUnit 等はあるが、√2 を生成する能力はない。spin-half の1/2だけの係数 helper では、この正規直交基底上の spin-1 を表せない。基底の重みを変えて√2を消すと内積規約も変わるため、本案では行わない。

状態への作用は `Act<Complex<R>> where R: num_traits::Float` として実装する方針に決定する。Float は複素係数 C 全体ではなく、その実数成分 R に付ける。標準の利用対象は `Complex<f32>` / `Complex<f64>` とし、√2 は R の精度で次のように生成する。

```rust
fn sqrt_two<R: num_traits::Float>() -> R {
    (R::one() + R::one()).sqrt()
}
```

この helper は private とする。旧案の公開 trait `SqrtTwo` は追加しない。実係数の複素数への埋め込みは `Complex::new(value, R::zero())` とし、定数を f64 で計算して R へキャストする実装は採らない。平方根の結果は R の精度に従う近似値である。

num-traits 0.2.19 の Float は Num・Copy・NumCast・PartialOrd・Neg を含むため、R の Copy / Clone も許容する。Float を満たさない厳密数・非 Copy の実数型や、Complex<R> 以外の独自係数型への作用は、初期の標準実装の対象外とする。一方、Product・基底変換・primitive の内積・随伴・trace は引き続き汎用の C と操作ごとの小さい bound を維持する。構築・保持には Float を要求しない。

| 能力 | 想定する C の能力 |
| --- | --- |
| Identity・Dimension・QuantumSystem、primitive 保持 | なし |
| InnerProduct・Trace | `Zero + One`。3は `one + one + one` |
| Adjoint・SzBasis | `One` |
| CartesianBasis | `One + Neg + ImaginaryUnit` |
| LadderBasis | `One + Add + Div + Neg + ImaginaryUnit` |
| Product | `One + Add + Div + Neg + ImaginaryUnit` |
| Act | `C=Complex<R>`、`R: num_traits::Float` |

`One` が含む乗算能力も利用する。係数は非零の係数による除算に閉じ、各規則で必要な1/2・1/3・iを数学的な意味どおり表せることを要求する。Rust の `Div<Output=C>` だけでは整数除算も通るため、この数学的な契約も併記する。浮動小数点ではその精度に応じた近似を許容する。局所 Product は match の各分岐で定数係数を直接生成するため、C の Clone / Copy / Float / Num / 型変換を追加要求しない。共通式の評価・座標集約には既存 core の別の bound が必要である。

例えば1/3は `C::one() / (C::one() + C::one() + C::one())` で作る。分子・分母を別の内部型に保持する必要はない。

上表の Add・Div・Neg はそれぞれ `Output=C` とする。`One` は `num_traits::One`、`Zero` は `num_traits::Zero` を指す。空の線型結合で zero を返せるため、Product 自体に Zero は要求しない。Act では Float が必要な算術・Zero・One を含む。構築・保持には演算能力を要求しない。

基底変換は1/2と i で書けるため、両演算子基底の expand に Float や平方根能力を要求しない。全15値は Cartesian 基底の有理数・i 係数の線型結合であり、局所積の係数も有理数・i の範囲に閉じる。この係数で積の規則を直接記述し、Product に Float や平方根能力を要求しない。Act の遷移振幅に必要な √2 のみ Float で生成する。

規則上 zero になる項は match の結果に含めず、返却係数の生成には C の算術を使う。浮動小数点の近似零判定で項を削除する設計にはしない。

現在の `Basis<P,C>` は expand と metadata を一つの impl に持つため、各基底の `index` / `elements` もその impl の係数 bound のもとで利用する。metadata だけを理由に core の Basis 契約を変更しない。全 variant を受け取る Product / Act は、呼び出した値によらず同じ impl の bound を要求する。

### 5.1. 既存の局所式 API で追加される係数能力

次の表は core の既存 API 自体が要求する bound である。上表の primitive 規則・利用する Basis の bound と合わせて適用する。

| 局所式の操作 | core が C に要求する能力 |
| --- | --- |
| 保持・和・遅延した積 / 作用 / 随伴の構築 | なし。scaled の入力は `Into<C>` |
| 状態の `evaluate` | `Clone + Conjugate<Output=C>`。primitive の Act / Adjoint も必要 |
| 内積・期待値・trace・基底展開・座標取得 | `Clone + Zero + Mul<Output=C> + Conjugate<Output=C>`。各操作で必要な primitive / Basis の能力も必要 |
| 状態の `norm_squared` | 内積の bound に `RealPart` を追加 |
| 状態の `norm` | norm_squared の bound と `C::Real: Float` |
| 状態の `normalized` | norm の bound に `From<C::Real>` を追加 |

Zero は加算能力、One は乗算能力を含む。式の共有・係数の再利用のため core の評価には Clone が必要になる。汎用 C 向けの primitive 規則には Clone を追加せず、Act の実数成分 R は Float の Copy / Clone を利用できる。Conjugate は複素共役であり、数値型が実数だけであることを意味しない。f32 / f64 自体は虚数単位を表せず、全 variant を扱う Product / Act の係数には使えない。

以上は core の局所式 API の bound であり、変更しない。spin-one の状態式の evaluate・内積等には Act も必要なため、標準実装では `C=Complex<R>`・`R: Float` の制約も適用される。primitive の内積だけを直接呼ぶ場合は上の小さい bound のままである。複合式の共通評価や正規化 API には別途の bound があり、例えば norm の平方根・状態の正規化は qmb-normalize の実数成分の Float 能力を使う。

## 6. 基底と既存式への接続

`SzBasis::expand` は状態をそのまま返す。`CartesianBasis::expand` / `LadderBasis::expand` は、選択した基底に含まれる primitive を自身一項へ、その他を以下の相互変換へ展開する。

```math
\begin{aligned}
S_\pm&=S_x\pm iS_y,&
S_x&=\tfrac12(S_++S_-),& S_y&=\tfrac1{2i}(S_+-S_-),\\
Q_{1\pm}&=Q_{zx}\pm iQ_{yz},&
Q_{zx}&=\tfrac12(Q_{1+}+Q_{1-}),& Q_{yz}&=\tfrac1{2i}(Q_{1+}-Q_{1-}),\\
Q_{2\pm}&=D_{xy}\pm iQ_{xy},&
D_{xy}&=\tfrac12(Q_{2+}+Q_{2-}),& Q_{xy}&=\tfrac1{2i}(Q_{2+}-Q_{2-}).
\end{aligned}
```

I・Sz・Q0 は両基底で共有する。どちらの基底でも全15 primitive を9基底要素の線型結合に展開できる。`{I,Sx,Sy,Sz}` や `{I,SPlus,SMinus,Sz}` の4個だけを spin-1 の完全基底として提供しない。通常の評価に利用者の選択基底を要求しない。

既存 `LocalState` / `LocalOperator` の評価・期待値・随伴・trace・展開・座標取得・正規化を利用する。演算子の九つの基底係数と、状態基底に対する3×3作用行列を区別する。作用行列は既存 export の列ごとの作用で生成し、primitive 側に backend 固有の行列化を実装しない。

既存の状態・演算子型に配置し、ordinary / graded 共通の複合評価 API に接続する。Core の `HomogeneousParity`・`ParityOperator` として、全状態・演算子を even、P=I という trivial grading を明示実装し、磁気量子数 m の整数としての偶奇を fermion parity と解釈しない。spin-half・fermion と混在する型消去や一般 spin-S は本書の対象にしない。

## 7. Crate 構成と実装順序

想定ファイルは `src/{lib,state,operator,basis,coefficient}.rs`、検証は `tests/{contracts,minimal_scalar,local_expressions}.rs` とする。本番依存候補は既存 workspace の `qmb-core`・`num-traits`・`num-complex` である。num-complex は Act の係数型を表すため、本番依存とする。num-traits はこの crate 側で std feature を有効にし、qmb-normalize と同様に Float を利用する。workspace 全体の既定 feature や core の公開契約は変更しない。geometry・数値 backend は必要としない。

想定する依存宣言は以下とする。

```toml
num-traits = { workspace = true, features = ["std"] }
num-complex.workspace = true
```

1. 規約・公開型、Dimension・Identity・QuantumSystem を追加する。
2. C の算術による定数係数 helper と Float による private な √2 helper を実装する。
3. Product・Act の規則を match で直接記述し、InnerProduct・Adjoint・Trace、SzBasis・CartesianBasis・LadderBasis を実装する。
4. 局所式・ordinary な複合式・既存 export / normalize の接続を検証し、利用例と公開 API 文書を追加する。

標準の spin-half の型・数値 bound を変更しない。Grading trait は core の実装済み契約へ接続し、共通機構を再実装しない。専用 Heisenberg / bilinear-biquadratic モデル utility は primitive の動作確認例とは分けて設計する。

## 8. 受入条件と検証計画

以下は実装時に満たす条件であり、現時点の Rust 実装のテスト結果ではない。

- 全15²=225組の積、15×3=45組の作用、15種の随伴・trace、9組の状態内積を、2節のスピン行列と3.1節の定義から独立に構成した四重極行列・状態ベクトルに照合する。参照行列を製品の match 規則や係数 helper から生成しない。Product の出力が4.2節で指定した基底・共通部分に閉じ、I との積が相手を保持することも確認する。
- 全15³=3,375組の三重積の結合則、675組の `(AB)ψ=A(Bψ)`、単位元、二回の随伴、随伴の反準同型性と内積整合性、trace の巡回性を確認する。
- `[Sx,Sy]=iSz` と巡回置換、`[Sz,S±]=±S±`、`[S+,S−]=2Sz`、`Sx²+Sy²+Sz²=2I` を検証する。spin-half の Casimir 3I/4 と混同しない。
- `Sz²=(2I+Q0)/3`、`SPlus²=Q2Plus`、`SMinus²=Q2Minus`、ladder の三乗 zero、`Tr(Sα²)=2`、端状態への zero 作用を確認する。
- 状態基底3要素・両演算子基底それぞれ9要素の列挙 / 番号 / 完全性を確認する。全15 primitive の直接展開と局所式の再展開、両基底間の往復、展開前後の作用・trace・座標を検証する。両基底の Hilbert–Schmidt Gram 行列が3.2節のノルムを持つ対角行列となることも確認する。
- `Complex<f32>` / `Complex<f64>` で局所式の作用・期待値・norm・正規化・3×3行列 export を照合する。ordinary 二因子では9状態・81演算子基底を持ち、積・作用・trace・mixed-radix 順序が既存 core と整合することも確認する。
- Vec・ndarray・sprs は同じ状態基底行列を返す。特に Sy の虚数符号、列ベクトル規約、m の順序を固定し、転置や逆順でも偶然通る対角行列だけを検証に使わない。
- Float・Clone / Copy のない独自係数で Product・両基底の変換・primitive の内積 / 随伴 / trace を呼べること、Identity / Dimension / 構築に数値 bound を要求しないことを capability テストにする。Act は `Complex<f32>` / `Complex<f64>` で全45組を照合し、各 R の精度で √2 を生成することを確認する。Act に非 Float・非 Copy の係数を要求する旧案のテストは初期の受入条件に含めない。
- 数学的に zero の積・作用が空の線型結合を返し、match 規則で必要な項だけを生成することを確認する。係数生成は C の算術を使い、近似零判定で項を削除しない。数値比較は既存要件の精度別の許容誤差を使う。

実装時の検証入口は `cargo test -p qmb-spin-one`、workspace 登録後の `cargo test --workspace --all-features`、`cargo check -p qmb-spin-one`、`cargo fmt --all -- --check` とする。no_std target への対応は初期要件に含めない。MSRV の実ビルド・性能計測は実施有無を別に記録する。

## 9. 確定方針と実装への移行

確定した方針は、正規直交な `Plus,Zero,Minus`、I・スピン5値・四重極9値の計15演算子 primitive、SzBasis と各9要素の CartesianBasis / LadderBasis、演算子基底を正規直交化しない規約、同じ基底の要素同士の積をその基底内で返す規則である。両演算子基底の和集合を primitive とし、公開の行列単位は初期提供しない。局所積の共通要素・混在時の扱いは4.2節のとおりとする。

std の利用を許可し、Act は `Complex<R>`・`R: Float` 向けに実装する。Product・基底変換等は操作ごとの小さい bound を維持し、構築・保持には数値 bound を付けない。√2 は R の sqrt で生成し、SqrtTwo の公開 trait は追加しない。

variant 名と基底の列挙順も3節の内容で確定する。係数は非零による除算に閉じることを要求し、定数係数は C の算術で生成する。局所規則は既存実装と同様に match で直接記述する。独自の有理数型・記号係数や実行時の行列単位分解を導入する旧案は採用しない。private helper の分け方や共通する match 分岐のまとめ方は、公開 API・数学的規則・基底内閉性を保ったうえで実装時に決める。

まず全225積と45作用を覆う規則・両完全基底・係数能力の検証を実装し、その後に追加基底・モデル・一般 spin-S を検討する。tensor 積の計画とは独立に完成を判定できる。

## 10. 本ドキュメント作成時の確認

初稿では現行の設計文書・core 契約・spin-half・三準位テスト代数を照合した。初稿の補助的な Python 標準ライブラリの複素行列計算は旧案（名前付き6値と行列単位9値）を対象としたものであり、今回採用した四重極を含む15値の検証結果としては扱わない。

今回、対話で合意した primitive と2種類の演算子基底の方針を3節へ記録し、局所規則・係数・検証計画もそれに合わせて更新した。変更は計画文書のみであり、Rust の spin-one 実装・crate 追加・ビルド・テストは行っていない。8節の受入条件は実装後に検証する。

今回の補助的な Python 標準ライブラリの複素3×3行列計算では、両基底の Gram 行列、各基底への全15値の展開、基底間変換、Sz²・ladder の積の例、対角 E の逆変換、q の選択則、全225積が Cartesian 基底の span に入ることを許容誤差1e-12で確認した。相対ファイルリンク・コードフェンス・差分の空白も確認した。この計算は Rust の局所規則や係数能力、f32 精度の検証を代替しない。

局所積の基底内閉性の方針決定時には、有理数・i・√2 の厳密な係数を使った Python 標準ライブラリの補助計算を行った。両基底の Gram 行列と全15値の展開に加え、全225組を共通要素9組・Cartesian の残り72組・ladder の残り72組・混在72組に分け、指定した基底への再構成と積の係数から √2 が消えることを確認した。これは Product に平方根能力が不要であるという実装案の根拠であり、Rust の capability テストは未実施である。

std / Float の方針決定では、利用中の num-traits 0.2.19 の Float の定義と std feature、num-complex 0.4.6 の Complex の実装、qmb-normalize の依存宣言を確認した。上の数理的な確認結果は維持するが、Float を使う Act の Rust 実装・ビルド・精度別テストは今後実施する。

係数の除算と match の方針決定に合わせ、記号係数・実行時の行列単位分解を前提とした実装案を、C の算術で係数を生成する直接の局所規則へ更新した。4.1節の行列単位分解は、規則の導出・検証用の参照として残す。
