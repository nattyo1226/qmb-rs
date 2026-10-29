# Spin-one の局所 primitive 実装案

作成日: 2026-10-29
更新日: 2026-10-29
状態: 検討案。記載する crate・型・追加 API は未実装。

## 1. 対象と参照した実装

spin-1 の三状態、全九次元演算子空間、局所規則と完全な基底を提供する案を記録する。現在の ordinary な core の契約で実装でき、graded tensor の追加を前提にしない。tensor 合成そのものは [別の計画](tensor-products.md)に委ねる。

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

## 3. 提案する型と基底

`qmb-spin-one` を `no_std + alloc`、`unsafe` 禁止の具体的代数 crate として追加する。

```rust
enum SpinOneState { Plus, Zero, Minus }

enum SpinOneOperator {
    I, Sx, Sy, Sz, SPlus, SMinus,
    MatrixUnit { row: SpinOneState, column: SpinOneState },
}
```

`SpinOneState::Zero` は磁気量子数0の正規化状態であり、zero vector ではない。`MatrixUnit { row:r, column:c }` は `E_rc=|r><c|` を表す。有限 enum を行・列に用いるため、範囲外の量子数や行列 index を受け入れる経路は不要である。Zero operator は共通の空の線型結合で表す。

状態は3値、演算子は名前付き6値と行列単位9値の計15値を持ち、`Clone, Copy, Debug, PartialEq, Eq, Hash` を実装する。`SpinOneState::dim()` は3、`SpinOneOperator::dim()` は9とする。variant 数を次元と混同しない。`Identity::identity()` は I を返し、数値能力を要求しない。

状態基底 `SzBasis` は `[Plus,Zero,Minus]` とする。既定の演算子基底 `MatrixUnitBasis` は9個の E のみを行優先で列挙し、番号を `3*row_index + column_index` とする。これは局所演算子の基底番号であり、export 行列の格納順序の指定ではない。

名前付き spin primitive と I はこの基底に含まれず、`index` は `None` を返す。`elements()` は15値全部を返してはいけない。I が基底要素でなくても `Basis::expand` が対角3項へ展開すれば、既存 core の Identity / 基底契約を満たす。

`SpinOne<C>` を `QuantumSystem` に接続し、共通の `LocalState` / `LocalOperator` を構築する。primitive や marker に C の bound を一律に要求しない。

### 表現の選択肢

| 案 | 判断 |
| --- | --- |
| 6種の spin primitive のみ | 積に閉じず、全演算子空間を張らないため不採用 |
| I と9個の E のみ。Sx 等は式 constructor | 最小の局所規則になるが、既存 spin-half と異なり spin 演算子自体を primitive として保持できない |
| 6種の自然な primitive と9個の E | 推奨。自然な入口と閉性を両立し、完全基底を単純に定義できる |
| spin と quadrupole の独立な9要素 | 物理的な解釈はよいが、積・基底変換の規約が増える。追加基底として後から検討 |

## 4. 全 primitive を覆う局所規則

### 4.1. 行列単位への有限分解

基本順序を `+,0,−` とし、次を private な有限分解の規則とする。これは基底をユーザに選ばせる操作ではなく、局所 Product / Act を定義する固定の実装手段である。

```math
\begin{aligned}
I &= E_{++}+E_{00}+E_{--}, & S_z &= E_{++}-E_{--},\\
S_+ &= r(E_{+0}+E_{0-}), & S_- &= r(E_{0+}+E_{-0}),\\
S_x &= \tfrac r2(E_{+0}+E_{0+}+E_{0-}+E_{-0}),\\
S_y &= \tfrac r2(-iE_{+0}+iE_{0+}-iE_{0-}+iE_{-0}).
\end{aligned}
```

E 自身の分解は係数1の自身一項である。この規則で全15値を処理できる。局所規則は sparse な有限項を扱い、密3×3行列や backend を生成しない。

### 4.2. Product

`Product<C, Term=SpinOneOperator>` は次を使う。

```math
E_{ab}E_{cd}=\delta_{bc}E_{ad}.
```

I との積は相手の primitive をそのまま返す。それ以外の組は4.1節の有限項を掛け合わせ、接続する中間 index が一致する項だけを集約し、E の有限線型結合を返す。返却項は最大9個、列挙順は行優先に揃える。同じ値を表す結果を名前付き spin primitive へ毎回戻すことは要求しない。

固定の局所規則のため、実装上は係数を小さな整数・二分の有理数・i・r の記号で保持し、`r²=2`・`i²=−1` を用いて厳密に集約してから C に変換する案を推奨する。これにより、入力 C の Clone や浮動小数点での近似零判定を局所 Product に要求せずに済む。汎用記号計算エンジンは追加せず、この有限表の係数だけを扱う。

確認すべき例は以下である。

```math
S_z^2=E_{++}+E_{--},\quad
S_+^2=2E_{+-},\quad S_-^2=2E_{-+},\quad S_+^3=S_-^3=0,
```

```math
S_+S_-=2(E_{++}+E_{00}),\qquad
S_-S_+=2(E_{00}+E_{--}).
```

spin-half の `Sα²=I/4`、`SPlus²=0` という分岐を流用しない。

### 4.3. Act と内積

`Act<C>` は `Operator=SpinOneOperator`、`Term=SpinOneState` とする。`E_ab|c>=δ_bc|a>` と4.1節の分解を用い、名前付き spin に対して次の表と一致させる。I は状態を保存する。

| 演算子 | Plus | Zero | Minus |
| --- | --- | --- | --- |
| Sz | Plus | 0 | −Minus |
| SPlus | 0 | r Plus | r Zero |
| SMinus | r Zero | r Minus | 0 |
| Sx | (r/2) Zero | (r/2)(Plus + Minus) | (r/2) Zero |
| Sy | i(r/2) Zero | −i(r/2) Plus + i(r/2) Minus | −i(r/2) Zero |

表中の `0` は zero vector、`Zero` は m=0 の状態である。InnerProduct は3状態の Kronecker delta とし、式の半双線型性と係数の共役は core に委ねる。

### 4.4. Adjoint と通常 Trace

Adjoint は `E_ab†=E_ba`、`SPlus†=SMinus`、I・Sx・Sy・Sz は自身とする。`Adjoint<C, Term=SpinOneOperator>` は名前付き primitive を保存し、余分な基底展開をしない。

Trace は `Tr(E_ab)=δ_ab`、`Tr(I)=3`、`Tr(Sx)=Tr(Sy)=Tr(Sz)=Tr(SPlus)=Tr(SMinus)=0` を直接返す。正規化や supertrace を採用しない。

## 5. √2 の生成と係数能力

現行 [coefficient.rs](../../../crates/qmb-core/src/coefficient.rs) には ImaginaryUnit 等はあるが、√2 を生成する能力はない。spin-half の1/2だけの係数 helper では、この正規直交基底上の spin-1 を表せない。基底の重みを変えて√2を消すと内積規約も変わるため、本案では行わない。

具体的代数側に独立した能力 `SqrtTwo` を提案する。

```rust
trait SqrtTwo {
    fn sqrt_two() -> Self;
}
```

返す r は実数の正の√2に対応し、数学的には `r*r=2`、`conjugate(r)=r` を満たす。f32 / f64 には各精度の `core::f32::consts::SQRT_2` / `core::f64::consts::SQRT_2`、`Complex<R>` には実部 `R::sqrt_two()`・虚部 `R::zero()` を使う。実装は `R: SqrtTwo + Zero` でよく、Float や実行時の平方根を必要としない。

この trait はまず `qmb-spin-one` が所有し、一般の平方根・norm 用 API を core に追加する前提にはしない。独自の厳密係数もこの trait を実装できる。定数を f64 で計算して C へキャストする実装は採らない。

| 能力 | 想定する C の能力 |
| --- | --- |
| Identity・Dimension・QuantumSystem、primitive 保持 | なし |
| InnerProduct・Trace | `Zero + One`。3は `one + one + one` |
| Adjoint・SzBasis | `One` |
| Product・Act・MatrixUnitBasis | `One + Add + Div + Neg + ImaginaryUnit + SqrtTwo`。各算術の Output は C |

`One` が含む乗算能力も利用する。局所 Product は記号的な有限係数を C へ出力する方式なら、C の Clone / Copy / Float / Num / 型変換を追加要求しない。共通式の評価・座標集約には既存 core の別の bound が必要である。

現在の `Basis<P,C>` は expand と metadata を一つの impl に持つため、同じ MatrixUnitBasis の `index` / `elements` もその impl の係数 bound のもとで利用する。metadata だけを理由に core の Basis 契約を変更しない。全 variant を受け取る Product / Act も、呼び出した値が E でも同じ bound を要求する。この制約を避けたい場合は「I と E だけの別 primitive 型」が必要であり、今回の推奨 API には追加しない。

## 6. 基底と既存式への接続

`SzBasis::expand` は状態をそのまま返す。`MatrixUnitBasis::expand` は4.1節の規則を利用し、E だけを返す。Product の固定分解 helper を共有できるが、通常の評価に利用者の選択基底を要求しない。

E の9要素は Hilbert–Schmidt 内積で正規直交であり、完全性は行列単位として保証できる。`{I,Sx,Sy,Sz}` や `{I,SPlus,SMinus,Sz}` を spin-1 の完全な演算子基底として提供しない。将来 quadrupole 基底を追加するなら、その要素を primitive の値として列挙できる設計と変換規則を別途検討する。

既存 `LocalState` / `LocalOperator` の評価・期待値・随伴・trace・展開・座標取得・正規化を利用する。演算子の九つの基底係数と、状態基底に対する3×3作用行列を区別する。作用行列は既存 export の列ごとの作用で生成し、primitive 側に backend 固有の行列化を実装しない。

ordinary な複合系では既存の状態・演算子型にそのまま配置できる。graded の能力を後から付与する場合は、全状態・演算子を even、P=I という trivial grading とし、磁気量子数 m の整数としての偶奇を fermion parity と解釈しない。spin-half・fermion と混在する型消去や一般 spin-S は本書の対象にしない。

## 7. Crate 構成と実装順序

想定ファイルは `src/{lib,state,operator,basis,coefficient}.rs`、検証は `tests/{contracts,minimal_scalar,local_expressions}.rs` とする。本番依存候補は既存 workspace の `qmb-core`・`num-traits` と、SqrtTwo の Complex 実装を提供する `num-complex` である。既存 spin-half の num-complex は dev-dependency だが、本案では trait の実装先に使うので直接の本番依存となる。feature は default 無効を維持し、std・geometry・数値 backend を必要としない。

1. 規約・公開型、Dimension・Identity・QuantumSystem を追加する。
2. 係数の有限タグと SqrtTwo、全 variant の行列単位への分解を実装する。
3. Product・Act・InnerProduct・Adjoint・Trace、状態 / 演算子の完全基底を実装する。
4. 局所式・ordinary な複合式・既存 export / normalize の接続を検証し、利用例と公開 API 文書を追加する。

標準の spin-half の型・数値 bound を変更しない。core の新しい grading trait の完成をこの順序に挟まない。専用 Heisenberg / bilinear-biquadratic モデル utility は primitive の動作確認例とは分けて設計する。

## 8. 受入条件と検証計画

以下は実装時に満たす条件であり、現時点の Rust 実装のテスト結果ではない。

- 全15²=225組の積、15×3=45組の作用、15種の随伴・trace、9組の状態内積を、2節の行列と独立な E 行列・状態ベクトルに照合する。参照行列を製品用分解 helper から生成しない。
- 全15³=3,375組の三重積の結合則、675組の `(AB)ψ=A(Bψ)`、単位元、二回の随伴、随伴の反準同型性と内積整合性、trace の巡回性を確認する。
- `[Sx,Sy]=iSz` と巡回置換、`[Sz,S±]=±S±`、`[S+,S−]=2Sz`、`Sx²+Sy²+Sz²=2I` を検証する。spin-half の Casimir 3I/4 と混同しない。
- `Sz²`、`SPlus²=2E_+-`、`SMinus²=2E_-+`、ladder の三乗 zero、`Tr(Sα²)=2`、端状態への zero 作用を確認する。
- 状態基底3要素・演算子基底9要素の列挙 / 番号 / 完全性、全15 primitive の直接展開と局所式の再展開、展開前後の作用・trace・座標を検証する。
- `Complex<f32>` / `Complex<f64>` で局所式の作用・期待値・norm・正規化・3×3行列 export を照合する。ordinary 二因子では9状態・81演算子基底を持ち、積・作用・trace・mixed-radix 順序が既存 core と整合することも確認する。
- Vec・ndarray・sprs は同じ状態基底行列を返す。特に Sy の虚数符号、列ベクトル規約、m の順序を固定し、転置や逆順でも偶然通る対角行列だけを検証に使わない。
- 独自係数で√2を供給できること、Float・Clone / Copy のない係数で primitive 規則を呼べること、Identity / Dimension / 構築に√2や ImaginaryUnit を要求しないことを capability テストにする。
- 有限係数の集約は記号的に zero を判定し、C に変換した後の近似零で項を削除しない。数値比較は既存要件の精度別の許容誤差を使う。

実装時の検証入口は `cargo test -p qmb-spin-one`、workspace 登録後の `cargo test --workspace --all-features`、`cargo check -p qmb-spin-one --no-default-features`、`cargo fmt --all -- --check` とする。専用 no_std target・MSRV の実ビルド・性能計測は実施有無を別に記録する。

## 9. 実装前に固定する詳細

推奨する規約は、正規直交な `Plus,Zero,Minus`、6種の自然な演算子と9個の E、全9要素の MatrixUnitBasis、係数能力 SqrtTwo、固定の有限分解による局所規則である。

残る詳細は、state / variant / 基底名、有限係数の内部表現、SqrtTwo を crate root からどう案内するかである。まず全225積と45作用を覆う規則・完全基底・最小係数の検証を実装し、その後に追加の物理的基底・モデル・一般 spin-S を検討する。tensor 積の計画とは独立に完成を判定できる。

## 10. 本ドキュメント作成時の確認

現行の設計文書・core 契約・spin-half・三準位テスト代数を照合した。補助的な Python 標準ライブラリの複素行列計算では、225組の積と45組の作用、15値の随伴・trace、spin の交換則、Casimir、ladder の二乗・三乗を確認した。有限分解による計算と直接記述した参照行列を比較し、許容誤差は1e-12とした。相対ファイルリンク・コードフェンス・差分の空白も確認している。

変更は計画文書のみであり、Rust の spin-one 実装・crate 追加・ビルド・テストは行っていない。この補助計算は C の各能力の実装や f32 の精度を検証するものではなく、8節の受入条件は実装後に検証する。
