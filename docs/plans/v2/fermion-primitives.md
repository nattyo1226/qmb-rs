# Fermion の局所 primitive 実装案

作成日: 2026-10-29
更新日: 2026-10-29
状態: 検討案。記載する crate・型・追加 API は未実装。

レビュー中。状態と基本4演算子、Majorana 系の3演算子の su(2) 規格化、独立した Parity 演算子の保持を合意済みとする。演算子基底は LadderBasis と MajoranaBasis の二種類を提供する。標準 Hubbard と symmetric Hubbard のモデル構築も実装範囲に含める。合意と残る提案を2.1節に記録する。

## 1. 対象と既存実装との関係

単一 fermion mode の二次元占有空間、その上の全演算子代数、および局所基底を実装する案を記録する。局所の積・作用・内積・随伴・通常 trace・基底展開は、現在の `qmb-core` の能力だけで実装・検証できる。因子間の符号、空間の種類、Jordan–Wigner string は [実装済みの共通評価機構](../../design/composite-expression.md#12-ordinary--graded-の共通評価と移行)の担当とし、本書の局所規則に埋め込まない。

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

Workspace に fermion crate はまだ存在しない。Core の `Parity`・`HomogeneousParity`・`ParityOperator` は実装済みであり、[局所契約](../../design/local-structure.md#8-grading-の局所能力と移行)に従って接続する。spin-half は通常のスピン行列を提供するが、別因子の反交換関係を持つ fermion と同一の型にはしない。

## 2. 提案する crate と公開型

`qmb-fermion` を追加し、`no_std + alloc`、`unsafe` 禁止、crate root からの主要型・基底の再公開という既存の具体的代数の構成に合わせる。

```rust
enum FermionState { Empty, Occupied }

enum FermionOperator {
    I,
    Annihilate,
    Create,
    Number,
    MajoranaX,
    MajoranaY,
    MajoranaZ,
    Parity,
}
```

全 variant は有効な局所値なので、整数ラベルや占有数を受ける検証付き constructor は不要である。`Clone, Copy, Debug, PartialEq, Eq, Hash` を実装する。Zero は primitive に追加せず、`LinearCombination::zero()` と共通式の zero を使う。

`Fermion<C>` を `PhantomData<fn() -> C>` の型選択 marker とし、`QuantumSystem` の State / Operator / Coeff を上記型と C に関連付ける。marker や primitive の保持・`state` / `operator` に評価・数値能力の bound を課さない。独自の式木・evaluator・演算子構文は作らない。

想定ファイルは `src/{lib,state,operator,basis,models}.rs`、必要なら小さな定数生成をまとめる `coefficient.rs`、`tests/{contracts,minimal_scalar,local_expressions,models}.rs` とする。本番依存は既存 workspace の `qmb-core`・`num-traits` とし、モデルの Chain 接続・エラー処理に `qmb-geometry`・`thiserror` を使う。`num-complex` は両精度の検証用 dev-dependency とする。局所規則・基底は geometry に依存せず、Index 発行・数値 backend を本番依存に加えない。

### 2.1. レビューでの合意と追加提案

2026-10-29 のレビューで、以下を確認した。

- 状態 primitive `FermionState { Empty, Occupied }` は確定する。
- 演算子 primitive の `I`・`Annihilate`・`Create`・`Number` は承認済みとする。
- 単一 mode の二つの Majorana 演算子と、その積に係数を掛けて Hermite にした第三の演算子を追加する。規格化は既存 spin-half の `Sα=σα/2`、`hbar=1`、`[Sx,Sy]=iSz` と合わせる。
- `Parity` は `MajoranaZ` と線型従属でも独立した variant として保持する。
- 演算子基底は既定の `LadderBasis` と、Majorana 系による `MajoranaBasis` の二種類を提供する。従来案の `ParityBasis` は今回の提供対象から外す。
- 係数型には各操作で必要な trait bound を個別に要求する。全 variant を扱う積・作用などに分数と虚数単位の能力が必要になることを受け入れ、式の保持・構築には数値制約を波及させない。
- モデル構築として標準 Hubbard と symmetric Hubbard を提供する。後者は前者の相互作用項 `U n↑n↓` だけを `U(n↑−1/2)(n↓−1/2)` に置き換える。

追加 variant の仮称を `MajoranaX`・`MajoranaY`・`MajoranaZ` とし、3節の X・Y・Z に対応させる。名称は引き続き検討対象とする。

生成・消滅演算子を使う演算子基底は、既存 spin-half の命名に合わせて `LadderBasis` と呼ぶ。状態基底の `OccupationBasis` と区別し、両者に Operator / State の長い接尾辞は付けない。

`Parity` の独立 variant は、core の `ParityOperator::parity_operator()` が係数を使わず P 自身を返す契約にも対応する。`MajoranaZ=P/2` をそのまま代用できない。以下の規則表・検証件数は、この `Parity` を含む8 variant を前提とする。

## 3. 規約・次元・parity

状態の順序を `(|0>, |1>) = (Empty, Occupied)` とし、内積を `⟨m|n⟩ = δ_mn` とする。`a = Annihilate`、`a† = Create`、`n = Number`、`P = Parity` と書く。追加演算子は `X = MajoranaX`、`Y = MajoranaY`、`Z = MajoranaZ` と書く。

```math
a=\begin{pmatrix}0&1\\0&0\end{pmatrix},\quad
a^\dagger=\begin{pmatrix}0&0\\1&0\end{pmatrix},\quad
n=\begin{pmatrix}0&0\\0&1\end{pmatrix},\quad
P=\begin{pmatrix}1&0\\0&-1\end{pmatrix}=I-2n.
```

標準的な Majorana の規約を `γ₁=a+a†`、`γ₂=−i(a−a†)` と置く。この定義と `{γj,γk}=2δjk I` は [Kitaev, Unpaired Majorana fermions in quantum wires, 式 (1)–(2)](https://arxiv.org/pdf/cond-mat/0010440) に合わせる。以下はこの規約から導出する。

```math
X=\frac{\gamma_1}{2}=\frac{a+a^\dagger}{2}
 =\frac12\begin{pmatrix}0&1\\1&0\end{pmatrix},\qquad
Y=\frac{\gamma_2}{2}=-\frac{i}{2}(a-a^\dagger)
 =\frac12\begin{pmatrix}0&-i\\i&0\end{pmatrix},
```

```math
Z=-\frac{i}{2}\gamma_1\gamma_2=-2iXY=\frac{I}{2}-n
 =\frac12\begin{pmatrix}1&0\\0&-1\end{pmatrix}=\frac{P}{2}.
```

X・Y は標準規格化の Majorana 自身の1/2倍である。Z はそれらの積を Hermite にした even な双線型演算子であり、第三の独立した odd Majorana ではない。三つすべてが Hermite で、既存 spin-half と同じ通常の交換関係を満たす。

```math
[X,Y]=iZ,\qquad [Y,Z]=iX,\qquad [Z,X]=iY,\qquad
X^2=Y^2=Z^2=\frac14 I,\qquad
X^2+Y^2+Z^2=\frac34 I.
```

`a=X+iY`、`a†=X−iY`、`n=I/2−Z`、`P=2Z` である。占有基底上の行列は spin-half と同じだが、X・Y は odd、Z は even であり、全演算子が even な spin-half とは異なる型・grading を持つ。局所の通常交換子による su(2) と、異なる mode 間の graded 合成則を区別する。

これらの行列で規約を固定し、局所の match 規則とは独立したテストの参照値にする。製品の局所規則を密行列の乗算として実装する指定ではない。

| 対象 | 値 |
| --- | --- |
| `FermionState::dim()` | 2 |
| `FermionOperator::dim()` | 4。8個の primitive は線型従属 |
| `Identity::identity()` | `FermionOperator::I`。係数型に依存しない |
| 状態の parity | Empty は even、Occupied は odd |
| 演算子の parity | I・Number・Parity・MajoranaZ は even、Annihilate・Create・MajoranaX・MajoranaY は odd |

P を primitive に持つ理由は、係数型に依存せず P 自身を生成して core の grading 契約に接続するためである。これは全 variant を扱う Act の係数 bound を軽くするものではない。`P = I - 2n` を保持しても全演算子空間は四次元のままで、Dimension や基底要素数を8にしない。

## 4. 局所能力の実装

### 4.1. Product の全規則

表の行を左 operand、列を右 operand とし、右端が先に作用する。結果は同じ O の `LinearCombination<O,C>` とする。

| 左 × 右 | I | a | a† | n | P | X | Y | Z |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| I | I | a | a† | n | P | X | Y | Z |
| a | a | 0 | I − n | a | −a | (I − n)/2 | i(I − n)/2 | −a/2 |
| a† | a† | n | 0 | 0 | a† | n/2 | −in/2 | a†/2 |
| n | n | 0 | a† | n | −n | a†/2 | ia†/2 | −n/2 |
| P | P | a | −a† | −n | I | iY | −iX | I/2 |
| X | X | n/2 | (I − n)/2 | a/2 | −iY | I/4 | iZ/2 | −iY/2 |
| Y | Y | in/2 | −i(I − n)/2 | −ia/2 | iX | −iZ/2 | I/4 | iX/2 |
| Z | Z | a/2 | −a†/2 | −n/2 | I/2 | iY/2 | −iX/2 | I/4 |

`Product<C, Term=FermionOperator>` は64組の直接 match とし、identity の積では元の自然な primitive を返す。例えば `a a†` は `I − n` の二項、nilpotent な積は空の線型結合とする。通常の Product に graded 交換則を追加しない。ここでの局所反交換関係そのものが具体的代数の規則である。

表から `a²=(a†)²=0`、`aa†+a†a=I`、`n=a†a`、`n²=n`、`P²=I`、`Pa=-aP`、`Pa†=-a†P` と、3節の su(2) の関係が得られる。P と a の掛け順、および `Z=−2iXY` の係数・符号を取り違えない。

### 4.2. Act・内積

`Act<C>` は `Operator=FermionOperator`、`Term=FermionState` とする。

| 演算子 | Empty への作用 | Occupied への作用 |
| --- | --- | --- |
| I | Empty | Occupied |
| a | 0 | Empty |
| a† | Occupied | 0 |
| n | 0 | Occupied |
| P | Empty | −Occupied |
| X | Occupied/2 | Empty/2 |
| Y | i Occupied/2 | −i Empty/2 |
| Z | Empty/2 | −Occupied/2 |

状態ごとの parity が明確なので、P の作用と even / odd の定義が一致する。混合 parity の状態は共通の局所状態式の和で表す。

### 4.3. Adjoint・Trace

`Adjoint<C, Term=FermionOperator>` は a と a† を交換し、I・n・P・X・Y・Z を保存する。primitive の係数は1とし、式の係数の共役は core が担当する。

| 演算子 | I | a | a† | n | P | X | Y | Z |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 通常の非正規化 trace | 2 | 0 | 0 | 1 | 0 | 0 | 0 | 0 |

`Tr(I)=2` とし、正規化 trace や supertrace に変更しない。局所 trace の実装に状態基底・行列を要求しない。

## 5. 基底と展開

状態基底 `OccupationBasis` は `[Empty, Occupied]` とする。演算子基底は次の二種類とし、順序も固定する。型を分けて `Basis` の対象を明確にする。

| 基底 | 順序付き要素 |
| --- | --- |
| `LadderBasis`（既定） | `[I, Annihilate, Create, Number]` |
| `MajoranaBasis` | `[I, MajoranaX, MajoranaY, MajoranaZ]` |

いずれも四要素が独立で全演算子空間を張り、各要素は homogeneous な parity を持つ。LadderBasis の I と n の Hilbert–Schmidt 内積は0ではなく、正規直交基底とは扱わない。

MajoranaBasis は Hermite 演算子だけからなり、Hilbert–Schmidt 内積で直交する。ただし `Tr(I†I)=2`、`Tr(X†X)=Tr(Y†Y)=Tr(Z†Z)=1/2` なので正規直交基底ではない。su(2) 規格化を保ち、基底のための追加の規格化はしない。Hermite な局所演算子はこの基底で実係数の線型結合になる。

全 primitive の展開は次のとおりとする。

| primitive | LadderBasis への展開 | MajoranaBasis への展開 |
| --- | --- | --- |
| I | I | I |
| a | a | X + iY |
| a† | a† | X − iY |
| n | n | I/2 − Z |
| P | I − 2n | 2Z |
| X | (a + a†)/2 | X |
| Y | −i(a − a†)/2 | Y |
| Z | I/2 − n | Z |

`index` は各基底の要素だけに上表の基底順の番号を返し、それ以外には `None` を返す。特に P は両基底で `None` となる。`elements()` は基底の四要素だけを一度ずつ返し、線型従属な全 primitive を列挙しない。

両基底間の直接展開・再展開と、局所式の作用・trace の一致を検証する。基底の順序を演算子の行列成分と同一視しない。ParityBasis は、具体的な用途が必要になった段階で追加を検討する。Parity 演算子自体は引き続き提供し、どちらの基底でも展開できる。

## 6. 係数型と能力の最小化

各操作の実装に必要な trait bound を要求する方針をレビューで合意済みとする。既存 `LinearCombination::primitive` が `C: One` を要求することと、`One`・`Zero` 自体が持つ算術上の制約を前提とする。次は具体的 impl の想定 bound であり、式の保持へは伝播させない。全能力をまとめて一律に要求する係数 trait は導入せず、実際の bound は各操作で使う能力に合わせて定める。

| 能力 | 想定する C の能力 |
| --- | --- |
| Identity、Dimension、parity 問い合わせ、QuantumSystem | 数値制約なし |
| InnerProduct、Trace | `Zero + One`。2は `one + one` で生成 |
| Adjoint、状態基底 | `One` |
| Product、Act、LadderBasis、MajoranaBasis | `One + Add<Output=C> + Div<Output=C> + Neg<Output=C> + ImaginaryUnit`。1/2・1/4・i を表せることも数学的契約。零は空の線型結合で返す |

局所規則には Float・平方根・C の Clone / Copy を要求しない。一方、Majorana 系の追加で虚数係数と分数が必要になるため、符号付き整数だけでは全局所規則を表せない。全 variant を対象にする trait impl の bound は、呼び出した値が a・n などでも同じである。独自の厳密な複素有理係数などで検証し、整数の切り捨て除算は使わない。共通式の評価・集約・export には、それぞれ core が要求する Clone・共役・算術能力が別途必要である。

精度は C が決め、`Complex<f32>` と `Complex<f64>` を暗黙に混ぜない。定数のために f64 経由のキャストやパースを追加しない。

## 7. Tensor 実装から独立した接続方針

局所 primitive の S・O は core の `HomogeneousParity` を実装する。O は `ParityOperator` も実装し、係数型に依存せず `FermionOperator::Parity` 自身を返す。局所の parity と P の契約は[局所能力](../../design/local-structure.md#8-grading-の局所能力と移行)を参照し、重複定義しない。

多 mode の順序・交換符号・先行区間の P は core の合成機構に委ねる。`FermionOperator` に Index・Site・mode 番号・parity string・global ordering を保持しない。既存 ordinary tensor へ複数の a を配置するだけでは fermion の異なる mode の CAR は得られないため、その構成を fermion 多体系の利用例として案内しない。

spinful の内部順序は `(up,down)`、二重占有の位相は `|up down>=f_up† f_down†|0>` とする。合成側との接続試験では、`f_up|up down>=|down>`、`f_down|up down>=-|up>`、総 parity `diag(1,-1,-1,1)` を確認する。表示順 `[empty, up, down, double]` は、二 mode の mixed-radix 順 `[empty, down, up, double]` と異なるので明示的に対応付ける。この表示順の違いを単一 mode primitive の符号へ埋め込まない。

### 7.1. Hubbard モデルの構築

標準 Hubbard と symmetric Hubbard の二種類を今回の実装範囲に含める。以下の Hamiltonian は合意した相互作用の違いを具体化したものであり、API と配置入力は既存 spin-half に合わせた提案とする。

#### Hamiltonian の規約

まず一様な実パラメータ t・U・μ、Chain の最近接辺、spin を保存するホッピングを対象とする。

```math
K(t)=-t\sum_{\langle i,j\rangle,\sigma}
  \left(f_{i\sigma}^\dagger f_{j\sigma}
       +f_{j\sigma}^\dagger f_{i\sigma}\right),\qquad
N=\sum_{i,\sigma}n_{i\sigma},
```

```math
H_{\mathrm{Hubbard}}(t,U,\mu)
 =K(t)+U\sum_i n_{i\uparrow}n_{i\downarrow}-\mu N,
```

```math
H_{\mathrm{symmetric}}(t,U,\mu)
 =K(t)+U\sum_i\left(n_{i\uparrow}-\frac I2\right)
                  \left(n_{i\downarrow}-\frac I2\right)-\mu N.
```

標準形は [ALPS の Hubbard 定義](https://alps.comp-phys.org/documentation/models/hubbard/)、symmetric 形は [TRIQS/TPRF の定義](https://triqs.github.io/tprf/3.1.x/user_guide/Linearized%20Eliashberg%20equation%20on%20the%20attractive%20Hubbard%20model.html) と同じ符号規約である。化学ポテンシャル μ を共通の引数に含める案とし、化学ポテンシャル項が不要なら0を渡す。

symmetric でもホッピング・化学ポテンシャル項は変更せず、定数項を落とさない。L Site に対して、展開から次の恒等式が得られる。

```math
H_{\mathrm{symmetric}}(t,U,\mu)
 =H_{\mathrm{Hubbard}}(t,U,\mu)
   -\frac U2 N+\frac{UL}{4}I
 =H_{\mathrm{Hubbard}}\left(t,U,\mu+\frac U2\right)+\frac{UL}{4}I.
```

μ の読み替えを API が自動で行うことはしない。`n−I/2=−Z` なので、symmetric な局所相互作用は `U Z↑Z↓` とも表せる。両表現の値の一致を検証する。

#### 公開入口と配置の案

既存の `SpinHalf<C>` と同様に `models.rs` で構築し、`Fermion<C>` の関連関数として公開する。QuantumSystem trait にモデル固有メソッドを追加しない。

```rust,ignore
Fermion::hubbard(&chain, &sites, t, u, mu)
Fermion::symmetric_hubbard(&chain, &sites, t, u, mu)
// いずれも Result<OperatorOf<Fermion<C>, I>, UtilityError> を返す想定。
```

- `sites: &[HilbertSpace<I>]` の一要素を一つの物理 Site とし、各 Site は葉順 `(up,down)` の二 mode の graded 空間とする。Index は利用側で発行し、関数内で再発行しない。
- 全作用空間は Site 順、各 Site 内の up・down 順で graded 合成し、入力の Site 空間を共有して構築する。返す演算子はこの全空間に bind する。
- ホッピングは隣接 Site の同じ spin 同士だけを結ぶ。spin-half の XYZ が使う全葉ペアの直積は用いない。各辺に両方向のホッピングを一度ずつ加え、符号は core の graded 合成に委ねる。
- open / periodic の辺は既存 Chain の契約に従う。0・1 Site は辺なし、2 Site の periodic でも辺は一つとし、ホッピングを重複させない。
- Site 数の一致、各 Site が二つの葉を持つ graded 空間であること、全 mode の Index 重複を、結合値が0でも検証する。不正な Site 構成もエラーとして返す。エラー型の詳細は実装時に定める。
- 0 Site は graded な単位空間に束縛した zero Hamiltonian を返す。1 Site はホッピングなしで相互作用・化学ポテンシャル項を保持する。
- モデル構築に必要な係数の変換・複製・符号反転・定数生成だけを各関数の bound とし、局所 Product・Act や Float を構築のために要求しない。t・U・μ が実数に対応することはモデルの数学的契約とする。

## 8. 実装順序と受入条件

1. 新 crate、enum、Dimension・Identity、型選択 marker と再公開を追加する。
2. Product・Act・InnerProduct・Adjoint・Trace と LadderBasis・MajoranaBasis を実装する。
3. core の局所式を通した評価・展開・再展開・座標取得、既存 export / normalize の単一 mode 利用を検証する。
4. 実装済みの core の grading trait と graded 空間の構築 API に接続し、多 mode / spinful の統合検証を別途行う。
5. 標準 Hubbard・symmetric Hubbard のモデル構築、利用例、独立な小規模 Fock 空間の参照行列との統合検証を追加する。

局所実装の完成条件は以下とする。これらは実装後の検証計画であり、現時点の Rust テスト結果ではない。

- 64組の Product、16組の Act、8種の Adjoint / Trace、4組の内積を、3節の独立な2×2行列・二成分ベクトルと照合する。
- primitive の全三重積512組の結合則、全 Product / 状態の128組の作用整合性、単位元、随伴の反準同型性・二回の随伴、trace の巡回性を検証する。構造の一致ではなく値を比較する。
- 局所 CAR、n の射影性、P の自己随伴性・平方・作用、各非零返却項の parity を確認する。X・Y・Z の Hermite 性、su(2) の巡回交換関係、各二乗・Casimir、`Z=−2iXY=P/2`、a・a†・n・P との変換も検証する。
- 状態基底2要素・両演算子基底それぞれ4要素の列挙 / 番号 / 完全性と、全 primitive の展開・両基底間の再展開を確認する。`P=I−2n=2Z` を経由した評価も同じ値になる。MajoranaBasis の Hilbert–Schmidt Gram 行列が `diag(2,1/2,1/2,1/2)` であることも確認する。
- `Complex<f32>` / `Complex<f64>` の局所式の作用・期待値・通常 trace・正規化・座標・2×2行列 export を照合する。Vec と任意 feature の ndarray / sprs を同じ基底順で比較する。
- Float のない厳密な複素有理係数で全局所規則を利用できること、Clone / Copy のない係数で primitive の規則を利用できること、評価能力のない型でも式を構築できることを compile / capability テストで確認する。
- spin-half と独立した型として動作し、core に fermion 固有分岐を追加しない。多 mode の CAR・spinful の接続は局所完成条件と別に記録する。

モデル実装の完成条件は局所代数と分けて記録し、今回の実装範囲として以下を検証する。

- 二 mode の CAR、spinful の二重占有からの消滅符号、および全空間の mode 順を確認する。
- 1 Site の表示順 `[empty,up,down,double]` で、標準形の対角成分が `[0,−μ,−μ,U−2μ]`、symmetric 形が `[U/4,−U/4−μ,−U/4−μ,U/4−2μ]` となることを確認する。実際の列挙順との対応は明示する。
- 2 Site（4 mode、16状態）の両モデルを、core の評価・export を使わない Fock 空間の符号規則から構築した参照行列と比較する。実パラメータでの Hermite 性、粒子数・各 spin の粒子数の保存、`t=0`・`U=0` の極限を確認する。
- 両モデルの差が `−(U/2)N+(UL/4)I` であることと、symmetric 相互作用の `U Z↑Z↓` による構築を、定数項を含めて比較する。
- 0・1・2 Site と3 Site 以上の open / periodic、Site 数不一致、mode 数の過不足、ordinary な Site、Index 重複、非 Clone Index を検証する。
- `Complex<f32>` / `Complex<f64>` と Vec・ndarray・sprs の同じ基底順で値を比較し、構築の係数 bound に評価能力が混入していないことも確認する。

実装時の検証入口は `cargo test -p qmb-fermion`、workspace 登録後の `cargo test --workspace --all-features`、`cargo check -p qmb-fermion --no-default-features`、`cargo fmt --all -- --check` とする。ホスト上の no-default-features check と専用 no_std target の検証は区別し、未実施事項を明記する。

## 9. 設計判断と残る詳細

レビューで合意した構成は、a・a†・n・I と、spin-half に規格化を合わせた X・Y・Z、および grading 用の P を直接保持する八つの線型従属な primitive とする。四つの行列単位だけにすると規則は単純になるが、現行 `Identity` に必要な I は別途 primitive の値として追加する必要がある。物理的な名前と P の直接保持を優先する。演算子基底には LadderBasis と MajoranaBasis を提供し、ParityBasis は今回の対象に含めない。

追加 variant・基底の最終命名と、7.1節のモデル API・配置入力の案は実装前に固定する。Core の grading trait と tensor API は実装済みであり、その契約に接続する。局所代数の規約・表・次元は多 mode の統合検証と独立に確定・検証できる。標準 Hubbard と symmetric Hubbard は今回の対象に含める。汎用の spinful 構築 utility、粒子数 sector、任意格子・Site ごとや辺ごとの可変結合、複素ホッピング・磁場・pairing は後続に分ける。

## 10. 本ドキュメント作成時の確認

初稿作成時に設計文書・core 契約・spin-half 実装を照合した。補助的な Python 標準ライブラリの行列計算では、Majorana 追加前の表を読み取って25組の積と10組の作用を参照行列と照合し、局所 CAR、P の平方・交換関係、trace の値も確認した。初稿の相対ファイルリンク・コードフェンス・差分の空白を確認した記録であり、追加後の確認は11節に記録する。

変更は計画文書のみであり、Rust の fermion 実装・crate 追加・ビルド・テストは行っていない。8節の受入条件は実装後に検証する。

## 11. Majorana 追加レビュー時の確認

Python 標準ライブラリの補助計算で、文書の表を読み取った64組の積・16組の作用を、独立に記述した2×2参照行列・二成分ベクトルと照合し、全件一致した。8演算子の随伴・trace、局所 CAR、parity、su(2) の巡回交換関係・Casimir、および追加時点の LadderBasis・ParityBasis への展開の恒等式も確認した。計算に使う整数・二進分数・虚数単位と中間値は二進浮動小数点で厳密に表せる範囲であり、許容誤差による比較は使っていない。差分の空白も確認した。

Parity の保持と MajoranaBasis の採用を反映した後、更新した基底展開表の全16展開を独立な参照行列と照合し、両方向の基底変換行列の積が単位行列になること、MajoranaBasis の Gram 行列が `diag(2,1/2,1/2,1/2)` になることを確認した。既存の64積・16作用の照合も再実行し、全件一致した。

変更は計画文書のみであり、追加後の Rust 実装・ビルド・テストは未実施である。

## 12. Hubbard モデル追加時の確認

Python 標準ライブラリの Fraction を使う独立な Fock 空間の参照計算で、0〜3 Site、open / periodic、非零結合・`t=0`・`U=0` の計24ケースを確認した。両モデルの差、μ のシフトと定数項を含む恒等式、1 Site の上記対角成分、Hermite 性、各 spin の粒子数保存、および `(n↑−1/2)(n↓−1/2)=Z↑Z↓` を厳密比較で確認した。

これは計画書の数式・検証条件を確認する補助計算であり、Rust のモデル構築・graded 評価・export の統合検証は未実施である。
