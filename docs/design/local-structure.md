# 局所 primitive の設計

更新日: 2026-11-01

## 1. 対象と検討状況

本書は局所状態・演算子 primitive の能力と整合性を定める。LocalStructure と追加の定義ラベルを廃止し、全 crate の式・空間・幾何・利用例を S・O の直接指定へ移行した。採用理由と契約は第7節に記録する。移行時には primitive 型の混在拒否、同じ空間の異なる primitive 型の式による共有、深い式の評価・破棄、および両精度と独自代数を検証した。

局所式は `LocalState<S, O, C>`・`LocalOperator<O, C>`、複合式は `State<I, S, O, C>`・`Operator<I, O, C>` とする。演算子式に状態型を要求しない。状態式は遅延した作用節点を保持するため O も指定する。保持・共有・係数計算を伴わない遅延構築に演算能力や primitive の Clone を要求しない。

同じ O が異なる S に作用する場合は、各 S の Act 実装で対応する。状態の構築時に O を推論できなければ、型注釈・型別名や後続の利用箇所から指定する。借用を含む primitive を保持でき、保持のために一律の Send・Sync・'static 制約を課さない。スレッド間の移動・共有の可否は保持する型と Arc の制約に従う。

## 2. 空間の情報と能力の関連付け

| 責任 | 実装した配置 |
| --- | --- |
| Primitive の意味・局所規則 | S・O と個別の能力実装。追加の定義 ID・ラベルは設けない |
| 局所状態空間の次元 | S の Dimension。primitive の値全体が張る空間の次元 d |
| 局所演算子空間の次元 | O の Dimension。全演算子空間では d²。variant の個数とは異なる |
| 恒等演算子 | `Identity::identity() -> Self` を O が実装する。係数型に依存しない |
| 作用の型対応 | 評価時の `S: Act<C, Operator = O, Term = S>` と `O: Adjoint<C, Term = O>` 等 |
| Scalar の整合性 | C と各操作の bound。異なる精度を暗黙に混ぜない |
| 基底 | `Basis<S, C>`・`Basis<O, C>`。列挙数は対象 primitive の span の次元と一致する |
| 因子構造・幾何 | `HilbertSpace<I>` が因子順序・grouping、Chain がサイト数と境界条件だけを保持する。配置は外部の局所空間列で指定する |

Dimension は `fn dim() -> Result<NonZeroUsize, DimensionOverflow>` として実装した。次元を整数化できない場合は取得時にエラーを返す。`dim_state<S>()`・`dim_op<O>()` は対象の次元をそのまま返す。Operator の d² を再度二乗しない。能力は次元を使う操作だけが要求する。

恒等元は primitive の値として表せることを要求するが、選択基底の要素である必要はない。基底に含まれなければ Basis が展開する。線型従属な primitive を許容するため、恒等元用の値を追加しても基底の次元は増えない。局所式の恒等元は `Primitive(O::identity())` として保持し、作用・積・随伴・trace・基底展開を通常の primitive と同じ経路で処理する。`LocalOperator::identity()` と空の `product` は `O: Identity` を要求するが、係数の数値能力は要求しない。複合式の記号的な `OpNode::Identity` は、未束縛の空の積や全作用空間の恒等演算子として維持する。

局所 identity と spectator 因子の trace は、O の Identity・Trace から得る。次元の整数値を係数へ変換する経路に置き換えず、局所 trace に係数の One・整数変換を追加要求しない。外部の O に Identity を追加実装できない場合は、提供元の実装か利用者が所有する wrapper で接続する。`Identity`・`Dimension` は `qmb_core::ops` と crate root から公開する。`Identity` は局所式（`O: Identity`）と複合式（追加の能力制約なし）にも実装する。

## 3. 数学的能力の整合性

以下の `A,B,C` は対応する primitive 型 O の演算子、`ψ,φ` はその状態、`I` は恒等演算子とする。primitive の線型結合も含め、次の契約を要求する。

- 閉性・所属の保存: Product・Act・随伴の結果は対応する S・O の型に閉じる。具体的な表現と各能力の実装が整合することを要求する。
- 線型性: 積と作用はそれぞれの引数について線型、内積は第一引数について共役線型・第二引数について線型、随伴は共役線型、trace は線型とする。
- 積の結合則と作用の整合性: `(AB)C = A(BC)` および `(AB)ψ = A(Bψ)` を満たす。
- 単位元: `IA = AI = A` および `Iψ = ψ` を満たす。
- 内積: 共役対称性と、状態空間上での正定値性を満たす。非独立な primitive では、非零の係数列でも zero state を表し得る。
- 随伴: 内積・作用との以下の関係、および `(A†)† = A`、`(AB)† = B†A†` を満たす。
- Trace: 局所 Hilbert 空間上の通常の非正規化 trace と一致し、`Tr(I) = d` を満たす。
- 線型関係の保存: 同じ状態・演算子を表す異なる primitive の線型結合に演算を適用したとき、数学的に同じ結果を得る。展開係数や式木の構造の一意性は要求しない。
- Grading〈対応する場合〉: homogeneous な primitive の parity と局所演算規則が整合する。因子間の符号処理は共通の複合代数機構に委ねる。

```math
\langle\phi,A\psi\rangle=\langle A^\dagger\phi,\psi\rangle.
```

これらは局所実装が満たす数学的な契約と検証対象である。任意の実装に対して結合則などを実行時に自動証明する機構は要求しない。初期版では S・O の型と能力実装が意味・規則を定め、追加の定義ラベルは設けない。実行時パラメータを持つ構造への拡張は後続へ延期する。

等式は数学的な意味について要求し、浮動小数点計算のビット単位の一致とは区別する。数値精度、特殊値、および近似比較は要件の 4.1 節・5.6 節に従う。

## 4. 局所構造の一致条件

S・O の型と個別の能力実装が primitive の意味と数学的規約を定める。同じ S・O を、別の L や定義 ID を付けて区別しない。異なる規則や表現が必要なら primitive の型と能力実装に反映する。既存の同じ primitive 型を別のラベルだけで隔離する契約は廃止した。

状態式の tensor・和・内積は同じ S・O・C、演算子式の和・積は同じ O・C の入力を受け取る。作用は対応する O を持つ状態式と演算子式を組み合わせ、評価時に Act の対応を要求する。型の一致は数学的な法則の証明ではなく、具体的代数の実装・検証がその責任を担う。

Index・Site・表示名・選択基底とその列挙順序・メモリ上の object identity・共有方法を理由に新しい定義ラベルを作らない。数学的同型の判定、異なる型間の暗黙の変換、primitive の自動対応付けは導入しない。

## 5. Index の契約との接続

[Index](index.md)と局所代数の意味は別の責任である。空間は因子の Index・順序・grouping を保持し、primitive 型は式が保持する。因子の重複や配置先・束縛先の存在は空間の構造から検証し、異種の式の不適切な合成は S・O・C の型制約で拒否する。Index と局所定義の全体登録や、一律の TypeId・'static 制約は導入しない。

空間の次元は `space.dim::<P>()` で指定した primitive の span を各因子に適用するか、`space.dim_with(dim_local)` で明示する。S を指定すると状態空間、O を指定すると演算子空間の次元を取得する。この取得操作は空間に primitive の定義を登録しない。空の tensor は次元1とし、通常の作用・内積・trace に総次元の整数化を要求しない。

## 6. 後続事項

一般の U→V の写像と State／Operator の型の統合は後続課題とする。Grading の情報保持・API は実装済みであり、局所契約を第8節、共通評価を[複合式の設計](composite-expression.md#12-ordinary--graded-の共通評価と移行)に記録する。ordinary な空間の責任分離は実装済みである。Chain はサイト数と境界条件だけを保持し、外部の局所空間列を借用する構成へ2026-10-28に移行した。契約は[Chain の設計](geometry.md)に従う。異種代数の型消去や実行時の表現パラメータを今回の移行には追加しない。

## 7. LocalStructure 廃止の採用方針（2026-10-26）

S・O の直接指定、primitive の span の次元と恒等元、追加の定義ラベルを設けない方針を採用し、実装へ反映した。API は第1〜2節、型と Index の責任分担は第4〜5節を現在の契約とする。旧 LocalStructure・Identity<C> と SpinHalf の定義用 marker は公開しない。

局所規則は S・O・C と個別の能力実装で決まり、定義用の L は規則を切り替える引数ではなかった。次元と恒等元を primitive 側へ移したため、型の組だけを保持する LocalStructure や空間側の定義ラベルは不要とした。Act の関連型から保持型の O を導く構成は、保持・遅延構築に評価能力を要求するため採用しない。全能力をまとめた StatePrimitive・OperatorPrimitive trait、必須の FactorLayout 層、実行時の Definition／provider も導入しない。

## 8. Grading の局所能力と移行

`qmb_core::ops` と crate root から、係数型に依存しない次の能力を公開する。

```rust,ignore
enum Parity { Even, Odd }
trait HomogeneousParity {
    fn parity(&self) -> Parity;
}
trait ParityOperator: Sized {
    fn parity_operator() -> Self;
}
```

Graded な局所代数の S・O は homogeneous な primitive とし、両方が `HomogeneousParity` を提供する。混合 parity は `LocalState` / `LocalOperator` の和で表し、混合した式に単一の parity を割り当てない。未対応を暗黙に even とする blanket impl は提供しない。

O は `ParityOperator` で局所 parity 演算子 P 自身を primitive の値として返す。Identity と同様に選択基底の要素である必要はなく、例えば基底展開で `P = I - 2n` と表せる。線型従属な P の値を追加しても span の次元は変わらない。

契約は `P† = P`、`P² = I`、`Pψ = (-1)^p(ψ) ψ`、`POP = (-1)^p(O) O` とする。Identity と P は even、随伴は parity を保存し、局所 Product / Act の非零な出力項の parity は入力の和 modulo 2 と一致する。状態の even・odd 部分は内積で直交する。局所基底の primitive も同じ契約に従う。これらの法則は具体的代数の実装・テストで保証し、core は任意の実装について実行時に自動証明しない。

Ordinary 専用の独自代数も、複合評価の共通 API を利用する場合は必要な grading trait を明示実装する。自明な grading を採る代数では全状態・演算子を even、P を I とする。

```rust,ignore
impl HomogeneousParity for MyState {
    fn parity(&self) -> Parity { Parity::Even }
}
impl HomogeneousParity for MyOperator {
    fn parity(&self) -> Parity { Parity::Even }
}
impl ParityOperator for MyOperator {
    fn parity_operator() -> Self { Self::identity() }
}
```

Spin-half はこの自明な grading を明示実装する。既存モデルの空間は ordinary を維持し、graded 空間でも同じ数値結果となる。局所だけで完結する評価や、構築・配置・束縛には grading 能力を一律に要求しない。具体的な操作の bound は[複合式の設計](composite-expression.md#12-ordinary--graded-の共通評価と移行)を参照。

## 9. Fermion の局所代数

`qmb-fermion` は `no_std + alloc`、unsafe 禁止とし、係数型から独立した `FermionState`・`FermionOperator` を提供する。状態は `Empty/Occupied`、演算子は `I/Annihilate/Create/Number/Mx/My/Mz` の7 variant とする。旧 `MajoranaX/Y/Z` と独立した `Parity` variant は公開しない。状態空間の次元は2、演算子空間の次元は4である。

`Fermion<C>` は `ExpressionFactory<FermionState,FermionOperator,C>` の型エイリアスである。公開入口は [lib.rs](../../crates/qmb-fermion/src/lib.rs)、局所規則は [state.rs](../../crates/qmb-fermion/src/state.rs)・[operator.rs](../../crates/qmb-fermion/src/operator.rs)、基底は [basis.rs](../../crates/qmb-fermion/src/basis.rs) に実装する。モデルは [Hubbard モデルの設計](geometry.md#6-fermion-の-hubbard-モデル)に従う。

### 9.1. 規約・次元・parity

状態順を `(|0>,|1>)=(Empty,Occupied)`、内積を Kronecker delta とする。`a=Annihilate`、`a†=Create`、`n=Number` と書く。

```math
a=\begin{pmatrix}0&1\\0&0\end{pmatrix},\quad
a^\dagger=\begin{pmatrix}0&0\\1&0\end{pmatrix},\quad
n=\begin{pmatrix}0&0\\0&1\end{pmatrix},
```

```math
M_x=a+a^\dagger=\begin{pmatrix}0&1\\1&0\end{pmatrix},\quad
M_y=\frac{a-a^\dagger}{i}=\begin{pmatrix}0&-i\\i&0\end{pmatrix},\quad
M_z=-iM_xM_y=I-2n=\begin{pmatrix}1&0\\0&-1\end{pmatrix}.
```

Mx・My は [Kitaev の式 (1)–(2)](https://arxiv.org/pdf/cond-mat/0010440) の標準規格化の Majorana に対応する。Mz は even な双線型演算子で parity 演算子 P 自身であり、第三の odd Majorana ではない。

```math
[M_x,M_y]=2iM_z,\quad [M_y,M_z]=2iM_x,\quad [M_z,M_x]=2iM_y,\quad
M_\alpha^2=I,\quad \sum_\alpha M_\alpha^2=3I.
```

Empty は even、Occupied は odd。a・a†・Mx・My は odd、I・n・Mz は even とする。`ParityOperator::parity_operator()` は Mz を返し、P の契約は第8節に従う。同一 mode 内の通常の局所積と、異なる mode 間の graded 合成則を区別する。

### 9.2. 局所能力

#### 9.2.1. Product の全規則

行を左 operand、列を右 operand とし、右端が先に作用する。

| 左 × 右 | I | a | a† | n | Mx | My | Mz |
| --- | --- | --- | --- | --- | --- | --- | --- |
| I | I | a | a† | n | Mx | My | Mz |
| a | a | 0 | I−n | a | I−n | i(I−n) | −a |
| a† | a† | n | 0 | 0 | n | −in | a† |
| n | n | 0 | a† | n | a† | ia† | −n |
| Mx | Mx | n | I−n | a | I | iMz | −iMy |
| My | My | in | −i(I−n) | −ia | −iMz | I | iMx |
| Mz | Mz | a | −a† | −n | iMy | −iMx | I |

49組を直接の match で返し、数学的に零の積は空の線型結合、恒等元との積は相手一項とする。局所の `a²=(a†)²=0`、`{a,a†}=I`、`n=a†a` と parity の法則を満たす。行列・基底変換を実行時の局所積に使わない。

#### 9.2.2. Act・内積

| 演算子 | Empty への作用 | Occupied への作用 |
| --- | --- | --- |
| I | Empty | Occupied |
| a | 0 | Empty |
| a† | Occupied | 0 |
| n | 0 | Occupied |
| Mx | Occupied | Empty |
| My | i Occupied | −i Empty |
| Mz | Empty | −Occupied |

状態式の線型性・半双線型性は core が提供する。

#### 9.2.3. Adjoint・Trace

随伴は a・a† を交換し、残る5要素を保存する。通常の非正規化 trace は I が2、n が1、他は0とする。HS 内積だけが `Tr(A†B)/2` による次元規格化を行う。

### 9.3. 基底と展開

| 基底 | 列挙順 |
| --- | --- |
| `OccupationBasis` | `[Empty,Occupied]` |
| `LadderBasis` | `[I,Annihilate,Create,Number]` |
| `MajoranaBasis` | `[I,Mx,My,Mz]` |

MajoranaBasis は `Tr(A†B)/2` に関して正規直交である。LadderBasis は I と n の内積が非零なので非直交であり、規格化を一律に適用しない。

```math
a=(M_x+iM_y)/2,\qquad a^\dagger=(M_x-iM_y)/2,\qquad n=(I-M_z)/2.
```

逆変換は第9.1節の定義を使う。index は選択基底のメンバーだけに列挙順の番号を返し、Mz は MajoranaBasis では3、LadderBasis では None とする。両基底間の再展開・座標・作用行列を照合する。

### 9.4. 係数型と能力

Identity・Dimension・parity と式の構築には数値制約を要求しない。InnerProduct・Trace は Zero + One、随伴・OccupationBasis は One。Product は One + Add + Neg + ImaginaryUnit、Act は One + Neg + ImaginaryUnit、LadderBasis は One + Add + Neg + ImaginaryUnit を要求する。MajoranaBasis はさらに Div を使う。算術の Output は C とする。

分数を使う操作では数学的な1/2を表せることを要求する。局所規則に Float・C の Clone / Copy を要求せず、共通式の評価・HS 内積・export に必要な能力は各操作で別途要求する。HS 内積は既存の Divide で `Tr(I)` による除算を行う。

### 9.5. Graded 合成・spinful への接続

mode の順序・交換符号・先行区間の parity string は core に委ね、primitive に Index・Site・mode 番号を保持しない。spinful の内部順序は `(up,down)`、二重占有の位相は `|up down>=f_up† f_down†|0>` とする。`f_up|up down>=|down>`、`f_down|up down>=−|up>`、総 parity `diag(1,−1,−1,1)` を検証する。

表示順 `[empty,up,down,double]` は二 mode の mixed-radix 順 `[empty,down,up,double]` と区別する。通常 tensor に配置するだけでは異なる mode の CAR は得られないため、fermion 多体系では graded 空間を使う。

### 9.6. 実装時の検証記録

名前・規格化の移行内容と受入検証は [v3 の修正記録](../plans/v3/primitive-normalization.md)に記録する。局所規則・基底の独立2×2行列との照合は [contracts.rs](../../crates/qmb-fermion/tests/contracts.rs)、非 Clone / Copy / Float の係数能力は [minimal_scalar.rs](../../crates/qmb-fermion/tests/minimal_scalar.rs)、式・全 backend は [local_expressions.rs](../../crates/qmb-fermion/tests/local_expressions.rs)、graded 接続・CAR・両 Hubbard の独立 Fock 空間参照は [models.rs](../../crates/qmb-fermion/tests/models.rs)で検証する。数値試験は f32 / f64 の両方を対象とする。

## 10. Spin-one の局所代数

`qmb-spin-one` は三状態の spin-1、15個の演算子 primitive、2種類の演算子基底と XYZ・Heisenberg・BLBQ を提供する。std を許可し、unsafe を禁止する。モデルの定義は [SpinOne のモデル構築](geometry.md#7-spinone-のモデル構築)に従う。

### 10.1. Primitive と規約

```rust
enum SpinOneState { Plus, Zero, Minus }
enum SpinOneOperator {
    I, Sx, Sy, Sz, Sp, Sm, Dxy, Q0, Qxy, Qyz, Qzx,
    Q1p, Q1m, Q2p, Q2m,
}
```

状態順は `(|+1>,|0>,|-1>)`、内積は Kronecker delta とする。Zero は磁気量子数0の正規化状態であり、zero vector は空の線型結合で表す。状態次元は3、演算子空間次元は9とする。

物理的な spin-1 演算子を J、`c=√(3/2)` とし、規格化した primitive は `S_alpha=c J_alpha` とする。J は ℏ=1 の通常の角運動量である。

```math
S_x=\frac{\sqrt3}{2}\begin{pmatrix}0&1&0\\1&0&1\\0&1&0\end{pmatrix},\quad
S_y=\frac{\sqrt3}{2}\begin{pmatrix}0&-i&0\\i&0&-i\\0&i&0\end{pmatrix},\quad
S_z=c\operatorname{diag}(1,0,-1).
```

四重極は次で定義する。

| primitive | 定義 |
| --- | --- |
| `Dxy` | `c (Jx²−Jy²)=(Sx²−Sy²)/c` |
| `Q0` | `(3Jz²−2I)/√2=√2(Sz²−I)` |
| `Qxy/Qyz/Qzx` | それぞれ `c {Jx,Jy}`、`c {Jy,Jz}`、`c {Jz,Jx}` |
| `Sp/Sm` | `Sx±iSy` |
| `Q1p/Q1m` | `Qzx±iQyz={Sz,Sp/m}/c` |
| `Q2p/Q2m` | `Dxy±iQxy=(Sp/m)²/c` |

Hermite な9要素の通常の HS ノルム二乗はすべて3で、`Tr(A†B)/3` に関して正規直交である。スピン3成分と四重極5成分の回転に関する部分空間は維持する。物理的なスピン J を得るには既存のスカラー倍で `√(2/3) S` を使い、補助 API は追加しない。

### 10.2. 基底と展開

| 基底 | 列挙順 |
| --- | --- |
| `SzBasis` | `[Plus,Zero,Minus]` |
| `CartesianBasis` | `[I,Sx,Sy,Sz,Dxy,Q0,Qxy,Qyz,Qzx]` |
| `LadderBasis` | `[I,Sp,Sm,Sz,Q0,Q1p,Q1m,Q2p,Q2m]` |

CartesianBasis は次元規格化した HS 内積に関して正規直交、LadderBasis は直交で I・Sz・Q0 のノルム二乗は1、昇降型6要素は2とする。昇降演算子はノルムの統一ではなく Cartesian 成分への直接の接続を優先する。通常の HS ノルム二乗はそれぞれ3・6となる。

各対の変換は `L_±=A±iB`、`A=(L_++L_-)/2`、`B=(L_+−L_-)/(2i)` とする。I・Sz・Q0 は両基底で共有し、非メンバーの index は None。9個の演算子座標と状態基底上の3×3作用行列は区別する。

### 10.3. 局所規則と返却表現

全225積・45作用を有限線型結合として直接返す。数学的に零の項は生成せず、恒等元との積は相手一項を返す。通常の評価は基底を要求しない。

| 入力の組 | 返却する primitive |
| --- | --- |
| I との積 | 相手をそのまま一項で返す |
| 両方が共通要素 I・Sz・Q0 | 共通部分 `[I,Sz,Q0]` に閉じる |
| 両方が CartesianBasis の要素 | CartesianBasis の要素で返す |
| 両方が LadderBasis の要素 | LadderBasis の要素で返す |
| Cartesian 専用要素と Ladder 専用要素の混在 | CartesianBasis の要素で返す |

返却順は対象基底の列挙順とし、局所積に行列単位への変換・密行列・記号計算を持ち込まない。係数は数学的に整理した最終形を選択した精度で生成する。

```math
[S_x,S_y]=icS_z,\quad \sum_\alpha S_\alpha^2=3I,\quad
S_z^2=I+Q_0/\sqrt2,\quad Q_0^2=I-Q_0/\sqrt2,
```

```math
S_\pm^2=cQ_{2\pm},\quad \{S_z,S_\pm\}=cQ_{1\pm},\quad S_\pm^3=0.
```

Sp・Sm の非零な遷移振幅は √3、Sx の隣接状態への振幅は √3/2。Sz の固有値は `c,0,−c` で、状態名の磁気量子数は物理的 Jz に対応する。随伴は Sp/Sm、Q1p/Q1m、Q2p/Q2m を交換し、他の9要素を保存する。通常 trace は I が3、他の14要素が0である。

### 10.4. 係数能力と共通 API

| 能力 | 係数に必要な能力 |
| --- | --- |
| Identity・Dimension・ExpressionFactory・保持 | なし |
| InnerProduct・Trace | `Zero + One` |
| Adjoint・SzBasis | `One` |
| CartesianBasis | `One + Neg + ImaginaryUnit` |
| LadderBasis | `One + Add + Div + Neg + ImaginaryUnit` |
| Product・Act | `Complex<R>`、`R: num_traits::Float` |

算術の Output は C とする。Product にも平方根係数が必要なので、当面は Float を用い、厳密な代数的係数のための新しい trait は追加しない。基底変換は有理係数と i だけを使い、保持・遅延構築・随伴・trace に Float を要求しない。平方根は R の One・算術・sqrt から生成し、f64 経由のキャストやパースは行わない。

共通の局所式・複合式・HS 内積・正規化・export に接続する。HS 内積は通常表現の局所 trace を Tr(I) で割ってから縮約し、総次元の整数化を要求しない。通常 trace と状態内積の意味は維持する。Float・std の既存 feature 条件は変えない。

### 10.5. 合成・export と grading

全状態・演算子は even、parity 演算子は I とする。ordinary / graded の両方で同じ値となり、磁気量子数の偶奇を fermion parity と解釈しない。二因子では9状態・81演算子基底、右端最速の mixed-radix 順を使う。作用行列は既存 export の列ごとの作用で生成し、primitive に backend 固有の行列化を実装しない。

### 10.6. 実装時の検証記録

規格化の移行と受入検証は [v3 の修正記録](../plans/v3/primitive-normalization.md)にまとめる。独立な物理的 J の行列と採用倍率から参照値を作り、[contracts.rs](../../crates/qmb-spin-one/tests/contracts.rs) で全225積・45作用・随伴・trace・結合則・作用整合性・両基底と Gram 行列を f32 / f64 で検証する。

[minimal_scalar.rs](../../crates/qmb-spin-one/tests/minimal_scalar.rs) は平方根を必要としない能力の最小 bound、[local_expressions.rs](../../crates/qmb-spin-one/tests/local_expressions.rs) は局所式と全 backend、[tensor_products.rs](../../crates/qmb-spin-one/tests/tensor_products.rs) は複合基底・ordinary / graded・全 backend、[models.rs](../../crates/qmb-spin-one/tests/models.rs) は物理的な3モデルとエネルギー・構築能力を検証する。

数値比較は複素数の絶対値による `|a−b| <= atol + rtol max(|a|,|b|)` とし、現在の小規模試験では f32 は atol=rtol=2e-5、f64 は atol=rtol=1e-12 を使う。検証用の閾値を通常評価の項削除に流用しない。局所規則の厳密な零と、一般の式の相殺残差は区別する。

### 10.7. 後続の範囲

公開の MatrixUnitBasis、標準的な球面テンソル位相を付けた基底、一般 spin-S、異種局所代数の混在は後続とする。Float を満たさない厳密係数による局所 Product / Act、専用 no_std target・MSRV の検証状況は別途扱う。
