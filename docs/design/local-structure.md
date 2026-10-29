# 局所 primitive の設計

更新日: 2026-10-30

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

Dimension は `fn dim() -> Result<NonZeroUsize, DimensionOverflow>` として実装した。次元を整数化できない場合は取得時にエラーを返す。`state_basis_len<S>()`・`operator_basis_len<O>()` は対象の次元をそのまま返す。Operator の d² を再度二乗しない。能力は次元を使う操作だけが要求する。

恒等元は primitive の値として表せることを要求するが、選択基底の要素である必要はない。基底に含まれなければ Basis が展開する。線型従属な primitive を許容するため、恒等元用の値を追加しても基底の次元は増えない。局所式の恒等元は `Primitive(O::identity())` として保持し、作用・積・随伴・trace・基底展開を通常の primitive と同じ経路で処理する。`LocalOperator::identity()` と空の `product_all` は `O: Identity` を要求するが、係数の数値能力は要求しない。複合式の記号的な `OpNode::Identity` は、未束縛の空の積や全作用空間の恒等演算子として維持する。

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

空間の次元は `space.dimension::<P>()` で指定した primitive の span を各因子に適用するか、`space.dimension_with(local_dimension)` で明示する。S を指定すると状態空間、O を指定すると演算子空間の次元を取得する。この取得操作は空間に primitive の定義を登録しない。空の tensor は次元1とし、通常の作用・内積・trace に総次元の整数化を要求しない。

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

`qmb-fermion` は `no_std + alloc`、`unsafe` 禁止とし、係数型から独立した `FermionState`・`FermionOperator` を提供する。状態は `Empty`・`Occupied`、演算子は `I`・`Annihilate`・`Create`・`Number`・`MajoranaX`・`MajoranaY`・`MajoranaZ`・`Parity` の8 variant とする。いずれも `Clone, Copy, Debug, PartialEq, Eq, Hash` を実装する。Zero は primitive に追加せず、共通の空の線型結合・式で表す。

`Fermion<C>` は `PhantomData<fn() -> C>` による型選択 marker とし、QuantumSystem の State・Operator・Coeff を関連付ける。構築は core の共通式を返し、係数精度の暗黙変換や独自 evaluator は追加しない。公開入口は [lib.rs](../../crates/qmb-fermion/src/lib.rs)、局所規則は [state.rs](../../crates/qmb-fermion/src/state.rs)・[operator.rs](../../crates/qmb-fermion/src/operator.rs)、基底は [basis.rs](../../crates/qmb-fermion/src/basis.rs) に実装する。Hubbard モデルの API と配置は [Chain・Site の設計](geometry.md#6-fermion-の-hubbard-モデル)に定める。

### 9.1. 規約・次元・parity

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

### 9.2. 局所能力

#### 9.2.1. Product の全規則

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

表から `a²=(a†)²=0`、`aa†+a†a=I`、`n=a†a`、`n²=n`、`P²=I`、`Pa=-aP`、`Pa†=-a†P` と、9.1節の su(2) の関係が得られる。P と a の掛け順、および `Z=−2iXY` の係数・符号を取り違えない。

#### 9.2.2. Act・内積

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

#### 9.2.3. Adjoint・Trace

`Adjoint<C, Term=FermionOperator>` は a と a† を交換し、I・n・P・X・Y・Z を保存する。primitive の係数は1とし、式の係数の共役は core が担当する。

| 演算子 | I | a | a† | n | P | X | Y | Z |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 通常の非正規化 trace | 2 | 0 | 0 | 1 | 0 | 0 | 0 | 0 |

`Tr(I)=2` とし、正規化 trace や supertrace に変更しない。局所 trace の実装に状態基底・行列を要求しない。

### 9.3. 基底と展開

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

### 9.4. 係数型と能力

各操作の実装に必要な trait bound を個別に要求する。既存 `LinearCombination::primitive` が `C: One` を要求することと、`One`・`Zero` 自体が持つ算術上の制約を前提とする。次は実装した各 impl の bound であり、式の保持へは伝播させない。全能力をまとめて一律に要求する係数 trait は導入せず、実際の bound は各操作で使う能力に合わせて定める。

| 能力 | C に要求する能力 |
| --- | --- |
| Identity、Dimension、parity 問い合わせ、QuantumSystem | 数値制約なし |
| InnerProduct、Trace | `Zero + One`。2は `one + one` で生成 |
| Adjoint、状態基底 | `One` |
| Product、Act、LadderBasis、MajoranaBasis | `One + Add<Output=C> + Div<Output=C> + Neg<Output=C> + ImaginaryUnit`。1/2・1/4・i を表せることも数学的契約。零は空の線型結合で返す |

局所規則には Float・平方根・C の Clone / Copy を要求しない。一方、Majorana 系の追加で虚数係数と分数が必要になるため、符号付き整数だけでは全局所規則を表せない。全 variant を対象にする trait impl の bound は、呼び出した値が a・n などでも同じである。独自の厳密な複素有理係数などで検証し、整数の切り捨て除算は使わない。共通式の評価・集約・export には、それぞれ core が要求する Clone・共役・算術能力が別途必要である。

精度は C が決め、`Complex<f32>` と `Complex<f64>` を暗黙に混ぜない。定数のために f64 経由のキャストやパースを追加しない。

### 9.5. Graded 合成・spinful への接続

局所 primitive の S・O は core の `HomogeneousParity` を実装する。O は `ParityOperator` も実装し、係数型に依存せず `FermionOperator::Parity` 自身を返す。局所の parity と P の契約は[局所能力](#8-grading-の局所能力と移行)を参照し、重複定義しない。

多 mode の順序・交換符号・先行区間の P は core の合成機構に委ねる。`FermionOperator` に Index・Site・mode 番号・parity string・global ordering を保持しない。既存 ordinary tensor へ複数の a を配置するだけでは fermion の異なる mode の CAR は得られないため、その構成を fermion 多体系の利用例として案内しない。

spinful の内部順序は `(up,down)`、二重占有の位相は `|up down>=f_up† f_down†|0>` とする。合成側との接続試験では、`f_up|up down>=|down>`、`f_down|up down>=-|up>`、総 parity `diag(1,-1,-1,1)` を確認する。表示順 `[empty, up, down, double]` は、二 mode の mixed-radix 順 `[empty, down, up, double]` と異なるので明示的に対応付ける。この表示順の違いを単一 mode primitive の符号へ埋め込まない。

### 9.6. 実装時の検証記録

2026-10-29、Rust 1.99.1 の direnv devShell で受入確認した。以下は実装時の検証記録であり、文書整理時の再実行結果ではない。局所代数と graded 接続・モデルは別々のテストで検証する。

| 検証ファイル | 確認した内容 |
| --- | --- |
| [contracts.rs](../../crates/qmb-fermion/tests/contracts.rs) | 独立な2×2行列による64積・16作用・8随伴/trace・4内積、512三重積の結合則、128組の作用整合性、随伴と内積、単位元・trace 巡回性、CAR・su(2)・Casimir・parity、両演算子基底の列挙・展開・再展開・Gram 行列 |
| [minimal_scalar.rs](../../crates/qmb-fermion/tests/minimal_scalar.rs) | 非 Clone / Copy / Float の厳密な複素有理係数で全局所積・作用・両基底を確認。One だけの係数による随伴・状態基底、整数による内積・trace、数値能力なしでの遅延構築、評価能力・算術能力のない係数と非 Clone Index による両モデル構築 |
| [local_expressions.rs](../../crates/qmb-fermion/tests/local_expressions.rs) | 複素係数の作用・期待値・共役・trace、norm・位相を保存する正規化・正規化期待値、展開・座標、全 primitive の行列 export と Vec / ndarray / sprs の一致 |
| [models.rs](../../crates/qmb-fermion/tests/models.rs) | 二 mode の CAR・spinful の消滅符号・総 parity、graded な基底再展開、独立な Fock 空間参照による0〜3 Site の両モデル、open / periodic、`t=0`・`U=0`、定数を含む両モデルの差、各 spin の粒子数保存、norm・期待値・全 backend、空間共有と不正入力 |

数値検証は `Complex<f32>` / `Complex<f64>` の両方で実行した。統合テスト13件と doctest 4件（異精度を混ぜる compile-fail を含む）が成功した。

- `cargo test -p qmb-fermion`
- `cargo test --workspace --all-features`（doctest を含む）
- `cargo check -p qmb-fermion --no-default-features`
- `cargo clippy --workspace --all-features --all-targets -- -D warnings`
- `cargo fmt --all -- --check`

利用例 [hubbard.rs](../../examples/src/bin/hubbard.rs) は `--all-features`、`--no-default-features`、`--no-default-features --features export` で実行し、2 Site の期待値0・−2と16×16の密疎行列の一致を確認した。

no-default-features の確認はホスト上のビルドである。専用 no_std target・MSRV 1.85 の実ビルド・性能ベンチマークは [未実施事項](../future-work.md#3-未実施の検証と性能)に残す。

## 10. Spin-one の局所代数

`qmb-spin-one` は三状態の spin-1 を実装する。`std` を許可し、unsafe を禁止する。公開入口は [lib.rs](../../crates/qmb-spin-one/src/lib.rs)、局所規則は [operator.rs](../../crates/qmb-spin-one/src/operator.rs)・[state.rs](../../crates/qmb-spin-one/src/state.rs)、基底は [basis.rs](../../crates/qmb-spin-one/src/basis.rs) に置く。XYZ・Heisenberg・BLBQ の定義と構築契約は [SpinOne のモデル構築](geometry.md#7-spinone-のモデル構築)にまとめる。

### 10.1. Primitive と規約

公開 variant とその列挙順は次のとおりである。

```rust
enum SpinOneState { Plus, Zero, Minus }

enum SpinOneOperator {
    I, Sx, Sy, Sz, SPlus, SMinus,
    Dxy, Q0, Qxy, Qyz, Qzx,
    Q1Plus, Q1Minus, Q2Plus, Q2Minus,
}
```

`hbar=1`、状態の順序を `(|+1>,|0>,|-1>)` とし、正規直交内積を使う。`SpinOneState::Zero` は磁気量子数0の正規化状態であり、zero vector は共通の空の線型結合で表す。

```math
S_z|m\rangle=m|m\rangle,\qquad
S_\pm|m\rangle=\sqrt{2-m(m\pm1)}|m\pm1\rangle,
\qquad S_\pm=S_x\pm iS_y.
```

端の状態への ladder 作用は zero とする。`r=√2` と書くと、参照行列は次である。

```math
S_z=\begin{pmatrix}1&0&0\\0&0&0\\0&0&-1\end{pmatrix},\quad
S_+=r\begin{pmatrix}0&1&0\\0&0&1\\0&0&0\end{pmatrix},\quad S_-=S_+^\dagger,
```

```math
S_x=\frac r2\begin{pmatrix}0&1&0\\1&0&1\\0&1&0\end{pmatrix},\quad
S_y=\frac r2\begin{pmatrix}0&-i&0\\i&0&-i\\0&i&0\end{pmatrix}.
```

反交換子を `{A,B}=AB+BA` とし、全演算子を以下で定義する。

| Primitive | 定義 |
| --- | --- |
| `I` | 恒等演算子 |
| `Sx`, `Sy`, `Sz` | スピンの Cartesian 3成分 |
| `SPlus`, `SMinus` | `S±=Sx±iSy` |
| `Dxy` | `Sx²−Sy²` |
| `Q0` | `3Sz²−2I` |
| `Qxy`, `Qyz`, `Qzx` | それぞれ `{Sx,Sy}`, `{Sy,Sz}`, `{Sz,Sx}` |
| `Q1Plus`, `Q1Minus` | `Q1±={Sz,S±}=Qzx±iQyz` |
| `Q2Plus`, `Q2Minus` | `Q2±=S±²=Dxy±iQxy` |

恒等元1個、スピン5値、四重極9値の計15値を持つ。Cartesian の四重極5値と昇降型4値は同じ5次元空間を張る。対応する構成は[双極子・四重極による spin-one の構成（式3–7）](https://arxiv.org/pdf/1305.2560)と[球面テンソルによる多重極展開（式29–30）](https://arxiv.org/pdf/1409.1106)を参照し、規格化・位相は上表の定義を使う。

両 primitive は `Clone, Copy, Debug, PartialEq, Eq, Hash` を実装し、構造的な等値性は variant を区別する。Dimension は状態が3、演算子が9であり、15という variant 数とは異なる。`Identity::identity()` は I を返し、係数の数値能力を要求しない。

### 10.2. 基底と展開

| 基底 | 列挙順 |
| --- | --- |
| `SzBasis` | `[Plus,Zero,Minus]` |
| `CartesianBasis` | `[I,Sx,Sy,Sz,Dxy,Q0,Qxy,Qyz,Qzx]` |
| `LadderBasis` | `[I,SPlus,SMinus,Sz,Q0,Q1Plus,Q1Minus,Q2Plus,Q2Minus]` |

両演算子基底は完全・線型独立で、Hilbert–Schmidt 内積 `Tr(A†B)` に関して直交する。正規直交化はしない。CartesianBasis のノルム二乗は I が3、Q0 が6、他が2、LadderBasis では I が3、Sz が2、Q0 が6、他が4となる。

CartesianBasis の全要素は Hermite である。LadderBasis の要素 B は `[Sz,B]=qB` を満たし、q は I・Sz・Q0 では0、S±・Q1± では±1、Q2± では±2となる。回転に対して I は rank 0、スピンは rank 1、四重極は rank 2に対応するが、標準的な球面テンソルの規格化・位相係数は付けない。

両演算子基底の和集合が15個の公開 primitive となる。各 `elements()` は選択した基底の要素だけを上表の順に列挙し、`index` は非メンバーに `None` を返す。SzBasis は状態をそのまま展開する。演算子基底のメンバーも自身一項に展開し、他の要素は次の相互変換で展開する。

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

I・Sz・Q0 は両基底で共有する。通常の局所評価は基底を要求せず、基底選択は展開・座標取得時に明示する。9個の演算子基底係数と SzBasis 上の3×3作用行列、および局所基底番号と export の格納順序を区別する。

### 10.3. 局所規則と返却表現

`Product<C,Term=SpinOneOperator>` の返却規則は次の優先順とする。

| 入力の組 | 返却する primitive |
| --- | --- |
| I との積 | 相手をそのまま一項で返す |
| 両方が共通要素 I・Sz・Q0 | 共通部分 `[I,Sz,Q0]` に閉じる |
| 両方が CartesianBasis の要素 | CartesianBasis の要素で返す |
| 両方が LadderBasis の要素 | LadderBasis の要素で返す |
| Cartesian 専用要素と Ladder 専用要素の混在 | CartesianBasis の要素で返す |

返却順は対象基底の列挙順、項数は最大9個とする。入力 primitive の値から決め、式へ基底情報を保持させたり Product に基底引数を追加したりしない。各基底内の線型結合同士の積も同じ基底内に閉じる。数学的に同じ結果の式構造の一致は要求しない。

局所積の全225組と作用の全45組を、既存代数と同様に `match` で直接記述する。各分岐で有限線型結合を生成し、数学的に zero の項を含めず、zero は空の線型結合で返す。浮動小数点の近似零判定で項を削除しない。行列単位 `E_ab=|a><b|` と `E_ab E_cd=δ_bc E_ad` は導出・参照に使えるが、製品の実行時には行列単位への分解・逆変換、密行列、独自の有理数型や記号係数の演算機構を使わない。

代表的な規則は以下である。

```math
\begin{aligned}
S_z^2&=\tfrac23 I+\tfrac13 Q_0,& S_zQ_0&=Q_0S_z=S_z,& Q_0^2&=2I-Q_0,\\
S_xS_y&=\tfrac12(Q_{xy}+iS_z),&
S_+^2&=Q_{2+},& S_-^2&=Q_{2-},\\
S_+S_-&=\tfrac43I-\tfrac13Q_0+S_z,&
S_-S_+&=\tfrac43I-\tfrac13Q_0-S_z,& S_+^3&=S_-^3=0.
\end{aligned}
```

`Act<Complex<R>>` は `Operator=SpinOneOperator,Term=SpinOneState` とする。`r=√2` として、I は状態を保存し、スピンの作用は次のとおりである。

| 演算子 | Plus | Zero | Minus |
| --- | --- | --- | --- |
| Sz | Plus | 0 | −Minus |
| SPlus | 0 | r Plus | r Zero |
| SMinus | r Zero | r Minus | 0 |
| Sx | (r/2) Zero | (r/2)(Plus + Minus) | (r/2) Zero |
| Sy | i(r/2) Zero | −i(r/2) Plus + i(r/2) Minus | −i(r/2) Zero |

表中の `0` は zero vector、`Zero` は m=0 の状態である。InnerProduct は3状態の Kronecker delta とし、式の半双線型性と係数の共役は core に委ねる。

Adjoint は SPlus / SMinus、Q1Plus / Q1Minus、Q2Plus / Q2Minus を交換し、追加の符号を付けない。他の9値は自身を返す。Trace は通常の非正規化 trace とし、Tr(I)=3、他の14値は0を直接返す。

### 10.4. 係数能力と共通 API

`SpinOne<C>` は QuantumSystem に状態・演算子・係数を関連付け、保持・係数計算を伴わない遅延構築には数値 bound を要求しない。primitive の能力は操作ごとに限定する。

| 能力 | C に要求する能力 |
| --- | --- |
| Identity・Dimension・QuantumSystem、primitive 保持 | なし |
| InnerProduct・Trace | `Zero + One`。3は `one + one + one` |
| Adjoint・SzBasis | `One` |
| CartesianBasis | `One + Neg + ImaginaryUnit` |
| LadderBasis・Product | `One + Add + Div + Neg + ImaginaryUnit` |
| Act | `C=Complex<R>`、`R: num_traits::Float` |

Add・Div・Neg はそれぞれ `Output=C` とする。One と Zero は num-traits の trait であり、それぞれ乗算・加算能力を含む。係数は非零による数学的な除算に閉じ、必要な1/2・1/3・iを表せることを要求する。Rust の Div だけでは整数の切り捨て除算も通るため、この数学的な契約も適用する。浮動小数点では精度に応じた近似を許容する。

1/2・1/3・1/6等の定数は C の One・加算・除算、虚係数は ImaginaryUnit・乗算・符号反転で生成し、private な [coefficient.rs](../../crates/qmb-spin-one/src/coefficient.rs) で共有する。局所 Product は C の Clone / Copy / Zero / Float / Num / 型変換を要求しない。全 variant を扱う impl の bound は、実際に渡す primitive の値によらず同じである。現在の Basis は expand と metadata が同じ impl にあるため、index / elements にもその基底の bound が適用される。

作用は `Complex<f32>` / `Complex<f64>` を標準の利用対象とし、√2 は R の精度で次のように生成する。

```rust
fn sqrt_two<R: num_traits::Float>() -> R {
    (R::one() + R::one()).sqrt()
}
```

helper は private とし、公開の SqrtTwo trait は追加しない。実係数の埋め込みは `Complex::new(value,R::zero())` を使い、f64 経由のキャストやパースを行わない。num-traits 0.2.19 の Float が要求する R の Copy / Clone 等は許容する。Product・基底変換・primitive の内積・随伴・trace には Float を要求しない。f32 / f64 自体は虚数単位を表せず、全 variant を扱う Product / Act の係数には使えない。

共通式の能力は [primitive / 式の設計](primitive-expression.md)に従う。局所状態式の evaluate は C の Clone・Conjugate と primitive の Act / Adjoint、内積・期待値・trace・基底展開・座標取得はさらに必要な Zero・乗算等を要求する。primitive の bound を一律に強めない。norm squared は RealPart、norm は実数成分の Float、normalized は実係数の From も使い、正規化 API は qmb-normalize が提供する。

本番依存は既存 workspace の qmb-core・qmb-geometry・thiserror・num-traits・num-complex とする。num-traits の std feature をこの crate で有効にする。数値 backend は本番依存に含めず、core の公開契約や workspace 全体の既定 feature は変更しない。no_std 対応は初期要件に含めない。

### 10.5. 合成・export と grading

共通の LocalState / LocalOperator・State / Operator を使い、評価・期待値・随伴・trace・展開・座標取得・正規化へ接続する。作用行列は既存 export の列ごとの作用で生成し、primitive に backend 固有の行列化を実装しない。ordinary 二因子は9状態・81演算子基底を持ち、既存の mixed-radix 順序を使う。

全状態・演算子を even、P=I として HomogeneousParity / ParityOperator を明示実装する。磁気量子数 m の整数としての偶奇を fermion parity と解釈しない。符号は[ordinary / graded の共通評価](composite-expression.md#12-ordinary--graded-の共通評価と移行)に委ね、自明な grading により両方で同じ数値結果を得る。

### 10.6. 実装時の検証記録

2026-10-30、Rust 1.99.1 の direnv devShell で受入確認した。以下は実装時の検証記録であり、文書整理時の再実行結果ではない。モデル追加後は spin-one の統合テスト25件と doctest 4件（異精度の混在を拒否する compile-fail を含む）、workspace の全231テストが成功した。

参照計算は直接記述したスピン行列と四重極の定義から作り、製品の match や係数 helper を使わない。両精度で、複素数の絶対値による `|a−b| <= atol + rtol * max(|a|,|b|)` を使い、f32 は atol=rtol=2e-5、f64 は atol=rtol=1e-12とする。

| 検証ファイル | 確認した内容 |
| --- | --- |
| [contracts.rs](../../crates/qmb-spin-one/tests/contracts.rs) | 全225積・45作用・15随伴 / trace・9内積、基底内閉性・返却順・零項、3,375三重積・675作用整合性、単位元・随伴の反準同型性と内積整合性・trace 巡回性、交換則・Casimir `Sx²+Sy²+Sz²=2I`・`Tr(Sα²)=2`・ladder 二乗 / 三乗、両基底の列挙・index・展開・往復・座標・Gram 行列、公開型・次元・trivial grading |
| [minimal_scalar.rs](../../crates/qmb-spin-one/tests/minimal_scalar.rs) | 非 Clone / Copy / Zero / Float / Num / NumCast / Conjugate の厳密な複素有理係数で全225積と全15値の両基底展開。加算・除算のない係数で CartesianBasis、One だけで随伴・SzBasis、整数で primitive 内積・trace、数値能力なしで恒等元・遅延構築・配置 |
| [local_expressions.rs](../../crates/qmb-spin-one/tests/local_expressions.rs) | 複素線型結合の作用・内積・期待値・trace・随伴・積、norm / 正規化 / 正規化期待値・zero-norm 拒否、全15値の3×3行列と Vec / ndarray / sprs、状態と演算子の座標の区別 |
| [tensor_products.rs](../../crates/qmb-spin-one/tests/tensor_products.rs) | 二因子の9状態・両81演算子基底、mixed-radix 順序、全列の作用・trace・座標、spectator の I、複素係数の積・随伴・作用・期待値・再展開・正規化、ordinary / graded の一致、全 backend |
| [models.rs](../../crates/qmb-spin-one/tests/models.rs) | 3モデルの独立な参照行列、全スピン別エネルギー、Site・境界・空間共有・不正入力・最小係数能力。詳細は[モデルの検証](geometry.md#74-検証と実装時の記録) |

以下の確認が成功した。

- `cargo test -p qmb-spin-one`（局所 primitive の初回受入時）
- `cargo test -p qmb-spin-one --test models`
- `cargo test --workspace --all-features`
- `cargo check -p qmb-spin-one`
- `cargo clippy --workspace --all-features --all-targets -- -D warnings`
- `cargo fmt --all -- --check`
- `RUSTDOCFLAGS='-D warnings' cargo doc -p qmb-spin-one --no-deps`

利用例 [spin_one.rs](../../examples/src/bin/spin_one.rs) は `--no-default-features` と `--all-features` で実行した。局所 Hamiltonian の未正規化期待値1.2・正規化期待値0.4、SPlus² の Q2Plus 座標、二因子の期待値−1、3×3とモデルの9×9の密疎行列の一致を確認した。MSRV 1.85 の実ビルド・性能ベンチマークは[未実施事項](../future-work.md#3-未実施の検証と性能)に残す。

### 10.7. 後続の範囲

Gell-Mann 基底、公開の MatrixUnitBasis、正規直交化した球面テンソル基底、一般 spin-S、異種局所代数の混在は後続とする。Float を満たさない厳密な実数や非 Copy の実数、Complex<R> 以外の係数への作用も、標準の Act 実装の対象外とする。これらを初期の局所 primitive・2演算子基底・3モデルの完成条件には含めない。
