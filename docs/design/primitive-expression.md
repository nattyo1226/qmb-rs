# 局所 primitive・式木の数学的能力

更新日: 2026-10-27

## 1. 対象と責任分担

本書は、状態・演算子の primitive と、それらを葉とする局所式の数学的能力を定める。[要件](../requirements.md)と[アーキテクチャ概要](../architecture.md)に基づく詳細設計であり、[局所構造](local-structure.md)が対応する primitive と能力を関連付け、その整合性を要求する。

具体的代数の crate が primitive の局所的な演算規則を実装し、`qmb-core` の共通の式木・評価機構が線型結合、積、作用、随伴などへ拡張する。独自の局所代数を追加する際は、primitive と規則を実装し、core の共通型を直接利用する。具体的な Rust の型、trait シグネチャ、および内部表現は本書では定めない。Product と Act は能力を区別する設計上の名前として用いる。

本書の局所式は、一つの局所構造の上の式を扱う。Index の付与と全作用空間への束縛は、[Index の契約](index.md)および要件の 5.1 節に従う。複合空間の tensor product、因子間の grading、partial trace は共通の合成・埋め込み・評価機構が担当する。

初期版では、複合式へ配置した後も式の primitive 型 S・O を揃える。独自代数は独自の S・O と各能力実装により別の系で単独利用できるが、一つの系への異種局所構造の混在・合成は後続へ延期する。局所式の取り込みに異種構造の型消去を要求しない。

複合空間は [Hilbert 空間の再帰構造](hilbert-space.md)で表す。局所複合代数の数学的能力は、この空間構成に沿って子の局所規則と共通の合成・評価機構から処理する。空間の木の構築・基本因子の列挙によって、状態・演算子の式を全展開することは要求しない。

局所式と Index の組を複合式の葉とし、複合 Hilbert 空間が子の局所契約を保存する能力の合成と整合性を担う。具体的な合成則・埋め込み・縮約の契約は、[埋め込み・複合空間の式と数学的能力](composite-expression.md)に定める。

## 2. Primitive の数学的能力

状態 primitive の線型結合は対象の全状態空間を、演算子 primitive の線型結合は対象の全演算子空間を張る。primitive 自体の線型独立性は要求せず、選択基底とは区別する。

| 実装先 | 能力 | 入力 | 結果 |
| --- | --- | --- | --- |
| 状態 primitive | 内積 | 状態 primitive 二つ | 複素 scalar |
| 演算子 primitive | Product | 演算子 primitive 二つ | 演算子 primitive の有限線型結合 |
| 状態 primitive | Act | 状態 primitive と演算子 primitive | 状態 primitive の有限線型結合 |
| 演算子 primitive | 随伴 | 演算子 primitive | 演算子 primitive の有限線型結合 |
| 演算子 primitive | trace | 演算子 primitive | 複素 scalar |
| 状態・演算子 primitive | parity の取得〈graded の場合〉 | primitive | even / odd |

Product は演算子 primitive 側、Act は状態 primitive 側の能力とする。Act の規則は、作用先の状態 primitive 型と演算子 primitive 型の組に対応し、`state.act(&operator)` は数学的には `operator * state` を計算する。対応する型と表現パラメータの整合性は局所構造の契約に従う。

Rust API は状態式を S・O・C、演算子式を O・C で指定する。恒等元は O が Identity として生成し、次元は S・O が SpanDimension で提供する。演算能力と入出力の対応は必要な評価メソッドで要求し、保持型に全演算能力を要求しない。LocalStructure と追加の定義ラベルは導入しない。

Product・Act・随伴の結果は、同じ局所構造の primitive の有限線型結合を表す薄い `LinearCombination` とする。単一 primitive への閉性は要求しない。演算 trait は係数全体 C を型引数、返却要素を関連型 `Term` で指定し、局所式や局所構造には依存しない。局所規則では `Term` を対応する primitive 型とする。式木側は `map` で要素を変換し、和・スカラー倍として取り込む。恒等元は `O: Identity` の値として生成し、選択基底への展開は Basis が担う。未評価の作用・積・随伴を返して局所計算を先送りしない。加算、scalar multiplication、および zero は共通処理で提供する。演算 trait の関連型への移行と式木への接続は Step 2 で行う。

trace は局所 Hilbert 空間上の通常の非正規化 trace を返す。正規直交基底、作用、内積から数学的に導出できるが、通常の評価経路では primitive の trace を直接利用できるようにし、基底展開や数値行列の生成を必須としない。

Graded な primitive は homogeneous なものとし、混合 parity の対象は式の和で表す。局所構造は、局所 parity 演算子を得る手段と能力間の grading の整合性を定める。

## 3. 共通の式木による構築と評価

| 実装先 | 構築・提供する操作 | 評価時の担当 |
| --- | --- | --- |
| 状態式 | 加算、scalar multiplication、zero | 線型結合の処理 |
| 状態式 | 他の状態式との内積 | 半双線型性による拡張と、状態 primitive の内積への委譲 |
| 演算子式 | 加算、scalar multiplication、zero | 線型結合の処理 |
| 演算子式 | 積 | 双線型性と積の順序を扱い、必要な primitive 間の Product へ委譲 |
| 状態式 | 演算子式の作用を受ける | `act` で遅延した作用を構築し、必要な評価時に両引数の線型性を扱って primitive の Act へ委譲 |
| 演算子式 | 随伴 | 係数の共役、積の順序反転、primitive の随伴への委譲 |
| 演算子式 | trace | 線型性による拡張と、primitive の trace への委譲 |

以下で、`A,B,C` は局所演算子式、`ψ` は局所状態式、`a,b` は同じ精度の複素 scalar とする。共通の式処理は、次の拡張規則を扱う。

```math
\left\langle\sum_i a_i s_i,\sum_j b_j t_j\right\rangle
=\sum_{i,j}\overline{a_i}b_j\langle s_i,t_j\rangle,
\qquad
(aA+bB)\psi=a(A\psi)+b(B\psi).
```

```math
(aA+bB)C=a(AC)+b(BC),\qquad
A(aB+bC)=a(AB)+b(AC),
\qquad
(AB)\psi=A(B\psi).
```

```math
(aA+bB)^\dagger=\overline a A^\dagger+\overline b B^\dagger,
\qquad
(AB)^\dagger=B^\dagger A^\dagger,
\qquad
\operatorname{Tr}(aA+bB)=a\operatorname{Tr}(A)+b\operatorname{Tr}(B).
```

内積は第一引数について共役線型、第二引数について線型とする。積は両引数について線型、作用は演算子・状態の両引数について線型、随伴は共役線型、trace は線型とする。`AB` は状態に右側から作用し、`B` が先に、`A` が後に作用する。

これらの規則は、式の構築時に全展開する要求ではない。積・作用・随伴などを未評価の式として保持し、内積、trace、作用の結果、または指定基底への展開・座標取得が必要な時点で、必要な部分だけを処理する。primitive の規則はその評価の過程で呼び出す。

状態式の `state.act(&operator)` は `Act` 節点で作用を表す遅延式を返し、`state.evaluate()` は状態式に含まれる作用を評価して通常の状態式を返す。primitive の Act はその局所計算を担当する。`ψ.act(&b).act(&a)` は `A(Bψ)` を表す。`evaluate` は基底を指定せず、局所規則と線型性・合成則だけで処理し、選択基底への展開や座標化・行列化は行わない。

初期版では局所規則による `evaluate` と、指定基底への展開・係数取得を提供する。分配・局所積などの特定の簡約ルールのみを適用する独立した操作は future work とする。評価や基底展開に必要な内部処理として、分配・局所積・係数の結合を行うことは維持する。

局所式は `ScaleExpr`・`SumExpr`・`ProductExpr`・`ActExpr`・`AdjointExpr` を各節点の内容を保持する共通パーツとして用いる。zero は空の和で表す。`LinearCombination<P, C>` は `SumExpr<ScaleExpr<P, C>>` で保持する薄い有限線型結合型として維持する。局所規則や局所基底展開の結果を、式木側が所有権付きの `map` と共通 constructor で取り込む。係数・項順序・重複・零係数を保ち、この変換自体に要素や係数の `Clone` は要求しない。

共通パーツは `expr::parts` で保持内容を含めて公開し、式ハンドルは crate root から再公開する。式の節点 enum と共有の内部表現は非公開とする。項型 `Term<P, C>` は `ScaleExpr<P, C>` の別名に統合し、要素名を `operand` とする。既存の線型結合の構築・走査 API を維持し、走査を `SumExpr` へ集約する。演算 trait の関連型 `Term` は係数を含まない要素型を指す。

評価器は状態評価・演算子の作用・演算子項の走査を分担する。状態の作用評価に演算子の Product 能力は要求せず、局所積が必要な trace などに限定する。型の保持・係数計算を伴わない遅延構築・map に数値制約を付けず、必要な制約は各メソッドに置く。全能力をまとめた係数 trait は導入しない。式同士の二項演算子は所有値・借用値の4通りを提供し、所有ハンドルを移動、借用ハンドルは共有参照を複製する。名前付きメソッドは共有借用を基本とする。

評価と最後の所有者になった式節点の破棄は明示的な作業スタックで反復的に行い、式の深さに比例した関数再帰を避ける。長い積・作用のチェーンと、共有先が残る部分式を検証する。性能上の合格条件は共有、不要な primitive 複製の抑制、空・単一項のコンテナ無確保とし、map 全般の無確保や特定の実行時間は保証しない。

式ハンドルは `Arc` で不変な節点を共有し、節点内の primitive もその共有に含まれる。同じ式ハンドルを複数箇所で利用しても木全体や primitive を複製しない。別々に構築した同形の式を自動検出して共有する機構は要求しない。操作は入力を変更せず、共通の完全展開形や評価結果のキャッシュを要求しない。source の自動保持、式の構造比較・表示は後続へ延期する。局所式のために局所代数ごとの独自式木を導入しない。

## 4. 共通機構から導出する操作

| 対象 | 操作 | 導出に用いる能力 |
| --- | --- | --- |
| 状態式 | norm | 内積から `sqrt(⟨ψ,ψ⟩)` を計算する |
| 状態式 | 正規化 | norm と scalar multiplication を用いる。zero-norm state の正規化は拒否する |
| 状態式（後続） | bra | Hilbert 空間の内積により ket の随伴として得る。ket から bra への対応は共役線型とする |
| 演算子式 | commutator / anticommutator | 積と加算から `AB−BA` / `AB+BA` を構築する |
| 演算子式 | Hilbert–Schmidt 内積 | `Tr(A†B)` を計算する。次元による除算は行わない |
| 演算子式・状態式 | 期待値 | 作用と内積から `⟨ψ,Aψ⟩` を計算する |
| 演算子式・状態式 | 正規化期待値 | 期待値を `⟨ψ,ψ⟩` で割る。zero-norm state では拒否する |

bra 型・構築 API は後続へ延期し、初期版の braket は内積で計算する。Norm・正規化・正規化期待値は共通機構から提供し、局所代数ごとの実装を必須としない。恒等演算子は遅延構築では `Identity` 節点で保持し、trace で必要になった時に `O::identity()` を取得する。

### 4.1. Norm・正規化の API と feature

係数全体を C、`C: RealPart` の実数型を R とする。内積・期待値・norm squared・スカラー除算は `qmb-core`、以下の拡張 trait と `NormalizationError` は std を利用する `qmb-normalize` が提供する。

| 操作 | 提供元 | 局所状態の返り値 | 複合状態の返り値 |
| --- | --- | --- | --- |
| `norm_squared()` | core の状態式 | R | `Result<R, ExpressionError>` |
| `norm()` | `Norm` | R | `Result<R, ExpressionError>` |
| `normalized()` | `Normalize` | `Result<Self, NormalizationError>` | 同左 |
| `normalized_expectation(&operator)` | `NormalizedExpectation` | `Result<C, NormalizationError>` | 同左 |
| `normalized_expectation_in(&operator, &ambient)` | `NormalizedExpectationIn` | 提供しない | `Result<C, NormalizationError>` |

norm squared は自己内積の実部、norm はその平方根とする。`norm()`・`normalized()` だけが R に `num_traits::Float` を要求し、正規化には実数を C に埋め込む `From<R>` も要求する。正規化は入力を変更せず、位相を保つ実数の逆ノルムで式をスカラー倍する。通常の作用・内積と同じ評価経路を使い、基底展開・行列化は要求しない。

正規化と正規化期待値は、自己内積全体の `Zero::is_zero()` による厳密なゼロを `NormalizationError::ZeroNorm` として拒否する。明示的な zero に加え、非独立な primitive や項の相殺による zero も対象とする。`norm(zero)` は0を返す。負値の補正・近似ゼロ判定は導入せず、NaN・infinity などは数値型に委ねる。複合式の空間・評価エラーは `NormalizationError::Expression(ExpressionError)` に変換し、元のエラーを source として保持する。

正規化期待値は期待値を自己内積全体で除し、`Divide` を要求する。平方根・実部取得・実数の埋め込みは要求しない。標準 Complex の除算では、純実の分母を不要に二乗しない。明示空間を渡す複合式では、その空間を演算子の既存束縛より優先する。

拡張 trait とエラーは `qmb-normalize` から直接利用する。core の libm feature は設けず、core の代数 API・norm squared・スカラー除算は正規化 crate がなくても利用できる。非零の自己内積による商を明示的に計算でき、平方根を使わない正規化期待値の拡張 trait も `qmb-normalize` に属する。

## 5. 線型従属性と基底との接続

Primitive が線型独立でないため、同じ数学的対象に複数の式表現を許す。たとえば `S+` と `Sx + i Sy` をどちらも記述できる。各能力は、primitive の線型関係を保存し、同じ対象を表す異なる入力表現から数学的に同じ結果を得ることを要求する。結果の係数列や式木の構造の一致は要求しない。

評価に primitive や式の構造比較を必須としない。指定基底への展開では、その基底の順序番号を用いて同じ基底要素の係数をまとめる。局所構造や Index の整合性を検証する契約は、式の構造比較を延期しても維持する。

指定基底への正確な展開と座標取得は、式を整理する明示的な操作と、基底との接続として提供する。そこでは primitive 間の線型関係も正しく扱い、同じ基底要素の係数をまとめる。基底選択を primitive に固定せず、通常の作用・内積・trace に基底展開を必須としない。数値的近似比較と作用としての等価性は、要件の 5.6 節に従う。

基底の定義を表す struct が対象の要素型を指定した `Basis<P, C>` を実装し、primitive の基底展開、基底要素の順序番号の取得、および基底順序での列挙を提供する。状態 primitive の展開は `LinearCombination<S, C>`、演算子 primitive の展開は `LinearCombination<O, C>` を返し、各要素を選択した基底の primitive に限る。式木側で `map` して通常の式へ取り込む接続も Step 2 に含める。基底要素は primitive の部分集合とし、要素数は対応する局所構造と対象種別から導出する。共通機構が線型性によって式全体へ拡張する。詳細は [基底の数学的能力](basis.md)に定める。

基底は展開操作への入力であり、結果は基底識別情報や特別な展開済み状態を要求しない通常の式木とする。入力の式は変更せず、結果は局所構造などを保持し、通常の演算と別基底への再展開を適用できる。式木のための独立した基底変更機構は設けない。

## 6. 整合性の検証対象と未決事項

各局所実装では、Product の結合則、Product と Act の整合性、恒等演算子、内積と随伴、通常の trace、および primitive の線型関係の保存を検証対象とする。数学的な契約は[局所構造の設計](local-structure.md)に定める。式側では、線型性・半双線型性・共役線型性の拡張、右側からの作用、随伴の積順序反転、および lazy 評価を検証する。

これらの等式は数学的な意味について要求し、浮動小数点計算のビット単位の一致とは区別する。数値精度、特殊値、および近似比較は既存要件の規約に従う。

局所構造と演算能力の分離、共通パーツの公開と項型の統合、薄い `LinearCombination` と関連型 `Term` による規則の出力、Step 2 での `map` と式木への接続は確定事項とする。`StateBasis`・`OperatorBasis` の削除と、利用側での `Basis<S, C>`・`Basis<O, C>` の直接指定、メソッドごとの制約、所有権、反復的な評価・破棄、および `Arc` による不変共有と `act`・`evaluate` の役割も採用済みとする。非公開の作業領域・補助関数と、具体的な数値演算の呼び出し位置に必要な制約は実装時に定める。線型結合の保持・変換と選択基底の条件、共有によって deep copy が生じないことも検証対象に含める。局所構造の一致契約は、[局所構造の一致条件](local-structure.md#4-局所構造の一致条件)に定める。

## 7. 式構築の入口と係数変換（2026-10-28）

構築・集約の基本の入口は qmb_core::QuantumSystem とし、構築済みの式の操作・計算はインスタンスメソッドで行う。LocalStateOf・LocalOperatorOf・StateOf・OperatorOf は主に型注釈で使う。既存の core constructor は維持し、QuantumSystem の既定メソッドと同じ基礎処理へ委譲する。具体的な配置・合成の契約は[複合式の設計](composite-expression.md#11-quantumsystem-の構築入口と-fallible-集約2026-10-28)に記録する。

局所・複合の4式型と LinearCombination の scaled は `impl Into<C>` を受け取り、入口で一度だけ C に変換する。式型は係数を保持するだけで数値能力・Clone を要求せず、LinearCombination は従来の `C: Clone + Mul<Output=C>` を維持する。変換元の Clone は不要で、空の LinearCombination でも変換は一度行う。C の精度は維持し、新たなキャストや精度変換は導入しない。Mul・scalar・term の引数は広げない。

旧 scaled(C) の期待型に依存していたコードでは、構築入口の係数型や Default の入力型を明示し、不要な into() を外す。推論の完全互換は条件とせず、選択済みの C と最小能力の契約を維持する。[移行例](../plans/expression-construction.md#23-型推論の変化と移行例)を参照する。

LocalState::try_sum・LocalOperator::try_sum は `Result<式,ExpressionError>` の列を集約し、最初の入力 Err 以降を取り出さずに、そのまま返す。全成功時は既存 sum と同じ項順・重複・零係数・空／単項の扱いを保ち、primitive を評価しない。QuantumSystem の局所式用転送メソッドは後続候補であり、今回追加していない。
