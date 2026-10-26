# 局所 primitive の設計

更新日: 2026-10-27

## 1. 対象と検討状況

本書は局所状態・演算子 primitive の能力と整合性を定める。LocalStructure と追加の定義ラベルを廃止し、全 crate の式・空間・幾何・利用例を S・O の直接指定へ移行した。採用理由と依存調査は[調査記録](../refactoring-local-structure.md)、移行と検証は[実装計画](../implementation-plan.md#8-localstructure-の廃止)に記録する。

局所式は `LocalState<S, O, C>`・`LocalOperator<O, C>`、複合式は `State<I, S, O, C>`・`Operator<I, O, C>` とする。演算子式に状態型を要求しない。状態式は遅延した作用節点を保持するため O も指定する。保持・共有・係数計算を伴わない遅延構築に演算能力や primitive の Clone を要求しない。

## 2. 空間の情報と能力の関連付け

| 責任 | 実装した配置 |
| --- | --- |
| Primitive の意味・局所規則 | S・O と個別の能力実装。追加の定義 ID・ラベルは設けない |
| 局所状態空間の次元 | S の SpanDimension。primitive の値全体が張る空間の次元 d |
| 局所演算子空間の次元 | O の SpanDimension。全演算子空間では d²。variant の個数とは異なる |
| 恒等演算子 | `Identity::identity() -> Self` を O が実装する。係数型に依存しない |
| 作用の型対応 | 評価時の `S: Act<C, Operator = O, Term = S>` と `O: Adjoint<C, Term = O>` 等 |
| Scalar の整合性 | C と各操作の bound。異なる精度を暗黙に混ぜない |
| 基底 | `Basis<S, C>`・`Basis<O, C>`。列挙数は対象 primitive の span の次元と一致する |
| 因子構造・幾何 | `HilbertSpace<I>` が因子順序・grouping、Chain が Site ごとの Index 列と境界条件を保持する。Chain も移行済み |

SpanDimension は `fn span_dimension() -> Result<NonZeroUsize, DimensionOverflow>` として実装した。次元を整数化できない場合は取得時にエラーを返す。`state_basis_len<S>()`・`operator_basis_len<O>()` は対象の次元をそのまま返す。Operator の d² を再度二乗しない。能力は次元を使う操作だけが要求する。

恒等元は primitive の値として表せることを要求するが、選択基底の要素である必要はない。基底に含まれなければ Basis が展開する。線型従属な primitive を許容するため、恒等元用の値を追加しても基底の次元は増えない。式木の記号的な Identity 節点は、数値能力も恒等 primitive の生成も不要な遅延構築と、単位元の直接処理のため維持した。trace で必要になった時に O の恒等元を取得する。

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

一般の U→V の写像と State／Operator の型の統合は初期版実装後、Grading の情報保持・API は fermion 系の実装時に検討する。ordinary な空間の責任分離は実装済みである。Chain は Site ごとの Index 列と境界条件を直接保持し、Geometry 構造体と空間保持を廃止する方針を2026-10-26に採用した。Chain の実装移行は2026-10-27に実施し、契約は[Chain の設計](geometry.md)に従う。異種代数の型消去や実行時の表現パラメータを今回の移行には追加しない。

## 7. LocalStructure 廃止の採用方針（2026-10-26）

S・O の直接指定、primitive の span の次元と恒等元、追加の定義ラベルを設けない方針を採用し、実装へ反映した。API は第1〜2節、型と Index の責任分担は第4〜5節を現在の契約とする。旧 LocalStructure・Identity<C> と SpinHalf の定義用 marker は公開しない。
