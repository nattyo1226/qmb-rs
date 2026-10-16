# 局所 primitive・式木の数学的能力

更新日: 2026-10-15

## 1. 対象と責任分担

本書は、状態・演算子の primitive と、それらを葉とする局所式の数学的能力を定める。[要件](../requirements.md)と[アーキテクチャ概要](../architecture.md)に基づく詳細設計であり、[局所構造](local-structure.md)が対応する primitive と能力を関連付け、その整合性を要求する。

Primitive は局所的な演算規則を実装し、core の共通の式木・評価機構が線型結合、積、作用、随伴などへ拡張する。独自の局所代数を追加する際は、primitive と規則を実装し、共通の式木を利用する。具体的な Rust の型、trait シグネチャ、および内部表現は本書では定めない。Product と ActOn は能力を区別する設計上の名前として用いる。

本書の局所式は、一つの局所構造の上の式を扱う。Index の付与と全作用空間への束縛は、[Index の契約](index.md)および要件の 5.1 節に従う。複合空間の tensor product、因子間の grading、partial trace は共通の合成・埋め込み・評価機構が担当する。

複合空間は [Hilbert 空間の再帰構造](hilbert-space.md)で表す。局所複合代数の数学的能力は、この空間構成に沿って子の局所規則と共通の合成・評価機構から処理する。空間の木の構築・基本因子の列挙によって、状態・演算子の式を全展開することは要求しない。

局所式と Index の組を複合式の葉とし、複合 Hilbert 空間が子の局所契約を保存する能力の合成と整合性を担う。具体的な合成則・埋め込み・縮約の契約は、[埋め込み・複合空間の式と数学的能力](composite-expression.md)に定める。

## 2. Primitive の数学的能力

状態 primitive の線型結合は対象の全状態空間を、演算子 primitive の線型結合は対象の全演算子空間を張る。primitive 自体の線型独立性は要求せず、選択基底とは区別する。

| 実装先 | 能力 | 入力 | 結果 |
| --- | --- | --- | --- |
| 状態 primitive | 内積 | 状態 primitive 二つ | 複素 scalar |
| 演算子 primitive | Product | 演算子 primitive 二つ | 演算子 primitive の有限線型結合 |
| 演算子 primitive | ActOn | 演算子 primitive と状態 primitive | 状態 primitive の有限線型結合 |
| 演算子 primitive | 随伴 | 演算子 primitive | 演算子 primitive の有限線型結合 |
| 演算子 primitive | trace | 演算子 primitive | 複素 scalar |
| 状態・演算子 primitive | parity の取得〈graded の場合〉 | primitive | even / odd |

Product と ActOn は演算子 primitive 側の能力とする。ActOn の規則は、演算子 primitive 型と作用先の状態 primitive 型の組に対応する。対応する型と表現パラメータの整合性は局所構造の契約に従う。

Product・ActOn・随伴の結果は、同じ局所構造の primitive の有限線型結合とする。単一 primitive への閉性は要求しない。この線型結合も core の共通の式で表す。加算、scalar multiplication、および zero は共通の式側で提供する。

trace は局所 Hilbert 空間上の通常の非正規化 trace を返す。正規直交基底、作用、内積から数学的に導出できるが、通常の評価経路では primitive の trace を直接利用できるようにし、基底展開や数値行列の生成を必須としない。

Graded な primitive は homogeneous なものとし、混合 parity の対象は式の和で表す。局所構造は、局所 parity 演算子を得る手段と能力間の grading の整合性を定める。

## 3. 共通の式木による構築と評価

| 実装先 | 構築・提供する操作 | 評価時の担当 |
| --- | --- | --- |
| 状態式 | 加算、scalar multiplication、zero | 線型結合の処理 |
| 状態式 | 他の状態式との内積 | 半双線型性による拡張と、状態 primitive の内積への委譲 |
| 演算子式 | 加算、scalar multiplication、zero | 線型結合の処理 |
| 演算子式 | 積 | 双線型性と積の順序を扱い、必要な primitive 間の Product へ委譲 |
| 演算子式 | 状態式への作用 | 両引数の線型性を扱い、必要な primitive 間の ActOn へ委譲 |
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

これらの規則は、式の構築時に全展開する要求ではない。積・作用・随伴などを未評価の式として保持し、内積、trace、作用の結果、座標取得、または明示的な簡約結果が必要な時点で、必要な部分だけを処理する。primitive の規則はその評価の過程で呼び出す。

式木の節点構成や所有権の詳細は未決とする。未簡約の source expression と不変な部分構造を共有し、共通の完全展開形や評価結果のキャッシュを要求しない。局所式のために局所代数ごとの独自式木を導入しない。

## 4. 共通機構から導出する操作

| 対象 | 操作 | 導出に用いる能力 |
| --- | --- | --- |
| 状態式 | norm | 内積から `sqrt(⟨ψ,ψ⟩)` を計算する |
| 状態式 | 正規化 | norm と scalar multiplication を用いる。zero-norm state の正規化は拒否する |
| 状態式 | bra | Hilbert 空間の内積により ket の随伴として得る。ket から bra への対応は共役線型とする |
| 演算子式 | commutator / anticommutator | 積と加算から `AB−BA` / `AB+BA` を構築する |
| 演算子式 | Hilbert–Schmidt 内積 | `Tr(A†B)` を計算する。次元による除算は行わない |
| 演算子式・状態式 | 期待値 | 作用と内積から `⟨ψ,Aψ⟩` を計算する |
| 演算子式・状態式 | 正規化期待値 | 期待値を `⟨ψ,ψ⟩` で割る。zero-norm state では拒否する |

これらについて局所代数ごとの実装を必須としない。恒等演算子は、局所構造から対応する演算子式として取得する。

## 5. 線型従属性、構造比較、および基底との接続

Primitive が線型独立でないため、同じ数学的対象に複数の式表現を許す。たとえば `S+` と `Sx + i Sy` をどちらも記述できる。各能力は、primitive の線型関係を保存し、同じ対象を表す異なる入力表現から数学的に同じ結果を得ることを要求する。結果の係数列や式木の構造の一致は要求しない。

Source の構造比較では、同じ局所構造のもとで primitive の種類・パラメータ、式の構文、および係数の厳密な一致を判定する。局所構造の一致は、[同じ定義・定義パラメータによる契約](local-structure.md#4-局所構造の一致条件)に従う。Index を付与した葉では Index も比較する。通常の同類項の結合はこの構造的一致を用い、primitive 間の一般の線型関係を構造比較に含めない。

指定基底への正確な展開と座標取得は、基底との接続として提供する。そこでは primitive 間の線型関係も正しく扱う。基底選択を primitive に固定せず、基底展開を通常の作用・内積・trace・簡約の一般経路に組み込まない。数値的近似比較と作用としての等価性は、要件の 5.6 節に従う。

基底の定義を表す struct が `Basis<Primitive>` を実装し、primitive の基底展開、基底要素の順序番号の取得、および基底順序での列挙を提供する。基底要素は primitive の部分集合とし、要素数は対応する局所構造と対象種別から導出する。共通機構が線型性によって式全体へ拡張する。詳細は [基底の数学的能力](basis.md)に定める。

## 6. 整合性の検証対象と未決事項

各局所実装では、Product の結合則、Product と ActOn の整合性、恒等演算子、内積と随伴、通常の trace、および primitive の線型関係の保存を検証対象とする。数学的な契約は[局所構造の設計](local-structure.md)に定める。式側では、線型性・半双線型性・共役線型性の拡張、右側からの作用、随伴の積順序反転、および lazy 評価を検証する。

これらの等式は数学的な意味について要求し、浮動小数点計算のビット単位の一致とは区別する。数値精度、特殊値、および近似比較は既存要件の規約に従う。

具体的な Rust の型・trait シグネチャ、線型結合の返却形式、および式木の内部表現・所有権設計は未決とする。局所構造の一致契約は、[局所構造の一致条件](local-structure.md#4-局所構造の一致条件)に定める。
