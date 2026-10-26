# 係数型の依存・能力の整理（一時）

更新日: 2026-10-25

## 1. 基本方針と調査範囲

係数型には **「必要な時に必要な trait bound を課す」**。これは合意した基本方針である。保持・構築・評価・導出操作を区別し、全数値能力をまとめた trait をすべての型や操作の必須条件にしない。

本書は [リファクタリングの整理メモ](refactoring-notes.md)の数値部分を具体化する。現行コードと、Cargo.lock で解決されている num-traits 0.2.19・num-complex 0.4.6 のローカルソースを調査した。第2〜8節は移行前の調査記録である。調査後に実装を移行し、[API・検証結果](refactoring-coefficient-proposal.md#10-実装した-api-と移行方法)に記録した。

**採用した方針**：2026-10-25 の利用者の指示により、core には num-complex・num-traits を含め、Zero・One は num-traits をそのまま使い、不足する trait を Complex<f32>・Complex<f64> に直接実装する。係数全体 C の保持と個別能力、std を利用する qmb-normalize の拡張 trait、標準 Complex の左スカラー倍の個別実装、正規化期待値の純実分母の処理の維持も採用済みである。以前の数値依存解消や Zero・One の細分化を目標とした比較は調査経緯として扱う。採用事項の詳細は [係数型の採用方針](refactoring-coefficient-proposal.md)、num-traits の使い分けは第14節に記録する。

数値ライブラリの再公開は、その後の整理で削除した。現在は実装・テスト・doctest・利用例のすべてで num-complex・num-traits を直接 import する。以下の再公開に関する記述は移行前の記録である。

現行の `T` は複素係数全体ではなく、`Complex<T>` の実部・虚部の型である。以下では現行の成分型を `T`、表現に依存しない係数全体を `C`、norm などの実数値の型を `R` と書き分ける。

## 2. 依存箇所

### 2.1. 設定と公開 API

| 箇所 | 現在の依存 |
| --- | --- |
| [workspace 設定](../Cargo.toml)、[core 設定](../crates/qmb-core/Cargo.toml) | num-traits・num-complex を通常依存とする。両者の default features は無効。num-complex 自体も num-traits に依存する |
| [core の再公開](../crates/qmb-core/src/lib.rs) | `num_complex::Complex` と `num_traits` を再公開していた。移行は import の変更だけでは済まず、公開型・型制約にも及ぶ |
| [qmb の設定](../crates/qmb/Cargo.toml)、[再公開](../crates/qmb/src/lib.rs) | qmb は Complex を再公開しており、既定で core の `libm` を有効にする |
| core の `libm` feature | `num-traits/libm` を有効にする。平方根を使う公開メソッド自体も feature により有効化される。平方根の能力契約と、その具体実装・feature の配置を見直す必要がある |

### 2.2. core の通常コード

| 箇所 | ライブラリ依存と用途 |
| --- | --- |
| [expr/parts.rs](../crates/qmb-core/src/expr/parts.rs) | `ScaleExpr` が `Complex<T>` を保持する。保持自体に数値制約はない。Clone・Debug・PartialEq の導出に伴う制約は、それらの能力を使う場合に限る |
| [linear_combination.rs](../crates/qmb-core/src/linear_combination.rs) | 係数の保持、primitive の係数1の生成、既存係数への乗算。`primitive` に Zero・One を要求する |
| [expr/local.rs](../crates/qmb-core/src/expr/local.rs)、[expr/composite.rs](../crates/qmb-core/src/expr/composite.rs) | スカラー倍の節点・scalar identity を `Complex<T>` で保持する。式ハンドルの Clone は係数の Clone を要求しない |
| [ops.rs](../crates/qmb-core/src/ops.rs) | InnerProduct・Trace の返却値が `Complex<T>`。Product・Act・Adjoint は係数を含む LinearCombination を返す。trait 自体に数値制約はない |
| [coefficient.rs](../crates/qmb-core/src/coefficient.rs) | 成分からの構築、加算、乗算、共役、Zero による複素零の生成。算術は num-complex の演算子実装を使わず補助関数で行う |
| [expr/operators.rs](../crates/qmb-core/src/expr/operators.rs) | 負号・減算・commutator は `(-1, 0)` を生成するため Zero・One・Neg を要求する。スカラー倍の構文は係数を格納するだけ |
| [expr/evaluate.rs](../crates/qmb-core/src/expr/evaluate.rs) | 状態の作用評価で係数の複製・共役。項の走査と局所演算子の評価で乗算。内積・trace の集計で零生成・加算・乗算・共役を使う |
| [expr/composite_eval.rs](../crates/qmb-core/src/expr/composite_eval.rs) | 作用評価で複製・共役。演算子項の走査で乗算・共役。Full trace で零・単位元生成、局所 trace の積と項の和を計算する |
| [expr/inner_product.rs](../crates/qmb-core/src/expr/inner_product.rs) | 零・単位元生成、加算・乗算・共役。空の tensor product の内積と因子の積の初期値に数値の1を使う |
| 旧 `expr/normalization.rs`（削除済み。正規化は [qmb-normalize](../crates/qmb-normalize/src/lib.rs) へ移行） | `.re`・`.im`、Zero の `is_zero()`、成分による除算式、Float の `sqrt()` を使う |
| [basis.rs](../crates/qmb-core/src/basis.rs)、[local_structure.rs](../crates/qmb-core/src/local_structure.rs) | 基底展開・恒等元の返却型を介して係数表現に依存するが、trait 自体は数値能力を要求しない。次元・要素数の取得は係数に依存しない |

### 2.3. spin-half と利用側

Spin-half は数値ライブラリを直接 import せず、core の再公開を利用している。そのため、core の依存を外す際には同時に移行する必要がある。

| 箇所 | 現在の用途 |
| --- | --- |
| [spin の係数補助](../crates/qmb-spin-half/src/coefficient.rs) | `Complex::new`、零・単位元、加算・除算による `1/2`・`1/4`、実係数・純虚係数の構築 |
| [状態規則](../crates/qmb-spin-half/src/state.rs) | 内積の0・1と、作用規則の実係数・純虚係数・正負の分数 |
| [演算子規則](../crates/qmb-spin-half/src/operator.rs) | 積規則の分数・純虚係数、随伴の係数1、trace の0・2 |
| [局所構造](../crates/qmb-spin-half/src/local_structure.rs) | 恒等元を係数1の線型結合として返す |
| [基底](../crates/qmb-spin-half/src/basis.rs) | Sz 基底の係数1、Cartesian 基底の `±i`、Ladder 基底の `±1/2`・`±i/2` |
| [モデル utility](../crates/qmb-spin-half/src/utility.rs) | XYZ・Heisenberg は実数結合の複製・負号・虚部の零生成。Neel の構築には係数の数値能力を要求しない |
| テスト・doc test・[利用例](../crates/qmb/examples/heisenberg.rs) | Complex の構築・演算・比較・表示と、再公開された num-traits を利用する。利用例での表示能力を core の係数全体の条件にしない |

qmb-basic・qmb-geometry の通常コードには、これらの数値型・trait の直接利用はない。

## 3. ライブラリ由来の制約と実際の演算を区別する

確認したバージョンでは、以下の追加条件が付く。

| 現行の trait・実装 | 付随する条件 | 今回の方針による整理案 |
| --- | --- | --- |
| `num_traits::Zero` | 加算を supertrait に持ち、零生成とゼロ判定を一緒に要求する | そのまま採用し、生成だけの箇所でも付随条件を受け入れる |
| `num_traits::One` | 乗算を supertrait に持つ | そのまま採用し、生成だけの箇所でも付随条件を受け入れる |
| `num_traits::Float` | Num・Copy・NumCast・PartialOrd・Neg と、多数の数学関数 | norm 系で追加する数学関数は実数の平方根だけに限定する。他に必要な四則演算等は別々に要求する |
| `num_traits::Num` | 零・単位元、四則演算・剰余、PartialEq、基数指定の文字列変換 | 一般の評価に除算・剰余・比較・文字列変換を波及させない |
| num-complex の Add・Sub・Mul | `T: Clone + Num` | 所有値の成分加算なら Add だけでよい。乗算は現行の補助関数のように Clone・Add・Sub・Mul で実装できる |
| num-complex の `conj()` | `T: Clone + Num + Neg` | 現行の所有値を受け取る共役補助関数は `T: Neg` だけでよい。係数全体の契約では共役を独立した能力として扱える |

独自の Zero・One へ置き換える案では暗黙の Add・Mul の移行が必要だったが、現在はその案を採用せず num-traits を利用する。各操作の bound を読む際には、Zero・One がこれらの能力も含むことを考慮する。

## 4. 操作別に必要な能力

### 4.1. 現行の成分型 T に対する制約

以下は具体的な新 trait の定義ではなく、現行 API の制約を読み解くための表である。`Add` 等の出力は同じ `T` とする。表中の Zero・One は現行の num-traits のものであり、前節の付随条件を含む。

| 操作 | 現在の係数側の制約 | 実際の処理と分離できる点 |
| --- | --- | --- |
| 線型結合の保持・map・連結、局所式の primitive・zero・identity・sum・product・act・adjoint の構築、配置・束縛、式ハンドルの Clone | 数値制約なし | 入力を保持・共有するだけ。primitive の演算 trait 等は評価時に要求する |
| 式の `scaled(c)`・所有係数との乗算構文・`Operator::scalar(c)` | 数値制約なし | スカラー倍節点を作るだけ。借用した係数との構文では係数の Clone が必要 |
| `LinearCombination::primitive` | Zero・One | `(1, 0)` の生成だけ。加算・乗算・ゼロ判定は呼ばない |
| 式の負号・減算・commutator の遅延構築 | Zero・One・Neg | `(-1, 0)` の生成だけ。係数の加算・乗算や比較は呼ばない |
| `LinearCombination::scaled` | Clone・Add・Sub・Mul | 既存係数への乗算と、乗数の再利用。零・単位元生成や比較は不要 |
| 局所・複合状態の `evaluate` | Clone・Neg | 係数を複製し、随伴側で共役する。係数の積をスカラー倍節点の入れ子として残すため、core のこの経路に Add・Mul は不要 |
| 局所の項走査 `LinearTerms` | Clone・Add・Sub・Mul | 経路上の係数の積を集める。暗黙の1を Option で表すので One は不要 |
| 複合の演算子項走査 `OpTerms`、局所演算子の内部評価 | Clone・Add・Sub・Mul・Neg | 係数の乗算と随伴の共役。項走査自体には零・単位元生成は不要 |
| 局所の内積・期待値・trace・Hilbert–Schmidt 内積 | Clone・Zero・Sub・Mul・Neg | Zero 経由で Add も要求する。零を初期値に加算し、係数を乗算・共役する。ゼロ判定・One・除算は core のこの経路では不要 |
| 複合の内積・期待値・full trace・Hilbert–Schmidt 内積 | Clone・Zero・One・Sub・Neg | Zero 経由で Add、One 経由で Mul も要求する。局所の場合に加えて、空の tensor product 等の数値の1が必要 |
| `norm_squared` | 対応する内積と同じ | 自己内積の `.re` を返す。平方根・除算・ゼロ判定は不要 |
| `normalized_expectation`・`normalized_expectation_in` | 対応する内積の制約と Div | 自己内積全体のゼロ判定、分母の虚部のゼロ判定、複素数の商。平方根は不要 |
| `norm` | 対応する内積の制約と Float、`libm` feature | 自己内積の実部の平方根。norm 自体に除算・順序比較・ゼロ判定は不要だが、現在は Float・Zero 経由で能力が付く |
| `normalized` | 対応する内積の制約と Float、`libm` feature | 自己内積全体のゼロ判定、実部の平方根、`1 / norm` の計算、実係数の構築。局所版の One・Div は Float から得る |

ここでいう制約は core の共通処理の条件である。例えば `evaluate` 自体が Add・Mul を使わなくても、具体的な primitive の Act 実装が分数生成などの能力を要求すれば、その条件は別途必要になる。また、随伴節点がない入力でも、同じ公開メソッドが随伴を含む式を処理する以上、現在の `evaluate` は共役できる型を要求する。

### 4.2. 表現に依存しない係数 C に対する能力の整理案

採用した係数全体 `C` への移行に向けて、能力の対応を示す。

| 場面 | 必要な能力 |
| --- | --- |
| 保持・map・構文の構築 | 数値能力なし |
| 係数1の明示生成 | 単位元の生成 |
| 係数−1の明示生成 | 単位元の生成と負号、または同等の生成能力 |
| 既存係数の乗算・項走査 | 係数の乗算と、再利用する場合の複製または借用演算 |
| 状態の作用評価 | 係数の複製と共役。係数の乗算を行わない現在の経路を維持できる |
| 局所の内積・trace の集計 | 零生成・加算・乗算・共役・係数の再利用 |
| 複合の内積・trace の集計 | 前行に加えて単位元生成 |
| `norm_squared` | 内積の能力と、係数から実数値 `R` を取り出す能力 |
| 正規化期待値 | 内積・期待値の能力、係数の厳密なゼロ判定、除算 |
| `norm` | `norm_squared` の能力と、`R` の平方根 |
| 状態の正規化 | 内積の能力、係数の厳密なゼロ判定、実数値の取得、`R` の平方根・単位元・除算、実数値を係数へ埋め込む能力 |

加減乗除・負号は標準の `core::ops` で表せる。零生成・単位元生成・ゼロ判定は num-traits の Zero・One を採用し、付随する能力を受け入れる。共役・実数値取得・実係数の埋め込み・純虚係数の生成は個別の能力とし、具体的な trait 名やシグネチャは実装時に定める。平方根は qmb-normalize の実数型 R に Float を要求して提供する。

実数値の取得・埋め込みには `C` と `R` の型の対応が必要だが、その対応をすべての式や四則演算に要求する必要はない。共役能力にも、実部・虚部の公開や実数型への変換を必須にする必要はない。

この表は計算に必要な能力を示すものであり、平方根の trait・具体実装・正規化 API のすべてを core に置くという結論ではない。tensory-regulated を参考にした責任の分け方を第11節に記録する。

## 5. 特に切り分ける必要がある操作

### 5.1. 構文上の zero・identity と数値の0・1

`LinearCombination::zero` と式の zero は空和、局所演算子の identity は専用節点、空の状態 tensor product は Unit 節点で表す。これらの構築時に係数の零・単位元は生成しない。

一方、scalar を返す空和・空積の評価では実際の0・1が必要になる。`LinearCombination::primitive` も係数1を明示するため、局所式の `primitive` とは条件が異なる。メソッド名の類似だけで同じ bound を付けない。

### 5.2. 零生成・ゼロ判定・一般の比較

内積・trace の集計では零を初期値にするが、値が零かどうかは調べない。厳密なゼロ判定を使うのは現在の正規化・正規化期待値と、商の実分母経路の選択である。能力としては区別できるが、API では num-traits の Zero にまとめたまま利用する。

ゼロ判定には Zero::is_zero を使い、独立した IsZero trait は追加しない。一般の PartialEq は一律に要求せず、近似ゼロ判定や許容誤差も持ち込まない。Zero を提供しない型でも、その能力を必要としない保持・構築などには利用できる。Default が数学的な零を表すとも仮定しない。

### 5.3. 共役と負号

現行の共役は `(a, b) -> (a, -b)` なので成分型 T の Neg を使う。係数全体 C の共役と、C 自体の負号は別の能力である。一般の表現に対して、共役の条件を自動的に Neg と同一視しない。

### 5.4. Norm と正規化期待値

`norm_squared` は状態の自己内積の実部であり、複素係数単体の絶対値の二乗を計算する操作ではない。`norm` はその実数値の平方根だけを追加する。正規化期待値には平方根を追加しない。

現在の商は分母の虚部が零なら各成分を実部で割り、それ以外は複素数の除算式を使う。前者は大きい／小さい実分母を二乗することを避ける。係数 C の Div に委譲する場合、この計算法や数値特性を具体型の責任としてどう維持するかも検討する。ゼロ除算の拒否は、商の前の自己内積の検査で行う。

Float を分解する際も、負値の補正、近似ゼロ判定、NaN・infinity の独自処理は追加しない。平方根や特殊値の挙動は具体的な数値実装に委ねる。

### 5.5. Clone と数値能力

Clone は演算能力とは別に扱う。式ハンドルの共有、所有係数の格納、map での係数の移動に係数の Clone は不要である。現在の評価器・項走査は共有節点から係数を取り出して再利用するため Clone を要求する。

借用演算で一部の Clone を減らせる可能性はあるが、借用から結果を作る際の実装まで無複製になる保証ではない。Copy・Send・Sync・'static を追加の共通条件にはしない。

## 6. 具体的な spin-half 規則に必要な能力

現在の分数生成方法と複素数の成分表現を維持する場合の、実際に使う能力は以下である。Zero・One に付随している演算能力を除いて記述する。

| 規則・構築 | 実際に必要な成分型 T の能力 |
| --- | --- |
| 状態 primitive の内積 | 零・単位元生成 |
| 演算子 primitive の随伴、局所構造の恒等元、Sz 基底への展開 | 零・単位元生成 |
| 演算子 primitive の trace | 零・単位元生成、`1 + 1` のための Add |
| Cartesian 基底への展開 | 零・単位元生成、Neg |
| Ladder 基底への展開、primitive の Product・Act | 零・単位元生成、分数生成の Add・Div、Neg |
| XYZ・Heisenberg の遅延構築 | 結合の Clone・Neg、虚部の零生成 |
| Neel の遅延構築 | 数値能力なし |

これらの具体的実装には成分型 T の Mul、Copy、Num、Float、文字列変換・数値 cast は不要である。`1/2`・`1/4` の生成方法を変えれば Div の代わりに分数定数を提供する能力も候補になるが、これは具体的な局所代数側の判断とする。core の全係数へ分数・虚数単位の生成を要求しない。

型が Div を実装するだけでは、spin の分数を正しく表せる保証にはならない。整数除算では契約を満たさない。必要な値と演算規則の数学的な整合性は、具体的実装の契約と検証で扱う。

基底の `index`・`elements` 自体は数値計算しないが、現在は展開と同じ `Basis<P, T>` 実装なので、その impl の数値制約も適用される。これは基底メタデータ自体に数値能力が必要という意味ではない。基底 API を分割するかどうかは、今回の能力調査だけでは決めない。

## 7. 能力の存在と数学的な契約

Rust の bound が保証するのは、呼び出す演算と入出力の型である。係数の単位元則、スカラーの可換性・分配則、共役の involution と演算との整合性、実数値の埋め込み、局所規則との精度の一致は、具体型・代数の数学的な契約として別に記述する。浮動小数点に厳密な結合則やビット単位の一致を要求しない。

除算を使わない経路に体としての全能力を要求せず、平方根を使わない経路に実数の平方根への閉性も要求しない。有理数などでも、表せる演算だけを利用できるようにする。

## 8. テスト側の依存と今後の確認

テストには、参照行列の演算・リテラルの変換・誤差比較のための Num・FromPrimitive・ToPrimitive・Copy・Debug や、複素数の演算子・`conj()` がある。これらは通常コードで要求する最小能力とは区別し、移行後の core の公開契約へ戻さない。

既存の以下のテストを、移行後も能力を増やしていないことの基準にする。

- [線型結合](../crates/qmb-core/tests/linear_combination.rs)：数値能力のない係数の map、Clone・Num のない係数1の生成、零・単位元・比較を持たない型でのスケーリング。
- [局所式](../crates/qmb-core/src/expr/tests.rs)、[複合式](../crates/qmb-core/src/expr/composite_tests.rs)：不透明な係数での構築、加算・乗算のない型での作用評価、Num・除算・剰余・Copy を持たない型での集計。
- [正規化期待値の移行先](../crates/qmb-normalize/tests/expectation.rs)：Float・Num・一般の PartialEq がない型での norm squared・正規化期待値、複素分母の商。
- [spin の最小 scalar](../crates/qmb-spin-half/tests/minimal_scalar.rs)：Clone・Copy・Num を持たない分数型での局所規則。
- [独立した三次元代数](../crates/qmb-core/tests/local_algebra.rs)、[spin の行列照合](../crates/qmb-spin-half/tests/contracts.rs)、[モデル](../crates/qmb/tests/chain.rs)：単精度・倍精度、非正規直交基底、局所・複合・縦断計算の一致。

新しい契約を実装した段階では、数値能力を持たない型での保持・構築、num-traits の Zero・One の利用、一般の PartialEq を持たない型での Zero::is_zero、実数型だけに Float を要求する qmb-normalize の利用を検証する。独自の零生成のみ／単位元生成のみの trait や平方根 provider の検証は不要とする。

今回行ったのはソースと解決バージョンの調査であり、上記テストの再実行や新しい契約のコンパイル検証はまだ行っていない。

## 9. 現時点の結論と次の設計判断

基本方針に沿って、まず次の能力を分けて考える。

1. 数値能力のない保持・共有・構文構築。
2. 零生成、単位元生成、厳密なゼロ判定。
3. 標準 trait による個別の算術、独立した共役。
4. 実数値との接続と、実数の平方根。
5. 具体的な代数に固有の分数・純虚係数の構築。

Zero・One は num-traits を再利用し、独自の細分化はしない。それ以外については、現在の core が独自の係数補助関数で保っている狭い制約を、移行時に Num・Float へ一律に広げない。

係数全体 C を保持し、個別の能力を使う方針を採用した。標準の Complex への追加 trait は core 内に置く。関連する実数型、実係数の埋め込み、純虚係数の生成は必要な操作に限定し、標準 Complex の左スカラー倍は個別実装、任意の C は scaled(c) を基本経路とする。利用者から挙げられた tensory の演算設計との比較を次節に記録する。

## 10. tensory の演算設計との比較

利用者が参考として挙げた [tensory の llm-work-in-use ブランチ](https://github.com/todo-group/tensory/tree/llm-work-in-use)を取得し、2026-10-25 時点のコミット `ac4bbcf175d1d89e0dab91a7f523454f151f1d59` のソースを確認した。以下の参照はそのコミットに固定する。GitHub の Web 取得は失敗したため、Git 経由で取得したソースを読んだ。ビルド・実行による検証は行っていない。

### 10.1. 確認した構成

| 対象 | ソースで確認した内容 | qmb への示唆 |
| --- | --- | --- |
| [加算](https://github.com/todo-group/tensory/blob/ac4bbcf175d1d89e0dab91a7f523454f151f1d59/tensory-core/src/arith/add.rs)・[負号](https://github.com/todo-group/tensory/blob/ac4bbcf175d1d89e0dab91a7f523454f151f1d59/tensory-core/src/arith/neg.rs) | AddOp・NegOp と、演算を表す遅延表現を定義する。Tensor の演算子実装は数値計算ではなく、その表現を構築する | 構文の構築時と係数を計算する時点の bound を分ける、という今回の基本方針に対応する |
| [スカラー乗算](https://github.com/todo-group/tensory/blob/ac4bbcf175d1d89e0dab91a7f523454f151f1d59/tensory-core/src/arith/scalar_mul.rs)・[除算](https://github.com/todo-group/tensory/blob/ac4bbcf175d1d89e0dab91a7f523454f151f1d59/tensory-core/src/arith/scalar_div.rs) | 個別の拡張 trait があり、`left_mul`・`right_mul` 等はスカラー E の数値能力を要求せず遅延表現を構築する | 演算ごとの能力分離と名前付きメソッドを参考にできる。qmb の `scaled` の構築にも係数の Mul は不要 |
| [Context・RuntimeFor](https://github.com/todo-group/tensory/blob/ac4bbcf175d1d89e0dab91a7f523454f151f1d59/tensory-core/src/concept/task.rs) | Context が task ごとの実行方法と結果型を提供し、RuntimeFor がその task の context を関連付ける | 必要な演算を実行する箇所に能力を要求する構成の具体例。qmb では既存の evaluator と係数の個別操作の境界に対応させて検討する |
| [ndarray 側の加算実装](https://github.com/todo-group/tensory/blob/ac4bbcf175d1d89e0dab91a7f523454f151f1d59/tensory-ndarray/src/arith.rs) | 実際の要素加算を行う Context 実装で `&'l E: Add<&'r E, Output = E>` を要求する。この加算実装は E: Clone を要求しない | 所有値・借用値のどちらを演算の入力にするかも、Clone と独立に検討する根拠になる |
| [arith::Scalar](https://github.com/todo-group/tensory/blob/ac4bbcf175d1d89e0dab91a7f523454f151f1d59/tensory-core/src/arith/mod.rs) | 空の marker trait。`Mul<E>`・`Div<E>` の E をスカラーとして扱う構文の条件に用いる。名前付きメソッドと `(E,)` を介した構文には、この marker を要求しない経路がある | 「スカラーとして扱う型」と「何の演算ができる型か」は別の役割として考えられる。marker 自体は零・単位元・四則演算・共役を保証しない |

このブランチの `arith` に零生成・ゼロ判定・平方根の契約はなく、共役 module の宣言はコメントアウトされている。今回不足している数値能力をそのまま満たす trait 一式として採用できるものではなく、能力を要求する場所と演算の組み立て方を参考にする。

### 10.2. 今回の設計候補へ反映する点

1. **構築と実行の分離を保つ**：数値能力を要求しない式の型・constructor と、能力を要求する評価メソッドを分ける。具体的な数値型を制限する marker だけで演算能力を代用しない。
2. **演算ごとの契約を使う**：標準の Add・Mul 等、num-traits の Zero・One、不足する共役等を操作に合わせて要求する。tensory の実行 context の役割は、qmb では既存の evaluator と各係数操作に対応付ける。今回 Context／Runtime は新設しない。
3. **借用演算も候補にする**：現在の Clone による係数の取り出しと、`&C` を入力とする演算を比較する。借用の Add があるだけで、評価器全体の Clone が不要になるとは判断しない。
4. **結果型の条件を操作側で決める**：一般の演算 trait は結果型を表せるようにし、同じ係数型・精度を維持する共通評価では、その操作で必要な `Output = C` を要求する。tensory の結果型の柔軟性を参考にしても、qmb で暗黙の精度変換を許可する結論にはしない。
5. **左・右の構文と数学的意味を区別する**：tensory は左・右のスカラー操作を別の演算表現にする。qmb の係数の左右のスカラー倍は現行契約では同じ意味であり、その意味を保ったうえで、係数 C をジェネリックにしたときの trait 実装と構文を検討する。

なお、tensory の NegOp は負号の構築自体に数値の−1を生成しない。qmb も独立した負号節点に変えれば構築時の One 等を外せる可能性がある。ただし、これは係数 trait の置換に加えて節点・評価器・破棄処理を変更する案なので、現行の−1のスカラー倍を維持する案との比較事項として扱う。

## 11. tensory-regulated の正規化と平方根の責任

同じコミット `ac4bbcf175d1d89e0dab91a7f523454f151f1d59` の tensory-regulated と、その利用先の tensory-linalg・tensory-ndarray を追加調査した。以下はソースによる確認であり、ビルド・実行による検証は行っていない。

### 11.1. 抽象契約・正規化方針・数値実装の分離

| 層 | 確認した責任 |
| --- | --- |
| [tensory-regulated の抽象契約](https://github.com/todo-group/tensory/blob/ac4bbcf175d1d89e0dab91a7f523454f151f1d59/tensory-regulated/src/lib.rs#L35-L98) | RegulatedRepr が本体・係数表現を分けて持つ。Regulator は本体を変換して係数を取り出し、Inflator は係数を本体へ戻す。Regulation が対応する scalar 型を指定する。これらの trait は Float・Sqrt を要求しない |
| 同 crate の CoefficientRepr・UnpackableCoefficientRepr | 取り出した係数の保持・結合・分解と、必要な場合の scalar 値への復元を担当する。平方根の契約は持たない。保持表現 C と、実際に本体へ掛ける C::Scalar が区別されている |
| [tensory-linalg の norm](https://github.com/todo-group/tensory/blob/ac4bbcf175d1d89e0dab91a7f523454f151f1d59/tensory-linalg/src/norm.rs) | NormRepr と構築用の拡張 trait を提供し、実際の norm 計算は Context に委譲する。構築時には係数の平方根能力を要求しない |
| [tensory-ndarray の L2Regulator](https://github.com/todo-group/tensory/blob/ac4bbcf175d1d89e0dab91a7f523454f151f1d59/tensory-ndarray/src/regulated.rs#L26-L64) | norm を取得し、本体をその値で割り、正規化した本体と実数係数を返す。実数から要素型への変換もこの具体実装で扱う |
| [tensory-ndarray の norm 実装](https://github.com/todo-group/tensory/blob/ac4bbcf175d1d89e0dab91a7f523454f151f1d59/tensory-ndarray/src/linalg.rs#L224-L230) | NormRepr の Context が `data.norm_l2()` を呼び、`E::Real` を返す。L2 norm の実際の数値計算はこの backend から ndarray-linalg に委譲される |
| [tensory-ndarray の LnCoeff](https://github.com/todo-group/tensory/blob/ac4bbcf175d1d89e0dab91a7f523454f151f1d59/tensory-ndarray/src/regulated.rs#L67-L109) | 係数を対数で保持し、積を対数の和、商を差、復元を exp で扱う。具体的な ln・exp と数値型の制約はここに置かれる |

したがって、正規化の抽象契約は tensory-regulated が担い、L2 norm を選ぶ方針と数値実装は tensory-ndarray 側にある。tensory-regulated が平方根を直接計算する構成ではない。

通常依存も [tensory-regulated の Cargo.toml](https://github.com/todo-group/tensory/blob/ac4bbcf175d1d89e0dab91a7f523454f151f1d59/tensory-regulated/Cargo.toml)では tensory-core・tensory-linalg のみで、num-traits・num-complex は dev-dependencies にある。

### 11.2. qmb の検討へ反映する点

前節までの「norm には実数の平方根が必要」という能力の分析と、その能力・処理をどこへ置くかは別の判断である。qmb では次の三つを分ける。

1. **代数的な共通処理**：式の保持・作用・内積・norm squared・スカラー除算。これらには平方根を要求しない。
2. **Norm・正規化の導出**：自己内積から norm と逆数係数を得て状態をスカラー倍する処理と、正規化期待値・零ノルムの検査。qmb-normalize の拡張 trait に移す。
3. **具体的な数値計算**：qmb-normalize で実数型 R に Float を要求し、std 対応の平方根を利用する。core の一般の係数には要求しない。

tensory の構成を参考に、調査時には正規化方針を渡す案も比較した。採用した構成では provider は導入しない。qmb の通常評価は式と局所規則から行列化せずに行うため、norm・正規化を提供する経路もその契約を維持する。数値コンテナへの export を必須の前提にはしない。

また、tensory-regulated の CoefficientRepr は mul・div・factorize 等をまとめた「規格化した本体に付随する係数表現」の契約である。qmb のすべての係数・保持操作へそのまま要求する案とは分ける。取り出したスケールを別表現で保持する設計は参考になるが、対数係数の導入や評価時の自動規格化は今回の採用事項ではない。

qmb の状態の正規化は zero-norm をエラーとし、結果には取り出したスケールを別途付随させない。一方、確認した L2Regulator は本体と係数の組を返し、norm が零の場合の明示的な分岐を持たない。この API とエラー契約の違いを踏まえて対応付ける。

調査後に qmb-normalize への分離と拡張 trait での公開を採用した。既存 libm feature と norm API の移行は、その方針に沿って実装時に整理する。

## 12. 設計案への接続

上記の調査を踏まえ、[係数型の採用方針](refactoring-coefficient-proposal.md)に比較と採用事項をまとめた。core の標準 Complex への直接実装、係数全体 C の保持と個別の能力契約、qmb-normalize への norm・正規化の分離を採用している。Rust の局所的な実装制約は独立した試作で確認したが、製品コードはその後に移行済みであり、現在の API は採用方針の第10節を参照する。

## 13. ndarray の要素型と数値ライブラリの接続

ndarray 0.17.2 の公開ドキュメントとソースを確認した。tensory の指定ブランチも、このバージョンを依存に指定している。ソース・依存情報の確認であり、ndarray のビルド・実行は行っていない。

| 論点 | ndarray の対応 | qmb との関係 |
| --- | --- | --- |
| 保持と演算 | [ArrayBase](https://docs.rs/ndarray/0.17.2/ndarray/struct.ArrayBase.html) は数値以外の要素も保持する一般のコンテナ。演算時に要素の能力を要求する | 数値能力を保持型へ一律に要求しない方針の参考になる |
| 加減乗除 | [演算実装](https://docs.rs/ndarray/latest/src/ndarray/impl_ops.rs.html) は標準 trait を使い、例えば配列間の加算の一経路では A: Clone + Add<B, Output = A> と B: Clone を要求する | Complex の提供側が実装済みの演算を利用するため、ndarray 独自の加算 trait を Complex に追加する必要はない |
| 数値用のまとまった能力 | [LinalgScalar](https://docs.rs/ndarray/0.17.2/ndarray/trait.LinalgScalar.html) は Copy・'static・Zero・One・四則演算をまとめ、条件を満たすすべての型へ blanket impl を提供する | 既存の能力をまとめるだけなら外部型へも適用できる。ただし qmb のすべての操作にこの広い条件を要求する必要はない |
| 右スカラー演算の識別 | [ScalarOperand](https://docs.rs/ndarray/0.17.2/ndarray/trait.ScalarOperand.html) を ndarray 自身が定義し、Complex<f32>・Complex<f64> に直接実装する。配列全般の要素型を制限する trait ではない | trait の定義元 crate が実装するため orphan rule に抵触しない。同 trait 自体は Clone・'static を要求する |
| 左スカラー演算 | ScalarOperand の文書と演算ソースでは、左側は具体的な scalar 型ごとに個別 impl を用意する | 任意の C に対する裸の左スカラーの blanket impl ができるわけではない |
| 数値ライブラリ依存 | [依存情報](https://docs.rs/crate/ndarray/0.17.2)で num-traits・num-complex を通常依存として持つ | qmb の旧依存解消目標とは異なるが、現在採用した core 内での標準 Complex への接続方針の参考になる |

ndarray が行っているのは、要素型の既存の標準演算の利用、既存能力の blanket impl による集約、trait 定義元での外部型への直接実装である。数値ライブラリへの依存を持ちながら一般の要素型も扱う設計であり、独自の新 trait を依存のない core から外部型へ接続する仕組みとは区別する。

qmb への示唆は、標準演算や Zero・One をそのまま使える部分には newtype・provider が不要であり、不足する共役等も、trait 定義元が数値ライブラリに依存すれば直接実装できるという点にある。qmb も標準の Complex についてこの接続方針を採用した。

## 14. num-complex を含める方針と num-traits の使い分け

**合意事項**：core の num-complex・num-traits 依存を維持し、零・単位元には num_traits::Zero・One をそのまま使う。付随する Add・Mul・ゼロ判定は受け入れ、独自の Zero・One・IsZero は作らない。不足する追加 trait は core 内で Complex<f32>・Complex<f64> に直接実装する。標準の複素数接続には newtype・provider を必須にしない。

両ライブラリは `default-features = false` の直接依存として維持する。One・Zero を必要なメソッドにだけ課し、保持や能力を使わない遅延構築に広げない。Num・Float 全体を一般の係数の必須条件にする方針ではない。

確認した num-complex 0.4.6 自体が num-traits 0.2 系を通常依存に持つため、同じ解決バージョンの num-traits を直接利用しても新しい数値ライブラリが増えるわけではない。Cargo.lock では 0.2.19 が解決されている。両 crate の default features を無効にし、libm は必要な場合だけ有効にする現在の構成も維持できる。

| 能力・用途 | 利用方法の案 |
| --- | --- |
| Add・Sub・Mul・Div・Neg | 標準の core::ops を使う。具体的な Complex<f32>・Complex<f64> は既存実装を利用できる |
| 零生成・単位元生成 | num-traits の Zero・One をそのまま使う。生成だけでも付随する能力を受け入れ、再定義しない |
| 加算・乗算 | 標準演算を使う。Zero が Add、One が Mul を含むため、重複する条件を考慮する |
| 厳密なゼロ判定 | Zero::is_zero を使う。一般の PartialEq・Num は追加しない |
| 共役・実数値との接続 | 必要な qmb の個別 trait を定義し、標準 Complex に直接実装する |
| 平方根 | qmb-normalize の平方根を使う実数型 R に Float を要求する。core の係数 C 全体には要求しない |
| リテラル変換・参照計算 | FromPrimitive・ToPrimitive 等は必要な具体接続・テストで利用する。全係数の共通条件にしない |

例えば、具体的な Complex<f64> の共役実装で num-complex の `conj()` を使っても、共通評価の C に Num を要求することにはならない。一方、Complex<T> の T までジェネリックな実装で `conj()` を使う場合は T: Num 等が必要になるため、成分の最小限の処理で実装する案と比較する。

方針採用後、能力契約と各メソッドの bound を実装し、qmb-normalize へ分離した。core の両数値ライブラリへの依存は維持し、core/libm は削除した。
