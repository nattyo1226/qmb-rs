# Ordinary / graded tensor product の実装修正案

作成日: 2026-10-29
更新日: 2026-10-29
状態: 検討案。以下の型名・メソッド名は提案であり、実装済み API ではない。

## 1. 目的と範囲

現在の ordinary 専用実装を、ordinary と Z₂-graded の二種類の tensor product に対応させるための変更案を記録する。依頼中の `ordinal` は、現行の要件・設計で用いる `ordinary`（通常のテンソル積）を指すものとして扱う。

推奨する構成は、空間が積の種類と因子順序を保持し、primitive が局所 parity と parity 演算子を提供し、core の共通評価機構が符号を処理する方式である。`HilbertSpace<I>`、`State<I,S,O,C>`、`Operator<I,O,C>` の型引数は維持する。利用者は既存 API を残さない breaking change を許容しており、互換性のための別入口は設けず、共通の評価 API が空間の種類を参照する方針とする。

対象は空間、状態の tensor 構築、配置・束縛、作用、随伴、内積、full trace、基底展開・座標取得、および既存の密・疎 export への接続とする。Fermion は規則を検証する具体例として用い、製品用 fermion crate の全 API・モデル設計は別計画に委ねる。Spin-one、異種 primitive 型の混在、ordinary 積と graded 積の実質的な混在、permutation、regrouping、partial trace、State / Operator の型統合は今回追加しない。これは [v2 の確定範囲](README.md)全体の実装完了を意味しない。

## 2. 参照した現行契約と実装

リポジトリの実際の設計ディレクトリは `docs/design/` である。既存文書には将来の数学的契約と実装済み機能の両方があるため、両者を区別する。

| 参照文書 | 本案で引き継ぐ契約 |
| --- | --- |
| [要件](../../requirements.md) 2.4、3.2、5.1〜5.5、6.2 節 | homogeneous な primitive、空間選択、graded の符号、通常の trace、基底順序 |
| [アーキテクチャ](../../architecture.md) 2、10.2 節 | core / 具体的代数 / geometry / export の責任分離、grading の保持・API の具体化 |
| [Hilbert 空間](../../design/hilbert-space.md) 4、6〜8 節 | 空・単因子の同一視、同種積の結合則、混在拒否、部分系、深い木の処理 |
| [複合式](../../design/composite-expression.md) 4〜9、11 節 | 埋め込みの準同型性、作用・随伴・内積・trace、遅延評価、fallible 集約 |
| [局所構造](../../design/local-structure.md)・[primitive / 式](../../design/primitive-expression.md) | LocalStructure を再導入しない、局所能力の分離、構築時に評価能力を要求しない |
| [基底](../../design/basis.md) 5、8 節・[export](../../design/export.md) 4 節 | 展開結果は共通式木、graded 演算子基底、右端が最速、列ごとの作用による行列化 |
| [幾何](../../design/geometry.md)・[今後の課題](../../future-work.md) | Chain は幾何のみ、grading と State / Operator 統合は独立、partial trace 等の延期 |

以下は調査時点のソースで確認した事実である。

| 現在の実装 | 二種類に対応する際の不足 |
| --- | --- |
| [space.rs](../../../crates/qmb-core/src/space.rs) の `Node::Tensor(SumExpr<...>)` | 積の種類を保持しない。互換性は葉列の一致、部分系は部分列だけを検査する |
| [ops.rs](../../../crates/qmb-core/src/ops.rs) | parity 取得・parity 演算子の能力がない |
| [expr/composite.rs](../../../crates/qmb-core/src/expr/composite.rs) の `State::tensor_collected` | ordinary 空間を構築する。`StateNode::Tensor` に種類がなく、`StateNode::Act` は選択済み空間を保持しない |
| [evaluate/composite.rs](../../../crates/qmb-core/src/evaluate/composite.rs) の `Task::Local` | 配置先だけを更新し、別因子には作用しない。`evaluate` は `Product` や係数の `One`・`Neg` を要求しない |
| [expr/traversal.rs](../../../crates/qmb-core/src/expr/traversal.rs) の `OpTerms` | 入力を左から読み、同じ Index の局所演算子をまとめる。異なる Index をまたぐ交換の符号は計算しない |
| [evaluate/trace.rs](../../../crates/qmb-core/src/evaluate/trace.rs) | 上記 `OpTerms` を使い、構文で参照されない因子を identity として局所 trace を掛ける |
| [evaluate/inner_product.rs](../../../crates/qmb-core/src/evaluate/inner_product.rs) | 評価済み状態を通常の正定値な積で縮約する。この縮約規則自体に graded 符号は不要 |
| [basis/expansion/composite.rs](../../../crates/qmb-core/src/basis/expansion/composite.rs) | 演算子展開にも ordinary な `OpTerms` を使用。局所因子の展開結果を因子順の配置積として再構築する |
| [basis/enumeration.rs](../../../crates/qmb-core/src/basis/enumeration.rs) | 空間の葉順で状態 tensor・演算子の配置積を生成する。graded 空間を保持する変更が必要 |
| [qmb-export の matrix.rs](../../../crates/qmb-export/src/matrix.rs) | 状態基底の各要素への作用と座標から列を作る。ここで独自に符号を足す必要はない |
| [system.rs](../../../crates/qmb-core/src/system.rs) | `state_tensor`・`try_state_tensor`・`state_tensor_at` は ordinary な構築入口のみ |

特に `OpTerms` は trace と演算子基底に共有されている。ここだけへ一律に parity string を追加すると、graded 基底の座標と通常の行列表示を混同するため、用途に応じた二つの内部表現を区別する必要がある。

## 3. 情報の配置と API 方針

### 3.1. 積の種類は空間が所有する

core に `TensorKind { Ordinary, Graded }` を追加し、非公開節点を概念的に `Tensor { kind, children }` とする。`Factor` に parity や局所代数の ID を追加しない。積の種類は合成則の情報であり、状態の even / odd とは異なる。

| 提案 API | 意味 |
| --- | --- |
| `HilbertSpace::tensor_with_kind(kind, children)` | 種類を明示して構築。保持する子の grouping は変更しない |
| `HilbertSpace::tensor(children)` | `tensor_with_kind(TensorKind::Ordinary, children)` の短縮形 |
| `HilbertSpace::tensor_g(children)` | `tensor_with_kind(TensorKind::Graded, children)` の短縮形 |
| `tensor_product_with_kind(kind, &rhs)`・`tensor_product(&rhs)`・`tensor_product_g(&rhs)` | 同じ三種類の二項版。因子順序は左、右の順 |
| `kind() -> Option<TensorKind>` | 保持する Tensor 節点の種類。Factor は `None` |
| `effective_kind() -> Option<TensorKind>` | 同一視後の種類。基本因子が0または1個なら `None` |

三種類の入口は共通の構築・検証処理へ委譲し、重複 Index・種類の混在・grouping の扱いを一致させる。`tensor` と `tensor_g` は呼び出しの便宜のための入口であり、primitive や子空間から種類を推測しない。種類の取得には、構築名と区別した `kind()`・`effective_kind()` を用いる。

`unit()` は ordinary の空 Tensor として保持するが、共通の単位空間として両種類と互換である。`tensor_g([])` で明示的な空の graded Tensor も構築でき、数学的意味は同じとする。

### 3.2. 空・単因子を正規化した後で検証する

今回の対応範囲では、一つの空間内で ordinary / graded の実質的な混在を許さない。これは利用者との合意事項とし、空・単因子の同一視後に残る積の種類について検証する。

積の種類の混在と異種局所代数の混在は別の条件である。同じ局所 primitive 型でも、ordinary / graded によって因子間の交換関係が変わる。混在を扱うには木の各境界に応じた合成則・互換性判定が必要になるため、今回は両者とも対応範囲から外す。

木は保持したまま、各空間に「葉数の区分（0 / 1 / 複数）」と「有効な積の種類」を構造 metadata として持たせる案を推奨する。局所次元・全基底・評価結果のキャッシュは追加しない。

新しい Tensor の子について、空の子を意味上除いた個数を `m` とする。

1. `m = 0` なら共通の単位空間。有効な種類はない。
2. `m = 1` なら、その唯一の非空の子の有効な種類を継承する。外側の種類は互換性に影響しない。
3. `m >= 2` なら、親の種類を有効な種類とする。有効な種類を持つ子はすべて親と同じ種類であることを要求する。

例えば `T_o(E,T_g(a,b))` は単位元除去後に `T_g(a,b)` と同じなので許可する。`T_o(a,T_g(b,c))` は二種類が残るので `SpaceError::MixedTensorKinds` 相当で拒否する。単純に「親子の kind が違えば拒否」してはいけない。

評価経路の選択には、外側の保持された kind だけでなく、正規化後の有効な種類を用いる。有効な種類がない0・1葉では両規則が一致し、交差符号は生じない。

互換性は従来の順序付き葉列に加え、有効な種類が一致することを検査する。部分系は、0・1葉では従来の規則、複数葉では相対順序と有効な種類の一致を要求する。ordinary と graded の二因子空間は、同じ Index 列でも非互換となる。

Index 重複の検証も従来どおり行う。深い単因子 wrapper は子の metadata を定数時間で継承し、構築のたびに全子孫を走査しない。`I: Eq` の重複検査の計算量は現状を引き継ぎ、Hash・Ord・Clone を新たに要求しない。走査・破棄は反復処理を維持する。

### 3.3. Primitive が提供する grading 能力

次の独立した能力を core に追加する。シグネチャは検討用である。

```rust,ignore
enum Parity { Even, Odd }

trait HomogeneousParity {
    fn parity(&self) -> Parity;
}

trait ParityOperator: Sized {
    fn parity_operator() -> Self;
}
```

graded な局所代数の S・O は homogeneous な primitive とし、両方が `HomogeneousParity` を提供する。混合 parity は `LocalState` / `LocalOperator` の和で表す。`parity() -> Option<Parity>` の `None` を even と扱う方式や、任意の型を自動的に even とする blanket impl は採用しない。

O は局所 parity 演算子 `P` 自身を primitive の値として `ParityOperator` で提供する。Identity と同様に係数型へ依存させない。P は選択基底の要素である必要はなく、基底展開では例えば `I - 2n` と表せる。必要なら P を表す線型従属な primitive を追加する。

契約は `P† = P`、`P² = I`、`Pψ = (-1)^p(ψ) ψ`、`POP = (-1)^p(O) O` とする。Identity と P は even、随伴は parity を保存し、局所 Product / Act の非零な出力項の parity は入力の和 modulo 2 と一致する。S の even・odd 部分は内積で直交する。局所基底の primitive もこの契約に従う。

これらの代数法則は具体的代数の契約・テストで保証する。一般のユーザ実装の法則を実行時に証明しない。式全体の parity を必要とする内部処理は even / odd 成分へ分け、混合した和に一つの parity を割り当てない。Zero に便宜的な parity を与えても結果を変えない。

### 3.4. 空間の種類に従う統一評価 API

既存 API を残さない breaking change を許容するという利用者の方針を受け、`evaluate`・`trace`・`expand_in`・`coordinates_in` 等を両種類に対応する共通入口へ変更する。`evaluate_graded`・`trace_graded`・`GradedExport` のような互換性維持用の別入口・adapter は追加しない。

Rust の trait bound は実行時の kind で切り替えられないため、複合評価にはその操作が graded 側で必要とする能力を ordinary 空間でも要求する。ordinary 専用の独自代数も、共通評価を利用する場合は明示的な grading 実装へ移行する。非対応を暗黙の even と扱う blanket impl は作らない。

実際の合成則は選択済み空間の有効な種類で決める。ordinary 空間では parity による交差符号・P の挿入を行わず、graded 空間ではそれらを処理する。primitive に grading 能力があることから空間の種類を推測・変更しない。

局所式の構築・`at`・`bind`・和・積・`adjoint`・遅延した `act` に評価能力を一律に要求しない。局所だけで完結する評価の bound も独立に保つ。互換性を前提にしないことは、全操作に全能力を要求する理由にはならない。

符号係数を出力する処理には `Neg<Output=C>` 等を要求する。一方、P の primitive への局所 Act と既存の Scale 節点を使う作用経路は、交差符号のためだけに係数の One・Neg を追加する必要がない。状態への作用に O の Product を要求しない契約も維持する。

変更する共通入口の範囲は次のとおりとする。明示空間版（`_in` / `_in_space`）も同じ評価方針へ接続し、旧束縛より明示空間を優先する契約を維持する。

| 対象 | graded 対応する操作 | 追加で必要となる主な能力 |
| --- | --- | --- |
| State | evaluate、inner_product、norm_squared、expectation と明示空間版 | O の parity と P。内積等の既存 bound はそのまま必要 |
| Operator | trace、Hilbert–Schmidt 内積と明示空間版 | O の parity と P。既存の局所 Product / Trace 等も必要 |
| State の基底処理 | expand_in、coordinates_in | graded の作用評価能力。状態 tensor の展開だけには交差符号を追加しない |
| Operator の基底処理 | expand_in、coordinates_in と明示空間版 | O の parity、係数の Neg。intrinsic な整列自体には P は不要 |
| export / normalize | 共通の Export・norm・正規化・期待値 API | 呼び出す core 評価経路の能力を引き継ぐ |

### 3.5. 利用側の移行

既存 API・trait bound の互換性は維持条件にしない。ordinary 専用の独自代数は、複合評価に必要な parity・P を明示実装する。自明な grading を採用する代数では全状態・演算子を even、P を I とする。graded を使う側は、空間・状態を `TensorKind::Graded` で構築し、同じ評価・export API を呼ぶ。既存空間の kind を書き換える変換は提供せず、構築時に種類を選ぶ。

`SpaceError`・`ExpressionError` への variant 追加は、外部の網羅的な match を壊し得る。0.2.0 の移行案内でこの変更を記録する。既存のエラー型を一律に `non_exhaustive` へ変えるかは別の API 判断とし、本案では前提にしない。

## 4. 式の構築と作用の空間

状態も `State::tensor_with_kind(kind, states)`・`tensor(states)`・`tensor_g(states)` の三種類に揃え、空間と状態構文を同じ kind で構築する。fallible 入力には `try_tensor_with_kind`・`try_tensor`・`try_tensor_g`、二項合成には空間と同じ `tensor_product_with_kind`・`tensor_product`・`tensor_product_g` を用意し、共通処理へ委譲する。`StateNode::Tensor` にも kind を保持するか、それと等価な内部 context を保持し、評価・再構築・基底列挙で失わないようにする。

`QuantumSystem` にも `state_tensor_with_kind`・`state_tensor`・`state_tensor_g` と対応する `try_` 版を用意する。配置も行う入口は `state_tensor_at_with_kind`・`state_tensor_at`・`state_tensor_g_at` とし、種類指定版は kind を第一引数に受け取る。すべて共通 constructor へ委譲する。`try_` の最初の入力エラー優先・以降の Iterator を読まない規則を維持する。全入力が成功した後に、重複と種類の混在を検証する。

演算子には独立した tensor 操作・節点を追加しない。因子順の配置積 `E_0(A_0) E_1(A_1) ...` が、その空間上の演算子 tensor を表す。未束縛演算子は kind を決めず、既存どおり明示空間、既存束縛、作用先状態の順で全空間を選択する。複数の束縛済み演算子を合成するときは、新しい空間互換性の検査で ordinary / graded の不一致も拒否する。

### 遅延した Act の scope を保存する

現在の `StateNode::Act(ActExpr { state, operator })` に、作用を構築した時点で選択した `HilbertSpace` を追加する。共通パーツ `ActExpr` を一律に変更せず、複合式側で包む案を推奨する。これは演算子の旧束縛を子節点へ残すこととは異なり、確定済みの状態への作用の意味を保存する変更である。

例えば odd な B と odd な状態 `s_a` に対し、次の二つは graded では異なる。

```text
tensor_g(s_a, s_b.act(B_b))      // B は H_b 内で既に作用を構築済み
tensor_g(s_a, s_b).act(B_b)      // B は H_a ⊗g H_b に埋め込んで作用
```

後者だけが `s_a` の parity による負号を持つ。評価時に外側の `State.space` だけを見る実装では、この区別を失う。評価タスクへ各 Act の選択済み空間を渡し、内側の Act ではその scope に切り替える。後から tensor 合成した作用は外側の空間を用いる。展開・再構築でも scope を保存する。

演算子の未束縛構文は選択した全空間に対して解釈できる従来契約を維持する。ただし、作用評価済み・基底展開済みの結果を単に再束縛しても、元の未評価式が別空間で持つ作用を再現するとは限らない。自動的な再評価・source 保持は追加しない。

## 5. 合成・評価の具体案

### 5.1. 二つの内部表現を区別する

| 内部表現 | 内容 | 用途 |
| --- | --- | --- |
| Intrinsic な因子順の homogeneous な項 | 符号係数と `A_0 ⊗g ... ⊗g A_n` | graded 演算子基底への展開・座標取得 |
| 通常 tensor 表現の局所因子列 | parity string を挿入済みの通常の作用 | 状態への作用、full trace、行列との照合 |

異なる private な型・走査入口で区別し、通常表現の項を再び graded の `OpTerms` へ流さない。これは一時的な評価表現であり、全式共通の完全展開形・公開の変換 API・永続キャッシュは要求しない。

homogeneous な intrinsic 項に対応する通常表現は次になる。因子は左から順に並べる。

```math
\rho(A_0\otimes_g\cdots\otimes_g A_{n-1})
=\bigotimes_{i=0}^{n-1}\left(A_i P_i^{q_i}\right),
\qquad q_i=\sum_{j>i}p(A_j)\pmod 2.
```

P は **A の右側** に置く。例えば `A ⊗g B` は `A P^p(B) ⊗ B` である。この順序を逆にすると odd 同士の符号が変わる。

### 5.2. 状態への作用

推奨する基本経路は、選択済み空間内の配置された homogeneous な primitive `A_j` を、評価タスクとして `P_0 ... P_{j-1} A_j` へ下ろし、既存の局所 Act を利用する方式である。even なら P は不要、ordinary 空間なら常に従来どおり局所作用だけとする。

局所演算子式が和・積・随伴を持つ場合も、和の分配、積の右からの順次作用、随伴の順序反転と局所 Adjoint を使い、primitive ごとに埋め込む。状態への作用のために O の `Product` を要求してはいけない。随伴の返す線型結合の各項は元の parity を保存する契約を利用する。

この経路が保存すべき homogeneous な多因子の作用則は次になる。通常 tensor 表現への変換後に、同じ指数の符号をさらに加えてはいけない。

```math
(\mathop{\otimes_g}_i A_i)(\mathop{\otimes_g}_i\psi_i)
=(-1)^{\sum_{i>j}p(A_i)p(\psi_j)}
\mathop{\otimes_g}_i(A_i\psi_i).
```

先行する状態が混合 parity・もつれた和でも、P の実際の作用を行うため、一つの定数符号へ誤って置換しない。先行区間を表すタスクまたは private な区間記述で P の挿入を遅延させ、全状態の積和展開や全基底列挙を必要にしない。先行区間と対象因子以外は共有できるが、P が作用する先行因子は spectator として共有するだけでは済まない。

S の `HomogeneousParity` は parity の問い合わせ・契約検証・homogeneous 状態に対する直接符号処理に利用できる。基本の P 作用経路自体は S の parity を読まずに実現できるため、共通の `evaluate` の S の bound は既存 Act を基本とし、直接符号処理を導入するときだけ追加する。最初から二つの独立した本番評価器を保守する必要はない。

### 5.3. 積と随伴

演算子の構築順序は変更しない。graded の intrinsic な項走査では、局所式を必要な範囲で homogeneous 成分へ分け、選択空間の葉順へ安定的に整列させる。異なる因子の odd 成分同士が交差するたびに負号を一度加え、同じ因子の局所積は元の順序で処理する。Index の値の大小では整列しない。

多因子の規則は、homogeneous な A・B に対して次になる。

```math
(\mathop{\otimes_g}_i A_i)(\mathop{\otimes_g}_i B_i)
=(-1)^{\sum_{i>j}p(A_i)p(B_j)}
\mathop{\otimes_g}_i(A_iB_i).
```

局所積の非零な返却項は入力 parity の和を持つ。和は枝ごとに扱う。符号は内部の bool などで保持し、係数を出力するときに必要なら `Neg` を使う。既存の `OpTerms` の「後から現れた同じ Index を符号なしに集約する」処理を graded に流用しない。

`Operator::adjoint` は既存の遅延節点を維持し、`(XY)† = Y†X†` を使う。intrinsic 項を再び因子順へ整列すると、tensor 随伴の符号が生じる。

```math
(\mathop{\otimes_g}_i A_i)^\dagger
=(-1)^{\sum_{i<j}p(A_i)p(A_j)}
\mathop{\otimes_g}_i A_i^\dagger.
```

この整列を使う場合、同じ tensor 随伴の符号を別に足さない。通常表現に下ろした側は、通常の積・随伴を使うだけとする。ordinary の場合は両方の指数を符号処理に使用しない。

### 5.4. 内積と full trace

状態の内積は既存の正定値な縮約を再利用し、追加の Koszul 符号を入れない。graded の遅延作用を正しく評価してから縮約する。期待値、norm squared、正規化にも同じ経路を用いる。

full trace は両種類とも通常の非正規化 trace とする。graded では配置式を通常表現へ下ろした因子列について局所積・局所 Trace を処理する。構文上参照されない因子も P が入る可能性があるため、下ろす前に identity spectator と判定しない。

例えば二 mode では `c_1 = P_0 ⊗ a_1`、`c_1†c_1 = I_0 ⊗ n_1` である。前者の trace だけでは誤った spectator 処理でもゼロになり得るため、後者の trace や随伴を含む積も検証する。`Tr(XY)=Tr(YX)` に parity 符号を追加しない。supertrace と partial trace は追加しない。

通常表現を一時的な局所因子列として処理すれば、trace に S・状態基底・`usize` で表せる総次元を要求せず、現在の局所 Trace による計算を維持できる。Hilbert–Schmidt 内積も同じ trace 経路へ接続する。

## 6. 基底・座標・export

状態基底の tensor 自体には追加の符号を入れない。`basis_states` と状態展開結果が元の graded 空間・kind を保持するように変更し、座標取得前の作用を graded 対応経路に接続する。右端が最速の mixed-radix 規約、空の基底 `{1}`、非正規直交基底に対する既存契約は維持する。

演算子基底の要素は `E_0(b_0) ... E_n(b_n)` という因子順の配置積と定義する。これは intrinsic な graded tensor 基底である。`basis_operators` の列挙はこの積を作るだけなので、新たな局所積の評価や全行列は不要である。

共通の `expand_in` / `coordinates_in` は、選択空間が graded の場合、5.3 節の intrinsic な項走査から局所基底へ展開する。P string を含む通常の因子列の係数を、そのまま graded 演算子基底の座標として返してはいけない。

例えば `c_1` の intrinsic な座標は `I ⊗g a` に対する係数1である。一方、状態基底上の通常の行列は `P ⊗ a` になる。`c_1` を graded 空間上の配置積 `P_0 a_1` として再構築すると P を二重に含み、異なる演算子になってしまう。展開結果は intrinsic な因子順の配置積で再構築し、元の束縛先を保持する。再展開後も作用・座標が一致することを検証する。

`qmb-export` の行列化は既存どおり「状態基底列挙 → 作用 → 状態座標 → 列格納」とする。空間の種類に対応した共通の作用・座標 API を呼び、backend 側には符号処理を追加しない。演算子基底の座標と状態基底上の行列は別の出力である。

`Export` 自体を両種類に対応させ、共通評価に必要な bound へ変更する。専用 adapter は設けず、Vec・ndarray・sprs の列構築・storage を共有する。vector 出力も同じ API から graded の遅延作用を評価する。metadata は `HilbertSpace` を保持するため kind を失わず、shape・格納順序を変更しない。

`qmb-normalize` も共通の内積・期待値を呼ぶ既存メソッドの bound を更新し、平方根・除算・零 norm の処理を共有する。core だけの対応で完了とせず、公開の導出操作まで接続する。

## 7. 変更箇所と実装順序

以下の段階を一つの機能として進める。途中段階では未対応の graded 評価を明示的に拒否し、ordinary として評価しない。

| 段階 | 主な変更箇所 | 到達条件 |
| --- | --- | --- |
| 1. 空間・局所能力 | `qmb-core/src/space.rs`、`ops.rs`、`lib.rs` | kind、正規化 metadata、混在拒否、互換性・部分系、parity / P の契約と再公開 |
| 2. 構築・scope | `expr/composite.rs`、`system.rs`、`basis/enumeration.rs` | kind 指定の状態構築、Act の選択済み空間、再構築・走査・破棄の更新 |
| 3. 作用と内積 | `evaluate/composite.rs`、`evaluate/inner_product.rs` | scope 内の P 挿入、混合 parity、期待値、必要な能力を要求する統一 API |
| 4. 演算子の項処理 | `expr/traversal.rs`、`evaluate/trace.rs` | intrinsic / 通常表現を区別し、積・随伴・通常 trace を検証 |
| 5. 基底と数値接続 | `basis/expansion/composite.rs`、`basis/coordinates.rs`、`qmb-export`、`qmb-normalize` | 展開・再展開・座標・密疎行列・正規化の graded 対応 |
| 6. 具体例と文書 | core のテスト用代数、`examples`、`docs/design/`、`docs/architecture.md`、`docs/future-work.md` | 独立した二・三 mode の参照計算と既存 ordinary 回帰、公開 API・移行例の記録 |

作用・trace・基底の全経路で空間検証を先に行い、Zero・零係数・identity でも不正な Index や空間を見逃さない。`BasisError`・`ExportError` は既存の複合式エラーの伝播を利用し、新しい種類の拒否も外側で失わないようにする。

`qmb-basic` の Index 発行と `qmb-geometry` の Chain に grading の所有・符号計算を追加しない。利用側が kind を指定して空間を構築し、従来の Site 順の外部空間列を接続する。`qmb-spin-half` の既存モデルは ordinary を明示する現在の挙動を維持する。共通の複合評価へ接続するため、S・O を even、P を I とする trivial grading を明示実装し、fermion の局所規則を流用しない。core / geometry の `no_std + alloc` を維持し、新規の本番依存を必要としない。

## 8. 受入条件と検証計画

本節は**将来の実装時に実施する条件**である。本ドキュメントの作成時点で graded の実装・回帰テストが通ったという記録ではない。

### 8.1. 空間と構築

- 両種類について、同じ葉順の再結合が互換、逆順が非互換、重複 Index がエラーとなる。
- `T_o()` と `T_g()`、`T_o(E,a)` と `T_g(a)`、`T_o(E,T_g(a,b))` と `T_g(a,b)` がそれぞれ互換となる。保持した子・kind の読み取りは元の木を返す。
- `T_o(a,b)` と `T_g(a,b)` が非互換、実質的な混在構築がエラーとなる。空・単因子 wrapper を挟んでも同じ判断になる。
- `T_g(a,c)` は `T_g(a,b,c)` の部分系、`T_o(a,c)` と逆順の `T_g(c,a)` は部分系でない。単因子・空の部分系規則も確認する。
- `try_` の入力エラー優先・消費停止、束縛の集約規則、Zero でも参照検証する規則を維持する。
- `tensor` と ordinary 指定の `tensor_with_kind`、`tensor_g` と graded 指定の `tensor_with_kind` が、それぞれ同じ空間・状態を構築する。三入口と fallible・二項・QuantumSystem の委譲経路で、同じ検証・エラー規則が適用されることを確認する。
- 既存と同じ2万段・128 KiB スタックの条件で、kind の取得・構築・比較・評価・最終所有者の破棄を検証する。単因子 wrapper の構築が繰り返し全葉走査を起こさないことも確認する。

### 8.2. 符号と評価

テスト用単一 mode を `|0>` even、`|1>` odd、`a=|0><1|`、`a†=|1><0|`、`P=diag(1,-1)` と定める。基底順序は二 mode なら `|00>,|01>,|10>,|11>` とする。参照行列は evaluator・export を使わずに小さな Kronecker 積から独立に組み立てる。

- `c_0=a⊗I`、`c_1=P⊗a` より `c_1|11>=-|10>`、`c_0 c_1|11>=-|00>`、`c_1 c_0|11>=|00>` を確認する。ordinary では異なる因子の a 同士が可換となる。
- 同一 mode の局所 CAR と異なる mode の CAR を区別し、`{c_i,c_j†}=δ_ij I`、`{c_i,c_j}=0` を確認する。
- 三 mode の離れた因子で中間の P が残り、空状態・占有状態・混合 parity の重ね合わせ・もつれた状態へ正しく作用する。
- `(AB)ψ=A(Bψ)`、`(AB)†=B†A†`、二回の随伴、`⟨φ,Aψ⟩=⟨A†φ,ψ⟩` を ordinary / graded の両方で確認する。積の基底要素に odd 成分が2個以上あるケースを含める。
- 4 節の Act scope の二経路に期待される符号差があり、先に評価してから tensor 合成する経路と遅延したまま合成する経路が一致する。
- 空・単因子・再結合で結果が変わらず、内積は正定値で余分な parity 符号を持たない。
- 通常 trace の巡回性、`Tr(I)=2^N`、`Tr(c_1†c_1)=2`（二 mode）、Hilbert–Schmidt 内積が独立な行列と一致する。partial trace はこの受入条件に含めない。

### 8.3. 基底と全利用経路

- 状態の直接展開と再展開、graded 演算子基底の直接座標と再展開後の座標、元の演算子と展開後の作用・trace が一致する。
- 特に `c_1` の intrinsic な `I ⊗g a` という座標と、通常の行列 `P ⊗ a` を別々に検証する。両方を同じ配列と比較するテストにしない。
- even / odd を含む複数の局所基底、混合 parity の式、非正規直交な状態基底を含める。基底列挙順序・要素数・空の基底を維持する。
- Vec・ndarray・sprs が同じ状態基底行列を出力し、kind・因子順序の metadata を保持する。export 側で再度符号を入れない。
- graded の期待値・正規化を行列計算と照合する。`Complex<f32>`・`Complex<f64>` の両方で既存要件の許容誤差を用い、近似的な零による成分の黙った削除をしない。
- 独自代数・テスト用代数・利用例を共通の grading 契約へ移行し、`Product` を持たない作用、操作ごとの最小係数能力、borrowed primitive、共有・疎性を検証する。spin-half の全 primitive が even で P=I となり、両種類の空間で同じ結果を返すこと、同じ fermion primitive でも ordinary / graded で因子間の符号が変わることを確認する。
- 通常の構築・作用・内積・trace が全基底・密行列を生成しないことを、既存の allocation / capability テストを拡張して確認する。parity string の長さに応じた作業量は必要だが、状態全体の積和展開を共通の前処理にしない。

実装時の基本検証は `cargo test --workspace --all-features`、`cargo check -p qmb-core -p qmb-geometry --no-default-features`、`cargo fmt --all -- --check` とする。後者の check はホスト上の確認であり、専用 no_std target の検証とは区別する。MSRV の実ビルド・専用 target・性能計測は、既存の未検証事項と区別して実施有無を記録する。

## 9. 実装着手前に固定する詳細

本案では「kind は空間」「parity / P は primitive」「符号は core」「混在は同一視後に拒否」「空間の種類に従う統一評価 API」「作用の scope 保存」「intrinsic 座標と通常表現の区別」を推奨する。数学的規約は既存設計を引き継ぎ、新たな規約へ変更しない。

利用者との合意により、構築入口は `tensor_with_kind`・`tensor`・`tensor_g` の三種類とし、空・単因子の同一視後にも残る ordinary / graded の混在は拒否する。残る詳細は、各共通操作の最小 trait bound、intrinsic 項走査のデータ構造、P の区間タスクの共有方法、エラーの具体名である。API を確定する際には、構築・作用・trace・基底・export・normalize を通した短い利用例と最小 trait bound を先にレビューする。Fermion の製品実装に進む前に、この共通機構をテスト用代数で検証する。

設計文書の「実装済み API」や future work の完了状態は実装後に更新する。今回の文書追加だけをもって、v2 の Graded tensor product または Fermion の実装を完了扱いにはしない。

## 10. 本ドキュメント作成時の確認

現行ソースと上記の設計文書を照合し、相対ファイルリンク、コードフェンス、差分の空白を確認した。補助的に Python 標準ライブラリの整数行列計算を使い、局所演算子 `{I,P,a,a†}` の二・三 mode の組合せについて、通常表現への変換と graded の積・随伴を4,432ケース照合した。二 mode の CAR、`Tr(c_1†c_1)=2`、8.2 節に記載した占有状態への作用の符号も確認した。

これは初回作成時の式の整合性の確認であり、Rust の graded 実装を検証したものではない。その後の API 方針の更新でも変更は計画文書のみで、Rust のビルド・テストは実行していない。
