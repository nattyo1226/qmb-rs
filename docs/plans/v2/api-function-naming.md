# 関数名・公開 API 名の整理候補

作成日: 2026-10-30
更新日: 2026-10-30

現在の実装を読み、同種の操作で名前の規則が揃っていない箇所と、公開 API として長くなっている名前を記録する。実装名・入力・動作はソースで確認した事実であり、使いにくさの評価・優先度・新しい名前は検討案である。第1〜3節は現行 API に対する候補と解消済み項目、第5節は採用状況と推奨案を記録する。N1 の通常の加算・乗算については、式型の名前付き構築を列の `sum` / `product` に絞り、二項は標準の `+` / `*` を使う方針を採用済み・未実装とする。N3 の基底処理名、N4 の次元問い合わせ名、N5 の Chain の命名方針、L1 の Hilbert–Schmidt 内積名、L4 の HilbertSpace のテンソル構築名、L5 の局所式の項列構築も採用済み・未実装とする。N6 のモデル名、L3 の QuantumSystem の演算子積名、L6 の正規化期待値名は現状維持に決定した。状態構築の変更は実装済みで、複合式の設計文書に記録する。内部関数 I1〜I4 も提示した推奨名を採用済み・未実装とする。primitive の独自 trait とメソッド `Product::product` は現行維持に決定した。改名に伴う公開 API の互換性維持は不要とする。N1〜N6・L1〜L7・I1〜I4 の採否と互換性方針は決定済みである。追加で採用した、型情報の QuantumSystem trait と構築用 ExpressionFactory struct の役割分担、型エイリアスとモデル構築入口の分離は第6節に記録する。型注釈・ジェネリック API の移行方針も同節で決定済みとする。いずれも実装時期は未定である。

確認対象は `crates/` 以下の `qmb-core`・`qmb-basic`・`qmb-geometry`・`qmb-spin-half`・`qmb-spin-one`・`qmb-fermion`・`qmb-normalize`・`qmb-export` の関数・メソッド・公開 trait と、`examples/src/bin/` の呼び出しである。非公開 module 内でも公開型の inherent method や再公開された trait は公開 API に含める。内部関数は後段に分ける。

優先度「高」は日常的な構築・計算で長さや名前の衝突が目立つ候補、「中」は API 間の規則を揃える候補、「低」は現行名を残す判断も有力な候補とする。識別子の長さは `_` を含み、型名・module path は含めない。文字数だけを改名の基準にはしない。

明示空間への移行後の実装を再確認した。以下の ID は従来の議論を追えるよう維持し、解消済みの項目は改名の実施対象から外す。

| 対象 | 現在の状態 |
| --- | --- |
| N1 の通常の加算・乗算 | 名前付き構築を列の `sum` / `product` に絞り、式型の二項 inherent `add` / `product` を廃止する方針を採用済み・未実装。標準の `+` / `*` は維持 |
| N3 の基底処理と空間指定 | 基底のみは `expand` / `coordinates`、明示空間版は `expand_in` / `coordinates_in` を採用済み、改名は未実装 |
| N4 の次元問い合わせ | `dim` / `dim_state` / `dim_op` / `dim_with` / `dim_local` を採用済み、改名は未実装。`Dimension::dim` は維持 |
| N5 の Chain の接続方法と空間指定 | `pairs_by_pos` → `zip`、`_on` → `_in` を採用済み、改名は未実装。公開名に `leaf` は追加しない |
| N6 のモデル名 | 現状維持を採用済み。`hubbard` / `symmetric_hubbard` は改名しない |
| L1 の Hilbert–Schmidt 内積 | `hs_inner_product` / `hs_inner_product_in` を採用済み、改名は未実装 |
| N2 の状態側の種類指定・修飾順序 | 明示空間から種類を取得する構築へ移行して解消済み |
| L2 の長い状態構築名 | `tensor_in` / `try_tensor_in` / `product_state` への移行で解消済み |
| L3 の QuantumSystem の演算子積 | 現状維持を採用済み。`operator_product` / `try_operator_product` / `operator_product_at` は改名しない |
| L4 の空間のテンソル構築名 | HilbertSpace の二項を `tensor*`、列を `tensor_all*`、種類指定の接尾辞を `_with` とする方針を採用済み、改名は未実装。State 側は現行維持 |
| L5 の局所式の項列構築 | 両局所式型の `from_terms(impl IntoIterator<Item = Term<Self, C>>)` を採用済み・未実装。改名と入力の一般化を行う |
| L6 の正規化期待値 | 現状維持を採用済み。`normalized_expectation` / `normalized_expectation_in` と二つの trait 名は改名しない |
| I1〜I4 の内部関数 | 提示した推奨名を採用済み、改名は未実装。公開範囲と処理内容は維持する |
| N1 に付随する trait の整理 | 現状維持を決定済み。独自の `Product<C>` と `Product::product` を残し、標準の `core::ops::Mul` には置き換えない |
| 改名時の互換性 | 互換性維持は不要。旧名の互換 alias・wrapper を設けず新名へ移行する |
| QuantumSystem とモデル構築入口 | 型情報の QuantumSystem trait、構築用 ExpressionFactory struct、具体的な型エイリアスと各 Models 型への分離を採用済み・未実装。既存の `*Of` 型注釈も維持。詳細は第6節 |

現在の状態構築契約は[複合式の設計](../../design/composite-expression.md#121-構築入口)を参照する。所属空間の必須化と、その空間からの種類取得を前提に API 整理を続ける。

## 1. 名前の規則が揃っていない箇所

### N1. 二項演算と列の集約の区別（採用済み・未実装）

| 操作 | 二項版 | 列の集約版 |
| --- | --- | --- |
| 複合状態・演算子の和 | `State::add`、`Operator::add` | `sum` / `try_sum` |
| 局所状態・演算子の和 | 二項の inherent method はなく、標準の `+` を使う | `sum` / `try_sum` |
| 演算子の積 | `LocalOperator::product`、`Operator::product` | `product_all`、複合式には `try_product_all` もある |
| 空間のテンソル積 | `HilbertSpace::tensor_product` | `HilbertSpace::tensor` |
| 状態のテンソル積 | 専用の二項版はなく、2要素の列で構築 | `State::tensor_in` / `try_tensor_in`、QuantumSystem にも同名の入口 |
| QuantumSystem の演算子積 | 二項専用の入口はない | `operator_product` / `try_operator_product` |

列を表す印が演算子の積では `_all`、空間のテンソル積では無印になり、QuantumSystem の演算子積では `_all` が消える。状態のテンソル構築は明示空間を受け取る列の入口へ統一済みで、`_in` は所属空間を表す。各シグネチャは明確だが、別の型で同じ操作を探す際に規則を推測しにくい。

**採用決定（2026-10-30、方針更新）**：通常の加算・乗算では、式型の名前付き構築メソッドを列の `sum` / `product` に絞る。二項演算には標準の `+` / `*` を使い、必要なら2要素の列を集約入口へ渡す。

| 操作 | 二項演算 | 名前付きの列集約 |
| --- | --- | --- |
| 加算 | 標準の `+` | `sum` / `try_sum` |
| 乗算 | 標準の `*` | `product` / 複合式の `try_product` |

`State` / `Operator` の inherent `add` と、`LocalOperator` / `Operator` の二項 inherent `product` は廃止する。代替の二項 inherent `add` / `mul` / `sum` / `product` は設けない。四つの式型の列集約 `sum` / `try_sum` は維持し、LocalOperator / Operator の `product_all` は `product`、Operator の `try_product_all` は `try_product` に改名する。通常の和・積に `sum_all` / `product_all` は残さない。以前の二項 `add` / `mul` を名前付きで残す案は、この決定で置き換える。

標準の `core::ops::Add::add` / `Mul::mul` と、既存の所有・借用の `+` / `*` 実装は維持する。二項の借用入力を保持したい呼び出しでは `&a + &b` / `&a * &b` を使える。標準 trait のメソッド自体を削除・改名するものではない。列集約の入力型・戻り値・能力制約・エラーの優先順位を維持し、新たな iterator 集約 trait の実装は追加しない。`State::sum(&space, states)` / `try_sum(&space, states)` は所属空間を明示する契約を維持する。

LocalOperator の列の積には、空入力を恒等演算子にするための `O: Identity` 制約が引き続き必要である。標準の二項 `*` は Identity 制約なしで使える契約を維持する。したがって、旧二項呼び出しをすべて列集約へ機械的に置き換えない。列集約には空入力・単項入力も含み、積の式木は既存の二項 ProductExpr を左から畳み込む構造を維持する。

primitive の独自 trait `Product<C>` と、その二項メソッド `Product::product` は現行維持とする。係数型 `C` を primitive 型から独立に指定でき、局所規則が線型結合を返す契約を維持する。標準の `core::ops::Mul` への置き換えや、primitive のメソッドの `mul` への改名は行わない。これで N1 に付随する trait 名・メソッド名の保留も解消する。採用は方針の確定であり、この時点では実装を変更していない。

HilbertSpace のテンソル積の二項 `tensor` / 列 `tensor_all` は、通常の乗算とは別に L4 として採用済み・未実装とする。空間側の現行 `tensor` は列を受け取るため、二項・列を同時に移行する。状態側の `tensor_in` / `try_tensor_in` は実装済みの名前を維持し、`tensor_all*` への改名は行わない。

根拠: [State の構築](../../../crates/qmb-core/src/expr/composite.rs#L262)、[Operator の積](../../../crates/qmb-core/src/expr/composite.rs#L514)、[局所演算子の積](../../../crates/qmb-core/src/expr/local.rs#L233)、[空間のテンソル積](../../../crates/qmb-core/src/space.rs#L181)、[QuantumSystem の積](../../../crates/qmb-core/src/system.rs#L150)。

### N2. 状態側の種類指定・修飾順序（解消済み）

状態側の `_g` / `_with_kind` と、配置を伴う `state_tensor*` の入口は撤去した。配置済み状態の合成は `State::tensor_in(&space, states)` / `try_tensor_in`、primitive の葉順の配置は `QuantumSystem::product_state(&space, primitives)` を使う。QuantumSystem にも `tensor_in` / `try_tensor_in` を提供している。

テンソル積の種類の情報源は `HilbertSpace` に集約した。State は所属空間の共有ハンドルを持ち、状態の Tensor 節点は子の列だけを保持する。遅延作用の StateAction は構築時に選択した対象空間を保持する。このため、状態 API の修飾順序を変えたり `_at_g` / `_with` を新設したりする必要はない。

HilbertSpace の `tensor_g` / `tensor_product_g` / `tensor_with_kind` は引き続き存在する。空間側には配置を表す `_at` はなく、今回問題にした修飾順序の揺れは残っていない。`_g` は L4 の採用内容に従って維持し、空間側の名前の長さ・二項／列の区別は N1・L4 で扱う。

根拠: [QuantumSystem の構築入口](../../../crates/qmb-core/src/system.rs#L62)、[State の構築](../../../crates/qmb-core/src/expr/composite.rs#L310)、[状態式と作用節点](../../../crates/qmb-core/src/expr/composite.rs#L95)、[空間の構築](../../../crates/qmb-core/src/space.rs#L181)。

### N3. `_in` が基底と空間の両方を表す（採用済み・未実装）

| 現行名 | `_in` が指す入力 |
| --- | --- |
| `LocalState` / `LocalOperator` / `State` / `Operator` の `expand_in(&basis)`、`coordinates_in(&basis)` | 基底 |
| `State` / `QuantumSystem` の `tensor_in(&space, states)` / `try_tensor_in(&space, states)` | 明示する所属空間 |
| `State::act_in` / `expectation_in`、`Operator::trace_in` / `hilbert_schmidt_inner_product_in`、`normalized_expectation_in` | 明示する作用空間 |
| `Operator::expand_in_space(&basis, &ambient)` / `coordinates_in_space(&basis, &ambient)` | 基底に加えて明示する作用空間 |
| `MatrixArgs::in_space(&ambient)` | 行列出力で優先する作用空間 |

新しい状態構築でも `_in` は空間を表すため、基底処理との用語の違いが残っている。基底処理は基本名の時点で `_in` を使い、明示空間を加えると `_in_space` になる。特に `coordinates_in_space` は20文字で、他の明示空間版から名前を推測しにくい。

**採用決定（2026-10-30）**：基底を受け取る基本名を `expand(&basis)` / `coordinates(&basis)`、明示空間版を `expand_in(&basis, &ambient)` / `coordinates_in(&basis, &ambient)` に揃える。`Basis::expand` とも語彙が揃う。実装済みの状態構築名 `tensor_in` / `try_tensor_in` と、葉順の配置を表す `product_state` は維持し、N3 では改名しない。同名メソッドの引数数が変わるため、局所式・複合式をまとめて移行する必要がある。`MatrixArgs::in_space` は builder の役割が明瞭であり、現行名を維持する。

対象は四つの式型（`LocalState`・`LocalOperator`・`State`・`Operator`）の基底のみの入口と、Operator の既存の明示空間版とする。引数順序・戻り値・能力制約・基底の借用・空間選択の契約は変更しない。`act_in`・`expectation_in`・`trace_in` と他の明示空間の計算入口も維持する。詳細な対応は第5.4節に記録する。採用は命名方針の確定であり、この時点では実装を変更していない。

根拠: [局所式の基底処理](../../../crates/qmb-core/src/basis/expansion/local.rs#L93)、[複合演算子の基底処理](../../../crates/qmb-core/src/basis/expansion/composite.rs#L241)、[作用・期待値](../../../crates/qmb-core/src/evaluate/inner_product.rs#L327)、[trace・内積](../../../crates/qmb-core/src/evaluate/trace.rs#L60)、[正規化期待値](../../../crates/qmb-normalize/src/expectation.rs#L23)、[MatrixArgs](../../../crates/qmb-export/src/args.rs#L24)。

### N4. 次元を表す `dim` / `dimension` / `basis_len`（採用済み・未実装）

primitive は `Dimension::dim()`、空間・座標は `dimension()` / `local_dimension()` を使う。また `state_basis_len::<S>()` / `operator_basis_len::<O>()` は基底を受け取らず、内部で primitive の `dim()?.get()` を返している。

同じ次元の問い合わせに略語・完全形・長さという三つの語彙がある。

**採用決定（2026-10-30）**：以下の次元問い合わせを `dim` に揃える。状態・演算子を選ぶ自由関数は `dim_state` / `dim_op`、明示する局所次元を受け取る入口は `dim_with`、局所基底の次元を返す入口は `dim_local` とする。局所次元の語順は `dim_local` を採用する。

| 提供型・入口 | 現行名 | 採用名 | 状態 |
| --- | --- | --- | --- |
| Dimension trait | `dim()` | `dim()` | 維持 |
| 自由関数 | `state_basis_len::<S>()` | `dim_state::<S>()` | 採用済み・未実装 |
| 自由関数 | `operator_basis_len::<O>()` | `dim_op::<O>()` | 採用済み・未実装 |
| HilbertSpace | `dimension::<P>()` | `dim::<P>()` | 採用済み・未実装 |
| Coordinates | `dimension()` | `dim()` | 採用済み・未実装 |
| HilbertSpace | `dimension_with(local)` | `dim_with(local)` | 採用済み・未実装 |
| Coordinates / VectorMetadata / MatrixMetadata | `local_dimension()` | `dim_local()` | 採用済み・未実装 |

`dim_state` / `dim_op` は基底を受け取らず、指定した primitive 型が張る線型空間の次元を返す。`dim_op` は演算子空間自体の次元であり、完全な演算子空間なら状態空間の次元 `d` に対して `d²` となる。二つの自由関数は `usize` を返す便宜を維持し、削除しない。primitive の戻り値の `NonZeroUsize` と、空間・自由関数の `usize` の違い、overflow の契約も維持する。

`dim_local()` は各因子で使う一様な局所基底の要素数を `Option<usize>` で返す。Coordinates とそのベクトル出力の VectorMetadata では使用する局所基底、MatrixMetadata では入出力に共通の状態基底の次元を表す。演算子基底の座標をベクトル出力する場合も、その基底の次元を維持する。ゼロ因子の単位空間では `None` を返し、座標空間全体の `dim()` は `Ok(1)` となる。

`Coordinates::len()` は保存された非零座標の数、改名後の `dim()` は座標空間全体の次元であり、この区別を維持する。primitive の次元は enum の variant 数ではなく張る線型空間の次元を表す。入力型・引数順序・戻り値・能力制約は変更しない。`Dimension` trait の名前と必須メソッド `dim` は維持する。完全形の `dimension` に統一する旧推奨案は採用しない。採用は命名方針の確定であり、この時点では実装を変更していない。

根拠: [Dimension](../../../crates/qmb-core/src/ops.rs#L81)、[HilbertSpace](../../../crates/qmb-core/src/space.rs#L149)、[基底長の自由関数](../../../crates/qmb-core/src/basis/enumeration.rs#L162)、[Coordinates](../../../crates/qmb-core/src/basis/coordinates.rs#L66)、[ベクトル出力の局所次元](../../../crates/qmb-export/src/metadata.rs#L51)、[行列出力の局所次元](../../../crates/qmb-export/src/metadata.rs#L115)。

### N5. Chain の接続方法と空間指定の命名（採用済み・未実装）

`nn_site_pairs` は Site 番号対、`nn_site_pairs_with` は与えた要素への参照対、`nn_pairs_with` は最近接 Site 間で局所空間の全組み合わせを返す。`nn_pairs_by_pos_with` は各 Site 内の局所空間を保存順に zip し、短い側で止まる。空間を入力する現行の `_on` 版にも同じ接続方法の違いがある。

**採用決定（2026-10-30）**：Site 単位は `nn_site_pairs*`、全組み合わせは `nn_pairs*` を維持し、位置対応を `nn_zip*` とする。Site 順の列を受け取る `_with` は維持し、HilbertSpace を受け取る `_on` は N3 の明示空間の命名に合わせて `_in` に変更する。`nn_site_pairs_on` も `nn_site_pairs_in` へ移行する。公開名に木構造由来の `leaf` を追加する旧案は採用しない。位置対応版は `nn_pairs_by_pos_with`（20文字）から `nn_zip_with`（11文字）へ短縮する。

Chain の `_in` は、指定空間に保存されたトップレベル分割を Site として読む意味とする。局所空間一つは一つの Site、空のテンソル積はゼロ Site として扱い、ネストした grouping は各 Site の内部に保持する。`_with` と `_in` の入力の解釈を統合せず、空間を組み直す変更も行わない。

入力型・引数順序・戻り値・能力制約は変更しない。Site 数を iterator 返却前に検証すること、全組み合わせ版・zip 版では空 Site が位置を保ったまま局所空間の対を生成しないこと、各 Site 内を保存順にたどること、全組み合わせ版は左側を外側として列挙すること、zip 版は短い側で止まることを維持し、実装時に rustdoc へ明記する。詳細な対応は第5.4節に記録する。採用は命名方針の確定であり、この時点では実装を変更していない。

根拠: [Chain の公開入口](../../../crates/qmb-geometry/src/geometry/chain.rs#L39)、[局所空間の対を列挙する実装](../../../crates/qmb-geometry/src/connection.rs#L24)。

### N6. モデル名の修飾語の位置（現状維持・決定済み）

`Fermion::hubbard` に対する派生名は `Fermion::symmetric_hubbard` で、操作名の後に修飾語を置く `_g` / `_in` 系と語順が異なる。

**採用決定（2026-10-30）**：現状維持とし、`hubbard` / `symmetric_hubbard` を残す。モデル名の形容詞は API の入力条件を示す接尾辞とは役割が異なるため、自然なモデルの呼び方を優先する。`hubbard_symmetric` への改名は行わず、N6 は改名の実施対象から外す。

根拠: [Hubbard の二つの入口](../../../crates/qmb-fermion/src/models.rs#L54)。

## 2. 長い公開名の短縮候補

N1〜N5 と同じ関数を含む。以下の候補名は相互に独立して採用できるとは限らず、二項／列の区別や修飾順序を決めてから確定する。

| ID・優先度 | 現行の関数名と長さ | 使いにくさと候補 |
| --- | --- | --- |
| L1・高 | `LocalOperator::hilbert_schmidt_inner_product`、`Operator::hilbert_schmidt_inner_product`（29）、`Operator::hilbert_schmidt_inner_product_in`（32） | 数値計算の呼び出しで長い数式名を毎回綴る。`hs_inner_product` / `hs_inner_product_in`（16 / 19）への改名を採用済み・未実装。`hs` の意味を rustdoc で示し、非正規化の `Tr(A†B)` という契約を残す。 |
| L2・解消済み | 旧 QuantumSystem の `state_tensor*` の種類指定付き構築名 | 長い旧名は撤去済み。現在は `tensor_in`（9）、`try_tensor_in`（13）、`product_state`（13）を提供し、所属空間を明示する。短縮用の別名や `tensor_all*` の追加は不要。 |
| L3・現状維持 | `QuantumSystem::operator_product_at`（19）、`try_operator_product`（20） | 現状維持を採用済み。`state_sum` / `operator_sum` と同様に対象を名前で明示する一貫性を優先し、`operator_product` を含む三入口の `operator_` を維持する。 |
| L4・高（空間側のみ） | `HilbertSpace::tensor_product_with_kind`（24） | 空間側の二項 `tensor_with(kind, &rhs)`、列 `tensor_all_with(kind, children)` への改名を採用済み・未実装。二項／列のファミリ全体を合わせて移行する。State の種類指定版は撤去済み。入力型 `TensorKind` と文書で `_with` の意味を明示する。 |
| L5・中 | `LocalState::from_linear_combination`、`LocalOperator::from_linear_combination`（23） | `from_terms`（10）への改名と、入力を `IntoIterator<Item = Term<Self, C>>` に広げる案3を採用済み・未実装。配列・iterator・既存の LinearCombination を受け取る。primitive の規則結果には `map(Self::primitive)` が引き続き必要。 |
| L6・現状維持 | `NormalizedExpectation::normalized_expectation`（22）、`NormalizedExpectationIn::normalized_expectation_in`（25） | 現状維持を採用済み。非正規化の `expectation` と正規化の有無を名前で区別し、既存の二つの trait 名も維持する。短縮案は採用しない。 |
| L7・中 | `Operator::coordinates_in_space`（20）、`Chain::nn_pairs_by_pos_with`（20） | 基底側の短縮は N3、Chain 側の `nn_zip_with` への改名は N5 として採用済み・未実装。名前だけを切り詰めず、それぞれのファミリを同時に揃える。 |

根拠:

- L1: [Hilbert–Schmidt 内積](../../../crates/qmb-core/src/evaluate/trace.rs#L33)。局所版と複合版の両方に同じ長い名前がある。
- L2: [実装済みの状態構築](../../../crates/qmb-core/src/system.rs#L62)と、現在の呼び出し例 [tensor_products](../../../examples/src/bin/tensor_products.rs#L17)。
- L3: [QuantumSystem の演算子積](../../../crates/qmb-core/src/system.rs#L105)。実際の呼び出しは [two_spins](../../../examples/src/bin/two_spins.rs) にもある。
- L4: [HilbertSpace](../../../crates/qmb-core/src/space.rs#L231)。
- L5: [LocalState](../../../crates/qmb-core/src/expr/local.rs#L137)、[LocalOperator](../../../crates/qmb-core/src/expr/local.rs#L204)、[LinearCombination の既存の from_terms](../../../crates/qmb-core/src/expr/parts/linear_combination.rs#L44)。
- L6: [正規化期待値の trait と実装](../../../crates/qmb-normalize/src/expectation.rs#L15)、[非正規化期待値](../../../crates/qmb-core/src/evaluate/inner_product.rs#L284)、利用例の [heisenberg](../../../examples/src/bin/heisenberg.rs#L22)。
- L7: [座標取得](../../../crates/qmb-core/src/basis/expansion/composite.rs#L321)、[Chain](../../../crates/qmb-geometry/src/geometry/chain.rs#L66)。

## 3. 内部関数の表記揺れ・処理内容が隠れる名前

| 現行名・公開範囲 | 確認した違いと整理候補 |
| --- | --- |
| `LocalOperator::evaluated`（`pub(crate)`）と公開の `LocalState::evaluate` / `State::evaluate` | いずれも元を変更せず評価済みの式を返すが、演算子側だけ過去分詞になっている。内部名も `evaluate` に揃える方針を採用済み・未実装。公開 API の追加は別の判断であり、現状の演算子に公開 `evaluate` があるとは扱わない。[局所評価](../../../crates/qmb-core/src/evaluate/local.rs#L138)、[状態評価](../../../crates/qmb-core/src/evaluate/composite.rs#L41)を参照。 |
| `StateExpr::same_node`（`pub(crate)`）と公開の `HilbertSpace::shares_node_with` | どちらも `Arc::ptr_eq` で共有節点の同一性を調べる。`shares_node_with` に揃えると数学的な同値判定と区別しやすい。[式](../../../crates/qmb-core/src/expr/composite.rs#L182)、[空間](../../../crates/qmb-core/src/space.rs#L171)を参照。 |
| `basis::enumeration::leaves`（`pub(super)`）と公開の `HilbertSpace::leaves` | 前者は葉をすべて `Vec` へ集め、後者は遅延 iterator を返す。内部の前者を `collect_leaves` にすると確保・集約を伴うことが伝わる。[集約](../../../crates/qmb-core/src/basis/enumeration.rs#L15)、[iterator](../../../crates/qmb-core/src/space.rs#L138)を参照。 |
| `qmb-export` の `coordinates::owned_entries`（`pub(crate)`） | 所有権の移動だけでなく、反復時に多桁の座標を flat index に変換し、検証・overflow 検出を行う。`into_flat_entries` とすれば、core の借用版 `Coordinates::flat_entries` と処理内容を対応づけられる。[変換](../../../crates/qmb-export/src/coordinates.rs#L21)、[借用版](../../../crates/qmb-core/src/basis/coordinates.rs#L100)を参照。 |

## 4. 整理時に維持する区別と検証

- **`try_` は単なる Result の印ではない。** `State::tensor_in` / `Operator::product_all` も `Result` を返すが、`try_tensor_in` / `try_product_all` は入力が `Result` の列である。最初の入力エラーで消費を止め、全入力成功後に構成を検証するという一貫した契約がある。名前を揃えるためにすべての fallible 操作へ `try_` を付けない。[契約](../../../crates/qmb-core/src/system.rs#L38)と[検証対象の既存テスト](../../../crates/qmb-core/tests/fallible_construction.rs)を参照。
- **名前が違っても、受け手や入力が違えば必要な区別がある。** `QuantumSystem::state` / `operator` は primitive と係数型を選択済みの構築入口、`LocalState::primitive` / `LocalOperator::primitive` は式型の入口である。`local.at(&factor)` と `State::place(&local, &factor)` / `Operator::place(&local, &factor)` はメソッドと関連関数の入口である。これらを長い `local_*` 名へ機械的に揃える優先度は低い。
- **短縮で数学的契約を変えない。** 左側の共役・随伴、積の右側が先に作用する順序、ordinary / graded の選択、正規化の有無、明示空間が古い束縛に優先する規則を維持する。局所式の恒等元と未束縛の複合恒等元、列集約の空入力・単項の扱いも維持する。
- **名前と能力を別に扱う。** primitive・Index の Clone 要求、数値能力、係数精度、no_std の境界を改名に伴って追加・変更しない。`normalized_expectation` の短縮のために `normalized().expectation()` へ置き換えると、平方根など現在不要な能力を要求するので同等の移行にはならない。
- **採用時に名前の衝突と利用例を確認する。** 改名時の互換性維持は不要とし、旧名の互換 alias は設けない。特に N1・N3 は既存名の引数形を変えるため、呼び出しの入力を確認して移行する。公開 trait の変更は独自代数の実装にも影響する。対応する crate のテスト・rustdoc・examples と、設計文書の API 表を合わせて更新・検証する。本書の調査ではソースと呼び出し例を確認し、実行テストは行っていない。

## 5. 各候補の推奨案

検討日: 2026-10-30

判断基準は、型をまたいで操作名から入力の形を推測できること、数学的な意味を失わないこと、頻出する呼び出しを短くすることの三点とする。既存名を残すだけで十分な項目も含めて結論を示す。第1〜3節の代替案のうち、以下で選ばなかったものは今回の推奨には含めない。根拠となる定義箇所は各候補のリンクを参照する。

### 5.1. 命名規則についての結論

| 対象 | 推奨 | 判断理由・選ばなかった案 |
| --- | --- | --- |
| N1 | **採用済み・未実装**。式型の名前付き構築は列の `sum` / `product` に絞り、二項の inherent `add` / `product` は廃止する。標準の `+` / `*` と primitive の `Product::product` は維持する。HilbertSpace は別途 L4 に従う。 | 名前付きの二項・列の区別を不要にし、列の積を `product_all` から `product` へ揃える。二項 inherent `mul` や `sum_all` は追加せず、独自の乗算規則 trait と標準 Mul の役割を維持する。 |
| N2 | **解消済み・現行維持**。状態側は明示空間を受け取る構築へ移行済み。HilbertSpace 側の `_g` は維持する。 | 種類の指定・保持を空間へ集約し、状態側の `_at_g` / `_with` の改名案は不要になった。状態の Tensor 節点に種類を戻さず、遅延作用の対象空間は維持する。 |
| N3 | **採用済み・未実装**。基底を指定する基本名を `expand` / `coordinates`、明示空間版を `expand_in` / `coordinates_in` とする。実装済みの `tensor_in` / `try_tensor_in` / `product_state` と `act_in`・`expectation_in`・`trace_in` と `MatrixArgs::in_space` は維持する。 | 式の構築・計算メソッドで `_in` の対象を所属・作用空間へ揃えられる。基底は必須引数なので無印でも情報を失わない。全計算名を `_in_space` に長くする案や、配置の `_at` を流用する案は採らない。builder は計算ではないため `in_space` を残す。 |
| N4 | **採用済み・未実装**。`dimension` → `dim`、`state_basis_len` → `dim_state`、`operator_basis_len` → `dim_op`、`dimension_with` → `dim_with`、`local_dimension` → `dim_local` とする。`Dimension::dim` は維持する。 | 次元の問い合わせを `dim` に揃え、状態・演算子・局所の指定を後ろに置く。二つの自由関数は削除せず、`Coordinates::len` と次元の区別、戻り値・overflow の契約も維持する。完全形への統一案は採用しない。 |
| N5 | **採用済み・未実装**。Site 単位は `nn_site_pairs*`、全組み合わせは `nn_pairs*`、位置対応は `nn_zip*` とする。列入力の `_with` は維持し、空間入力の `_on` は `_in` に揃える。 | 公開名に木構造由来の `leaf` を追加せず、`zip` で短い側までの位置対応を示す。`_in` は明示空間を表す N3 の方針に揃える。`nn_` と、列／空間の入力の区別は維持する。 |
| N6 | **現状維持・決定済み**。`hubbard` / `symmetric_hubbard` を残し、改名の実施対象から外す。 | モデル名の形容詞と API の入力条件を示す接尾辞は役割が異なる。自然なモデル名を変えて得られるのは一覧上の並びだけであり、短縮や意味の改善がない。 |

### 5.2. 長い公開名についての結論

| 対象 | 推奨する公開名 | 判断理由・移行範囲 |
| --- | --- | --- |
| L1 | **採用済み・未実装**。`hs_inner_product` / `hs_inner_product_in` | Hilbert–Schmidt の部分だけを `hs` とし、内積であることは省略しない。局所版・複合版を同時に変える。`inner_product` だけにする案は内積の種類を隠し、`dot` / `hs_inner` はさらに語彙を減らす必要がないため採らない。rustdoc の冒頭に正式名と `Tr(A†B)`、非正規化であることを残す。 |
| L2 | **解消済み・現行維持**。`tensor_in` / `try_tensor_in` / `product_state` | 型と種類指定を重ねた旧名は撤去済み。第5.3節の実装済み入口を使い、`tensor_all_at_with` などの旧推奨名は追加しない。 |
| L3 | **現状維持・決定済み**。`operator_product` / `try_operator_product` / `operator_product_at` を残す。 | `state_sum` / `operator_sum` と同様に、QuantumSystem の集約入口では対象を名前で明示する一貫性を優先する。`product` / `try_product` / `product_at` への短縮案は採用せず、L3 は改名の実施対象から外す。 |
| L4 | **採用済み・未実装（HilbertSpace 側のみ）**。二項は `tensor` / `tensor_g` / `tensor_with`、列は `tensor_all` / `tensor_all_g` / `tensor_all_with` | 空間の二項／列の区別と合わせて短縮する。状態側の `tensor_in` / `try_tensor_in` は維持する。空間の種類の暗黙推論は導入しない。 |
| L5 | **採用済み・未実装（案3）**。両局所式型で `from_terms` | 入力を `IntoIterator<Item = Term<Self, C>>` に一般化し、配列・iterator・既存の LinearCombination を一つの入口で受け取る。構築時の非評価、入力順・重複項・ゼロ係数の保持、数値演算・Clone の能力制約を要求しない契約を維持する。`from_combination` への単純改名と From 実装の追加は今回採用しない。 |
| L6 | **現状維持・決定済み**。`normalized_expectation` / `normalized_expectation_in` と二つの trait 名を残し、改名の実施対象から外す。 | 非正規化の `expectation` と式上で明確に区別できることを優先する。`mean` は短いが正規化規約が名前に表れず、`expectation` への単独改名は既存メソッドと競合する。名前の短縮のために平方根を要する正規化処理へ置き換えない。 |
| L7 | **採用済み・未実装**。`coordinates_in_space` → `coordinates_in` は N3、`nn_pairs_by_pos_with` → `nn_zip_with` は N5 として採用した。 | 基底側は基本版と明示空間版を同時に移行する。Chain 側は位置対応名と `_on` → `_in` を合わせて移行する。 |

**L1 の採用決定（2026-10-30）**：次の三入口を改名する。引数順序・戻り値・能力制約は変更しない。左側の随伴と非正規化の `Tr(A†B)`、明示空間が既存の束縛に優先する契約を維持する。採用は命名方針の確定であり、この時点では実装を変更していない。

| 提供型 | 現行名 | 採用名 | 状態 |
| --- | --- | --- | --- |
| LocalOperator | `hilbert_schmidt_inner_product(&rhs)` | `hs_inner_product(&rhs)` | 採用済み・未実装 |
| Operator | `hilbert_schmidt_inner_product(&rhs)` | `hs_inner_product(&rhs)` | 採用済み・未実装 |
| Operator | `hilbert_schmidt_inner_product_in(&rhs, &ambient)` | `hs_inner_product_in(&rhs, &ambient)` | 採用済み・未実装 |

**L5 の採用決定（2026-10-30）**：案3を採用し、LocalState / LocalOperator の `from_linear_combination` を次の入口に変更する。名前の短縮に加えて入力型を一般化する変更であり、採用済み・未実装とする。

```rust
pub fn from_terms(terms: impl IntoIterator<Item = Term<Self, C>>) -> Self
```

`Term<Self, C>` は係数 `C` と局所式 `Self` を所有する項を表す。既存の `LinearCombination<Self, C>` も IntoIterator を実装しているため、そのまま渡せる。`sum` は式の列、`from_terms` は係数付きの式の列を受け取る入口として区別する。既存の `LinearCombination::from_terms` は維持する。

| 提供型 | 現行 | 採用する入口 | 状態 |
| --- | --- | --- | --- |
| LocalState / LocalOperator | `from_linear_combination(LinearCombination<Self, C>)` | `from_terms(impl IntoIterator<Item = Term<Self, C>>)` | 採用済み・未実装 |

以下は採用後の使用イメージで、`L` は局所式の型、`a`・`b` はその式、`c1`・`c2` は係数とする。それぞれ独立した呼び出し例である。

```rust
let expr = L::from_terms([
    Term { coefficient: c1, operand: a },
    Term { coefficient: c2, operand: b },
]);
let expr = L::from_terms(combination); // LinearCombination<L, C>
let expr = L::from_terms(primitive_combination.map(L::primitive));
```

入力を順に消費し、各項をスケール節点として和へ取り込む。構築時に評価・同類項の集約・ゼロ係数の除去を行わず、入力順・重複項・ゼロ係数を保持する。空入力は空の和、単項入力はその項のスケール式となる。係数と operand は所有権を移し、数値演算・Clone・比較などの能力制約を追加しない。primitive を直接受け取る別の入口や From 実装の追加は今回の採用範囲に含めない。

**L6 の採用決定（2026-10-30）**：現状維持とし、`normalized_expectation` / `normalized_expectation_in`、`NormalizedExpectation` / `NormalizedExpectationIn` を維持する。非正規化の core の `expectation` との区別を保ち、L6 は改名の実施対象から外す。

`hs` と `_g` は特定の数学的概念の短縮形として限定して使う。`coordinates` を `coords` にするなど、一般の単語まで一律に略す規則にはしない。また QuantumSystem の状態側は `tensor_in` / `try_tensor_in` / `product_state` として実装済みである。L3 の演算子側も現行名を維持する。将来別の合成操作を追加する場合も、実装済み入口の入力型や意味を暗黙に変更しない。

### 5.3. 現行の状態構築と、残るテンソル・積の対応表

実装済みの状態構築と、未実装の命名案を分けて示す。未実装の改名では入力型・引数順序・戻り値・trait bound を変更しない。ここでの QuantumSystem は現行の構築入口を指す。第6節の移行後は同じメソッド名と契約を ExpressionFactory の inherent method に置き、SpinHalf などの型エイリアスから呼ぶ。

#### State / QuantumSystem の状態構築（実装済み・維持）

| 入口 | 現行の入力と役割 |
| --- | --- |
| `State::tensor_in` / `QuantumSystem::tensor_in` | `(&space, states)`。所有する配置済み状態の列を明示空間上で合成 |
| `State::try_tensor_in` / `QuantumSystem::try_tensor_in` | `(&space, states)`。Result の列を集約後、明示空間との整合性を検証 |
| `QuantumSystem::product_state` | `(&space, primitives)`。葉の順序に primitive を一つずつ配置して合成 |

二項も `State::tensor_in(&space, [left, right])` を使う。状態側に二項専用版・種類指定版・`_g` 版は設けない。指定空間のハンドルと grouping を保持し、入力の因子順序・被覆と部分系の種類の互換性を検証する。種類は HilbertSpace から取得し、Tensor 節点には保持しない。`State::unit()` と `State::zero(&space)` も維持する。

#### HilbertSpace（L4・採用済み・未実装）

**採用決定（2026-10-30）**：種類指定名の短縮と二項／列の区別を合わせ、次の六入口を一緒に改名する。

| 現行名 | 採用名 |
| --- | --- |
| `tensor_product` | `tensor` |
| `tensor_product_g` | `tensor_g` |
| `tensor_product_with_kind` | `tensor_with` |
| `tensor` | `tensor_all` |
| `tensor_g` | `tensor_all_g` |
| `tensor_with_kind` | `tensor_all_with` |

二項版は `lhs.tensor_with(kind, &rhs)`、列版は `HilbertSpace::tensor_all_with(kind, children)` とし、`_with` の先頭引数は現在と同じ `TensorKind` とする。入力型・引数順序・戻り値・能力制約、ordinary / graded の選択、空入力・単項入力・grouping の保持、重複 Index と種類の不整合の検証は維持する。HilbertSpace に `try_` 版を新設する判断は含めない。この表を State に適用しない。採用は命名方針の確定であり、この時点では実装を変更していない。

#### 通常の加算・乗算（採用済みの命名方針）

| 提供型 | 現行名 | 変更後の名前 |
| --- | --- | --- |
| State / Operator | inherent `add`（二項） | 廃止。標準の `+` または列の `sum` を使用 |
| 四つの式型 | `sum` / `try_sum` | `sum` / `try_sum`（維持） |
| LocalOperator / Operator | inherent `product`（二項） | 廃止。標準の `*` または列の `product` を使用 |
| LocalOperator / Operator | `product_all`（列） | `product` |
| Operator | `try_product_all` | `try_product` |
| primitive の Product trait | `Product::product`（二項規則） | trait 名・メソッド名とも維持 |

式型に二項専用の inherent method は追加しない。LocalOperator に `try_product` を新設する判断も含めない。既存の `core::ops::Add` / `Mul` の実装と、primitive の `Product<C>` / `Product::product` は維持する。

#### QuantumSystem の演算子積（現状維持・決定済み）

**採用決定（2026-10-30）**：三入口とも現行名を維持する。

| 現行名 | 決定後の名前 |
| --- | --- |
| `operator_product` | `operator_product`（維持） |
| `try_operator_product` | `try_operator_product`（維持） |
| `operator_product_at` | `operator_product_at`（維持） |

`state_sum` / `try_state_sum` と `operator_sum` / `try_operator_sum` は区別を維持する。現在の通常の積の対象は演算子だけだが、和と同様に対象を名前で明示する一貫性を優先する。式型の集約名を `product` にする N1 の採用決定は維持し、QuantumSystem の構築名では `operator_` を残す。`product_state` は明示空間上の積状態の構築として維持する。

#### 呼び出しの比較

| 現行 | 整理後 | 状態 |
| --- | --- | --- |
| `State::tensor_in(&h, [left, right])?` | 同じ | 実装済み・維持 |
| `Spin::product_state(&h, primitives)?` | 同じ | 実装済み・維持 |
| `h_a.tensor_product(&h_b)?` | `h_a.tensor(&h_b)?` | L4 の採用済み・未実装 |
| `Spin::try_operator_product(terms)?` | 同じ | L3 の現状維持・決定済み |
| `a.product(&b)?`（複合演算子） | `(&a * &b)?` | N1 の採用済み・未実装 |
| `a.add(&b)?`（複合状態・演算子の inherent method） | `(&a + &b)?` | N1 の採用済み・未実装 |
| `Operator::product_all(terms)?` | `Operator::product(terms)?` | N1 の採用済み・未実装 |

### 5.4. 基底・幾何の統合対応表

| 提供型 | 現行 | 整理後 | 状態 |
| --- | --- | --- | --- |
| 四つの式型 | `expand_in(&basis)` | `expand(&basis)` | N3・採用済み・未実装 |
| 四つの式型 | `coordinates_in(&basis)` | `coordinates(&basis)` | N3・採用済み・未実装 |
| Operator | `expand_in_space(&basis, &ambient)` | `expand_in(&basis, &ambient)` | N3・採用済み・未実装 |
| Operator | `coordinates_in_space(&basis, &ambient)` | `coordinates_in(&basis, &ambient)` | N3・採用済み・未実装 |
| Chain | `nn_site_pairs()` | 同じ | N5・維持 |
| Chain | `nn_site_pairs_with(&sites)` | 同じ | N5・維持 |
| Chain | `nn_site_pairs_on(&space)` | `nn_site_pairs_in(&space)` | N5・採用済み・未実装 |
| Chain | `nn_pairs_with(&sites)` | 同じ | N5・維持 |
| Chain | `nn_pairs_on(&space)` | `nn_pairs_in(&space)` | N5・採用済み・未実装 |
| Chain | `nn_pairs_by_pos_with(&sites)` | `nn_zip_with(&sites)` | N5・採用済み・未実装 |
| Chain | `nn_pairs_by_pos_on(&space)` | `nn_zip_in(&space)` | N5・採用済み・未実装 |

四つの式型は `LocalState`・`LocalOperator`・`State`・`Operator` を指す。局所式や State に明示空間版を追加する提案ではない。Chain の `_in` も、現行 `_on` と同じトップレベル分割を Site として使う。局所空間の展開順序や接続方法は改名で変更しない。

### 5.5. 内部関数についての結論

**採用決定（2026-10-30）**：第3節の四候補は、以下の推奨名を採用する。I1〜I4 はいずれも採用済み・未実装とし、公開範囲と処理内容は現行のままとする。

| ID | 現行名 | 採用名 | 理由・変更範囲 |
| --- | --- | --- | --- |
| I1 | `LocalOperator::evaluated` | `evaluate` | 状態側と同じ動詞にする。公開する必要はなく、`pub(crate)` のまま呼び出しを置換する。`evaluated_inner_product` は「既に評価済みの入力を要求する」という前提を表すので、この規則で機械的に変更しない。 |
| I2 | `StateExpr::same_node` | `shares_node_with` | 空間側と共有節点の同一性の語彙を揃える。`Arc::ptr_eq` を使う実装を維持し、式の数学的等価性を比較する名前にはしない。 |
| I3 | `basis::enumeration::leaves` | `collect_leaves` | `Vec` への集約と、公開 `HilbertSpace::leaves` の遅延列挙を区別する。既存走査を公開 iterator の利用へ置き換えるかどうかは改名とは別に扱う。 |
| I4 | `qmb-export::coordinates::owned_entries` | `into_flat_entries` | 所有権の移動と flat index 化の両方を表す。借用版の `flat_entries` と対応し、反復時にエラーを返す契約も維持できる。 |

### 5.6. 採用時の移行単位

採用済みの N1 では式型の二項 inherent `add` / `product` の廃止と、列の `product_all` → `product` / `try_product_all` → `try_product` を同時に適用する。primitive の `Product<C>` / `Product::product` と標準の Add / Mul は維持する。L3 の QuantumSystem の構築名は現状維持とし、N1 の式型の改名に伴って接頭辞を削除しない。N2・L2 の状態側の移行は完了している。L4 は HilbertSpace 側の六入口を一緒に移行する。採用済みの N3（L7 の基底側を含む）も、四つの式型の基本名と Operator の明示空間版を一緒に変更する。N5（L7 の Chain 側を含む）も採用済みとし、位置対応名と空間指定の接尾辞を合わせて移行する。旧名の互換 alias は設けず、v2 の破壊的変更として対応表に従って移行する。

通常の積では旧 inherent `product` が二項、新 inherent `product` が列を受け取るため、受け手と引数の形を見て移行する。旧二項の借用呼び出しは標準の `&a * &b`、旧 inherent `add` は標準の `&a + &b` などへ変更する。LocalOperator の二項に Identity 制約を追加しないようにする。`Product::product` や標準 trait のメソッド、および新しい列集約 `product` まで一括置換しない。

実装済みの状態構築名と契約は維持する。旧状態 API や状態側の `tensor_all*` の改名案は再導入しない。採用済みの L4 でも、HilbertSpace の旧 `tensor*`（列）を `tensor_all*` へ、旧 `tensor_product*`（二項）を `tensor*` へ対応づけ、引数の形を見て移行する。別の既存名の一部や置換済みの新名まで一括置換しないようにし、呼び出し・テスト・rustdoc・examples・設計文書を合わせて更新する。`expand_in` / `coordinates_in` も同様に、基底だけの呼び出しと明示空間付きの呼び出しを区別する。

N3 では Operator の旧 `expand_in(&basis)` / `coordinates_in(&basis)` と、新しい同名の明示空間版を併設できない。旧基本名を `expand` / `coordinates` へ、旧 `_in_space` 版を `_in` へ対応づけ、引数の役割・個数を見て移行する。明示空間の優先順位や基底展開・座標取得の意味は維持する。

採用済みの N4 では、第1節の対応表に従い、qmb-core の空間・座標・自由関数と qmb-export の局所次元メタデータを一緒に移行する。呼び出し・再公開・テスト・rustdoc・examples・設計文書も更新する。`Dimension::dim` は既存の必須メソッドのまま維持し、primitive の独自実装にメソッド改名を要求しない。単位空間の局所次元が `None` で全体次元が `1` であること、次元計算の overflow、非零座標数と次元の違い、演算子基底と行列出力の状態基底の次元の違いを維持する。

N5 では `nn_site_pairs_on` → `nn_site_pairs_in`、`nn_pairs_on` → `nn_pairs_in`、`nn_pairs_by_pos_with` → `nn_zip_with`、`nn_pairs_by_pos_on` → `nn_zip_in` を一緒に移行する。`nn_site_pairs`・`nn_site_pairs_with`・`nn_pairs_with` は維持する。対応する crate の呼び出し・テスト・rustdoc・examples・設計文書も更新し、Site 数の検証、保存されたトップレベル分割、列挙順序と zip の打ち切りを維持する。

L1 は LocalOperator・Operator の基本版と Operator の明示空間版を一緒に改名し、呼び出し・テスト・rustdoc・examples・設計文書も更新する。非正規化の `Tr(A†B)` と明示空間の優先を検証する。N6 のモデル名と L6 の正規化期待値名は現状維持のため、改名作業は行わない。

L5 は LocalState / LocalOperator の二入口を `from_terms` に変更し、項の列を受け取る入力へ一般化する。旧入口に渡していた LinearCombination は新入口にもそのまま渡せ、primitive の `map(Self::primitive)` も維持できる。呼び出し・テスト・rustdoc・examples・設計文書を更新し、配列・iterator・既存の LinearCombination、空入力・単項入力、重複項・ゼロ係数・入力順の保持、構築時に評価しないこと、および数値演算や Clone の制約が不要であることを検証する。

**互換性方針の採用決定（2026-10-30）**：今回の整理では既存 API との互換性維持を要求しない。名前を再利用しない改名も含め、旧名の deprecated alias や互換 wrapper は設けず、新名を正式な入口に揃える。対応する呼び出し・テスト・rustdoc・examples・設計文書を更新する。

I1〜I4 は内部の定義と呼び出しを採用名へ揃え、公開範囲・既存アルゴリズム・エラーの契約を維持する。I1 で LocalOperator の評価を新たに公開したり、I3 の集約処理を別の iterator 設計へ変更したりする作業は含めない。

採用後は該当 crate の既存テスト・rustdoc・examples を更新して実行し、第4節の契約を確認する。特に fallible 集約のエラー優先順位、graded の符号、明示空間の優先、基底座標と行列の区別を検証対象にする。本書の未実装の推奨名は、この文書更新では実装・コンパイル検証していない。実装済みの状態構築と未実装の命名案を区別し、解消済み項目を次回の改名対象に含めない。

## 6. QuantumSystem・ExpressionFactory とモデル構築入口の分離

**採用決定（2026-10-30、名称と役割を更新）**：型情報を関連付ける公開 trait は `QuantumSystem`、共通の式木を構築する公開 struct は `ExpressionFactory<S, O, C>` とする。構築メソッドは struct の inherent method に置き、具体的な primitive と係数型を選ぶ名前は型エイリアスにする。モデル固有構築は各 primitive crate の別の入口へ移す。以下は採用済み・未実装の設計であり、この文書更新ではコードを変更しない。

### 6.1. 型情報と共通の構築入口

qmb-core の `QuantumSystem` は trait として名前を維持し、状態 primitive・演算子 primitive・係数型の三つの関連型だけを定義する。構築メソッドや数値能力の制約は持たせない。構築用 struct は `ExpressionFactory<S, O, C>` とし、qmb-core から両方を公開する。

```rust
pub trait QuantumSystem: Sized {
    type State;
    type Operator;
    type Coeff;
}

pub struct ExpressionFactory<S, O, C>(
    PhantomData<fn() -> (S, O, C)>,
);

impl<S, O, C> QuantumSystem for ExpressionFactory<S, O, C> {
    type State = S;
    type Operator = O;
    type Coeff = C;
}
```

ExpressionFactory は値を保持せず、インスタンス生成を必要としない型選択・構築の入口とする。QuantumSystem は全ての ExpressionFactory に一括実装し、通常の利用者に trait 実装を要求しない。

既存の `state` / `operator` / `tensor_in` / `try_tensor_in` / `product_state` / `operator_product_at` / `state_sum` / `try_state_sum` / `operator_sum` / `try_operator_sum` / `operator_product` / `try_operator_product` は ExpressionFactory の inherent method へ移し、名前と構築契約を維持する。L3 で採用した接頭辞の維持も変わらない。空間・Index の登録や保持、構築時の評価・数値能力・Clone の一律要求は追加しない。

具体的な入口は次の型エイリアスとする。

```rust
pub type SpinHalf<C> = ExpressionFactory<SpinHalfState, SpinHalfOperator, C>;
pub type SpinOne<C> = ExpressionFactory<SpinOneState, SpinOneOperator, C>;
pub type Fermion<C> = ExpressionFactory<FermionState, FermionOperator, C>;
```

これにより `SpinHalf::<C>::state(...)` などを QuantumSystem trait の import なしで呼べる。独自の primitive を利用する場合も、利用側で `type Custom = ExpressionFactory<Ket, Op, Coeff>` と定義できる。同じ primitive と係数型からは同じ共通式型を構築し、システム名を式の同一性へ追加しない。

### 6.2. モデル固有構築の入口

各 primitive crate の公開 `models` module に、以下の構築専用 struct を置く。前の `ChainModels<Q>` 案の型選択方法を引き継ぎ、`Q` には係数型を選択済みのシステム型を渡す。これらの struct は空間や結合定数などの値を保持せず、関連関数を提供する。

| 定義する crate | 採用する入口 | 移動する既存のメソッド |
| --- | --- | --- |
| qmb-spin-half | `models::SpinHalfModels<Q>` | `xyz` / `heisenberg` / `neel` |
| qmb-spin-one | `models::SpinOneModels<Q>` | `xyz` / `heisenberg` / `blbq` |
| qmb-fermion | `models::FermionModels<Q>` | `hubbard` / `symmetric_hubbard` |

各型で対応するシステム型への inherent impl を定義する。例えば qmb-spin-half では `impl<C> SpinHalfModels<SpinHalf<C>>` にモデル構築メソッドを置く。型自体を各 primitive crate が所有するため、qmb-core の型エイリアスへ inherent method を追加する制約を回避できる。モデル固有メソッドを qmb-core へ移したり、呼び出すための拡張 trait を要求したりしない。

`ChainModels` という名前は採用しない。現在のメソッドが Chain を受け取ることは各シグネチャで表し、構築入口全体を Chain に限定しない。将来、別の幾何に対応する構築関数も同じ入口へ追加できる。今回の変更で既存モデルの対応幾何自体を増やすことはしない。

以下は採用後の使用イメージである。

```rust
use num_complex::Complex;
use qmb_spin_half::{
    SpinHalf, SpinHalfState as Ket, SpinHalfOperator as Op,
    models::SpinHalfModels,
};

type Spin = SpinHalf<Complex<f64>>;
type Models = SpinHalfModels<Spin>;

let up = Spin::state(Ket::Up);
let sx = Spin::operator(Op::Sx);
let h = Models::heisenberg(&chain, &sites, j)?;
let psi = Models::neel(&chain, &sites, Ket::Up)?;
```

係数型は Spin で一度選び、モデル構築でも同じ選択を使う。モデル固有の能力制約は各メソッドへ限定する。N6 はモデルのメソッド名の維持として引き継ぎ、`symmetric_hubbard` の名前は変えず、所属する入口を `FermionModels` へ移す。

### 6.3. 型注釈とジェネリック API

**採用決定（2026-10-30）**：`LocalStateOf<Q>` / `LocalOperatorOf<Q>` / `StateOf<Q, I>` / `OperatorOf<Q, I>` は公開型エイリアスとして維持する。型情報を取得する trait には QuantumSystem を使う。

システム単位で抽象化する公開のジェネリック API は `Q: QuantumSystem` を使う。構築を行う箇所は、型が具体化された式型の入口、または `ExpressionFactory<Q::State, Q::Operator, Q::Coeff>` の inherent method を使う。旧 trait の構築メソッドに依存する `Q::state(...)` などは書き換える。個々の式型を処理する低レベルの関数では、必要な `S, O, C` を直接型引数に取れる。

通常の `Spin::state(...)` には trait の import を要求しない。型情報を扱うジェネリック API の定義側だけが、必要に応じて QuantumSystem を利用する。型情報専用の別名の補助 trait や、例示した `make_state` 関数を新たな公開 API として追加することはしない。

### 6.4. 移行範囲

QuantumSystem の構築メソッドの移動と関連型への役割の限定、ExpressionFactory の追加、具体的なシステム型の型エイリアス化、共通構築の呼び出しと import、モデル構築の定義と呼び出し、テスト・rustdoc・examples・設計文書を合わせて移行する。既存の計算結果、明示空間の選択と共有、入力エラーの優先、係数精度と能力制約は維持する。旧システム型のモデル構築メソッドや QuantumSystem の構築メソッドを互換 wrapper として残さない。

旧 QuantumSystem trait の実装は、基本構築に関しては ExpressionFactory の型エイリアスへ置き換えられる。ジェネリック API では関連型を維持して構築処理を書き換え、既存の `*Of` 型注釈が新しい型エイリアスでも利用できることを検証する。

根拠: [現在の QuantumSystem と型エイリアス](../../../crates/qmb-core/src/system.rs#L43)、[SpinHalf の型選択](../../../crates/qmb-spin-half/src/lib.rs#L123)、[SpinHalf のモデル構築](../../../crates/qmb-spin-half/src/models.rs#L34)、[SpinOne のモデル構築](../../../crates/qmb-spin-one/src/models.rs#L31)、[Fermion のモデル構築](../../../crates/qmb-fermion/src/models.rs#L29)、[ジェネリック構築と型注釈の既存テスト](../../../crates/qmb-core/tests/systems.rs#L25)。

## 関連文書

- [明示空間による状態構築](../../design/composite-expression.md#121-構築入口)：状態側のテンソル構築に関する実装済みの契約。
- [API の使いやすさに関する残る検討事項](../api-usability/README.md)：補助 API、型推論、利用開始の導線など。
- [局所式の設計](../../design/primitive-expression.md)：構築・評価・数値能力と正規化の境界。
- [複合式の設計](../../design/composite-expression.md)：集約・束縛・ordinary / graded の契約。
- [Hilbert 空間の設計](../../design/hilbert-space.md)：因子・順序・合成規則。
- [幾何の設計](../../design/geometry.md)：Site 対と葉対、モデル構築。
