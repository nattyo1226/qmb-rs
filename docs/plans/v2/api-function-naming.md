# 関数名・公開 API 名の整理候補

作成日: 2026-10-30
更新日: 2026-10-30

現在の実装を読み、同種の操作で名前の規則が揃っていない箇所と、公開 API として長くなっている名前を記録する。実装名・入力・動作はソースで確認した事実であり、使いにくさの評価・優先度・新しい名前は検討案である。第1〜3節は現行 API に対する候補と解消済み項目、第5節は採用状況と推奨案を記録する。N1 の通常の加算・乗算については、二項を `add` / `mul`、列の集約を `sum` / `product` とする方針を採用済みとする。状態構築の変更は実装済みで、複合式の設計文書に記録する。それ以外は推奨段階であり、実装時期は未決である。

確認対象は `crates/` 以下の `qmb-core`・`qmb-basic`・`qmb-geometry`・`qmb-spin-half`・`qmb-spin-one`・`qmb-fermion`・`qmb-normalize`・`qmb-export` の関数・メソッド・公開 trait と、`examples/src/bin/` の呼び出しである。非公開 module 内でも公開型の inherent method や再公開された trait は公開 API に含める。内部関数は後段に分ける。

優先度「高」は日常的な構築・計算で長さや名前の衝突が目立つ候補、「中」は API 間の規則を揃える候補、「低」は現行名を残す判断も有力な候補とする。識別子の長さは `_` を含み、型名・module path は含めない。文字数だけを改名の基準にはしない。

明示空間への移行後の実装を再確認した。以下の ID は従来の議論を追えるよう維持し、解消済みの項目は改名の実施対象から外す。

| 対象 | 現在の状態 |
| --- | --- |
| N1 の通常の加算・乗算 | `add` / `mul` と `sum` / `product` の方針は採用済み、改名は未実装 |
| N2 の状態側の種類指定・修飾順序 | 明示空間から種類を取得する構築へ移行して解消済み |
| L2 の長い状態構築名 | `tensor_in` / `try_tensor_in` / `product_state` への移行で解消済み |
| L4 の種類指定付き tensor 名 | State 側は撤去済み。HilbertSpace 側だけを引き続き検討する |
| その他の候補 | 現行実装に残っているため、従来の採用決定・推奨案を維持する |

現在の状態構築契約は[複合式の設計](../../design/composite-expression.md#121-構築入口)を参照する。所属空間の必須化と、その空間からの種類取得を前提に API 整理を続ける。

## 1. 名前の規則が揃っていない箇所

### N1. 二項演算と列の集約の区別（中）

| 操作 | 二項版 | 列の集約版 |
| --- | --- | --- |
| 状態の和 | `State::add` | `State::sum` / `try_sum` |
| 演算子の積 | `LocalOperator::product`、`Operator::product` | `product_all`、複合式には `try_product_all` もある |
| 空間のテンソル積 | `HilbertSpace::tensor_product` | `HilbertSpace::tensor` |
| 状態のテンソル積 | 専用の二項版はなく、2要素の列で構築 | `State::tensor_in` / `try_tensor_in`、QuantumSystem にも同名の入口 |
| QuantumSystem の演算子積 | 二項専用の入口はない | `operator_product` / `try_operator_product` |

列を表す印が演算子の積では `_all`、空間のテンソル積では無印になり、QuantumSystem の演算子積では `_all` が消える。状態のテンソル構築は明示空間を受け取る列の入口へ統一済みで、`_in` は所属空間を表す。各シグネチャは明確だが、別の型で同じ操作を探す際に規則を推測しにくい。

**採用決定（2026-10-30）**：通常の加算・乗算は、二項演算の動詞と列の集約名で区別する。

| 操作 | 二項演算 | 列の集約 |
| --- | --- | --- |
| 加算 | `add` | `sum` |
| 乗算 | `mul` | `product` |

`LocalOperator` / `Operator` の二項 `product` は `mul`、列の `product_all` は `product`、Operator の `try_product_all` は `try_product` に変更する方針とする。列の集約には空入力・単項入力も含む。`sum` / `sum_all`、`product` / `product_all` という `_all` による区別案は通常の加算・乗算には採用しない。

ここでの名前付き二項メソッドは式型の inherent method を指す。`core::ops::Add::add` / `Mul::mul` による `+` / `*` の実装も存在するが、名前を揃えることで借用・所有権や戻り値の契約を変更するものではない。`Iterator::sum` / `product` と語彙は揃うが、新たな iterator 集約 trait の実装を決めるものでもない。

primitive の二項規則 `Product::product` もメソッド名を `mul` に揃える対象とする。独自の `Product` trait 自体の名前、および標準の `core::ops::Mul` との呼び分けは別途整理する。局所規則が線型結合を返す契約は維持する。

HilbertSpace のテンソル積の `tensor` / `tensor_all` は通常の乗算とは別の推奨案として残し、N1 の採用決定には含めない。空間側の現行 `tensor` は列を受け取るため、この案を採用する場合も二項・列を同時に移行する必要がある。状態側の `tensor_in` / `try_tensor_in` は実装済みの名前を維持し、`tensor_all*` への改名は行わない。

根拠: [State の構築](../../../crates/qmb-core/src/expr/composite.rs#L262)、[Operator の積](../../../crates/qmb-core/src/expr/composite.rs#L514)、[局所演算子の積](../../../crates/qmb-core/src/expr/local.rs#L233)、[空間のテンソル積](../../../crates/qmb-core/src/space.rs#L181)、[QuantumSystem の積](../../../crates/qmb-core/src/system.rs#L150)。

### N2. 状態側の種類指定・修飾順序（解消済み）

状態側の `_g` / `_with_kind` と、配置を伴う `state_tensor*` の入口は撤去した。配置済み状態の合成は `State::tensor_in(&space, states)` / `try_tensor_in`、primitive の葉順の配置は `QuantumSystem::product_state(&space, primitives)` を使う。QuantumSystem にも `tensor_in` / `try_tensor_in` を提供している。

テンソル積の種類の情報源は `HilbertSpace` に集約した。State は所属空間の共有ハンドルを持ち、状態の Tensor 節点は子の列だけを保持する。遅延作用の StateAction は構築時に選択した対象空間を保持する。このため、状態 API の修飾順序を変えたり `_at_g` / `_with` を新設したりする必要はない。

HilbertSpace の `tensor_g` / `tensor_product_g` / `tensor_with_kind` は引き続き存在する。空間側には配置を表す `_at` はなく、今回問題にした修飾順序の揺れは残っていない。`_g` は維持する推奨とし、空間側の名前の長さ・二項／列の区別は N1・L4 で扱う。

根拠: [QuantumSystem の構築入口](../../../crates/qmb-core/src/system.rs#L62)、[State の構築](../../../crates/qmb-core/src/expr/composite.rs#L310)、[状態式と作用節点](../../../crates/qmb-core/src/expr/composite.rs#L95)、[空間の構築](../../../crates/qmb-core/src/space.rs#L181)。

### N3. `_in` が基底と空間の両方を表す（中）

| 現行名 | `_in` が指す入力 |
| --- | --- |
| `LocalState` / `LocalOperator` / `State` / `Operator` の `expand_in(&basis)`、`coordinates_in(&basis)` | 基底 |
| `State` / `QuantumSystem` の `tensor_in(&space, states)` / `try_tensor_in(&space, states)` | 明示する所属空間 |
| `State::act_in` / `expectation_in`、`Operator::trace_in` / `hilbert_schmidt_inner_product_in`、`normalized_expectation_in` | 明示する作用空間 |
| `Operator::expand_in_space(&basis, &ambient)` / `coordinates_in_space(&basis, &ambient)` | 基底に加えて明示する作用空間 |
| `MatrixArgs::in_space(&ambient)` | 行列出力で優先する作用空間 |

新しい状態構築でも `_in` は空間を表すため、基底処理との用語の違いが残っている。基底処理は基本名の時点で `_in` を使い、明示空間を加えると `_in_space` になる。特に `coordinates_in_space` は20文字で、他の明示空間版から名前を推測しにくい。

候補は基底を受け取る基本名を `expand(&basis)` / `coordinates(&basis)`、明示空間版を `expand_in(&basis, &ambient)` / `coordinates_in(&basis, &ambient)` に揃えること。`Basis::expand` とも語彙が揃う。実装済みの状態構築名 `tensor_in` / `try_tensor_in` と、葉順の配置を表す `product_state` は維持し、N3 では改名しない。同名メソッドの引数数が変わるため、局所式・複合式をまとめて移行する必要がある。`MatrixArgs::in_space` は builder の役割が明瞭であり、計算メソッドと同じ形へ無理に短縮しなくてもよい。

根拠: [局所式の基底処理](../../../crates/qmb-core/src/basis/expansion/local.rs#L93)、[複合演算子の基底処理](../../../crates/qmb-core/src/basis/expansion/composite.rs#L241)、[作用・期待値](../../../crates/qmb-core/src/evaluate/inner_product.rs#L327)、[trace・内積](../../../crates/qmb-core/src/evaluate/trace.rs#L60)、[正規化期待値](../../../crates/qmb-normalize/src/expectation.rs#L23)、[MatrixArgs](../../../crates/qmb-export/src/args.rs#L24)。

### N4. 次元を表す `dim` / `dimension` / `basis_len`（中）

primitive は `Dimension::dim()`、空間・座標は `dimension()` / `local_dimension()` を使う。また `state_basis_len::<S>()` / `operator_basis_len::<O>()` は基底を受け取らず、内部で primitive の `dim()?.get()` を返している。

同じ次元の問い合わせに略語・完全形・長さという三つの語彙がある。候補は `dimension` に統一し、二つの自由関数も `state_dimension` / `operator_dimension` にするか、primitive の次元問い合わせへ集約すること。primitive の戻り値の `NonZeroUsize` と、空間・自由関数の `usize` の違いは維持する。

`Coordinates::len()` は保存された非零座標の数、`dimension()` は座標空間全体の次元であり、こちらの区別は必要である。primitive の次元も enum の variant 数ではなく張る線型空間の次元なので、単なる `len` にはしない。

根拠: [Dimension](../../../crates/qmb-core/src/ops.rs#L81)、[HilbertSpace](../../../crates/qmb-core/src/space.rs#L149)、[基底長の自由関数](../../../crates/qmb-core/src/basis/enumeration.rs#L162)、[Coordinates](../../../crates/qmb-core/src/basis/coordinates.rs#L66)。

### N5. Chain の `pairs` が Site 対か葉対か分かりにくい（中）

`nn_site_pairs` は Site 番号対、`nn_site_pairs_with` は与えた要素への参照対、`nn_pairs_with` は最近接 Site 間の全葉の直積を返す。`nn_pairs_by_pos_with` は葉を位置順に zip し、短い側で止まる。後二者だけ対象の `leaf` が名前に現れない。`_on` 版にも同じ違いがある。

候補は Site の入口を維持し、葉の入口を `nn_leaf_pairs_with` / `nn_leaf_pairs_on`、位置対応を `nn_leaf_zip_with` / `nn_leaf_zip_on` のようにすること。全組み合わせ版は少し長くなるが、位置対応版は `nn_pairs_by_pos_with`（20文字）から短くできる。`by_pos` を `zip` とするなら、短い側で止まる契約も説明しやすい。

`_with` は外部の Site 順の列、`_on` は一つの空間の保存されたトップレベル分割を使うという意味の違いがある。この区別や、列長を iterator 返却前に検証する契約は維持する。全組み合わせと zip を単一の `pairs` へまとめない。

根拠: [Chain の公開入口](../../../crates/qmb-geometry/src/geometry/chain.rs#L39)、[葉対の実装](../../../crates/qmb-geometry/src/connection.rs#L24)。

### N6. モデル名の修飾語の位置（低）

`Fermion::hubbard` に対する派生名は `Fermion::symmetric_hubbard` で、操作名の後に修飾語を置く `_g` / `_in` 系と語順が異なる。`hubbard_symmetric` にすれば名前順の一覧で同じモデル族が並ぶが、文字数は変わらない。一般的なモデルの呼び方を優先して現行名を残す判断もある。

根拠: [Hubbard の二つの入口](../../../crates/qmb-fermion/src/models.rs#L54)。

## 2. 長い公開名の短縮候補

N1〜N5 と同じ関数を含む。以下の候補名は相互に独立して採用できるとは限らず、二項／列の区別や修飾順序を決めてから確定する。

| ID・優先度 | 現行の関数名と長さ | 使いにくさと候補 |
| --- | --- | --- |
| L1・高 | `LocalOperator::hilbert_schmidt_inner_product`、`Operator::hilbert_schmidt_inner_product`（29）、`Operator::hilbert_schmidt_inner_product_in`（32） | 数値計算の呼び出しで長い数式名を毎回綴る。`hs_inner_product` / `hs_inner_product_in`（16 / 19）を候補とする。`hs` の意味を rustdoc で示し、非正規化の `Tr(A†B)` という契約を残す。単に `dot` とすると左側の随伴が読み取りにくい。 |
| L2・解消済み | 旧 QuantumSystem の `state_tensor*` の種類指定付き構築名 | 長い旧名は撤去済み。現在は `tensor_in`（9）、`try_tensor_in`（13）、`product_state`（13）を提供し、所属空間を明示する。短縮用の別名や `tensor_all*` の追加は不要。 |
| L3・中 | `QuantumSystem::operator_product_at`（19）、`try_operator_product`（20） | 基本構築では繰り返し現れる。現在の積構築は演算子だけなので、採用済みの N1 に合わせて `product_at` / `try_product`（10 / 11）などへ寄せる余地がある。`state_sum` と `operator_sum` は両方存在するため、全入口から型の接頭辞を取り除くことはできない。 |
| L4・高（空間側のみ） | `HilbertSpace::tensor_product_with_kind`（24） | State の種類指定版は撤去済み。空間側は、N1 の空間の二項／列の区別案に合わせ、二項を `tensor_with(kind, &rhs)`、列を `tensor_all_with(kind, children)` とする候補が残る。入力型 `TensorKind` と文書で `_with` の意味を明示する。 |
| L5・中 | `LocalState::from_linear_combination`、`LocalOperator::from_linear_combination`（23） | 変換元の型名をすべて綴る。`from_combination`（16）、または `From<LinearCombination<Self, C>>` を実装して `.into()` で変換する案。現在の入力は primitive ではなく局所式の線型結合であり、primitive の規則結果には `map(Self::primitive)` が必要な点を隠さない。`from_terms` に寄せて iterator を受け取るなら、単なる改名とは別に入力型も検討する。 |
| L6・低 | `NormalizedExpectation::normalized_expectation`（22）、`NormalizedExpectationIn::normalized_expectation_in`（25） | `mean` / `mean_in` のような短い名前を候補にできるが、計算の意味を明確にする利点が現行名にはある。core の `expectation` は非正規化なので、正規化版だけを `expectation` へ改名すると既存の inherent method と衝突する。単独の改名は避け、必要なら期待値の名前全体を再検討する。 |
| L7・中 | `Operator::coordinates_in_space`（20）、`Chain::nn_pairs_by_pos_with`（20） | N3 の基底／空間の語尾と、N5 の葉の zip の整理で短縮できる。名前だけを切り詰めず、それぞれのファミリを同時に揃える。 |

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
| `LocalOperator::evaluated`（`pub(crate)`）と公開の `LocalState::evaluate` / `State::evaluate` | いずれも元を変更せず評価済みの式を返すが、演算子側だけ過去分詞になっている。内部名も `evaluate` に揃える案。公開 API の追加は別の判断であり、現状の演算子に公開 `evaluate` があるとは扱わない。[局所評価](../../../crates/qmb-core/src/evaluate/local.rs#L138)、[状態評価](../../../crates/qmb-core/src/evaluate/composite.rs#L41)を参照。 |
| `StateExpr::same_node`（`pub(crate)`）と公開の `HilbertSpace::shares_node_with` | どちらも `Arc::ptr_eq` で共有節点の同一性を調べる。`shares_node_with` に揃えると数学的な同値判定と区別しやすい。[式](../../../crates/qmb-core/src/expr/composite.rs#L182)、[空間](../../../crates/qmb-core/src/space.rs#L171)を参照。 |
| `basis::enumeration::leaves`（`pub(super)`）と公開の `HilbertSpace::leaves` | 前者は葉をすべて `Vec` へ集め、後者は遅延 iterator を返す。内部の前者を `collect_leaves` にすると確保・集約を伴うことが伝わる。[集約](../../../crates/qmb-core/src/basis/enumeration.rs#L15)、[iterator](../../../crates/qmb-core/src/space.rs#L138)を参照。 |
| `qmb-export` の `coordinates::owned_entries`（`pub(crate)`） | 所有権の移動だけでなく、反復時に多桁の座標を flat index に変換し、検証・overflow 検出を行う。`into_flat_entries` とすれば、core の借用版 `Coordinates::flat_entries` と処理内容を対応づけられる。[変換](../../../crates/qmb-export/src/coordinates.rs#L21)、[借用版](../../../crates/qmb-core/src/basis/coordinates.rs#L100)を参照。 |

## 4. 整理時に維持する区別と検証

- **`try_` は単なる Result の印ではない。** `State::tensor_in` / `Operator::product_all` も `Result` を返すが、`try_tensor_in` / `try_product_all` は入力が `Result` の列である。最初の入力エラーで消費を止め、全入力成功後に構成を検証するという一貫した契約がある。名前を揃えるためにすべての fallible 操作へ `try_` を付けない。[契約](../../../crates/qmb-core/src/system.rs#L38)と[検証対象の既存テスト](../../../crates/qmb-core/tests/fallible_construction.rs)を参照。
- **名前が違っても、受け手や入力が違えば必要な区別がある。** `QuantumSystem::state` / `operator` は primitive と係数型を選択済みの構築入口、`LocalState::primitive` / `LocalOperator::primitive` は式型の入口である。`local.at(&factor)` と `State::place(&local, &factor)` / `Operator::place(&local, &factor)` はメソッドと関連関数の入口である。これらを長い `local_*` 名へ機械的に揃える優先度は低い。
- **短縮で数学的契約を変えない。** 左側の共役・随伴、積の右側が先に作用する順序、ordinary / graded の選択、正規化の有無、明示空間が古い束縛に優先する規則を維持する。局所式の恒等元と未束縛の複合恒等元、列集約の空入力・単項の扱いも維持する。
- **名前と能力を別に扱う。** primitive・Index の Clone 要求、数値能力、係数精度、no_std の境界を改名に伴って追加・変更しない。`normalized_expectation` の短縮のために `normalized().expectation()` へ置き換えると、平方根など現在不要な能力を要求するので同等の移行にはならない。
- **採用時に名前の衝突と利用例を確認する。** 特に N1・N3 は既存名の引数形を変えるため、同じ型に旧名の互換 alias をそのまま残せるとは限らない。公開 trait の変更は独自代数の実装にも影響する。対応する crate のテスト・rustdoc・examples と、設計文書の API 表を合わせて更新・検証する。本書の調査ではソースと呼び出し例を確認し、実行テストは行っていない。

## 5. 各候補の推奨案

検討日: 2026-10-30

判断基準は、型をまたいで操作名から入力の形を推測できること、数学的な意味を失わないこと、頻出する呼び出しを短くすることの三点とする。既存名を残すだけで十分な項目も含めて結論を示す。第1〜3節の代替案のうち、以下で選ばなかったものは今回の推奨には含めない。根拠となる定義箇所は各候補のリンクを参照する。

### 5.1. 命名規則についての結論

| 対象 | 推奨 | 判断理由・選ばなかった案 |
| --- | --- | --- |
| N1 | **採用済み**。通常の二項演算は `add` / `mul`、列の集約は `sum` / `product` とする。HilbertSpace の `tensor` / `tensor_all` は推奨段階、状態の `tensor_in` は実装済みとして維持する。 | 標準の二項演算と iterator 集約の語彙に揃い、加算・乗算の区別方法が対称になる。通常の積の `product_all` は不要となる。`_all` を通常の和・積にも共通適用する案は採らない。 |
| N2 | **解消済み・現行維持**。状態側は明示空間を受け取る構築へ移行済み。HilbertSpace 側の `_g` は維持する。 | 種類の指定・保持を空間へ集約し、状態側の `_at_g` / `_with` の改名案は不要になった。状態の Tensor 節点に種類を戻さず、遅延作用の対象空間は維持する。 |
| N3 | **改名**。基底を指定する基本名を `expand` / `coordinates`、明示空間版を `expand_in` / `coordinates_in` とする。実装済みの `tensor_in` / `try_tensor_in` / `product_state` と `act_in`・`expectation_in`・`trace_in` と `MatrixArgs::in_space` は維持する。 | 式の構築・計算メソッドで `_in` の対象を所属・作用空間へ揃えられる。基底は必須引数なので無印でも情報を失わない。全計算名を `_in_space` に長くする案や、配置の `_at` を流用する案は採らない。builder は計算ではないため `in_space` を残す。 |
| N4 | **`dimension` に統一**。`Dimension::dim` → `dimension`、`state_basis_len` → `state_dimension`、`operator_basis_len` → `operator_dimension` とする。 | 基底を受け取らない自由関数の意味も名前に反映できる。二つの自由関数は `usize` を返す便宜を残すため削除しない。既存の `dimension_with`・`local_dimension`・`Coordinates::len` は維持し、戻り値・overflow の契約も変更しない。 |
| N5 | **葉を明示して改名**。全組み合わせは `nn_leaf_pairs_with` / `nn_leaf_pairs_on`、位置対応は `nn_leaf_zip_with` / `nn_leaf_zip_on` とする。Site の三入口は維持する。 | `leaf` で返す対象、`zip` で短い側までの位置対応を示す。`nn_` の削除は最近接という情報を失い、`_with` と `_on` の統合は入力の解釈を隠すため採らない。全組み合わせ版が少し長くなることは意味の明確化のために受け入れる。 |
| N6 | **現状維持**。`hubbard` / `symmetric_hubbard` を残す。 | モデル名の形容詞と API の入力条件を示す接尾辞は役割が異なる。自然なモデル名を変えて得られるのは一覧上の並びだけであり、短縮や意味の改善がない。関連する二関数の rustdoc から相互に参照できれば十分である。 |

### 5.2. 長い公開名についての結論

| 対象 | 推奨する公開名 | 判断理由・移行範囲 |
| --- | --- | --- |
| L1 | **改名**。`hs_inner_product` / `hs_inner_product_in` | Hilbert–Schmidt の部分だけを `hs` とし、内積であることは省略しない。局所版・複合版を同時に変える。`inner_product` だけにする案は内積の種類を隠し、`dot` / `hs_inner` はさらに語彙を減らす必要がないため採らない。rustdoc の冒頭に正式名と `Tr(A†B)`、非正規化であることを残す。 |
| L2 | **解消済み・現行維持**。`tensor_in` / `try_tensor_in` / `product_state` | 型と種類指定を重ねた旧名は撤去済み。第5.3節の実装済み入口を使い、`tensor_all_at_with` などの旧推奨名は追加しない。 |
| L3 | **改名を推奨**。`product_at`、`try_product`。通常の `operator_product` も `product` とする。 | QuantumSystem の積は演算子だけを受け取るので `operator_` を省ける。式型の列集約と同名になり、配置も `_at` で表せる。`op_product` のような新しい略語は導入しない。状態・演算子の双方がある `state_sum` / `operator_sum` とそれぞれの `try_` 版は維持する。 |
| L4 | **HilbertSpace 側のみ改名を推奨**。二項は `tensor_with`、列は `tensor_all_with` | 空間の二項／列の区別と合わせて短縮する。状態側の種類指定 API は存在せず、状態に `tensor_with` / `try_tensor_all_with` を追加する提案ではない。空間の種類の暗黙推論は導入しない。 |
| L5 | **改名**。両局所式型で `from_combination` | `LinearCombination<Self, C>` という既存の入力型と型推論を保ち、名前だけを短縮できる。主入口を `From` / `.into()` に置き換える案は変換先の推論に依存し、`from_terms` は入力を iterator へ広げる別の API 変更になるため今回は採らない。primitive の線型結合を受け取る別 overload も追加しない。 |
| L6 | **現状維持**。`normalized_expectation` / `normalized_expectation_in` と二つの trait 名を残す。 | 非正規化の `expectation` と式上で明確に区別できることを優先する。`mean` は短いが正規化規約が名前に表れず、`expectation` への単独改名は既存メソッドと競合する。名前の短縮のために平方根を要する正規化処理へ置き換えない。 |
| L7 | **N3・N5 に統合して改名**。`coordinates_in_space` → `coordinates_in`、`nn_pairs_by_pos_with` → `nn_leaf_zip_with` | それぞれのファミリ全体に対する決定を適用し、この二関数だけ別の略称を追加しない。対応する基本版・`_on` 版も同時に移行する。 |

`hs` と `_g` は特定の数学的概念の短縮形として限定して使う。`coordinates` を `coords` にするなど、一般の単語まで一律に略す規則にはしない。また QuantumSystem の状態側は `tensor_in` / `try_tensor_in` / `product_state` として実装済みである。L3 の演算子側の接頭辞省略は引き続き推奨段階とする。将来別の合成操作を追加する場合も、実装済み入口の入力型や意味を暗黙に変更しない。

### 5.3. 現行の状態構築と、残るテンソル・積の対応表

実装済みの状態構築と、未実装の命名案を分けて示す。未実装の改名では入力型・引数順序・戻り値・trait bound を変更しない。

#### State / QuantumSystem の状態構築（実装済み・維持）

| 入口 | 現行の入力と役割 |
| --- | --- |
| `State::tensor_in` / `QuantumSystem::tensor_in` | `(&space, states)`。所有する配置済み状態の列を明示空間上で合成 |
| `State::try_tensor_in` / `QuantumSystem::try_tensor_in` | `(&space, states)`。Result の列を集約後、明示空間との整合性を検証 |
| `QuantumSystem::product_state` | `(&space, primitives)`。葉の順序に primitive を一つずつ配置して合成 |

二項も `State::tensor_in(&space, [left, right])` を使う。状態側に二項専用版・種類指定版・`_g` 版は設けない。指定空間のハンドルと grouping を保持し、入力の因子順序・被覆と部分系の種類の互換性を検証する。種類は HilbertSpace から取得し、Tensor 節点には保持しない。`State::unit()` と `State::zero(&space)` も維持する。

#### HilbertSpace（改名案・未採用）

| 現行名 | 推奨名 |
| --- | --- |
| `tensor_product` | `tensor` |
| `tensor_product_g` | `tensor_g` |
| `tensor_product_with_kind` | `tensor_with` |
| `tensor` | `tensor_all` |
| `tensor_g` | `tensor_all_g` |
| `tensor_with_kind` | `tensor_all_with` |

空間側の `_with` の先頭引数は現在と同じ `TensorKind` とする。HilbertSpace に `try_` 版を新設する提案は含まない。この表を State に適用しない。

#### 通常の加算・乗算（採用済みの命名方針）

| 提供型 | 現行名 | 変更後の名前 |
| --- | --- | --- |
| State / Operator | `add` | `add`（維持） |
| 四つの式型 | `sum` / `try_sum` | `sum` / `try_sum`（維持） |
| LocalOperator / Operator | `product`（二項） | `mul` |
| LocalOperator / Operator | `product_all`（列） | `product` |
| Operator | `try_product_all` | `try_product` |
| primitive の Product trait | `product`（二項） | `mul`（trait 自体の名前は別途検討） |

局所式に新しい inherent `add` や `try_product` を追加する判断は含めない。既存の `core::ops::Add` / `Mul` の実装は維持する。

#### QuantumSystem の演算子積（改名案・未採用）

| 現行名 | 推奨名 |
| --- | --- |
| `operator_product` | `product` |
| `try_operator_product` | `try_product` |
| `operator_product_at` | `product_at` |

`state_sum` / `try_state_sum` と `operator_sum` / `try_operator_sum` は区別を維持する。`product_state` は積状態の構築、提案する `product` は演算子の積の集約であり、型名の接頭辞を全て削除する規則にはしない。

#### 呼び出しの比較

| 現行 | 整理後 | 状態 |
| --- | --- | --- |
| `State::tensor_in(&h, [left, right])?` | 同じ | 実装済み・維持 |
| `Spin::product_state(&h, primitives)?` | 同じ | 実装済み・維持 |
| `h_a.tensor_product(&h_b)?` | `h_a.tensor(&h_b)?` | HilbertSpace の未採用案 |
| `Spin::try_operator_product(terms)?` | `Spin::try_product(terms)?` | L3 の未採用案 |
| `a.product(&b)?`（複合演算子） | `a.mul(&b)?` | N1 の採用済み・未実装 |
| `Operator::product_all(terms)?` | `Operator::product(terms)?` | N1 の採用済み・未実装 |

### 5.4. 基底・幾何の統合対応表

| 提供型 | 現行 | 推奨案 |
| --- | --- | --- |
| 四つの式型 | `expand_in(&basis)` | `expand(&basis)` |
| 四つの式型 | `coordinates_in(&basis)` | `coordinates(&basis)` |
| Operator | `expand_in_space(&basis, &ambient)` | `expand_in(&basis, &ambient)` |
| Operator | `coordinates_in_space(&basis, &ambient)` | `coordinates_in(&basis, &ambient)` |
| Chain | `nn_pairs_with(&sites)` | `nn_leaf_pairs_with(&sites)` |
| Chain | `nn_pairs_on(&space)` | `nn_leaf_pairs_on(&space)` |
| Chain | `nn_pairs_by_pos_with(&sites)` | `nn_leaf_zip_with(&sites)` |
| Chain | `nn_pairs_by_pos_on(&space)` | `nn_leaf_zip_on(&space)` |

四つの式型は `LocalState`・`LocalOperator`・`State`・`Operator` を指す。局所式や State に明示空間版を追加する提案ではない。Chain の `nn_site_pairs` / `nn_site_pairs_with` / `nn_site_pairs_on` は変更しない。

### 5.5. 内部関数についての結論

第3節の四候補に、追跡用の ID を付ける。公開範囲は現行のままとする。

| ID | 現行名 | 推奨名 | 理由・変更範囲 |
| --- | --- | --- | --- |
| I1 | `LocalOperator::evaluated` | `evaluate` | 状態側と同じ動詞にする。公開する必要はなく、`pub(crate)` のまま呼び出しを置換する。`evaluated_inner_product` は「既に評価済みの入力を要求する」という前提を表すので、この規則で機械的に変更しない。 |
| I2 | `StateExpr::same_node` | `shares_node_with` | 空間側と共有節点の同一性の語彙を揃える。`Arc::ptr_eq` を使う実装を維持し、式の数学的等価性を比較する名前にはしない。 |
| I3 | `basis::enumeration::leaves` | `collect_leaves` | `Vec` への集約と、公開 `HilbertSpace::leaves` の遅延列挙を区別する。既存走査を公開 iterator の利用へ置き換えるかどうかは改名とは別に扱う。 |
| I4 | `qmb-export::coordinates::owned_entries` | `into_flat_entries` | 所有権の移動と flat index 化の両方を表す。借用版の `flat_entries` と対応し、反復時にエラーを返す契約も維持できる。 |

### 5.6. 採用時の移行単位

採用済みの N1 では通常の積の二項 `product` → `mul` と列の `product_all` → `product` を同時に適用し、primitive の二項規則も合わせて整理する。L3 の QuantumSystem の接頭辞省略を採用する場合も、この集約名に揃える。N2・L2 の状態側の移行は完了している。L4 は HilbertSpace 側の短縮名だけを検討する。N3・L7 の基底側も、基本名と明示空間版を一緒に変更する。古い名前を新しい意味で再利用する箇所は、互換 alias を残すだけでは移行できないため、v2 の破壊的変更として対応表を公開する方針を推奨する。

通常の積では旧 `product` が二項、新 `product` が列を受け取るため、同じ型に旧名の互換 alias を残せない。`product_all` → `product` の後に `product` → `mul` を一括置換して、新しい集約名まで変えないようにする。

実装済みの状態構築名と契約は維持する。旧状態 API や状態側の `tensor_all*` の改名案は再導入しない。HilbertSpace 側で `tensor` → `tensor_all` と `tensor_product` → `tensor` の短縮案を採用する場合も、別の既存名の一部や置換済みの新名まで一括置換しないようにする。`expand_in` / `coordinates_in` も同様に、基底だけの呼び出しと明示空間付きの呼び出しを区別する。

名前を再利用しない改名も、新名を正式な入口に揃える。互換期間が必要になった場合に限り、同じシグネチャで旧名から新名へ委譲する deprecated alias を検討し、両名を恒久的な選択肢にはしない。`Dimension` の必須メソッドの改名は独自実装側の更新も必要であり、旧名と新名のデフォルト実装を相互に委譲させて対応しない。

採用後は該当 crate の既存テスト・rustdoc・examples を更新して実行し、第4節の契約を確認する。特に fallible 集約のエラー優先順位、graded の符号、明示空間の優先、基底座標と行列の区別を検証対象にする。本書の未実装の推奨名は、この文書更新では実装・コンパイル検証していない。実装済みの状態構築と未実装の命名案を区別し、解消済み項目を次回の改名対象に含めない。

## 関連文書

- [明示空間による状態構築](../../design/composite-expression.md#121-構築入口)：状態側のテンソル構築に関する実装済みの契約。
- [API の使いやすさに関する残る検討事項](../api-usability/README.md)：補助 API、型推論、利用開始の導線など。
- [局所式の設計](../../design/primitive-expression.md)：構築・評価・数値能力と正規化の境界。
- [複合式の設計](../../design/composite-expression.md)：集約・束縛・ordinary / graded の契約。
- [Hilbert 空間の設計](../../design/hilbert-space.md)：因子・順序・合成規則。
- [幾何の設計](../../design/geometry.md)：Site 対と葉対、モデル構築。
