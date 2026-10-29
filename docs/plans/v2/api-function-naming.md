# 関数名・公開 API 名の整理候補

作成日: 2026-10-30
更新日: 2026-10-30

現在の実装を読み、同種の操作で名前の規則が揃っていない箇所と、公開 API として長くなっている名前を記録する。実装名・入力・動作はソースで確認した事実であり、使いにくさの評価・優先度・新しい名前は検討案である。改名の採用や v2 での実施は未決とする。

確認対象は `crates/` 以下の `qmb-core`・`qmb-basic`・`qmb-geometry`・`qmb-spin-half`・`qmb-spin-one`・`qmb-fermion`・`qmb-normalize`・`qmb-export` の関数・メソッド・公開 trait と、`examples/src/bin/` の呼び出しである。非公開 module 内でも公開型の inherent method や再公開された trait は公開 API に含める。内部関数は後段に分ける。

優先度「高」は日常的な構築・計算で長さや名前の衝突が目立つ候補、「中」は API 間の規則を揃える候補、「低」は現行名を残す判断も有力な候補とする。識別子の長さは `_` を含み、型名・module path は含めない。文字数だけを改名の基準にはしない。

## 1. 名前の規則が揃っていない箇所

### N1. 二項演算と列の集約の区別（中）

| 操作 | 二項版 | 列の集約版 |
| --- | --- | --- |
| 状態の和 | `State::add` | `State::sum` / `try_sum` |
| 演算子の積 | `LocalOperator::product`、`Operator::product` | `product_all`、複合式には `try_product_all` もある |
| 状態・空間のテンソル積 | `State::tensor_product`、`HilbertSpace::tensor_product` | `tensor`、状態には `try_tensor` もある |
| QuantumSystem の演算子積 | 二項専用の入口はない | `operator_product` / `try_operator_product` |

列を表す印が積では `_all`、テンソル積では無印になり、QuantumSystem に移ると `_all` が消える。各シグネチャは明確だが、別の型で同じ操作を探す際に規則を推測しにくい。

候補は二項版を `product` / `tensor`、列を `product_all` / `tensor_all` に揃え、QuantumSystem 側も列であることを同じ規則で示すこと。ただし現行の `tensor` は列を受け取るため、`tensor_product` だけを `tensor` に変えることはできない。和の `add` / `sum` は既に二項と列を区別できており、同じ語尾へ機械的に揃える必要はない。

根拠: [State の構築](../../../crates/qmb-core/src/expr/composite.rs#L261)、[Operator の積](../../../crates/qmb-core/src/expr/composite.rs#L529)、[局所演算子の積](../../../crates/qmb-core/src/expr/local.rs#L233)、[空間のテンソル積](../../../crates/qmb-core/src/space.rs#L181)、[QuantumSystem の積](../../../crates/qmb-core/src/system.rs#L182)。

### N2. graded・配置・種類指定の修飾順序（中）

`tensor_g` / `tensor_product_g` / `state_tensor_g` は末尾の `_g` が graded を表す。一方、配置も行う入口は `state_tensor_g_at`、種類指定版は `state_tensor_at_with_kind` で、`_at` と合成規則の修飾順序が異なる。

`_g` は短いものの `TensorKind::Graded` との対応を名前だけでは読み取りにくい。略記を維持するなら「基本操作 → 配置 → 合成規則」の順にし、例えば `state_tensor_at_g` と `state_tensor_at_with_kind` を対応させる案がある。`graded_*` のように綴る案は分かりやすさと名前の長さを比較する。ordinary / graded を名前から消して、入力の primitive や既存空間から暗黙に推測する変更は行わない。

根拠: [QuantumSystem の各入口](../../../crates/qmb-core/src/system.rs#L65)、[State の各入口](../../../crates/qmb-core/src/expr/composite.rs#L305)、[空間の各入口](../../../crates/qmb-core/src/space.rs#L181)。

### N3. `_in` が基底と作用空間の両方を表す（中）

| 現行名 | `_in` が指す入力 |
| --- | --- |
| `LocalState` / `LocalOperator` / `State` / `Operator` の `expand_in(&basis)`、`coordinates_in(&basis)` | 基底 |
| `State::act_in` / `expectation_in`、`Operator::trace_in` / `hilbert_schmidt_inner_product_in`、`normalized_expectation_in` | 明示する作用空間 |
| `Operator::expand_in_space(&basis, &ambient)` / `coordinates_in_space(&basis, &ambient)` | 基底に加えて明示する作用空間 |
| `MatrixArgs::in_space(&ambient)` | 行列出力で優先する作用空間 |

基底処理だけは基本名の時点で `_in` を使い、明示空間を加えると `_in_space` になる。特に `coordinates_in_space` は20文字で、他の明示空間版から名前を推測しにくい。

候補は基底を受け取る基本名を `expand(&basis)` / `coordinates(&basis)`、明示空間版を `expand_in(&basis, &ambient)` / `coordinates_in(&basis, &ambient)` に揃えること。`Basis::expand` とも語彙が揃う。同名メソッドの引数数が変わるため、局所式・複合式をまとめて移行する必要がある。`MatrixArgs::in_space` は builder の役割が明瞭であり、計算メソッドと同じ形へ無理に短縮しなくてもよい。

根拠: [局所式の基底処理](../../../crates/qmb-core/src/basis/expansion/local.rs#L93)、[複合演算子の基底処理](../../../crates/qmb-core/src/basis/expansion/composite.rs#L246)、[作用・期待値](../../../crates/qmb-core/src/evaluate/inner_product.rs#L327)、[trace・内積](../../../crates/qmb-core/src/evaluate/trace.rs#L60)、[正規化期待値](../../../crates/qmb-normalize/src/expectation.rs#L23)、[MatrixArgs](../../../crates/qmb-export/src/args.rs#L24)。

### N4. 次元を表す `dim` / `dimension` / `basis_len`（中）

primitive は `Dimension::dim()`、空間・座標は `dimension()` / `local_dimension()` を使う。また `state_basis_len::<S>()` / `operator_basis_len::<O>()` は基底を受け取らず、内部で primitive の `dim()?.get()` を返している。

同じ次元の問い合わせに略語・完全形・長さという三つの語彙がある。候補は `dimension` に統一し、二つの自由関数も `state_dimension` / `operator_dimension` にするか、primitive の次元問い合わせへ集約すること。primitive の戻り値の `NonZeroUsize` と、空間・自由関数の `usize` の違いは維持する。

`Coordinates::len()` は保存された非零座標の数、`dimension()` は座標空間全体の次元であり、こちらの区別は必要である。primitive の次元も enum の variant 数ではなく張る線型空間の次元なので、単なる `len` にはしない。

根拠: [Dimension](../../../crates/qmb-core/src/ops.rs#L81)、[HilbertSpace](../../../crates/qmb-core/src/space.rs#L149)、[基底長の自由関数](../../../crates/qmb-core/src/basis/enumeration.rs#L164)、[Coordinates](../../../crates/qmb-core/src/basis/coordinates.rs#L66)。

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
| L2・高 | `QuantumSystem::state_tensor_at_with_kind`（25）、`state_tensor_with_kind`（22）、`try_state_tensor_with_kind`（26） | `state`・`tensor`・配置・種類指定・fallible 入力が積み重なる。現在この trait の tensor 構築は状態だけなので、`state_` を省く余地がある。N1 の列名と組み合わせて検討し、例えば現行の種類指定付き列集約を `tensor_all_with` / `try_tensor_all_with`（15 / 19）へ整理する案。配置版も同じ語順で揃える。 |
| L3・中 | `QuantumSystem::operator_product_at`（19）、`try_operator_product`（20） | 基本構築では繰り返し現れる。現在の積構築は演算子だけなので、N1 と合わせて `product_all_at` / `try_product_all`（14 / 15）などへ寄せる余地がある。`state_sum` と `operator_sum` は両方存在するため、全入口から型の接頭辞を取り除くことはできない。 |
| L4・高 | `State::tensor_product_with_kind`、`HilbertSpace::tensor_product_with_kind`（24）、`State::try_tensor_with_kind`（20） | 二項／列の区別と種類指定が長さを増やす。N1 の整理後なら二項の `tensor_with(kind, &rhs)`（11）、列の `try_tensor_all_with(kind, states)`（19）などが候補。引数の `TensorKind` で明示性を保てるか確認する。列版と二項版の両方を同じ `tensor_with` にしない。 |
| L5・中 | `LocalState::from_linear_combination`、`LocalOperator::from_linear_combination`（23） | 変換元の型名をすべて綴る。`from_combination`（16）、または `From<LinearCombination<Self, C>>` を実装して `.into()` で変換する案。現在の入力は primitive ではなく局所式の線型結合であり、primitive の規則結果には `map(Self::primitive)` が必要な点を隠さない。`from_terms` に寄せて iterator を受け取るなら、単なる改名とは別に入力型も検討する。 |
| L6・低 | `NormalizedExpectation::normalized_expectation`（22）、`NormalizedExpectationIn::normalized_expectation_in`（25） | `mean` / `mean_in` のような短い名前を候補にできるが、計算の意味を明確にする利点が現行名にはある。core の `expectation` は非正規化なので、正規化版だけを `expectation` へ改名すると既存の inherent method と衝突する。単独の改名は避け、必要なら期待値の名前全体を再検討する。 |
| L7・中 | `Operator::coordinates_in_space`（20）、`Chain::nn_pairs_by_pos_with`（20） | N3 の基底／空間の語尾と、N5 の葉の zip の整理で短縮できる。名前だけを切り詰めず、それぞれのファミリを同時に揃える。 |

根拠:

- L1: [Hilbert–Schmidt 内積](../../../crates/qmb-core/src/evaluate/trace.rs#L33)。局所版と複合版の両方に同じ長い名前がある。
- L2・L3: [QuantumSystem](../../../crates/qmb-core/src/system.rs#L65)。実際の呼び出しは [tensor_products](../../../examples/src/bin/tensor_products.rs#L17)、[two_spins](../../../examples/src/bin/two_spins.rs) にもある。
- L4: [State](../../../crates/qmb-core/src/expr/composite.rs#L334)、[HilbertSpace](../../../crates/qmb-core/src/space.rs#L231)。
- L5: [LocalState](../../../crates/qmb-core/src/expr/local.rs#L137)、[LocalOperator](../../../crates/qmb-core/src/expr/local.rs#L204)、[LinearCombination の既存の from_terms](../../../crates/qmb-core/src/expr/parts/linear_combination.rs#L44)。
- L6: [正規化期待値の trait と実装](../../../crates/qmb-normalize/src/expectation.rs#L15)、[非正規化期待値](../../../crates/qmb-core/src/evaluate/inner_product.rs#L284)、利用例の [heisenberg](../../../examples/src/bin/heisenberg.rs#L22)。
- L7: [座標取得](../../../crates/qmb-core/src/basis/expansion/composite.rs#L326)、[Chain](../../../crates/qmb-geometry/src/geometry/chain.rs#L66)。

## 3. 内部関数の表記揺れ・処理内容が隠れる名前

| 現行名・公開範囲 | 確認した違いと整理候補 |
| --- | --- |
| `LocalOperator::evaluated`（`pub(crate)`）と公開の `LocalState::evaluate` / `State::evaluate` | いずれも元を変更せず評価済みの式を返すが、演算子側だけ過去分詞になっている。内部名も `evaluate` に揃える案。公開 API の追加は別の判断であり、現状の演算子に公開 `evaluate` があるとは扱わない。[局所評価](../../../crates/qmb-core/src/evaluate/local.rs#L138)、[状態評価](../../../crates/qmb-core/src/evaluate/composite.rs#L41)を参照。 |
| `StateExpr::same_node`（`pub(crate)`）と公開の `HilbertSpace::shares_node_with` | どちらも `Arc::ptr_eq` で共有節点の同一性を調べる。`shares_node_with` に揃えると数学的な同値判定と区別しやすい。[式](../../../crates/qmb-core/src/expr/composite.rs#L181)、[空間](../../../crates/qmb-core/src/space.rs#L171)を参照。 |
| `basis::enumeration::leaves`（`pub(super)`）と公開の `HilbertSpace::leaves` | 前者は葉をすべて `Vec` へ集め、後者は遅延 iterator を返す。内部の前者を `collect_leaves` にすると確保・集約を伴うことが伝わる。[集約](../../../crates/qmb-core/src/basis/enumeration.rs#L15)、[iterator](../../../crates/qmb-core/src/space.rs#L138)を参照。 |
| `qmb-export` の `coordinates::owned_entries`（`pub(crate)`） | 所有権の移動だけでなく、反復時に多桁の座標を flat index に変換し、検証・overflow 検出を行う。`into_flat_entries` とすれば、core の借用版 `Coordinates::flat_entries` と処理内容を対応づけられる。[変換](../../../crates/qmb-export/src/coordinates.rs#L21)、[借用版](../../../crates/qmb-core/src/basis/coordinates.rs#L100)を参照。 |

## 4. 整理時に維持する区別と検証

- **`try_` は単なる Result の印ではない。** `State::tensor` / `Operator::product_all` も `Result` を返すが、`try_tensor` / `try_product_all` は入力が `Result` の列である。最初の入力エラーで消費を止め、全入力成功後に構成を検証するという一貫した契約がある。名前を揃えるためにすべての fallible 操作へ `try_` を付けない。[契約](../../../crates/qmb-core/src/system.rs#L37)と[検証対象の既存テスト](../../../crates/qmb-core/tests/fallible_construction.rs)を参照。
- **名前が違っても、受け手や入力が違えば必要な区別がある。** `QuantumSystem::state` / `operator` は primitive と係数型を選択済みの構築入口、`LocalState::primitive` / `LocalOperator::primitive` は式型の入口である。`local.at(&factor)` と `State::place(&local, &factor)` / `Operator::place(&local, &factor)` はメソッドと関連関数の入口である。これらを長い `local_*` 名へ機械的に揃える優先度は低い。
- **短縮で数学的契約を変えない。** 左側の共役・随伴、積の右側が先に作用する順序、ordinary / graded の選択、正規化の有無、明示空間が古い束縛に優先する規則を維持する。局所式の恒等元と未束縛の複合恒等元、列集約の空入力・単項の扱いも維持する。
- **名前と能力を別に扱う。** primitive・Index の Clone 要求、数値能力、係数精度、no_std の境界を改名に伴って追加・変更しない。`normalized_expectation` の短縮のために `normalized().expectation()` へ置き換えると、平方根など現在不要な能力を要求するので同等の移行にはならない。
- **採用時に名前の衝突と利用例を確認する。** 特に N1・N3 は既存名の引数形を変えるため、同じ型に旧名の互換 alias をそのまま残せるとは限らない。公開 trait の変更は独自代数の実装にも影響する。対応する crate のテスト・rustdoc・examples と、設計文書の API 表を合わせて更新・検証する。本書の調査ではソースと呼び出し例を確認し、実行テストは行っていない。

## 関連文書

- [API の使いやすさに関する残る検討事項](../api-usability/README.md)：補助 API、型推論、利用開始の導線など。
- [局所式の設計](../../design/primitive-expression.md)：構築・評価・数値能力と正規化の境界。
- [複合式の設計](../../design/composite-expression.md)：集約・束縛・ordinary / graded の契約。
- [Hilbert 空間の設計](../../design/hilbert-space.md)：因子・順序・合成規則。
- [幾何の設計](../../design/geometry.md)：Site 対と葉対、モデル構築。
