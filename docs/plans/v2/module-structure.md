# モジュール構造整理：現状調査

作成日: 2026-10-29
更新日: 2026-10-29
調査時の HEAD: `33049f1`

## 1. 対象と読み方

[v2 の先行する構造整理](README.md)に向け、6つのライブラリ crate の現在のモジュール構成、公開範囲、責任分担と接続を記録する。根拠は各 crate の `src/`、`Cargo.toml`、テストの実装である。設計文書にある将来の契約と、現在の実装を区別する。

以下の「公開」「非公開」はモジュールの外部からの可視性を示す。非公開モジュールでも、公開型の再公開や、その型への公開メソッド・trait 実装を通じて外部 API を提供している。ファイルの配置と API の利用パスは必ずしも一致しない。

本書は構造変更前の調査記録とする。変更後のモジュール名、分割単位、公開範囲は今後決定する。

## 2. crate 間の責任と依存

通常の workspace 内依存は次のとおり。dev-dependency と examples はこの表の依存欄に含めない。宣言は各 crate の `Cargo.toml` にある。

| crate | 現在の責任 | workspace 内の通常依存 | 実行環境の宣言 |
| --- | --- | --- | --- |
| [qmb-core](../../../crates/qmb-core/Cargo.toml) | 抽象契約、有限線型結合、式・空間、評価、基底処理 | なし | `no_std`、`alloc` 使用 |
| [qmb-spin-half](../../../crates/qmb-spin-half/Cargo.toml) | spin-1/2 の primitive・局所規則・基底、Chain 上のモデル構築 | core、geometry | `no_std`、`alloc` 使用 |
| [qmb-basic](../../../crates/qmb-basic/Cargo.toml) | UUID Index の発行、標準 Index を使う空間の型別名 | core | `no_std` 宣言なし |
| [qmb-geometry](../../../crates/qmb-geometry/Cargo.toml) | Chain の幾何と、外部の Site 列・空間への借用接続 | core | `no_std` |
| [qmb-normalize](../../../crates/qmb-normalize/Cargo.toml) | norm、状態の正規化、正規化期待値 | core | std を利用する構成 |
| [qmb-export](../../../crates/qmb-export/Cargo.toml) | 基底・空間の解釈を伴う数値コンテナへの変換 | core | std 使用 |

Core は具体的代数、UUID 発行、幾何、正規化拡張、数値 backend に依存しない。spin-half はモデル構築で geometry を使う。export の spin-half 依存はテスト用である。`no_std` の記載はソースと依存宣言の確認に基づき、本調査で専用 target のビルドは行っていない。

## 3. qmb-core

### 3.1. モジュール構成

[lib.rs](../../../crates/qmb-core/src/lib.rs) と [expr/mod.rs](../../../crates/qmb-core/src/expr/mod.rs) がモジュール宣言と再公開を行う。

```text
qmb_core                          lib.rs
├── coefficient                  pub mod
├── primitive                    pub mod
├── ops                          pub mod
├── error                        pub mod
├── linear_combination           pub mod
├── space                        pub mod
├── system                       pub mod
├── basis                        pub mod（basis.rs）
│   └── coordinates              非公開
├── expr                         pub mod
│   ├── parts                    pub mod
│   ├── composite                pub mod
│   ├── local                    非公開
│   ├── operators                非公開
│   ├── evaluate                 非公開
│   ├── composite_eval           非公開
│   │   └── tests                cfg(test)、非公開、ファイル内定義
│   ├── inner_product            非公開
│   ├── expansion                非公開
│   ├── tests                    cfg(test)、非公開
│   └── composite_tests          cfg(test)、非公開
└── test_coefficient              cfg(test)、非公開、tests/support から読込
```

### 3.2. トップレベルと basis の責任

| モジュール／ファイル | 主な型・関数 | 責任範囲 |
| --- | --- | --- |
| [crate root](../../../crates/qmb-core/src/lib.rs) | 各モジュールの主要型・trait の再公開 | crate の案内、`no_std`・`alloc` の設定、公開窓口。具体的な代数や評価の実装は置かない |
| [coefficient](../../../crates/qmb-core/src/coefficient.rs) | `Conjugate`、`RealPart`、`ImaginaryUnit`、`Divide` | 独立した係数能力の契約と `Complex<T>`・実数への実装。内部向けの共役・加算・乗算・零の補助関数も保持する |
| [primitive](../../../crates/qmb-core/src/primitive.rs) | `Identity`、`SpanDimension` | primitive の恒等演算子生成と、その型全体が張る線型空間の次元の契約 |
| [ops](../../../crates/qmb-core/src/ops.rs) | `InnerProduct`、`Product`、`Act`、`Adjoint`、`Trace` | primitive が提供する局所的な演算の契約。積・作用・随伴は `LinearCombination` を返し、具体的な規則は利用側の代数が実装する |
| [error](../../../crates/qmb-core/src/error.rs) | `DimensionOverflow` | 共通の次元 overflow エラー。このモジュールに core の全エラーが集まっているわけではない |
| [linear_combination](../../../crates/qmb-core/src/linear_combination.rs) | `LinearCombination`、`Term`、`IntoTerms` | 順序付き有限線型結合の保持、構築、走査、map、加算、係数の乗算。重複・零係数を保存する。`Term` は `ScaleExpr` の別名で、格納には `SumExpr` を使う |
| [space](../../../crates/qmb-core/src/space.rs) | `HilbertSpace`、`SpaceError`、`Indices`、`Leaves` | Index を持つ葉と順序付き Tensor の共有木構造、構築時の重複検証、互換性・部分系判定、次元取得、借用走査、反復的な破棄 |
| [system](../../../crates/qmb-core/src/system.rs) | `QuantumSystem`、`LocalStateOf`、`LocalOperatorOf`、`StateOf`、`OperatorOf` | primitive・係数型の選択と、局所式・配置・和・積・tensor の構築入口。既存の式 constructor に委譲する。空間や Index の保持・発行、評価は担当しない |
| [basis](../../../crates/qmb-core/src/basis.rs) | `Basis`、`state_basis_len`、`operator_basis_len` | 局所基底の展開・Index・順序付き列挙の契約、primitive の span に基づく基底要素数の取得。`coordinates` の公開型も再公開する |
| [basis::coordinates](../../../crates/qmb-core/src/basis/coordinates.rs) | `Coordinates`、`BasisError` | 集約済み疎座標と基底・空間情報の保持、座標走査、平坦 Index と総次元の取得。展開・集約の計算本体は `expr::expansion` が担当する |

現在の `HilbertSpace<I>` は Index と木構造だけを保持する。局所代数、局所次元、grading は保持せず、次元取得では `SpanDimension` または明示した局所次元を受け取る。木の括り方は保持し、互換性は葉の Index 列と順序で判定する。

### 3.3. expr の責任

| モジュール／ファイル | 主な型・処理 | 責任範囲 |
| --- | --- | --- |
| [expr](../../../crates/qmb-core/src/expr/mod.rs) | 子モジュール宣言、式型・エラーの再公開 | `LocalState`・`LocalOperator` と `State`・`Operator` を同じ窓口から公開する |
| [expr::parts](../../../crates/qmb-core/src/expr/parts.rs) | `ScaleExpr`、`SumExpr`、`IntoOperands`、`ProductExpr`、`ActExpr`、`AdjointExpr` | 式の内容を入れる汎用の構文コンテナ。`SumExpr` は空・単項を inline に保持し、順序付き列の走査・map・連結を提供する。代数評価は行わない |
| [expr::local](../../../crates/qmb-core/src/expr/local.rs) | `LocalState`、`LocalOperator`、内部の節点 enum | 局所式の共有ハンドル、遅延した和・係数・積・作用・随伴の構築、有限線型結合の取り込み、Clone・Default・反復的な Drop。状態は未評価の作用の有無も保持する |
| [expr::composite](../../../crates/qmb-core/src/expr/composite.rs) | `State`、`Operator`、`ExpressionError`、内部の式ハンドル・節点 enum | 局所式の葉への配置、状態の空間保持、tensor・和・積・作用の遅延構築、演算子の束縛と作用空間の選択・検証、反復的な Drop。局所式の `.at()` もここで実装する |
| [expr::operators](../../../crates/qmb-core/src/expr/operators.rs) | `Add`・`Sub`・`Mul`・`Neg`、借用対応 macro、commutator | 4種類の式型の演算子構文。所有／借用の組合せ、係数によるスカラー倍、標準複素係数の左スカラー倍を提供し、式の構築処理に委譲する。commutator の名前付きメソッドもここに置く |
| [expr::evaluate](../../../crates/qmb-core/src/expr/evaluate.rs) | 局所状態の `evaluate`・内積・期待値、局所演算子の内部評価・trace、`LinearTerms` | 局所規則に基づく反復的な評価とスカラーへの縮約。状態作用、演算子の積・随伴の評価、線型項の共通走査を同じファイルが担当する |
| [expr::composite_eval](../../../crates/qmb-core/src/expr/composite_eval.rs) | 複合状態の `evaluate`・内積・期待値、演算子の trace、`OpTerms` | 配置済み演算子の作用を局所評価につなぐ。未変更の状態部分を共有し、内積は `inner_product::contract` に委譲する。演算子項の分配・同一 Index 上の局所積の組立て、spectator の恒等演算子を含む full trace も担当する |
| [expr::inner_product](../../../crates/qmb-core/src/expr/inner_product.rs) | `contract`、`FactorCounts`、因子境界の照合 | 評価済み複合状態の縮約。一致する tensor 境界ごとに内積を計算し、必要な和だけを分配する。局所内積は `evaluate` 側の処理を使う |
| [expr::expansion](../../../crates/qmb-core/src/expr/expansion.rs) | 4式型の `expand_in`・`coordinates_in`、明示空間版、`HilbertSpace::basis_states`・`basis_operators` | 局所／複合の基底展開、係数集約、厳密な零の除去、疎座標取得、通常の式木への再構築、tensor product 基底の遅延列挙と列挙結果の検証。局所・複合の項走査を再利用する |

`evaluate.rs` は局所式向け、`composite_eval.rs` は配置済み複合式向けである。両方に状態の評価と演算子側の処理があり、状態／演算子だけで完全にファイルを分けた構成ではない。

状態への積の作用は右側の演算子から順に適用し、primitive の `Product` を要求しない。通常の評価は基底を選ばず、重複項の集約は基底展開・座標取得が担当する。`LocalOperator::evaluated` は内部メソッドであり、公開の状態 `evaluate` と同じ公開範囲ではない。

### 3.4. 公開パスとモジュール間の接続

主要な trait・型・型別名は crate root から再公開される。`LocalState`・`LocalOperator` の定義先 `expr::local` は非公開で、利用パスは `qmb_core::LocalState` または `qmb_core::expr::LocalState` などとなる。`State`・`Operator`・`ExpressionError` は公開モジュール `expr::composite` にもパスを持つ。`Coordinates`・`BasisError` は非公開の `basis::coordinates` から `basis` と crate root に再公開する。構文パーツは `expr::parts` で公開し、crate root には再公開していない。

実装上の接続として、次がある。

- `space` と `linear_combination` は、`expr::parts::SumExpr` を格納・走査に使う。`linear_combination::Term` も同じ場所の `ScaleExpr` を使う。
- `ops` が返す `LinearCombination` を、式の `from_linear_combination` で通常の構文へ取り込む。
- `system` は局所・複合の式 constructor を呼び、`ExpressionError` を使う。
- `local` の fallible 集約も `composite` で定義された `ExpressionError` を使う。
- `basis::coordinates` は `ExpressionError` を包み、`expr::expansion` が `Basis`・`Coordinates`・`BasisError` を使う。
- `expr::expansion` は `evaluate::LinearTerms` と `composite_eval::OpTerms` を再利用する。`HilbertSpace` の基底列挙の実装もこのモジュールにある。
- 式の節点・内部ハンドルは主に `pub(super)` で `expr` 配下に共有され、評価・演算子構文・基底処理から参照する。

このため、現在の module 間の参照は「契約 → 空間 → 式 → 基底」の一方向だけでは表せない。

## 4. qmb-spin-half

```text
qmb_spin_half                    lib.rs
├── state                        pub mod
├── operator                     pub mod
├── basis                        pub mod
├── coefficient                  非公開
└── utility                      非公開
```

| モジュール／ファイル | 主な型・処理 | 責任範囲 |
| --- | --- | --- |
| [crate root](../../../crates/qmb-spin-half/src/lib.rs) | `SpinHalf<C>`、`QuantumSystem` 実装、主要型の再公開 | 状態・演算子 primitive と係数型の選択。精度の既定値はなく、共通式型は core が提供する |
| [state](../../../crates/qmb-spin-half/src/state.rs) | `SpinHalfState::{Up, Down}` | 状態 primitive、局所内積、演算子の作用、状態 span の次元2 |
| [operator](../../../crates/qmb-spin-half/src/operator.rs) | `SpinHalfOperator::{I, Sx, Sy, Sz, SPlus, SMinus}` | 演算子 primitive、局所積・随伴・通常の trace・恒等演算子、演算子 span の次元4 |
| [basis](../../../crates/qmb-spin-half/src/basis.rs) | `SzBasis`、`CartesianBasis`、`LadderBasis` | 状態基底と2種類の演算子基底の展開・Index・順序付き列挙 |
| [coefficient](../../../crates/qmb-spin-half/src/coefficient.rs) | `half`、`quarter`、`imaginary` | 局所規則で使う分数・虚係数の内部生成補助。core の係数能力の定義とは別の責任 |
| [utility](../../../crates/qmb-spin-half/src/utility.rs) | `SpinHalf<C>::xyz`・`heisenberg`・`neel`、`UtilityError`、`model_space` | Chain と外部の Site 順空間列からモデル・積状態を構築する。Site 数・Index 重複を検証し、隣接 Site の全葉ペアに結合を配置する |

`state`・`operator`・`basis` の主要型と `UtilityError` は crate root に再公開する。モデル構築は非公開の `utility` に実装した `SpinHalf<C>` の関連関数として公開される。

状態・演算子 primitive は係数精度を保持しない。局所規則は core の能力 trait を実装して `LinearCombination` を返し、共通式木や評価器はこの crate に重ねて実装していない。geometry との接続は `utility` にある。

## 5. qmb-basic

子モジュールはなく、すべて [lib.rs](../../../crates/qmb-basic/src/lib.rs) にある。

| 対象 | 責任範囲 |
| --- | --- |
| `Id128` | UUID v4 による Index。同一性の比較・コピー・hash・順序・表示と、`new`／`Default` による発行 |
| `HilbertSpace` 型別名 | `qmb_core::HilbertSpace<Id128>` を標準 Index の空間として提供する |

乱数源と UUID の具体実装はこの crate とその依存が担当する。Index と局所代数を登録・対応付ける機構や、式の実装はない。

## 6. qmb-geometry

子モジュールはなく、すべて [lib.rs](../../../crates/qmb-geometry/src/lib.rs) にある。

| 対象 | 責任範囲 |
| --- | --- |
| `Boundary`、`Chain` | Open／Periodic と Site 数を保持し、自己辺・重複を除いた最近接 Site 対を列挙する |
| `ChainError` | 外部の列の長さと Site 数の不一致を表す |
| `nn_site_pairs_with` | 任意の外部 Site 列を検証し、隣接要素への参照を返す |
| `nn_pairs_with`、`nn_pairs_by_pos_with` | 外部の Site 順 `HilbertSpace` 列を受け、隣接 Site 間の全葉ペアまたは同じ位置の葉ペアを列挙する |
| `nn_site_pairs_on`、`nn_pairs_on`、`nn_pairs_by_pos_on` と内部の `sites_on` | 空間の直下の子を Site 列として読み、上記の列向け処理へ接続する。葉自体は1 Site とする |

幾何の定義と core の空間への接続補助を同じ root に置く。`Chain` の保持情報は Site 数・境界条件だけであり、Index・空間・代数は保持しない。入力空間の Index 重複検証は空間・モデルの構築側に任せる。

## 7. qmb-normalize

子モジュールはなく、すべて [lib.rs](../../../crates/qmb-normalize/src/lib.rs) にある。

| 対象 | 責任範囲 |
| --- | --- |
| `Norm` | 局所・複合状態の自己内積の実部から平方根を取り、norm を取得する |
| `Normalize` | 自己内積が厳密な零なら拒否し、逆 norm による共有式のスカラー倍で正規化する |
| `NormalizedExpectation` | 局所・複合状態の期待値を自己内積全体で除算する。平方根・実部取得は要求しない |
| `NormalizedExpectationIn` | 複合状態について明示した全作用空間を優先して正規化期待値を計算する |
| `NormalizationError`、内部の `nonzero` | 厳密な零 norm と core の `ExpressionError` を扱う |

公開 trait、その局所・複合状態向け実装、エラー、内部補助が同じファイルにある。内積・非正規化期待値・norm squared は core に委譲する。平方根を使う処理だけが係数の実数成分に `Float` を要求する。

## 8. qmb-export

```text
qmb_export                       lib.rs
├── data                         非公開
├── error                        非公開
├── vector                       非公開
├── matrix                       非公開
├── backend_ndarray              非公開、feature = "ndarray"
└── backend_sprs                 非公開、feature = "sprs"
```

| モジュール／ファイル | 主な型・処理 | 責任範囲 |
| --- | --- | --- |
| [crate root](../../../crates/qmb-export/src/lib.rs) | `Export<T, Args>`、型・エラー・外部 backend crate の再公開 | 変換の共通 trait と公開窓口、feature による backend の有効化 |
| [data](../../../crates/qmb-export/src/data.rs) | `Exported`、`Layout`、`VectorMetadata`、`MatrixMetadata`、`MatrixArgs`、内部の `Ambient` | 数値データと基底・空間・shape の解釈、行列出力の状態型・基底・明示空間の指定。結果型、引数型、空間指定の内部抽象を同居させる |
| [error](../../../crates/qmb-export/src/error.rs) | `ExportError`、内部の `reserve` | core・backend のエラーと shape・サイズ・割り当て失敗の変換。Vec の要素数・byte size の検証と確保補助も担当する |
| [vector](../../../crates/qmb-export/src/vector.rs) | `VectorStorage`、`flat_index`、`owned_entries`、ベクトル向け `Export` 実装 | 局所・複合状態から疎座標を取得し、メタデータを組み立て、出力先へ渡す。所有／借用 `Coordinates` の出力、平坦 Index 化、Vec の密ベクトル生成も担当する |
| [matrix](../../../crates/qmb-export/src/matrix.rs) | `MatrixStorage`、`MatrixBuilder`、`DenseBuilder`、行列向け `Export` 実装 | 作用空間・状態基底を選び、各基底状態への作用と座標取得で行列の列を生成する。出力先の内部抽象と、Vec の row-major 密行列生成も担当する |
| [backend_ndarray](../../../crates/qmb-export/src/backend_ndarray.rs) | `Array1`・`Array2` 向け実装、`ArrayBuilder` | 取得済み Vec 出力から ndarray への所有変換、および式からの直接出力。密データの生成には `vector`・`matrix` の Vec 向け処理を使う |
| [backend_sprs](../../../crates/qmb-export/src/backend_sprs.rs) | `CsVec` 向け実装、`SparseBuilder` | 疎座標から疎ベクトルを構築し、列ごとの疎座標から CSC 行列を構築する。密データを中間生成しない |

自前の子モジュールはすべて非公開で、公開 trait・結果型・メタデータ・引数型・エラーは root に集約する。`ndarray`・`sprs` の root での公開は、条件付きの外部 crate の再公開である。既定 feature は空で、Vec 出力は backend feature なしで利用できる。

`matrix` は座標を所有して取り出す補助を `vector::owned_entries` から利用する。backend は `vector`／`matrix` にある内部 storage trait を実装するため、共通処理と backend の接続はこれらのモジュールを境界とする。行列の列生成に演算子 primitive の `Product` は要求しない。

演算子基底における座標取得は core の `Operator::coordinates_in` が提供する。export の行列生成は状態基底への作用を用い、演算子基底座標をそのまま行列として出す処理ではない。

## 9. テストと examples の配置

製品コードのモジュールを動かす際に影響するため、テストの配置・共有も調査した。以下は既存の配置の記録であり、本調査では Rust のテストを再実行していない。

| 対象 | 配置と主な検証範囲 |
| --- | --- |
| core の unit test | [expr/tests.rs](../../../crates/qmb-core/src/expr/tests.rs) が局所式、[expr/composite_tests.rs](../../../crates/qmb-core/src/expr/composite_tests.rs) が複合式の内部共有・遅延構築・評価・破棄を検証する。`composite_eval.rs` 内にも項走査の能力制約のテストがある |
| core の integration test | [tests/](../../../crates/qmb-core/tests) に係数・primitive・有限線型結合・空間・構築 API・fallible 集約・局所代数・基底・割り当てのテストがある |
| spin-half | [tests/](../../../crates/qmb-spin-half/tests) に局所規則・代数則・基底・精度・最小の係数能力・モデル構築・複合式と独立な行列との照合がある |
| basic／geometry | [basic/tests/](../../../crates/qmb-basic/tests) が Index 発行と core 接続、[geometry/tests/](../../../crates/qmb-geometry/tests) が境界条件と外部空間・Site 列への接続を検証する |
| normalize | [tests/](../../../crates/qmb-normalize/tests) に能力制約・正規化・正規化期待値・零 norm・極端な実分母の検証がある |
| export | [tests/](../../../crates/qmb-export/tests) に公開契約・密疎出力・shape／サイズ overflow・割り当て・非正規直交基底での照合があり、[tests/support/mod.rs](../../../crates/qmb-export/tests/support/mod.rs) に独立な行列演算などの補助を置く |
| doctest | 各 crate の rustdoc に利用例を置く。core・spin-half・export には型や係数精度の不一致などを拒否する compile-fail の例もある |

[core/tests/support/coefficient.rs](../../../crates/qmb-core/tests/support/coefficient.rs) は integration test と、`lib.rs` の `#[path]` を使った unit test の両方から読み込む。[core/tests/support/three_level.rs](../../../crates/qmb-core/tests/support/three_level.rs) は core のテストに加え、normalize・export のテストからも相対パスで読み込む。テスト補助専用の共有 crate はない。

Workspace の [examples](../../../examples/Cargo.toml) はライブラリではなく、`two_spins`・`basis_expansion`・`heisenberg` の3つの bin と [tests/chain.rs](../../../examples/tests/chain.rs) を持つ。手作業の構築、基底展開、Chain モデルと正規化・密疎出力の接続を扱う。`normalize` は既定 feature、`export` は任意 feature である。

## 10. 構造整理で検討する境界

次は実装から確認した責任の同居・横断と、それに対応する検討事項である。具体的な分割方針の採用を意味しない。

| 確認した現状 | 構造整理での検討事項 |
| --- | --- |
| core の `local`・`composite` は型・内部節点・構築・破棄を持ち、`composite` は配置・束縛・検証も持つ | 型と構文の保持、構築・空間接続、所有・破棄をどこまで同じ単位に置くか |
| `evaluate`・`composite_eval` に状態評価・演算子処理・縮約 API・再利用される項走査が同居する | 局所／複合という分割軸と、作用・項走査・縮約という責任の対応 |
| `basis` が基底契約・座標型、`expr::expansion` が4式型の展開・座標取得と空間の基底列挙を担当する | 基底契約、座標データ、式への接続、列挙の配置 |
| `space`・`linear_combination` が `expr::parts` を使い、基底側のエラーも式側のエラーを参照する | 共有コンテナと各責任の依存方向、エラーの配置 |
| basic・geometry・normalize は子モジュールを持たず、root に公開型・実装を置く | 各 crate の責任に応じた分割と root の役割 |
| spin-half の root が系の型を定義し、非公開 `utility` がその関連関数としてモデル構築を実装する | 具体的代数・系の構築入口・幾何を使うモデル構築の配置 |
| export の `data` は結果・メタデータ・引数を持ち、`vector`・`matrix` は変換経路と密 storage 実装を持つ | 公開データ型、出力引数、座標変換、storage・backend の責任分担 |
| 公開モジュールのパス、root への再公開、非公開 module 内の公開メソッドが併存する | 内部配置の変更と公開 API の扱いを区別して決める必要がある |
| unit test が内部節点を参照し、テスト補助の相対パス共有が crate をまたぐ | 製品 module の移動に伴うテスト配置と共有方法 |

Graded tensor・fermion・spin-one は [v2 の確定範囲](README.md)に含まれるが、現在の製品コードには grading／parity の能力・節点、fermion・spin-one の具体的代数はない。構造整理後にこれらをどの module・crate に置くかは、別途設計する。
