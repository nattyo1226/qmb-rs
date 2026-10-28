# モジュール構造整理：決定事項

作成日: 2026-10-29
更新日: 2026-10-29

[現状調査](module-structure.md)を踏まえた構造整理の決定事項と、それを具体化する構成案を記録する。以下は変更後の方針であり、現時点の実装とは区別する。構成案の詳細に残る未決事項は、決定済みの配置方針と分けて示す。

## 1. qmb-core

### 1.1. 次元情報を primitive の能力として維持し、Dimension に改名する

次元情報は従来どおり状態・演算子 primitive 型の能力として提供する。`SpanDimension` trait を `Dimension`、取得メソッドを `dim()` に改名する。

その型の値全体が張る線型空間の次元を返す契約は維持する。状態 primitive では d、全演算子空間の primitive では d² を返す。個々の primitive 値に次元を格納するものではない。

返却型は従来どおり `Result<NonZeroUsize, DimensionOverflow>` とし、整数化できない次元は取得時にエラーとする。次元を使う空間・基底・export の処理は、引き続き primitive 型の能力を参照する。

`Dimension` は空間の次元取得、基底の列挙・展開・座標検証、および数値出力のサイズ決定に必要な操作で要求する。通常の代数的な作用評価・内積・trace に次元取得能力は追加要求しない。

次元情報を `QuantumSystem` に移す方針は取り消す。`QuantumSystem` は従来の型選択・構築入口・型別名の役割を維持し、この整理のために式型を Q でパラメータ化する変更は行わない。

### 1.2. Identity trait と operator primitive への実装を維持する

`Identity` trait は維持し、`SpinHalfOperator` などの operator primitive 型に実装する。恒等演算子を表す primitive の取得は、引き続き `Identity::identity()` が担当する。この責任は `QuantumSystem` に移さない。

### 1.3. 局所式の恒等元表現を primitive に統一する

局所式の `OperatorNode::Identity` variant は削除し、恒等演算子は `OperatorNode::Primitive(O::identity())` で表す。

- `O: Identity` のもとで、`OperatorNode<O, C>` と `LocalOperator<O, C>` にも `Identity` trait を実装する。
- `OperatorNode` の実装は `Primitive(O::identity())` を返し、`LocalOperator` の実装はその節点を保持する式を返す。
- 局所式の `identity()` と、空入力で恒等元を生成する `product_all` は `O: Identity` を要求する。係数型 C への数値能力の要求は追加しない。
- 局所式の恒等元の作用・積・随伴・trace・基底展開は、他の primitive と同じ経路で処理する。

複合式の `OpNode::Identity` variant は維持する。これは未束縛の空の積や、選択した全作用空間の恒等演算子を表す。複合式の恒等元構築に `O: Identity` の要求は追加しない。

### 1.4. 複合式にも Identity trait を実装する

複合式の `Operator<I, O, C>` にも `Identity` trait を実装し、既存の `OpNode::Identity` を保持する未束縛の演算子を返す。Index 型 I、operator primitive 型 O、係数型 C に追加の trait bound は要求しない。作用空間は既存の `bind` で指定する。

既存の関連関数 `Operator::identity()` も維持し、trait 経由の生成と同じ表現を返す。これにより、primitive・局所式・複合式を `T: Identity` という共通の契約で扱えるようにする。

`Identity` trait の説明は「恒等演算子を表す値の生成」とし、primitive と式型の両方を対象とする契約へ改める。

### 1.5. Identity・Dimension を ops に移し、primitive モジュールを廃止する

`Identity` と、`SpanDimension` から改名する `Dimension` の定義は、`qmb_core::primitive` から `qmb_core::ops` へ移す。恒等演算子の生成・次元取得を、`Product`・`Act`・`Adjoint`・`Trace` などと同じ数学的能力の領域で扱う。

- `qmb_core::ops::Identity` を公開し、既存の crate root の `qmb_core::Identity` 再公開も維持する。
- `qmb_core::ops::Dimension` を公開し、crate root からも `qmb_core::Dimension` として再公開する。
- `ops` の説明は、次元取得と代数演算の能力を扱うものへ改める。各 trait の契約に応じて、primitive と式型の両方を対象に含める。
- `Identity`・`Dimension` の移動が完了した段階で、現在の `primitive` モジュールを廃止する。

### 1.6. 配置方針

- `LinearCombination` とその項・イテレータを `expr::parts` に置く。独立したトップレベルの `linear_combination` モジュールは、この配置へ移行する。
- 空・単項・複数項を保持する共通の列コンテナと、その走査・連結・map も `expr::parts` に置く。現在の `SumExpr`・`IntoOperands` が該当する。独立した `collections` モジュールは設けない。
- 式の型・共有構造・構築・配置・束縛・破棄・構造の走査は `expr` が担当する。
- 評価・作用・内積の縮約・trace は、`expr` の外にあるトップレベルの `evaluate` にまとめる。
- 基底契約・座標データ・基底列挙・展開・座標取得は `basis` にまとめる。現在の `expr::expansion` の処理もこの領域へ移す。

モジュールの名前空間と実際の依存関係は区別する。`expr::parts` は具体的な式ハンドルや評価器に依存しない共通データ構造とし、`ops` が局所規則の返却型として `LinearCombination` を使い、`space` が子列として共通の列コンテナを使える構成にする。全モジュール間の参照を完全に一方向へ揃えることは、本案の前提にしない。

### 1.7. 構成案と責任範囲

以下は上記の方針を具体化した構成案である。子モジュールの細分化・名前は、実装設計時に調整する。公開範囲は第1.9節に従って現行を維持する。

```text
qmb_core
├── coefficient           係数の能力と標準実装
├── ops                   Dimension・Identity・積・作用・随伴・内積・trace の契約
├── system                QuantumSystem・型選択・構築入口・型別名
├── space                 空間の構造・走査・互換性・部分系判定
├── expr
│   ├── parts             共通データ構造と構文パーツ
│   │   └── linear_combination
│   │                     有限線型結合・項・イテレータ
│   ├── local             局所式の保持・構築・破棄
│   ├── composite         複合式の保持・構築・配置・束縛・破棄
│   ├── arithmetic        Rust の演算子構文
│   └── traversal         再利用する線型項・演算子項の走査
├── evaluate
│   ├── local             局所規則による評価
│   ├── composite         複合式の作用評価
│   ├── inner_product     内積の縮約
│   └── trace             trace・Hilbert–Schmidt 内積
├── basis
│   ├── coordinates       座標データ
│   ├── enumeration       局所・tensor product 基底の列挙
│   └── expansion
│       ├── local         局所式の展開・座標取得
│       └── composite     複合式の展開・座標取得
└── error                 複数領域で使う共通エラー
```

`parts` には有限線型結合・共通列に加え、現在の `ScaleExpr`・`ProductExpr`・`ActExpr`・`AdjointExpr` などの構文パーツも置く。列コンテナの空・単項での無確保、項の順序・重複・明示した零係数を保持する契約は維持する。

### 1.8. 現在の実装からの整理

| 現在の配置 | 移行先・整理方針 |
| --- | --- |
| `primitive::SpanDimension` | primitive の能力として維持し、`ops::Dimension` へ改名・移動する |
| `primitive::Identity` | `ops::Identity` へ移す |
| トップレベルの `linear_combination` | `expr::parts` 配下へ移す |
| `expr::parts` の共通列・構文パーツ | 同じ領域に保持し、有限線型結合と併せて共通データ構造を管理する |
| `expr::operators` | `expr::arithmetic` として演算子構文の役割を明示する |
| `expr::evaluate` の局所評価 | トップレベルの `evaluate` 配下へ移す |
| `expr::composite_eval` の作用・縮約・trace | `evaluate` 配下で処理の責任に応じて分ける |
| `expr::inner_product` | `evaluate::inner_product` へ移す |
| `evaluate::LinearTerms`・`composite_eval::OpTerms` など | `expr::traversal` に集め、評価と基底処理から再利用する |
| `expr::expansion` | `basis` 配下で列挙・局所展開・複合展開に分ける |

基底処理は共通の評価・走査・構築処理を利用する。基底処理専用の評価器を重ねて実装せず、式の内部節点への直接参照は必要な範囲に限定する。

`ExpressionError` は `expr`、`BasisError` は `basis`、`SpaceError` は `space` が担当する。トップレベルの `error` は `DimensionOverflow` など複数領域で使う共通エラーを扱う。

### 1.9. 公開範囲

公開範囲は現行を維持する。既存の公開型・trait・メソッドと crate root の再公開を維持し、内部節点・評価器・内部走査は引き続き外部へ公開しない。合意済みの改名・モジュール移動・廃止は本書の方針に従う。

トップレベルへ移す `evaluate` も非公開モジュールとし、既存の公開型へのメソッド実装を担当する。型の定義先とメソッド実装の配置を一致させることは要求しない。

## 2. qmb-geometry

### 2.1. 承認済みの構成

境界条件・幾何・関係・外部の列や空間への接続を、それぞれ別の責任として分離する。

```text
qmb_geometry
├── boundary
│                         共有する境界条件の定義
├── geometry
│   └── chain             Chain の構造・構築・情報取得
├── relation
│   └── nearest_neighbor  各幾何に対する最近接 Site 関係
└── connection
                          Site 関係を外部の列・HilbertSpace に接続
```

| モジュール | 責任範囲 |
| --- | --- |
| `boundary` | `Open`・`Periodic` など、幾何をまたいで共有する境界条件を定義する。具体的な Site 関係は、幾何と関係の実装が条件を解釈して生成する |
| `geometry::chain` | Site 数・境界条件など Chain 自体の情報と、その構築・情報取得を担当する。今後追加する幾何は `geometry` 配下で `chain` と並べる |
| `relation::nearest_neighbor` | 幾何に対する最近接 Site 関係を生成する。今後追加する関係は `relation` 配下に置く |
| `connection` | Site 対を外部の列の要素への参照へ変換し、Site に対応する Hilbert 空間の葉を全組合せ・位置対応で列挙する。空間の直下の子を Site 列として読む接続も担当する |

最近接 Site 対の生成と、Site に配置された自由度同士の組合せを分離する。接続処理は、別の幾何や関係でも再利用できる責任として扱う。

### 2.2. 既存 API と共通化

既存の `Chain` の関連メソッドは、対応する担当モジュールの処理へ委譲する形で維持する。`Boundary`・`Chain`・`ChainError` の crate root からの公開も維持する。

この承認はモジュール構成と責任分担を対象とする。幾何・関係の共通 trait の具体的な形は、別途検討する。

## 3. qmb-export

### 3.1. 承認済みの構成

公開データ型・出力引数・変換手順・数値コンテナへの格納を分離する。

```text
qmb_export
├── lib.rs
│                         Export trait・公開型の再公開・feature 設定
├── output.rs
│                         Exported：数値データとメタデータの組
├── metadata.rs
│                         VectorMetadata・MatrixMetadata・Layout
├── args.rs
│                         MatrixArgs・作用空間の指定
├── vector.rs
│                         状態・Coordinates のベクトル出力
├── matrix.rs
│                         基底状態への作用による行列の列生成
├── coordinates.rs
│                         共通の平坦 Index 化・座標取り出し
├── backend
│   ├── mod.rs
│   │                     内部 storage trait・構築契約・確保処理
│   ├── vec.rs
│   │                     Vec による密ベクトル・密行列の構築
│   ├── ndarray.rs
│   │                     ndarray への変換［feature: ndarray］
│   └── sprs.rs
│                         疎ベクトル・CSC 行列の構築［feature: sprs］
└── error.rs
                          ExportError・外部エラーの変換
```

### 3.2. 現在の実装からの整理

| 現在の配置 | 移行先・整理方針 |
| --- | --- |
| `data::Exported` | `output` へ移す |
| `data` のメタデータ型・`Layout` | `metadata` へ移す |
| `data::MatrixArgs`・内部の `Ambient` | `args` へ移す |
| `vector`・`matrix` の出力処理 | 同じモジュールで、座標取得・作用空間の選択・列生成などの変換手順を担当する |
| `vector::flat_index`・`owned_entries` | 共通の `coordinates` へ移し、ベクトル・行列の両方から利用する |
| `VectorStorage`・`MatrixStorage`・`MatrixBuilder` | `backend` に集め、コンテナへの格納と構築の内部契約とする |
| Vec 向けの実装・`DenseBuilder` | `backend::vec` へ移す |
| `backend_ndarray`・`backend_sprs` | それぞれ `backend::ndarray`・`backend::sprs` へ移す |
| `error::reserve` | `backend` の共通確保処理へ移す |

`vector`・`matrix` は出力対象と生成手順を、`backend` は数値コンテナへの格納方法を担当する。共通の座標変換を `coordinates` に置き、行列処理がベクトル処理の補助関数を参照する依存を整理する。

### 3.3. 公開範囲と既存の出力契約

子モジュールは非公開とし、`Export` trait、結果型、メタデータ、引数型、エラーは従来どおり crate root から公開する。feature に応じた外部 crate の `ndarray`・`sprs` の再公開も維持する。

既定 feature は空とし、Vec 出力は backend feature なしで利用できる構成を維持する。ndarray は Vec の密データ生成を再利用し、取得済み Vec 出力からの変換も提供する。sprs は疎座標から疎ベクトル・CSC 行列を直接構築し、密データを中間生成しない。

## 4. qmb-spin-half

### 4.1. 承認済みの構成

非公開の `utility` モジュールを `models` に改名し、XYZ・Heisenberg・Neel の構築を担当することを明示する。既存の primitive・基底の分離は維持する。

```text
qmb_spin_half
├── lib.rs                SpinHalf<C>・QuantumSystem 実装・再公開
├── state                 状態 primitive・内積・作用・次元
├── operator              演算子 primitive・局所積・随伴・trace・恒等元・次元
├── basis                 状態基底・演算子基底
├── coefficient           局所規則で使う係数の内部生成補助
└── models                XYZ・Heisenberg・Neel の構築と入力検証
```

### 4.2. 公開範囲と既存 API

`models` は引き続き非公開とし、`SpinHalf<C>::xyz`・`heisenberg`・`neel` の関連関数を維持する。公開エラー名 `UtilityError` と crate root からの再公開も維持する。

`state`・`operator`・`basis` の公開、主要型の crate root からの再公開、および非公開の `coefficient` の配置は変更しない。

## 5. qmb-normalize

### 5.1. 承認済みの構成

局所／複合の区別ではなく、操作の責任に沿ってモジュールを分ける。

```text
qmb_normalize
├── lib.rs
│                         crate の案内・公開 trait とエラーの再公開
├── norm.rs
│                         Norm・Normalize と局所／複合状態への実装
├── expectation.rs
│                         NormalizedExpectation・NormalizedExpectationIn
│                         と局所／複合状態への実装
└── error.rs
                          NormalizationError・内部の零 norm 検証
```

`norm` は平方根を使う norm と状態の正規化を担当する。`expectation` は期待値を自己内積で割る処理を担当し、平方根・実部取得を要求しない契約を維持する。共通エラーと、自己内積が厳密な零であることを検出して拒否する内部処理は `error` にまとめる。

### 5.2. 公開範囲と既存 API

子モジュールは非公開とし、`Norm`・`Normalize`・`NormalizedExpectation`・`NormalizedExpectationIn` と `NormalizationError` は従来どおり crate root から公開する。

既存のメソッド・返却型・操作ごとの能力制約を維持する。内積・非正規化期待値・norm squared は引き続き core に委譲する。

本書は実装前の計画であり、コードへの反映は別途行う。
