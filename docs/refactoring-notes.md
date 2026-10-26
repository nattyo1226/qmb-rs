# リファクタリングの整理メモ（一時）

更新日: 2026-10-26

## 1. 目的と位置づけ

[アーキテクチャ第10節](architecture.md#10-将来の必須課題)の5項目について、現状、検討事項、依存関係、着手順序を整理する作業用文書である。State・Operator・Structure・Space の所有・共有・関連付けと、Geometry と空間の依存関係の見直しを全体の対象とする。

必須課題そのものは既存の合意事項である。係数型は3.1節のとおり移行済みであり、LocalStructure の廃止方針は3.2節で採用した。その他の着手順序と設計候補は提案であり、State と Operator の型の統合、HilbertSpace・Chain の内部表現の詳細は未決である。

結論を得た項目には、採否、理由、公開 API と実装への影響、移行方法、検証結果を追記する。確定した内容は対応する要件・設計文書・実装計画へ反映し、反映後に本書を整理または削除する。

## 2. 棚卸し時点の実装（係数移行前）

Step 4 までの機能を実装済みである。Partial trace、式全体の基底展開・複合基底・座標取得・写像の表現係数、数値 export は後続段階にある。今回の見直しだけを理由に、それらの機能実装を追加しない。

| 対象 | 現在の構成 | 主な実装 |
| --- | --- | --- |
| 数値 | 係数は `num_complex::Complex<T>`。`Zero`・`One` と、平方根を使う操作の `Float` は num-traits を利用する。両ライブラリを core から再公開していた | [core の公開 API](../crates/qmb-core/src/lib.rs)、[係数演算](../crates/qmb-core/src/coefficient.rs)、[正規化の移行先](../crates/qmb-normalize/src/lib.rs) |
| 局所構造 | `LocalStructure` が状態・演算子 primitive の関連型と固定次元を定義する。恒等元は `Identity<T>`、他の演算能力は primitive 側の trait に分離する | [局所構造](../crates/qmb-core/src/primitive.rs)、[演算能力](../crates/qmb-core/src/ops.rs) |
| 局所式 | `LocalState<L, T>`・`LocalOperator<L, T>` が別々の節点を `Arc` で共有する。和・スカラー倍などの内容を共通パーツで保持する | [局所式](../crates/qmb-core/src/expr/local.rs)、[共通パーツ](../crates/qmb-core/src/expr/parts.rs) |
| 複合式 | `State<I, L, T>` は常に所属空間を持つ。`Operator<I, L, T>` は外側のハンドルに任意の束縛先を持ち、子節点に旧束縛を残さない | [複合式](../crates/qmb-core/src/expr/composite.rs) |
| Hilbert 空間 | `HilbertSpace<I, L>` は Factor(Index) と順序付きの入れ子 Tensor を共有する。`L` の値は保持せず、型によって局所定義を関連付ける。実装は ordinary のみ | [空間](../crates/qmb-core/src/space.rs) |
| Geometry・Chain | Geometry が Site 順の空間列と全空間を持ち、Chain が Geometry と境界条件を持つ。両型が `L` に依存する | [幾何](../crates/qmb-geometry/src/lib.rs)、[モデル utility](../crates/qmb-spin-half/src/utility.rs) |

状態への作用、内積、full trace は目的別の評価経路を持つ。通常の評価に基底展開や数値行列を要求せず、状態への作用評価には primitive の Product 能力を要求しない。

## 3. 5項目の整理

### 3.1. core の係数型・数値能力の整理

**目的**：core に num-complex を含めたまま、係数型の能力契約を最小限に整理する。2026-10-25 の利用者の指示により、両数値ライブラリへの依存解消を目標とした旧方針から変更した。

**移行前**：複素数型は公開 API、線型結合、式節点、局所規則の返却値まで使われている。spin-half も core の再公開を利用する。複素数の基本演算には独自の補助関数があり、num-complex の乗算が要求する Num より狭い制約で処理している。core の任意の `libm` feature は `num-traits/libm` を有効にする。

**合意した基本方針**：係数型には「必要な時に必要な trait bound を課す」。保持・構築のために全数値能力を要求せず、各操作で実際に使う能力だけを要求する。全能力をまとめた係数 trait を共通の必須条件にしない。

零・単位元には num-traits の Zero・One をそのまま採用する。Zero の Add・ゼロ判定、One の Mul は付随する条件として受け入れる。厳密に分割した独自 Zero・One は作らず、ゼロ判定も Zero::is_zero を使う。

**合意した接続方針**：core は num-complex・num-traits の直接依存を維持する。不足する trait は core 自身が定義し、Complex<f32>・Complex<f64> へ直接実装する。標準の複素数接続のために newtype や演算 provider を必須にしない。

**追加の採用事項**：

- 成分型 T から係数全体 C を保持する API へ移行する。
- 共役・実数値取得・実係数の埋め込み・純虚係数の生成を個別の能力として扱う。
- norm・normalized・正規化期待値は std を利用できる qmb-normalize の拡張 trait で提供し、NormalizationError も同 crate に置く。実数型 R に Float を要求し、provider・独自 Sqrt trait は導入しない。
- 標準 Complex の左スカラー倍は個別実装し、任意の C は scaled(c) を基本経路とする。
- 正規化期待値では純実分母を不要に二乗しない既存の数値的な処理を維持する。

core に内積・期待値・norm squared・スカラー除算を残し、正規化期待値と零ノルムの検査・NormalizationError は qmb-normalize に移す。core/expr の正規化 module は削除した。qmb から新 crate と拡張 trait を再公開する。実装は Conjugate・RealPart・ImaginaryUnit・Divide、From<R>、Norm／Normalize に具体化した。数値ライブラリの再公開は削除し、必要な箇所では num-complex・num-traits を直接 import する。

この新 crate は std の利用を許容し、no_std 対応を要件にしない。core 自体の no_std + alloc は維持する。

**実施した変更範囲**：数値の契約、係数表現、具体型との接続、およびその利用箇所に限定した。LocalStructure、空間の構造、State／Operator の分離は維持した。

**完了条件**：操作ごとの必要能力と具体的な Complex の接続を実装し、保持・遅延構築へ不要な数値能力を要求しないことを確認する。単精度・倍精度、独自 scalar、数値能力のない型での保持・遅延構築、平方根を使う／使わない操作を検証する。依存と feature は実際の利用に合わせて確認する。

#### 3.1.1. 完了確認（2026-10-26）

コミット `3d6d2af` の実装・テストを採用方針と照合し、本項の完了条件を満たすことを確認した。確認した範囲で完了を妨げる問題はなく、係数型のリファクタリングを完了とする。

- 線型結合・式節点・演算 trait・局所／複合評価・spin-half・モデル utility・facade が係数全体 C を使う。保持・map・共有・係数計算を伴わない遅延構築に数値能力を要求せず、評価時に個別の能力を要求する。
- 標準 Complex への直接実装、num-traits の Zero／One、左スカラー倍、任意の C の scaled(c)、qmb-normalize の拡張 trait と feature の接続が採用方針に沿っている。純実分母を二乗しない除算と、自己内積全体による零ノルム検査も維持している。
- 両精度・独自係数・数値能力のない係数・最小限の型制約・非正規直交な三次元代数・独立行列との一致・深い式の共有／評価／破棄を既存テストで確認した。

Rust 1.98.0 で次を再実行し、すべて成功した。

- `cargo test --workspace --all-features -j 2`（doc test を含む）
- `cargo test --workspace --no-default-features -j 2`（doc test を含む）
- `cargo fmt --all -- --check`
- `cargo clippy --workspace --all-targets --all-features -j 2 -- -D warnings`
- `RUSTDOCFLAGS='-D warnings' cargo doc --workspace --all-features --no-deps -j 2`
- `cargo check -p qmb-core --no-default-features -j 2`
- `cargo run -p qmb --example heisenberg -j 2` と同コマンドの `--no-default-features` 構成

`cargo tree -p qmb-core --no-default-features -e normal,features` で、core 単独の通常依存に std・libm・qmb-normalize がないことを確認した。facade の default feature 無効時にも qmb-normalize を含まない。実行例は両構成で期待値 `1.25+0i`、既定構成で norm `1` を返した。

MSRV 1.85 と専用 no_std target の実ビルド、性能ベンチマークは未確認である。core の no_std はソース・依存 feature・ホスト上の単独ビルドによる確認である。正式文書への反映は[実装計画の第7節](implementation-plan.md#7-係数型のリファクタリング)を参照する。後続項目の採否は各節に記録する。

### 3.2. LocalStructure の必要性

**目的**：維持・変更・廃止のいずれかを決め、理由を記録する。

**移行前**：関連型 State・Operator と `DIMENSION` を提供し、一つの型 `L` が一つの固定定義・固定パラメータを表す。異なる局所定義の接続は型で拒否する。式・空間・Geometry にこの型が渡るが、局所構造の実行時の値を共有・照合する仕組みではない。

**採用した結論（2026-10-26）**：LocalStructure を廃止し、追加の定義ラベルを設けない。primitive の意味・規則は S・O の型と各能力の実装で定める。詳細な依存調査と採用理由は[別紙](refactoring-local-structure.md)、数学的契約は[局所構造の設計第7節](design/local-structure.md#7-localstructure-廃止の採用方針2026-10-26)に記録する。

| 責任 | 採用した配置 |
| --- | --- |
| primitive の型対応 | 演算子式に O、状態式に S・O を直接指定する。Act 等の対応は評価時に要求する |
| 次元 | primitive 型全体が張る線型空間の次元を提供する。現在の契約では S が d、O が d² |
| 恒等元 | 演算子 primitive O の値として、係数型に依存せず生成する |
| 定義の区別 | S・O と能力の実装が担う。同じ S・O を別の L や定義 ID で区別しない |
| 因子構造・幾何 | 空間・Geometry／Chain から L を外す。空間比較・bind は Index・順序等を確認する |

**理由**：S・O を直接指定すれば、保持・遅延構築に演算能力を要求せず、型の対応を利用箇所で要求できる。異種の式の混在拒否も式側の型制約で維持できる。現在も積・作用・内積等は S・O・C で決まり、L はそれらの規則を切り替える引数ではない。次元と恒等元を移した後に名目的なラベルを残す必要はない。

**完了（2026-10-26）**：全6 crate の製品コード・テスト・利用例を移行した。Identity・SpanDimension と次元取得 API を実装し、記号的な Identity 節点は遅延構築と単位元の直接処理のため維持した。[実装記録と検証結果](refactoring-local-structure.md#8-実装した-api-と移行結果2026-10-26)に公開 API の移行、108／107テスト、Clippy・format・ドキュメント生成・core 単独ビルド・実行例の結果を記録する。3.3〜3.5の追加の設計判断は後続とする。

### 3.3. Operator と State の整理

**目的**：操作的量子論の枠組みで両者を一般の演算子の特別な形として扱えるかを検討し、数学的意味と API の対応を決める。

**現状**：State は ket であり、`ℂ→H` として始域・終域を参照できる。束縛済み Operator は `H→H` として参照できる。ただし、局所・複合それぞれに状態用と演算子用の式木・評価処理があり、一般の `U→V` の写像を表す型はない。

**決めること**：

- 「一般の演算子」が表す数学的対象と、今回扱う状態の意味。Ket、密度演算子、量子操作を同じ対象として扱う前提を置かない。
- 合成、tensor product、随伴、内積、trace、zero・identity を始域・終域とどう対応付けるか。
- 数学的な共通化を公開型、内部表現、共通処理のどこまで反映するか。数学的に整理できることと、公開型を統合することは別に判断する。
- 未束縛の配置式と、始域・終域が確定した写像の関係をどう表すか。

**依存関係**：数学的な整理は3.2・3.4の設計前または同時に行う。型と評価器の統合を伴う実装は、責任分担が定まってから範囲を決める。

**完了条件**：採用する数学的対象、対応する API、共通化する範囲、別に維持する範囲を理由とともに記録する。統合を採用しない結論も含む。

### 3.4. HilbertSpace の最小構成

**目的**：Index と Grading だけで空間を議論・表現できるかを決め、追加情報の責任の置き場所を明らかにする。

**現状**：LocalStructure 廃止により空間は Index と子空間を保持し、primitive 型は式側が持つ。次元は取得時の primitive 型または明示的な局所次元から得る。ordinary の互換性は葉の Index 列、部分系判定は順序を保存する部分列の照合で行う。Grading を含む最小構成の詳細は未決である。

**決めること**：

- 因子の同一性・順序・grouping・積の種類を表す構造と、Hilbert 空間の次元・局所定義を分離できるか。
- Index と Grading だけの構造が提供する能力と、別の情報を必要とする能力は何か。
- 次元取得、配置の整合性検証、恒等的に作用する因子の trace、将来の基底処理に必要な情報をどこから得るか。
- Grading が表す情報の範囲。Tensor の積の種類と、状態・演算子成分の parity の責任を区別する。
- 木の共有、空・単因子の同一視、互換性と保持する grouping の区別をどう維持するか。

**依存関係**：3.2 と一体で設計し、結果を3.5へ渡す。表現の検討に際して graded 演算の実装完了までは要求しない。

**完了条件**：最小構成の採否、空間だけで提供する能力、外部の情報を要する能力とその受け渡し方法を記録する。情報を外へ移した場合も、整合性検証の責任を残す。

### 3.5. Chain の内部表現

**目的**：3.4で検討する構成から Grading を除いた構造を、Chain の内部表現として利用できるかを決める。

**現状**：Chain と Geometry から局所構造型への依存は除去済みであり、Index 型だけを指定して HilbertSpace の因子構造を保持する。近接関係の列挙は Site 順序・境界条件・葉の Index を使う。一方、spin-half の utility は Factor への配置と全空間への束縛に現在の空間 API を利用する。

**決めること**：

- Site・Index・順序・grouping・近接関係に必要な構造を、局所代数と空間の数学的情報から分離できるか。
- 3.4の構造を共用するか、共通の内部部品を持つか、別の表現にするか。
- Site ごとの配置と全体の因子順序をどう保持するか。Site 順序と明示された全空間の因子順序は異なり得る。
- Hilbert 空間の構築・参照、局所式の配置・束縛、モデル utility との接続をどこへ置くか。
- 既存 Index の再利用、新規発行、重複・未配置の検証、空 Site をどう扱うか。

**依存関係**：3.4の結論を受けて設計・実装する。Chain だけを先に独自の Index コンテナへ変更し、後から再移行する順序は推奨しない。

**完了条件**：構成の採否と理由、Geometry・HilbertSpace との境界、モデル構築までの利用方法を記録する。変更する場合は既存の配置・近接関係・モデル utility の契約を検証する。

## 4. 着手順序の提案

まず数値部分に着手する。利用箇所は広いが、空間・局所構造・状態の数学的な整理から比較的独立しているためである。これは公開 API の変更が小さいという意味ではない。係数型を抽象化する場合は型引数や演算子構文にも影響し得るため、候補の比較時に変更範囲を確認する。

| 順序 | 作業 | 成果物 |
| --- | --- | --- |
| 1 | 数値能力と利用箇所の棚卸し、係数表現・具体型の配置の比較 | 3.1の設計判断と移行範囲 |
| 2 | 係数の能力契約と Complex の具体実装、利用側の移行 | 空間・式の構造を維持した実装と検証結果 |
| 3 | State／Operator の数学的整理と、LocalStructure／HilbertSpace の責任分担の検討 | 3.2〜3.4の採否・根拠・API 案 |
| 4 | 決定した責任分担に沿った移行、Chain の表現の検討 | 変更単位ごとの実装・検証、3.5の設計判断 |
| 5 | Chain とモデル構築の接続の移行、採用した式の共通化の実装 | 縦断的な検証と正式文書への反映 |

後半の実装順序は3.2〜3.4の結論に応じて更新する。設計をまとめて検討することと、複数の責任を一度の変更で移行することは区別する。

## 5. 移行時の検証

各変更は、次の既存契約を基準に確認する。契約自体を変更する結論を採る場合は、変更理由と代わりの保証を設計文書に明記する。

- core の `no_std + alloc`、具体的代数・幾何・UUID・数値 backend からの独立性。
- 操作ごとの最小限の型制約。保持・map・遅延構築へ不要な数値能力や primitive の Clone を要求しない。
- 行列化・基底展開を必須にしない評価、入力の不変性、不変実体の共有、深い式・空間の反復的な走査・評価・破棄。
- Index の同一性、因子順序と grouping、重複の拒否、非連続な部分系の判定、局所定義の整合性。
- Zero・identity の空間保持、空の tensor product、明示空間の優先、配置と束縛の区別。
- 単精度・倍精度、非正規直交基底、spin と独立した三次元代数、独立な小規模行列との一致。
- Chain の空・小規模・open／periodic、複数自由度・空 Site、Index の発行と再利用、XYZ・Heisenberg・Neel の縦断例。

実装変更時は対象テストから確認し、workspace の test、format check、Clippy と必要なドキュメント生成へ進む。数値契約の移行では依存グラフと feature 構成も確認する。

事前確認として、2026-10-25 に以下のテストを実行し、doc test を含め成功している。これは移行前の確認であり、リファクタリング後の検証結果ではない。

- `cargo test --workspace --all-features`
- `cargo test --workspace --no-default-features`

## 6. 状態と記録上の注意

| 項目 | 現在の状態 |
| --- | --- |
| 1. 係数型・数値能力 | **完了（2026-10-26）**。実装を採用方針と照合し、完了条件と検証の成功を確認。範囲・結果・未検証事項は3.1.1節に記録 |
| 2. LocalStructure | **完了（2026-10-26）**。S・O の直接指定、Identity・SpanDimension、空間・幾何からの L 除去を実装・検証。追加の定義ラベルは設けない |
| 3. Operator／State | 現状と論点を整理済み。数学的対象・共通化範囲は未決 |
| 4. HilbertSpace | 現状と論点を整理済み。最小構成と責任分担は未決 |
| 5. Chain | 現状と論点を整理済み。内部表現と接続方法は未決 |
| 着手順序 | 順序1〜2の係数型の整理・移行は完了し、実装計画の第7節へ反映済み。順序3以降は本書の提案 |

参照されている `docs/step-3-design.md` は、現時点の作業ツリーには存在しない。現状の把握には、現行コード・テストと実装計画の段階4の実装記録を用いた。実装計画の冒頭に残っていた古い進捗・crate 数は、係数型の採用方針の反映時に更新した。
