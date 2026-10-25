# Step 3: 空間・配置・複合評価の実装方針

更新日: 2026-10-23

## 1. 位置付けと実装範囲

本書は [実装計画](implementation-plan.md)の Step 3 で採用した設計判断と実装結果をまとめる。Step 1〜3 は実装済みである。[要件](requirements.md)、[Index](design/index.md)、[局所構造](design/local-structure.md)、[Hilbert 空間](design/hilbert-space.md)、[複合式](design/composite-expression.md)の数学的契約に従う。

Step 3 は同じ局所構造の ordinary な複合系を対象とし、Geometry なしで空間・状態・演算子を構築・評価できるところまで実装する。Index の抽象化と Dof の葉への統合、およびその他の推奨案を採用した記録として、以下の方針を確定事項とする。メソッド名・型引数の順序・非公開の補助型などは実装時に整える。

| 対象 | Step 3 の提供範囲 |
| --- | --- |
| `qmb-core` の空間 | ジェネリックな Index の照合、Factor・入れ子 Tensor、構成要素・Index の列挙、互換性・部分系判定、次元の取得 |
| `qmb-core` の複合式 | 配置、加算・係数乗算・zero・identity、状態の tensor product、演算子積、束縛、随伴、遅延作用と評価、内積、期待値、full trace |
| 導出操作 | Commutator・anticommutator、Hilbert–Schmidt 内積 |
| `qmb-basic` | UUID による具体的な Index 型と seed 指定不要の発行 API |

Norm・正規化・bra・正規化期待値は、局所式と複合式に共通する別の作業として Step 4 の受入確認までに設計・実装する。Partial trace・基底展開・座標取得は Step 5 とする。これらは初期版の必須機能として維持する。

数値制約は実際に使う操作に限定する。特に norm の平方根のために、基本的な式の構築・作用・内積・trace へ一律な `Float` 制約や数学関数 backend の依存を追加しない。

## 2. Index の抽象化と具体実装

core の空間・配置済み状態・演算子は Index 型 `I` をジェネリックに受け取る。同じ Index であることを照合でき、同一性が安定した等価関係であることを契約とする。Rust では `Eq` 相当とし、UUID の表現・発行・乱数源・seed・大小関係・`Hash`・`Ord` を一律に要求しない。保持・共有に必要な制約は実際に使う箇所に限定する。

具体的な表現と一意性の確保は Index の提供側の責任とする。標準の具体型と発行 API は `qmb-basic` に置き、UUID を用い、利用者の seed 指定を要求しない。core は `qmb-basic`・UUID・乱数源に依存せず、`no_std + alloc` を維持する。Index の独立した発行と、空間の構築・複製・配置で暗黙に再発行しない契約を守る。

Geometry の発行補助も外部の発行処理を利用し、UUID の具象化を必須依存にしない。core のテストは UUID に依存しない独自 Index 型で実施し、basic の発行テストと分離する。

## 3. 局所構造の同一性

初期版では一つの型 `L: LocalStructure` が一つの定義と固定パラメータの組を表す。現在の固定次元・静的な primitive の関連付けを維持し、同じ `L` の空間・式だけを接続する。型の実装者は、その型の内部で局所構造の定義や表現パラメータを変えない契約を負う。別の `L` との接続は型で拒否する。

これは [局所構造の数学的な一致条件](design/local-structure.md#4-局所構造の一致条件)を初期版の制約下で実現する方式である。独立した定義 ID、異なる実装型間の定義照合、実行時パラメータの保持・照合は導入しない。`TypeId` のための一律な `'static` 制約も追加しない。

## 4. HilbertSpace と基本自由度

公開の空間型は `HilbertSpace` に統一する。基本自由度は Index と局所構造を持つその葉（Factor）として表し、独立した `Dof` 型は設けない。Tensor は順序付きの子空間を持つ。構築時に保持した葉を、そのまま空間の合成・局所式の配置先指定・部分系判定に再利用する。

`HilbertSpace` は不変な `Arc` の共有ハンドルとし、Factor・Tensor の内部節点は非公開にする。利用者が指定した入れ子を保持し、複製・合成で変更されない子空間を共有する。共有ハンドルの複製に局所構造や Index 本体の不要な複製を要求しない。

`children()` は保持している直下の子空間を構築順に、`indices()` は全葉の Index を左から順に遅延列挙する読み取り API とする。借用を基本とし、列挙のための Index の再発行や式・基底の展開を行わない。Factor の直下の子はなく、所属 Index は一つとする。

互換性判定は、保持している木を変更せず、ordinary な葉の反復走査によって Index・順序・局所構造を照合する。空・単因子の同一視と ordinary の結合則を適用し、括り方に依存させない。部分系判定は同じ葉列の部分列照合とする。Index の値で因子を並べ替えず、次元やポインタ一致だけでは判定しない。

Tensor の構築時には全階層の葉に対する Index の重複を検証する。空間の走査と最後の所有者になった節点の破棄にも明示的な作業スタックを使い、深さに比例した関数再帰を避ける。Graded な構成は今回実装しない。

## 5. 次元と overflow

総次元が `usize` に収まることを空間構築の条件にしない。各葉の有限次元と因子構成を保持できれば空間を構築でき、互換性判定・作用・内積・trace に総次元の整数化を必須としない。

次元を整数として取得する API は `Result<usize, DimensionOverflow>` 相当とし、子の次元の積を checked arithmetic で求め、収まらなければエラーを返す。空の Tensor の次元は `1` とする。数値 export などが要素数を要求する段階でも overflow を検証し、暗黙の切り詰めは行わない。

Full trace における spectator の寄与は局所恒等演算子の trace の積から処理できるようにし、総次元を `usize` にしてから係数型へ変換することを必須にしない。係数の丸め・overflow・特殊値は既存の数値規約に従い、空間次元の整数 overflow と区別する。

## 6. 複合状態と空間の保持

複合状態型 `State` は、Index 型・局所構造・係数精度を型引数で関連付け、常に確定した所属空間を共有ハンドルとして持つ。式本体は `Arc` で不変に共有する。局所状態式の配置で葉の空間に属する状態を作り、tensor product で空間と状態を指定順序に合成する。葉に保持する局所式は単一 primitive に限定しない。

状態の加算は互換な所属空間上で行い、tensor product は因子の重複を拒否する。局所式の配置先には一つの基本因子を指定していることを検証する。状態への因子の追加はその因子の状態との tensor product で行い、identity spectator では補わない。

Zero state は `State::zero(&space)` 相当として所属空間を明示する。空の和にも空間を要求する。空の tensor product は0因子空間上の状態 `1` とする。構築時に数値能力を広げないため、単位状態を暗黙の `1` として保持できる内部表現を用いる。

評価で zero になっても所属空間を残す。状態式の括り方が異なっていても、互換な空間上の状態には内積・加算を適用できる。状態は所属空間を参照でき、写像としての始域は0因子の単位空間、終域は所属空間とする。共通の汎用写像 trait は今回追加しない。

## 7. 複合演算子と束縛

複合演算子型 `Operator` は、Index 型・局所構造・係数精度を型引数で関連付け、不変な式本体と任意の束縛先空間を持つ。一つの型で未束縛と束縛済みを扱い、束縛情報は外側のハンドルに集約する。葉は局所演算子式全体と配置先 Index を保持し、演算子用の tensor product 節点は設けない。

演算子の和・積を構築する際は以下に従う。

- 全 operand が束縛済みなら空間の互換性を検証し、結果にも束縛先を保持する。
- 一つでも未束縛なら作用空間を比較せず、結果を未束縛とする。
- 子節点には operand の旧束縛を取り込まず、式本体を共有する。元の operand は変更しない。
- 積の構築から全作用空間を自動合成せず、参照 Index だけから spectator を含む全空間を推論しない。

束縛や評価で用いる空間の選択は以下の順序で一か所に実装する。

| 条件 | 使用する空間 |
| --- | --- |
| 明示空間あり | 指定先を優先し、旧束縛先との一致は要求しない |
| 明示空間なし・束縛済み | 既存の束縛先 |
| 明示空間なし・未束縛の作用 | 状態の所属空間 |
| 明示空間なし・未束縛の trace | 空間未確定のエラー |

選択した空間に全参照 Index が存在することを検証する。作用では状態空間との互換性も検証し、状態の所属空間を暗黙に拡張しない。束縛・作用で式を返す場合は選択した空間を結果に保持する。束縛済み演算子は同じ空間を始域・終域として参照できる。

演算子の空和・空積は未束縛の zero・identity として構築でき、束縛後はその空間を保持する。係数だけの演算子は選択した空間上の scalar identity として扱う。評価後の zero・identity でも束縛先を失わない。

## 8. エラーと遅延評価の境界

名前付きメソッドを基本 API とする。空間構築、状態の加算・tensor product、演算子の加算・積・束縛・作用、内積・期待値・trace など、整合性検証が必要な操作は `Result` を返す。係数乗算や随伴など検証が不要な遅延構築は、式を直接返す。

複合式の `+`・`*` など検証が必要な演算子構文も `Output = Result<式, Error>` とし、不整合を panic にしない。局所式の infallible な API は維持する。式同士の二項演算子は所有値・借用値の4通りを提供し、名前付きメソッドは共有借用を基本とする。

空間構築のエラーと式の接続・評価のエラーを区別し、少なくとも重複 Index、局所配置先の不正、参照 Index の不存在、空間不整合、作用空間未確定、および次元 overflow を判別できるようにする。局所構造・精度の不整合は今回の型構成で拒否できるものを型で拒否する。具体的なエラー型・variant の構成は実装時に定める。

構築時に分かる空間の不整合は早期に検出する。状態の `act` は選択した空間を検証して遅延作用を構築し、`evaluate` が局所規則を実行する。明示的な zero・identity や零係数の項でも、配置された Index の不正を計算の省略によって見逃さない。局所規則は既存どおり infallible とする。

## 9. 評価器と中間表現

Step 2 の共通パーツ、局所式、局所規則を再利用し、複合系専用の primitive や局所規則を追加しない。状態評価・演算子の作用・trace に必要な項走査を分け、全経路に共通の完全展開形を強制しない。

| 経路 | 処理方針 |
| --- | --- |
| 作用 | `ABψ` を `B`、`A` の順に処理し、対象因子の局所状態式へ作用する。状態への作用評価に primitive の `Product` 能力を要求しない |
| 内積 | 線型性・共役線型性で必要な項をたどり、積状態間では対応する局所内積の積を用いる。入れ子の括り方を照合条件にしない |
| Full trace | 必要な演算子項を順次走査し、同一因子上の積順序を保って局所 trace の積を求める。局所積が必要な経路で `Product` を要求する |
| 随伴 | 遅延節点で保持し、評価時に係数の共役と積順序の反転を処理する |

必要な項・因子だけを扱い、総 Hilbert 次元に比例する密な中間領域や全大域演算子の完全展開結果を一括確保しない。入力と再利用できる部分式・空間を共有し、結果を通常の式として返す。評価・走査・最後の所有者となった節点の破棄は反復的に行う。

特定の評価済み表現、恒久的な support metadata、キャッシュ、同類項の自動結合、一般の入力に対する最適な計算量は要求しない。Step 2 で確認した局所作用チェーンの係数節点の再走査・再構築の最適化も、今回の完了条件に追加しない。

## 10. 実装順序と検証条件

| 順序 | 実装・検証する内容 |
| --- | --- |
| 1 | core の Index のジェネリック化と空間の葉・Tensor。独自 Index 型、共有、列挙、重複検出、空・単因子・入れ子、互換性・部分系判定、次元 overflow |
| 2 | 複合状態・演算子の構築、配置、束縛、空間選択、エラーと演算子構文。構築時に局所規則を呼ばないこと |
| 3 | 作用・評価・内積・期待値・随伴・full trace と導出操作。spin と独立した三次元代数のそれぞれで単精度・倍精度の行列参照計算と照合 |
| 4 | `qmb-basic` の具体型・発行 API と core への接続。seed 指定不要の発行、複製での同一性保持、core の通常依存への UUID・乱数源の混入がないこと |

数学的結果に加え、以下を受入条件とする。

- 非連続な部分系、並べ替えの拒否、一次元 Factor と空の Tensor の区別、括り方の違いによらない演算結果。
- 束縛済み同士の不整合、未束縛を含む合成、明示空間による旧束縛の上書き、状態空間の暗黙拡張の拒否、未束縛 trace の拒否。
- 異なる因子への演算子和と積、同じ因子上の積順序、恒等作用、spectator を含む full trace、`Tr(I_∅) = 1`。
- 0因子の状態 `1`、所属空間を持つ zero state、評価後の zero・identity の空間保持、不正な参照を含む zero・identity の拒否。
- 総次元が `usize` を超える空間を構築でき、次元取得だけが overflow を返し、基本的な行列を作らない作用・内積・trace が利用できること。
- 数値能力を持たない型での構築、不要な primitive・Index・局所構造の複製の抑制、共有先が残る場合の保持、長い式と深い空間の評価・走査・破棄。

対象テストの後、workspace の test・format check・Clippy・警告をエラーとするドキュメント生成を行う。新しい依存を導入する段階で、core の `no_std + alloc` と MSRV への影響を確認する。環境上実行できない確認は未確認として記録する。

## 11. 実装した API と後続事項

Step 3 の公開型の役割、Index の抽象化、Dof の統合、局所構造の静的な同一性、`Arc` による共有、次元 overflow の検出時点、空間の保持と選択、fallible な複合操作、反復的な評価・破棄は採用済みである。

Step 3 では以下を実装した。型の保持・共有に不要な `Clone`・`Hash`・`Ord`・数値能力は要求しない。Index の照合には専用 trait を追加せず、必要なメソッドで標準の `Eq` を要求する。

| 公開型・対象 | 実装した API |
| --- | --- |
| `HilbertSpace<I, L>` | `factor(index)`、`unit()`、`tensor(iter)`、`tensor_product(&rhs)`、`children()`、`indices()`、`factor_index()`、`dimension()`、`is_compatible_with(&rhs)`、`is_subsystem_of(&ambient)`、`contains_index(&index)`。`shares_node_with` は共有実体の確認用 |
| 局所式の配置 | `LocalState::at(&factor)`・`LocalOperator::at(&factor)`。式全体を共有して葉へ配置する |
| `State<I, L, T>` | `place`、`zero(&space)`、`unit()`、`sum(&space, iter)`、`add`、`tensor`・`tensor_product`、`scaled`、`space`・`domain`・`codomain`、`act`・`act_in`、`evaluate`、`inner_product`、`expectation`・`expectation_in` |
| `Operator<I, L, T>` | `place`、`zero()`、`identity()`、`scalar`、`sum`、`add`、`product_all`・`product`、`scaled`、`adjoint`、`bind`、`space`・`domain`・`codomain`、`commutator`・`anticommutator`、`trace`・`trace_in`、`hilbert_schmidt_inner_product` とその `_in` 版 |
| エラー | `SpaceError::DuplicateIndex`、`ExpressionError::{NotAFactor, IncompatibleSpaces, MissingIndex, UnboundOperator, Space}`、次元取得の `DimensionOverflow` |
| `qmb_basic::Index` | UUID v4 を包む不透明な具体型。`new()` と `Default` は seed 指定不要の新規発行、`Clone`・`Copy` は同一性を保持する。`Display` と比較・hash にも対応する。`qmb_basic::HilbertSpace<L>` はこの Index を使う別名 |

名前付きの二項メソッドは共有借用を基本とし、式同士の `+`・`-`・演算子積・作用の演算子構文は所有値・借用値の4通りを実装した。検証を伴う演算子構文は `Result` を返す。単項負号と複素係数による左右の乗算も提供する。

UUID の依存指定は `uuid = { version = "1.18.1", features = ["v4"] }` とし、検証環境では `uuid 1.26.1`・`getrandom 0.4.3` が解決された。UUID の発行は OS の乱数源を使い、乱数取得失敗時の挙動は `uuid::Uuid::new_v4` に従う。UUID・乱数源は core の通常依存に含めない。

Norm・正規化・bra・正規化期待値の具体的な数値 API は Step 4 の受入確認までに決める。Partial trace・基底展開・座標取得は Step 5、数値 backend は Step 6 で具体化する。異種局所構造、実行時の定義パラメータ、異なる実装型間の定義照合、graded、保存・復元、キャッシュ、並列化などは既存の延期方針に従う。

## 12. 実装・検証状況

2026-10-23 に Step 3 の実装を確認した。

- Spin-half の三因子系と、core の公開 API だけで定義した三次元代数の二因子系を、単精度・倍精度の独立した Kronecker 行列計算と照合した。積・順次作用・随伴・内積・期待値・full trace・Hilbert–Schmidt 内積、非正規直交基底から取り込んだ局所式を含む。
- UUID に依存しない `Eq` だけの Index、非 Clone な primitive と借用、数値能力のない型での構築、Product を持たない演算子の作用、Num・Float・除算・剰余・Copy を要求しない係数計算を検証した。空間・精度の混在と非 Send の Index は compile-fail テストで確認した。
- 所有値・借用値の演算子構文、入力・spectator の共有、束縛の合成と明示空間による上書き、zero・identity の空間と参照の検証を確認した。総次元が `usize` を超える系でも、作用・内積・期待値・trace が行列を作らず利用できる。
- 2万段の入れ子空間、状態の tensor・和・係数乗算・作用、演算子の積・随伴と trace を128 KiB のスタック上で走査・評価・破棄した。
- 既存の direnv 環境（Rust 1.98.0）で全70テスト（doc test を含む）、format check、Clippy、警告をエラーとするドキュメント生成を確認した。core の依存に UUID・乱数源・`std`・`libm` が入らないことも確認した。
- 解決された UUID・getrandom の package metadata は Rust 1.85 を要求する。宣言した MSRV（Rust 1.85）での実行と `std` のない target 向けビルドは、対応 toolchain・target がないため未確認である。
