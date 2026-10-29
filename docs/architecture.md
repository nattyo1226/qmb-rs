# qmb-rs アーキテクチャ概要

更新日: 2026-10-30

## 1. 文書の位置づけ

本書は、[要件](requirements.md)に基づく設計の概要と、合意済みの設計判断を記録する。各対象の詳細な契約は `docs/design/` 以下の文書に記録する。具体的な Rust の型、公開 API、および内部実装は、対応する設計を決めた時点で追記する。

2026-10-28 に QuantumSystem と SpinHalf<C> の構築 API を導入し、Chain をサイト数・境界条件だけの幾何へ移行した。本書は移行後の構成を表す。

## 2. 責任分離

| 対象 | 担当すること |
| --- | --- |
| 局所構造 | 局所空間の情報、代数・表現、primitive と数学的能力の関連付け・整合性 |
| 状態・演算子 primitive | 局所的な内積・積・作用・随伴・trace、および対応する場合の parity |
| Index | 局所構造から独立した因子の識別子の発行と同一性 |
| 自由度・Hilbert 空間 | Index と局所構造の組、Factor・Tensor の再帰構造、構成要素・Index の列挙、因子順序、積の種類、grouping・互換性・部分系判定、複合系の数学的能力の合成と整合性 |
| 状態式・演算子式 | core の共通構文、線型性・共役線型性による拡張、不変な部分構造の共有 |
| 埋め込み・評価 | 局所式の配置と作用空間の確定、局所契約を保存する合成則、整合性検証、grading、必要な部分の評価 |
| 基底・export | 局所基底と tensor product 基底、通常の式木への展開・再展開、順序・列挙・要素数の導出、座標取得、外部の数値表現への変換。写像の表現行列は future work |
| 幾何 | サイト数・境界条件、外部の Site 順列への関係の借用列挙 |
| モデル構築・利用側 | Index と Factor の対応、空間の構築・保持・因子順序・grouping、局所式の配置・束縛、状態や Hamiltonian の構築 |

代数 core は幾何構造を必須とせず、自由度と式を直接構築できる。完全な密・疎行列を生成せずに代数オブジェクトを構築・合成・評価し、必要に応じて数値 backend に export する。

### 2.1. crate の構成

Chain はサイト数と境界条件だけを保持する。`qmb-geometry` は外部の局所空間列への接続補助のため `qmb-core` に依存するが、Chain 自体は Index・Hilbert 空間を保持しない。core は geometry に依存しない。

代数系の共通処理は `qmb-core` に集約する。抽象契約・有限線型結合に加え、局所・複合式木、Index・自由度・Hilbert 空間、局所式の配置・束縛、評価、基底展開・座標取得を同じ crate が所有する。写像の表現係数は将来の同 crate の担当とする。本文書と詳細設計でいう「代数 core」は `qmb-core` を指す。

具体的代数は専用 crate とし、`qmb-spin-half`・`qmb-spin-one`・`qmb-fermion` は primitive・局所規則・具体的基底を実装し、core の共通型を直接利用できる。UUID による具体的 Index と発行は `qmb-basic`、Site・近接関係と外部列への借用接続は `qmb-geometry`、外部数値コンテナへの変換は `qmb-export`、Heisenberg・Neel の構築 utility は `qmb-spin-half`、標準 Hubbard・symmetric Hubbard の構築は `qmb-fermion`、norm・状態の正規化・正規化期待値は `qmb-normalize` が担当する。workspace は8つのライブラリ crate と [examples パッケージ](../examples/README.md)で構成する。basic・geometry・spin-half・spin-one・fermion・export・normalize は core に依存し、spin-half と fermion のモデル構築は geometry の公開契約を利用する。normalize と spin-one は std を利用できる。core と geometry は `no_std + alloc` を維持し、core は basic・UUID・乱数源・具体的代数・geometry・normalize・数値 backend に依存しない。

core 内では局所契約、式、空間、評価、基底処理を module で分離する。共通表現と基本操作・名前付きメソッドを同じ crate で設計し、具体的代数ごとに別の式木や evaluator を作らない。式の内容を保持する `ScaleExpr`・`SumExpr`・`ProductExpr`・`ActExpr`・`AdjointExpr` を共通パーツとし、局所式・複合式で再利用する。`LinearCombination<P, C>` は薄い有限線型結合型として維持し、共通パーツで保持する。局所規則・局所基底展開は線型結合を返し、式木側が `map` で要素を変換して共通 constructor で取り込む。Step 2 で core・spin-half・テストの演算 trait を関連型に揃え、式木との接続を実装した。将来の追加 crate は、独立した機能・依存関係・実行基盤の境界が必要になった段階で検討する。

共通パーツは `expr::parts` で保持内容も含めて公開し、式ハンドルは crate root から再公開する。節点 enum と共有の内部表現は非公開とする。項型 `Term<P, C>` は `ScaleExpr<P, C>` の別名に統合し、要素名を `operand` に揃える。演算 trait の関連型 `Term` は係数を含まない要素型として区別する。既存の線型結合の構築・走査 API を維持し、走査を共通パーツへ集約する。

この分割は、共通型と基本演算を所有する [tensory-core](https://github.com/todo-group/tensory/blob/c84f091e5c65f1408dee91db3730bdb8729e966e/tensory-core/src/lib.rs)と、公開テンソル表現の所有先を集約しつつ実行系・演算群を分ける [tenferro-rs の crate 設計](https://github.com/tensor4all/tenferro-rs/blob/f1d4d8d961ff6b937768111888ad2a493314892a/docs/architecture/tenferro-crates.md)を参考にした。qmb-rs では primitive の局所規則と共通合成則に基づく評価を core に含め、外部数値 backend との接続を分離する。

### 2.2. crate 内のモジュール構造

core は `ops` に `Dimension::dim()`・`Identity` と代数契約、`expr::parts` に有限線型結合と共通構文パーツ、`expr::traversal` に内部走査を置く。式の保持・構築は `expr`、評価・内積・trace は非公開のトップレベル `evaluate`、基底列挙・展開・座標取得は `basis` が担当する。旧 `primitive` とトップレベル `linear_combination` は廃止し、主要型・trait の crate root からの再公開は維持する。

geometry は `boundary`・`geometry::chain`・`relation::nearest_neighbor`・`connection` に分け、Chain の既存メソッドから委譲する。export は `output`・`metadata`・`args`・共通の `coordinates` と `vector`・`matrix`、storage 契約と Vec・ndarray・sprs の `backend` に分ける。spin-half は `models`、normalize は `norm`・`expectation`・`error` に責任をまとめる。これらの crate の公開窓口は引き続き crate root とし、basic の構成は維持する。

### 2.3. 数値 export（実装済み）

2026-10-27に Step 6 の数値 export を完了した。qmb-export の `Export<T, Args = ()>` と、数値データ・メタデータを持つ `Exported<D, M>` を共通の利用窓口にする。密表現は Vec と shape・格納順序のメタデータを出力し、ndarray の Array1／Array2 への直接出力・取得済み Vec 出力からの変換も提供する。密行列専用の数値コンテナは設けない。

疎表現は core の Coordinates から sprs の疎ベクトル・CSC 行列を構築する。行列の列は状態基底への作用と座標取得から生成し、密・疎で共有する。空間・基底・因子順序・係数型の解釈を保持し、通常の式木に backend や基底識別情報を追加しない。

qmb-export は既定 feature なしで Vec 出力を提供し、ndarray・sprs は個別の feature で有効化する。examples の任意の export feature は両方を有効にする。公開 API・所有権・エラー・依存指定・未実施の検証は[数値 export の設計](design/export.md)に定める。Step 6 の実装と Rust 1.98.0 での受入確認は完了した。MSRV の実ビルドは利用者の依頼により省略している。

### 2.4. Fermion（実装済み）

`qmb-fermion` は `no_std + alloc` で、二状態の単一 mode、I・a・a†・n・su(2) 規格化した MajoranaX/Y/Z・Parity の八つの primitive を提供する。状態基底は OccupationBasis、演算子基底は LadderBasis と MajoranaBasis とする。`Fermion<C>` は既存 QuantumSystem の構築入口を使い、局所規則に Float・係数の Clone / Copy を要求しない。

複数 mode の交換符号は primitive の parity と core の graded 評価に委ねる。モデル構築は `Fermion::hubbard`・`symmetric_hubbard` に、Chain と各 Site が `(up,down)` の2 mode を持つ graded 空間列を渡す。両モデルの係数 bound は `Clone + From<R> + Neg` だけであり、構築に局所 Product・Act を要求しない。規約・基底・検証記録は [Fermion の局所代数](design/local-structure.md#9-fermion-の局所代数)、モデル API は [Hubbard モデルの設計](design/geometry.md#6-fermion-の-hubbard-モデル)、利用例は [hubbard.rs](../examples/src/bin/hubbard.rs) に記録する。

### 2.5. Spin-one（実装済み）

`qmb-spin-one` は三状態と I・スピン5値・四重極9値の計15演算子 primitive を提供し、状態空間3次元・演算子空間9次元を表す。状態基底は SzBasis、演算子基底は CartesianBasis / LadderBasis。同じ基底の primitive 同士の積はその基底内に閉じる。std を許可し、作用は `Complex<R>`・`R: Float`、他の局所規則は操作ごとの小さい bound とする。規則は直接の match で実装し、全要素 even・P=I で共通の ordinary / graded 評価へ接続する。[Spin-one の局所代数](design/local-structure.md#10-spin-one-の局所代数)と[実装仕様・検証](plans/v2/spin-one-primitives.md)を参照。

## 3. Index・自由度・局所構造・空間

Index は識別子のみを表し、局所構造や次元を持たない。局所構造は、局所 Hilbert 空間と、その上で実現される代数・表現を定める。基本自由度は Index を持つ `HilbertSpace` の葉（Factor）として表し、独立した `Dof` 型は設けない。葉も複合空間も同じ空間型で保持・合成する。同じ局所構造を持つ複数の因子を Index で区別する。

初期版では一つの系の全基本因子を同じ primitive 型 S・O で表す一種類の局所代数に揃える。異種局所構造の混在・合成は後続へ延期し、そのための型消去を初期版の要件にしない。ユーザ定義の局所構造は、別の系で単独利用できる。同種因子の入れ子や一つの Site への複数自由度の配置は維持する。

局所状態・演算子に Index を付与することで、所属先・作用先の因子を識別する。元の局所構造は局所状態・演算子が保持する。この付与は、spectator を含む全 Hilbert 空間への束縛とは別の操作である。

局所 primitive の意味・能力と数学的な整合性は、[局所 primitive の設計](design/local-structure.md)に定める。S・O の型と能力実装で定義を表し、追加の定義ラベル・定義 ID は導入しない。異なる規則や表現は primitive の型と実装に反映する。型消去による異種代数の合成と実行時パラメータは後続へ延期する。

Hilbert 空間自体を、基本自由度の Index を持つ Factor、および ordinary / graded の種類と順序付きの子 Hilbert 空間を持つ Tensor の再帰構造として表す。Tensor の子に Tensor を含められる。Factor と Tensor は概念上の名前とする。基本因子の順序は子を左からたどった葉の Index 列で表し、Index の値の大小から決めない。基底は空間とは別に選ぶ。

Grouping は空間の木の括り方で表し、そのための独立した Group 型は追加しない。Chain はサイト数と境界条件だけを保持し、Site 順の局所空間列と grouping は利用側が管理する。Site 順序と別途存在する全空間の因子順序を区別し、既存状態・行列を暗黙に並べ替えない。内部ラベルは初期版では導入しない。局所複合代数の数学的能力は、子の局所規則と共通の合成・評価機構から処理する。

空間の互換性は、基本因子の Index・順序・局所構造と、空・単因子 tensor product の同一視および同種の積の結合則を適用した後の tensor 構造に基づく。保持する木構造と数学的な互換性を分け、括り方だけの違いでは互換性を変えない。構築時の平坦化は要求せず、不変な子空間を共有できるようにする。

部分系判定は、基本因子の Index の包含、対応する局所構造の一致、相対順序の保存、および積の種類の整合性に基づく。木の部分木や連続した因子列であることは要求しない。任意の因子集合を指定した部分系の取得と Regrouping の操作は提供しない。Tensor の直下の子空間と、全葉の Index は共有借用を基本として順序を保って列挙し、既存の構成を読み取る。空間は `Arc` で共有し、互換性・部分系判定には有効な積の種類と反復的な葉列の照合を用いる。深い空間の走査・破棄も反復的に行う。総次元の `usize` overflow は次元取得時に検出し、空間構築の拒否理由にしない。詳細は [Hilbert 空間の再帰構造](design/hilbert-space.md)に定める。

## 4. Index の契約

Index は局所構造を指定せずに独立に発行し、自由度はあらかじめ発行した Index と局所構造から構築する。局所状態・演算子には発行済みの Index を明示的に付与する。自由度の構築や局所式への付与で Index を自動発行する経路は設けない。

core は Index 型をジェネリックに扱い、同じ Index であることを照合できることを契約とする。具体的な表現・発行・一意性の確保は Index の提供側が担い、標準の具体実装は `qmb-basic` に UUID と seed 指定を要求しない発行 API を置く。同一ランタイム内および複数ランタイム間で、独立に発行した Index の一意性を仮定する。Index と局所構造の対応を全体で登録・固定する機構は設けず、合成・埋め込み・評価に必要な整合性は、その操作の対象から検証する。

詳細は [Index の契約](design/index.md)に定める。同一性の照合は標準の `Eq` で行う。Step 4 では UUID v4 の具体型を `qmb_basic::Id128` とし、`new()` と `Default` による seed 指定不要の新規発行を提供する。発行が必要な箇所だけ `I: Default` を要求し、専用の Index trait は設けない。

保存・復元と複数ランタイムの契約は以下のとおりとする。詳細は [Index の契約の 7 節](design/index.md#7-保存復元と複数ランタイム)に定める。

| 操作 | 契約 |
| --- | --- |
| 保存した対象の復元 | 元の Index を保持し、同じ自由度の記述を復元する。共有所有や同一のメモリ上のオブジェクトへの復元は要求しない |
| 保存した対象から独立した系の作成 | 新しい Index を発行し、空間・状態・演算子の参照を一貫して付け直す |
| 複数ランタイムからの持ち寄り | Index で因子を照合し、局所構造と空間の互換性を検証する |

保存形式と具体的な公開 API は未決とする。これらの契約は、保存機能や並列評価を初期版の必須実装に追加するものではない。

## 5. Primitive・式木の数学的能力

状態 primitive は内積と演算子の作用を受ける Act、演算子 primitive は Product・随伴・trace を提供する。Product・Act・随伴の trait は係数全体 `C` を型引数、返却要素を関連型 `Term` とし、`LinearCombination<Self::Term, C>` を返す。Act の演算子型は関連型 `Operator` で指定する。局所規則と恒等演算子は、同じ局所構造の primitive の有限線型結合を返す。数値制約は各実装に置き、演算 trait は局所式や局所構造に依存させない。内積・trace は従来どおり複素 scalar を返す。局所構造はこれらの能力を関連付け、次元、所属、表現パラメータ、恒等演算子と、結合則・作用・内積・随伴などの整合性を定める。

Rust の契約では、LocalStructure と追加の定義ラベルを廃止し、S・O と個別の能力実装で局所規則を定める。恒等元は `O: Identity` の値、次元は `Dimension` が返す primitive 型全体の span の次元とする。能力は利用するメソッドで要求し、保持・係数計算を伴わない遅延構築・map に数値制約を付けない。全能力をまとめた係数・primitive trait は導入しない。状態への作用評価に Product を要求せず、局所積を必要とする trace などで要求する。

Core の共通の式木は、加算・scalar multiplication・zero と、積・作用・随伴などの構築を提供する。局所式は `LocalState<S, O, C>`・`LocalOperator<O, C>` のハンドルが `Arc` で節点 enum を共有し、和・スカラー倍・積・作用・随伴の variant は対応する共通パーツを保持する。zero は空の和で表す。同じ式ハンドルを再利用した箇所で不変な部分式を共有し、別々に構築した同形の式を自動検出して共有する機構は設けない。

必要な評価時に線型性・半双線型性・共役線型性を用いて primitive の規則へ委譲し、構築時には式全体を展開しない。状態式の `act` は `Act` 節点による遅延した作用を構築し、`evaluate` は基底を指定せず局所規則に基づいて状態式を評価する。`ψ.act(&b).act(&a)` は `A(Bψ)` を表す。Norm、正規化、期待値、正規化期待値、commutator、Hilbert–Schmidt 内積などは共通機構から導出する。bra 型は初期版では導入せず、braket は内積で計算する。

評価器は状態評価・演算子の作用・演算子項の走査を分担し、必要な範囲を処理する。評価と、最後の所有者になった式節点の破棄は明示的な作業スタックで反復的に行う。式同士の二項演算子は所有値・借用値の4通りに対応し、所有ハンドルを移動、借用ハンドルは共有参照を複製する。名前付きメソッドは共有借用を基本とする。

Primitive の線型独立性は要求せず、数学的に同じ対象を表す結果の式構造が一致することも要求しない。指定基底への座標取得は別の接続として扱う。式の操作は入力を変更しない。評価前の source の自動保持、式の構造比較・表示は後続へ延期し、初期版の必須機能にしない。詳細は [局所 primitive・式木の数学的能力](design/primitive-expression.md)に定める。ここで定めた公開範囲・能力分離・評価経路は採用済みとし、非公開の作業領域や補助関数は実装時に定める。局所式・評価と API の移行は Step 2、空間・配置・束縛・複合作用・内積・full trace は Step 3 で実装済みである。

## 6. 基底の数学的能力

基底の定義を表す struct が対象の要素型を指定した `Basis<P, C>` を実装し、primitive の基底展開、基底要素の順序番号の取得、および順序付きの遅延列挙を提供する。状態基底の展開は `LinearCombination<S, C>`、演算子基底の展開は `LinearCombination<O, C>` を返す。各要素は選択した基底の primitive に限り、式木側で `map` して通常の式へ取り込む。結果に基底情報を付随させない。局所基底要素は primitive の部分集合とし、対象空間全体に線型独立かつ完全であることを要求する。

基底要素数は、対応する局所構造と対象種別から導出する。局所 Hilbert 空間の次元が `d` なら状態基底は `d`、全演算子基底は `d²` 要素となる。Basis に独立した要素数の提供を要求しない。

Core の共通機構が展開を線型に式全体へ拡張し、順序番号から座標を得る。内積・Gram 行列・双対基底は展開の実装に利用できるが、Basis の必須の独立能力とはしない。Step 2 では `Basis<P, C>` を維持し、`StateBasis`・`OperatorBasis` を削除した。利用側で状態には `Basis<S, C>`、演算子には `Basis<O, C>` を直接要求する。次元からの要素数取得と基底の数学的契約は維持する。Step 5 で `expand_in`・`coordinates_in`・Coordinates・BasisError と HilbertSpace 上の複合基底の遅延列挙を実装した。詳細な契約と公開 API は [基底の数学的能力](design/basis.md)に定める。

初期版では局所規則に基づく `evaluate` に加え、指定基底への展開と係数取得を提供する。分配だけ、局所積だけなど、特定の簡約ルールのみを選択して適用する独立した操作は future work とする。必要な分配・局所積・係数の結合は、基底展開や作用・内積・trace の内部処理として維持する。入力の式は変更せず、通常の評価に基底展開を必須としない。

複合基底は基本因子ごとの局所基底の tensor product に限定し、Hilbert 空間の因子順序と既存の mixed-radix に従う。Grouping は基底を変えない。展開結果は通常の式木とし、基底識別情報や展開済み状態を保持させない。別基底への変更は再展開で扱い、Index・局所構造・確定した空間を保持する。数値座標を export する場合は、その基底・順序などの解釈を別途確定する。

写像の表現行列・表現係数は future work とする。将来の数学的契約として、再利用可能な係数データは写像の入力・出力基底に対する表現行列として整理する。各入力の線型性／共役線型性を明示し、複数入力の係数列を指定順序で tensor product にまとめて行列を適用する。基底の列挙・通常の演算・出力基底への展開と係数取得から構成する共通操作とし、専用 trait の実装は要求しない。数式・共役指定・非正規直交基底・検証対象の詳細は [写像の表現行列](design/representation-matrix.md)に保持する。

## 7. 埋め込み・複合空間の数学的能力

局所演算子式と Index の組を複合演算子式の葉とし、和・積などで合成する。状態側も局所状態式と Index の組から tensor product と線型結合で構成する。Index の付与を配置、spectator を含む全作用空間の確定を束縛と呼び分ける。埋め込みは局所式の線型結合・積・随伴・恒等演算子を保存する。状態に新しい因子を追加する場合は、その状態または明示的な埋め込み写像を必要とする。

複合 Hilbert 空間が、局所構造に相当する整合性の責任を持つ。子の構造・因子順序・ordinary / graded の種類から積・作用・内積・随伴・trace などを合成し、共通の合成・埋め込み・評価機構がその契約を実現する。独立した複合構造型は要求しない。

同じ合成則を入れ子にも適用し、同種の積の括り方と空・単因子の同一視に対して数学的結果を一致させる。Graded では積・作用・tensor の随伴に符号規則を適用し、通常の正定値な内積と非正規化 trace を維持する。Partial trace は future work とし、埋め込み後の通常の作用を縮約して残存空間を保持する数学的契約のみを将来の実装向けに残す。

局所契約を満たすことは局所構造の提供者、契約を保存する合成則は共通機構の責任とする。操作時には適用の前提を検証し、数学的法則を任意の実装について実行時に自動証明することは要求しない。能力一覧・数式・検証対象の詳細は [埋め込み・複合空間の式と数学的能力](design/composite-expression.md)に定める。Ordinary / graded の共通評価を実装済み。空間が積の種類、primitive が parity / P、core が符号を処理する。

Step 3 の複合状態 `State` は常に所属空間を持ち、複合演算子 `Operator` は一つの型で未束縛・束縛済みを扱う。束縛情報は外側のハンドルに集約し、明示空間・既存束縛・未束縛の作用時の状態空間の順に選択する。検証が必要な複合操作とその演算子構文は `Result` を返し、空間不整合を早期に検出して局所規則の実行を遅延させる。評価は作用・内積・full trace の経路を分担し、全大域演算子の完全展開を必須にしない。採用方針と公開 API は[複合式の設計](design/composite-expression.md)に記録する。

## 8. 型選択・幾何・モデル構築

`qmb_core::QuantumSystem` は primitive の `State`・`Operator` と係数型 `Coeff` を関連付ける。`state(primitive)`・`operator(primitive)` は既存の局所式を返し、構築時の数値能力や Clone を要求しない。評価の能力は引き続き各操作で要求する。`LocalStateOf<Q>`・`LocalOperatorOf<Q>`・`StateOf<Q,I>`・`OperatorOf<Q,I>` は既存の式型への投影別名であり、系を表す識別子や新たな式ラッパーは保持しない。

`qmb_spin_half::SpinHalf<C>` で代数と係数型を一度選択する。既定精度・個別 primitive 用メソッド・異精度の暗黙変換は設けない。既存の core constructor は維持する。明示的な戻り値型には投影別名を使え、Q は構築入口で選ぶ。構築・集約は QuantumSystem の state_tensor_at・operator_product_at・state_sum・operator_sum・state_tensor・operator_product と集約の try_ 版に揃え、構築済み式の操作はインスタンスメソッドで行う。scaled は `impl Into<C>` を一度変換して保持する。try_ は最初の ExpressionError で入力消費を止め、全成功時だけ合成を検証する。契約は[複合式の設計](design/composite-expression.md#11-quantumsystem-の構築入口と-fallible-集約2026-10-28)に記録する。

Chain はサイト数と境界条件のみを持ち、Index の入力・発行・ToSpace による空間生成を行わない。`Chain::new(n, boundary)` は Self を返す。近接 Site 対は辞書順で、0・1 Site は辺なし、2 Site は open／periodic とも一辺とする。

外部の Site 順列を借用する `_with` API は、列長を列挙前に検証する。`nn_site_pairs_with` は入力要素対、`nn_pairs_with` は Site 間の全葉の直積、`nn_pairs_by_pos_with` は短い方までの葉の位置対応を返す。`_on` は保存された全空間のトップレベル分割を Site と読み `_with` へ委譲する。葉の列挙には core の `leaves()` を使い、要素・Index の Clone・Eq を要求しない。

モデル構築は SpinHalf<C> の関連関数 xyz・heisenberg・neel で提供し、Chain と局所空間列を受け取る。Index の重複と列長を辺の有無によらず検証する。Hamiltonian は局所空間ハンドルを共有して列順に tensor した空間へ束縛する。標準モデルは複数葉と空 Site を受け入れ、XYZ・Heisenberg は Site 間の全葉ペアへ結合を与え、Neel は Site ごとに向きを決めて全葉へ割り当てる。空 Site は位置と偶奇を保ち、前後を直接結び直さない。0 Site では単位空間上の zero Hamiltonian と単位状態を返す。

Index は引き続き Eq だけで照合し、共有に Clone・再発行を要求しない。局所 grouping は維持するが、別途保持する全空間の上位 grouping・異なる因子順序を復元する契約は持たない。契約と検証は [Chain の設計](design/geometry.md)に記録する。

## 9. 初期版の提供範囲と完成条件

**段階1〜6完了（2026-10-27）**。初期版の機能実装と受入確認を完了した。数値 export の検証結果は[数値 export の設計](design/export.md#7-受入確認)、MSRV・専用 no_std target・性能ベンチマークの未実施事項は[今後の課題](future-work.md#3-未実施の検証と性能)に記録する。

各設計文書の Step または段階は、1 が局所契約と具体例、2 が局所式と評価、3 が空間と複合評価、4 が Chain と導出操作、5 が基底展開と座標取得、6 が数値 export と受入確認を指す。

初期版は、spin-1/2 の 1D Chain を代数的に構築・評価して密・疎の数値表現へ変換する版とする。ユーザ定義の有限次元局所構造の単独利用、同種の基本因子による ordinary な入れ子空間、状態・演算子の演算と full trace、指定基底への展開・再展開・係数取得、および単精度・倍精度の複素数を提供する。

2026-10-27に、partial trace と写像の表現行列・表現係数を future work とし、Step 5 の担当範囲と初期版の必須機能・検証項目から外した。基底処理と数値接続の実装は、これらの機能を前提条件にしない。

初期版の標準の状態 primitive・基底は `up,down`、演算子 primitive は `I,Sx,Sy,Sz,S+,S-`、演算子基底は `{I,Sx,Sy,Sz}` と `{I,S+,S-,Sz}` とする。初期版から延期した fermion は第2.4節のとおり実装済みである。異種局所構造の混在・合成、source の自動保持・式の構造比較・表示、保存・復元、個別の簡約ルール、キャッシュ、並列化、追加の ML・tensor network 接続は引き続き後続とする。

完成条件は、Heisenberg chain と Neel 状態の代表シナリオ、独立な小規模行列との一致、基底展開・再展開と座標取得の一貫性、別の局所代数を独立した系で使う抽象化と非正規直交基底の検証、lazy 評価・共有、および実行可能な利用例と不正入力の検証とする。数値 export の API・backend・依存指定は[数値 export の設計](design/export.md)として実装済みである。後続へ延期した機能は[今後の課題](future-work.md)にまとめる。

## 10. 将来の必須課題

State・Operator・Structure・Space の管理方法は、将来必ず再検討し、より適切な設計へリファクタリングする。所有・共有・関連付け、および Geometry と空間の依存関係を含めて見直す。現行設計は当面の実装方針であり、最終形として固定しない。具体的な設計と実施時期は後続で決定する。

このリファクタリングと併せて、以下を必須課題として扱う。

1. **core の係数型・数値能力の整理（実装移行完了）**：`qmb-core` に `num-complex`・`num-traits` を含め、係数型には必要な操作で必要な trait bound を課す。零・単位元には `num_traits::Zero`・`One` をそのまま使い、付随する能力も受け入れる。不足する能力 trait は core で定義し、`Complex<f32>`・`Complex<f64>` に直接実装できる構成にする。保持・構築へ不要な数値能力を要求せず、係数表現・演算・norm・正規化の責任を整理する。
2. **LocalStructure の必要性（実装移行完了）**：primitive 型 S・O を式へ直接指定し、次元と恒等元を primitive 側へ移す。追加の定義ラベルは設けず、意味・規則は S・O と各能力の実装で定める。
3. **Operator と State の整理（初期版実装後）**：操作的量子論の枠組みで、Operator と State をいずれも一般の演算子の特別な形として整理できるかを議論し、解決する。その整理が数学的意味と API にどう対応するかを検討し、結論と根拠を記録する。2026-10-26に、検討事項が多いため初期版実装後へ延期した。初期版では現行の分離を維持し、この整理の結論を実装の前提条件にしない。
4. **HilbertSpace の最小構成（初期版の範囲は完了）**：LocalStructure 廃止により、Index と順序付き Tensor の構造だけを保持する構成へ移行した。primitive 型・次元・局所規則は式・取得時の入力・個別能力が担う。Grading の具体的な保持方法・API は 2026-10-29 に実装済み。第10.2節を参照。
5. **Chain の内部表現（実装移行済み）**：サイト数と境界条件だけを保持し、Index・Geometry 構造体・HilbertSpace の保持を廃止する。空間の構築・保持・grouping・全体の因子順序・配置と束縛は利用側へ移す。理由と契約は[Chain の設計](design/geometry.md)に記録する。Grading の設計完了を前提にしない。

初項は、2026-10-25 の方針変更により、数値ライブラリ依存の解消から、`num-complex` を含めた最小限の能力契約の設計・実装へ変更し、2026-10-26 に完了を確認した。未検証事項は[今後の課題](future-work.md#3-未実施の検証と性能)に記録する。LocalStructure は同日に廃止し、ordinary な空間の責任分離も実装済みである。Chain の Index 列保持・ToSpace への2026-10-27の移行を経て、2026-10-28に純粋な Chain と外部局所空間列・SpinHalf<C> の関連関数へ移行した。Operator と State の統合、Grading の具体的な保持方法・API のうち後者は実装済みであり、前者の Operator / State 統合は引き続き後続課題とする。

係数型については、成分型 T から係数全体 C への移行と、共役・実数値取得・実係数の埋め込み・純虚係数の生成を個別の能力として扱う方針を採用した。norm・normalized・正規化期待値と NormalizationError は std を利用する qmb-normalize へ移し、Float は平方根を使う実数型 R に要求する。内積・期待値・norm squared・スカラー除算は core に残し、provider は導入しない。標準 Complex の左スカラー倍は個別実装、任意の C は scaled(c) を基本とする。正規化期待値の純実分母を不要に二乗しない処理も維持する。これらの実装移行は完了し、追加能力は Conjugate・RealPart・ImaginaryUnit・Divide、実係数の埋め込みは From<R> で提供する。拡張 trait は qmb_normalize::{Norm, Normalize, NormalizedExpectation, NormalizedExpectationIn} とする。core に正規化専用 module・エラーは置かない。詳細は[正規化 API の契約](design/primitive-expression.md#41-norm正規化の-api-と-feature)に記録する。

LocalStructure 廃止後は、演算子式に O・C、状態式に S・O・C を直接指定する。次元は primitive 型全体の span の次元（S は d、O は d²）、恒等元は係数型に依存せず生成する O の値とする。空間・Geometry／Chain から L を外し、因子構造と局所規則を分離する。採用した契約と移行前の仕様の扱いは[局所構造の設計第7節](design/local-structure.md#7-localstructure-廃止の採用方針2026-10-26)に記録する。全 crate の製品コード・テスト・利用例を移行済みである。

この見直しは Step 4 の完了条件には追加しなかった。

### 10.1. Operator と State の整理

2026-10-26に初期版実装後へ延期した設計検討である。再開時期と実装範囲は別途決める。

操作的量子論の枠組みで両者を一般の演算子の特別な形として扱えるかを検討し、数学的意味と API の対応を決める。現在の State は ket であり、`ℂ→H` として始域・終域を参照できる。束縛済み Operator は `H→H` として参照できる。ただし、局所・複合それぞれに状態用と演算子用の式木・評価処理があり、一般の `U→V` の写像を表す型はない。

再開時には次の論点を検討する。

- 「一般の演算子」が表す数学的対象と、扱う状態の意味。Ket、密度演算子、量子操作を同じ対象として扱う前提を置かない。
- 合成、tensor product、随伴、内積、trace、zero・identity と始域・終域の対応。
- 数学的な共通化を公開型、内部表現、共通処理のどこまで反映するか。公開型を統合するかは別に判断する。
- 未束縛の配置式と、始域・終域が確定した写像の関係。

空間・幾何との境界を再確認し、型と評価器の統合を伴う実装は責任分担が定まってから範囲を決める。Grading の情報保持・API の具体化とは独立に検討する。完了条件は、採用する数学的対象、対応する API、共通化する範囲、別に維持する範囲を理由とともに記録することとし、統合を採用しない結論も含む。

### 10.2. Grading の情報保持と API

2026-10-29 に共通機構を実装した。`HilbertSpace<I>` が `TensorKind` と順序を保持し、primitive が `HomogeneousParity` / `ParityOperator` を提供する。共通の作用評価では Act ごとの選択済み空間に従って P を挿入する。演算子の intrinsic な基底展開と通常表現による trace は private な走査型を分ける。export / normalize は core の共通 API に接続する。

空・単因子の同一視後にも残る異種 tensor kind の混在は拒否する。幾何・Index に grading の責任は追加しない。構築 API・最小能力と breaking change は [複合式の共通評価と移行](design/composite-expression.md#12-ordinary--graded-の共通評価と移行)を参照。Fermion と Spin-one の製品実装は第2.4〜2.5節に記録する。Operator / State 統合は後続に残す。
