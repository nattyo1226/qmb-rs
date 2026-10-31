# 多体相互作用の関係列挙

作成日: 2026-10-30
更新日: 2026-10-31

基礎検討・初期試作の記録。本書の各案をそのまま正式 API とするものではない。

2026-10-31、num_sites の名称指定を含め、[幾何・Relation・space 列挙の統合設計](geometry-refactoring.md)を採用した。本書の単独 Relation と Pauli string サンプルは初期試作であり、採用した設計では Relation<G> による規則の適用と SiteMapper による接続へ発展させる。ライブラリ本体への実装は未着手である。

## 設計の起点となる利用例

利用者の問題意識は、隣接する HilbertSpace の tuple を返すだけの機能では、「高々 K Site に作用する Pauli string をすべて作る」処理を記述する際に大量の trait bound が必要になる、というものである。関係列挙の返却型だけでなく、局所演算子の選択から式構築までを通した使いやすさを設計対象にする。

先に提示した連続 window と固定長配列は、その一部を扱う候補に留める。連続するちょうど K Site の配置だけでは、高々 K Site という要求を満たせない。体数ごとに処理を手書きすることも前提にしない。

将来の2次元正方格子での方向付きの列や中心付き近傍も、拡張性の検討例とする。関係の意味は幾何ごとに定義し、すべてを連続 window に合わせない。共通の Site 列に加え、中心・方向などの役割を保持できる必要がある。具体的な整理は[中心・方向を持つ関係](geometry-refactoring.md#中心方向を持つ関係の扱い)に記録する。

まず1 Site が1 spin の場合を考え、「高々 K」を非恒等演算子が作用する Site 数（support の大きさ）と解釈して検討する。支持集合を任意に選べるか、幾何の関係で制限するか、恒等演算子を含めるかは、利用要件として確認する。複数葉の Site については[統合設計の選択方針](geometry-refactoring.md#複数葉を持つ-site-への接続)を参照する。Site 数と実際に作用する葉の数を区別し、Site 内の複数葉へ同時に作用する Pauli string も含める。

任意の支持集合を対象に恒等演算子も含めるなら、構築の流れは次になる。

1. `r = 0..=min(K, N)` の各大きさについて、Site の支持集合を重複なく列挙する。
2. 選んだ r Site に対する `{X, Y, Z}` の全割り当てを列挙する。
3. Site と局所演算子の対応を保って、既存の式構築 API へ渡す。支持集合の外には明示的な I を配置しない。
4. 得られた式を、相互作用に現れない Site も含む全作用空間へ bind する。

この定義では項数は `Σ(r=0..min(K,N)) binomial(N,r) × 3^r` となる。`r=0` は恒等演算子一つであり、window の長さ0を認めるかとは別の契約である。全件を Vec に集めず順に構築・消費できる経路を検討する。

geometry は支持集合・配置を担当し、`X/Y/Z` の選択や Pauli の規約は具体的代数、式の配置・積は core が担当する。利用者にこれらの接続を毎回手書きさせない補助の配置と必要性を、具体的なコード例で検討する。

なお現行の `SpinHalfOperator::Sx/Sy/Sz` は Pauli 行列の1/2である。Pauli string と spin 演算子の積を区別し、係数の規格化を暗黙に混同しない。

## trait bound を検討する際の基準

現行の `ExpressionFactory::operator_product_at` は、Index の Eq と寿命の条件のみで任意長の入力を扱い、係数型や primitive に評価・数値演算の bound を要求しない。tuple を返すこと自体が bound 増加の原因とはまだ断定しない。

- 関係の列挙、局所演算子の選択・規格化、配置・積構築、評価を分け、どこで bound が増えるかを利用例で特定する。
- 構築だけの利用に、`Product`・`Act` や係数の数値演算など、評価のための能力を要求しない。
- 実行時の上限 K に対しても、同じ構築処理を使えるか検討する。固定長配列で返す入口だけを共通 API と決めない。
- 必要な能力を一つの包括 trait に隠すだけの解決にはせず、要求する能力そのものを最小化する。局所ラベルの再利用や規格化に必要な能力は個別に説明する。

## 現行の制約

- `Chain::nn_site_pairs()` は `(usize, usize)` を返し、各ペアを `i < j` に正規化して辞書順に列挙する。periodic の2 Site も一辺だけとする。
- 接続処理も2要素に固定されている。Site の要素を借用する処理、Site 間の葉の全組み合わせ、葉の位置対応が別の処理として実装されている。
- core の `ExpressionFactory::operator_product_at` は既に任意個数の配置を受け取り、入力順に積を構築できる。多体相互作用のために core の演算子積を新設する必要はない。

したがって、まず幾何による多体の組の列挙と、その組から既存の式構築へ渡す接続を検討する。

参照: [Chain](../../../crates/qmb-geometry/src/geometry/chain.rs)、[最近接ペアの列挙](../../../crates/qmb-geometry/src/relation/nearest_neighbor.rs)、[借用接続](../../../crates/qmb-geometry/src/connection.rs)、[式構築](../../../crates/qmb-core/src/system.rs)。

## 順序付きの配置と、順序を持たない集合

最初に決めるべき点は、列挙する組の各位置に役割があるかどうかである。例えば3体の項 `X_i Z_j X_k` では、中央の位置にどの Site を置くかが意味を持つ。

3 Site の periodic な連続配置は、次のように列挙する案がある。

```text
[0, 1, 2] → X_0 Z_1 X_2
[1, 2, 0] → X_1 Z_2 X_0
[2, 0, 1] → X_2 Z_0 X_1
```

Site 集合は同じでも異なる項になるため、組を昇順にそろえて一つにまとめることはできない。graded な積では、演算子の順序を変えると符号にも影響する。

この連続配置の API については、順序付きの配置を返し、組内の並べ替えや Site 集合による重複除去をしない案である。`Z_i Z_j Z_k` のように配置を回転しても同じ項になるモデルでは、並進した項をすべて加えるか、同じ支持集合を一度だけ数えるかをモデル側で定める。相互作用の対称性を geometry が推測しない。

一方、全 Pauli string の列挙では、支持集合を一度ずつ選び、その上の局所ラベル割り当てをすべて生成する方が重複を避けやすい。支持集合を正規順で表しても、各 Site と演算子の対応を保てば役割は失われない。順序付きの配置と支持集合は異なる列挙契約として扱う。

## 組の表現

| 候補 | 利点 | 注意点 |
| --- | --- | --- |
| `[usize; K]` / `[&T; K]` | 体数を型で表し、`[a, b, c]` と分解できる。組ごとのヒープ確保が不要 | 体数を実行時だけに決める用途は別途検討が必要 |
| `Vec<usize>` / `Vec<&T>` | 実行時の体数や、異なる長さの組を同じ列で扱える | 組ごとの確保と長さの実行時確認が必要 |

ちょうど K Site の配置には const generics による固定長配列が候補になる。一方、高々 K Site の列挙では異なる長さを同じ処理へ流せる表現も必要になる。可変長列、借用 slice を受けるコールバックなどの利用形と比較し、共通 API を決める。借用データを次の列挙時に再利用する方式では、その寿命と所有値として保持する経路も明示する。

## 有力案: 可変長の組と Relation の共通インターフェース

2026-10-31 の検討では、組の格納形式と、利用者が組を処理するインターフェースを分ける案を有力候補とする。

### 固定長配列を「最大長」として使う

長さごとに異なる配列型を返す代わりに、非公開の固定長配列と実際の長さを保持する。以下は API 案の抜粋で、constructor の検証などは省略している。

```rust,ignore
pub struct SiteGroup<const CAP: usize> {
    sites: [usize; CAP],
    len: usize,
}

impl<const CAP: usize> AsRef<[usize]> for SiteGroup<CAP> {
    fn as_ref(&self) -> &[usize] {
        &self.sites[..self.len]
    }
}
```

空の組、`[0]`、`[0, 2]`、`[0, 2, 4]` を、すべて `SiteGroup<3>` として同じ Iterator で返せる。CAP は格納できる最大長であり、その組の体数ではない。未使用の要素は公開せず、構築時に `len <= CAP` を保証する。Site 番号を保持するため、Index の Default・Clone などは必要ない。

この型なら組ごとのヒープ確保は不要である。ただし容量 CAP に応じた領域と移動コストがあり、大きな上限に常に適するわけではない。上限が実行時にしか決まらない場合には Vec を使う経路も必要になる。

### 格納形式を Relation の利用側へ露出させない

```rust,ignore
pub trait Relation {
    type Group<'a>: AsRef<[usize]>
    where
        Self: 'a;

    fn num_sites(&self) -> usize;

    fn groups(&self)
        -> impl Iterator<Item = Self::Group<'_>> + '_;
}
```

組の型が slice として読めることを trait 自体で保証する。共通の利用側は `R: Relation` とだけ記述でき、各体数の tuple・配列用 trait や、組に対する追加の AsRef bound を書く必要がない。これは評価能力を包括 trait にまとめる案ではなく、列挙結果を読むための最小の契約である。

| Relation の用途 | Group の候補 |
| --- | --- |
| 最近接ペア | `[usize; 2]` |
| 高々 CAP Site の支持集合 | `SiteGroup<CAP>` |
| 実行時の上限に従う支持集合 | `Vec<usize>` |
| 保存済みの任意の関係 | `&'a [usize]` |

GAT の寿命は、保存済みの関係をコピーせず借用する用途に使う。生成中の Iterator 自身が持つ可変バッファへの参照を返し、その次の next で書き換えることは、この通常の Iterator の契約ではできない。生成型 Relation は組を所有値として返すか、別途コールバック方式を検討する。

Relation は Site 番号を列挙し、接続処理が呼び出し元の Site 列を借用する。列長と番号範囲を検証してからアクセスする。異なる体数の組も同じ接続処理と `operator_product_at` へ流す。結果の式は SiteGroup や一時バッファを借用し続けない。

trait の GAT や Iterator の具体型は実装側が受け持つ。動的 dispatch が必要になれば別の検討が必要だが、今回の利用例では前提にしない。保存済み関係の借用が不要なら、GAT のない単純な関連型でも同じ方針を実現できる。

### 列挙規則は具体的な Relation で選ぶ

共通 trait は体数・連続性・順序の対称性を決めない。それらは具体的な列挙規則として明示する。以下は役割を説明する仮称である。

- `SupportCombinations` は、指定された Site の支持集合を、指定サイズの範囲で一度ずつ列挙する。高々 K Site の全 Pauli string の用途に対応する。
- `NearestNeighbors` は、従来の正規化された最近接ペアを列挙する。
- `SiteWindows` は、連続する Site の順序付き配置を列挙する。
- `ExplicitRelations` は、利用者が指定した Site の組を、その順序で読む。

その上に、局所ラベルの全割り当てを行う代数側の補助と、既存の式構築 API を接続する。geometry に Pauli 固有の規則を持たせず、利用者が体数別の処理を手書きしなくて済む構成を検討する。

### 試作で確認した範囲

2026-10-31、リポジトリ外の一時パッケージで上記の Relation インターフェースを試作し、既存の qmb-core に接続した。Rust 1.99.1 の direnv 環境で offline の cargo test を実行し、4件が成功した。

- 固定長配列、可変長の SiteGroup、保存済み slice を同じ `R: Relation` の式構築関数で扱えること。
- Index の Eq 以外に、state・operator primitive・係数の Clone や数値演算・評価能力を要求せず、既存の `operator_product_at` で構築できること。局所ラベルは生成関数から新しい値を受け取った。
- 3 Site・高々2 Site で、恒等演算子を含む局所ラベル列が37通りになること。空系、容量0、Site 数より大きな指定、保存済みの異なる長さの組も確認した。
- Site 列長の不一致と、範囲外の Site 番号を検出できること。

この試作は no_std のライブラリとしてコンパイルしたが、テストはホストで実行した。Pauli の規格化・数値評価、全作用空間への bind、性能、MSRV 1.85、専用 no_std target は未検証である。正式な実装・公開名の確定には進んでいない。

### Pauli string を構築する実行例

後続の試作を[pauli_strings.rs](examples/pauli_strings.rs)に保存した。Relation・SiteGroup・SupportCombinations は、このサンプル内だけにある提案 API の実装である。既存のライブラリへ公開 API を追加していない。1 Site は1 spin の葉であり、空の支持集合を恒等演算子として含める。

利用側のプログラムは次の形になる。

```rust,ignore
let sites: Vec<_> = (0..3).map(HilbertSpace::factor).collect();
let relation = SupportCombinations::<2>::new(sites.len());
let mut strings = Vec::new();

for_each_pauli_string::<Complex64, _, _>(&relation, &sites, |word| {
    strings.push(word);
})?;

assert_eq!(strings.len(), 37);
```

構築関数は `C: From<f64>, I: Eq, R: Relation` を要求する。前処理で Site 列長、各要素が葉であること、Index の一意性を確認し、Site 列から全作用空間を構築する。各組の番号範囲・重複を確認した後の本体は次のようになる。

```rust,ignore
for group in relation.groups() {
    let support = group.as_ref();

    for axes in axis_assignments(support.len()) {
        let factors = support.iter().zip(axes).map(|(&i, axis)| {
            SpinHalf::<C>::operator(axis).scaled(2.0).at(&sites[i])
        });
        let word = SpinHalf::<C>::try_operator_product(factors)?.bind(&space)?;
        emit(word);
    }
}
```

`axis_assignments(r)` は長さ r の `Sx/Sy/Sz` の全割り当てを一つずつ返す。例えば支持集合 `[0, 2]` なら、`XX, XY, XZ, YX, YY, YZ, ZX, ZY, ZZ` の9個を作る。`Sx/Sy/Sz` をそれぞれ2倍するため、Site 1 が spectator の `X_0 I_1 Y_2` など、Pauli の規格化で構築される。空の支持集合は空の積となり、全作用空間上の恒等演算子を一つ生成する。

体数は slice の長さで扱うため、1体・2体などで型や構築コードを分けない。`C: From<f64>` は2倍の係数を作るための制約であり、構築自体に係数の Clone・数値演算や局所演算子の評価能力は必要ない。全件の保存も必須ではなく、emit へ順に渡せる。サンプルの局所ラベル列には、分かりやすさを優先して Vec を使用している。

同じ支持集合を複数回返す Relation は同じ項を複数回生成する。全 Pauli string を一度ずつ作る場合は SupportCombinations のような、支持集合を一度ずつ返す Relation を選ぶ。構築関数は関係全体の重複除去を行わない。

リポジトリルートから以下で実行できる。サンプルは[独立した Cargo パッケージ](examples/Cargo.toml)とし、本体 workspace に組み込まない。

```sh
direnv exec . cargo run --offline --manifest-path docs/plans/v3/examples/Cargo.toml --bin pauli-strings
direnv exec . cargo test --offline --manifest-path docs/plans/v3/examples/Cargo.toml
```

2026-10-31、Rust 1.99.1 で実行し、37個の出力を確認した。テスト4件、Clippy の全 target・警告をエラーとする確認が成功した。3 Site・高々2 Site の37個について、独立な Pauli 行列の tensor product と全行列成分を比較し、重複がないことを確認した。空系、恒等演算子だけの列挙、CAP が Site 数より大きい場合、不正な入力、非 Clone の Index・係数での構築も確認した。このサンプルの性能、MSRV 1.85、専用 no_std target は未検証である。

## 連続する Site の配置の規約案

連続配置も共通の Relation として提供する方向で検討する。公開名は仮であり、利用形は次を候補とする。これは実行可能な現行 API の例ではない。

```rust,ignore
let windows = chain.windows::<3>()?; // Relation を実装する型

for group in windows.groups() {
    let positions = group.as_ref(); // &[usize] として読める。
    // 共通の接続・式構築処理へ渡す。
}
```

この例の `K` は参加する Site 数である。複数葉を持つ Site では、直ちに局所自由度の体数と同じになるわけではない。

連続配置の境界規約として、以下を推奨する。

| 条件 | 列挙する配置 |
| --- | --- |
| open、`1 <= K <= N` | 開始位置 `0..=N-K` ごとに `[i, i+1, ..., i+K-1]` |
| periodic、`1 <= K <= N` | 開始位置 `0..N` ごとに、各位置を N で折り返した配置 |
| `K > N` | 配置なし。同一 Site を繰り返して K 要素を埋めない |
| `K == 0` | 不正な指定としてエラー。零体の定数項は式構築側で扱う |
| `K == 1` | 各 Site を一度ずつ列挙 |

`N == 0`・`K >= 1` は配置なしとする。`_with` の入力列長は、配置が空になる場合も Iterator を返す前に検証する。演算は `usize` の overflow や0による剰余を避けて実装する。

例えば `N=5, K=3` では次のようになる。

```text
open:     [0,1,2], [1,2,3], [2,3,4]
periodic: [0,1,2], [1,2,3], [2,3,4], [3,4,0], [4,0,1]
```

既存の `nn_site_pairs()` は現在の正規化・重複除去・列挙順を維持する。新案の `windows::<2>()` による配置は periodic の向きと項数が異なるため、既存 API の単純な置き換えには使わない。特に2 Site の periodic では、新案は `[0,1]`・`[1,0]` を返すが、従来 API は `(0,1)` の一辺だけである。

## Site から局所自由度への接続

幾何は Site 番号の slice として読める Group を返し、外部列への接続は可変長の組をそのまま読める借用ビューを検討する。要素・Index の Clone・Eq を要求しない。列挙は遅延処理を基本とし、関係の全件収集を前提にしない。`qmb-geometry` の `no_std` を維持する。

複数葉の Site を扱う場合は、次の選択を明示する。

- Site 自体の借用: 利用側が各 Site 内の葉の役割を選ぶ。Hubbard の spin 選択のような処理もここで行える。
- 葉の全組み合わせ: K Site の葉集合の直積を列挙する。最初の Site の葉を最も外側のループとし、各 Site 内の葉順を保つ。
- 葉の位置対応: K Site の同じ葉位置を組にし、最も短い Site の葉数まで列挙する。

空 Site の位置は保持し、それを含む葉の組は生成しない。空 Site を取り除いて前後をつなぎ直さない。`_in` を提供する場合は現行同様、トップレベルの分割を Site と読む。

全組み合わせ・位置対応の具体的な公開名と、共通の接続補助として公開するかは引き続き検討する。関係の種類ごとに `_with`・`_in`・全組み合わせ・位置対応を増やす場合の API 数にも注意する。

crate 内の責務分担と具体的な移行手順は[qmb-geometry のリファクタリング案](geometry-refactoring.md)に記録する。

## 次に詰める点

- 高々 K Site の Pauli string を構築するコード例を先に具体化する。支持集合の範囲、恒等演算子の扱い、Site と葉の対応を確認し、必要な trait bound を列挙する。
- 順序付き配置と、小さい periodic 系での項数を上記の規約でよいか確認する。順序を持たない多体の集合が必要なら、対称性と重複除去を別の契約として定める。
- 最初の対象を連続配置とするか、`[0, 2, 4]` のような相対位置パターンの並進も含めるか。後者では、境界で異なる相対位置が同一 Site に重なる場合の扱いも必要になる。
- 体数が実行時に決まる用途、任意の組の明示入力、接続補助の共通化が必要か。一般の Relation trait を導入するかは、具体的な利用形を確認してから決める。

## 実装時の検証案

高々 K Site の Pauli string は、小規模な系で独立に定めた全列挙と比較し、支持集合・局所ラベル・件数・重複の有無を確認する。構築だけの利用例で評価用の trait bound を必要としないことも確認する。

3体・4体について open／periodic の期待列を手作業で定め、0 Site、`K < N`・`K == N`・`K > N`、不正な `K == 0` を確認する。列長不一致は列挙件数によらず拒否し、借用先の同一性と非 Clone・非 Eq 要素も確認する。

非対称な演算子配置を用いて、周期境界で組内の役割が保たれることを確認する。複数葉・空 Site・入れ子の全組み合わせと位置対応、多体の項の手作業構築との一致、graded な順序の維持、従来の2体モデルの回帰確認を行う。
