# 幾何・Relation・space 列挙の統合設計

作成日: 2026-10-31
更新日: 2026-10-31

**採用済み・本体未実装（2026-10-31）**。利用者の承認により本設計を採用した。対応付けの型名は SiteMapper、Site 数の取得名は num_sites とする。エラー型などの実装詳細と追加の補助 API は実装時に確定する。[多体相互作用の基礎検討](multi-site-relations.md)と[実行可能な統合試作](examples/geometry_relations.rs)を参照する。

採用した構成は、最小の Geometry、geometry ごとの Relation<G>、可変長の Group を扱う SiteGroup、外部列への SiteMapper、結果を借用する MappedGroup である。葉の標準操作は leaves_cartesian と leaves_zip とし、共通入口に combine_leaves(FnOnce) を設ける。専用の選択 trait は追加しない。

## 起点となる二つの要求

1. 幾何を利用して、連続・隣接などの条件を満たす Site の HilbertSpace を列挙する。Chain の連続 n Site に限らず、将来の正方格子の方向付きの列や中心付き近傍にも対応できるようにする。
2. 同じ geometry に最近接・連続配置・支持集合など複数の relation を適用し、共通の接続・式構築処理へ渡せる trait 設計にする。

前の単独 Relation 案では geometry との対応が型の契約に表れていなかった。今回は **Relation<G> を geometry G に適用する規則**とし、geometry と外部の space 列の対応を借用ビューにまとめる。

## 採用する三つの役割

| 部分 | 責務 |
| --- | --- |
| Geometry | Site 数と、具体的な幾何が持つ構造・境界条件 |
| Relation<G> | G に適用する配置・近傍・支持集合の列挙規則 |
| SiteMapper<G, T> | 一つの geometry と外部の Site 順の &[T] を対応付け、関係を T の借用へ写す |

SiteMapper の T を HilbertSpace<I> にすれば space 列挙になる。元の空間・Index は利用側が所有し、geometry や relation が発行・複製しない。

### 実装する型と入口の整理

| 型・trait | 主な入口・読み取り |
| --- | --- |
| Geometry | num_sites()。既存 Chain に trait を実装する |
| Relation<G> | groups(&geometry)。G に対して有効な Site 番号の組を列挙する |
| SiteGroup<CAP> | 最大容量の配列と有効長。AsRef<[usize]> で有効な部分を読む |
| SiteMapper<'g, 's, G, T> | new(&geometry, &values)、groups(&relation)。外部列との接続を検証する |
| MappedGroup<'s, T, Group> | iter() で &T、source() で元の Group を読む |
| HilbertSpace 用の MappedGroup | combine_leaves(FnOnce)、leaves_cartesian()、leaves_zip() |

全体の利用形は以下になる。これは設計用サンプルで確認した構成であり、ライブラリ本体への import はまだできない。

```rust,ignore
let chain = Chain::new(4, Boundary::Open);
let mapper = SiteMapper::new(&chain, &site_spaces)?;

for group in mapper.groups(&Consecutive::<3>)? {
    let group = group?; // MappedGroup

    for factors in group.leaves_cartesian() {
        // Vec<&HilbertSpace<I>>。ここで局所演算子を配置して積を構築する。
    }
}

// 同じ mapper に別の relation と葉の組み合わせを指定できる。
for group in mapper.groups(&NearestNeighbors)? {
    let group = group?;
    for factors in group.leaves_zip() {
        // 同じ葉位置にある局所自由度の組。
    }
}
```

MappedGroup は Site 番号の組を外部列へ写した結果であり、Relation が返す Group と区別する。leaves の三つの入口は T=HilbertSpace<I> の場合にだけ提供する。独自の選択 trait は追加しない。

### Geometry は最小の共通能力にする

```rust,ignore
pub trait Geometry {
    fn num_sites(&self) -> usize;
}
```

Chain はフィールドと幾何の責務を保ち、この trait を実装する。Site 数の取得メソッドは、正式実装時に現行の site_count() から num_sites() へ統一する。方向への移動・座標などは具体的な幾何が提供する。全幾何に共通の「連続 n Site」や方向型を要求しない。

現時点の統合試作では、Geometry::num_sites() から現行 Chain::site_count() へ委譲するアダプターを使う。これは未変更のライブラリ本体に接続するためのもので、採用 API の名称は num_sites() である。

将来の正方格子では、Site 番号を座標へ対応させる処理と方向への移動を、その具体型が持つ。Site 順の外部列への番号付けと、幾何上の隣接・方向は区別する。

### Relation は geometry を内部に保持する規則型に固定しない

```rust,ignore
pub trait Relation<G: Geometry> {
    type Group<'a>: AsRef<[usize]>
    where
        Self: 'a,
        G: 'a;

    fn groups<'a>(
        &'a self,
        geometry: &'a G,
    ) -> Result<impl Iterator<Item = Self::Group<'a>> + 'a, RelationError>;
}
```

relation 自体は設定を保持し、適用時に geometry を受け取る。同じ geometry に複数の規則を適用でき、同じ規則型を複数の geometry に実装することもできる。

```rust,ignore
impl Relation<Chain> for NearestNeighbors { /* ... */ }
impl Relation<Chain> for Consecutive<3> { /* ... */ }

impl Relation<SquareLattice> for NearestNeighbors { /* ... */ }
impl Relation<SquareLattice> for AlongAxis<3> { /* ... */ }
impl Relation<SquareLattice> for ClosedNeighborhood { /* ... */ }
```

Consecutive と AlongAxis の上記 impl は説明用の一例であり、実装は const generics で長さを一般化する。Geometry に単一の関連 Relation 型を割り当てる設計にはしない。適用できない組み合わせは Relation<G> の未実装として型で区別できる。

幾何の隣接を使わない SupportCombinations<CAP> は、Geometry の Site 数だけを使う blanket 実装を候補とする。高々 K の実行時上限を CAP 内で指定する経路と、実行時容量の Vec を使う経路も、この利用インターフェースへ接続できる。

GAT は Group が geometry・規則内の保存済みデータを借用する場合にも対応する。生成 Iterator の可変バッファへの参照を next ごとに使い回す lending iterator にはしない。生成型では SiteGroup などの所有値を返す。

### SiteMapper から HilbertSpace を借用する

```rust,ignore
let mapped = SiteMapper::new(&chain, &site_spaces)?;

for group in mapped.groups(&NearestNeighbors)? {
    let group = group?;
    // group.iter(): 隣接する Site の &HilbertSpace<I>
}

for group in mapped.groups(&Consecutive::<3>)? {
    let group = group?;
    // 同じ geometry・同じ space 列から、連続3 Site を借用する。
}
```

SiteMapper は geometry と外部列への参照だけを持つ。新しい空間の所有者や Index 管理者にはしない。SiteMapper を使わずに呼び出せる接続用の補助関数を提供するかは、実装時に判断する。

各 MappedGroup は元の Group と外部列への参照を保持する。共通の iter() で組内順の &T を列挙し、source() で関係固有の情報へアクセスできる。要素・Index の Clone・Eq は接続のためには必要ない。

geometry と対応付ける外部列は、この SiteMapper に一度渡す。列長の一致は確認できるが、利用側が意図した座標と列順の意味まで自動検証できるわけではない。

## 中心・方向を持つ関係の扱い

Group をすべて裸の配列・slice へ変換しない。Relation ごとの具体的な Group 型に、中心・方向などの情報を持たせる。AsRef<[usize]> は共通走査のための最小能力であり、情報の全体ではない。

```rust,ignore
let mapped = SiteMapper::new(&square, &site_spaces)?;

for group in mapped.groups(&AlongAxis::<3> { axis: Axis::X })? {
    let group = group?;
    // x 方向に並ぶ3 Site の space。y 方向の行折り返しを混ぜない。
}

for neighborhood in mapped.groups(&ClosedNeighborhood)? {
    let neighborhood = neighborhood?;

    let center = neighborhood.center(); // &HilbertSpace<I>
    for (direction, neighbor) in neighborhood.neighbors() {
        // 中心との相互作用などを構築できる。
    }
}
```

center()・neighbors() は、近傍用の Group を保持する接続ビューへの具体的な補助として提供できる。試作ではこの方法を確認した。将来、複数の近傍 Group を汎用に扱う必要が出た場合に、その能力だけの追加 trait を検討する。すべての消費側へ中心・方向の bound を要求しない。

通常の十分大きい open 正方格子では、中心込みの近傍長は角・辺・内部で3・4・5と変わる。最大長5の配列と有効長で格納でき、各隣接先に方向を対応付ければ、境界で欠けた方向も取り違えない。

中心と各近傍の2体項を作るか、近傍全体で一つの多体項を作るかはモデル側が決める。近傍の列挙と、無向辺を一度ずつ数える最近接列挙は別の規約である。

## SiteGroup と列挙の契約

| 関係 | Group の例 | 主な規約 |
| --- | --- | --- |
| NearestNeighbors | [usize; 2] | 無向辺の正規化・重複除去を具体的な幾何ごとに定める |
| Consecutive<K> / AlongAxis<K> | [usize; K] | 各位置の役割を保つ順序付きの列 |
| SupportCombinations<CAP> | SiteGroup<CAP> | 指定サイズ範囲の支持集合を一度ずつ列挙 |
| ClosedNeighborhood | 有効長・中心・方向を持つ具体的 Group | 中心と各近傍の役割を保つ |
| 明示指定の関係 | 借用 slice など | 指定された組内順と列挙順を保つ |

SiteGroup<CAP> は最大容量の配列と有効長を持つ型として維持する。異なる有効長を同じ型で返せる。すべての Group をこの型へ統一する必要はない。

支持集合は昇順・組内重複なしにできるが、順序付きの列を勝手にソートしない。小さい周期系での折り返し、同一 Site の再出現、方向の違いによる重複は、それぞれの relation が規約を定める。

特に2 Site の periodic な最近接辺は従来通り一辺とし、順序付きの列 [0,1]・[1,0] の二配置と区別する。空の支持集合は一つの恒等項に対応するが、連続列の長さ0を認める必要はない。

## 検証とエラーの責任

- Geometry の構築: Site 数の overflow など、幾何自体の不整合を検出する。
- Relation の適用: 長さ0など、不正な設定を Iterator の返却前に検出する。
- SiteMapper の構築: geometry の Site 数と外部列長の一致を、関係が空でも検証する。
- 各 Group の接続: 独自 Relation の範囲外番号を、要素へのアクセス前に検出して Result で返す。全件の事前収集はしない。
- モデル・core の空間構築: Index の一意性、配置可能性、全作用空間への bind を検証する。

空間の組は借用して列挙し、暗黙に新しい tensor 空間へ合成しない。元の空間内の grouping、順序、空の Site を保つ。接続時のエラーまでに消費済みの結果を巻き戻す契約にはしない。具体的な公開エラー型は未確定である。

## Pauli string の構築へつなぐ

幾何・関係・接続を扱う汎用関数の制約は G: Geometry, R: Relation<G> でよい。Group に関する条件は Relation の契約内に含める。

その後は MappedGroup::iter() の space と局所ラベルを zip して、既存の try_operator_product へ渡せる。構築段階の追加制約は Index の Eq と、必要な係数変換だけに留める。中心・方向を使わない全割り当てに、それらの能力を要求しない。

X/Y/Z の全割り当て、spin と Pauli の規格化、全作用空間の選択はモデル・spin-half 側、配置・積・bind・評価は core が担当する。geometry はどの演算子を配置するかを知らない。

## 複数葉を持つ Site への接続

SiteMapper<G, HilbertSpace<I>> は各 Site の空間をそのまま借用する。一つの Site が複数葉の tensor 空間でもよい。例えば Site 0 が A⊗B、Site 1 が C⊗D なら、最近接の Group は &[A⊗B, C⊗D] に相当する借用ビューになる。

Site 内の局所自由度は既存の HilbertSpace::leaves() でテンソル積順に列挙できる。直下の grouping を読みたい場合は children() がある。これらの列挙 API を置き換えず、複数の leaves Iterator の組み合わせだけを補助する。

### 標準の操作は直積と zip の二つ

採用した標準の葉の組み合わせは次の二つである。現在は設計用サンプルにのみ実装している。

```rust,ignore
for group in mapped.groups(&relation)? {
    let group = group?;

    for leaves in group.leaves_cartesian() {
        // 各 Site から葉を一つずつ選ぶ全組み合わせ。
    }

    for leaves in group.leaves_zip() {
        // 各 Site の同じ葉位置を、最も短い葉列まで組にする。
    }
}
```

例えば [A⊗B, C⊗D] に対し、直積は [A,C], [A,D], [B,C], [B,D]、zip は [A,C], [B,D] を返す。Site が3個以上、または Group ごとに Site 数が異なる場合も同じ API を使う。

結果は Vec<&HilbertSpace<I>> など、長さを実行時に扱える列とする。先頭の Site を最も外側のループとし、各 Site 内の leaves() の順序を維持する。結果の位置 k は元の Group の Site 位置 k に対応するため、外側の Group が持つ中心・方向などの情報と対応付けられる。

### 共通の入口は Iterator を組み合わせる FnOnce

二つの補助を実装する共通の入口を、独自の組み合わせにも使えるようにする。

```rust,ignore
impl<'s, I, Group: AsRef<[usize]>> MappedGroup<'s, HilbertSpace<I>, Group> {
    pub fn combine_leaves<Out: Iterator>(
        &self,
        combine: impl FnOnce(Vec<Leaves<'s, I>>) -> Out,
    ) -> Out {
        combine(self.iter().map(HilbertSpace::leaves).collect())
    }

    pub fn leaves_cartesian(&self)
        -> impl Iterator<Item = Vec<&'s HilbertSpace<I>>>
    {
        self.combine_leaves(cartesian_leaves)
    }

    pub fn leaves_zip(&self)
        -> impl Iterator<Item = Vec<&'s HilbertSpace<I>>>
    {
        self.combine_leaves(zip_leaves)
    }
}
```

Leaves は現在の qmb_core::space::Leaves を使う。FnOnce へ渡すのは Site ごとの Iterator の Vec であり、Site 数に対応する tuple trait を必要としない。各 Group で一度呼ぶだけなので FnOnce で足り、Fn・FnMut として使えるクロージャも渡せる。

```rust,ignore
let leaves = group.combine_leaves(|streams| {
    // Vec<Leaves<'s, I>> を受け取り、一つの Iterator を返す。
    my_combiner(streams)
});
```

独自の出力には、固定長配列・Vec・Site ごとの区切りを保持する型などを選べる。共通入口では Out::Item の形を固定しない。独自 Fn が順序・対応を変える場合、その意味はその Fn と利用側の契約にする。

新しい FactorSelection trait、NonEmptySubsetsPerSite、ExplicitFactorChoices などは標準 API に追加しない。必要な処理は既存の Iterator 操作か、この共通入口へ渡す関数で記述する。

### 直積の再走査と空の扱い

現行 Leaves は Clone を実装していない。設計用サンプルの直積では、読み取った葉への参照をキャッシュして再走査する。入力 Iterator や HilbertSpace・Index の Clone は要求しない。

直積の全結果を事前生成することはなく、葉への参照のキャッシュは、全件消費時で入力の葉数の和に比例する。zip は入力 Iterator を順に進め、最短の列で終了する。正式実装時にキャッシュと leaves() の再生成のコストを比較する余地はあるが、Iterator の Clone を利用側へ要求する前提にはしない。

空 Site を一つでも含む場合、直積・zip とも結果は空とする。空 Site を取り除いて前後の位置を詰めない。

Site の組自体が空の場合は、両方とも空の葉列を一つ返す案とする。恒等演算子の構築へつなげられ、特に zip で空列を無限に生成しないようにする。この0入力の規約は明示しておく。

葉の組み合わせは既存の有効な参照を扱うため、標準の二つの操作に新たなエラー型を追加しない。Site 番号の範囲は Group の接続時に検証済みである。独自 Fn が必要とする追加の検証は、その Fn が扱う。

### 全 Pauli string との境界

高々 K Site の全 Pauli string を作る際、Site に複数の spin-1/2 の葉があるなら、標準の直積・zip だけで全項を作れるわけではない。1 Site に m 葉あれば、その Site の非恒等な Pauli string は 4^m−1 種類である。

この用途では、spin-half 側または利用側の関数で、各 Site の葉の非空な部分集合を選び、その選んだ葉へ X/Y/Z を全割り当てする。既存 Iterator 操作や combine_leaves へ渡す独自 Fn で実装する方向とし、geometry の標準操作や選択 trait を増やさない。

Site 数と配置する葉の数を区別する。K/CAP は作用する Site 数であり、配置する葉はそれを超えてよい。既存 Pauli サンプルの複数葉への全割り当ては未実装である。

全空間のトップレベルを Site と読む現行 _in の経路も、共通の接続への補助にできる。この部分では qmb-core の HilbertSpace を参照するため、geometry から core への現行依存を維持する。正式実装では no_std と必要な alloc の範囲を保つ。

## リファクタリングの範囲と手順

1. Geometry と Relation<G> の共通契約、SiteGroup、SiteMapper と接続ビューを整える。Chain のフィールド・構築入口は維持し、Site 数の取得を num_sites() に統一する。
2. 現行の最近接列挙を NearestNeighbors の Relation<Chain> 実装へまとめる。初期の nn_* API は互換ラッパーとして委譲し、小さい周期系・列挙順・検証タイミングを保つ。
3. Consecutive と SupportCombinations を追加し、同じ geometry と space 列を複数の規則で利用できることを確認する。
4. HilbertSpace のトップレベル Site の読み取り・葉の直積・位置対応を共通の接続補助へ整理する。Site 関係そのものと葉の選択を混同しない。
5. 各モデルの内部処理と設計用 Pauli string サンプルを移行する。既存モデルが任意の多項 relation をそのまま受け付けるとは限らないため、公開シグネチャまで無条件に一般化しない。
6. 参照行列・既存モデルの検証と rustdoc・examples の更新を行い、旧 API の最終整理を判断する。

主な内部区分は geometry、relation、group、connection とする。今回の v3 では Chain と多体関係の正式実装を対象とし、正方格子は拡張性の検証に使う設計用の型に留める。

## 試作で確認した結果と限界

[geometry_relations.rs](examples/geometry_relations.rs) は上記 trait と接続ビューを実装した、独立パッケージ内の設計用サンプルである。正方格子の型は、open 境界の例を試すためのものであり、製品版の格子実装ではない。

```sh
direnv exec . cargo run --offline --manifest-path docs/plans/v3/examples/Cargo.toml --bin geometry-relations
direnv exec . cargo test --offline --manifest-path docs/plans/v3/examples/Cargo.toml
```

2026-10-31、Rust 1.99.1 で統合試作8件と既存の Pauli string の4件、計12件のテストが成功した。全 target の Clippy（警告をエラー）も成功した。

新規5件では、同じ Chain・space 列に複数の規則を適用すること、同じ NearestNeighbors 型を Chain と正方格子で使うこと、x/y 方向の3 Site 列、角・辺・内部の近傍長、中心と方向の保持、借用先の同一性、不正な入力、空系、非 Clone・非 Eq の Index での space 接続を確認した。式構築では Eq のみの Index と評価能力・Clone を持たない係数型で、接続した space から既存 core の積・bind を利用できることも確認した。

葉の組み合わせを確認する追加3件では、入れ子で長さの異なる3 Site の直積・zip の期待列、空 Site・0 Site・単一 Site・終了後の next、実際の借用先、非 Clone・非 Eq の Index、独自の FnOnce を確認した。サンプルの実行でも、二つの2葉 Site から直積4組・zip 2組を出力した。

従来の37個の参照行列比較は初期 Pauli サンプルのテストであり、新しい Relation<G> と SiteMapper を用いた全 Pauli string の数値比較はまだ実施していない。将来の周期正方格子、役割を持つ関係の汎用能力 trait、性能、MSRV 1.85、専用 no_std target は未検証である。葉の二つの補助と combine_leaves も、正式ライブラリの公開 API にはまだ追加していない。
