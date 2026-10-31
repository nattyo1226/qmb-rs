# qmb-geometry のリファクタリング計画

作成日: 2026-10-31
更新日: 2026-10-31

**設計採用済み・本体未実装**。geometry と relation の採用設計、移行手順、検証条件を本書に集約する。採用した名称は `SiteMapper` と `num_sites`。旧 API の互換性は維持しない。

現行の契約は[Chain・Site と自由度の配置](../../design/geometry.md)、試作は[geometry_relations.rs](examples/geometry_relations.rs)を参照する。本書の API は試作内で確認したもので、ライブラリ本体にはまだ追加していない。

## 1. 目的と範囲

現在の幾何・空間接続は最近接の2要素ペアに固定されている。次の二つを満たす構成へ移行する。

1. 幾何を利用し、連続・隣接などの条件を満たす Site の HilbertSpace を列挙する。
2. 同じ geometry と外部の space 列に複数の relation を適用し、体数ごとの tuple trait や不要な数値能力を要求せずに式構築へ渡す。

今回の正式実装は Chain と、その最近接・連続配置、共通の借用接続を対象とする。将来の正方格子の方向付き列・中心付き近傍も追加できる設計にするが、正方格子自体の正式実装は今回に含めない。

Relation は具体的な幾何の隣接・方向・境界条件などから配置を定める。候補数だけで任意の部分集合を作る一般の組み合わせ列挙は、geometry の Relation に含めない。

## 2. 型と責務

| 型・trait | 責務 |
| --- | --- |
| `Geometry` | 共通能力として `num_sites()` を保証する。構造・境界条件は具体型が持つ |
| `Relation<G>` | geometry G に適用する規則として、Site 番号の Group を列挙する |
| `SiteGroup<CAP>` | 最大容量の配列と有効長で、異なる長さの組を同じ型で保持する |
| `SiteMapper<'g, 's, G, T>` | 一つの geometry と外部の Site 順の `&[T]` を借用し、対応を検証する |
| `MappedGroup<'s, T, Group>` | 元の Group を保持し、Site 番号を組内順の `&T` へ写す |

`T = HilbertSpace<I>` なら Site の空間を借用列挙できる。Index の発行・複製、空間の並べ替え・暗黙の tensor 合成は行わない。空間の所有と全作用空間の選択は利用側が担当する。

### Geometry と Relation

```rust,ignore
pub trait Geometry {
    fn num_sites(&self) -> usize;
}

pub trait Relation<G: Geometry> {
    type Error;

    type Group<'a>: AsRef<[usize]>
    where
        Self: 'a,
        G: 'a;

    fn groups<'a>(
        &'a self,
        geometry: &'a G,
    ) -> Result<impl Iterator<Item = Self::Group<'a>> + 'a, Self::Error>;
}
```

Chain のフィールド・幾何の責務・構築入口を維持し、Site 数の取得を `site_count()` から `num_sites()` に変更して Geometry を実装する。試作では未変更の本体に接続するため、`Geometry::num_sites()` から現行の `Chain::site_count()` に委譲している。

Relation は設定を保持し、適用時に geometry を受け取る。同じ G に複数の規則を実装でき、同じ規則型を複数の G に実装することもできる。Geometry に単一の関連 Relation 型は割り当てず、全幾何へ共通の連続性や方向型も要求しない。適用可能性は `R: Relation<G>` で表す。

`Error` は関連型とし、独自 Relation のエラー型をそのまま保持する。Debug・Display・数値能力などは一律に要求せず、診断などで必要な時点に限定する。

GAT は geometry・規則内の保存済みデータの借用にも対応する。生成中の可変バッファへの参照を next ごとに再利用する lending iterator にはせず、生成型は Group を所有値として返す。

### Group の表現

固定長の配置には `[usize; K]`、可変長には `SiteGroup<CAP>`、中心・方向が必要な場合には役割情報を持つ具体型を使える。すべての Group を一つの格納型へ統一する必要はない。

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

CAP は実際の Site 数ではなく格納上限である。構築時に `len <= CAP` を保証し、未使用部分は公開しない。組ごとのヒープ確保は不要だが、容量に比例する領域・移動コストはある。

番号範囲、順序、重複の規則は各 Relation と接続の契約で定める。SiteGroup 自体に昇順・重複除去を強制せず、順序付き配置の役割を維持する。

## 3. 具体的な Relation の列挙契約

### NearestNeighbors

Chain の無向の最近接辺を `[i, j]`（i < j）へ正規化し、辞書順・重複なしで列挙する。open は隣接位置、periodic は先頭と末尾の辺も含む。

0・1 Site は辺なし、2 Site はどちらの境界条件でも一辺とする。現在の最近接列挙の数学的な契約を維持し、自己ループを作らない。

### Consecutive<K>

各位置の役割を保つ、長さ K の順序付き連続配置を列挙する。N を `num_sites()` とする。

| 条件 | 列挙 |
| --- | --- |
| K == 0 | 不正な設定としてエラー |
| 1 <= K <= N、open | 開始位置 `0..=N-K` ごとの連続配置 |
| 1 <= K <= N、periodic | 開始位置 `0..N` ごとに N で折り返す配置 |
| K > N | 配置なし。同じ Site を繰り返して埋めない |

K == 1 は各 Site を一度ずつ返す。N == 0 かつ K >= 1 は配置なしとする。算術 overflow と0による剰余を避けて実装する。

```text
N=5, K=3
open:     [0,1,2], [1,2,3], [2,3,4]
periodic: [0,1,2], [1,2,3], [2,3,4], [3,4,0], [4,0,1]
```

組内をソートしたり、Site 集合が同じという理由で統合したりしない。例えば3 Site の周期系の `X_i Z_j X_k` は、回転ごとに中央の役割が変わる。graded な積でも順序を保つ必要がある。

NearestNeighbors と Consecutive<2> は別の規約である。2 Site の periodic では、前者は一辺、後者は `[0,1]`・`[1,0]` の二配置を返す。既存モデルの辺を後者へ置き換えて二重に数えない。

### 将来の幾何への拡張

正方格子では x/y 方向の列、中心付き近傍、plaquette などを具体的な Relation として追加できる。外部列への Site 番号付けと、幾何上の座標・方向を区別する。

中心付き近傍の Group は中心と隣接先の方向を保持する。十分大きい open 正方格子では中心込みの長さが角・辺・内部で3・4・5となるため、有効長付きの格納が使える。境界で欠けた方向を、残った要素の位置だけで推測させない。

接続後も元の Group 型を保持し、近傍用の MappedGroup に `center()`・`neighbors()` を提供できる。一般の消費側には中心・方向の能力を要求しない。近傍全体で多体項を作るか、中心と各近傍の2体項を作るかはモデル側が決める。

## 4. SiteMapper と借用接続

```rust,ignore
let chain = Chain::new(4, Boundary::Open);
let mapper = SiteMapper::new(&chain, &site_spaces)?;

for group in mapper.groups(&Consecutive::<3>)? {
    let group = group?; // MappedGroup
    for factors in group.leaves_cartesian() {
        // Vec<&HilbertSpace<I>> を使い、局所演算子を配置する。
    }
}

for group in mapper.groups(&NearestNeighbors)? {
    let group = group?;
    for factors in group.leaves_zip() {
        // 同じ geometry・space 列に別の規則を適用する。
    }
}
```

SiteMapper は geometry と外部列への参照だけを保持する。列長の一致は検証するが、利用側が意図した座標と列順の意味までは自動検証しない。

MappedGroup の `iter()` は組内順の `&T`、`source()` は元の Group を返す。要素・Index の Clone・Eq は不要であり、全件の事前収集も要求しない。結果から式を構築した後に、式が一時的な Group を借用し続ける前提にはしない。

全空間のトップレベルを Site と読む現行 `_in` の経路も、共通接続への補助として整理できる。その場合は葉入力を一 Site、空 Tensor を0 Site、空の空間を一つ包んだ Tensor を一つの空 Site と読む。全葉数から Site 数を決めず、保存された grouping を保持する。

## 5. 複数の局所自由度と leaves の組み合わせ

一つの Site に複数の自由度があれば、その合成 HilbertSpace をそのまま対応付ける。Site 0 が A⊗B、Site 1 が C⊗D なら、MappedGroup は両方の合成空間を借用する。

Site 内の列挙は既存の `children()` と `leaves()` を利用する。例えば (A⊗B)⊗C では、children は A⊗B と C、leaves は A・B・C をテンソル積順に返す。

### 標準の二つのメソッド

`MappedGroup<'s, HilbertSpace<I>, Group>` の場合だけ、次を提供する。

| メソッド | [A⊗B, C⊗D] に対する結果 |
| --- | --- |
| `leaves_cartesian()` | [A,C], [A,D], [B,C], [B,D] |
| `leaves_zip()` | [A,C], [B,D] |

Site が3個以上でも同じ入口を使う。結果は `Vec<&HilbertSpace<I>>` などの可変長列とし、位置 k を元の Group の Site 位置 k に対応させる。直積は先頭の Site を最も外側のループとし、各 Site の葉順を維持する。zip は最も短い葉列までとする。

### 共通入口 combine_leaves

```rust,ignore
pub fn combine_leaves<Out: Iterator>(
    &self,
    combine: impl FnOnce(Vec<Leaves<'s, I>>) -> Out,
) -> Out {
    combine(self.iter().map(HilbertSpace::leaves).collect())
}
```

入力は `qmb_core::space::Leaves` を Site ごとに並べた Vec であり、Site 数ごとの tuple trait は不要である。各 Group で一度呼ぶだけなので FnOnce を受け取る。標準の二つのメソッドも、この入口に直積・zip の関数を渡して実装する。

独自の組み合わせは `group.combine_leaves(my_combiner)` として指定する。Out::Item は固定せず、独自 Fn が順序や Site ごとの区切りを変える場合は、その意味も利用側で定める。専用の FactorSelection trait や追加の選択規則型は導入しない。

現行 Leaves は Clone を実装していない。試作の直積は読み取った葉への参照をキャッシュして再走査し、Index・HilbertSpace・入力 Iterator の Clone を要求しない。全結果を事前生成せず、全件消費時のキャッシュは入力の葉数の和に比例する。正式実装時のキャッシュと再生成の比較は実装判断とする。

| 入力 | 直積・zip の扱い |
| --- | --- |
| 空 Site を含む | 結果なし。位置を詰めて接続し直さない |
| Site の組自体が空 | 空の葉列を一つ返す。zip で無限生成しない |
| 一つの Site | その Site の葉を一つずつ返す |

標準の二つは有効な参照の組み合わせであり、新しいエラー型を加えない。独自 Fn が必要とする検証はその Fn が扱う。

## 6. エラーと能力制約

| 時点 | 検証とエラー |
| --- | --- |
| geometry の構築 | Site 数の overflow など、幾何自体の不整合 |
| Relation の適用 | 長さ0などの不正設定を列挙前に検出し、R::Error を返す |
| SiteMapper の構築 | Relation が空でも、geometry の Site 数と外部列長を検証 |
| 各 Group の接続 | 独自 Relation の範囲外番号を、要素アクセス前に検出 |
| モデル・core の空間構築 | Index の一意性、配置可能性、全作用空間への bind |

mapper の `groups()` の外側の Result は R::Error、各 Group の接続結果は mapping 側のエラーを返す。独自エラーへの Debug・Display・Error の実装は必須にしない。公開エラー型名・診断項目は実装時に整理する。

列挙は遅延処理を基本とし、エラーまでに消費済みの結果を巻き戻さない。関係の読み取り・借用接続へ数値演算・評価能力を要求せず、式構築に必要な Index の Eq などは、その処理の段階で要求する。

## 7. 式構築・Pauli string との責務境界

既存の `ExpressionFactory::operator_product_at` と `try_operator_product` は任意個数の配置・積を構築できる。体数ごとの構築処理や、geometry 固有の評価器を追加しない。

X/Y/Z の選択・規格化と全作用空間の選択は具体的代数・モデル側、配置・積・bind・評価は qmb-core が担当する。core から geometry への依存は追加せず、geometry は HilbertSpace 接続のための現行 core 依存と `no_std + alloc` を維持する。

一般の支持集合の全列挙は geometry の Relation ではない。[pauli_strings.rs](examples/pauli_strings.rs) の SupportCombinations は通常の Iterator の補助であり、標準 Relation・今回の正式実装には含めない。統合試作の TestGroups も、不正な番号などを流す test-only fixture である。

Pauli サンプルは、1 Site が1 spin の場合に、全支持集合への非恒等な X/Y/Z の割り当てを列挙する。3 Site・高々2 Site なら恒等項を含む37個となる。Sx/Sy/Sz は Pauli 行列の半分なので各因子を2倍し、全作用空間へ bind する。構築関数の制約は `C: From<f64>, I: Eq, Group: AsRef<[usize]>` と、Group を要素に持つ IntoIterator であり、評価能力・係数の Clone は要求しない。

「高々 K Site」の K と、実際に作用する葉の数は区別する。1 Site に m 個の spin-1/2 があれば、その Site 上の非恒等 Pauli string は `4^m - 1` 種類である。2葉では XI・IX・XX などを含む15種類となり、各 Site から葉を一つずつ選ぶだけでは全項を作れない。

この全割り当ては利用側の Iterator や独自 combiner で扱い、geometry の標準操作を増やさない。既存 Pauli サンプルの複数葉への全割り当ては未実装である。

## 8. 実装・移行手順

1. Geometry・Relation<G>・SiteGroup・SiteMapper・MappedGroup を正式実装する。Chain の構築入口と幾何の責務を維持し、Site 数の取得を num_sites() に統一する。
2. NearestNeighbors と Consecutive を実装し、同じ geometry・space 列へ複数の規則を適用する経路を完成させる。
3. 葉の直積・zip・combine_leaves と、トップレベルを Site と読む接続補助を整理する。
4. spin-half・spin-one・fermion のモデル、examples、テストを移行する。旧 nn_* と site_count() は削除し、非推奨エイリアスは残さない。
5. 設計用サンプルの重複した実装を正式 API の利用へ置き換え、rustdoc・現行設計文書を更新する。
6. 以下の受入確認を行い、旧 API と重複した列挙処理の残存を確認する。

既存モデルの公開シグネチャまで無条件に任意の多項 Relation へ一般化しない。モデル固有の葉の役割、結合係数、全作用空間、符号規約、数値結果を維持する。辺なし・結合0の場合も入力検証を省略しない。

内部は geometry・relation・group・connection に責務を分け、公開窓口は crate root を基本とする。ファイル分割、constructor の細部、Iterator の具体型、キャッシュ、追加の読み取り補助は、採用済みの契約を保つ範囲で実装時に判断する。追加の合意を要する大きな未決事項はない。

## 9. 検証計画と試作の記録

### 正式実装の受入確認

- 最近接は0・1・2 Site とそれ以上の両境界条件、連続配置は3体・4体を含む K<N・K=N・K>N、open・periodic、N=0・K=0 の期待列で確認する。順序付き配置と無向辺の規約を区別する。
- 複数葉・入れ子・空 Site・可変長 Group について、直積・zip の期待列、終了後の next、借用先の同一性を確認する。非 Clone・非 Eq の Index に不要な制約を加えない。
- 中心・方向などの役割情報、独自 Relation のエラー、列長の先行検証と範囲外番号の逐次検証、空入力を確認する。
- 非対称な演算子配置と graded な入力を使い、葉順・周期境界の役割・符号を維持する。多体項の手作業構築と比較する。
- 各既存モデルの独立な参照計算・行列、Pauli の参照行列比較を回帰確認に使う。構築だけの例に評価用 bound を持ち込まない。
- workspace の該当テスト・format・Clippy・ドキュメント生成、MSRV 1.85、専用 no_std target と必要な feature 構成を確認する。

### 保存した試作と既存の検証結果

[独立した Cargo パッケージ](examples/Cargo.toml)に、次のサンプルを保持する。

| サンプル | 確認対象 |
| --- | --- |
| [geometry_relations.rs](examples/geometry_relations.rs) | 採用した trait・mapper、Chain と設計用の open 正方格子、葉の組み合わせ |
| [pauli_strings.rs](examples/pauli_strings.rs) | geometry と独立した支持集合 Iterator からの Pauli string 構築 |

```sh
direnv exec . cargo run --offline --manifest-path docs/plans/v3/examples/Cargo.toml --bin geometry-relations
direnv exec . cargo run --offline --manifest-path docs/plans/v3/examples/Cargo.toml --bin pauli-strings
direnv exec . cargo test --offline --manifest-path docs/plans/v3/examples/Cargo.toml
direnv exec . cargo clippy --offline --manifest-path docs/plans/v3/examples/Cargo.toml --all-targets -- -D warnings
```

2026-10-31、Rust 1.99.1 で統合試作9件と Pauli サンプル4件、計13件のテスト、format、全 target の Clippy が成功した。これは既存の試作検証記録であり、今回の文書統合時にテストを再実行したものではない。

統合試作では複数 Relation の適用、方向付き列、中心・近傍長の保持、直積・zip、空入力、借用、独自 FnOnce と独自エラーを確認した。Pauli サンプルでは37個すべてを独立な Pauli 行列の tensor product と比較し、重複・規格化・全作用空間を確認した。

新しい Relation<G>・SiteMapper を通した全 Pauli string の数値比較、複数葉への全割り当て、周期正方格子、性能、MSRV 1.85、専用 no_std target は未検証である。正式ライブラリへの API 追加もまだ行っていない。
