# Chain・Site と自由度の配置

更新日: 2026-11-01

## 1. 対象と責任分担

`qmb-geometry` は幾何、Site の配置規則、外部列への借用接続を分離する。公開型は crate root から利用する。core は geometry に依存せず、geometry は HilbertSpace 接続のため core に依存する。`no_std + alloc` を維持する。

| 型・trait | 責務 |
| --- | --- |
| `Geometry` | `num_sites()` を保証する最小の共通契約 |
| `Chain` | Site 数と境界条件を保持する不変な幾何 |
| `Relation<G>` | G の隣接・方向・境界条件などから Site 番号の Group を列挙 |
| `SiteGroup<CAP>` | 最大容量の配列と有効長。`AsRef<[usize]>` は使用部分だけを公開 |
| `SiteMapper<G,T>` | geometry と Site 順の `&[T]` を借用し、対応付けを検証 |
| `MappedGroup<T,Group>` | 元の Group と外部列を保持し、組内順の `&T` を列挙 |

`Chain::new(num_sites, boundary)` は検証エラーなしで Self を返し、0 Site も認める。`num_sites()` と `boundary()` で保持内容を取得する。Index・HilbertSpace・局所代数・係数型は保持・発行しない。Site は `0..num_sites()` の位置番号であり、専用の SiteId や所属トークンは設けない。

同じ geometry へ複数の Relation を適用でき、同じ規則型を複数の geometry に実装できる。共通の Group 契約は `AsRef<[usize]>` だけとし、固定長配列・可変長・保存済み slice・中心や方向を持つ独自型を使える。Group の順序や重複の意味は各 Relation が定める。SiteGroup 自体はソート・重複除去を行わない。

`SiteGroup<CAP>` の CAP は格納上限であり、実際の Site 数は有効長で表す。`new` は容量超過を `SiteGroupError::CapacityExceeded { capacity, len }` として拒否し、未使用部分を公開しない。組ごとのヒープ確保は不要だが、容量に比例した領域・移動コストはある。

v3 では旧 `nn_*`・`site_count()`・`ChainError` を撤去し、各モデルの接続エラーを `UtilityError::Mapping` に統一した。互換エイリアスは残さない。実装範囲は共通の trait・借用接続と Chain の最近接・連続配置である。正方格子は [geometry_relations.rs](../../examples/src/bin/geometry_relations.rs) 内の使用例に留め、一般の支持集合の組み合わせ列挙は標準 Relation に含めない。

## 2. 関係の列挙と借用接続

`Relation<G>::groups(&geometry)` は独自の関連型 `Error` を返し、成功時は Group を遅延列挙する。設定の検証は Iterator の返却前に行う。エラー型に Debug・Display・Error の実装を一律には要求しない。GAT により規則や geometry 内の保存済みデータも借用できるが、生成 Iterator の可変バッファを貸し出す形式ではない。

| Chain の Relation | 列挙規約 |
| --- | --- |
| `NearestNeighbors` | `[i,j]`（i<j）の辞書順・一意な無向辺。0・1 Site は辺なし、2 Site は両境界で一辺 |
| `Consecutive<K>` | 役割順を保持する K Site の連続配置。open は開始位置 `0..=N-K`、periodic は `0..N` で折り返す |

`Consecutive<0>` は `ConsecutiveError::ZeroLength`、K>N は配置なしとする。K=1 は各 Site を一度ずつ返す。配置内の並べ替えや、同じ支持集合になる回転の統合は行わない。2 Site の periodic で `Consecutive<2>` は `[0,1]`・`[1,0]` を返し、最近接の一辺とは区別する。NearestNeighbors の Chain 実装のエラー型は Infallible である。

N=0・K>=1 は配置なしとし、K>N の場合も同じ Site を繰り返して埋めない。算術 overflow と0による剰余を避ける。例えば N=5・K=3 は次の順序になる。

```text
open:     [0,1,2], [1,2,3], [2,3,4]
periodic: [0,1,2], [1,2,3], [2,3,4], [3,4,0], [4,0,1]
```

組内の位置は局所演算子の役割に対応する。回転した支持集合が同じでも、`X_i Z_j X_k` の中央の役割や graded な積の符号を保つために統合しない。相互作用の対称性による重複の扱いはモデル側が決める。

```rust,ignore
use qmb_geometry::{Boundary, Chain, Consecutive, NearestNeighbors, SiteMapper};
let chain = Chain::new(4, Boundary::Open);
let mapper = SiteMapper::new(&chain, &site_spaces)?;
for group in mapper.groups(&Consecutive::<3>)? {
    let group = group?;
    for leaves in group.leaves_cartesian() {
        // 各 Site から葉を一つずつ選び、既存 core で配置・積を構築する。
    }
}
for group in mapper.groups(&NearestNeighbors)? {
    let group = group?;
    // 同じ geometry と外部列に別の Relation を適用できる。
}
```

`SiteMapper::new` は列長の完全一致を Relation の適用前に検証し、空の列挙でも省略しない。不一致は `MappingError::SiteCountMismatch { expected, actual }`。`groups()` の外側の Result は `R::Error`、各 Group の接続結果は `Result<MappedGroup<...>, MappingError>` とする。番号範囲は各 Group を消費する時点で検証し、不正なら `SiteOutOfRange { site, num_sites }` を返す。全件の事前収集や消費済みの結果の巻き戻しは行わない。

geometry の Site 数・番号の意味は mapper の借用中に変えない。列長が一致しても、利用側が意図した座標と列順の意味までは自動検証しない。

MappedGroup の `iter()` は組内順の `&T`、`get(position)` は組内位置の要素、`source()` は元の Group を返す。中心・方向なども保持されるため、独自 Group 用の拡張 trait から読み取れる。要素・Index の Clone・Eq は接続に不要であり、Index 重複は空間・モデル構築側で検証する。

## 3. 空間と局所自由度

局所空間列 `&[HilbertSpace<I>]` の一要素を一 Site に対応させる。葉・複数葉・入れ子・単項 Tensor・空の空間を使える。利用側は既存の空間ハンドルを共有して列を構成でき、Index の Clone・再発行は不要である。

`SiteMapper::from_space` は保存されたトップレベル分割を Site と読む。葉入力は一 Site、空 Tensor は0 Site、`Tensor([unit])` は一つの空 Site、単項 Tensor は一 Site とする。全葉数から Site 数を決めず、局所 grouping を保持する。

HilbertSpace 用 MappedGroup は次の三つの入口を持つ。

| メソッド | `[A⊗B,C⊗D]` に対する結果 |
| --- | --- |
| `leaves_cartesian()` | `[A,C], [A,D], [B,C], [B,D]` |
| `leaves_zip()` | `[A,C], [B,D]` |
| `combine_leaves(FnOnce)` | Site ごとの `Vec<Leaves>` を一つの Iterator にまとめる独自処理 |

任意の Site 数で同じ入口を使い、体数ごとの tuple trait を要求しない。標準の二つは `combine_leaves` を利用する。直積は先頭 Site を最も外側とし、読み取った葉への参照をキャッシュして再走査する。全結果の事前生成や Leaves・Index の Clone は不要である。zip は最短の葉列までとし、いずれも Site 内のテンソル積順を保つ。

`combine_leaves` は Site ごとの `Vec<qmb_core::space::Leaves<'s, I>>` を FnOnce に一度だけ渡し、返された Iterator をそのまま返す。出力 Item は固定せず、専用の選択 trait は導入しない。直積を全件消費した場合の参照キャッシュは、入力の葉数の和に比例する。

空 Site を含めば結果なしとし、位置を詰めてつなぎ直さない。Site の組自体が空なら両方とも空の葉列を一つ返す。一つの Site はその葉を一つずつ返す。独自 Fn が順序や対応を変える場合の契約は利用側が定める。

Site 順序と別途存在する全空間の因子順序を区別する。接続は空間・状態・行列を暗黙に並べ替えず、tensor 合成もしない。全作用空間と X/Y/Z などの選択・規格化は利用側、配置・積・bind は既存 core が担当する。各 Site から葉を一つ選ぶ操作は、複数葉に対する全 Pauli string の列挙とは区別する。

既存の `ExpressionFactory::operator_product_at`・`try_operator_product` は任意長の配置・積を受け取れる。体数別の構築 API や geometry 固有の評価器は追加せず、列挙・接続に評価用の数値能力を要求しない。

全 Pauli string の使用例は [pauli_strings.rs](../../examples/src/bin/pauli_strings.rs) に置く。1 Site が1 spin-1/2 の場合、通常の Iterator で大きさ `r = 0..=min(K,N)` の支持集合を一度ずつ選び、各集合へ `3^r` 通りの X/Y/Z を割り当てる。項数は `Σ(r=0..min(K,N)) binomial(N,r) × 3^r` となり、3 Site・高々2 Site では恒等項を含む37項になる。空の支持集合は全作用空間上の恒等項一つであり、`Consecutive<0>` とは別の契約である。Sx/Sy/Sz は Pauli 行列をそのまま表すため、追加の倍率を掛けずに全作用空間へ bind する。構築に必要な制約は `I: Eq, Group: AsRef<[usize]>` と Group の IntoIterator に留める。

K は作用する Site 数であり、実際に作用する葉の数とは区別する。1 Site に m 個の spin-1/2 があれば、その Site 上の非恒等 Pauli string は `4^m - 1` 種類となる。複数葉への全割り当ては利用側の Iterator や独自 combiner で扱い、geometry の標準操作は増やさない。この全割り当ては現在の使用例には含めていない。

## 4. SpinHalf のモデル構築

モデル構築の公開入口は `qmb_spin_half::models::SpinHalfModels<SpinHalf<C>>` の関連関数である。以下では `type Models = SpinHalfModels<Spin>` とする。旧 `utility` の自由関数は残さず、エラー型 `UtilityError` を crate root から公開する。

| 呼び出し | 成功時の戻り値 |
| --- | --- |
| `Models::xyz(&chain, &factors, jx, jy, jz)` | `OperatorOf<Spin, I>` |
| `Models::heisenberg(&chain, &factors, coupling)` | `OperatorOf<Spin, I>` |
| `Models::neel(&chain, &factors, first)` | `StateOf<Spin, I>` |

全関数は Result を返し、列長と全葉の Index 重複を、辺がなくても検証する。列長の不一致は `UtilityError::Mapping(MappingError::SiteCountMismatch { expected, actual })`、重複は `UtilityError::Space(SpaceError::DuplicateIndex)`、配置・束縛エラーは `UtilityError::Expression` とし、元エラーを source として保持する。Index には Eq のみを要求する。

XYZ は `−(1/4)Σ_(i,j) Σ_(a∈Li,b∈Lj) (Jx Sx_a Sx_b + Jy Sy_a Sy_b + Jz Sz_a Sz_b)` を構築する。Li は Site i の葉集合であり、`NearestNeighbors` と `leaves_cartesian()` の全組み合わせを使う。Site 内の葉間結合は追加しない。Heisenberg は三結合を J に揃え、正の J を強磁性的とする。標準モデルは位置対応を暗黙に選ばない。

XYZ・Heisenberg の全作用空間は、入力の局所空間ハンドルを共有し、列順に tensor した空間である。局所 grouping と相互作用に現れない Site も保持し、返す演算子をその空間に bind する。別途存在する全空間の root・上位 grouping・異なる因子順序は保存しない。

Neel は Site 位置と first から向きを交互に決め、その Site の全葉に同じ向きを置く。空 Site も偶奇に含める。状態の因子順は Site 順に局所空間の葉をたどった順とし、奇数長 periodic も認める。例えば `[(a,b), c]`、first=Up は `Up_a ⊗ Up_b ⊗ Down_c` となる。

0 Site は単位空間に束縛した zero Hamiltonian と単位状態を返す。全 Site が空でも同様であり、Hamiltonian の空間には空 Site の grouping を残す。

数値能力は個別関数へ限定する。xyz は `C: Clone + From<R> + Neg<Output=C> + One + Add<Output=C> + Div<Output=C>`、heisenberg はさらに `R: Clone` を要求する。neel の C に Clone・数値演算は不要である。モデル固有関数は QuantumSystem trait に追加しない。

## 5. 検証対象

geometry のテストは 0・1・2・3 以上の open／periodic、列長過不足、非 Clone・非 Eq 要素、借用の同一性、入れ子・空・可変長局所空間の順序、全組み合わせと位置対応、トップレベル分割の解釈と `new`／`from_space` の一致、3体・4体の連続配置、独自エラー・役割情報と遅延範囲検証を検証する。

spin-half のテストは複数葉の全組み合わせ結合、Site ごとの向き、空 Site の偶奇と辺、全空・0 Site、重複と列長拒否、非 Clone Index、深い grouping、最小の係数能力を検証する。examples の縦断テストは両精度の独立な Pauli 行列との一致、符号・norm・期待値、および密疎 export を検証する。

### 5.1. テストと使用例

| 対象 | テスト・使用例 |
| --- | --- |
| 最近接・多体配置の列挙規約、容量制限、保存済み Group の借用、独自エラーと遅延範囲検証 | [geometry.rs](../../crates/qmb-geometry/tests/geometry.rs) |
| 最大 Site 数での折り返し算術 | [consecutive.rs の単体テスト](../../crates/qmb-geometry/src/relation/consecutive.rs) |
| 入れ子・複数葉・可変長 Group、空入力、終了後の next、借用の同一性、FnOnce、トップレベル解釈 | [space.rs](../../crates/qmb-geometry/tests/space.rs) |
| 3体・4体の Pauli 配置、Pauli string 37項、周期境界の fermion 積を手作業構築・独立な行列と比較 | [geometry_relations.rs](../../examples/tests/geometry_relations.rs) |
| 複数 Relation、異なる geometry、方向付き列、中心・方向を持つ可変長近傍、最小限の型制約 | [geometry_relations.rs のサンプル](../../examples/src/bin/geometry_relations.rs) |
| 支持集合の Iterator と Pauli の規格化・独立行列 | [pauli_strings.rs のサンプル](../../examples/src/bin/pauli_strings.rs) |

サンプルは workspace の examples パッケージに含める。実行方法と feature の選択は [examples の README](../../examples/README.md)を参照する。

### 5.2. 実装・移設時の検証記録

2026-10-31、Rust 1.99.1 の direnv devShell で geometry・relation の実装と受入確認を完了した。workspace 全 feature・default feature 無効のテスト、最終状態の geometry テスト13件と doc test、多体項・Pauli string・graded 積の数値比較3件、format、全 target の Clippy、警告をエラーとする rustdoc の生成が成功した。全 feature の workspace 一括実行後に追加した境界値テストは、geometry 単独と default feature 無効の workspace 全体でも確認した。

examples への移設時には両サンプルを実行し、feature 無効時のテスト11件・全 feature 有効時のテスト13件、format、両 feature 構成の全 target の Clippy が成功した。Pauli の行列比較を含む2件は `export` feature を有効にする。これらは実装・移設時の記録であり、今回の文書整理でテストを再実行したものではない。

```sh
direnv exec . cargo test -p qmb-geometry --all-features --offline
direnv exec . cargo test -p examples --all-features --offline --test geometry_relations
direnv exec . cargo test -p examples --all-features --offline --bin geometry_relations --bin pauli_strings
direnv exec . cargo test --workspace --all-features --offline -j2
direnv exec . cargo test --workspace --no-default-features --offline -j2
direnv exec . cargo fmt --all --check
direnv exec . cargo clippy --workspace --all-targets --all-features --offline -- -D warnings
direnv exec . sh -c 'RUSTDOCFLAGS="-D warnings" cargo doc --workspace --all-features --no-deps --offline'
```

Rust 1.85 と専用 no_std target の実ビルドは利用者の指定により省略した。専用 target の試行は依存解決とネットワーク制限で完了しておらず、成功した検証には数えない。性能測定・キャッシュと再生成の比較、周期正方格子の正式実装も未実施である。共通の未実施事項は[今後の課題](../future-work.md#3-未実施の検証と性能)にまとめる。

## 6. Fermion の Hubbard モデル

標準 Hubbard と symmetric Hubbard の二種類を実装した。以下の Hamiltonian は合意した相互作用の違いを具体化したものであり、API と配置入力は既存 spin-half に合わせる。

### 6.1. Hamiltonian の規約

一様な実パラメータ t・U、Chain の最近接辺、spin を保存するホッピングを対象とする。両モデルとも化学ポテンシャル項を含めない。

```math
K(t)=-t\sum_{\langle i,j\rangle,\sigma}
  \left(f_{i\sigma}^\dagger f_{j\sigma}
       +f_{j\sigma}^\dagger f_{i\sigma}\right),\qquad
N=\sum_{i,\sigma}n_{i\sigma},
```

```math
H_{\mathrm{Hubbard}}(t,U)
 =K(t)+U\sum_i n_{i\uparrow}n_{i\downarrow},
```

```math
H_{\mathrm{symmetric}}(t,U)
 =K(t)+U\sum_i\left(n_{i\uparrow}-\frac I2\right)
                  \left(n_{i\downarrow}-\frac I2\right).
```

標準形は [ALPS の Hubbard 定義](https://alps.comp-phys.org/documentation/models/hubbard/)、symmetric 形は [TRIQS/TPRF の定義](https://triqs.github.io/tprf/3.1.x/user_guide/Linearized%20Eliashberg%20equation%20on%20the%20attractive%20Hubbard%20model.html) のホッピング・相互作用と同じ符号規約である。化学ポテンシャル引数は設けない。

symmetric ではホッピングを変更せず、相互作用の展開で生じる定数項を落とさない。L Site に対して、展開から次の恒等式が得られる。

```math
H_{\mathrm{symmetric}}(t,U)
 =H_{\mathrm{Hubbard}}(t,U)
   -\frac U2 N+\frac{UL}{4}I.
```

`n−I/2=−Mz/2` なので、symmetric な局所相互作用は `(U/4) Mz↑Mz↓` とも表せる。両表現の値の一致を検証する。

### 6.2. 公開入口と配置

既存の `SpinHalf<C>` と同様に `models.rs` で構築し、`models::FermionModels<Fermion<C>>` の関連関数として公開する。QuantumSystem trait にモデル固有メソッドを追加しない。

```rust,ignore
FermionModels::<Fermion<C>>::hubbard(&chain, &sites, t, u)
FermionModels::<Fermion<C>>::symmetric_hubbard(&chain, &sites, t, u)
// いずれも Result<OperatorOf<Fermion<C>, I>, UtilityError> を返す。
```

- `sites: &[HilbertSpace<I>]` の一要素を一つの物理 Site とし、各 Site は葉順 `(up,down)` の二 mode の graded 空間とする。Index は利用側で発行し、関数内で再発行しない。
- 全作用空間は Site 順、各 Site 内の up・down 順で graded 合成し、入力の Site 空間を共有して構築する。返す演算子はこの全空間に bind する。
- ホッピングは隣接 Site の同じ spin 同士だけを結ぶ。spin-half の XYZ が使う全葉ペアの直積は用いない。各辺に両方向のホッピングを一度ずつ加え、符号は core の graded 合成に委ねる。
- open / periodic の辺は既存 Chain の契約に従う。0・1 Site は辺なし、2 Site の periodic でも辺は一つとし、ホッピングを重複させない。
- Site 数の一致、各 Site が二つの葉を持つ graded 空間であること、全 mode の Index 重複を、結合値が0でも検証する。不正な Site 構成は `UtilityError::InvalidSite { site }` とし、列長・空間・式のエラーはそれぞれ `Mapping`・`Space`・`Expression` variant に元エラーを保持する。
- 0 Site は graded な単位空間に束縛した zero Hamiltonian を返す。1 Site はホッピングなしで相互作用を保持する。
- standard Hubbard の係数 bound は `C: Clone + From<R> + Neg<Output=C>`、symmetric Hubbard はさらに `One + Add<Output=C> + Div<Output=C>` を要求する。入力 t・U は同じ R とする。symmetric 相互作用は `(U/4) Mz↑Mz↓` を使い、1/4を係数型の中で生成する。虚数単位・局所 Product・Act・Float は要求しない。t・U が実数に対応すること、分数を数学的に表せることを契約とする。

局所演算子・Majorana の規格化・基底・spinful の位相規約は [Fermion の局所代数](local-structure.md#9-fermion-の局所代数)に従う。係数の中で1/4を生成して `(U/4) Mz↑Mz↓` を構築する実装は [models.rs](../../crates/qmb-fermion/src/models.rs) に置く。

### 6.3. 検証

- 二 mode の CAR、spinful の二重占有からの消滅符号、および全空間の mode 順を確認する。
- 1 Site の表示順 `[empty,up,down,double]` で、標準形の対角成分が `[0,0,0,U]`、symmetric 形が `[U/4,−U/4,−U/4,U/4]` となることを確認する。実際の列挙順との対応は明示する。
- 2 Site（4 mode、16状態）の両モデルを、core の評価・export を使わない Fock 空間の符号規則から構築した参照行列と比較する。実パラメータでの Hermite 性、粒子数・各 spin の粒子数の保存、`t=0`・`U=0` の極限を確認する。
- 両モデルの差が `−(U/2)N+(UL/4)I` であることと、symmetric 相互作用の `(U/4) Mz↑Mz↓` による構築を、定数項を含めて比較する。
- 0・1・2 Site と3 Site 以上の open / periodic、Site 数不一致、mode 数の過不足、ordinary な Site、Index 重複、非 Clone Index を検証する。
- `Complex<f32>` / `Complex<f64>` と Vec・ndarray・sprs の同じ基底順で値を比較し、構築の係数 bound に評価能力が混入していないことも確認する。

[モデルのテスト](../../crates/qmb-fermion/tests/models.rs)で局所規則と独立した Fock 空間の参照行列を使う。各辺に一度だけホッピングを追加すること、同じ spin のみが接続すること、symmetric 相互作用の定数項を保存することを確認する。実装時の実行結果・コマンド・未実施事項は [受入確認の記録](local-structure.md#96-実装時の検証記録)にまとめる。

## 7. SpinOne のモデル構築

`qmb_spin_one::models::SpinOneModels<SpinOne<C>>` は XYZ・Heisenberg・BLBQ の3モデルを関連関数として提供する。局所 primitive・基底・数値規約は [Spin-one の局所代数](local-structure.md#10-spin-one-の局所代数)に従い、構築は [models.rs](../../crates/qmb-spin-one/src/models.rs) に置く。

### 7.1. Hamiltonian の規約

`hbar=1` とし、まず1 Site が1 spin の場合を示す。和は Chain の最近接結合を各1回数える。primitive は `S=√(3/2) J_physical` なので、係数を2/3・4/9で補正し、結合定数の従来の物理的な意味を維持する。

```math
H_{\mathrm{XYZ}}=-\frac23\sum_{\langle ij\rangle}
\left(J_xS_i^xS_j^x+J_yS_i^yS_j^y+J_zS_i^zS_j^z\right),
```

```math
H_{\mathrm{Heisenberg}}=-\frac{2J}{3}\sum_{\langle ij\rangle}\mathbf S_i\cdot\mathbf S_j,
```

```math
H_{\mathrm{BLBQ}}=-\sum_{\langle ij\rangle}
\left[\frac{2J}{3}\,\mathbf S_i\cdot\mathbf S_j
+\frac{4K}{9}\,(\mathbf S_i\cdot\mathbf S_j)^2\right].
```

符号は既存 SpinHalf のモデルにそろえる。J>0 の bilinear 結合は強磁性的、J<0 は反強磁性的となる。`heisenberg(j)` は `xyz(j,j,j)`、`blbq(j,0)` は `heisenberg(j)` と数学的に同じである。単イオン異方性・磁場項・追加の定数シフトは設けず、上の式のエネルギー原点を保持する。

BLBQ は bilinear 結合 J と biquadratic 結合 K を独立に受け取る。角度表示 `H=sum[cos(theta) A_ij+sin(theta) A_ij²]`、`A_ij=J_physical,i·J_physical,j=(2/3)S_i·S_j` は `blbq(-cos(theta),-sin(theta))` で表せる。角度専用の入口は設けない。Hermite Hamiltonian に対して実数値の結合を渡す契約とし、型による実数制限や実行時の判定は追加しない。

### 7.2. 公開入口と空間

| 公開入口 | 成功時の戻り値 |
| --- | --- |
| `SpinOneModels::<SpinOne<C>>::xyz(&chain, &sites, jx, jy, jz)` | `OperatorOf<SpinOne<C>,I>` |
| `SpinOneModels::<SpinOne<C>>::heisenberg(&chain, &sites, j)` | `OperatorOf<SpinOne<C>,I>` |
| `SpinOneModels::<SpinOne<C>>::blbq(&chain, &sites, j, k)` | `OperatorOf<SpinOne<C>,I>` |

全関数は `&Chain` と `&[HilbertSpace<I>]` を受け取り、`Result<_,qmb_spin_one::UtilityError>` を返す。Site 列・全葉ペア・束縛の基本契約は第4節の SpinHalf の XYZ / Heisenberg とそろえる。

- 一つの入力空間を一 Site とし、隣接 Site 間の全葉ペア `(a,b)` に相互作用を加える。BLBQ の二乗は各葉ペアの `A_ab=S_a·S_b` に個別に適用し、Site 内の結合や異なる葉ペア間の交差項を加えない。
- 列長と全 Site にわたる Index の一意性を、辺の有無・結合定数によらず検証する。Index には Eq のみを要求し、Clone・Hash・Ord は要求しない。
- 入力 Site の空間ハンドルを列順に tensor した全空間へ Hamiltonian を bind し、内部 grouping・葉順・相互作用に現れない Site を保持する。既存の別の全空間の上位 grouping や因子順は復元しない。
- 空 Site は位置を保ち、前後を直接つなぎ直さない。0 Site は単位空間上の zero Hamiltonian、1 Site と全空 Site は入力の空間構造に束縛した zero Hamiltonian を返す。
- 境界条件と辺は第2節の Chain に従う。0・1 Site は辺なし、2 Site は open / periodic とも一辺とする。

エラーは `Mapping(MappingError)`・`Space(SpaceError)`・`Expression(ExpressionError)` の3 variant とし、元エラーを source として保持する。列長不一致は SiteCountMismatch、重複 Index は DuplicateIndex とし、UtilityError を crate root から公開する。

### 7.3. 遅延構築と係数能力

XYZ は配置した Cartesian 成分の積と和、Heisenberg は XYZ の等方的な場合として構築する。BLBQ は各葉ペアの `A_ab=Sx_a Sx_b+Sy_a Sy_b+Sz_a Sz_b` を一度構築し、式を共有して `(&A_ab * &A_ab)` を作る。J・K は2/3・4/9で補正して最後に各項へ適用し、J の符号や大きさを二乗へ混入させない。

| 関数 | 構築に必要な bound |
| --- | --- |
| `xyz<I: Eq,R>` | `C: Clone + From<R> + Neg<Output=C> + One + Add<Output=C> + Div<Output=C>` |
| `heisenberg<I: Eq,R: Clone>` | `C: Clone + From<R> + Neg<Output=C> + One + Add<Output=C> + Div<Output=C>` |
| `blbq<I: Eq,R>` | `C: Clone + From<R> + Neg<Output=C> + One + Add<Output=C> + Div<Output=C>` |

構築時に係数 C の中で有理数の補正を行う。局所規則の評価・基底展開・密行列・数値 backend は使わず、C の Float・Zero・共役を要求しない。評価・trace・export 時には既存 core と局所 primitive の操作別の bound が適用される。モデル固有メソッドを QuantumSystem trait に追加しない。

### 7.4. 検証と実装時の記録

[models.rs](../../crates/qmb-spin-one/tests/models.rs) では、直接記述した3×3スピン行列から参照行列を構成し、モデル実装や Chain の辺列挙を参照計算に使わない。

- `Complex<f32>` / `Complex<f64>` の両方で0〜3 Site・両境界の全行列成分・Hermite 性・trace を比較する。異方的結合、正負の J・K、zero 結合、pure biquadratic 項、等方的・bilinear 極限を含める。
- 二 spin の全スピン0・1・2に対する `A_ij` の固有値 -2・-1・1と、BLBQ のエネルギー `2J-4K`・`J-K`・`-J-K` を、独立した状態の期待値と作用の残差で検証する。
- 複数葉 Site の全葉ペアとペアごとの二乗、空 Site、0 / 1 Site、2 Site の周期辺の重複防止、非連番の因子順、grouping の共有と spectator の保持を確認する。
- 辺なし・zero 結合でも列長不一致と重複 Index を拒否すること、エラーの source、Clone・From・Neg だけの係数と非 Clone の結合入力・Index で構築できることを確認する。

2026-10-30、Rust 1.99.1 の direnv devShell で `cargo test -p qmb-spin-one --test models` の9件が成功した。これは実装時の記録であり、文書整理時の再実行結果ではない。workspace 全体を含む確認コマンドと未実施事項は [Spin-one の検証記録](local-structure.md#106-実装時の検証記録)を参照する。

利用例 [spin_one.rs](../../examples/src/bin/spin_one.rs) は `--no-default-features` と `--all-features` で実行した。`|+,->` の期待値は J=1 の等方的 XYZ / Heisenberg で1、J=1・K=0.25 の BLBQ で0.5となり、export 有効時には9×9の Vec / sprs 行列が一致した。
