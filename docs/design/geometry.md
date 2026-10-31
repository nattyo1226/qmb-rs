# Chain・Site と自由度の配置

更新日: 2026-10-31

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

v3 では旧 `nn_*`・`site_count()`・`ChainError` を撤去した。移行設計と検証記録は[リファクタリング計画](../plans/v3/geometry-refactoring.md)に記録する。正方格子は設計用サンプルに留め、一般の支持集合の組み合わせ列挙は標準 Relation に含めない。

## 2. 関係の列挙と借用接続

`Relation<G>::groups(&geometry)` は独自の関連型 `Error` を返し、成功時は Group を遅延列挙する。設定の検証は Iterator の返却前に行う。エラー型に Debug・Display・Error の実装を一律には要求しない。GAT により規則や geometry 内の保存済みデータも借用できるが、生成 Iterator の可変バッファを貸し出す形式ではない。

| Chain の Relation | 列挙規約 |
| --- | --- |
| `NearestNeighbors` | `[i,j]`（i<j）の辞書順・一意な無向辺。0・1 Site は辺なし、2 Site は両境界で一辺 |
| `Consecutive<K>` | 役割順を保持する K Site の連続配置。open は開始位置 `0..=N-K`、periodic は `0..N` で折り返す |

`Consecutive<0>` は `ConsecutiveError::ZeroLength`、K>N は配置なしとする。K=1 は各 Site を一度ずつ返す。配置内の並べ替えや、同じ支持集合になる回転の統合は行わない。2 Site の periodic で `Consecutive<2>` は `[0,1]`・`[1,0]` を返し、最近接の一辺とは区別する。NearestNeighbors の Chain 実装のエラー型は Infallible である。

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

空 Site を含めば結果なしとし、位置を詰めてつなぎ直さない。Site の組自体が空なら両方とも空の葉列を一つ返す。一つの Site はその葉を一つずつ返す。独自 Fn の出力 Item は制限せず、順序や対応を変える場合の契約は利用側が定める。

Site 順序と別途存在する全空間の因子順序を区別する。接続は空間・状態・行列を暗黙に並べ替えず、tensor 合成もしない。全作用空間と X/Y/Z などの選択・規格化は利用側、配置・積・bind は既存 core が担当する。各 Site から葉を一つ選ぶ操作は、複数葉に対する全 Pauli string の列挙とは区別する。

## 4. SpinHalf のモデル構築

モデル構築の公開入口は `qmb_spin_half::models::SpinHalfModels<SpinHalf<C>>` の関連関数である。以下では `type Models = SpinHalfModels<Spin>` とする。旧 `utility` の自由関数は残さず、エラー型 `UtilityError` を crate root から公開する。

| 呼び出し | 成功時の戻り値 |
| --- | --- |
| `Models::xyz(&chain, &factors, jx, jy, jz)` | `OperatorOf<Spin, I>` |
| `Models::heisenberg(&chain, &factors, coupling)` | `OperatorOf<Spin, I>` |
| `Models::neel(&chain, &factors, first)` | `StateOf<Spin, I>` |

全関数は Result を返し、列長と全葉の Index 重複を、辺がなくても検証する。列長の不一致は `UtilityError::Mapping(MappingError::SiteCountMismatch { expected, actual })`、重複は `UtilityError::Space(SpaceError::DuplicateIndex)`、配置・束縛エラーは `UtilityError::Expression` とし、元エラーを source として保持する。Index には Eq のみを要求する。

XYZ は `−Σ_(i,j) Σ_(a∈Li,b∈Lj) (Jx Sx_a Sx_b + Jy Sy_a Sy_b + Jz Sz_a Sz_b)` を構築する。Li は Site i の葉集合であり、`NearestNeighbors` と `leaves_cartesian()` の全組み合わせを使う。Site 内の葉間結合は追加しない。Heisenberg は三結合を J に揃え、正の J を強磁性的とする。標準モデルは位置対応を暗黙に選ばない。

XYZ・Heisenberg の全作用空間は、入力の局所空間ハンドルを共有し、列順に tensor した空間である。局所 grouping と相互作用に現れない Site も保持し、返す演算子をその空間に bind する。別途存在する全空間の root・上位 grouping・異なる因子順序は保存しない。

Neel は Site 位置と first から向きを交互に決め、その Site の全葉に同じ向きを置く。空 Site も偶奇に含める。状態の因子順は Site 順に局所空間の葉をたどった順とし、奇数長 periodic も認める。例えば `[(a,b), c]`、first=Up は `Up_a ⊗ Up_b ⊗ Down_c` となる。

0 Site は単位空間に束縛した zero Hamiltonian と単位状態を返す。全 Site が空でも同様であり、Hamiltonian の空間には空 Site の grouping を残す。

数値能力は個別関数へ限定する。xyz は `C: Clone + From<R> + Neg<Output=C>`、heisenberg はさらに `R: Clone` を要求する。neel の C に Clone・数値演算は不要である。モデル固有関数は QuantumSystem trait に追加しない。

## 5. 検証対象

geometry のテストは 0・1・2・3 以上の open／periodic、列長過不足、非 Clone・非 Eq 要素、借用の同一性、入れ子・空・可変長局所空間の順序、全組み合わせと位置対応、トップレベル分割の解釈と `new`／`from_space` の一致、3体・4体の連続配置、独自エラー・役割情報と遅延範囲検証を検証する。

spin-half のテストは複数葉の全組み合わせ結合、Site ごとの向き、空 Site の偶奇と辺、全空・0 Site、重複と列長拒否、非 Clone Index、深い grouping、最小の係数能力を検証する。examples の縦断テストは両精度の独立な Pauli 行列との一致、符号・norm・期待値、および密疎 export を検証する。

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

`n−I/2=−Z` なので、symmetric な局所相互作用は `U Z↑Z↓` とも表せる。両表現の値の一致を検証する。

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
- 両モデルの係数 bound は `C: Clone + From<R> + Neg<Output=C>`、入力 t・U は同じ R とする。symmetric 相互作用は `U Z↑Z↓` を使うため、構築に定数生成・除算・虚数単位や局所 Product・Act・Float を要求しない。t・U が実数に対応することはモデルの数学的契約とする。

局所演算子・Majorana の規格化・基底・spinful の位相規約は [Fermion の局所代数](local-structure.md#9-fermion-の局所代数)に従う。係数の恒等元・分数生成を避けて `U Z↑Z↓` を構築する実装は [models.rs](../../crates/qmb-fermion/src/models.rs) に置く。

### 6.3. 検証

- 二 mode の CAR、spinful の二重占有からの消滅符号、および全空間の mode 順を確認する。
- 1 Site の表示順 `[empty,up,down,double]` で、標準形の対角成分が `[0,0,0,U]`、symmetric 形が `[U/4,−U/4,−U/4,U/4]` となることを確認する。実際の列挙順との対応は明示する。
- 2 Site（4 mode、16状態）の両モデルを、core の評価・export を使わない Fock 空間の符号規則から構築した参照行列と比較する。実パラメータでの Hermite 性、粒子数・各 spin の粒子数の保存、`t=0`・`U=0` の極限を確認する。
- 両モデルの差が `−(U/2)N+(UL/4)I` であることと、symmetric 相互作用の `U Z↑Z↓` による構築を、定数項を含めて比較する。
- 0・1・2 Site と3 Site 以上の open / periodic、Site 数不一致、mode 数の過不足、ordinary な Site、Index 重複、非 Clone Index を検証する。
- `Complex<f32>` / `Complex<f64>` と Vec・ndarray・sprs の同じ基底順で値を比較し、構築の係数 bound に評価能力が混入していないことも確認する。

[モデルのテスト](../../crates/qmb-fermion/tests/models.rs)で局所規則と独立した Fock 空間の参照行列を使う。各辺に一度だけホッピングを追加すること、同じ spin のみが接続すること、symmetric 相互作用の定数項を保存することを確認する。実装時の実行結果・コマンド・未実施事項は [受入確認の記録](local-structure.md#96-実装時の検証記録)にまとめる。

## 7. SpinOne のモデル構築

`qmb_spin_one::models::SpinOneModels<SpinOne<C>>` は XYZ・Heisenberg・BLBQ の3モデルを関連関数として提供する。局所 primitive・基底・数値規約は [Spin-one の局所代数](local-structure.md#10-spin-one-の局所代数)に従い、構築は [models.rs](../../crates/qmb-spin-one/src/models.rs) に置く。

### 7.1. Hamiltonian の規約

`hbar=1` とし、まず1 Site が1 spin の場合を示す。和は Chain の最近接結合を各1回数える。

```math
H_{\mathrm{XYZ}}=-\sum_{\langle ij\rangle}
\left(J_xS_i^xS_j^x+J_yS_i^yS_j^y+J_zS_i^zS_j^z\right),
```

```math
H_{\mathrm{Heisenberg}}=-J\sum_{\langle ij\rangle}\mathbf S_i\cdot\mathbf S_j,
```

```math
H_{\mathrm{BLBQ}}=-\sum_{\langle ij\rangle}
\left[J\,\mathbf S_i\cdot\mathbf S_j
+K\,(\mathbf S_i\cdot\mathbf S_j)^2\right].
```

符号は既存 SpinHalf のモデルにそろえる。J>0 の bilinear 結合は強磁性的、J<0 は反強磁性的となる。`heisenberg(j)` は `xyz(j,j,j)`、`blbq(j,0)` は `heisenberg(j)` と数学的に同じである。単イオン異方性・磁場項・追加の定数シフトは設けず、上の式のエネルギー原点を保持する。

BLBQ は bilinear 結合 J と biquadratic 結合 K を独立に受け取る。角度表示 `H=sum[cos(theta) A_ij+sin(theta) A_ij²]`、`A_ij=S_i·S_j` は `blbq(-cos(theta),-sin(theta))` で表せる。角度専用の入口は設けない。Hermite Hamiltonian に対して実数値の結合を渡す契約とし、型による実数制限や実行時の判定は追加しない。

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

XYZ は配置した Cartesian 成分の積と和、Heisenberg は XYZ の等方的な場合として構築する。BLBQ は各葉ペアの `A_ab=Sx_a Sx_b+Sy_a Sy_b+Sz_a Sz_b` を一度構築し、式を共有して `(&A_ab * &A_ab)` を作る。J・K は最後に各項へ適用し、J の符号や大きさを二乗へ混入させない。

| 関数 | 構築に必要な bound |
| --- | --- |
| `xyz<I: Eq,R>` | `C: Clone + From<R> + Neg<Output=C>` |
| `heisenberg<I: Eq,R: Clone>` | `C: Clone + From<R> + Neg<Output=C>` |
| `blbq<I: Eq,R>` | `C: Clone + From<R> + Neg<Output=C>` |

構築時には局所規則の評価・基底展開・係数同士の加算や乗算を行わず、密行列・数値 backend・追加の局所規則は使わない。C の Float・Zero・One・除算・共役等も要求しない。評価・trace・export 時には既存 core と局所 primitive の操作別の bound が適用される。モデル固有メソッドを QuantumSystem trait に追加しない。

### 7.4. 検証と実装時の記録

[models.rs](../../crates/qmb-spin-one/tests/models.rs) では、直接記述した3×3スピン行列から参照行列を構成し、モデル実装や Chain の辺列挙を参照計算に使わない。

- `Complex<f32>` / `Complex<f64>` の両方で0〜3 Site・両境界の全行列成分・Hermite 性・trace を比較する。異方的結合、正負の J・K、zero 結合、pure biquadratic 項、等方的・bilinear 極限を含める。
- 二 spin の全スピン0・1・2に対する `A_ij` の固有値 -2・-1・1と、BLBQ のエネルギー `2J-4K`・`J-K`・`-J-K` を、独立した状態の期待値と作用の残差で検証する。
- 複数葉 Site の全葉ペアとペアごとの二乗、空 Site、0 / 1 Site、2 Site の周期辺の重複防止、非連番の因子順、grouping の共有と spectator の保持を確認する。
- 辺なし・zero 結合でも列長不一致と重複 Index を拒否すること、エラーの source、Clone・From・Neg だけの係数と非 Clone の結合入力・Index で構築できることを確認する。

2026-10-30、Rust 1.99.1 の direnv devShell で `cargo test -p qmb-spin-one --test models` の9件が成功した。これは実装時の記録であり、文書整理時の再実行結果ではない。workspace 全体を含む確認コマンドと未実施事項は [Spin-one の検証記録](local-structure.md#106-実装時の検証記録)を参照する。

利用例 [spin_one.rs](../../examples/src/bin/spin_one.rs) は `--no-default-features` と `--all-features` で実行した。`|+,->` の期待値は J=1 の等方的 XYZ / Heisenberg で1、J=1・K=0.25 の BLBQ で0.5となり、export 有効時には9×9の Vec / sprs 行列が一致した。
