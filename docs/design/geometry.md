# Chain・Site と自由度の配置

更新日: 2026-10-30

## 1. 対象と責任分担

`qmb-geometry` の `Chain` はサイト数と境界条件だけを保持する純粋な幾何である。Index・Factor・HilbertSpace・局所代数・係数型を保持せず、Index を発行しない。Site は `0..site_count` の位置番号であり、専用の SiteId や所属トークンは設けない。

```rust
struct Chain {
    site_count: usize,
    boundary: Boundary,
}
```

`Chain::new(site_count, boundary)` は検証エラーなしで Self を返し、0 Site も認める。`site_count()` と `boundary()` で保持内容を取得する。構築後に変更しない。

2026-10-28 に Index 所有型の `Chain<I>`、`from_sites`・`with_sites`・Index 発行付き `new`、`ToSpace`、Index ペア用 `NearestNeighborRelation`・`index_pairs` を撤去し、以下の借用 API へ移行した。空間の構築・保持と Index 発行は利用側が担当する。core は geometry に依存せず、geometry は接続補助のため core に依存する。

## 2. 関係の列挙

`nn_site_pairs()` は正規化された番号対 `(i,j)`（i < j）を辞書順で遅延列挙する。open は隣接位置、periodic はさらに先頭と末尾を接続する。self-loop・重複辺は作らず、0・1 Site は辺なし、2 Site はどちらの境界条件でも一辺だけとする。

| メソッド | 入力と返すペア |
| --- | --- |
| `nn_site_pairs_with(&values)` | Site 順の任意の要素列から、要素をそのまま借用 |
| `nn_pairs_with(&factors)` | Site 順の局所空間列から、最近接 Site 間の葉の全組み合わせ |
| `nn_pairs_by_pos_with(&factors)` | 同じ葉位置同士を短い方まで zip |
| `nn_site_pairs_on(&space)` | トップレベル因子を Site として、その局所空間を借用 |
| `nn_pairs_on(&space)` | トップレベル Site 間の葉の全組み合わせ |
| `nn_pairs_by_pos_on(&space)` | トップレベル Site 間の葉の位置対応 |

`_with` は入力列の長さと site_count の一致を Iterator を返す前に検証し、過不足を `ChainError::SiteCountMismatch` とする。0・1 Site でも省略しない。葉ペアの入口は `nn_site_pairs_with` へ委譲し、`leaves()` で実際の葉を借用する。全組み合わせは左の葉を外側・右を内側とし、Site 対順を保つ。空 Site は位置を保ち、葉ペアを生成せず前後を直接接続し直さない。

`_on` は保存されたトップレベルの分割を読み取り、収集や Index の複製なしで `_with` へ委譲する。葉入力は一 Site、空 Tensor は0 Site、`Tensor([unit])` は一つの空 Site、単項 Tensor は一 Site と読む。局所空間内の grouping は保持する。全葉数から Site 数を決めない。

任意の要素 T に Clone・Eq・Default などを要求せず、葉ペアの列挙も Index の Clone・Eq を必要としない。Chain は入力列を保持しない。Index 重複は空間・モデル構築側の検証であり、幾何の列挙は重複を検証しない。

## 3. 空間と Index

局所空間列 `&[HilbertSpace<I>]` の一要素を一 Site に対応させる。葉・複数葉・入れ子・単項 Tensor・空の空間を使える。利用側は既存の空間ハンドルを共有して列を構成でき、Index の Clone・再発行は不要である。

Site 順序と別途存在する全空間の因子順序を区別する。接続補助は空間・状態・行列を暗黙に並べ替えない。任意の因子集合から部分系を抽出する API や、必須の配置ラッパーは導入しない。

## 4. SpinHalf のモデル構築

モデル構築の公開入口は `qmb_spin_half::SpinHalf<C>` の関連関数である。旧 `utility` の自由関数は残さず、エラー型 `UtilityError` を crate root から公開する。

| 呼び出し | 成功時の戻り値 |
| --- | --- |
| `Spin::xyz(&chain, &factors, jx, jy, jz)` | `OperatorOf<Spin, I>` |
| `Spin::heisenberg(&chain, &factors, coupling)` | `OperatorOf<Spin, I>` |
| `Spin::neel(&chain, &factors, first)` | `StateOf<Spin, I>` |

全関数は Result を返し、列長と全葉の Index 重複を、辺がなくても検証する。列長の不一致は `UtilityError::Chain(ChainError::SiteCountMismatch)`、重複は `UtilityError::Space(SpaceError::DuplicateIndex)`、配置・束縛エラーは `UtilityError::Expression` とし、元エラーを source として保持する。Index には Eq のみを要求する。

XYZ は `−Σ_(i,j) Σ_(a∈Li,b∈Lj) (Jx Sx_a Sx_b + Jy Sy_a Sy_b + Jz Sz_a Sz_b)` を構築する。Li は Site i の葉集合であり、`nn_pairs_with` の全組み合わせを使う。Site 内の葉間結合は追加しない。Heisenberg は三結合を J に揃え、正の J を強磁性的とする。標準モデルは位置対応を暗黙に選ばない。

XYZ・Heisenberg の全作用空間は、入力の局所空間ハンドルを共有し、列順に tensor した空間である。局所 grouping と相互作用に現れない Site も保持し、返す演算子をその空間に bind する。別途存在する全空間の root・上位 grouping・異なる因子順序は保存しない。

Neel は Site 位置と first から向きを交互に決め、その Site の全葉に同じ向きを置く。空 Site も偶奇に含める。状態の因子順は Site 順に局所空間の葉をたどった順とし、奇数長 periodic も認める。例えば `[(a,b), c]`、first=Up は `Up_a ⊗ Up_b ⊗ Down_c` となる。

0 Site は単位空間に束縛した zero Hamiltonian と単位状態を返す。全 Site が空でも同様であり、Hamiltonian の空間には空 Site の grouping を残す。

数値能力は個別関数へ限定する。xyz は `C: Clone + From<R> + Neg<Output=C>`、heisenberg はさらに `R: Clone` を要求する。neel の C に Clone・数値演算は不要である。モデル固有関数は QuantumSystem trait に追加しない。

## 5. 検証対象

geometry のテストは 0・1・2・3 以上の open／periodic、列長過不足、非 Clone・非 Eq 要素、借用の同一性、入れ子・空・可変長局所空間の順序、全組み合わせと位置対応、トップレベル分割の解釈と `_with`／`_on` の一致を検証する。

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

既存の `SpinHalf<C>` と同様に `models.rs` で構築し、`Fermion<C>` の関連関数として公開する。QuantumSystem trait にモデル固有メソッドを追加しない。

```rust,ignore
Fermion::hubbard(&chain, &sites, t, u)
Fermion::symmetric_hubbard(&chain, &sites, t, u)
// いずれも Result<OperatorOf<Fermion<C>, I>, UtilityError> を返す。
```

- `sites: &[HilbertSpace<I>]` の一要素を一つの物理 Site とし、各 Site は葉順 `(up,down)` の二 mode の graded 空間とする。Index は利用側で発行し、関数内で再発行しない。
- 全作用空間は Site 順、各 Site 内の up・down 順で graded 合成し、入力の Site 空間を共有して構築する。返す演算子はこの全空間に bind する。
- ホッピングは隣接 Site の同じ spin 同士だけを結ぶ。spin-half の XYZ が使う全葉ペアの直積は用いない。各辺に両方向のホッピングを一度ずつ加え、符号は core の graded 合成に委ねる。
- open / periodic の辺は既存 Chain の契約に従う。0・1 Site は辺なし、2 Site の periodic でも辺は一つとし、ホッピングを重複させない。
- Site 数の一致、各 Site が二つの葉を持つ graded 空間であること、全 mode の Index 重複を、結合値が0でも検証する。不正な Site 構成は `UtilityError::InvalidSite { site }` とし、列長・空間・式のエラーはそれぞれ `Chain`・`Space`・`Expression` variant に元エラーを保持する。
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

`qmb_spin_one::SpinOne<C>` は XYZ・Heisenberg・BLBQ の3モデルを関連関数として提供する。局所 primitive・基底・数値規約は [Spin-one の局所代数](local-structure.md#10-spin-one-の局所代数)に従い、構築は [models.rs](../../crates/qmb-spin-one/src/models.rs) に置く。

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
| `SpinOne::<C>::xyz(&chain, &sites, jx, jy, jz)` | `OperatorOf<SpinOne<C>,I>` |
| `SpinOne::<C>::heisenberg(&chain, &sites, j)` | `OperatorOf<SpinOne<C>,I>` |
| `SpinOne::<C>::blbq(&chain, &sites, j, k)` | `OperatorOf<SpinOne<C>,I>` |

全関数は `&Chain` と `&[HilbertSpace<I>]` を受け取り、`Result<_,qmb_spin_one::UtilityError>` を返す。Site 列・全葉ペア・束縛の基本契約は第4節の SpinHalf の XYZ / Heisenberg とそろえる。

- 一つの入力空間を一 Site とし、隣接 Site 間の全葉ペア `(a,b)` に相互作用を加える。BLBQ の二乗は各葉ペアの `A_ab=S_a·S_b` に個別に適用し、Site 内の結合や異なる葉ペア間の交差項を加えない。
- 列長と全 Site にわたる Index の一意性を、辺の有無・結合定数によらず検証する。Index には Eq のみを要求し、Clone・Hash・Ord は要求しない。
- 入力 Site の空間ハンドルを列順に tensor した全空間へ Hamiltonian を bind し、内部 grouping・葉順・相互作用に現れない Site を保持する。既存の別の全空間の上位 grouping や因子順は復元しない。
- 空 Site は位置を保ち、前後を直接つなぎ直さない。0 Site は単位空間上の zero Hamiltonian、1 Site と全空 Site は入力の空間構造に束縛した zero Hamiltonian を返す。
- 境界条件と辺は第2節の Chain に従う。0・1 Site は辺なし、2 Site は open / periodic とも一辺とする。

エラーは `Chain(ChainError)`・`Space(SpaceError)`・`Expression(ExpressionError)` の3 variant とし、元エラーを source として保持する。列長不一致は SiteCountMismatch、重複 Index は DuplicateIndex とし、UtilityError を crate root から公開する。

### 7.3. 遅延構築と係数能力

XYZ は配置した Cartesian 成分の積と和、Heisenberg は XYZ の等方的な場合として構築する。BLBQ は各葉ペアの `A_ab=Sx_a Sx_b+Sy_a Sy_b+Sz_a Sz_b` を一度構築し、式を共有して `A_ab.product(&A_ab)` を作る。J・K は最後に各項へ適用し、J の符号や大きさを二乗へ混入させない。

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
