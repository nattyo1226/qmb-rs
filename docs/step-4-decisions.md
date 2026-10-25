# Step 4 の実装方針と合意事項

更新日: 2026-10-25

## 1. 位置付けと実装範囲

本書は、[実装計画](implementation-plan.md)の Step 4 に適用する合意済みの設計と受入条件をまとめる。実装と検証の状況は[実装計画の段階4](implementation-plan.md#44-段階4の実装検証状況)に記録する。数学的契約は[要件](requirements.md)、提供範囲は[初期版の設計](design/initial-release.md)に従う。

| crate | Step 4 の担当 |
| --- | --- |
| `qmb-core` | 局所式・複合式の norm、正規化、正規化期待値、共通エラー、任意の `libm` feature |
| `qmb-basic` | UUID の具体型を `Index` から `Id128` へ改名し、関連する型別名・利用例・テストを移行 |
| `qmb-geometry` | Chain、空間の保持と構築、Site・Index の近接ペア列挙。`no_std + alloc` を維持 |
| `qmb-spin-half` | XYZ・Heisenberg Hamiltonian・Neel 状態の構築 utility |
| `qmb` | 各 crate の再公開、利用例、縦断的な結合テスト。既定で core の `libm` を有効化 |

Partial trace・複合基底・座標取得は Step 5、数値 backend・export は Step 6 とする。bra 型・内部ラベルは初期版では導入しない。

## 2. Geometry・空間の保持と構築

- Geometry の生成時に `L: LocalStructure` を指定し、Site 順の `Vec<HilbertSpace<I, L>>` と全 Hilbert 空間を保持する。既存空間からは、その型により `L` が定まる。
- HilbertSpace は内部の `Arc` で節点を共有するハンドルであり、追加の `Arc<HilbertSpace<…>>` は設けない。Site の空間と全空間で節点を共有する。
- Site は `Vec` の位置を表す `usize` とする。専用の SiteId・所属トークンは設けず、Geometry を跨ぐ同一性・混用を管理しない。範囲外は検出する。
- 配置情報は Geometry と一体で保持し、独立した Placement 型や対応表を要求しない。構築後の配置・順序は不変とする。
- Grouping は HilbertSpace の構造をそのまま使う。Site 内の自由度は構成要素・葉の順序と既存 Factor のハンドルで参照し、内部ラベルの型・対応表・検証処理は追加しない。

構築経路は次の二つとする。

1. **既存空間から構築**：Site ごとの Factor・Tensor を Site 順に受け取り、Index・因子順序・grouping を保持する。異なる Site に重複する Index は拒否する。既存の全空間を受け取る場合も、その因子順序を暗黙に変更しない。
2. **個数から構築**：Site 数と各 Site 共通の局所自由度数を受け取り、指定した Index 型と局所構造型で空間を構築する。自由度数は基本因子の個数を意味し、空間の次元とは区別する。

全空間を新しく構成する場合は Site 順・Site 内の順に tensor product を取る。Site 数ゼロでは空の Geometry と単位空間を作る。各 Site の自由度数ゼロでは Site 数を保ち、各 Site に空の Tensor を対応付ける。保持済みの全 Hilbert 空間を取得するメソッドを提供し、取得時に Index を再発行したり空間を再構築したりしない。

## 3. Index の発行と具体型

個数からの構築など、新規 Index の発行が必要な箇所でのみ `I: Default` を要求し、`I::default()` を呼ぶ。この用途に指定する型は、呼び出すたびに新しい独立した Index を発行する契約を満たす。一般の Default 自体が一意性を保証するわけではなく、型の提供側が発行の意味と一意性の前提を満たす。

照合は標準の `Eq` を用いる。既存空間の保持・参照・構築には Default・Clone を一律に要求しない。専用の Index 生成 trait や生成関数を渡す引数は追加せず、UUID・乱数源への依存は具体型の実装側に置く。

`qmb_basic::Index` は `qmb_basic::Id128` に改名する。UUID v4 を包む不透明な型、seed 指定不要の `new()`、Clone・Copy による同一性保持、表示・比較・hash の既存契約を維持する。`Id128::default()` は `Id128::new()` と同じく新しい UUID v4 を発行する。Geometry は `qmb-basic` に依存しない。詳細は [Index の契約](design/index.md)に従う。

## 4. 幾何関係と公開 API

Chain の境界条件は open / periodic とする。近接 Site ペアは `(min(i,j), max(i,j))` に正規化して辞書順で列挙し、self-loop・重複辺を作らない。長さ `0,1` では辺なし、長さ `2` では open / periodic とも `(0,1)` 一つとする。

Geometry 側には次を公開する。

- `sites(&self) -> &[HilbertSpace<I, L>]`
- `nn_site_pairs(&self) -> impl Iterator<Item = (usize, usize)> + '_`

返す Site 番号は `sites()` の範囲内とする。任意の番号の参照には `sites().get(site)` を利用できる。

幾何関係ごとに能力 trait を設ける。近接関係は以下とし、異なる Geometry 型にも同じ trait を実装できるようにする。

```rust
pub trait NearestNeighborRelation {
    type Index;

    fn nn_pairs(&self)
        -> impl Iterator<Item = (&Self::Index, &Self::Index)> + '_;

    fn nn_pairs_by_pos(&self)
        -> impl Iterator<Item = (&Self::Index, &Self::Index)> + '_;
}
```

`qmb-geometry` は二つの HilbertSpace を借用する共通関数を公開する。葉の Index 順序を用い、返す Index は元の空間から借用する。

| 共通関数 | 関係 trait のメソッド | 規則 | `[a,b]` と `[c,d,e]` の例 |
| --- | --- | --- | --- |
| `index_pairs(left, right)` | `nn_pairs()` | 全組合せ。左側を外側、右側を内側として列挙 | `(a,c),(a,d),(a,e),(b,c),(b,d),(b,e)` |
| `index_pairs_by_pos(left, right)` | `nn_pairs_by_pos()` | 同じ位置同士を先頭から短い方の長さまで列挙 | `(a,c),(b,d)` |

trait の実装は、`nn_site_pairs()` の各ペアに対応する空間を `sites()` から参照し、共通関数の結果を `flat_map` する。いずれの方式も Site ペアの列挙順を維持し、片側が空ならペアを返さない。自由度数の違い自体はエラーにしない。共通関数は Index の Clone・Default を要求せず、ペア全体を Vec に格納しない。追加の変換用 trait・Result は要求しない。

これらは幾何関係の列挙であり、Hamiltonian の項を自動生成するものではない。具体的な追加 Geometry・関係の実装は初期版の対象外とする。

## 5. XYZ・Heisenberg・Neel の構築 utility

構築ロジックは `qmb-spin-half` に置き、`qmb` は再公開と利用例を担当する。幾何関係は `qmb-geometry` の公開契約を利用する。

- 対象は各 Site がちょうど一つの spin Factor を持つ配置。空 Site・複数自由度の Site はこの utility ではエラーにする。Geometry 自体はこれらの Site を扱える。
- XYZ は実数の結合 `Jx, Jy, Jz: T` を受け取り、`−Σ(Jx Sx_i Sx_j + Jy Sy_i Sy_j + Jz Sz_i Sz_j)` を通常の演算子積と和で構築し、全空間に bind して返す。Heisenberg は `Jx = Jy = Jz = J` として XYZ に委譲し、同じ負符号の規約を用いる。`J > 0` は強磁性的結合とする。
- Neel は先頭の Up / Down を指定でき、同じ配置から交互の積状態を構築する。奇数長 periodic でも構築を許可し、端点間は同じ向きになる。
- Index は再発行しない。0-site では単位状態と、全空間へ束縛した zero Hamiltonian を構築する。
- 複数自由度から作用対象を選ぶ utility は後続で明示的な選択方法を設計する。

## 6. Norm・正規化・内積

局所式と複合式に以下を追加する。

| 操作 | 契約 | `libm` feature |
| --- | --- | --- |
| `norm_squared()` | 自己内積の実部を `T` で返す | 不要 |
| `norm()` | norm squared の平方根を `T` で返す | 必要 |
| `normalized()` | norm で除した係数の状態式を返す | 必要 |
| `normalized_expectation()` | 期待値を複素数の自己内積で除算する | 不要 |

複合式には `normalized_expectation_in()` も設け、既存の明示空間優先の選択に従う。通常の作用・内積を利用し、行列化・基底展開は要求しない。自己内積が実数になる数学的契約に基づき norm の実部を扱い、負値の補正や近似ゼロ判定は導入しない。

局所式・複合式で共通の `NormalizationError` を設け、`ZeroNorm` と既存の `ExpressionError` のラップを持たせる。正規化・正規化期待値は計算された自己内積の厳密なゼロを拒否する。明示的な zero に加え、相殺による zero も対象にする。`norm(zero)` は `0` を返す。NaN・infinity などは数値型の挙動に委ねる。

bra 型・ラッパー・構築メソッドは初期版では導入しない。braket は `inner_product` に、演算子を含む braket は作用と内積に委ねる。左引数の共役線型性・右引数の線型性は内積の契約として扱う。

## 7. Feature・依存関係・再公開

core は `no_std + alloc` を維持し、以下の任意 feature を追加する。

```toml
[features]
default = []
libm = ["num-traits/libm"]
```

`norm()`・`normalized()` の平方根には `num_traits::Float` を利用し、Float 制約はそのメソッドに限定する。通常の演算、norm squared、正規化期待値へ制約を広げない。`libm` の production 依存の追加は了承済みとする。

`qmb` は core の `libm` を既定で有効にし、無効化も可能にする。数値 export は Step 6 で導入する。各 crate に対応する module を設け、crate root の再公開は主要型に絞る。

Cargo の feature は加算される。確認した `num-traits 0.2.19` は std と libm が両方有効なら std 側の Float 実装を使うため、この feature は libm の実装を強制しない。libm 有効化自体は std を要求しない。

## 8. 受入条件と検証

- 両精度・両方の Neel 開始向きで、行列化せず open chain の期待値 `J(L−1)/4`（`L >= 1`）と norm `1` を得る。実数 J に対する Hermiticity を独立した小規模行列とも照合する。
- 0-site の状態は単位状態、norm は `1`、Hamiltonian と期待値は `0`。長さ `0,1,2` の open / periodic、辺の辞書順・正規化・self-loop と重複辺の不在を確認する。
- Factor・入れ子 Tensor・空 Site・複数自由度の配置、既存 Index・因子順序・grouping・共有の保持、範囲外・Index 重複・入力数不一致の拒否を確認する。内部ラベルと Geometry を跨ぐ所属チェックは対象外とする。
- 個数からの構築で Default の呼び出し回数・発行順序と Id128 の発行・複製契約を確認する。保持・参照・Clone で再発行せず、Default・Clone を持たない Index でも既存空間から構築できることを検証する。
- Chain と非公開の別 Geometry を同じ関係 trait で扱う。二つのペア方式について、可変長・空・入れ子の葉順・Site ペア順を検証する。
- 内積の共役線型性・線型性、作用を含む内積、非正規直交な例の norm、zero と相殺による zero-norm の拒否を確認する。既存の非公開の別代数による検証も拡張する。
- libm 有効・無効、qmb の既定構成での norm 利用、無効時の平方根を使わない操作、core・geometry の no_std と通常依存を確認する。
- 対象テストの後、workspace の test・format check・Clippy・警告をエラーとするドキュメント生成を行う。許容誤差は精度別に明示する。MSRV・no_std target の確認が環境上できない場合は未確認として報告する。

## 9. 後続の必須課題

State・Operator・Structure・Space の管理方法と Geometry との依存関係は将来必ず見直す。core の数値ライブラリ依存解消と四つの設計上の問いは、[アーキテクチャの将来の必須課題](architecture.md#10-将来の必須課題)を正として管理する。今回の Step 4 の完了条件には追加せず、現行の型構成・libm 方針で実装する。

## 10. 実装した構築 API と利用例

- `Geometry::<I, L>::new(site_count, degrees_per_site)`：`I: Default + Eq` で個数から生成する。1自由度の Site は Factor、0自由度と複数自由度の Site は Tensor とする。
- `Geometry::from_sites(sites)`：既存空間を Site 順に保持し、その順序で全空間を構成する。
- `Geometry::with_space(sites, &space)`：指定した全空間を共有し、各 Site が内部の葉順を保って全 Index を重複なく分割することを確認する。Site 順と全空間の因子順は一致しなくてもよい。
- `Chain::new(site_count, degrees_per_site, boundary)`、`Chain::from_sites(sites, boundary)`、`Chain::from_geometry(geometry, boundary)`：`Boundary::{Open, Periodic}` を指定して Chain を作る。`Chain::with_sites(site_count, sites, boundary)` は指定した Site 数との一致も確認する。
- `Geometry::space()`・`Chain::space()` は保持済み全空間を借用する。`Geometry::site(site)` は範囲外をエラーとし、`sites().get(site)` も利用できる。
- spin-half の構築補助は `qmb_spin_half::utility::{xyz, heisenberg, neel}`、エラーは `utility::UtilityError` とする。`xyz(&chain, jx, jy, jz)` に対して `heisenberg(&chain, j)` は `xyz(&chain, j, j, j)` に相当する。
- エラーの定義・表示・標準の変換は原則として `thiserror::Error` と属性で記述する。Geometry と utility の既存メッセージ・変換元の `source()` は保持する。
- 複合式の `norm_squared()`・`norm()` は既存の fallible な内積と同じく `Result<T, ExpressionError>`、局所式は `T` を返す。正規化と正規化期待値は共通の `NormalizationError` を返す。

再公開は `qmb::{core, basic, geometry, spin_half}` にまとめ、主要型を root にも公開した。[実行可能な例](../crates/qmb/examples/heisenberg.rs)は `cargo run -p qmb --example heisenberg` で実行できる。6-site open chain、`J = 1` の Neel 期待値は `1.25`、norm は `1` となる。`--no-default-features` でも、平方根を使わない正規化期待値を計算できる。
