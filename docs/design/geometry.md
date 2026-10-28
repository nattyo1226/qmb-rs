# Chain・Site と自由度の配置

更新日: 2026-10-28

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
