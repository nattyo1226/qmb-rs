# HilbertSpace と Geometry の接続・保持構造の見直し

作成日・更新日: 2026-10-28
状態: 2026-10-28 実装済み。以下は採用した設計の記録であり、実装予定と表記したシグネチャも製品 API に反映した。正式な契約は[Chain の設計](../design/geometry.md)と[Hilbert 空間の設計](../design/hilbert-space.md)を参照する。Index 所有 API・ToSpace・Index ペア trait は外部列を借用する Chain のメソッドへ移行した。find_factor は任意候補として未追加。

## 1. 採用する方針と API

**Chain は純粋な幾何としてサイト数と境界条件だけを保持する。Site の内容・Factor・Index・HilbertSpace は保持しない。** Factor 列を受け取る近接ペアの補助を Chain に設け、幾何と自由度の配置を借用で接続する。

### 1.1. 保持と構築

```rust
// 実装予定の API。Chain 自体に型引数はない。
struct Chain {
    site_count: usize,
    boundary: Boundary,
}

impl Chain {
    pub fn new(site_count: usize, boundary: Boundary) -> Self;
    pub fn site_count(&self) -> usize;
    pub fn boundary(&self) -> Boundary;
    pub fn nn_site_pairs(&self) -> impl Iterator<Item = (usize, usize)> + '_;
}
```

0 Site を認め、構築では Index を発行しない。サイト数と境界条件だけの構築に検証エラーを必要としないため、`new` は Self を返す。近接 Site ペアの辞書順、open／periodic、小規模系の規約は現行を維持する。

### 1.2. Site 要素のペアと葉ペアを分ける

名前は次の規則に揃える。位置番号、入力要素、葉の三段階はメソッドの基本名で区別し、入力の形を末尾で示す。

- `_with`：Site 順の要素列・局所空間列を明示的に渡す。
- `_on`：全体の HilbertSpace を渡し、トップレベルの因子を Site と読む。
- 引数なしの `nn_site_pairs()`：Site の位置番号だけを返す。

`by_pos` は組み合わせ方を表す基本名の一部とし、その後ろに `_with`／`_on` を付ける。

| メソッド | 戻り値の各ペア | 葉への展開 |
| --- | --- | --- |
| `nn_site_pairs()` | Site の位置番号 `(i, j)` | 行わない |
| `nn_site_pairs_with(&values)` | `(&values[i], &values[j])` | 行わない |
| `nn_pairs_with(&factors)` | 最近接 Site 間の葉の全組み合わせ | 両方の局所空間の葉を直積で組み合わせる |
| `nn_pairs_by_pos_with(&factors)` | 最近接 Site 間の同じ葉位置同士 | 両方の葉列を短い方まで zip する |

`nn_site_pairs_with` は、入力の一要素を一 Site に対応させ、その要素をそのまま返す。HilbertSpace に依存しない汎用の補助とする。

```rust
// 実装予定の API。Iterator の具体型や内部ヘルパーの配置は実装時に決める。
impl Chain {
    pub fn nn_site_pairs_with<'a, T>(
        &'a self,
        values: &'a [T],
    ) -> Result<impl Iterator<Item = (&'a T, &'a T)> + 'a, ChainError>;

    pub fn nn_pairs_with<'a, I>(
        &'a self,
        factors: &'a [HilbertSpace<I>],
    ) -> Result<
        impl Iterator<Item = (&'a HilbertSpace<I>, &'a HilbertSpace<I>)> + 'a,
        ChainError,
    >;

    pub fn nn_pairs_by_pos_with<'a, I>(
        &'a self,
        factors: &'a [HilbertSpace<I>],
    ) -> Result<
        impl Iterator<Item = (&'a HilbertSpace<I>, &'a HilbertSpace<I>)> + 'a,
        ChainError,
    >;
}
```

処理の委譲関係は次のとおりとする。

```text
nn_site_pairs()
  └─ nn_site_pairs_with(&values): 位置対を入力要素への参照対へ写す
       ├─ nn_pairs_with(&factors): 双方の leaves() の全組み合わせ
       └─ nn_pairs_by_pos_with(&factors): 双方の leaves() の位置対応
```

上の二つの葉ペア用メソッドは、まず `nn_site_pairs_with(&factors)` を呼ぶ。列長の検証と最近接 Site 関係への接続を共通化し、その後の組み合わせ方だけを変える。葉を得て二空間を組み合わせる処理も、必要に応じて共通の補助へ分けられる。

| 項目 | 契約 |
| --- | --- |
| 列の長さ | `values.len() == site_count` を呼び出し時に検証し、不一致は `ChainError::SiteCountMismatch` 相当。余剰も不足も拒否する |
| 検証の時点 | Iterator を返す前。0・1 Site などペアが空の場合も省略しない |
| 列挙方法 | 入力への借用ペアを遅延列挙。要素・ペア全体の複製や収集を要求しない |
| Site 対の順序 | 三つの補助とも `nn_site_pairs()` の辞書順・境界条件を引き継ぐ |
| 葉ペアの順序 | Site 対ごとに、全組み合わせは左の葉を外側・右の葉を内側とする。位置対応は各空間の葉順を使う |
| 空の局所空間 | どちらかの葉列が空なら、その Site 対から葉ペアを生成しない |
| 型制約 | 汎用の Site 要素ペアは T に Clone・Eq・Hash・Default を要求しない。葉の走査・ペア列挙も I の Clone・Eq を要求しない |
| 保持期間 | 列を Chain に格納しない。Iterator が存在する間だけ Chain と列を借用する |
| Index 検証 | 各局所空間の葉にわたる重複 Index はモデル構築・空間側で検証する。標準モデルも葉数では入力を拒否しない |

**`at()` への接続には Index ではなく葉そのものを使う。一方、入力の `factors[i]` は葉に限定しない。** 各要素は Site i の局所空間を表す `HilbertSpace<I>` とし、葉・複合 Tensor・空 Tensor を扱える。

`leaves()` は、葉自身ならその一要素、Tensor なら子孫の葉を左から順、空 Tensor ならゼロ要素を返す借用 Iterator とする。葉ペア用メソッドはその実際の葉への参照を返すため、Index から配置先を検索し直したり、Factor を再構築したりする手順は不要である。

例えば二つの最近接 Site の局所空間の葉が `[a,b]` と `[c,d,e]` なら、結果は次のようになる。

```text
nn_site_pairs_with:   局所空間そのものの参照ペア
nn_pairs_with:        (a,c), (a,d), (a,e), (b,c), (b,d), (b,e)
nn_pairs_by_pos_with: (a,c), (b,d)
```

```rust
// 提案 API。factors の各要素は Site の局所空間。
let chain = Chain::new(factors.len(), Boundary::Open);
let sx = Spin::operator(SpinHalfOperator::Sx);

for (left, right) in chain.nn_pairs_with(&factors)? {
    // left と right は葉なので、そのまま at() に渡せる。
    let term = sx.at(left)?.product(&sx.at(right)?)?;
}
```

periodic な 3 Site では位置対の順は `(0,1), (0,2), (1,2)`、2 Site では open／periodic とも `(0,1)` 一つである。各位置対に対して、上記の葉ペアを順に展開する。

Chain が HilbertSpace を保持しないことと、葉ペア用メソッドが HilbertSpace の型・走査 API を利用することは区別する。この構成でも Chain の保持はサイト数・境界条件だけであり、core は geometry に依存しない。

#### 1.2.1. 全空間のトップレベルを Site と読む入口

局所空間列を渡す基本 API に加え、全体の `&HilbertSpace<I>` を受け取り、そのトップレベルの tensor 因子を Site の局所空間として解釈する入口を設ける。名前と委譲先は次のとおりとする。

| メソッド | 委譲先 | 戻り値の各ペア |
| --- | --- | --- |
| `nn_site_pairs_on(&space)` | `nn_site_pairs_with(top_level)` | 直下の局所空間への参照 |
| `nn_pairs_on(&space)` | `nn_pairs_with(top_level)` | 最近接 Site 間の葉の全組み合わせ |
| `nn_pairs_by_pos_on(&space)` | `nn_pairs_by_pos_with(top_level)` | 最近接 Site 間の同じ葉位置同士 |

`top_level` は、全空間の直下の子をそのまま並べた借用列である。全葉に平坦化してから Site を決めることはしない。各子の内部 grouping を保持し、葉への展開は委譲先の役割とする。Chain は受け取った全空間を保持せず、元の空間の木も変更しない。

```rust
// 提案 API。space = Tensor([Tensor([a,b]), c, Tensor([d,e])])。
// Site 0 は (a,b)、Site 1 は c、Site 2 は (d,e)。
let chain = Chain::new(3, Boundary::Open);
for (left, right) in chain.nn_pairs_on(&space)? {
    // 順に (a,c), (b,c), (c,d), (c,e) の実際の葉を返す。
    let term = sx.at(left)?.product(&sx.at(right)?)?;
}
```

端点は次の規約とする。葉だけの入力は一つの tensor 因子として、その葉自身を一 Site に対応させる。

| 入力 | Site の局所空間列 |
| --- | --- |
| 葉 a | `[a]`：1 Site |
| 空 Tensor `unit` | `[]`：0 Site |
| `Tensor([unit])` | `[unit]`：1 個の空 Site |
| `Tensor([a])` | `[a]`：1 Site |
| `Tensor([Tensor([a])])` | `[Tensor([a])]`：1 Site、局所 grouping を保持 |

解釈した Site 数と `chain.site_count()` の一致は、委譲先で列挙前に検証する。空・単項 Tensor と葉が数学的に互換でも、この入口は保存されたトップレベルの分割を明示的に Site と読むため、Site 数や配置が異なり得る。境界条件は引き続き Chain が持つ明示情報であり、空間から推定しない。

現行の `children()` は直下の子を借用する slice Iterator を返すので、その `as_slice()` を利用できる。葉入力は `core::slice::from_ref(space)` により一要素の借用列として扱える。したがって、この入口のための Factor の Clone・収集用 Vec・新しい Index の発行は不要である。既存の借用列入力に委譲し、入力全空間内の実際の葉への参照を返す。

### 1.3. 空間・配置・既存 API の責任

Factor 列は Site 順の対応を与えるもので、全空間の因子順序を書き換える指定ではない。同じ Chain に異なる列を渡せる。全空間と Factor の共有は利用側が維持し、Geometry 自体の共有データにはしない。

Site ごとの複数自由度・空の配置は、外側の長さが site_count の局所空間列として呼び出し側に置ける。`nn_site_pairs_with` は各局所空間をそのまま返し、`nn_pairs_with`／`nn_pairs_by_pos_with` がその葉を展開する。従来の全組み合わせ／位置対応という規約は維持し、戻り値を Index ペアから葉の空間への参照ペアに変更する。

Chain 単独から Space は生成できない。`ToSpace for Chain`、Index 入力の `from_sites`、Index 発行を行う従来の `new(site_count, degrees_per_site, boundary)` は、この設計では維持しない。空間・配置の構築は利用側へ移す。必須の配置ラッパーや集約型を新たに導入しない。

モデル関数は `&chain` と Site 順の局所空間列 `&factors` を受け取り、`&space` は不要とする。一要素を一 Site に対応させ、局所空間の共有ハンドルを列順に tensor して全空間を構築する。基本自由度の順序は、Site 順に各局所空間の葉順を連結したものになる。

### 1.4. モデルへの入力と実装上の整理

現行ソースと照合した結果、接続そのものに技術的な障害は見つかっていない。ただし新 API 全体を実装・コンパイル検証した結果ではない。

**モデルは Site 順の局所空間列を明示入力とし、全空間を引数にしない。** 具体的な呼び出しは次のとおりである。

| 関数 | 呼び出し |
| --- | --- |
| XYZ | `Spin::xyz(&chain, &factors, jx, jy, jz)` |
| Heisenberg | `Spin::heisenberg(&chain, &factors, coupling)` |
| Neel | `Spin::neel(&chain, &factors, first)` |

モデルの引数は局所空間のハンドル列 `&[HilbertSpace<I>]` とする。各要素の内部 grouping は共有した局所空間の木として保持できる。列順は Site 順であり、生成空間の葉順は各要素の葉順を連結して決まる。別途存在する全空間の root・上位 grouping や異なる因子順序を復元する契約は持たない。core の grouping・互換性・因子順序の契約は変更せず、既存の状態や行列を並べ替える機能も追加しない。

XYZ・Heisenberg は局所空間を共有して Site 順の全空間を構築し、その空間へ束縛した演算子を返す。これにより単独の trace 等に必要な作用空間と、相互作用がなく式に現れない Site の情報も保持できる。Neel の向きは Site 順で決め、対応する葉上の状態から積状態を構築する。0 Site では単位空間に束縛した zero Hamiltonian と単位状態を返す。

一方、次の点は現在の合意と API のまま進められる。

- **`at()` は変更不要**：現行の LocalState／LocalOperator は `&HilbertSpace<I>` を受け取り、葉でなければ `ExpressionError::NotAFactor` を返す。配置時にハンドルを共有するので、入力スライスを式の生存期間全体にわたって借用し続ける必要はない。新たな Factor 専用型や `at()` の対応範囲拡大は今回の必須変更にしない。
- **葉の列挙は借用 Iterator にする**：`leaves()` 相当は入れ子をたどって実際の葉の `&HilbertSpace<I>` を返す。直下の子とは異なり、既存の木のまま全葉の連続スライスを返すことはできない。各 Site の局所空間から必要時に葉をたどればよく、入力列の要素を葉に限定する理由にはならない。
- **型による葉保証は導入しない**：新たな Factor 専用型を必須にせず、`leaves()` が葉だけを返し、`at()` の既存の実行時検証も維持する。
- **検証をペア列挙に依存させない**：Chain は列長を検証し、モデル側は全局所空間にわたる葉の重複 Index を検証する。指定 Space との照合や、一 Site 一葉の検証は不要になる。0・1 Site や相互作用項のない場合も検証を省略しない。
- **標準モデルでも複数葉を受け入れる**：局所空間の葉数を理由に入力を拒否しない。XYZ・Heisenberg は最近接 Site 間の葉の全組み合わせへ結合を与え、Neel は Site ごとの向きをその Site の全葉へ割り当てる定義とする。詳細は第3.2節に記録する。

移行作業としては、ToSpace の Chain 実装の撤去、既存の Index ベースの関係 trait の見直し、Index 発行の利用側への移動、例・要件・回帰テストの更新が必要になる。これらは今回の責任分担から生じる作業範囲であり、それ自体を新たな幾何の意味の未決事項にはしない。

## 2. 設計時に確認した移行前実装と要件

### 2.1. 実装から確認できること

- [HilbertSpace](../../crates/qmb-core/src/space.rs) は `Arc<Node<I>>` を共有する不変な木であり、葉の `Factor(I)` と、順序付きの子空間を持つ `Tensor` がある。公開型はどちらも `HilbertSpace<I>` である。
- `children()` は直下の子、`indices()` は全葉の Index を左から順に列挙する。当時は葉ハンドルを列挙する公開 API はなく、今回 `leaves()` を追加した。
- ordinary な空間の互換性は順序付きの Index 列で決まり、grouping、空 Tensor、単項 Tensor の違いを同一視する。構築した木そのものは保存する。
- [Chain](../../crates/qmb-geometry/src/lib.rs) は `Vec<Vec<I>>` と境界条件を所有する。Site 順と Site 内の Index 順を保持し、全空間は保持しない。
- `ToSpace`（撤去済み） は `I: Clone + Eq` で葉を作り直し、Site ごとの空間を Site 順に tensor する。既存の全空間へ戻す逆変換ではない。
- [spin utility](../../crates/qmb-spin-half/src/utility.rs) は Chain と明示 Space を受け取り、Index の Eq で指定空間の葉に接続する。Neel の向きは Site 順、状態の因子順は指定空間に従う。

したがって、今回問題となる `space → geometry` は既存の自動変換 API の不具合ではなく、新しい構築 API にどの契約を与えるかという論点である。

### 2.2. 維持する条件と、変更し得る条件

| 条件 | 扱い |
| --- | --- |
| Index は代数・係数型・次元から独立 | 維持する |
| 同一性は Index の Eq、ポインタ一致は必須でない | 共有は所有・最適化の手段とする |
| 因子順序を暗黙に変更しない | 維持する。順序を変える空間の新規構築を既存状態の permutation と扱わない |
| Site 順と全空間の因子順は異なり得る | core の空間と幾何は独立のまま維持。新しい標準モデルの生成空間は入力 Factor 列順とし、現行 utility の別順序 Space 入力とは契約が変わるため、移行時に要件へ反映する |
| grouping の保存と数学的互換性の区別 | 維持する。各 Site の局所空間内の grouping を保持する |
| 0 Site、空 Site、可変長・複数自由度の Site | 維持する。空間の次元だけから Site 数を決めない |
| 重複 Index 拒否、明示生成時だけの Index 発行 | 系全体として維持する。合意した Chain は Index を扱わず、検証・発行を空間・配置構築側へ移す |
| 最小の型制約、no_std + alloc、深い木、次元 overflow の遅延検出 | 維持する |
| core は geometry に依存せず、空間構築に Geometry は不要 | 維持する |
| Chain は Index 列だけを持つ | サイト数・境界条件だけの保持へ変更する |

参照: [要件第3節・第7節](../requirements.md)、[Hilbert 空間](../design/hilbert-space.md)、[Chain](../design/geometry.md)、[Index](../design/index.md)。製品実装への移行時に、現行仕様の文書を新しい責任分担へ更新する。

## 3. API・モデル構築への影響

採用方針の Chain は Space も配置も持たず、モデルの配置入力は Site 順の Factor 列とする。

### 3.1. 子空間と葉の名前を分ける

- `space.children()`：保持された直下の tensor 因子を列挙する既存 API。
- `space.leaves()`：全葉の Factor ハンドルを tensor 順に遅延列挙する実装済み API。入力自身が葉の場合も一要素として列挙する。
- `space.find_factor(&index)` 相当：Index の Eq で実際の葉を検索する提案。

全葉列挙・検索は geometry と独立に core に置ける。任意の因子集合から新しい部分系を抽出する API は今回の前提にしない。

### 3.2. モデルは全空間ではなく Factor 列を受け取る

[SpinHalf の関連関数化の合意](api-type-parameters.md#43-モデル構築の関連関数合意済みの方針)を、Chain と Site 順の Factor 列を受け取る形に更新する。全空間の引数は不要である。

列長の一致と、全局所空間にわたる葉の Index の重複がないことを検証する。入力要素が葉であることや、各要素の葉数が1であることは要求しない。局所空間を共有する tensor 構築には I の Clone や新しい Index の発行は不要で、Eq による重複検証を維持する。

Neel の向きは Site 順、状態の tensor 順は Site 順に局所空間の葉をたどった順とする。Hamiltonian は局所空間列から構築した全空間に bind する。局所空間内の grouping はそのハンドルを共有して保存できるが、別途保持する全空間の元 root・上位 grouping の再現は保証しない。

標準モデルも複数葉・空の局所空間を受け入れる。一 Site 一 spin の場合を含む定義として、次の規約を使う。

- **XYZ・Heisenberg**：`nn_pairs_with(&factors)` の返す全ての葉ペアへ従来の結合を与える。Site i の葉集合を Li とすると、XYZ は `−Σ_(i,j) Σ_(a∈Li,b∈Lj) (Jx Sx_a Sx_b + Jy Sy_a Sy_b + Jz Sz_a Sz_b)`、Heisenberg は Jx・Jy・Jz を J に揃える。Site 内の葉同士の相互作用は自動追加しない。
- **Neel**：向きは Site 位置 i と first から従来どおり交互に決め、Site i の全葉にその同じ向きを置く。葉ごとに向きを交互にする処理にはしない。
- **空 Site**：空の局所空間からは相互作用項も局所状態も生成しない。Site 自体は位置を保ち、前後の Site を直接接続し直さない。Neel の偶奇にもその Site の位置を含める。

この Heisenberg の定義は、Site の総スピン `T_i = Σ_(a∈Li) S_a` を用いると `H = −J Σ_(i,j) T_i · T_j` と書ける。複数葉を持つ Site の総スピン同士の結合として扱い、固定の総スピン部分空間への射影は行わない。

例えば局所空間列が `[(a,b), c]` で first が Up なら、Neel は葉順に `Up_a ⊗ Up_b ⊗ Down_c` となる。open Chain の Heisenberg 結合は葉ペア `(a,c)` と `(b,c)` に与える。各 Site が一葉なら従来の式に一致する。全ての Site が空なら単位状態と単位空間に束縛した zero Hamiltonian を返す。

葉の位置対応は `nn_pairs_by_pos_with` が提供する別の組み合わせ規約として維持し、標準の XYZ・Heisenberg が暗黙に位置対応を選ぶことはしない。0 Site、norm、結合の符号、周期境界の小規模規約は維持する。

## 4. 実装へ進む場合の範囲と確認事項

合意した方向は第1節の純粋な Chain と借用列への補助である。製品実装への移行時には次を行う。

1. Chain の保持をサイト数・境界条件だけに変更し、型引数と Index 入力・発行機能を除く。
2. 既存の `nn_site_pairs()` の規約を維持し、外部列を受け取る `nn_site_pairs_with(&values)` を位置対への写像として実装する。`nn_pairs_with`／`nn_pairs_by_pos_with` はこれを使い、それぞれ全組み合わせ／位置対応の葉ペアを返す。
   全空間を受け取る `_on` の入口はトップレベルの借用列を取り出し、各基本 API へ委譲する。
3. 長さの一致を列挙前に検証し、0・1・2・3 以上の Site、open／periodic、過不足、非 Clone・非 Eq の要素、借用の同一性を確認する。葉の全組み合わせと位置対応の順序、空・可変長・入れ子の局所空間も確認する。全空間入口ではトップレベルと全葉の区別、葉入力の一 Site、空 Tensor と一つの空 Site の区別、列入力経路との結果の一致を確認する。
4. 現行の Index ペア用 trait・共通関数・ToSpace を移行する。spin モデルは全空間の代わりに Site 順の局所空間列を受け取り、列順の空間を構築する。Site 対と入力列の長さの不一致は拒否する。葉の位置対応だけは従来どおり短い方まで列挙する。
5. 空間構築・Factor の共有・Index の重複検証を利用側へ移し、必要な core の読み取り補助は独立に検討する。
6. 要件・geometry／HilbertSpace の設計・architecture・examples を新しい責任分担へ揃え、既存のモデルの符号・norm・因子順序・数値結果を回帰検証する。標準モデルの複数葉・空 Site の受理、全組み合わせの相互作用、Site ごとの Neel の向き、一葉の場合との一致も確認する。

本書は採用する設計の記録である。製品実装と回帰テストは geometry・spin-half・examples に反映した。Iterator は借用の impl Iterator とし、_on はトップレベル slice を取り出して _with へ委譲する。
