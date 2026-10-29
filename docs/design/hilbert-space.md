# Hilbert 空間の再帰構造

更新日: 2026-10-29

## 1. 対象と合意した契約

本書は、[要件](../requirements.md)と[アーキテクチャ概要](../architecture.md)に基づき、具体的な自由度を持つ Hilbert 空間の構成、grouping、および互換性を定める。

ジェネリックな Index の照合・Hilbert 空間と、その上の共通の式・評価は `qmb-core` が担当する。具体的な UUID Index の発行は `qmb-basic`、Site への対応付けは `qmb-geometry` に置き、空間の構築・合成に Geometry を必須としない。現在の公開 API は第10節に記録する。

Hilbert 空間自体を再帰的な構造とし、基本自由度の空間と、順序付きの子空間の tensor product を表す。複合空間をさらに複合させられるようにし、grouping はその木の括り方で表す。提供する能力は構築・合成、構成要素・Index の列挙、互換性判定、および部分系判定とする。任意の因子集合を指定した部分系の取得と Regrouping の操作は提供しない。Grouping のための独立した Group 型は追加しない。Factor と Tensor は概念上の名前であり、具体的な Rust の型・公開 API を指定するものではない。

初期版の空間は局所代数の型を保持しない。式側の S・O・C の型制約が同種の式の合成を保証する。空間構造を共有しても、異なる primitive 型の式を暗黙に混在させる機構は導入しない。

## 2. Factor と Tensor

公開の空間型は `HilbertSpace` に統一し、自由度（Dof）は Index を持つその葉（Factor）として表す。独立した `Dof` 型や、Dof から単因子空間への変換は設けない。葉の空間をそのまま保持・共有し、tensor product の入力や、局所式の配置先を指定するために利用できる。配置対象は一つの基本因子であることを検証し、Tensor から局所式の配置先の Index を推測しない。具体的な constructor・参照 API と型引数は別途定める。

| 対象 | 情報・役割 |
| --- | --- |
| Factor | 基本因子の Index を持つ葉。局所次元は取得時の入力から得る |
| Tensor | ordinary / graded の積の種類と、順序付きの子 Hilbert 空間を持つ合成対象。子に Factor と Tensor を含められる |

```text
Hilbert 空間
├─ Factor(Index)
└─ Tensor(積の種類, 順序付きの子 Hilbert 空間列)
```

次元は `dim::<P>()` で指定した primitive の span を各葉に適用するか、`dim_with(dim_local)` で局所次元を与える。Tensor の次元は子空間の次元の積、空の Tensor は一次元とする。

総次元が `usize` に収まることを空間構築の条件にしない。次元取得 API は `Result<usize, DimensionOverflow>` 相当とし、checked arithmetic により overflow を検出する。互換性・部分系判定や行列を作らない作用・内積・trace では、総次元の整数化を必須にしない。

```math
\dim \operatorname{Tensor}(k;\mathcal H_1,\ldots,\mathcal H_n)
=\prod_{j=1}^{n}\dim\mathcal H_j,
\qquad
\dim\mathcal H_{\varnothing}=1.
```

基底は空間とは別に選ぶ。空間の木の構築や基本因子の列挙によって、全基底、状態・演算子の式、または行列を全展開することは要求しない。

## 3. Index と因子順序

自由度の Index は Factor に置き、Tensor に追加の自由度の Index を発行しない。Factor は、あらかじめ発行した Index によって指定する。[Index の契約](index.md)に従い、空間の構築・複製・既存因子の再利用で新しい Index を暗黙に発行しない。

基本因子の順序は、Tensor の子を左から順にたどった葉の順序とする。Factor が持つ Index の順序付き列として確認でき、Index の値の大小から決めない。空間の tensor product は左 operand の葉を左側、右 operand の葉を右側に置く。

重複は全階層の基本因子の Index について検証する。親となる Tensor とその内部因子を独立した因子として追加する場合や、共有した部分空間を二度追加する場合も、葉の Index が重複すれば拒否する。既存因子の記述を複製・保持することとは区別する。

## 4. 保持する木構造と空間の互換性

空間の木は、利用者が指定した括り方を保持できる。互換性判定や評価のために構築時から平坦化することは要求しない。必要な時点で基本因子をたどり、以下の同一視・結合則を用いて判定する。

- 空の Tensor は、積の種類によらず共通の一次元の単位空間と同一視する。graded な構成ではその状態 `1` を even とする。
- 単因子 Tensor は、その子空間と同一視する。
- 同じ種類の tensor product に対する結合則を、互換性判定の規則として適用する。括り方の違いだけなら同一視する。
- Index を持つ一次元の Factor は、基本因子の存在を持たない空の Tensor と区別する。次元だけを理由に Factor を取り除かない。

二つの空間構造は、基本因子の Index・順序と、空・単因子の同一視・結合則を適用した tensor 構造が一致するときに互換とする。primitive の型と規則は式側が持ち、空間比較は追加の定義ラベルを確認しない。木の括り方だけの違いは互換性を変えず、object identity や次元の一致だけでは判定しない。

初期版では、空・単因子の同一視を適用した後にも残る ordinary / graded の混在を一つの空間内で拒否する。互換性判定で積の種類を変更しない。利用者向けの Regrouping 操作は設けず、因子の並べ替えと permutation もサポートしない。

同一視は空間の数学的な互換性について定める。保持・表示する階層を必ず消すことは要求しない。数値精度・scalar・grading などの合成条件も、既存要件に従って別途検証する。

## 5. Grouping、幾何、および局所複合代数

Site 上の複数の自由度は、利用側が構成する Tensor の内部因子として表せる。Core は空間の木を扱い、Chain はサイト数と境界条件だけを保持する。Site 順の局所空間列・全空間・grouping は利用側が管理する。初期版は内部ラベルを導入せず、幾何を空間構築の必須情報にしない。

```text
大域空間: Tensor
├─ Tensor                         ← 配置層で site 0 に対応付ける
│  ├─ Factor(Index u0, mode)      ← Site 内の因子 0
│  └─ Factor(Index d0, mode)      ← Site 内の因子 1
└─ Tensor                         ← 配置層で site 1 に対応付ける
   ├─ Factor(Index u1, mode)      ← Site 内の因子 0
   └─ Factor(Index d1, mode)      ← Site 内の因子 1
```

この spinful fermion の例では、複数因子を合成する Tensor を graded とする。基本因子の順序は `u0, d0, u1, d1` となる。Site ordering と内部順序からこの並びを構築する処理は配置・モデル構築側が担当し、対応を明示する。Core は指定された子空間の順序を用いる。

局所複合代数の空間構成にも、この再帰構造を用いる。積・作用・内積・随伴・trace・grading などの数学的能力は、子の局所規則と共通の合成・評価機構から処理する。[局所 primitive・式木の数学的能力](primitive-expression.md)との責任分担を維持する。Grouping のために基本局所構造の primitive や規則を作り直すことは要求しない。

複合 Hilbert 空間は、局所構造に相当する整合性の責任を持ち、子の構造・因子順序・積の種類から数学的能力を合成する。共通機構は局所契約を保存する合成則を実現し、同種の積の括り方や空・単因子の同一視によって結果を変えない。独立した複合構造型は要求しない。能力一覧と合成・埋め込み・縮約の契約は、[埋め込み・複合空間の式と数学的能力](composite-expression.md)に定める。

## 6. 部分系の判定

ここでの部分系は、テンソル因子からなる空間を指す。一般の線形部分空間や、部分系の状態を大域状態へ埋め込む写像を定義するものではない。

構築条件を満たす二つの空間 `K,H` について、`K` が `H` の部分系であることを、以下のすべてを満たす場合に認める。

1. `K` のすべての基本因子の Index が `H` に存在する。
3. `K` の基本因子列が、`H` の基本因子列の部分列になっている。因子の相対順序を保存し、並べ替えによる対応付けは行わない。
4. 空・単因子の同一視後の tensor product の種類が整合する。複数の基本因子を含む `K` は、`H` と同じ積の種類を持つ必要がある。空・単因子については 4 節の同一視を適用する。

`K` が `H` の木の部分木として存在することや、基本因子が連続して並ぶことは要求しない。利用者が構築時に保持した基本因子から別途 `K` を構築し、判定対象として渡せるようにする。空間自身はその部分系とし、空の Tensor はすべての空間の部分系とする。Index を持つ一次元の Factor は、その Index が対応する場合にのみ部分系とする。

数値精度の互換性は、部分系の数学的な判定とは別に、演算の接続時に既存要件の規約に従って検証する。部分系判定は、partial trace の対象などに対する空間の整合性検証にも用いる。判定のために部分系の取得操作や Regrouping を行うことは要求しない。

## 7. 構成要素の列挙と不変構造の共有

利用者は、構築時に用意した葉（基本自由度）や Tensor の `HilbertSpace` を保持し、作用先や partial trace の対象として指定できる。さらに、配置などで既存の構成を読み取れるよう、次の三つの遅延列挙を提供する。

| 能力 | 列挙する対象 | 順序・契約 |
| --- | --- | --- |
| Tensor の構成要素の列挙 | 直下の子 Hilbert 空間 | 子の構築順序を保持する。子が Tensor ならその括りを保持する |
| 葉の列挙 | 入れ子内の実際の Factor ハンドル | `leaves()` で因子順に借用。葉入力は自身を一要素として返し、空 Tensor はゼロ要素 |
| 所属 Index の列挙 | 入れ子を再帰的にたどった全 Factor の Index | 葉を左からたどる因子順序。Factor は一要素、空の Tensor はゼロ要素 |

列挙は共有借用を基本とする Iterator とし、`children()` 相当で直下の子、`leaves()` で全葉のハンドル、`indices()` で全葉の Index を読み取る。Factor の直下の子と空の Tensor の子はゼロ要素とする。全要素の一括生成を要求しない。空・単因子の数学的な同一視にかかわらず、構成要素の列挙は保持している子構造を読む。元の空間・Index・局所構造・grouping を変更せず、式や基底や行列の展開も行わない。具体的な返却型・シグネチャは実装時に定める。

これは既存の構成の読み取りであり、任意の因子集合から新しい部分系を取得する操作とは区別する。利用側は Site 順の局所空間列を Chain に借用して渡す。Chain は列長を検証し、モデル・空間構築側が全葉の Index 重複を検証する。責任分担は [Chain・Site と自由度の配置](geometry.md)に定める。

大域空間へ演算子を埋め込むときは、指定された Index と、その空間の基本因子との整合性を検証する。部分系の取得元や木の位置を記録する特別な対象は設けない。埋め込み・作用・partial trace の契約は要件の 5 節に従い、core 内部での基本因子の走査は検証・評価のために行える。

`HilbertSpace` は `Arc` による不変な共有ハンドルとし、Factor・Tensor の節点 enum は非公開とする。複製で Index を変更せず、変更されない部分空間の不要な deep copy を避ける。走査と最後の所有者になった節点の破棄には明示的な作業スタックを使い、深さに比例した関数再帰を避ける。初期版の ordinary な空間の互換性は葉列の照合、部分系判定は順序を保つ部分列の照合で処理し、保持している木を変更しない。

## 8. 例と検証条件

以下で `a,b,c` は重複のない基本因子で、同じ局所構造を持ち、scalar などの合成条件は互換とする。`T_o` と `T_g` はそれぞれ ordinary と graded の Tensor、`E` は空の Tensor とする。Graded の例も実装済みの契約であり、異種局所構造の混在を認める例ではない。

| 対象・操作 | 契約上の結果 |
| --- | --- |
| `T_o(T_o(a,b),c)` と `T_o(a,T_o(b,c))` の比較 | 同じ基本因子・順序・積の種類を持つため互換 |
| `T_o(a,b)` と `T_o(b,a)` の比較 | 基本因子の順序が異なるため非互換 |
| `T_o(T_o(a,b),a)` の構築 | 親と内部因子の葉の Index が重複するため拒否 |
| 同じ部分空間 `T_o(a,b)` を二つ合成 | 基本因子の Index が重複するため拒否 |
| `T_o()` と `T_g()` の比較 | 共通の一次元の単位空間として互換 |
| `T_o(E,a)`、`T_g(a)`、`a` の比較 | 空・単因子の同一視後に同じ基本因子を持つため互換 |
| `T_o(a,T_g(b,c))` の構築 | 同一視後にも積の種類の混在が残るため拒否 |
| 一次元の Factor と `E` の比較 | Factor の Index を保存するため非互換 |
| 構築時の因子を配置・演算・空間構成に再利用 | 元の Index と grouping を保持する |

部分系判定の例では、`H = T_o(T_o(a,b),c)` とする。

| 判定対象 `K` | `H` の部分系か |
| --- | --- |
| `T_o(a,b)` | 部分系 |
| `T_o(a,c)` | 木の部分木でなくても、相対順序を保存するため部分系 |
| `T_o(c,a)` | 相対順序が異なるため部分系ではない |
| `H` 自身 | 部分系 |
| 空の Tensor `E` | 部分系 |
| `H` に存在しない Index を持つ Factor | 部分系ではない |
| `T_g(a,b)` | 複数因子の積の種類が異なるため部分系ではない |
| `T_g(a)` | 単因子の同一視により `a` と同じ部分系 |

空の Tensor に対しては、空の Tensor と同一視される空間だけが部分系となる。一次元でも Index を持つ Factor は、その Index を持たない空の Tensor の部分系ではない。

階層的な評価と基本因子へ展開した評価は、積の種類と因子順序を保持した場合に一致する。基本因子の列挙で全基底・状態式・行列の生成が起きないことも検証対象とする。

`H = T_o(T_o(a,b),c)` の直下の子の列挙は `[T_o(a,b),c]`、所属 Index の列挙は `[Index(a),Index(b),Index(c)]` とする。単因子 Tensor は保持する子を一つ列挙し、空の Tensor は子も Index もゼロ個とする。列挙で既存の Index・順序・括り方が変わらないことを確認する。

## 9. 未決事項

構築・合成、構成要素・Index の列挙、互換性判定、および部分系判定を提供する。任意の因子集合を指定した部分系の取得と Regrouping は提供しない。これらを本書の合意済みの契約とし、Hilbert 空間の数学的な基本契約について、現時点で追加の未決事項はない。

公開型 `HilbertSpace` への統一、Dof の葉への統合、`Arc` による共有、非公開の節点 enum、借用による列挙、ordinary な葉列の照合、反復的な走査・破棄、および次元取得時の overflow 検出は採用済みとする。段階3の検証では、空間・束縛・エラー、総次元の overflow、および2万段の空間の128 KiB スタック上での処理を確認した。

Step 3 で空間・constructor・列挙・判定・次元取得の API、非公開の節点と作業スタックを実装した。LocalStructure 廃止後の公開型は `HilbertSpace<I>` とする。Index は照合するメソッドで `Eq` を要求し、重複検出には等値比較を用いる。現在の API は第10節に記録する。

Chain と、利用側が保持する Site 順の局所空間列・基本自由度の対応は、[Chain の配置契約](geometry.md)に従う。

これらの決定では、再帰的な空間構成、基本因子の Index・順序の保存、および保持する木構造と数学的な互換性の区別を維持する。

## 10. 現在の公開 API

| API | 結果・制約 |
| --- | --- |
| `HilbertSpace::factor(index)`・`unit()` | 既存 Index の Factor・空 Tensor を構築する。Index の Clone・Default は不要 |
| `tensor_all_with(kind, children)`・`tensor_all(children)`・`tensor_all_g(children)` | 種類指定・ordinary・graded の構築。`I: Eq` で重複と混在を検査し、`Result<Self, SpaceError>` を返す |
| `tensor_with(kind, &rhs)`・`tensor(&rhs)`・`tensor_g(&rhs)` | 同じ三種類の二項版。因子順序は左、右の順 |
| `kind()`・`effective_kind()` | 保持する Tensor 節点と、空・単因子の同一視後の `Option<TensorKind>` を返す |
| `factor_index()`・`children()`・`indices()` | Factor の `Option<&I>`・直下の子・全葉の Index を借用する |
| `is_compatible_with(&rhs)`・`is_subsystem_of(&ambient)` | `I: Eq` で葉順・相対順序と有効な積の種類を照合する |
| `dim::<P>()` | `P: Dimension` の次元を因子ごとに適用し、`Result<usize, DimensionOverflow>` を返す |
| `dim_with(dim_local)` | `NonZeroUsize` の局所次元を明示し、同じ Result を返す |
| `clone()`・`shares_node_with(&rhs)` | Index を複製せず共有ハンドルを複製・比較する |

`dim::<S>()` は状態空間、`dim::<O>()` は演算子空間の次元を返す。Operator は S を保持しないため、その作用先の次元 d が必要なら `dim_with` に局所次元を指定する。空 Tensor の次元は1とし、P の局所次元が overflow する場合も空 Tensor の取得は成功する。空間比較・配置・束縛には primitive 型を登録せず、異なる primitive 型の式が同じ空間構造を利用できる。

`TensorKind::{Ordinary, Graded}` は構築時に明示する。`tensor_all` / `tensor` は ordinary、`tensor_all_g` / `tensor_g` は graded の短縮形で、共通の種類指定 constructor へ委譲する。子や primitive から種類を推測せず、既存空間の kind を書き換える変換は提供しない。`kind()` は Factor なら `None`、Tensor なら空・単因子 wrapper も含めた保持値を返す。`unit()` は ordinary の空 Tensor として保持し、`tensor_all_g([])` と数学的に互換である。

### 10.1. 有効な種類と混在の検証

各ハンドルは葉数の区分（0 / 1 / 複数）と有効な種類だけを構造 metadata として保持する。局所次元・全基底・評価結果のキャッシュではない。空の子を意味上除いた個数を `m` とすると、次のように判定する。

1. `m = 0` は単位空間で、有効な種類は `None`。
2. `m = 1` は唯一の非空の子の有効な種類を継承する。外側の kind は互換性を変えない。
3. `m >= 2` は親の kind が有効な種類となり、有効な種類を持つ子はすべて親と一致する必要がある。

このため `T_o(E,T_g(a,b))` は許可し、`T_o(a,T_g(b,c))` は `SpaceError::MixedTensorKinds` で拒否する。0・1葉の有効な種類はなく、両規則の評価は一致する。複数葉の互換性・部分系判定と評価経路の選択には `effective_kind()` を用いる。

木と子の grouping は保持し、単因子 wrapper は子の metadata を定数時間で継承する。単因子 wrapper の構築のたびに全子孫を走査せず、重複検証にも `I: Eq` 以外の Hash・Ord・Clone を追加要求しない。走査・破棄は反復処理を維持する。

`SpaceError::MixedTensorKinds` の追加により、0.2.0 への移行時には外部の網羅的 match を更新する必要がある。複合式・基底・export は既存の `ExpressionError::Space`・`BasisError::Expression`・`ExportError` の経路でエラーを伝播する。
