# Chain・Site と自由度の配置

更新日: 2026-10-26

## 1. 対象と責任分担

本書は[要件](../requirements.md)に基づき、Chain と Site への自由度の配置を定める。2026-10-26 に、Chain が Site ごとの Index 列と境界条件だけを保持する方針を採用した。以下は移行先の設計であり、製品コード・テスト・利用例の移行は未着手である。[Step 4 の合意事項](../step-4-decisions.md)にある Geometry・HilbertSpace の保持は移行前の記録とする。

`qmb-geometry` は Site の順序、Site ごとの自由度の Index 列、境界条件と幾何関係の列挙を担当する。Hilbert 空間の構築・保持・grouping・全体の因子順序は利用側が管理する。局所式の配置・束縛・評価は `qmb-core`、UUID の具体型 `Id128` は `qmb-basic`、XYZ・Heisenberg・Neel の構築 utility は `qmb-spin-half`、再公開・縦断例は `qmb` が担当する。

保持構成は次のとおりとする。Geometry 構造体は廃止し、Chain が直接データを所有する。代わりの配置用ラッパー構造体は設けない。`qmb-geometry` crate と関係ごとの能力 trait は維持する。

```rust
struct Chain<I> {
    sites: Vec<Vec<I>>,
    boundary: Boundary,
}
```

Chain は局所代数・HilbertSpace・Grading の情報を保持しない。Grading の情報保持・API は fermion 系の実装時に具体化する。core は Chain を必須としない。

## 2. Site の指定と順序

Site は外側の Vec の位置 `0,1,…,N−1` を表す `usize` で指定し、内側の Vec は Site 内の自由度の順序を表す。SiteId・共有識別トークンは設けず、Chain を跨ぐ Site の同一性・混用を管理しない。範囲外は検出する。構築後の配置と順序は不変とする。

Site ordering と tensor ordering は区別する。Site 内の Index 順序は関係の列挙に使う順序であり、既存の全空間の因子順序や grouping を暗黙に変更する指定ではない。

## 3. Site と Hilbert 空間・Index の対応

各 Site は0個以上の Index を持つ。既存の空間を利用する場合、利用側がその空間から必要な Index を読み取り、Site ごとの列として渡す。Tensor の grouping や Factor ハンドルそのものを Chain に保持させない。

| 入力の例 | Site に属する Index |
| --- | --- |
| 一つの自由度 | `[a]` |
| 複数の自由度 | `[a, b]` |
| 空 Site | `[]` |

一つの Chain 内では各 Index をちょうど一つの位置に割り当て、Site 内・Site 間の重複を構築時に `I: Eq` で拒否する。別の Chain や空間で同じ Index を再利用することは認める。Index と Site の対応を全体で登録する機構は設けない。

## 4. 二つの構築経路と既定の配置規則

| 経路 | 契約 |
| --- | --- |
| 既存 Index からの配置 | Site ごとの所有 Index 列を Site 順に受け取り、その順序を保持。Index は移動して取り込める |
| 個数からの生成 | Site 数と各 Site 共通の自由度数を受け取り、Site 順・Site 内の順に Index を発行 |

新規 Index の発行が必要な箇所でのみ `I: Default` を要求し、`I::default()` を呼ぶ。この用途の Default は新しい独立した Index を発行する契約とする。UUID・乱数源は具体型側に閉じ、Chain は basic に依存しない。専用の Index 生成 trait や発行関数の引数は設けない。保持・借用列挙には Default・Clone を要求しない。

所有 Index 列の複製や、借用した Index から別途 Factor を構築する操作では `I: Clone` を要求できる。Clone／Copy は同じ識別値の複製であり、新規発行ではない。標準の `Id128` は Copy を提供する。Clone できない Index を複数の所有対象で共有するための追加の仕組みは今回必須にしない。自由度数は Index の個数であり、Hilbert 空間の次元とは区別する。

### 4.1. 初期版の 1D Chain の配置規則

標準の幾何は open / periodic の Chain のみとする。外側の Vec が Site 順序、各 Index 列が Site 内の自由度順序を表す。

```text
入力: [[a,b], [c]]
site 0 → [a,b]
site 1 → [c]
```

個数からの経路では例えば Site 数2・各 Site の自由度数2から `[a,b]`, `[c,d]` の順に Index を発行する。利用側がこの列順に全空間を新規構築すると基本因子順序は `[a,b,c,d]` となる。既存の全空間を利用する経路では、その因子順序を維持して配置との対応を検証する。

## 5. 配置の完全性と空の場合

- Chain は Index の重複、範囲外の Site、Site 数を別途指定した経路での入力数不一致を検出する。
- 外部の全空間との対応・未配置・過不足は、その空間と Chain を接続するモデル構築側・利用側が検証する。
- 自由度数ゼロの Site は空の Index 列として保持し、Site や幾何関係を削除しない。
- Site 数ゼロでは Site 列が空の Chain を作る。単位空間は必要な利用側が構築する。
- 入力不足を理由に空の Site を補わない。

## 6. Grouping と自由度の参照

HilbertSpace の grouping と全体の因子順序は利用側が保持する。Chain は Site と Site 内の位置から Index を参照できればよく、内部ラベルの型・名前による指定は導入しない。空間との接続は Index の等値比較で行い、同じ Arc 節点の共有を必須条件にしない。core の式・空間内部の不変な共有は維持する。

### 6.1. 空間への変換能力 ToSpace

2026-10-26 に、空間への変換を独立した拡張 trait `ToSpace`、メソッド `to_space` で提供する方針を採用した。Chain と HilbertSpace の両方に依存する `qmb` の接続用 module に置き、qmb から公開する。Chain の保持構造や `NearestNeighborRelation` の必須能力には含めず、qmb-geometry に空間への依存を戻さない。実装は未着手である。

ordinary な空間に対するシグネチャは次の形を想定する。

```rust
pub trait ToSpace {
    type Index;

    fn to_space(&self) -> Result<HilbertSpace<Self::Index>, SpaceError>;
}
```

Chain に対しては `I: Clone + Eq` のときに実装する。Clone はこの変換実装で要求し、trait 全体や Chain の保持・借用列挙の条件にしない。既存の空間を共有保持する別の型は、その Index の Clone を要求せずに実装できる。

変換は元の対象を借用し、Index の識別値を維持する。Chain の実装は各 Site の Index から Factor と Site ごとの空間を構築し、Site 順・Site 内の Index 順で全空間を構成する。空の Chain は単位空間になる。同じ Chain から繰り返し変換した空間は同じ Index・因子順序を持つが、同じ Arc 節点を返す保証は設けない。生成した空間の保持は呼び出し側が担当する。

Chain の Index 列に含まれない任意の grouping や、Site 順と異なる既存全空間の因子順序を復元する API ではない。既存空間を利用する場合は、その空間を利用側で保持し、モデル構築時に Chain と併せて渡して対応を検証する。

変換契約は標準の Into／From に固定せず、独立した能力として定める。将来 Grading のために必要となる情報の渡し方と変換 API の拡張は、fermion 系の実装時に具体化する。現時点で Grading の型や変換コンテキストを先行導入しない。

## 7. 関係の列挙とモデル構築の境界

Chain は Site ごとの Index 列の借用と、`nn_site_pairs()` による近接 Site ペアの列挙を提供する。近接 Site ペアは正規化した番号の辞書順とし、self-loop・重複辺を生成しない。長さ0・1では辺なし、長さ2では open / periodic とも `(0,1)` 一つとする。

関係ごとの trait を設け、`NearestNeighborRelation` の `nn_pairs()` は近接 Site 間の Index 全組合せ、`nn_pairs_by_pos()` は同じ位置同士を短い方まで列挙する。Index は Chain から借用し、Clone・Default やペア全体の Vec は要求しない。片側が空ならペアなしとし、Site ペア順と Site 内の Index 順序を保持する。共通関数も Index 列の組合せ・位置対応として整理し、HilbertSpace を引数とする契約を移行する。

演算子・係数の選択、Factor・全空間の構築または既存空間との接続、局所式の配置・束縛、状態の構築はモデル utility・利用側が担当する。既存空間との接続では必要な Factor の存在を検証する。Neel を全空間上に構築する経路では、Site ごとの Index が全因子に対応することを確認し、Up／Down は Site 順で決め、状態は全空間の因子順で構築する。

XYZ・Heisenberg・Neel utility の一 Site 一 spin の制約、0-site の zero Hamiltonian・単位状態、解析的期待値は維持する。空間の構築や接続で新しい Index を発行し直さない。既存空間を明示する経路と Index から空間を構築する経路の具体的な関数シグネチャは、実装移行時に定める。

次近接などの追加の関係 trait と追加の幾何型は、その機能を導入するときに設計する。Geometry 構造体の保持は拡張の前提条件にしない。

## 8. 構築後の不変性

Chain の Site 順序・関係・配置は構築後に変更しない。異なる幾何や配置が必要なら新しい対象を作り、既存の空間・式を変更しない。新規発行と既存 Index の再利用を区別し、参照や Clone だけで再発行しない。

## 9. 検証対象

以下は移行後に実施する検証であり、成功済みの記録ではない。

- Clone・Default のない Index の所有入力・借用列挙、Site 内・Site 間の重複拒否、範囲外・Site 数不一致、空・可変長・複数自由度の列。
- 新規発行の回数と順序、参照・複製・空間接続による識別値の保持。
- open／periodic の小規模規約、二つの Index ペア方式の順序と空の扱い、別の幾何型での関係 trait の再利用。
- 利用側での空間との対応検証、Site 順と異なる全空間の因子順序、利用側で保持する grouping、既存 Index の再利用。
- ToSpace の変換前後での Index の同一性と順序、空の Chain、繰り返し変換した空間の互換性、元の Chain の不変性。Clone の要求が変換実装に限定されること。
- 一 Site 一 spin の制約、0-site、両精度の Neel の norm、XYZ・Heisenberg の Hermiticity と独立行列・解析的期待値との一致。

## 10. 実装移行の境界

Geometry とその再公開、`Chain::from_geometry`・`geometry()`・`space()`、HilbertSpace を受け渡す Site API は新しい責任分担へ移行する。空間の構築・参照・完全性検証を利用側へ移し、テスト・facade・実行例を追従させる。エラーは Chain の入力検証と空間接続側の検証に分け、具体名・constructor 名・借用 API の詳細は実装時に定める。

移行前の「Chain の Clone に I: Clone は不要」「Site と全空間の Arc 節点を Chain が共有する」という保証は、所有 Index 列の契約へ置き換える。空間接続時にだけ必要な Clone を、Chain の保持・借用列挙へ波及させない。`qmb-geometry` の core 依存は空間 API の移動に合わせて棚卸しし、不要になった依存を除く。
