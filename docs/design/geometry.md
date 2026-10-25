# Geometry・Site と自由度の配置

更新日: 2026-10-25

## 1. 対象と責任分担

本書は[要件](../requirements.md)と[Step 4 の合意事項](../step-4-decisions.md)に基づき、Geometry と Site への自由度の配置を定める。公開 API と受入条件は Step 4 の合意事項を参照する。

`qmb-geometry` は Site の順序・幾何関係・配置・空間を保持し、局所式の配置・束縛・評価は `qmb-core` に委譲する。UUID の具体型 `Id128` は `qmb-basic`、Heisenberg・Neel の構築 utility は `qmb-spin-half`、再公開・縦断例は `qmb` が担当する。

Geometry は生成時に利用者が指定した `L: LocalStructure` を用い、Site 順の `Vec<HilbertSpace<I, L>>` と全空間を保持する。独立した Placement 型は設けない。空間の既存の Arc 共有を使い、Index や木全体を複製しない。core は Geometry を必須としない。異種局所構造の混在・graded は初期版の対象外とする。

## 2. Site の指定と順序

Site は保持する Vec の位置 `0,1,…,N−1` を表す `usize` で指定する。SiteId・共有識別トークンは設けず、Geometry を跨ぐ Site の同一性・混用を管理しない。範囲外は検出する。構築後の配置と順序は不変とする。

Site ordering と tensor ordering は区別し、既存空間の順序を配置によって暗黙に変更しない。

## 3. Site と Hilbert 空間・Index の対応

各 Site は Factor・Tensor のどちらも保持でき、その空間に属する全葉の Index を対応付ける。入れ子の順序・grouping と既存の Index を保持する。

| 空間 | Site に属する Index |
| --- | --- |
| Factor | その Index 一つ |
| Tensor | 葉を左からたどった Index 列 |
| 空の Tensor | ゼロ個 |

一つの配置では各 Index をちょうど一つの Site に割り当て、Site 間で重複する Index は拒否する。別の配置での再利用は認める。Index と Site の対応を全体で登録する機構は設けない。

## 4. 二つの構築経路と既定の配置規則

| 経路 | 契約 |
| --- | --- |
| 既存空間からの配置 | Site ごとの HilbertSpace を Site 順に受け取り、既存の Index・順序・grouping を保持 |
| 個数からの生成 | Site 数と各 Site 共通の自由度数を受け取り、生成時に指定した Index 型と局所構造型で空間を構築 |

新規 Index の発行が必要な箇所でのみ `I: Default` を要求し、`I::default()` を呼ぶ。この用途の Default は新しい独立した Index を発行する契約とする。UUID・乱数源は具体型側に閉じ、Geometry は basic に依存しない。専用の Index 生成 trait や発行関数の引数は設けない。既存空間の保持や列挙には Default・Clone を一律に要求しない。

全空間の取得で Index を再発行したり空間を再構築したりしない。自由度数は葉の個数であり、Hilbert 空間の次元や Tensor 直下の子の個数とは区別する。

### 4.1. 初期版の 1D Chain の配置規則

標準の Geometry は open / periodic の Chain のみとする。Site ごとの空間を保持する Vec の外側が Site 順序、各要素の葉の順序が Site 内の自由度順序となる。

```text
入力: [Tensor(a,b), c]
site 0 → [a,b]
site 1 → [c]
```

個数からの経路では例えば Site 数2・各 Site の自由度数2から `[a,b]`, `[c,d]` の順に Index を発行する。新しい全空間の基本因子順序は `[a,b,c,d]` とする。既存全空間を受け取る場合は、その因子順序を変更しない。

## 5. 配置の完全性と空の場合

- 既存空間の全 Index を対応付け、重複・未配置を許さない。必要数を別途指定した経路では過不足をエラーとする。
- 自由度数ゼロの Site を空の Tensor として保持し、Site や幾何関係を削除しない。
- Site 数ゼロでは空の Geometry と単位空間を作る。
- 入力不足を理由に空の Site を補わない。

## 6. Grouping と自由度の参照

初期版は HilbertSpace をそのまま grouping に用い、構成要素・葉の順序と保持した Factor を参照する。内部ラベルの型・対応表・名前による指定は導入しない。

## 7. 関係の列挙とモデル構築の境界

Geometry は `sites()` と `nn_site_pairs()` を公開する。近接 Site ペアは正規化した番号の辞書順とし、self-loop・重複辺を生成しない。長さ0・1では辺なし、長さ2では open / periodic とも `(0,1)` 一つとする。

関係ごとの trait を設け、`NearestNeighborRelation` の `nn_pairs()` は近接 Site 間の Index 全組合せ、`nn_pairs_by_pos()` は同じ位置同士を短い方まで列挙する。二つの空間を借用する共通関数 `index_pairs`・`index_pairs_by_pos` を利用し、Index の Clone・Default やペア全体の Vec は要求しない。片側が空ならペアなしとし、Site ペア順と葉の順序を保持する。

関係の列挙は Hamiltonian の項の生成とは区別する。演算子・係数の選択と状態の構築はモデル utility が担当する。次近接などの追加の関係 trait と追加 Geometry は、その機能を導入するときに設計する。

## 8. 構築後の不変性

Geometry の Site 順序・関係・配置は構築後に変更しない。異なる幾何や配置が必要なら新しい対象を作り、既存の空間・式を変更しない。新規発行と既存 Index の再利用を区別し、参照や Clone だけで再発行しない。

## 9. 検証対象

Factor・入れ子・空・複数自由度の配置、全 Index の対応、重複・未配置・過不足・範囲外の拒否、発行と共有の保持を確認する。内部ラベルと Geometry を跨ぐ所属確認は対象外とする。

Chain の小規模規約・順序、二つの Index ペア方式の可変長・空・入れ子の葉順、別 Geometry のテスト実装での関係 trait の再利用を検証する。モデル utility では一 Site 一 spin Factor の制約、Neel の norm、Hermiticity、解析的期待値を確認する。具体的な受入条件は [Step 4 の検証](../step-4-decisions.md#8-受入条件と検証)に従う。
