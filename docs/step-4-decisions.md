# Step 4 着手前の判断事項と推奨案

作成日: 2026-10-25

本書は、Step 3 まで実装済みの qmb-rs における Step 4 の事前検討を一時的に記録する。第1・2節、第3.1節および第4〜8節は対話で合意済みであり、第3節の配置の構築経路・Index 発行・内部ラベルは採用前の推奨案である。いずれも実装済みの API を示すものではない。

第3節の推奨案と、合意済みの方針で未決としている具体的な API は引き続き検討する。

## 1. Geometry・Site の型と同一性（合意済み）

- Site は Geometry 内の `Vec` の位置を表す `usize` で指定する。専用の `SiteId` 型や Site ハンドルは設けない。
- Geometry を跨ぐ Site の同一性を定義・管理せず、別 Geometry 由来の位置番号の混用も検出・拒否しない。位置番号は利用する Geometry の文脈で解釈し、混用を避ける責任はユーザに置く。
- 所属を追跡する共有トークン、UUID、グローバル採番器は追加しない。
- 範囲外の位置番号など、対象の Geometry の情報から検出できる不正は検証する。
- 構築後は配置と順序を変更しない。

## 2. Geometry と配置情報の保持・自動配置（合意済み）

- Geometry が配置情報も保持し、独立した `Placement` 型は設けない。初期版の Chain は境界条件、Site 順の配置、全 Hilbert 空間をまとめて保持する。
- 入力された `Vec` の各要素を一つの Site に対応付け、並び順に自由度を自動配置する。ユーザに別の Site 対応表を要求しない。
- 各要素は Factor・Tensor を扱え、要素内の因子順序・grouping と既存の Index を保持する。例えば `[a,b,c]` は Site ごとに一因子、`[Tensor(a,b),c]` は Site 0 に二因子、Site 1 に一因子を配置する。
- 具体的には `Vec<HilbertSpace<I, L>>` として保持する。既存の `HilbertSpace` が内部の節点を `Arc` で共有するハンドルなので、追加の `Arc<HilbertSpace<…>>` は設けない。Site の空間と全空間で節点を共有し、Index や空間の木全体を複製しない。ただしハンドルの格納領域、節点の割り当て、参照カウントのコストは残る。
- 全空間を新規構築するときは、この Site 順で tensor product を取る。
- 既存空間を受け取る経路では、元の全空間と因子順序を保持する。Site 順に並べ直して再構築しない。
- 局所式の `at(&factor)` に必要な葉を、保持したハンドルから参照できるようにする。
- Chain の近接関係は、配置の長さと境界条件から位置番号の組として列挙する。

幾何・配置・モデル構築の責任は区別するが、それぞれを独立した公開型にすることは要求しない。

## 3. 配置の構築経路・Index 発行・内部ラベル

既存契約に対応する次の三つの経路を提供する。

1. Site ごとの既存 Factor・Tensor を配置する。
2. Site ごとの Index 列から、明示した局所構造で空間を構築する。
3. 自由度数の列と `FnMut() -> I` を受け取り、Index 列を先に発行する。

発行済み Index から Factor への変換は所有権を移動して行い、`I: Clone` を一律に要求しない。検証も基本は `Eq` に限定する。

内部ラベルはジェネリックな型とし、「ラベル→Site 内の自由度順序番号」の対応を保持する。Index を複製せず、ラベルなし・一部だけラベルありの配置も扱う。

### 3.1. 幾何関係の trait 化（方針・近接関係の公開 API は合意済み）

Geometry が列挙する Site の幾何関係を利用し、関係ごとの trait が対応する Index の組を列挙する。複数の Geometry 実装に同じ関係 trait を実装し、利用側が Geometry の具体型によらず同じ API を使えるようにする。Geometry と配置を別の保持型に分離する必要はない。

近接関係では、次の流れとする。

1. Geometry が近接 Site の位置番号 `(usize, usize)` を列挙する。位置番号は Site を保持する `Vec` の位置に対応する。
2. 近接関係を表す `NearestNeighborRelation` trait を設け、`nn_pairs()` で近接 Site 間の全組合せの Index ペアを列挙する。同じ位置同士の列挙は `nn_pairs_by_pos()` で提供する。
3. Chain などでの trait の実装は、Geometry が返す Site ペアと各 Site が保持する `HilbertSpace` を利用し、葉の Index を参照して Index ペアを構築する。

`NearestNeighborRelation` と、そのメソッド `nn_pairs()`・`nn_pairs_by_pos()` を採用する。公開 API は以下の形とする（コードへの実装は未着手）。

```rust
pub trait NearestNeighborRelation {
    type Index;

    /// 近接 Site 間の Index の全組合せを列挙する。
    fn nn_pairs(
        &self,
    ) -> impl Iterator<Item = (&Self::Index, &Self::Index)> + '_;

    /// 同じ位置の Index 同士を、短い方の長さまで列挙する。
    fn nn_pairs_by_pos(
        &self,
    ) -> impl Iterator<Item = (&Self::Index, &Self::Index)> + '_;
}
```

- 関係ごとに trait を設け、次近接なども対応する関係 trait として追加できるようにする。未対応の関係をすべての Geometry に必須化しない。
- 一つの Site に複数の自由度がある場合も、Site ペアから Index ペアへの変換側で対応する。trait の出力は Site ごとの Index 列の組ではなく、Index ペアとする。
- 複数自由度間のペア選択規則として、全組合せを取る方式と、並び順が同じ位置の要素のみを取る方式の二つを提供する。具体的な規則は下記のとおりとする。
- Index は保持した空間から借用でき、列挙のために複製・再発行しない。
- 関係の向き、順序、重複の扱いは関係ごとに定める。Chain の近接 Site 関係には既存の辞書順・self-loop と重複辺を作らない規約を適用する。
- Site 番号は対象の Geometry 内で解釈する。Geometry を跨ぐ Site の同一性や混用の検出は追加しない。
- 具体的な追加 Geometry・関係の実装は初期版に要求しない。

複数自由度間のペア選択規則は次の二方式を採用する。Site 内の順序には、保持した `HilbertSpace` の葉を左からたどった Index の順序を用いる。

| 方式 | 規則 | 左 Site が `[a,b]`、右 Site が `[c,d,e]` の例 |
| --- | --- | --- |
| 全組合せ | 両 Site の Index 列の直積を列挙する。左側を外側、右側を内側のループとして順序を定める | `(a,c), (a,d), (a,e), (b,c), (b,d), (b,e)` |
| 同じ位置同士 | 同じ位置の Index 同士を先頭から組にする。長さが異なる場合は短い方の長さで打ち切り、余った要素は使用しない | `(a,c), (b,d)` |

どちらの方式も片側の Site の自由度がゼロなら Index ペアを返さない。自由度数が異なること自体はエラーにしない。これは関係の列挙規則であり、配置時に別途指定された必要自由度数の検証とは区別する。いずれの方式も近接 Site ペアの列挙順を維持する。

近接関係の trait 名、両メソッドの名称・シグネチャ、および `nn_pairs()` が全組合せ、`nn_pairs_by_pos()` が同じ位置同士を返す挙動は確定する。Geometry 側の Site ペア列挙や変換処理の具体的な API は未決とする。

## 4. Heisenberg・Neel 構築補助の配置先と対象範囲（合意済み）

Heisenberg Hamiltonian と Neel 状態の構築補助は、`qmb-spin-half` の utility として実装する。`qmb` に構築ロジックを置く従来の案から変更し、`qmb` は利用窓口として再公開や縦断例を担当する。

- 最初の便利 API は、各 Site がちょうど一つの spin Factor を持つ配置を対象とする。
- 空 Site・複数自由度の Site は Geometry の機能として扱うが、このモデル API では明確なエラーにする。
- 複数自由度から対象を選ぶ API は、明示的な選択方法を設けて追加できるようにする。
- Heisenberg は実数 `J: T` を受け取り、通常の演算子積と和から構築し、配置の全空間へ bind して返す。
- Neel は先頭の `Up`／`Down` を指定可能にし、同じ配置から状態を構築する。
- どちらの構築時にも Index を再発行しない。
- 奇数長の periodic chain でも Neel 積状態の構築を許可し、端点間が同じ向きになるものとする。
- 0-site chain は許可し、空の積である単位状態と、全空間へ束縛した zero Hamiltonian を構築する。

## 5. norm・正規化の数値 API と平方根（合意済み）

局所式・複合式に共通して、次の API を設ける。

- `norm_squared`
- `norm`
- `normalized`
- `normalized_expectation`

複合式には、既存の空間選択に合わせて `normalized_expectation_in` も設ける。

`norm_squared` と `norm` は実数型 `T` を返す。自己内積が実数になる数学的契約に基づいて実部を使う方針を明文化する。負値のゼロへの補正や、近似許容誤差によるゼロ判定は導入しない。

正規化期待値は定義どおり自己内積で除算し、平方根を必要としない実装にする。

平方根の依存・feature 方針は以下を採用する。

- `qmb-core` に任意の `libm` feature を設け、`num-traits/libm` を有効にする。core 自体の既定 feature では有効にしない。
- 利用窓口の `qmb` では既定で core の `libm` feature を有効にし、依存を最小限にしたい利用者は無効にできるようにする。
- `norm()`・`normalized()` は `libm` feature 有効時に提供し、平方根には `num_traits::Float` を利用する。`Float` 制約はこれらのメソッドに限定する。
- 平方根を使わない `norm_squared()`・`normalized_expectation()` とその `_in` 版は、`libm` feature 無効時にも提供する。既存の作用・内積やこれらの操作に `Float` 制約を広げない。
- `no_std + alloc` を維持する。`libm` は Rust 製の数学関数ライブラリであり、OS の数学共有ライブラリの別途インストールは要求しない。
- feature 有効時に `num-traits` 経由の production 依存として `libm` が追加されることを了承済みとする。

core の feature 設定は次の形とする（manifest への反映は未着手）。

```toml
[features]
default = []
libm = ["num-traits/libm"]
```

確認した現在の依存設定では `num-traits` の `std`・`libm` は無効である。Cargo の feature は依存グラフ内で加算されるため、別の依存元が `num-traits/std` を有効にした場合はそれも有効になる。確認した `num-traits 0.2.19` では、`std` と `libm` が両方有効なら `Float` は `std` 側の実装を使う。今回の `libm` feature は数学関数を利用可能にするためのもので、常に libm の実装を強制する指定ではない。

追加後のビルド時間・バイナリサイズ・実行速度は未計測であり、採用する依存バージョンでの MSRV と `no_std` 向けビルドも実装時に確認する。

## 6. zero-norm とエラー型（合意済み）

- `normalized` と `normalized_expectation` は、計算された自己内積の厳密なゼロを拒否する。
- 明示的な zero に加え、`ψ − ψ` のような相殺による zero も検証する。
- 局所式・複合式で共通の `NormalizationError` を設け、`ZeroNorm` と既存の `ExpressionError` のラップを持たせる。
- `norm(zero)` 自体は `0` を返す。
- NaN・infinity などは、既存要件どおり数値型の挙動に委ねる。

## 7. bra の延期と braket の計算（合意済み）

- bra 型・ラッパー・構築メソッドは初期版では導入しない。Step 4 および初期版の完成条件にも含めない。
- braket の計算はすべて既存の内積に委ねる。`⟨φ|ψ⟩` は `phi.inner_product(&psi)` で計算する。
- `⟨φ|A|ψ⟩` は状態への作用と内積を組み合わせて計算する。期待値・正規化期待値も、作用と内積から導出する既存方針に従う。
- 左引数についての共役線型性、右引数についての線型性、複合状態の空間互換性は内積の契約として扱う。
- bra を独立した対象として提供する API は後続へ延期し、具体的な型や操作は今回設計しない。

## 8. Step 4 の受入テストと公開 API（合意済み）

`qmb-geometry` は `no_std + alloc` で実装し、`qmb` は各 crate と spin crate の構築 utility の再公開、および利用例・縦断的な結合テストを担当する。公開名の衝突を避けるため、各 crate に対応する module を設け、crate root の再公開は主要型に絞る。

完了条件には次を明記する。

- 両精度・両方の Neel 開始向きで、open chain の期待値 `−J(L−1)/4` と norm `1` を確認する。
- 0-site の状態は単位状態とし、norm `1`、Hamiltonian と期待値は `0` とする。
- periodic の長さ `0,1,2`、辺の辞書順、重複辺・self-loop の不在を確認する。
- 複数自由度・空 Site・入れ子空間・ラベル、範囲外の Site 指定、Index の重複や入力数の不一致を検証する。別 Geometry 由来の Site 番号の拒否は検証対象にしない。
- Chain と非公開の別 Geometry のテスト実装で、同じ関係 trait を通じて Index ペアを列挙できることを検証する。Geometry による近接 Site ペアの列挙と、Site ペアから Index ペアへの変換をそれぞれ確認する。全組合せ・同じ位置同士の二方式について、複数自由度、左右で異なる自由度数、片側が空、入れ子 Tensor の葉の順序、およびペアの列挙順を検証する。
- braket は内積の左引数に関する共役線型性・右引数に関する線型性と、作用を含む内積の結果で検証する。非正規直交な例の norm、zero-norm の拒否も確認する。
- `libm` feature の有効・無効の両構成を確認する。無効時にも平方根を使わない操作が利用できること、有効化が `std` を要求しないこと、および `qmb` の既定構成で norm・正規化を利用できることを検証する。

## 決定済みの事項・後続段階との境界

既存の `docs/design/geometry.md` と `docs/requirements.md` にある、Geometry を跨ぐ Site の同一性や別 Geometry の Site の拒否という契約は、第1節の合意に従って見直す。本メモへの記録時点では、これらの既存文書は未更新である。

既存の Geometry 設計で述べている Site ごとの Index 列の組の取得と、今回合意した関係 trait による Index ペアの列挙との関係も、第3.1節の方針および全組合せ・同じ位置同士のペア選択規則に合わせて整理する。

既存の実装計画・アーキテクチャ・Geometry 設計にある、Heisenberg・Neel の構築補助を `qmb` に置く責任分担は、第4節の合意に合わせて見直す。幾何関係の trait を `qmb-spin-half` の utility から利用するための依存関係も合わせて整理する。これらの既存文書・manifest は本メモへの記録時点では未更新である。

既存の実装計画・初期版の提供範囲・関連設計文書にある、bra を Step 4 および初期版で提供するという要求は、第7節の合意に従って後続へ延期する。本メモへの記録時点では、これらの既存文書は未更新である。

Chain の配置順序と periodic の小規模規約は決定済みであり、今回改めて決める対象ではない。

Partial trace・基底処理・数値 backend は Step 5・6 の担当であり、今回の事前決定には含めない。

## 確認した資料と検証範囲

- `docs/implementation-plan.md`（特に Step 4 の実装内容と完了条件）
- `docs/step-3-design.md`
- `docs/design/geometry.md`
- `docs/design/initial-release.md`
- `docs/design/primitive-expression.md`
- `docs/requirements.md`
- 現在の workspace・core の依存設定、局所式・複合式・空間・Index の実装
- ローカルにある `num-traits 0.2.19` のソースと feature 定義

今回の検討では実装コードを変更しておらず、テストも実行していない。本ファイルは合意済みの方針と未採用の推奨案の記録であり、実装完了を意味しない。
