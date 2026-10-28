# QuantumSystem を入口とする式構築 API

作成日・更新日: 2026-10-28
状態: 構築入口の統一方針と、初回は scaled の変更および優先順位1〜3を実装する範囲を合意済み。既存の state・operator 以外は未実装。以下の公開名は入力形式を区別するための再検討後の推奨案であり、まだ確定していない。

## 1. 方針と責任分担

構築補助が StateOf・OperatorOf などの型ごとに増えることで、利用者が処理の所在を覚える負担を増やさないよう、次の使い分けを基本にする。

- 新しい式の構築・複数の式の集約は `QuantumSystem` の関連関数から行う。
- 構築済みの式への操作・計算は、その式のインスタンスメソッドから行う。
- `LocalStateOf`・`LocalOperatorOf`・`StateOf`・`OperatorOf` は、主に戻り値や変数の型注釈で使う。

例えば、構築は `Spin::state_tensor_from_primitives`・`Spin::try_operator_sum`、構築後は `.scaled(...)`・`.at(...)`・`.bind(...)`・`.act(...)`・`.expectation(...)` とする。式の演算子構文も引き続き利用できる。複数の式を受け取る集約用関連関数と、既存式を起点にした `.add()`・`.product()`・`.tensor_product()` は併存する。

共通の構築入口は QuantumSystem の既定メソッドとし、既存 core constructor と共通の検証・構築処理へ委譲する。既存 constructor は維持する。QuantumSystem に空間・Index・幾何を登録・保持せず、系の識別子や新たな式ラッパーを返り値に追加しない。必要な能力は各メソッドに限定し、trait 全体に Clone・数値演算・評価能力を要求しない。

Heisenberg・Neel などの代数固有モデルは引き続き SpinHalf<C> の関連関数とし、QuantumSystem の共通メソッドに含めない。macro は後続で検討し、現段階では通常の関連関数を基本 API とする。`Mul` の変換入力・実数入力への拡張は保留とし、[scaled と try_sum の計画](expression-construction.md)の `scaled(impl Into<C>)` は維持する。

## 2. 命名と入力の基本規則

- 基本名は `state_...`・`operator_...` に揃える。既存の `state(primitive)`・`operator(primitive)` は維持する。
- 局所式の追加 constructor・集約は `local_state_...`・`local_operator_...` で区別する。
- `Result` の列を受け取る集約には `try_` を先頭に付ける。通常入力でも空間検証などで Result を返す場合があるため、try の有無は返り値だけで区別しない。
- 列入力は配列固定にせず IntoIterator とする。構築済みの式は所有ハンドルとして受け取り、借用した式の共有が必要なら利用側でハンドルを clone できる。
- 返り値の primitive・係数型は Self の関連型で固定し、Index 型 I は配置先・入力式・明示空間から決める。

### 2.1. primitive 入力と構築済みの式入力の命名再検討

同じ数学的な操作には同じ基本名を使い、primitive を配置してから合成する入口に `_from_primitives` を付ける形を推奨する。構築済みの式を受け取る基本名は `state_tensor`・`operator_product` とする。以下の表・コードはこの推奨案に揃えるが、命名はまだ未確定である。

| 候補 | 状態の primitive 配置／構築済みの式 | 演算子の primitive 配置／構築済みの式 | 評価 |
| --- | --- | --- | --- |
| 従来の仮称 | `state_product`／`state_tensor` | `operator_product`／`operator_product_all` | 名前から入力形式を予測しにくい。状態の二つはどちらも tensor 合成であり、all も要素数の違いに見えるが、両方とも列入力である |
| `_at` で区別 | `state_tensor_at`／`state_tensor` | `operator_product_at`／`operator_product` | 短く、配置を表せる。一方 primitive 限定と分からず、局所式を配置する入口にも読める。既存 at は単一の葉を受け取り、こちらは配置の組の列を受け取る点も異なる |
| `_from_primitives` で区別 | `state_tensor_from_primitives`／`state_tensor` | `operator_product_from_primitives`／`operator_product` | 長いが、型・操作・入力形式を明示でき、両者に同じ命名規則を適用できる。推奨 |
| `_from_assignments` で区別 | `state_tensor_from_assignments`／`state_tensor` | `operator_product_from_assignments`／`operator_product` | 配置の組を受け取ることは示せるが、primitive か局所式かを示せず、assignment という追加の用語の説明も必要になる |

推奨案では、基本名が数学的な操作、接尾辞が変換元を表す。state_product という別の基本名を覚えず、演算子の product にも all を付けて入力形式を区別する必要がなくなる。既存の from_linear_combination と同様に、入力から式を構築することを from で表せる。

`_from_primitives` の入力は primitive の裸の列ではなく、`(&factor, primitive)` の列である。rustdoc の冒頭と引数の型・利用例で配置先も必須と示す。構築済みの式は既に配置情報を持つため、基本名の入口には配置先を重ねて渡さない。

try 版にも同じ基本名を使い、`try_state_tensor`・`try_operator_product` とする。primitive の配置込み入口は通常の組の列を受け取り、配置・合成の失敗を Result で返す。今回、Result 化された配置の組を受け取る別の try 版までは追加しない。

QuantumSystem の関連関数には二項 operator_product を設けないため、列入力であることを示す all は不要と判断する。基礎式型の既存 `Operator::product_all` は改名せず、QuantumSystem の operator_product が委譲する。後続候補の局所式の積も `local_operator_product`・`try_local_operator_product` に揃える。

仮称 product_state・state_product・operator_product_all からの整理は、未実装の QuantumSystem メソッドの命名検討であり、既存公開 API の改名ではない。長い名前を許容する代わりに、入力形式ごとの別名を重ねて公開せず、入口の選択を明確にする。

## 3. primitive から構築する入口

下表の Q は選択した QuantumSystem を表す。

| メソッド | 入力 | 成功時の返り値 | 状態 |
| --- | --- | --- | --- |
| `state(primitive)` | `Self::State` | `LocalStateOf<Q>` | 実装済み。Result を返さない |
| `operator(primitive)` | `Self::Operator` | `LocalOperatorOf<Q>` | 実装済み。Result を返さない |
| `state_at(primitive, &factor)` | 状態 primitive と葉 | `StateOf<Q,I>` | 追加候補。Result を返す |
| `operator_at(primitive, &factor)` | 演算子 primitive と葉 | `OperatorOf<Q,I>` | 追加候補。Result を返す |
| `state_tensor_from_primitives(assignments)` | `(&factor, 状態 primitive)` の列 | `StateOf<Q,I>` | 追加候補。Result を返す |
| `operator_product_from_primitives(assignments)` | `(&factor, 演算子 primitive)` の列 | `OperatorOf<Q,I>` | 追加候補。Result を返す |

```rust
// 未実装の API 案。
let up_down = Spin::state_tensor_from_primitives([
    (&a, Ket::Up),
    (&b, Ket::Down),
])?;

let xx = Spin::operator_product_from_primitives([
    (&a, Op::Sx),
    (&b, Op::Sx),
])?;
```

単独配置は既存の state／operator と at、積構築は各要素を配置して State::tensor／Operator::product_all へ渡す処理に相当する。

| 契約 | state_tensor_from_primitives | operator_product_from_primitives |
| --- | --- | --- |
| 合成・順序 | 入力順の tensor product。状態の因子順も入力順 | 入力順の演算子積。右側が先に作用 |
| 同じ Index の再登場 | Eq によって拒否する | 許可し、積順を保持する |
| 空入力 | 単位空間上の単位状態 | 未束縛の恒等演算子 |
| 全作用空間 | 入力から状態の所属空間を構築 | 未束縛。必要に応じて利用側で bind |

配置先は既存 at と同じ葉に限定し、Tensor は NotAFactor とする。モデル構築が Site の局所空間から全葉を展開する契約とは区別し、暗黙の葉展開・並べ替え・spectator の補完・Index の発行は行わない。operator_product_from_primitives の重複 Index は、同じ因子上の演算子積を記述するために必要である。

primitive は移動して受け取り、配置先のハンドルは共有する。primitive・係数・Index の Clone や数値能力は要求しない。単独配置は現行 at と同様に I: Eq を要求せず、複数式の合成は既存 constructor と同じ I: Eq の範囲で実装する。構築時に primitive の評価や正規化はしない。

配置エラーを入力順に確認し、最初のエラー以降は入力を消費しない。全配置が成功した後に合成の検証を行う。状態の Index 重複と後続入力の配置エラーが共存する場合には、この順序により配置エラーが優先する。

## 4. 構築済みの複合式を集約する入口

通常入力は所有式の列、try 版は同じ式の `Result<_, ExpressionError>` 列を受け取る。返り値はいずれも `Result<式, ExpressionError>` とする。

| 通常入力のメソッド | Result 入力のメソッド | 処理 |
| --- | --- | --- |
| `state_sum(&space, states)` | `try_state_sum(&space, states)` | 同じ空間の状態の和 |
| `operator_sum(operators)` | `try_operator_sum(operators)` | 演算子の和 |
| `state_tensor(states)` | `try_state_tensor(states)` | 状態の入力順 tensor product |
| `operator_product(operators)` | `try_operator_product(operators)` | 演算子の入力順の積 |

state_tensor_from_primitives／operator_product_from_primitives は primitive と配置先の組を受け取り、state_tensor／operator_product は構築済みの式を受け取る。後者では局所重ね合わせを配置した式や、複数因子にまたがる式も合成できる。

空入力・単項・空間互換性・束縛の選択は既存 constructor に合わせる。state_sum の明示空間は、空の和の所属空間を確定するためにも維持する。try 版は、最初の入力 Err で停止し、全入力が Ok の場合に合成を検証する。[try_sum のエラー順序](expression-construction.md#32-入力の集約とエラーの順序)を tensor・積にも適用する。

2026-10-28 に、try 版の入力エラー型は ExpressionError に固定し、最初の入力 Err をそのまま返し、全入力が成功した後に合成を検証する契約を利用者と合意した。入力エラーは合成時の空間不整合などに優先する。短絡は Err より後の要素を Iterator から取り出さない保証であり、先に作った配列・Vec の生成処理の停止や副作用の巻き戻しは保証しない。具体例・任意のエラー型を受け取らない方針・短絡範囲は[try_sum の契約](expression-construction.md#3-try_sum-を追加する)に記録する。これらを実装基準とし、製品 API の検証結果とは区別する。

Operator の和・積では、全項が束縛済みかどうかを見て束縛の整合性を判断する。二項演算の fold へ単純に置き換えてこの契約を変えない。式の構築は遅延し、項の簡約や暗黙の因子順序変更は行わない。

```rust
// 未実装の API 案。構築の入口を Spin に揃える。
let h = Spin::try_operator_sum(
    [Op::Sx, Op::Sy, Op::Sz].into_iter().map(|axis| {
        Spin::operator_product_from_primitives([(&a, axis), (&b, axis)])
    }),
)?.bind(&space)?;

let energy = psi.expectation(&h)?;
```

core の各式型の sum・tensor・product_all、および try 版の共通処理を利用する。StateOf／OperatorOf は型別名であり、独立した実装先ではない。構築処理を QuantumSystem 側へ複製せず、利用者向けの説明と例を QuantumSystem の入口に揃える。

## 5. 零・単位元と局所式の補助

| メソッド | 返り値・契約 |
| --- | --- |
| `state_zero(&space)` | 指定空間上の零状態 |
| `state_unit()` | 単位空間上の単位状態 |
| `operator_zero()` | 未束縛の零演算子 |
| `operator_identity()` | 未束縛の恒等演算子 |
| `local_state_zero()` | 局所零状態 |
| `local_operator_zero()` | 局所零演算子 |
| `local_operator_identity()` | 局所恒等演算子 |

既存の象徴的な zero・unit・identity の constructor に委譲し、係数型の Zero・One や primitive の恒等元生成・評価能力を構築時に要求しない。返り値は各式型であり、Result は返さない。空間を引数に取らない複合式では、Index 型 I を後段の利用や戻り値注釈、明示型引数で指定する必要がある。

局所式の集約を追加する場合は、次を候補とする。

| 通常入力のメソッド | Result 入力のメソッド | 返り値 |
| --- | --- | --- |
| `local_state_sum(states)` | `try_local_state_sum(states)` | 局所状態 |
| `local_operator_sum(operators)` | `try_local_operator_sum(operators)` | 局所演算子 |
| `local_operator_product(operators)` | `try_local_operator_product(operators)` | 局所演算子 |

通常入力版は式を直接返し、try 版だけが Result を返す。局所式では空間を検証する必要がない。入力は構築済みの局所式であり、primitive の列を受け取る入口とは区別する。

## 6. 追加の優先順位と確認事項

1. 計画中の state_sum・operator_sum と try 版を QuantumSystem から呼べるようにする。scaled の変換対応と try_sum の既存計画にもこの入口を反映する。
2. primitive の配置をまとめる state_tensor_from_primitives・operator_product_from_primitives を追加する（名前は推奨案）。
3. 構築済みの式を受け取る state_tensor・operator_product と try 版を追加する（名前は推奨案）。
4. 零・単位元、state_at・operator_at、局所式の集約補助を追加候補として検討する。

初回は scaled の変換対応と上記1〜3までを実装する範囲を2026-10-28に合意した。4は初回の対象外である。各段階で既存の基礎 API に委譲できることを確認し、実装済みと未実装の候補を分けて記録する。

初回に含める QuantumSystem の集約入口は state_sum・try_state_sum・operator_sum・try_operator_sum、配置込み入口は上記2の2メソッド、構築済みの式の tensor・積は state_tensor・try_state_tensor・operator_product・try_operator_product である。これらに必要な core の try_sum・tensor／積の try 版の共通処理も実装する。既存の計画にある4式型の try_sum は基礎処理として含めるが、QuantumSystem の局所式集約用の転送メソッドは4の後続候補のままとする。
新 API の公開名・引数・戻り値は、通常の Rust の型推論とメソッド解決を確認して確定する。

scaled の入力を impl Into<C> に広げる際の推論差は許容し、必要な型注釈や余分な into() の除去で対応することを2026-10-28に合意した。従来の推論との完全互換は初回実装の条件にしない。具体的な移行例と小規模な検証結果は[scaled の型推論の計画](expression-construction.md#23-型推論の変化と移行例)を参照する。

実装時には、独自 QuantumSystem でも既定メソッドが使えること、最小能力、非 Clone の保持要素、空・単項入力、状態の重複拒否と因子順、演算子の重複受理と作用順、葉以外の配置拒否、束縛と未束縛の混在、try 版の短絡・エラー順序を検証する。入力の借用期間を超えて返り値を保持できること、既存 core constructor と同じ型・数学的結果になることも確認する。

rustdoc・examples の基本の案内を QuantumSystem に揃え、型別名には型注釈の用途と基礎式型へのリンクを記載する。現時点では本書の追加候補のコンパイルや回帰テストは行っていない。
