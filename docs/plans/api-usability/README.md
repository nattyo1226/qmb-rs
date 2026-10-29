# API の使いやすさに関する残る検討事項

作成日: 2026-10-27
更新日: 2026-10-30

未解決の使いやすさの論点と、今回のリファクタリングでは実装しなかった補助 API の候補を記録する。以下は追加を決めた機能の一覧ではなく、公開名・シグネチャ・実施時期は未確定である。

現在の構築 API は [ExpressionFactory と型情報の実装・rustdoc](../../../crates/qmb-core/src/system.rs)、[局所式の設計](../../design/primitive-expression.md#7-式構築の入口と係数変換2026-10-28)、[複合式の設計](../../design/composite-expression.md#11-expressionfactory-の構築入口と-fallible-集約2026-10-28)、幾何と配置は [Chain の設計](../../design/geometry.md)を参照する。実装済みの API・移行経緯を未完了の計画として重ねて管理しない。実装済みの Grading は[ordinary / graded の共通評価](../../design/composite-expression.md#12-ordinary--graded-の共通評価と移行)、State／Operator の統合などの延期した機能と未実施の検証は [今後の課題](../../future-work.md)で管理する。

## 1. ExpressionFactory の追加補助

次は初回の実装範囲から外した候補である。既存 constructor と式メソッドで同じ構築ができるため、実際の利用例で繰り返しが目立つ場合に追加を検討する。

| 候補 | 役割 |
| --- | --- |
| `state_at(primitive, &factor)` | 単独の状態 primitive の配置 |
| `operator_at(primitive, &factor)` | 単独の演算子 primitive の配置 |
| `state_zero(&space)` | 指定空間上の零状態 |
| `state_unit()` | 単位空間上の単位状態 |
| `operator_zero()` | 未束縛の零演算子 |
| `operator_identity()` | 未束縛の恒等演算子 |
| `local_state_zero()` | 局所零状態 |
| `local_operator_zero()` | 局所零演算子 |
| `local_operator_identity()` | 局所恒等演算子 |
| `local_state_sum` / `try_local_state_sum` | 構築済みの局所状態の集約 |
| `local_operator_sum` / `try_local_operator_sum` | 構築済みの局所演算子の集約 |
| `local_operator_product` / `try_local_operator_product` | 構築済みの局所演算子の積 |

局所式の集約用メソッドは ExpressionFactory からの入口の候補であり、既にある式型の sum・try_sum・product と区別する。追加する場合も、既存 constructor の空入力・単項・順序・遅延構築の契約を維持し、必要な能力だけを各メソッドに課す。零・単位元の構築に数値演算を要求せず、空間入力のない複合式では Index 型を注釈などで指定する。

## 2. 空間の読み取りと数値出力の型推論

### 2.1. Index から実際の葉を探す補助

`HilbertSpace::find_factor(&index)` 相当は未追加である。現在は leaves() から検索できるため、呼び出しの繰り返しを減らす価値があるか検討する。追加するなら Eq で実際の葉への参照を返し、Index の Clone・発行、因子順序の変更、新しい部分系の構築を行わない。

### 2.2. 状態から MatrixArgs の型を推論する補助

行列 export では、係数型とは別に `MatrixArgs::<SpinHalfState, _>` などの状態 primitive 型の指定が残る。`MatrixArgs::for_state(&psi, &basis)` 相当でこの指定を省く候補がある。

- 状態を型推論だけに使い、作用空間は従来の束縛・明示指定に従わせるか。
- 状態から作用空間も取得するなら、その選択と寿命をどう明示するか。
- 局所状態にも提供するか、状態が存在しない行列出力では何を案内するか。
- 戻り値で選ぶ出力コンテナと、状態型の推論補助をどう説明するか。
- 戻り値型によるコンテナ選択に加え、名前付きの出力補助が必要か。

実装は qmb-export 側で検討し、この補助のために spin-half から export への依存を追加しない。既存の export 契約は [数値 export の設計](../../design/export.md)に従う。

## 3. 作用空間・評価の案内

配置先と spectator を含む全作用空間、bind、明示空間の優先、合成時に未束縛になる条件は、計算結果の解釈に関わる。現在の単一 Operator 型と説明で十分に理解できるか、利用例を通じて確認する。束縛状態を型で区別する案を採用するかは未決である。

明示的な evaluate が必要な場面と、内積・期待値へそのまま渡せる場面、通常の評価と基底展開による集約の違いを、利用者向けの説明で整理する余地がある。これは評価器の変更を決めるものではない。構造表示や計算量・メモリ使用の診断を追加する場合は、[延期した式の管理・性能の課題](../../future-work.md)との範囲を整理する。

## 4. エラーの診断情報

`ExpressionError::MissingIndex`・`IncompatibleSpaces` は対象の Index や期待する順序を保持しない。修正箇所を特定しやすくする情報を、Index に一律の Debug・Display・Clone を課さずに提供できるか検討する。

エラーを `?` で扱う基本経路を保ちつつ、必要な場合に不正な配置・空間の過不足・順序の不一致を診断できることを目指す。既存の try_ 集約のエラー型・優先順位・短絡契約は、診断改善のために暗黙に変更しない。

## 5. 利用開始と独自代数への導線

[examples](../../../examples/README.md) には手作業の式構築、基底展開、Chain モデル、数値 export の実行例がある。一方、次の説明の配置・充実は引き続き検討する。

- 典型的な用途の依存設定・feature・trait import を最初に確認できる案内。trait の import 不足と feature 不足によるエラーを区別できるか。
- crate root の案内と各 rustdoc・設計文書・実行例のつなぎ方。リポジトリ直下の README は未作成。
- テスト内の [3準位系](../../../crates/qmb-core/tests/support/three_level.rs) を、独自の ExpressionFactory の型エイリアス・能力 trait・Basis を定義する利用者向けの最小例へ整理すること。
- 因子・基底・コンテナの順序、メタデータの共有・借用、非正規直交基底での内積の解釈を、数値出力の利用者へどう案内するか。

import の行数削減自体は優先課題とせず、prelude や統合窓口の追加は採用しない。共通の構築入口と、能力ごとに拡張する API の関係を説明する。

## 6. 保留中の構文補助

- **macro**：通常の関連関数を基本 API とする。必要になった段階で、同じ構築処理へ展開する構文補助を検討する。
- **Mul の入力拡張**：任意の `Into<C>` を受け取る blanket 実装は、演算子積・状態への作用の Mul と重複する。標準の複素係数に対応する実数入力の個別実装も含め、拡張は保留とする。変換可能な入力には既存の `scaled(impl Into<C>)` を使う。

## 7. 候補を具体化する際の確認

利用例の変更前後で、記述量だけでなく意図の読み取りやすさとエラー時の修正方法を比較する。Index の同一性・因子順序・作用空間・係数精度を保ち、独自代数・非 Clone 型・no_std の能力を必要以上に制限しない。採用時には公開名と契約を設計文書へ反映し、未解決事項としての記載を整理する。
