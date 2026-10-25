# 係数型の採用方針（一時）

更新日: 2026-10-25

## 1. 位置づけと採用事項

[係数型の依存・能力の調査](refactoring-coefficients.md)と、tensory・ndarray の設計を踏まえた作業用文書である。2026-10-25 に、利用者が未決事項に対する推奨案をすべて採用した。本書は採用事項と実装・検証の基準を記録する。製品コードの移行はまだ行っていない。

**式が係数全体 C を保持し、個別の演算能力で評価する構成**を採用する。core に num-complex・num-traits を含め、Zero・One は num-traits を利用し、付随する能力を受け入れる。不足する能力は個別の trait に分け、標準の Complex へ core 内で直接実装する。

| 項目 | 採用した方針 |
| --- | --- |
| 型引数 | 成分型 T から係数全体 C へ移行する |
| 追加能力 | 共役・実数値取得・実係数の埋め込み・純虚係数の生成を必要な箇所に限定して提供する |
| Norm | std を利用する qmb-normalize の拡張 trait で norm・normalized を提供し、平方根を使う処理の実数型に Float を要求する。provider は導入しない |
| 左スカラー倍 | 標準 Complex は個別実装し、任意の C は scaled(c) を基本経路とする |
| 除算 | 正規化期待値の純実分母を不要に二乗しない、既存の数値的な処理を維持する |

- 式の保持・共通評価は C の内部表現に依存させない。core 内の具体的な trait 実装では Complex の実部・虚部を利用できる。
- Add・Mul 等は標準 trait、不足する能力は小さな独立した trait で表す。
- Norm・正規化は現在の内積から導出し、必要な操作だけに追加の数値能力を要求する。
- 標準の複素数接続は core に置き、平方根を伴う norm・正規化は qmb-normalize に置く。

tensory の構築と実行の分離、Regulator と具体 backend の分離を参考にする。qmb の共通 evaluator と不変共有を維持し、最初から全演算を task・runtime に置き換える要求にはしない。

## 2. 比較した構成

| 案 | 構成 | 利点 | 影響・制約 | 評価 |
| --- | --- | --- | --- | --- |
| A | core が依存のない独自 Complex<T> を所有し、成分型 T の個別能力で計算する | 現行の型引数と成分演算に近く、移行範囲を抑えやすい | core が複素数の二成分表現と商の計算法を所有し続ける。係数表現の差し替えは制限される | 不採用 |
| B | core は数値依存を外し、係数全体 C と個別 trait を使い、具体的な数値処理を外側で接続する | 保持と演算を分離でき、core から数値表現への依存を外せる | 外部型の接続に newtype 等が必要 | 依存解消を目標としていた時点の旧案 |
| C | 零生成・加算・乗算・共役等をすべて演算 provider／context に委譲する | 外部型を直接接続しやすく、同じ型に複数の実行方法を提供できる | evaluator・primitive の規則・基底・モデル構築へ provider の受け渡し方法を設計する必要がある | 今回は不採用 |
| D | core が num-complex に依存し、係数全体 C の個別 trait と、Complex32・Complex64 への直接実装を提供する | 既存の複素数型を直接使え、追加の数値型 wrapper や演算 provider を必須にしない | 係数全体 C への型引数の移行と、操作ごとの能力契約を設計する必要がある | **採用** |

以後は D に従って実装する。数値ライブラリへの依存を持つことと、全係数へ広い数値 bound を要求することは別である。標準型への接続を core に置きつつ、保持・構築・共通評価は必要な能力だけで記述する。

## 3. 係数全体 C を直接保持する

現在の `T` は Complex<T> の成分型である。採用した方針では、その場所を係数全体 C に置き換える。

```rust
pub struct ScaleExpr<E, C> {
    pub coefficient: C,
    pub operand: E,
}

pub trait InnerProduct<C> {
    fn inner_product(&self, rhs: &Self) -> C;
}
```

LinearCombination<P, C>、LocalState<L, C>、LocalOperator<L, C>、State<I, L, C>、Operator<I, L, C> も同じ C を使う。これらの保持型に共通の Scalar・Number・Field 等を必須にしない。LocalStructure の状態・演算子の関連型と次元は、この変更では維持する。

Product・Act・Adjoint・Identity・Basis の係数型も C に揃える。Primitive の規則は LinearCombination<対応する primitive, C> を返し、式木へ map で取り込む現在の接続を保つ。

この変更は型引数の意味を変える。例えば現行の LocalState<SpinHalf, f64> は、新 API では LocalState<SpinHalf, Complex64> 相当へ明示的に移行する。単に T を C へ改名して、f64 を複素係数として扱う旧来の意味を失わせないよう、全利用箇所を確認する。必要なら facade 側に成分型を引数とする便利な型別名を置く。

共通評価の演算結果は同じ C に閉じることを要求する。異なる精度の係数や、異なる C の式を暗黙に混ぜない。明示的な型変換の API は別の判断とする。

## 4. 個別の演算能力

### 4.1. 標準 trait と小さな追加契約

加減乗除・負号には core::ops の trait を用いる。標準 trait の一般性は保ち、共通評価メソッドで必要な場合に `Output = C` を課す。

Zero・One は num-traits を使う。共役は借用から係数を返す能力、実数値の取得は関連型 Real への射影として個別に定義する。具体的な名前の案とシグネチャは以下のとおりである。

```rust
use num_traits::{One, Zero};

pub trait Conjugate {
    type Output;
    fn conjugate(&self) -> Self::Output;
}
pub trait RealPart {
    type Real;
    fn real_part(self) -> Self::Real;
}
```

Zero の Add・ゼロ判定、One の Mul は受け入れる。ゼロ判定は Zero::is_zero を使い、独立した IsZero は追加しない。Conjugate に Neg・RealPart を必須にせず、RealPart も norm squared などで使う時にだけ要求する。Default を零の代用にせず、一般の PartialEq も追加しない。

数学的な整合性は実装の契約として記述する。これらの trait に unsafe は必要としない。Scalar の可換性や共役と乗算の関係が必要な操作については、その法則も契約・検証対象に含める。

Conjugate・RealPart 等の追加 trait の定義元が core なので、同じ crate 内で Complex<f32>・Complex<f64> に直接実装できる。Zero・One と標準演算は既存実装を利用し、再実装しない。成分型までジェネリックな Complex<T> の実装は、必要な条件を確認したうえで追加を検討する。

num-traits の Zero・One は、生成だけの操作でも必要な場所に bound として課す。能力を使わない保持・構築には要求しない。Num・Float の一律な導入は避け、Float は qmb-normalize の平方根を使う実数型に限定する。

### 4.2. 評価経路ごとの条件

| 経路 | 係数側の条件の案 | 保持する特性 |
| --- | --- | --- |
| 保持・map・式の遅延構築・配置・束縛 | 数値 bound なし | 式ハンドルの Clone に C: Clone を追加しない |
| 係数1を明示する constructor | One | One の Mul は受け入れる。Zero を別途要求しない |
| 係数−1による負号の構築 | One + Neg<Output = C> | 初回移行は現在のスカラー倍による表現を保つ |
| LinearCombination のスケーリング | Clone + Mul<Output = C> | 係数全体の乗算が、成分型の Add・Sub を隠蔽する |
| 状態の作用評価 | Clone + Conjugate<Output = C> | 係数の加算・乗算を要求しない現在の評価経路を保つ |
| 項走査・演算子の内部評価 | Clone + Mul<Output = C>、随伴を処理する場合は Conjugate<Output = C> | 暗黙の1を保持する Option も維持できる |
| 局所の内積・trace 等の集計 | Clone + Zero + Mul<Output = C> + Conjugate<Output = C> | Add は Zero が含む。One・Div を追加しない |
| 複合の内積・full trace 等の集計 | 前行と One | 空の tensor product 等に数値の1を使う |
| norm squared | 内積の条件と RealPart | 平方根・除算・ゼロ判定を追加しない |
| 正規化期待値 | 内積・期待値の条件と係数の除算能力 | Zero::is_zero で検査し、自己内積全体で割る。標準 Complex の純実分母の処理を維持し、平方根は要求しない |

実際には primitive 側の Act・Product・Adjoint 等の能力も各経路で要求する。具体的な規則が必要とする分数等の能力は、その実装側の条件であり、上表の共通処理の条件とは分ける。

標準 Complex の Div をそのまま呼ぶだけで既存の数値特性を維持できるとは限らない。正規化期待値の除算は、純実分母を二乗せずに実部・虚部を直接割る経路を残す。係数 C の抽象化との接続方法は実装時に定め、この経路の回帰検証を行う。

Clone は現在の評価器が共有節点の係数を再利用するための条件として残す。tensory の借用加算を参考に、借用演算で扱える箇所は検討するが、全係数に所有値・借用値の全組合せを要求することはしない。最初の移行で評価器全体の Clone 除去を完了条件にしない。

## 5. Norm・正規化の数値方針

### 5.1. qmb-normalize への分離

平方根を伴う norm・状態の正規化を `qmb-normalize` に分ける。新 crate は std を利用でき、no_std 対応は要件としない。norm 用 provider・独自 Sqrt trait・core の `_with` API は導入しない。

| 担当 | 機能・責任 |
| --- | --- |
| qmb-core | 係数の基本契約と標準 Complex への不足する trait の実装、式の構築・共有・評価、作用・内積・期待値、norm squared、平方根を使わない正規化期待値 |
| qmb-normalize | norm・状態の正規化、平方根と実数の逆数の計算、実係数の構築。std の利用を許容する |
| qmb | 新 crate と拡張 trait を再公開し、既定の利用例と feature を接続する |

依存方向は `qmb-normalize -> qmb-core` とし、core は新 crate に依存しない。新 crate は core の公開された inner_product・scaled 等を利用し、非公開節点へのアクセスや evaluator の複製は必要としない。core の既存の num-complex・num-traits 依存は維持する。

現行 core の norm・normalized の inherent method を移行し、新 crate が定義する拡張 trait で提供する。別 crate から core の型へ inherent method は追加できないため、利用側には拡張 trait の import が必要になる。qmb から再公開し、prelude からも利用できるようにする。trait の具体名は実装時に定める。

norm squared は自己内積の実数値への接続、正規化期待値は期待値と自己内積の商なので core に残す。既存の NormalizationError も、正規化期待値と新 crate で共用するため core に残す。

### 5.2. 演算と能力の範囲

自己内積 s に対し、norm は `sqrt(Re(s))`、normalized は `1 / sqrt(Re(s))` を実係数として埋め込み、既存の scaled で状態を返す。正規化で入力の位相を変えない。式を行列へ export する必要はない。

平方根を使うメソッドの実数型 R に num_traits::Float を要求し、係数 C 全体には Float を要求しない。norm は実数値取得と平方根まで、normalized はさらに逆数・実係数の埋め込みを要求する。Complex<f32> の norm は f32、Complex<f64> の norm は f64 とし、暗黙の精度昇格は行わない。

normalized は自己内積全体の厳密なゼロを Zero::is_zero で検査し、零なら既存の ZeroNorm を返す。norm(zero) は零を返す。負値の補正・近似ゼロ判定・特殊値の独自処理は追加しない。自己内積の計算を不要に繰り返さない。

### 5.3. feature と移行範囲

新 crate の標準経路では num-traits の std 対応による Float の平方根実装を利用できるため、no_std 用の libm 接続を必須にする必要はない。既存の qmb・core の libm feature の互換性と再公開方法は、移行時に整理する。Cargo の feature は依存間で統合されるため、新 crate と同時に使う構成では num-traits の std も統合され得る。core 単独で no_std + alloc を使えることと、core の機能・公開 bound が平方根を要求しないことを維持する。

現行の normalized は正規化した状態だけを返す。元の値を「正規化した本体と取り出した係数」の組として保持する regulate、自動規格化、RegulatedState の導入は今回の実装範囲に含めない。

## 6. 外部の数値型と Rust の実装制約

### 6.1. 具体的な複素数の接続

採用した接続方針では、qmb-core が定義した Conjugate 等を qmb-core 自身で num_complex::Complex<f32>・Complex<f64> に直接実装する。trait の定義元による実装なので orphan rule に抵触しない。num-traits の Zero・One は既存の実装を利用する。

共通処理は C の個別能力を要求する。標準の Complex<f32>・Complex<f64> は既存の四則演算を利用できる。Complex<T> の成分型 T までジェネリックにする処理では、num-complex の演算が要求する Num 等を不用意に広げないよう、従来の成分演算の経路も比較する。

標準接続の複素除算では、現行の正規化期待値にある純実分母の経路も引き継ぐ。大きい／小さい実分母を不要に二乗しない処理を、単なる Div への置換で失わせない。具体型へ責任を移す際も数値的な挙動を確認する。

利用者が独自の係数型を所有していれば、その型に直接契約を実装できる。core に実装のない別の外部型を利用者が接続する場合には、newtype 等の対応が必要になる。標準の Complex の接続とは区別する。

### 6.2. 左スカラー倍

`impl<C> Mul<State<..., C>> for C` のような blanket implementation は、標準 trait の orphan rule に抵触する。C に空の Scalar marker を付けるだけでは解決しない。

名前付きの `state.scaled(c)` を任意の C の基本経路とし、所有係数 C を保持する右側の演算子構文は式同士の積との整合性を確認して提供する。標準の Complex<f32>・Complex<f64> には core 内で左スカラー倍を個別に実装する。任意の C 向けの構文用 wrapper は今回導入しない。左右のスカラー倍が同じ意味である現在の数学的契約は保つ。

## 7. 具体的な局所代数との接続

Spin-half の規則では C の0・1・純虚係数・分数を生成する必要がある。一方、一般の式の保持にはそれらを生成する能力は不要である。

既存の規則を C に移す際は、具体的な impl ごとに次を整理する。

- 内積・随伴・恒等元・Sz 基底の係数生成には零・単位元を要求する。
- Trace の係数2は加算を使うか、具体的な定数生成能力を使う。
- Cartesian 基底・Sy を含む規則には、純虚係数の生成方法を用意する。
- 分数は具体的な規則が必要とする能力だけで生成する。全係数へ Float・cast・文字列変換を要求しない。
- XYZ の実数結合は、必要な場所だけ実数から C へ埋め込む。Neel の遅延構築は数値能力を要求しない。

共役・実数値取得・実係数の埋め込み・純虚係数の生成を独立した能力として扱う方針は採用済みである。標準の From 等で表せる能力は既存 trait を利用し、不足する生成能力の trait 名と具体的な bound は移行時に定める。共通処理で C の .re・.im に依存せず、具体代数が要求する能力として接続する。

## 8. 段階的な移行と確認

1. Zero・One は num-traits を使い、係数 C と共役・実数値取得等の契約、Complex<f32>・Complex<f64> への core 内の追加実装を用意する。
2. 薄い線型結合・共通パーツ・演算 trait・局所式・局所評価を C へ移し、spin と独立代数を接続する。
3. 複合式・複合評価・モデル utility・再公開を C へ移す。空間・LocalStructure の構造は維持する。
4. qmb-normalize と拡張 trait を実装し、既存の norm API・feature・例を移行する。
5. core の通常依存・公開契約・テストの依存を確認し、正式な設計文書へ反映する。

各段階は作業順の目安であり、公開 API が不整合な途中状態を完成として扱わない。両精度・非正規直交基底・独立代数・行列照合に加え、次を確認する。

- 数値能力のない C を保持でき、map が C を移動し、式ハンドルの Clone が C の Clone を要求しない。
- num-traits の One で係数1を生成し、Zero・One が不要な構築へ bound を追加しない。ゼロ判定は既存の Zero::is_zero で確認する。
- 加算・乗算を持たない C で、局所規則が対応していれば作用評価できる。
- 正規化期待値に平方根・RealPart を要求しない。零の norm は零を返し、正規化は相殺による zero-norm も拒否する。
- qmb-normalize の拡張 trait を介して norm・normalized を呼べる。core 単独で no_std + alloc を維持し、Float が通常の代数処理へ波及しない。
- 標準の複素数接続で純実分母・複素分母の除算を確認し、現行の正規化期待値の実分母経路による数値特性を維持する。
- 深い式の共有・評価・破棄、最小限の型制約、零係数・項順序・重複の保持を維持する。

## 9. 試作で確認した点と実装上の確認事項

Zero・One の採用を確定する前に、Rust 1.98.0 の独立した一時ファイルで以下をコンパイル・一部実行して確認した。分離した Zero・One・IsZero を使う試作は旧案の確認であり、現在の採用仕様ではない。製品コードは変更していない。

- 数値能力のない型を係数とする構築と、同じ式型の積・右スカラー倍の実装が共存する最小例。
- Add・Mul 等を要求しない One、および外部の数値型を包むローカル newtype への追加 trait の実装。
- no_std の契約側と、外側の平方根を使う norm 方針の分離。正規化係数の計算と、共通処理での零の拒否。
- 外部 trait と外部数値型への接続 crate からの直接 impl が E0117、ジェネリックな裸の左係数への Mul impl が E0210 となること。

この試作は trait と構文の局所的な成立性の確認であり、qmb の全型・借用構文・評価器へ接続した検証ではない。MSRV 1.85 での確認、製品コードのテスト・性能比較も未実施である。

上記の orphan rule の失敗例は、trait と数値型の両方が外部にある別 crate からの接続の確認である。採用した core 内の標準 Complex への不足する trait の実装には当てはまらない。

第1節の5項目に関する設計方針はすべて採用済みである。追加 trait の具体名・シグネチャ、純実分母の除算経路の接続、既存 libm feature と呼び出し箇所の移行は、採用方針に沿って実装・検証する事項として残る。LocalStructure の廃止や State／Operator の統合、対数係数・自動規格化の導入は、本方針の前提にも実装目標にも追加しない。
