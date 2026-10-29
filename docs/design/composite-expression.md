# 埋め込み・複合空間の式と数学的能力

更新日: 2026-10-30

## 1. 対象と責任分担

本書は、[要件](../requirements.md)と[アーキテクチャ概要](../architecture.md)に基づき、局所式の配置、全 Hilbert 空間への束縛、および複合系の状態・演算子式の数学的能力と整合性を定める。[局所 primitive・式木の数学的能力](primitive-expression.md)を、[Hilbert 空間の再帰構造](hilbert-space.md)に沿って複合系へ拡張する契約である。

2026-10-27に、partial trace を一旦 future work とし、Step 5 の担当範囲と初期版の必須機能・検証項目から外した。本書の partial trace に関する契約と検証条件は将来の実装向けに保持する。

初期版で合成する全基本因子は、同じ primitive 型 S・O で表す一種類の局所代数を持つ。配置・束縛・作用・縮約も同じ局所代数の系の内部で行う。異種局所構造の混在・合成と、そのための型消去・保持機構は後続へ延期する。異なる Index や入れ子の grouping は、異種局所構造を意味しない。

| 対象 | 責任 |
| --- | --- |
| 局所構造 | 状態・演算子 primitive と局所的な数学的能力を関連付け、局所契約の整合性を担う |
| 複合 Hilbert 空間 | 子の構造・基本因子の Index と順序・ordinary / graded の種類から数学的能力を合成し、複合系の整合性を担う |
| 状態・演算子式 | 対応する空間の規則に従って構築・評価される |
| 共通の合成・埋め込み・評価機構 | 局所契約を保存する共通の合成則を提供し、複合 Hilbert 空間の数学的能力を実現する |

複合 Hilbert 空間は、局所構造に相当する整合性の責任を持つ。これは概念上の責任分担であり、独立した「複合構造」型、追加の空間の木、または個々の演算を空間のメソッドとすることを要求しない。複合式の実装方針は第9節、現在の公開 API は第10節にまとめる。

局所構造の提供者が局所契約を満たすことを前提に、共通機構はその契約を保存する合成則を提供する。子が基本因子でも複合空間でも、同じ数学的契約に従って再帰的に合成する。複合系ごとに積・作用・随伴などを独立に再定義することは要求しない。

## 2. 配置・束縛と式の構築

局所演算子式と Index の組を、複合演算子式の葉とする。状態側も、局所状態式と Index の組を葉とし、tensor product と有限線型結合によって複合状態式を構成する。配置された葉は元の局所構造を保持する。葉が保持する局所式は、単一 primitive に限らない。

Index の付与を「配置」、spectator を含む全 Hilbert 空間の確定を「束縛」と呼び分ける。配置時に全作用空間を要求せず、Index は新しく発行しない。局所式全体に同じ Index を付与する。

| 対象 | 構築・提供する操作 | 結果・契約 |
| --- | --- | --- |
| 状態・演算子式 | 局所式への Index の付与 | 指定因子に属する、または作用する葉を作る |
| 状態・演算子式 | 加算、scalar multiplication、zero | 有限線型結合を構築する。状態の加算は互換な所属空間上で行う |
| 状態式 | Tensor product | 独立した因子上の状態を指定順序で合成する。基本因子の重複と合成条件の不整合を拒否する |
| 演算子式 | 積 | 同じ因子・異なる因子への演算子を合成し、入力順序を保持する |
| 演算子式 | 全 Hilbert 空間への束縛 | 参照因子・局所構造・合成条件を検証し、全作用空間を確定する |
| 状態式 | 所属空間の確定・確認 | 各因子の状態と tensor 構造から所属空間を定める |
| 状態・演算子式 | 空間情報の保持 | 空間確定後は、評価・基底展開で zero・identity となっても所属空間・始域・終域を保持する |

演算子用の独立した tensor product 操作・節点は設けず、異なる因子上の演算子にも通常の積を用いる。本書の演算子の tensor 記法は、合成則を数学的に記述するためのものである。積の構築だけから全作用空間を自動合成しない。

状態の配置は、その Index の因子に属する状態を作る。新しい因子を状態に追加するには、その因子の状態または明示的な埋め込み写像が必要であり、状態を identity spectator で補わない。

明示空間・既存の束縛先・状態空間の選択、束縛済み operand の合成、およびエラーの検出時点は、要件の 4.3 節・5.1 節・5.7 節に従う。特に、明示空間があれば旧束縛先より優先し、元の式を変更せず指定先の整合性を検証する。

## 3. 数学的能力と導出する操作

| 対象 | 能力 | 入力 | 結果・評価時の担当 |
| --- | --- | --- | --- |
| 状態式 | 内積 | 互換な空間上の二状態 | 複素 scalar。半双線型性と tensor product の内積則を用い、局所式の内積へ委譲する |
| 演算子式 | Product | 二つの演算子式 | 演算子式。同一因子では局所積、異なる因子では空間の合成則を適用する |
| 状態式 | Act | 状態式と演算子式 | 状態式。線型性・因子順序・grading に従い、状態 primitive の Act による局所作用へ委譲する |
| 演算子式 | 随伴 | 演算子式 | 演算子式。係数の共役・積の順序反転・局所式の随伴を扱う |
| 演算子式 | Full trace | 全作用空間が確定した演算子式 | 複素 scalar。通常の非正規化 trace を計算する |
| 演算子式（future work） | Partial trace | 演算子式・全作用空間・縮約対象の部分系 | 残存空間に束縛された演算子式。残存因子の Index・局所構造・順序を保持する |
| 状態・演算子式 | Grading の処理〈graded の場合〉 | homogeneous な成分 | 成分の parity から符号を処理する。混合 parity の式に単一の parity を要求しない |

積と作用は各引数について線型、内積は第一引数について共役線型・第二引数について線型、随伴は共役線型、full/partial trace は線型とする。`AB` は右側の `B` から状態に作用する。局所式に含まれる primitive 間の線型関係を保存し、結果の係数列や式構造の一意性は要求しない。

| 対象 | 導出する操作 | 導出に用いる能力 |
| --- | --- | --- |
| 状態式 | Norm・正規化 | 内積による `sqrt(⟨ψ,ψ⟩)` と scalar multiplication。zero-norm state の正規化は拒否する |
| 状態式（後続） | Bra | Hilbert 空間の内積により ket の随伴として得る。ket から bra への対応は共役線型とする |
| 演算子式 | Commutator / anticommutator | 積と加算による `AB−BA` / `AB+BA` |
| 演算子式 | Hilbert–Schmidt 内積 | 非正規化の `Tr(A†B)` |
| 状態・演算子式 | 期待値・正規化期待値 | 作用と内積による `⟨ψ,Aψ⟩`、および `⟨ψ,ψ⟩` による除算。後者は zero-norm state を拒否する |

共通の複合式・配置・束縛・評価は `qmb-core` が担当し、具体的代数の局所規則を合成する。source の自動保持・式の構造比較・表示は、要件の 4.5 節に従って後続へ延期する。局所式と Index・局所構造・作用空間の整合性判定は維持し、数学的な等価性に式構造の一致を要求しない。

指定基底への展開・座標取得は、局所規則による評価とは別の操作として core が提供する。[局所基底の能力](basis.md)と複合基底の順序規約から処理し、同じ基底要素の係数をまとめる。通常の積・作用・内積・trace の評価に基底展開を必須としない。複合演算子の基底展開には、要件の 5.1 節に従って全作用空間を確定する。

複合基底は基本因子ごとの局所基底の tensor product に限定し、空間の因子順序と合成則に従う。Grouping によって基底を変えない。展開結果は通常の式木とし、基底識別情報・展開済み状態を保持させず、Index・局所構造・確定済みの空間を保持する。別基底への変更は通常の再展開で行い、数値係数の解釈に必要な基底情報は export 側の契約で扱う。

構築時に式全体を展開・評価しない。必要な結果を要求した時点で必要な部分を評価し、入力は変更せず不変な部分構造を共有する。作用・内積・trace に完全な数値行列を要求せず、共通の完全展開形や評価結果のキャッシュも要求しない。

分配・局所積などの特定の簡約ルールのみを選んで適用する独立した操作は future work とする。これらを基底展開・作用・内積・trace などの内部処理として用いることは維持する。

## 4. 配置・埋め込みが保存する関係

局所演算子式 `A` と Index `i` を、全空間 `H` の対応する因子へ埋め込んだ演算子を `E_i^H(A)` と書く。この埋め込みは、Index の存在、局所構造の一致、および合成条件の整合性を前提とする。同じ局所構造の式 `A,B` と複素 scalar `a,b` に対して、以下を要求する。

```math
\begin{aligned}
E_i^H(aA+bB)&=aE_i^H(A)+bE_i^H(B),\\
E_i^H(AB)&=E_i^H(A)E_i^H(B),\\
E_i^H(A^\dagger)&=E_i^H(A)^\dagger,\\
E_i^H(I_i)&=I_H.
\end{aligned}
```

局所式を計算してから埋め込む経路と、埋め込んでから計算する経路は同じ演算子を表す。同じ Index 上の積には局所代数を用いる。異なる Index `i,j` 上の演算子は ordinary では可換、graded では homogeneous な `A,B` に対して次の関係を満たす。

```math
E_i^H(A)E_j^H(B)
=(-1)^{p(A)p(B)}E_j^H(B)E_i^H(A),\qquad i\ne j.
```

Ordinary な埋め込みは他の因子を identity とする。Graded な埋め込みは、次節の作用則に従い因子順序に応じた parity を含む。参照 Index の付与だけから、その演算子が通常の tensor 表現でも他因子上で identity であるとは判断しない。

## 5. Ordinary / graded の合成則

以下の二因子の規則を、有限線型結合と多因子へ拡張する。`A,C` は第一因子、`B,D` は第二因子の演算子、`ψ,ψ'` は第一因子、`φ,φ'` は第二因子の状態とする。

Ordinary では以下に従う。

```math
\begin{aligned}
(A\otimes B)(C\otimes D)&=AC\otimes BD,\\
(A\otimes B)(\psi\otimes\phi)&=A\psi\otimes B\phi,\\
(A\otimes B)^\dagger&=A^\dagger\otimes B^\dagger,\\
I_{H\otimes K}&=I_H\otimes I_K.
\end{aligned}
```

Graded では、符号に用いる成分を homogeneous とし、`p` を even `0` / odd `1` とする。

```math
\begin{aligned}
(A\otimes_g B)(C\otimes_g D)
&=(-1)^{p(B)p(C)}AC\otimes_g BD,\\
(A\otimes_g B)(\psi\otimes_g\phi)
&=(-1)^{p(B)p(\psi)}A\psi\otimes_g B\phi,\\
(A\otimes_g B)^\dagger
&=(-1)^{p(A)p(B)}A^\dagger\otimes_g B^\dagger,\\
I_{H\otimes_g K}&=I_H\otimes_g I_K.
\end{aligned}
```

Graded でも状態の内積は通常の正定値な Hilbert 空間の内積とする。以下では `⊗_*` を ordinary / graded のいずれかとし、どちらでも同じ規則を用いる。

```math
\langle\psi\otimes_*\phi,\psi'\otimes_*\phi'\rangle
=\langle\psi,\psi'\rangle\langle\phi,\phi'\rangle.
```

内積は第一引数について共役線型とし、graded であることだけを理由に符号を追加しない。Graded な空間の even・odd 部分は直交させる。Homogeneous な状態・演算子の tensor 積の parity は各成分の parity の和 modulo 2、恒等演算子は even、随伴は parity を保存する。混合 parity の式は成分ごとに処理する。Zero の parity の選択は結果に影響させない。

通常の演算子積では、graded でも `(XY)† = Y†X†` とする。上記の graded な tensor の随伴の符号は、積の順序反転後に因子順序を元に戻す際の交換則と整合する。

これらの graded な積・作用・随伴は、[Crismale・Rossi・Zurlo, On C*-norms on Z₂-graded tensor products, 2–3 節](https://link.springer.com/article/10.1007/s43037-021-00168-0)の規約に対応する。qmb-rs ではこの規約を、通常の Hilbert 空間上の作用・随伴と整合する合成則として採用する。

通常の tensor 表現へ変換する場合は、要件の 5.5 節に従って parity 演算子を挿入する。Intrinsic な graded 規則と parity string による経路は同じ作用を表し、変換済みの符号を再度適用しない。

多因子の homogeneous な成分については、同じ規則を次の形で用いる。

```math
\begin{aligned}
(\mathop{\otimes_g}_i A_i)(\mathop{\otimes_g}_i B_i)
&=(-1)^{\sum_{i>j}p(A_i)p(B_j)}\mathop{\otimes_g}_i(A_iB_i),\\
(\mathop{\otimes_g}_i A_i)(\mathop{\otimes_g}_i\psi_i)
&=(-1)^{\sum_{i>j}p(A_i)p(\psi_j)}\mathop{\otimes_g}_i(A_i\psi_i),\\
(\mathop{\otimes_g}_i A_i)^\dagger
&=(-1)^{\sum_{i<j}p(A_i)p(A_j)}\mathop{\otimes_g}_i A_i^\dagger.
\end{aligned}
```

## 6. 入れ子・空・単因子との整合性

子が複合空間の場合も、その子の数学的能力と parity を用いて同じ合成則を再帰的に適用する。同じ因子順序・同じ積の種類について、以下の括り方の対応で結果を変えない。

```math
(\psi_1\otimes_*\psi_2)\otimes_*\psi_3
\longleftrightarrow
\psi_1\otimes_*(\psi_2\otimes_*\psi_3).
```

この対応は因子の交換を含まないので、graded でも符号を追加しない。積・作用・内積・随伴・trace は、この対応のもとで一致する。数学的な整合性を要求するものであり、利用者向けの Regrouping 操作、permutation、または構築時の平坦化は追加しない。

- 空の Tensor は一次元の単位空間とし、その状態 `1` と恒等演算子は even、恒等演算子の trace は `1` とする。
- 単因子 Tensor の能力は、その因子の能力に一致する。
- Index を持つ一次元 Factor は空の Tensor と区別し、その Index を保持する。
- 空・単因子の同一視後も ordinary / graded の混在が残る構成は、既存の初期版の制限に従って拒否する。

## 7. Full trace・partial trace と埋め込み

Full trace は、ordinary / graded ともに通常の非正規化 trace とする。全作用空間 `H` の任意の完全正規直交基底 `{e_n}` に対して次で定義し、基底選択に依存しない。

```math
\operatorname{Tr}_H(X)=\sum_n\langle e_n,Xe_n\rangle.
```

`Tr(XY) = Tr(YX)` と `Tr(I_H) = dim H` を要求する。巡回性に parity の符号を追加せず、supertrace は採用しない。

以下の partial trace は future work の数学的契約であり、今回の実装・検証対象には含めない。

Partial trace は、確定した因子順序と通常の tensor 表現に対して定義する。`ambient = H` の部分系 `T` を縮約し、残存空間を `R` とする。残存因子の基底 index を `r,r'`、縮約因子の完全正規直交基底 index を `t` としたとき、`B = Tr_T^H(X)` の行列要素は以下とする。

```math
B_{r,r'}=\sum_t X_{r,t;r',t}.
```

複合 index は元の因子順序に配置する。この式は意味の定義であり、計算時に全基底の列挙や行列生成を要求するものではない。結果は `R→R` の演算子式とし、残存因子の Index・局所構造・相対順序を保存する。次を満たす。

- 線型性と随伴を保存する: `Tr_T^H(X†) = Tr_T^H(X)†`。
- 残存空間で full trace を取ると元の full trace に一致する: `Tr_R(Tr_T^H(X)) = Tr_H(X)`。
- 互いに重ならない因子集合の逐次縮約は、まとめた縮約と一致する。それぞれ残存空間の通常の作用に対して縮約する。
- 通常の tensor 表現で恒等的に作用する縮約因子は、その次元を係数へ掛ける。
- 空の対象の縮約は元の演算子と等価とし、全因子の縮約は 0-site 空間の scalar operator expression とする。
- Zero・scalar identity の結果でも残存空間を保持する。

`ambient` の選択と `trace_out` の部分系判定は要件の 5.4 節に従う。残存空間を結果として保持することは、Hilbert 空間から部分系を取得する利用者向け操作を追加するものではない。縮約前の空間と対象を保持した遅延節点は、結果の束縛先を変更しても元の縮約の意味を保持する。

Graded な場合も埋め込み後の実際の作用を縮約する。二 mode の例では、局所消滅演算子 `a_1` の配置から得られる通常の作用は `c_1 = P_0 ⊗ a_1` であり、

```math
\operatorname{Tr}_0(c_1)=\operatorname{Tr}(P_0)a_1=0.
```

結果は残存空間 `H_1→H_1` の zero operator である。構文上で参照されない mode 0 を identity spectator と判断したり、縮約対象の Index を消して元の式を再埋め込みしたりしない。

## 8. 整合性の保証と検証対象

局所契約を満たす子から、複合 Hilbert 空間の合成則によって以下を保存する。

- Product の結合則、`(XY)ψ = X(Yψ)`、恒等演算子の単位元則
- 内積の共役対称性・正定値性、`⟨φ,Xψ⟩ = ⟨X†φ,ψ⟩`
- `(XY)† = Y†X†` と `(X†)† = X`
- 線型性・半双線型性・共役線型性、および primitive 間の線型関係
- 通常の trace、partial trace、基底座標による計算との整合性
- 同種の積の括り方によらない数学的結果と、graded / 通常 tensor 表現の間での作用の一致

任意のユーザ定義規則について、結合則などを実行時に自動証明することは要求しない。操作時には、Index の存在、局所構造の一致、状態・作用空間の互換性、因子の重複、scalar の精度・grading・積の種類、および縮約対象の部分系条件を検証する。必要な情報が構築時にあれば早期に検証し、空間依存の条件は束縛・評価時に検証する。Index と局所構造の対応を全体で登録する機構は設けない。

検証では、積を先に計算してから作用する経路と順次作用する経路、階層的な合成と基本因子に沿う計算、代数計算と独立に構成した小規模行列による参照計算を比較する。随伴・full/partial trace・基底座標についても確認する。Graded な場合は odd・even・混合 parity の和、離れた因子、括り方、parity string、および上記の縮約例を含める。

等式は数学的な意味について要求し、浮動小数点計算のビット単位の一致とは区別する。数値精度・特殊値・許容誤差は要件の 4.1 節・5.6 節に従う。

## 9. 未決事項と実装範囲

Step 3 の設計判断は以下を採用済みとする。Ordinary / graded な配置・束縛・複合評価と基底処理は実装済みであり、後続の partial trace とは区別する。段階3の検証では、spin-half の三因子系と三次元代数の二因子系を両精度の独立した行列計算と照合し、空間・束縛・エラー・共有、および深い式の処理を確認した。

| 対象 | 採用した方針 |
| --- | --- |
| 複合状態 `State` | 常に確定した所属空間を持つ。Zero と空和には空間を明示し、空の tensor product は0因子空間上の状態 `1` とする |
| 複合演算子 `Operator` | 式本体と任意の束縛先を持つ一つの型とし、束縛は外側のハンドルに集約する。合成時には全 operand が束縛済みの場合だけ空間を比較し、子節点に旧束縛を残さない |
| 空間の選択 | 明示空間を最優先し、次に既存の束縛先を使う。未束縛の作用は状態空間を使い、空間未指定の未束縛 trace はエラーとする |
| 共有 | 式本体と空間を不変な `Arc` ハンドルで共有し、節点 enum は非公開とする。状態・演算子は所属空間または確定した始域・終域を参照できる |
| エラー | 名前付きメソッドを基本とし、検証が必要な構築・接続・評価は `Result` を返す。検証を伴う演算子構文も `Output = Result<式, Error>` とし、空間不整合を panic にしない |
| 遅延評価 | 空間に関する検証は情報が揃った時点で行い、局所規則の実行を遅延させる。Zero・identity・零係数でも不正な参照は見逃さない |
| 評価経路 | 作用・内積・full trace を分担し、全経路共通の完全展開形を強制しない。作用は primitive の Product を要求せず、full trace は必要な項の局所積・局所 trace へ委譲する |
| 深い式 | 評価・走査・最後の所有者になった節点の破棄には明示的な作業スタックを使う |

Step 3 は期待値、commutator・anticommutator、Hilbert–Schmidt 内積までを含む。局所式・複合式の Norm・正規化・正規化期待値は qmb-normalize に実装済みであり、[正規化 API の契約](primitive-expression.md#41-norm正規化の-api-と-feature)に従う。bra 型・構築 API は初期版では導入せず、braket は内積で計算する。基底展開・座標取得は Step 5 で実装済みであり、初期版の必須範囲を維持する。Partial trace は future work とし、Step 5 と初期版の必須範囲に含めない。

Step 3 で非公開の節点・作業スタック、`ExpressionError` と必要なシグネチャを実装した。LocalStructure 廃止後の公開型は `State<I, S, O, C>`・`Operator<I, O, C>` とする。現在の API は第10節に記録する。Partial trace の節点と数値 API など、後続段階の詳細はその段階で設計する。

本書は数学的な契約を定める。[初期版の提供範囲](../architecture.md#9-初期版の提供範囲と完成条件)は同種局所構造の ordinary な演算とした。現在は第12節の ordinary / graded 共通評価まで実装済みであり、異種局所構造の合成は引き続き[後続へ延期](../future-work.md#2-延期した機能と検討事項)する。独立した複合構造型や、複合系専用の primitive・局所規則の追加実装も要求しない。

## 10. 現在の公開 API

| 対象 | 主な API と契約 |
| --- | --- |
| 局所式の配置 | `LocalState::at(&factor)`・`LocalOperator::at(&factor)` は因子への配置を検証し、局所式を共有する。状態は所属空間を持ち、演算子は未束縛とする |
| `State<I, S, O, C>` の構築・空間 | `place`・`zero(&space)`・`unit()`・`sum`・`try_sum`・`tensor_in(&space, states)`・`try_tensor_in(&space, states)`・`scaled` と、`space`・`domain`・`codomain` を提供する |
| `State<I, S, O, C>` の演算・評価 | `act`・`act_in`・`evaluate`・`inner_product`・`norm_squared`・`expectation`・`expectation_in` を提供する。作用は遅延構築する |
| `Operator<I, O, C>` の構築・空間 | `place`・`zero()`・`identity()`・`scalar`・`sum`・`try_sum`・`product`・`try_product`・`scaled`・`adjoint`・`bind` と、`space`・`domain`・`codomain` を提供する。未束縛の場合、空間の取得は `None` を返す |
| `Operator<I, O, C>` の演算・評価 | `commutator`・`anticommutator`・`trace`・`trace_in`・`hs_inner_product`・`hs_inner_product_in` を提供する。交換子は遅延構築する |
| エラー | 空間の重複には `SpaceError`、複合式の接続・評価には `ExpressionError`、次元取得には `DimensionOverflow` を使う |

通常の和・積の名前付き構築は列を受け取る `sum` / `product` とし、二項には標準の `+` / `*` を使う。Index の照合と数値能力は必要な操作で要求する。検証が必要な複合操作は `Result` を返し、状態への作用評価は演算子 primitive の `Product` を要求しない。基底展開・座標取得は[基底の設計](basis.md)、norm・正規化・正規化期待値は[正規化 API の契約](primitive-expression.md#41-norm正規化の-api-と-feature)に従う。Hilbert 空間の API は[空間の設計](hilbert-space.md#10-現在の公開-api)に記録する。

## 11. ExpressionFactory の構築入口と fallible 集約（2026-10-28）

QuantumSystem は状態・演算子 primitive と係数型だけを関連付ける。構築は ExpressionFactory の以下の inherent method から既存の core 式を生成する。空間や Index を登録・保持せず、係数・primitive・Index の Clone や評価能力を一律に要求しない。複数式の合成は既存と同じ I: Eq を要求する。

`ExpressionFactory<S,O,C>` は値を保持せず、インスタンス生成を必要としない構築入口であり、QuantumSystem を一括実装する。`SpinHalf<C>`・`SpinOne<C>`・`Fermion<C>` は具体的な primitive を選ぶ型エイリアスで、構築に trait の import は不要である。独自の primitive も ExpressionFactory の型エイリアスから利用できる。同じ primitive と係数型からは同じ式型を構築し、システム名を式の同一性へ追加しない。

`LocalStateOf<Q>`・`LocalOperatorOf<Q>`・`StateOf<Q,I>`・`OperatorOf<Q,I>` は型注釈として維持する。システム単位のジェネリック API では `Q: QuantumSystem` から型情報を取得し、構築には `LocalStateOf::<Q>::primitive(...)` などの式型の入口、または `ExpressionFactory<Q::State,Q::Operator,Q::Coeff>` を使う。QuantumSystem 自体に構築メソッドは置かない。モデル固有の構築は各 primitive crate の公開 `models` module に分離し、[幾何とモデルの設計](geometry.md)に従う。

| 入口 | 入力・契約 |
| --- | --- |
| product_state | 明示空間と状態 primitive の列。空間の葉の順に一つずつ配置し、個数の完全一致を検証する |
| operator_product_at | 葉と演算子 primitive の組の列。同一 Index の再登場を許可し、入力順に積を構築。右端が先に作用。空入力は未束縛 identity |
| state_sum / try_state_sum | 明示空間と状態列。空入力はその空間の zero |
| operator_sum / try_operator_sum | 演算子列。空入力は未束縛 zero |
| tensor_in / try_tensor_in | 明示空間と配置済み状態列。葉の連結が指定空間を順序も含めて完全に覆い、部分系の種類も互換であることを検証する |
| operator_product / try_operator_product | 配置済み演算子列の入力順の積。空入力は未束縛 identity |

operator_product_at の _at は配置も行う入口を意味し、裸の primitive の列ではなく (&factor, primitive) の IntoIterator を受け取る。配置先は葉のみで、Tensor の自動展開はしない。演算子は明示的に bind するまで未束縛であり、spectator を含む全作用空間を積の配置先から自動構築しない。入力スライスや配置先の変数の生存期間を越えて、共有した式・空間を保持できる。

try_ は Result 入力を表し、通常版も複合式では Result を返す。エラー型は ExpressionError に固定する。最初の入力 Err をそのまま返し、その後は Iterator の next を呼ばない。全入力が Ok の場合だけ合成を検証するため、入力エラーは先の項の空間不整合・重複より優先する。事前に配列・Vec を作った場合の生成処理や副作用は巻き戻さない。全入力成功後の合成検証が失敗する場合も、入力生成は最後まで行う。

State::try_sum・try_tensor_in、Operator::try_sum・try_product を基礎 API として公開し、ExpressionFactory はこれらへ委譲する。演算子の束縛は全入力を見て決め、後続に未束縛項があれば、それ以前の束縛が互換でなくても未束縛として構築する。逐次二項 fold で早期に拒否する実装にはしない。

各入口は構築のみを行い、評価・簡約・基底展開・自動正規化・暗黙の因子順序変更は行わない。構築の共有処理を使い、try 版のためだけの再収集や式の複製を避ける。単独 primitive 配置、zero／identity の ExpressionFactory 補助、局所式の集約用転送メソッドは後続候補として維持する。

## 12. Ordinary / graded の共通評価と移行

空間の `TensorKind` と有効な種類は[空間の API](hilbert-space.md#10-現在の公開-api)、primitive の parity / P は[局所能力](local-structure.md#8-grading-の局所能力と移行)に従う。`HilbertSpace<I>`・`State<I,S,O,C>`・`Operator<I,O,C>` の型引数は維持する。幾何・Index に grading の所有や符号処理を追加しない。

### 12.1. 構築入口

| 対象 | ordinary / graded 共通の入口 |
| --- | --- |
| State | `tensor_in(&space, states)` |
| fallible な State | `try_tensor_in(&space, states)` |
| ExpressionFactory | `tensor_in(&space, states)` |
| fallible な ExpressionFactory | `try_tensor_in(&space, states)` |
| primitive の配置・合成 | `product_state(&space, primitives)` |

所属空間を先に構築し、状態構築ではその空間を共有する。状態の Tensor 節点は子の列だけを保持し、種類を重複して保存しない。状態側の種類指定版・`_g` 版・二項専用版は撤去した。HilbertSpace 側の種類指定 API は維持する。

入力状態の葉の連結は指定空間の葉と順序も含めて完全一致し、複数因子の部分系の有効な種類も互換でなければならない。grouping は異なっていてよい。不足・過剰・重複・順序違い・種類不整合、primitive 数の不一致は `ExpressionError::IncompatibleSpaces` とする。第11節の入力エラー優先・消費停止を try_ 入口で維持する。

空入力は零因子の単位空間だけで受理し、指定した空間のハンドルと grouping を保持する。単項も検証後に指定空間を保持する。`State::unit()` は正準的な単位空間を作る便宜として維持し、`State::zero(&space)` は指定空間上の零状態を作る。構築に Index・primitive の Clone や係数演算を追加要求しない。

### 12.2. 共通 API と最小能力

`evaluate`・`trace`・`expand`・`coordinates` と各明示空間版を ordinary / graded 共通の入口とし、種類別の評価 API や adapter を追加しない。Rust の trait bound は実行時の kind で切り替えられないため、ordinary 空間でも操作に必要な grading 能力を要求する。実際の符号処理は選択空間の `effective_kind()` で決め、ordinary では P の挿入や交差符号を加えない。

| 操作 | 従来の能力に追加するもの |
| --- | --- |
| 複合状態の evaluate / 内積 / norm squared / 期待値 / 展開 / 座標と各明示空間版 | O の `HomogeneousParity + ParityOperator` |
| 複合演算子の trace / Hilbert–Schmidt 内積と明示空間版 | O の `HomogeneousParity + ParityOperator` |
| 複合演算子の展開 / 座標と明示空間版 | O の `HomogeneousParity`、C の `Neg<Output=C>`。Intrinsic な整列に P は不要 |
| export / norm / 正規化 / 正規化期待値 | 呼び出す core の共通評価経路の能力を引き継ぐ |

構築・配置・束縛・和・積・遅延した作用や随伴には評価能力を一律に要求せず、局所だけの API の bound も維持する。状態への作用は O の Product を要求しない。P を局所 Act で作用させて既存の Scale 節点を使うため、交差符号のためだけに C の One / Neg を追加要求せず、S の parity を直接問い合わせる必要もない。

### 12.3. 遅延した Act の scope

複合 Act は構築時に選択した HilbertSpace を保存する。未束縛の演算子構文は、明示空間、既存束縛、作用先状態の順で全空間を選ぶ。明示空間は旧束縛に優先する。束縛済み演算子同士の合成では、有効な種類の不一致も空間不整合となる。

Odd な B と odd な状態 `s_a` に対し、次の二つを区別する。

```text
tensor_in(&H_ab, [s_a, s_b.act(B_b)])      // B の作用は H_b 内で構築済み
tensor_in(&H_ab, [s_a, s_b]).act(B_b)      // B は H_a ⊗g H_b 上で作用
```

後者だけが `s_a` の parity による負号を持つ。評価タスクは Act ごとの scope に切り替わり、先に単因子上で構築した作用を外側の tensor に組み込んでも scope を拡大しない。先に評価してから tensor 合成する経路と遅延したまま合成する経路は一致し、展開・再構築でも意味を保存する。

未評価の演算子構文は選択した全空間に対して解釈する。一方、作用評価済み・基底展開済みの結果を再束縛しても、元の未評価式が別空間で持つ作用を再現するとは限らない。自動的な再評価・source 保持は追加しない。

### 12.4. Intrinsic な項と通常表現

| 内部表現 | 内容・用途 |
| --- | --- |
| `IntrinsicOpTerms` | 演算子式を homogeneous な成分へ分け、選択空間の葉順に整列した項。graded 演算子基底の展開・座標に用いる |
| `OrdinaryOpTerms` | 配置ごとに parity string を挿入した通常表現の局所因子列。full trace に用いる |

Intrinsic な整列は異なる因子の odd 成分同士が交差するたびに負号を一度加え、同じ因子の局所積は元の順序を保つ。Index の値の大小では整列しない。積の随伴による順序反転後にも同じ整列を使い、tensor 随伴の符号を別に重ねない。

Homogeneous な intrinsic 項の通常表現は次となる。

```math
\rho(A_0\otimes_g\cdots\otimes_g A_{n-1})
=\bigotimes_{i=0}^{n-1}\left(A_i P_i^{q_i}\right),
\qquad q_i=\sum_{j>i}p(A_j)\pmod 2.
```

P は A の右側に置く。例えば `A ⊗g B` は `A P^p(B) ⊗ B` となる。通常表現の因子列を再び intrinsic な graded 項として解釈せず、P や符号を二重に適用しない。これらは一時的な private 表現であり、公開の変換 API・全式共通の完全展開形・永続キャッシュではない。

状態への作用は、選択した Act の scope 内で odd primitive の先行因子へ P を挿入する。混合 parity・もつれた状態にも P の実際の局所作用を使い、単一の定数符号に置き換えない。局所式の和・積・随伴も primitive ごとに処理する。状態全体の積和展開・全基底・密行列を共通の前処理にせず、作用しない枝を共有する。

状態内積には追加符号がなく、full trace は通常表現で局所積・Trace を縮約する。P が入る可能性のある因子を下ろす前に identity spectator と判定しない。二 mode の `c_1 = P_0 ⊗ a_1` に対し、`c_1†c_1 = I_0 ⊗ n_1` の trace も検証する。通常の trace の巡回性に parity 符号を加えず、S・状態基底・usize で表せる総次元も要求しない。Hilbert–Schmidt 内積は同じ trace 経路を使う。Supertrace / partial trace は追加しない。

空間検証は Zero・零係数・identity でも先に行い、不正な Index や空間を見逃さない。利用例は [tensor_products.rs](../../examples/src/bin/tensor_products.rs)。この共通機構を使う製品実装は [Fermion](local-structure.md#9-fermion-の局所代数)・[Spin-one](local-structure.md#10-spin-one-の局所代数)に記録する。異種 primitive、実質的な ordinary / graded 混在、permutation、regrouping、State / Operator の統合は後続の範囲とする。

## 13. Ordinary / graded の受入確認（2026-10-29）

Ordinary / graded の共通機構の実装時に行った受入確認を記録する。Fermion / spin-one の製品実装の検証は、それぞれ[局所代数の検証記録](local-structure.md#96-実装時の検証記録)・[Spin-one の検証記録](local-structure.md#106-実装時の検証記録)を参照する。共通機構の実装では新しい本番依存は追加していない。

| 対象 | 実装時の検証根拠 |
| --- | --- |
| kind・同一視・混在拒否・部分系・構築入口・fallible 優先順位 | [qmb-core/tests/tensor_products.rs](../../crates/qmb-core/tests/tensor_products.rs)、既存の spaces / fallible_construction / systems テスト |
| 符号・CAR・作用・随伴・内積・trace・Hilbert–Schmidt 内積 | テスト専用二状態代数と独立 Kronecker 行列による二・三因子の照合（f32 / f64、ordinary / graded） |
| Act scope・local embedding の準同型性・三つの odd 成分 | 同テストの scope / local expression / three odd factor ケース |
| intrinsic 座標と通常行列・再展開・列挙・空の基底 | 同テストと [qmb-core/tests/basis_expansion.rs](../../crates/qmb-core/tests/basis_expansion.rs)。非正規直交な三状態基底も両種類で検証 |
| 密疎 export と metadata | [qmb-export/tests/tensor_products.rs](../../crates/qmb-export/tests/tensor_products.rs) の Vec / ndarray / sprs、両精度の行列・遅延作用 vector |
| 正規化・期待値 | [qmb-normalize/tests/tensor_products.rs](../../crates/qmb-normalize/tests/tensor_products.rs) の独立行列との照合（両精度） |
| 最小能力・借用・共有・疎性 | 既存 capability / contract テストの移行、Product / 係数演算なしの graded 作用、80因子の分離状態と局所 trace、32因子の割り当て上限 |
| 深い木 | 既存 ordinary 回帰と graded の2万段 wrapper の構築・比較・評価・破棄（128 KiB stack） |
| spin-half・利用例 | 全 primitive が even、P=I、両種類での作用と trace の一致。[tensor_products.rs](../../examples/src/bin/tensor_products.rs) を実行 |

非自明な grading の参照代数は[テスト専用二状態代数](../../crates/qmb-core/tests/support/graded.rs)であり、製品用 fermion API ではない。`|0>` even、`|1>` odd、`a=|0><1|`、`P=diag(1,-1)` とし、二 mode の基底順は `|00>,|01>,|10>,|11>` とする。独立した Kronecker 行列で `c_1|11>=-|10>`、`c_0 c_1|11>=-|00>`、`c_1 c_0|11>=|00>`、局所・異なる mode の CAR、`Tr(c_1†c_1)=2` を確認した。参照行列は評価器・export を使わずに構築し、intrinsic 座標と通常行列を別々に照合した。近似的な零成分の削除は行わない。

Rust 1.98.0 の direnv devShell で次を実行し、成功した。

- `cargo test --workspace --all-features`（doctest を含む）
- `cargo check -p qmb-core -p qmb-geometry --no-default-features`
- `cargo clippy --workspace --all-features --all-targets -- -D warnings`
- `cargo fmt --all -- --check`
- `cargo run -p examples --bin tensor_products`

最終調整した構築入口と export の contract テストも再実行した。no_std の確認はホスト上のビルドであり、専用 target の検証ではない。MSRV 1.85 の実ビルド・専用 no_std target・性能ベンチマークは今回も未実施で、既存の [未実施事項](../future-work.md#3-未実施の検証と性能)として残す。
