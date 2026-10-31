# 演算子 primitive の名前・規格化の修正案

作成日: 2026-11-01
更新日: 2026-11-01

## 1. 採用範囲と実装状況

以下の名前と Hermite 演算子の規格化への移行を実装した。本書は採用した規約と移行内容を記録し、第10節に検証結果をまとめる。現行の規約は[要件の spin-half の定義](../../requirements.md#23-spin-12-の規約)、[局所 primitive の設計](../../design/local-structure.md)、各 crate の実装・rustdoc にも反映した。

| 対象 | 採用する変更 |
| --- | --- |
| `qmb-spin-half` | `Sx/Sy/Sz` を Pauli 行列として再定義し、`SPlus/SMinus` を `Sp/Sm` に変更する |
| `qmb-fermion` | `MajoranaX/Y/Z` を変更前の値の2倍に変更し、`Mx/My/Mz` とする。`Parity` は `Mz` に統合する |
| `qmb-spin-one` | Hermite な9要素 `I/Sx/Sy/Sz/Dxy/Q0/Qxy/Qyz/Qzx` を第5節の規約に揃え、演算子名の `Plus/Minus` を `p/m` に変更する |

spin-half の Cartesian 成分の公開名は `Sx/Sy/Sz` とし、変更前の同名 primitive の2倍に再定義する。旧規格化との互換性は維持しない。spin-one の Hermite primitive の名前も現行名を使い、意味を本書の規格化へ更新する。両 spin の昇降演算子名は `Plus/Minus` を `p/m` に短縮する。

状態 primitive・状態の内積・局所次元は維持する。fermion の `Annihilate/Create/Number` は現行の消滅・生成・占有数演算子を表す。これらまで一律に単位ノルムへ変更する案ではない。

モデルの結合定数は従来の物理的な意味を保ち、spin-one の平方根係数は当面 `Float` を要求して扱う方針も採用する。昇降演算子はノルムの統一を選択基準にせず、新しい演算子との相互変換・積・作用の簡潔さから第6節の定義を採用する。第7.3節の数値的相殺への対応も採用済みとする。

### 1.1. 採用する primitive の一覧

昇降演算子の公開名は `Sp/Sm`、`Q1p/Q1m`、`Q2p/Q2m` とする。spin-one の状態名 `Plus/Zero/Minus` は維持する。以下の variant 数と、それらが張る空間の次元は区別する。

| 代数 | 状態 primitive | 演算子 primitive |
| --- | --- | --- |
| spin-half | `Up, Down`（2個） | `I, Sx, Sy, Sz, Sp, Sm`（6個） |
| fermion | `Empty, Occupied`（2個） | `I, Annihilate, Create, Number, Mx, My, Mz`（7個） |
| spin-one | `Plus, Zero, Minus`（3個） | `I, Sx, Sy, Sz, Sp, Sm, Dxy, Q0, Qxy, Qyz, Qzx, Q1p, Q1m, Q2p, Q2m`（15個） |

演算子空間の次元はそれぞれ4・4・9である。zero vector・zero operator 用の variant は追加せず、既存の空の線型結合・式で表す。

## 2. 次元で割った HS 内積（採用）

局所 Hilbert 空間の次元を d とし、次で規格化を定義する。d は状態空間の次元であり、演算子空間の次元 d² や primitive の variant 数ではない。

```math
\langle A,B\rangle_d=\frac{1}{d}\operatorname{Tr}(A^\dagger B),\qquad
\|A\|_d=\sqrt{\frac{1}{d}\operatorname{Tr}(A^\dagger A)}.
```

対象の primitive に要求する条件は次である。

```math
\|A\|_d=1
\quad\Longleftrightarrow\quad
\operatorname{Tr}(A^\dagger A)=d.
```

I はそのまま条件を満たす。spin-half・fermion は d=2、spin-one は d=3 とする。固有値の絶対値を一律に1へ揃えることは共通条件にしない。

次元規格化は HS 内積とそこから導かれるノルムに適用する。通常の `Trace` は引き続き非正規化 trace を返し、局所空間では `Tr(I)=d`、全次元 D の複合空間では `Tr(I)=D` とする。複合空間上の HS 内積は全次元 D で割る。

期待値 `Tr(rho A)` や密度行列の純度 `Tr(rho²)` に使う trace も通常の trace とし、追加の次元除算は行わない。状態の内積・norm と `qmb-normalize` の状態正規化も既存の意味を維持する。任意の独自代数・全 `Basis` 実装に正規直交性を要求する共通契約は追加しない。

この規格化では、どの次元でも `||I||_d=1` となり、作用しない自由度を追加しても同じ演算子の内積を維持できる。

```math
\langle A\otimes I_m,B\otimes I_m\rangle_{dm}
=\frac{\operatorname{Tr}(A^\dagger B)\,m}{dm}
=\langle A,B\rangle_d.
```

Hermite 演算子 A の固有値を a_j とすれば、`||A||_d²=(1/d)Σ_j a_j²=Tr((I/d) A²)` であり、最大混合状態での二乗平均として解釈できる。作用しない自由度の次元による増大を除く規格化であり、相互作用項そのものを増やした場合のノルムまで一定になることは要求しない。

Hermite な完全基底が `Tr(B_mu B_nu)=d delta_mu,nu` を満たす場合、演算子の展開は次となる。

```math
A=\sum_{\mu=0}^{d^2-1}a_\mu B_\mu,\qquad
a_\mu=\frac1d\operatorname{Tr}(B_\mu A).
```

同じ局所基底を N 因子へ tensor 合成した基底は、全次元 D=d^N で割った HS 内積に関して正規直交となる。通常の HS ノルム二乗は全要素で D に揃う。fermion の graded 配置では、core が既存の parity string と符号を処理し、通常の状態基底上の作用行列でもこの直交性を確認する。

### 2.1. 局所 trace ごとの次元規格化（採用）

HS 内積は、現行の trace 評価が局所 trace の積を取る構造を利用し、各局所 trace を `Tr(I)` で割ってから掛け合わせる。全次元 D を `usize` として取得したり、非正規化した full trace を先に計算して D で割ったりする必要はない。

```math
A^\dagger B=\sum_t\gamma_t\bigotimes_f K_{tf},\qquad
\langle A,B\rangle_D
=\sum_t\gamma_t\prod_f\frac{\operatorname{Tr}(K_{tf})}{\operatorname{Tr}(I_f)}.
```

局所次元は、既存の `O: Identity + Trace<C>` から `O::identity()` の trace を係数型 C の値として得る。この計算に状態 primitive 型、`Dimension`、総次元の整数化、整数から C への変換を追加要求しない。係数型は局所次元を非零として表し、その値による数学的な除算を表現できることを前提とする。

- `LocalOperator::hs_inner_product` は、局所の `Tr(A†B)` を `Tr(I)` で割る。
- `Operator::hs_inner_product` と `hs_inner_product_in` は、既存の所属空間選択・積・graded から通常表現への変換を維持し、通常表現の各局所 trace を規格化する。parity string に現れる Mz を恒等演算子の spectator と取り違えない。
- 作用しない因子では `Tr(I)/Tr(I)=1` となる。空の tensor 空間は D=1 とし、係数の縮約だけで計算する。
- 係数の除算には既存の `qmb_core::Divide` を利用する。`Complex<R>` では実数の分母を不要に二乗しない既存実装を使える。必要な bound は HS 内積側に限定し、通常の `trace` / `trace_in` に除算能力を追加要求しない。

現行の入口と局所 trace の走査は [evaluate/trace.rs](../../../crates/qmb-core/src/evaluate/trace.rs)、除算能力は [coefficient.rs](../../../crates/qmb-core/src/coefficient.rs)を参照する。通常 trace と HS 内積で内部の走査を共有する場合も、公開の trace の意味と能力要件を維持する。公開の正規化 trace API は追加しない。

これにより総次元の `usize` overflow と、作用しない因子の次元を掛け重ねることによる中間値の増大を避ける。一般の係数の桁あふれ・丸め誤差まで防ぐ保証ではない。

## 3. Spin-half

状態順を現行どおり `(Up, Down)` とし、定義を次とする。

```math
S_x=\begin{pmatrix}0&1\\1&0\end{pmatrix},\qquad
S_y=\begin{pmatrix}0&-i\\i&0\end{pmatrix},\qquad
S_z=\begin{pmatrix}1&0\\0&-1\end{pmatrix}.
```

公開名はそれぞれ `Sx/Sy/Sz`。変更前の同名 primitive の2倍であり、固有値は ±1 となる。物理的な角運動量を J と書けば、ℏ=1 で `J_alpha=S_alpha/2` である。

```math
S_\alpha^2=I,\qquad
S_xS_y=iS_z,\qquad
[S_x,S_y]=2iS_z,\qquad
S_x^2+S_y^2+S_z^2=3I.
```

循環置換した関係も同様とする。`CartesianBasis` の要素は `[I,Sx,Sy,Sz]` の順を維持し、正規化 HS Gram 行列を単位行列にする。`Up/Down` は引き続き物理的な磁気量子数 ±1/2 の状態で、`Sz` の固有値は ±1 となる。状態基底名は `SzBasis` を維持する。

`Sp/Sm` は第6.1節の `(Sx±iSy)/2` を採用する。`LadderBasis` の要素は `[I,Sp,Sm,Sz]` の順とし、昇降演算子との展開係数を更新する。

## 4. Fermion

状態順を `(Empty, Occupied)`、消滅演算子を a、占有数を `n=a†a` とする。

```math
M_x=a+a^\dagger,\qquad
M_y=\frac{a-a^\dagger}{i},\qquad
M_z=-iM_xM_y=I-2n.
```

公開名は `Mx/My/Mz`。占有基底での行列は第3節の Pauli 行列と同じで、変更前の `MajoranaX/Y/Z` の2倍となる。`Mx/My` は odd、`Mz` は even とし、Z 成分を第三の odd Majorana として扱わない。

```math
M_\alpha^2=I,\qquad
[M_x,M_y]=2iM_z,\qquad
\{M_\alpha,M_\beta\}=2\delta_{\alpha\beta}I.
```

ここで3成分の関係は同一 mode 内の通常の局所積について述べている。異なる mode 間の符号は既存の grading に従う。X・Y の規格化は標準的な Majorana の規約に対応する。

基底変換の関係は次となる。

```math
a=\frac{M_x+iM_y}{2},\qquad
a^\dagger=\frac{M_x-iM_y}{2},\qquad
n=\frac{I-M_z}{2},\qquad
P=M_z.
```

`MajoranaBasis` は `[I,Mx,My,Mz]` の順とし、正規化 HS Gram 行列を単位行列にする。`LadderBasis` は `[I,Annihilate,Create,Number]` で、従来どおり非直交な基底である。

`Parity` は `Mz` の名前で同じ値へ統合し、独立 variant・互換用の別名は設けない。`ParityOperator::parity_operator()` は `Mz` を返す。`Mz²=I` と既存の graded 合成契約を満たし、`MajoranaBasis` では `Mz` の index 3を使う。`LadderBasis` では `Mz` 自身の index はなく、`I−2n` に展開する。variant 数は8から7へ減るが、演算子空間の次元は4のままとする。

## 5. Spin-one の Hermite primitive

### 5.1. 物理的なスピンとの対応

ℏ=1 の物理的な spin-1 演算子を J と書く。状態順は現行どおり `(Plus, Zero, Minus)` とする。

```math
J_x=\frac1{\sqrt2}\begin{pmatrix}0&1&0\\1&0&1\\0&1&0\end{pmatrix},\qquad
J_y=\frac1{\sqrt2}\begin{pmatrix}0&-i&0\\i&0&-i\\0&i&0\end{pmatrix},\qquad
J_z=\begin{pmatrix}1&0&0\\0&0&0\\0&0&-1\end{pmatrix}.
```

`c=√(3/2)` と置き、v3 の Hermite primitive を次で定義する。右辺の J は物理的なスピンであり、左辺の規格化した S と区別する。反交換子は `{A,B}=AB+BA` とする。

| 公開 primitive | v3 の定義 | 変更前の値に掛ける係数 |
| --- | --- | --- |
| `I` | 恒等演算子 | 1 |
| `Sx/Sy/Sz` | `S_alpha=c J_alpha` | `√(3/2)` |
| `Dxy` | `c (Jx²−Jy²)` | `√(3/2)` |
| `Qxy/Qyz/Qzx` | それぞれ `c {Jx,Jy}`、`c {Jy,Jz}`、`c {Jz,Jx}` | `√(3/2)` |
| `Q0` | `(3 Jz²−2I)/√2` | `1/√2` |

特に、採用する `Sx` の行列は次である。

```math
S_x=\frac{\sqrt3}{2}\begin{pmatrix}0&1&0\\1&0&1\\0&1&0\end{pmatrix},\qquad
\operatorname{Tr}(S_x^\dagger S_x)=3,\qquad \|S_x\|_3=1.
```

`Sx/Sy/Sz` の固有値は `−c,0,c`、`Q0` の対角成分は `(1,−2,1)/√2` となる。状態名の `Plus/Zero/Minus` は物理的な `Jz` の磁気量子数 `+1,0,−1` に対応し、変更後の `Sz` 固有値をそのまま磁気量子数と呼ばない。

### 5.2. 基底と局所規則

`CartesianBasis` の順序は `[I,Sx,Sy,Sz,Dxy,Q0,Qxy,Qyz,Qzx]` を維持する。全9要素が Hermite で、通常の HS Gram 行列は `3 I_9`、正規化 HS Gram 行列は `I_9` となる。

変更前のトレースレスな8要素のうち Q0 をまず `Q0/√3` とすれば、すべて通常の HS ノルム二乗が2になる。本案はこの8要素全体を c 倍する規格化である。スピン3成分と四重極5成分の回転に関する部分空間を保つ。

変更後の代表的な関係は次となる。

```math
[S_x,S_y]=icS_z,\qquad
S_x^2+S_y^2+S_z^2=3I,\qquad
S_z^2=I+\frac{Q_0}{\sqrt2},\qquad
Q_0^2=I-\frac{Q_0}{\sqrt2}.
```

```math
D_{xy}=\frac1c(S_x^2-S_y^2),\qquad
Q_{xy}=\frac1c\{S_x,S_y\},\qquad
Q_0=\sqrt2(S_z^2-I).
```

Qyz・Qzx も同じ係数で定義される。変更前の `SPlus²=Q2Plus` は、v3 の primitive では `Sp²=c Q2p` となる。各倍率から多項式関係を再導出する。随伴は Hermite な9要素を保存し、通常の trace は I が3、他の8要素が0となる。

## 6. 非 Hermite な昇降演算子（採用）

新しい演算子との線型関係、積による関係、状態間の遷移振幅を基準に以下の定義を採用する。昇降演算子のノルムを1にすることは選択基準にしない。

### 6.1. Spin-half

`Sp/Sm` を数式では `S_±` と書く。

```math
S_\pm=\frac{S_x\pm iS_y}{2},\qquad
S_x=S_++S_-,\qquad
S_y=\frac{S_+-S_-}{i}.
```

`Sp Sm=(I+Sz)/2`、`Sm Sp=(I−Sz)/2` が射影演算子となる。`Sp`・`Sm` 自体は行列単位で、非零な遷移振幅は1であり、物理的な spin-half の昇降演算子とも一致する。行列としては fermion の a・a† と対応するが、grading と型は引き続き区別する。`[Sz,Sp]=2Sp`、`[Sz,Sm]=−2Sm`、`Sp²=Sm²=0` を満たす。

### 6.2. Spin-one

Cartesian 基底との直接の相互変換を優先し、公開名 `Sp/Sm`、`Q1p/Q1m`、`Q2p/Q2m` で次を定義する。数式ではそれぞれ `S_±`、`Q_1±`、`Q_2±` と書く。

```math
S_\pm=S_x\pm iS_y,\qquad
Q_{1\pm}=Q_{zx}\pm iQ_{yz},\qquad
Q_{2\pm}=D_{xy}\pm iQ_{xy}.
```

この案では全6要素が変更前の値の c 倍となり、Cartesian 基底との変換は和差・1/2・i だけで表せる。`S_±²=c Q_2±`、`{Sz,S_±}=c Q_1±` と、多項式からの定義には c が付く。状態への `S_±` の非零な遷移振幅は √3 となる。

基底展開を重視し、単純な線型変換を優先する。多項式関係の係数は文書・テストで明示する。通常の HS ノルム二乗は全6要素で6、正規化 HS ノルムは √2 となるが、その値は定義の結果として扱う。`LadderBasis` の順序は `[I,Sp,Sm,Sz,Q0,Q1p,Q1m,Q2p,Q2m]` とし、通常の HS Gram 行列は `diag(3,6,6,3,3,6,6,6,6)` となる。

実装時には、これらの定義に合わせて基底変換・随伴・積・作用・LadderBasis の Gram 行列を一括で更新する。

## 7. モデル・係数能力への影響

### 7.1. 物理的な結合定数の維持（採用）

既存モデルの結合定数と物理的な Hamiltonian の対応を保つ。結合定数を j・k・U、新しい spin-one の `A_ij=Σ_alpha S_i^alpha S_j^alpha` と書き、モデル構築側で以下の係数補正を行う。

| モデル | 新 primitive での表現 |
| --- | --- |
| spin-half XYZ / Heisenberg | 変更前の各 `−j_alpha J_i^alpha J_j^alpha` を `−(j_alpha/4) S_i^alpha S_j^alpha` とする |
| spin-one XYZ / Heisenberg | 変更前の各 `−j_alpha J_i^alpha J_j^alpha` を `−(2j_alpha/3) S_i^alpha S_j^alpha` とする |
| spin-one BLBQ | 各葉ペアで `−[(2j/3) A_ij+(4k/9) A_ij²]` とする |
| fermion symmetric Hubbard | 変更前の `U (n_up−I/2)(n_down−I/2)` を `(U/4) Mz_up Mz_down` とする |
| fermion standard Hubbard | a・a†・n の定義を保つため、Hamiltonian の係数は維持する |

現行のモデルと配置の契約は[Chain・Site の設計](../../design/geometry.md)、symmetric Hubbard の実装は [models.rs](../../../crates/qmb-fermion/src/models.rs)を参照する。モデルの符号、境界条件、葉ペアの列挙と定数項も合わせて検証する。

モデル構築で数値係数を直接補正する場合、変更前の `Clone + From<R> + Neg` だけでは分数の生成を保証できない。係数補正に必要な One・算術・除算等の trait bound は、実装時に演算ごとに適宜要求し、公開関数の能力要件とテストを更新する。std 等の feature 要件に影響する変更は第9節に従って別途検討する。

### 7.2. Spin-one の平方根係数（Float を採用）

変更前の spin-one は局所積と基底変換を複素有理係数で扱え、作用だけが `Complex<R>, R: Float` を使う。本案では Hermite primitive 同士の積にも c・√2 等が現れるため、変更前の有理係数だけの検証をそのまま維持することはできない。

当面、平方根係数を必要とする局所積・作用は `Complex<R>`、`R: num_traits::Float` を対象とする。厳密な代数的係数のための新しい能力 trait は今回追加しない。平方根を必要としない基底変換・随伴・trace 等は操作に必要な bound を維持する。式の保持・遅延構築へ一律の `Float` 制約を追加しない。定数は R の One・算術・sqrt から選択した精度で生成し、f64 経由の変換やパースを前提にしない。

### 7.3. 浮動小数点の相殺への対応（採用）

通常の評価・基底展開・座標取得は、現行どおり厳密な `Zero::is_zero()` のみで零係数を除く。微小な物理的結合を保存し、ユーザーが指定していない許容誤差で式の値を変更しない。この方針は[数値比較の要件](../../requirements.md#56-数値的近似等価性と検証の規約)と、[現行の座標集約](../../../crates/qmb-core/src/basis/expansion/mod.rs)に沿う。以下の対応を採用する。

- 局所 Product・Act の係数は数学的に導出・整理した最終形を直接返す。数学的に零の結果は空の線型結合とし、既知の零項を出力しない。数学的に1となる係数には One を使い、例えば c² を計算してから割り戻す経路で単位係数を作らない。
- 定数の生成式を private helper で揃える。旧 primitive への数値変換・積・新 primitive への逆変換を実行時に繰り返す実装は採らず、局所規則を直接記述する。
- 一般の式に対して、数学的な零が浮動小数点でも常に厳密零になるとは保証しない。基底の往復変換・結合則・一般の式の相殺は、行列・作用・座標の数値で比較する。近似的に同じ式で疎な項数まで一致することは要求しない。
- 検証は複素数の絶対値に対する `|a−b| <= atol + rtol max(|a|,|b|)` を基本とする。現行の小規模 spin-one テストの f32 の `atol=rtol=2e-5`、f64 の `atol=rtol=1e-12` を出発点に、計算規模・係数の大きさに合わせて確認する。積の連鎖・基底の往復、大小の係数を含む相殺、微小だが非零な項を保存するケースを加える。既知の零を返す局所規則の空性は厳密に検証する。
- 相殺残差に対する状態の正規化も既存の厳密零判定に従う。数学的には零の式でも計算結果が厳密零でなければ `ZeroNorm` を保証できず、残差を正規化するとその成分が拡大され得る。この制約を通常の浮動小数点評価の限界として扱う。

近似的な成分の切り捨てが必要な場合は、将来の明示的な後処理として検討する。基底・許容誤差・適用時点を利用者が選ぶ操作とし、通常の評価や export に暗黙に組み込まない。本変更では新しい切り捨て API の実装を必須にしない。

正規化 HS 内積で正規直交な演算子基底において、捨てる成分を `delta A=Σ_mu delta a_mu B_mu` とすると `||delta A||_D²=Σ_mu |delta a_mu|²` である。後処理を導入する場合は、単一係数の閾値だけでなく、削除した成分全体の誤差をこの式で評価できる。非正規直交な LadderBasis にこの式をそのまま適用しない。

## 8. 実装範囲と検証方針

修正対象は qmb-core の局所・複合 `hs_inner_product` と `hs_inner_product_in`、各代数 crate の enum・局所 `Product/Act/Adjoint/Trace`・基底変換、モデル、テスト、rustdoc・examples・現行設計文書である。基底の列挙順と座標の意味を対応付け、名前だけの機械的な置換で終わらせない。`Basis` が primitive の部分集合を基底にする契約は維持する。

現行の [spin-half](../../../crates/qmb-spin-half/tests/contracts.rs)、[fermion](../../../crates/qmb-fermion/tests/contracts.rs)、[spin-one](../../../crates/qmb-spin-one/tests/contracts.rs) の独立行列参照と能力別テストを更新する。受入時には以下を確認する。

- 参照行列を物理的な J・a と採用倍率から独立に構成し、全 primitive の積・作用・随伴・trace を両精度で照合する。
- 対象の Hermite 基底の正規化 HS Gram 行列が単位行列であること、特に spin-one の Sx が第5.1節の行列であることを確認する。
- 積の結合則、作用との整合性、交換・反交換関係、Casimir、parity の法則と mode 間の CAR を確認する。
- 全 primitive の両基底への展開・往復・座標を確認し、第6節の昇降演算子の定義との整合性を検証する。
- tensor 基底の正規化 HS Gram 行列と、graded な通常作用行列の符号を確認する。演算子基底係数と状態基底上の行列成分を区別する。
- 第2.1節の局所規格化を、小規模系の明示行列による `Tr(A†B)/D` と照合する。空の空間、identity spectator の追加、明示空間による束縛の上書き、fermion の parity string、総次元が `usize` に収まらない系での恒等演算子の HS 内積1を確認する。通常 trace の値・除算不要の能力要件も維持する。
- モデルは独立な物理的 Hamiltonian と照合し、結合定数・エネルギー・定数シフトの維持を確認する。
- 第7.3節の相殺・微小非零係数の保存を検証し、局所規則の厳密零と一般の式の数値的な零を区別する。
- export の Vec / ndarray / sprs、最小係数能力、構築時の能力要件、公開名の doctest と examples を更新する。Pauli string の既存例にある因子ごとの2倍を取り除く。

workspace の関連テスト・doctest、feature 構成、format、Clippy、ドキュメント生成の検証結果を第10節に記録する。

## 9. 補助 API と trait bound の実装方針（採用）

物理的なスピンを構築する補助 API は追加しない。spin-half の `J_alpha=S_alpha/2`、spin-one の `J_alpha=√(2/3) S_alpha` は既存のスカラー倍で表し、規格化した primitive と物理量の対応をドキュメント・examples で示す。

モデルの係数補正を含め、各操作の実装に必要な trait bound は実装時に適宜要求する。係数型の能力はそれを使う操作に限定し、式の保持・遅延構築へ一律の数値制約を追加しない。std / no_std や依存 crate の feature 要件に影響する変更は、この trait bound の具体化とは分けて検討する。

本書の primitive と規格化、HS 内積の局所 trace ごとの次元規格化と通常の Trace の維持、`Mz` への統合、モデルの結合定数の物理的な意味、平方根係数の Float による扱い、第7.3節の相殺への対応、および本節の実装方針は実装済みとする。現時点で本件の設計上の未決事項はない。

## 10. 実装と検証の記録

2026-11-01、Rust 1.99.1 の direnv devShell で実装と受入検証を行った。以下のコマンドはリポジトリ直下で `direnv exec .` を付けて実行し、すべて成功した。

```sh
cargo test --workspace --all-features --no-fail-fast --quiet
cargo test --workspace --no-default-features --offline --no-fail-fast --quiet
cargo test -p qmb-spin-one -p qmb-fermion --all-features --offline --quiet
cargo test -p qmb-spin-one --test tensor_products -p qmb-fermion --test models --offline --quiet
cargo check -p qmb-core -p qmb-geometry --no-default-features --offline
cargo fmt --all --check
cargo clippy --workspace --all-targets --all-features --offline -- -D warnings
cargo clippy --workspace --all-targets --no-default-features --offline -- -D warnings
sh -c 'RUSTDOCFLAGS="-D warnings" cargo doc --workspace --all-features --no-deps --offline'
```

qmb-export は単独でも `cargo test -p qmb-export --no-default-features --offline --quiet` を実行し、追加 feature なし・`--features ndarray`・`--features sprs` の3構成を確認した。両 backend 有効の構成は workspace の全 feature テストで確認した。

確認した内容は次のとおり。

- 各代数の全 primitive の積・作用・随伴・通常 trace、両演算子基底の変換と座標を、独立行列と f32 / f64 で照合した。spin-one の全225積・45作用、結合則・作用との整合性も検証した。
- 局所の正規化 HS Gram 行列に加え、spin-one の Cartesian / Ladder tensor 基底と fermion の graded Majorana tensor 基底の Gram 行列を検証した。graded 評価の符号・CAR と通常作用行列も照合した。
- [HS 内積の回帰テスト](../../../crates/qmb-core/tests/hs_inner_product.rs)で、明示行列による `Tr(A†B)/D`、空空間、spectator の追加、明示空間による束縛の上書き、および2048因子での `⟨I,I⟩=1` を確認した。状態型・Dimension・整数変換を要求しないこと、Divide だけを持つ係数でも計算できること、通常 trace に除算能力を要求しないことも検証した。
- 物理的なモデルの独立参照行列とエネルギー・定数シフトを維持した。微小係数の保存と厳密な相殺、解析的に零の局所規則の空性も確認した。暗黙の近似的な切り捨ては追加していない。

examples の全8バイナリを `cargo run -p examples --bin <name> --offline --quiet` により、既定 feature・`--all-features`・`--no-default-features` の3構成で実行した。一重項・三重項のエネルギー −0.75 / 0.25、symmetric Hubbard のエネルギー −2、spin-one BLBQ の期待値0.5、Pauli string の37項と密疎 export の一致が保たれた。

依存関係と feature 要件は変更していない。MSRV 1.85 と専用 no_std target での実ビルドは今回実行しておらず、[未実施の検証](../../future-work.md#3-未実施の検証と性能)として引き続き扱う。
