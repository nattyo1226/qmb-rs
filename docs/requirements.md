# qmb-rs 初期版要件

## 1. 文書の目的

本文書は、qmb-rs の目的、基本概念、初期版の開発範囲、採用する規約、および完成条件を定める。

## 2. 中心命題

qmb-rs は、物理的な自由度とその関係から、表現非依存の Hilbert 空間、状態、および演算子を構築する。それらに対して代数演算を行い、必要に応じて具体的な数値表現へ変換できるようにする。

初期版は、以下の三層を明確に分離する。

1. 物理空間と幾何構造
2. Hilbert 空間上の代数的な状態と演算子
3. 密ベクトル、疎ベクトル、密行列、および疎行列への変換

上位層の概念や操作は、下位層の具体的な数値表現に依存してはならない。特に、代数的な状態や演算子を扱う際に、密行列や疎行列の生成を必要としないことを基本方針とする。

## 3. 初期版の代表ユースケース

初期版は、次の処理を一貫して実行できることを代表的な縦断シナリオとする。

1. 長さ `L` の一次元 spin chain を構築する。
2. 各 site に spin-1/2 の局所自由度を配置し、全 Hilbert 空間を構築する。
3. nearest-neighbor 関係から Heisenberg Hamiltonian を構築する。
4. 各 site の局所状態を交互に割り当て、Neel 積状態を構築する。
5. Hamiltonian の状態への代数的な作用を通じて、行列を生成せずにエネルギー期待値を計算する。
6. 同じ状態と演算子を密および疎の数値表現へ変換する。
7. 小規模系において、代数計算、密表現、および疎表現から得られる結果が一致することを検証する。

等方的 Heisenberg Hamiltonian を

```math
H = J \sum_{\langle i,j\rangle}
    \left(S_i^x S_j^x + S_i^y S_j^y + S_i^z S_j^z\right)
```

とする。open chain の Neel 状態では、標準規約の下でエネルギー期待値が

```math
\langle \mathrm{Neel}|H|\mathrm{Neel}\rangle
= -\frac{J}{4}(L-1)
```

となることを基準ケースとする。

この基準ケースでは、結合定数 `J` は実数とする。代数オブジェクト自体は複素数係数を許容するが、複素数の `J` を与えた一般の場合に Hamiltonian が Hermitian であることは要件としない。

## 4. アーキテクチャ上の責務

### 4.1. 物理空間と幾何構造

物理空間は、自由度を置く場所と、それらの間に存在する幾何学的関係を表す。物理空間そのものは、各 site が spin、fermion、または別種の自由度を持つかを決定しない。

初期版では以下を扱う。

- 有限個の site
- 安定した一意な site 識別子
- site の明示的な順序
- 一次元 chain
- open boundary condition
- periodic boundary condition
- nearest-neighbor 関係
- 関係の列挙

一次元 chain は、一般的な site 集合と関係からなる物理空間を構築するための便利な生成手段として位置づける。幾何学的な関係と Hamiltonian の結合は同一視しない。物理空間は関係を提供し、モデル構築処理がその関係に演算子と係数を割り当てる。

向きのない二 site 関係は、明示された site ordering において先に現れる site を `min(i,j)`、後に現れる site を `max(i,j)` として、常に `(min(i,j), max(i,j))` の形へ正規化する。関係を列挙するときは、正規化した第一要素、第二要素の順に site ordering に基づく辞書順で並べる。

将来、三体相互作用、plaquette、および hyperedge を表現できるように、関係の概念を単純な最近接ペアだけに不必要に固定しない。ただし、それらの実装は初期版の対象外とする。

### 4.2. Hilbert 空間

初期版では以下を扱う。

- 有限次元の局所 Hilbert 空間
- spin-1/2 の局所 Hilbert 空間
- 局所 Hilbert 空間の tensor product
- site と局所 Hilbert 空間の対応
- tensor product を構成する因子の明示的な順序
- 局所基底 `|up>` および `|down>`

一つの系を構成するすべての局所 Hilbert 空間は、同一種類の局所代数を使用する。初期版で builtin algebra として公開する局所代数は spin-1/2 のみとし、異なる種類の局所代数を同一系に混在させない。spin-1/2 以外の有限次元局所空間を将来追加できるように、系全体で使用する局所代数を差し替え可能とする。局所代数を差し替え可能であることを検証するため、公開対象ではない test 用の有限次元局所代数を実装してよい。

### 4.3. 状態

状態は、数値ベクトルではなく、自由度の識別子を付与した局所状態からなる代数的なオブジェクトとして表現する。空間依存の検証または評価に必要な全 Hilbert 空間は、必要な操作の引数として与える。

初期版では、状態式の表現要素として以下を扱う。

- 局所基底状態
- 基底積状態
- 一般の積状態
- 状態の有限線形結合

状態式に対して初期版で公開する代数演算は、4.5.1 節に定める以下の操作とする。

- 状態の加算
- 状態と複素 scalar の scalar multiplication
- 状態間の tensor product

内積および norm は、状態式を消費する評価操作として扱う。演算子の状態への作用、期待値、および正規化期待値も評価操作とし、その意味論は 7 節に定める。

Neel 状態は専用の状態型とせず、site ごとに局所状態を交互に割り当てた積状態として表現する。Neel 状態を作る便利な構築処理を提供することは妨げない。先頭 site を `|up>` とするか `|down>` とするかは、構築時に明示できることが望ましい。

### 4.4. 演算子

演算子は、数値行列ではなく、自由度の識別子を付与した primitive な局所演算子からなる代数式として表現する。演算子式自身は全作用 Hilbert 空間を保持しない。全作用空間、site 順序、および恒等的に作用する自由度は、状態への作用、trace、または数値化などで必要になった時点で、状態または明示的な評価 context から与える。

初期版では、演算子式の表現要素として以下を扱う。

- 恒等演算子
- 局所 spin 演算子 `Sx`, `Sy`, `Sz`, `S+`, `S-`
- 異なる site に作用する局所演算子の積
- 演算子項の有限和

演算子式に対して初期版で公開する代数演算は、4.5.2 節に定める以下の操作とする。

- 演算子の加算
- 演算子の乗算
- 演算子と複素 scalar の scalar multiplication
- 演算子間の tensor product
- full trace
- partial trace

状態への作用は 7 節に定める評価操作として扱う。同類項の結合および局所 spin 代数に基づく基本的な簡約は、9.1 節に定める lazy simplification の一部として扱う。adjoint および commutator を公開操作として提供することは初期版の対象外とする。

Heisenberg Hamiltonian は特別な内部表現を持つ演算子型とせず、一般的な局所演算子積の有限和として表現する。nearest-neighbor 関係に演算子項を割り当てるモデル構築用の便利な処理を、その一般表現の上に提供できるものとする。

演算子式の support は、式中の primitive な局所演算子に付与された自由度識別子から導出できるものとする。support または最小作用空間を、すべての複合式へ必須 metadata として保存することは要件としない。

演算子積を常に完全展開すると項数が増大し得るため、代数式の展開と簡約は初期版から遅延して行う。演算子の構築や合成の時点で式全体を展開または正規化してはならない。ただし、あらゆる式に対する最適な簡約結果や計算量を保証することは、初期版の要件としない。

### 4.5. 代数演算の機能要件

#### 4.5.1. 状態に対する演算

初期版は、状態に対して以下の演算を提供する。

- 同一 Hilbert 空間上の二状態の加算
- 状態と複素 scalar の scalar multiplication
- 二状態の tensor product

状態の加算結果は、同じ Hilbert 空間上の状態の有限線形結合とする。同類項の結合および zero term の除去は lazy simplification により行う。

tensor product の結果は、左 operand の自由度を左側、右 operand の自由度を右側に配置した複合状態とする。二 operand が同じ自由度識別子を含む場合、または異なる種類の局所代数を使用する場合はエラーとする。

状態式だけから対象空間を一意に決められない場合は、検証または評価時に Hilbert 空間を明示的に与える。特に、簡約後の zero state から所属空間を推論できない場合は、Hilbert 空間 context を必須とする。

#### 4.5.2. 演算子に対する演算

初期版は、演算子に対して以下の演算を提供する。

- 演算子式同士の加算
- 演算子と複素 scalar の scalar multiplication
- 演算子の通常の積
- 二演算子の tensor product
- 演算子の full trace
- 演算子の partial trace

演算子は全作用 Hilbert 空間を保持しないため、異なる自由度を参照する演算子式も形式的に加算できる。加算結果の support は各 operand の support の和集合として導出できる。評価時には、各 primitive が参照する自由度を含む Hilbert 空間上へ、恒等演算子を用いて暗黙に拡張する。

同一の局所代数と数値精度を使用する演算子式だけを加算できる。式が参照する自由度の存在および局所代数との整合性は、状態への作用または Hilbert 空間上での評価時に検証する。

演算子の tensor product は通常の演算子積とは異なる演算として提供し、左 operand を左側、右 operand を右側に配置する。二 operand が同じ自由度識別子を含む場合、または異なる種類の局所代数を使用する場合はエラーとする。同一自由度上の演算子を続けて作用させる場合は、tensor product ではなく通常の積を使用する。

#### 4.5.3. Zero、identity、および空の tensor product

zero state、zero operator、および identity operator を明示的に構築するときは Hilbert 空間を指定せず、site、site 順序、または局所 Hilbert 空間の情報を紐付けない。zero term または identity factor は lazy simplification によって式から除去できる。これらは、評価時に与えられる Hilbert 空間上の zero または identity として解釈する。zero state の所属空間を式から推論できない場合は、評価時に Hilbert 空間 context を必須とする。

primitive な局所演算子を含まない scalar `c` は、評価対象の Hilbert 空間 `H` 上では

```math
c I_{\mathcal H}
```

として解釈する。scalar `0` は zero operator、scalar `1` は identity operator に対応する。空の tensor product は scalar `1` とする。

0-site Hilbert 空間は一次元とする。

```math
\dim \mathcal H_{\varnothing}=1, \qquad
\operatorname{Tr}(I_{\varnothing})=1
```

#### 4.5.4. Full trace と partial trace

演算子 `A` の full trace は、引数として与えた有限次元 Hilbert 空間 `H` 上で

```math
\operatorname{Tr}_{\mathcal H}(A)
= \sum_k \langle e_k|A|e_k\rangle
```

と定義する。ここで `{|e_k>}` は `H` の任意の完全正規直交基底である。full trace は非正規化 trace とし、複素 scalar を返す。trace の値は基底の選択に依存してはならない。

trace は代数的な演算として提供し、完全な密行列または疎行列の生成を必要としてはならない。局所演算子の tensor product に対しては、可能な場合に

```math
\operatorname{Tr}(A \otimes B)
= \operatorname{Tr}(A)\operatorname{Tr}(B)
```

を利用できるものとする。

partial trace では、次の二つを明示的に与える。

- `ambient`: 演算子を解釈する全 Hilbert 空間
- `trace_out`: trace out する自由度の集合

`trace_out` は `ambient` の部分集合でなければならない。演算子式が参照するすべての自由度も `ambient` に存在しなければならない。結果は `ambient` から `trace_out` を除いた自由度上で解釈される演算子式とし、残存自由度の順序は `ambient` の順序を保存する。結果の演算子式自身には、残存 Hilbert 空間を必須 metadata として保存しない。

演算子が trace out する自由度上で恒等的に作用する場合、その局所空間の次元を係数へ掛ける。特に、

```math
\operatorname{Tr}_{\mathcal H_1}
\left(A_0;\mathcal H_0\otimes\mathcal H_1\right)
= \dim(\mathcal H_1)A_0
```

とする。空集合を trace out する場合は元の演算子と等価な式を返す。すべての自由度を trace out する場合、partial trace の結果は 0-site Hilbert 空間上の scalar operator expression とする。full trace は同じ数学的結果を複素 scalar として返す。

状態ベクトルから reduced density operator を構築する操作は初期版に含めない。

## 5. 数学的および数値的規約

### 5.1. Spin 演算子

`hbar = 1` とし、Pauli 行列との関係を

```math
S^\alpha = \frac{1}{2}\sigma^\alpha
```

と定める。ladder operator は

```math
S^+ = S^x + iS^y, \qquad
S^- = S^x - iS^y
```

とする。

局所基底の順序は

```math
(|\uparrow\rangle, |\downarrow\rangle)
```

とし、それぞれ `Sz` の固有値 `+1/2` および `-1/2` に対応させる。この基底における局所行列表現は次のとおりとする。

```math
S^x = \frac{1}{2}
\begin{pmatrix}0&1\\1&0\end{pmatrix}, \qquad
S^y = \frac{1}{2}
\begin{pmatrix}0&-i\\i&0\end{pmatrix}, \qquad
S^z = \frac{1}{2}
\begin{pmatrix}1&0\\0&-1\end{pmatrix}
```

```math
S^+ =
\begin{pmatrix}0&1\\0&0\end{pmatrix}, \qquad
S^- =
\begin{pmatrix}0&0\\1&0\end{pmatrix}
```

### 5.2. Tensor product と基底順序

`N` site 系の Hilbert 空間を

```math
\mathcal H = \mathcal H_0 \otimes \mathcal H_1
\otimes \cdots \otimes \mathcal H_{N-1}
```

とし、site 0 を tensor product の最も左側に置く。複合基底状態は

```math
|s_0,s_1,\ldots,s_{N-1}\rangle
= |s_0\rangle\otimes|s_1\rangle\otimes\cdots
\otimes|s_{N-1}\rangle
```

とする。

数値表現では、右端の site の局所基底 index が最も速く変化する。site `i` の局所次元を `d_i`、局所基底 index を `s_i` としたとき、global basis index を

```math
\operatorname{index}(s_0,\ldots,s_{N-1})
= \sum_{i=0}^{N-1} s_i \prod_{j=i+1}^{N-1} d_j
```

と定める。

spin-1/2 の二 site 系では、基底順序は次のとおりとなる。

1. `|up, up>`
2. `|up, down>`
3. `|down, up>`
4. `|down, down>`

ここで実際の zero-based global basis index は順に `0`, `1`, `2`, `3` である。密表現と疎表現は同じ基底順序を使用する。

### 5.3. 演算子積の作用順

`A * B` は数学的な積 `AB` を表し、状態には右側から作用する。

```math
(AB)|\psi\rangle = A(B|\psi\rangle)
```

したがって、`A * B` では `B` が先に、`A` が後に作用する。

異なる site に作用する spin 演算子は可換として正規化できる。同一 site 上の演算子については順序を保存し、局所 spin 代数に従って処理する。

### 5.4. スカラーと精度

代数オブジェクトのスカラーは複素数に統一し、以下を扱う。

- 単精度複素数
- 倍精度複素数

実数値は、対応する複素数値に変換して受け入れてよい。記号係数は初期版では扱わない。Heisenberg Hamiltonian の結合定数 `J` も、具体的な複素数値として与える。

左 scalar multiplication と右 scalar multiplication は同じ意味とする。単精度と倍精度のオブジェクトを同じ演算で混在させず、暗黙の精度昇格も行わない。精度を変更する場合は明示的な変換を必要とする。

使用する数値型が表現できる場合、NaN と正負の infinity を scalar として許容する。これらを含む scalar の算術、厳密な比較、および特殊値の伝播は、使用する数値型の定義に従う。ただし、5.5 節で明示的に定義する近似比較は、その規約に従う。qmb-rs は NaN または infinity に対する独自の正規化、拒否、例外化、または置換規則を設けない。代数的な簡約においても、数値型による演算結果を変更してはならず、演算を省略することでその数値型とは異なる結果を生じさせてはならない。

### 5.5. 等価性と許容誤差

構造的な等価性と数値的な近似等価性を区別する。

- 構造的等価性: 項、site、局所状態または局所演算子、および係数の表現が厳密に一致すること。評価に使用する Hilbert 空間は比較対象に含めない
- 数値的近似等価性: 対応する数値が指定された許容誤差の範囲で一致すること
- 作用としての等価性: 同一の Hilbert 空間上で同じ線形作用を表すこと

近似比較には、絶対許容誤差 `atol` と相対許容誤差 `rtol` を用い、複素数 `a`, `b` に対して

```math
|a-b| \leq \mathrm{atol}
       + \mathrm{rtol}\max(|a|,|b|)
```

を満たす場合に近似的に等しいと判定する。

近似等価性は推移律を満たさないため、通常の厳密な等値判定とは別の明示的な操作とする。NaN を含む値は、自身を含めて近似的に等しいと判定しない。作用としての一般的な等価性判定は初期版の要件とせず、小規模系を数値表現へ変換した検証で代用できるものとする。

## 6. Boundary condition と小規模 chain

### 6.1. Open chain

長さ `N` の open chain の nearest-neighbor 関係は

```text
(0, 1), (1, 2), ..., (N - 2, N - 1)
```

とする。

### 6.2. Periodic chain

通常は open chain の関係に `(N - 1, 0)` を加える。ただし、nearest-neighbor 関係は 4.1 節の規約に従って正規化された、向きのない一意な site 対として扱い、self-loop および重複辺を生成しない。

小規模系の規約を次のように定める。

- `N = 0`: 空の chain を許可し、nearest-neighbor 関係は存在しない。
- `N = 1`: self-loop を生成せず、nearest-neighbor 関係は存在しない。
- `N = 2`: `(0, 1)` の一つだけを nearest-neighbor 関係とする。

したがって、periodic boundary condition は多重辺を持つ周期格子としてではなく、単純な無向グラフとして解釈する。

## 7. 行列化しない計算

演算子の状態への作用を、行列化しない計算の基本操作とする。

```math
|\phi\rangle = A|\psi\rangle
```

期待値は、この作用と内積を用いて

```math
\operatorname{expectation}(\psi,A)
= \langle\psi|A|\psi\rangle
= \langle\psi|\phi\rangle
```

と定義する。これは正規化されていない状態に対しても、式どおりの値を返す。

正規化期待値は別の操作として

```math
\operatorname{normalized\_expectation}(\psi,A)
= \frac{\langle\psi|A|\psi\rangle}
       {\langle\psi|\psi\rangle}
```

と定義する。zero-norm state に対する正規化期待値の計算はエラーとする。

初期版では、以下を満たすことを要件とする。

- 代数的な演算子の状態への作用に、完全な演算子行列を必要としない。
- 代数的な内積と期待値の計算に、密行列または疎行列への変換を必要としない。
- 積状態や局所演算子に対する局所的な因子分解は、上記の意味論を維持した最適化として追加できる。

一般の線形結合状態への作用では、状態を構成する項数が増大する可能性がある。あらゆる入力に対して Hilbert 空間の次元によらない計算量を保証することは、初期版の要件としない。

## 8. 数値表現への変換

初期版では以下の変換を提供する。

- 代数的な状態から密ベクトルへの変換
- 代数的な状態から疎ベクトルへの変換
- 代数的な演算子から密行列への変換
- 代数的な演算子から疎行列への変換

変換処理は代数層から分離し、代数オブジェクトの意味を変更してはならない。すべての数値表現は、5.2 節で定めた同一の基底順序を使用する。

数値表現とともに、使用した Hilbert 空間および基底順序を確認できる手段を設けることが望ましい。

tensor network への変換は future work とする。ただし、将来の MPO 変換が局所演算子積の有限和を利用できるように、代数情報を不要に失わない設計とする。

## 9. Non-functional Requirements

### 9.1. Performance: lazy simplification

代数的な展開、同類項の結合、および局所代数に基づく簡約は、常に lazy に行う。

- 状態または演算子の加算、scalar multiplication、通常の積、および tensor product の構築時に、式全体を eager に展開または簡約しない。
- 未評価の代数構造を保持し、状態への作用、内積、trace、数値表現への変換、または明示的な簡約結果の要求など、結果を必要とする時点で必要な部分だけを評価する。
- 評価に必要のない部分式を展開しない。
- 同一の部分式が繰り返し必要になる場合は、意味論と memory usage を損なわない範囲で評価結果を再利用できる設計とする。
- lazy evaluation の有無によって、代数オブジェクトの数学的な意味や計算結果が変化してはならない。

未簡約の source expression は評価後も保持する。簡約または展開した結果は、source expression から再計算できる破棄可能な cache として扱う。cache の有無によって数学的な結果が変化してはならない。どの単位で cache するか、および明示的な汎用 simplify 操作を提供するかは、具体的なプログラム設計で決定する。

式構造の表示と簡約結果の表示は別の操作として提供する。式構造の表示は、暗黙に式全体を展開または簡約してはならない。

### 9.2. Memory usage

状態、演算子、および Hilbert 空間の操作では、データのコピーを最小限に抑え、memory friendly な実装とする。

- 不変な Hilbert 空間情報、局所空間、基底情報、および代数式の部分構造を可能な限り共有する。
- 加算、積、scalar multiplication、および tensor product のたびに、変更されない operand 全体を複製しない。
- API の所有権設計では、参照、共有所有、move、view、または同等の仕組みを適切に用い、不要な deep copy を避ける。
- 密表現への変換を明示的に要求されるまで、Hilbert 空間の次元に比例する密な領域を確保しない。
- 疎表現への変換では、可能な限り密表現を中間生成せず、代数表現から直接構築する。
- 一時的な中間状態および中間演算子の peak memory usage を抑える。

具体的な内部表現は実装設計に委ねるが、structural sharing や同等の仕組みを採用できる構造が望ましい。性能向上のための cache は、無制限に保持して memory usage を悪化させてはならず、必要に応じて解放できなければならない。

### 9.3. Extensibility: replaceable local algebra

局所 Hilbert 空間が持つ代数構造は、適切な抽象化の下で差し替え可能とする。spin-1/2 固有の規則を、状態、演算子、または数値変換の共通処理へ直接埋め込んではならない。

局所代数の抽象化は、少なくとも以下の情報または操作を提供できるものとする。

- 局所 Hilbert 空間の次元
- 順序付けられた局所基底
- 局所状態および局所演算子の識別
- 恒等演算子
- 局所演算子の積
- 局所演算子の adjoint
- 局所演算子の局所状態への作用
- 局所状態の内積
- 局所演算子の trace
- 局所代数に固有の簡約規則
- 密および疎の局所数値表現を得るために必要な情報

spin-1/2 はこの抽象化の最初の公開実装とする。将来、spin-S、qudit、boson、または fermion を追加する際に、物理空間、複合 Hilbert 空間、代数式、および数値変換の共通機構を全面的に変更する必要がないことを目標とする。

新しい局所代数の追加は、既存の局所代数を列挙した closed な分岐を共通処理へ追加する方法ではなく、抽象化の新しい実装を提供する方法で行えるものとする。一つの系では単一種類の局所代数だけを使用し、異なる局所代数を持つ状態または演算子間の tensor product は初期版では許可しない。何を trait、type class、interface、または具象型として表現するかは、具体的なプログラム設計で決定する。

fermion の符号規則など、局所空間だけでは完結しない代数構造の具体的な対応は future work とする。

## 10. エラー条件

初期版では、少なくとも以下を検出してエラーとして扱う。

- 存在しない site への局所状態または演算子の配置
- 重複または不正な site 識別子
- 一つの積状態における、同一 site への矛盾した局所状態の割り当て
- 演算子式が参照する自由度を含まない状態または評価 context 上での作用、trace、または数値化
- 状態、演算子式、および評価 context の局所代数が一致しない場合の演算または評価
- 互換でない Hilbert 空間に属する状態同士の加算または内積
- 同じ自由度識別子を含む状態同士または演算子同士の tensor product
- 異なる種類の局所代数を使用する状態同士または演算子同士の tensor product
- `trace_out` が `ambient` の部分集合でない partial trace
- 基底順序を確定できない状態または演算子の数値表現への変換
- zero-norm state に対する正規化期待値の計算
- 単精度と倍精度を混在させた代数演算

二つの Hilbert 空間は、少なくとも自由度識別子、site 順序、および各 site に割り当てられた局所 Hilbert 空間が一致するときに互換とみなす。内部的な object identity だけを互換性の条件とはしない。

operand および構築時に与えられた情報だけで判定できるエラーは、原則として式の構築時に検出する。代表例として、互換でない Hilbert 空間に属する状態同士の加算、および同じ自由度識別子を含む演算子同士の tensor product がある。同一自由度上の通常の演算子積はエラーではなく、4.5.2 節の規約に従って許可する。

一方、lazy な式では、構築時点で未評価の部分式や、後から与えられる Hilbert 空間または評価 context に依存しなければ判定できない条件がある。このため、簡約時または評価時にもエラーを検出し得るものとする。同種の不正な式について、構築時に判定可能であれば早期にエラーを返し、構築時に判定できなければ簡約時または評価時にエラーを返す。

## 11. 初期版の対象外

以下は初期版の対象外とする。

- spin-1/2 以外の builtin algebra の提供。ただし、局所代数を差し替え可能であることを検証する非公開の test 用 algebra は除く
- boson などの無限次元局所 Hilbert 空間
- fermion と反交換関係
- 粒子数固定などの symmetry sector
- density operator および mixed state
- 熱状態
- 自動対角化
- 基底状態探索
- 時間発展
- 測定と確率 sampling
- 大規模な記号計算
- 記号係数
- 一般的な式変形 system
- tensor network への変換
- MPO の生成または最小化
- tensor network の truncation、gauge、および canonical form
- GPU および分散計算

## 12. 検証要件

### 12.1. 物理空間

- 長さ `N` の open chain の辺数が、`N > 0` では `N - 1` となる。
- 通常の長さの periodic chain が末端と先頭の関係を持つ。
- `N = 0`, `N = 1`, `N = 2` の periodic chain が 6.2 節の規約に従う。
- nearest-neighbor 関係に self-loop や重複辺が存在しない。
- 向きのない nearest-neighbor 関係が site ordering に従って `(min(i,j), max(i,j))` に正規化され、関係の列挙結果が同じ ordering に基づく辞書順となる。

### 12.2. 代数

- `Sx`, `Sy`, `Sz` が Hermitian である。
- `S+` と `S-` が互いに adjoint である。
- spin 演算子が採用した規約における交換関係を満たす。
- 演算子積が右側から作用する規約と一致する。
- 状態の加算、scalar multiplication、および tensor product が線形性と定められた因子順序を満たす。
- 演算子の加算、scalar multiplication、および tensor product が線形性と定められた因子順序を満たす。
- 異なる自由度を参照する演算子の和が、評価 context 上で適切に identity extension される。
- 通常の演算子積と演算子の tensor product が明確に区別される。
- 恒等演算子の trace が Hilbert 空間の次元と一致する。
- 局所 spin 演算子 `Sx`, `Sy`, `Sz` の trace が 0 となる。
- tensor product の trace が各因子の trace の積と一致する。
- 代数的に求めた trace が、密表現および疎表現から求めた trace と一致する。
- support に含まれない自由度を trace out したとき、その局所次元が係数へ掛けられる。
- 空集合、真部分集合、および全自由度に対する partial trace が定められた結果を返す。
- Neel 積状態の norm が 1 である。
- 構築した Heisenberg Hamiltonian が Hermitian である。
- 代数的な演算子作用および内積から得た期待値が解析解と一致する。
- operand だけで判定できる不正な状態和および不正な tensor product が構築時にエラーとなる。
- 構築時には判定できない不正な式が、必要な情報を得た簡約時または評価時にエラーとなる。

### 12.3. 数値変換

- 二 site spin-1/2 系の基底順序が 5.2 節の規約と一致する。
- 状態の密表現と疎表現が近似的に一致する。
- 演算子の密表現と疎表現が近似的に一致する。
- 密表現から求めた期待値が代数的に求めた期待値と一致する。
- 疎表現から求めた期待値が代数的に求めた期待値と一致する。
- 単精度および倍精度のそれぞれで、適切な許容誤差を用いた検証が成功する。

テストでは、検証に使用した `atol` と `rtol` を明示する。

Hermiticity と解析的なエネルギー期待値の検証には、実数の結合定数 `J` を使用する。

### 12.4. Non-functional Requirements

- 長い加算列や積の構築時に、式全体の展開または簡約が eager に実行されないことを確認する。
- 結果の要求時に、必要な部分式だけが評価されることを確認する。
- lazy な代数計算と、小規模系を完全に具体化した計算の結果が一致することを確認する。
- 式構造の表示が source expression を維持し、簡約を暗黙に実行しないことを確認する。
- cache を破棄して再評価しても同じ結果が得られることを確認する。
- 代表ユースケースにおいて、演算子の合成だけでは Hilbert 空間の次元に比例する密な領域が確保されないことを確認する。
- 疎表現への変換に密表現の中間生成が不要であることを確認する。
- 代数式の合成時に、変更されない大きな部分構造が不要に deep copy されないことを確認する。
- spin-1/2 の実装が局所代数の抽象化を通じて利用され、共通処理が spin-1/2 固有の分岐を必要としないことを確認する。
- core 部分の差し替え可能性を検証するため、公開対象ではない test 用の別の有限次元局所代数を同じ抽象化で実装し、別の系として状態への作用、trace、および数値変換の共通機構を利用できることを確認する。この test 用 algebra は、初期版が提供する builtin algebra には含めない。
- 異なる種類の局所代数を同一系で混在させる操作が拒否されることを確認する。

## 13. 初期版の完成条件

初期版は、以下をすべて満たした時点で完成とする。

> 有限一次元 spin-1/2 chain について、定められた小規模系の規約を含む open または periodic な物理空間を構築できる。その nearest-neighbor 関係から、一般的な局所演算子積の有限和として Heisenberg Hamiltonian を構築できる。状態間および演算子間の加算、scalar multiplication、および tensor product を行え、演算子の trace を行列化せずに計算できる。Neel 状態を積状態として構築し、完全な行列を生成せず、演算子の状態への代数的な作用を通じて期待値を評価できる。また、同じ状態と演算子を一貫した基底規約の下で密および疎の数値表現へ変換でき、代数計算を含む各経路の結果が指定した許容誤差内で一致することを検証できる。これらの処理は lazy simplification、不要な copy の抑制、および差し替え可能な局所代数の抽象化という Non-functional Requirements を満たす。

## 14. Future work

初期版の後に検討する主な拡張候補を以下に示す。

- spin-S および一般の qudit
- boson および局所 Hilbert 空間の truncation
- fermion と符号規則
- symmetry sector と保存量
- 任意 graph、座標、距離、plaquette、および hyperedge
- parameterized operator および記号係数
- 高度な簡約戦略と簡約結果の最適化
- tensor network、MPS、および MPO への変換
- 時間発展、固有値問題、および基底状態探索
- 大規模計算向け backend
