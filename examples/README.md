# qmb-rs examples

qmb-rs の利用例を、workspace に属する Cargo パッケージの bin としてまとめる。各利用例は `qmb-core`、`qmb-basic`、`qmb-geometry`、`qmb-spin-half` などを直接利用し、依存設定は workspace から継承する。

リポジトリ直下から実行する。

```sh
cargo run -p examples --bin two_spins
cargo run -p examples --bin tensor_products
cargo run -p examples --bin basis_expansion
cargo run -p examples --bin heisenberg
cargo run -p examples --bin heisenberg --features export
cargo run -p examples --bin hubbard
cargo run -p examples --bin hubbard --features export
cargo run -p examples --bin spin_one
cargo run -p examples --bin spin_one --features export
cargo run -p examples --bin geometry_relations
cargo run -p examples --bin pauli_strings
```

`two_spins` は `Spin::product_state`・`Spin::operator_product_at`・`Spin::try_operator_sum` から2スピンの一重項・三重項と交換 Hamiltonian `H = J_a · J_b = (S_a · S_b)/4` を手で組み立てる。norm squared、エネルギー期待値（一重項 −0.75、三重項 0.25）、固有状態の残差、状態間の内積を計算・検証する。追加 feature は不要で、`--no-default-features` でも実行できる。

`hubbard` は `Fermion<Complex<f64>>` と、各 Site が `(up,down)` の2 mode を持つ graded 空間を使う。2 Site の標準 Hubbard・symmetric Hubbard を `t=1,U=4` で構築し、`|up,down>` の期待値がそれぞれ0・−2になることを確認する。parity 演算子 `Mz=I−2n` の基底座標も取得する。既定の normalize feature では正規化・正規化期待値、任意の export feature では16×16の Vec・sprs 行列の一致を確認する。`--no-default-features` でも実行でき、Index は利用側で発行する。

`spin_one` は `SpinOne<Complex<f64>>` の三状態と四重極を使い、単イオン異方性と横磁場を持つ局所 Hamiltonian を構築する。`Sp²=√(3/2) Q2p` の ladder 座標、Cartesian / ladder の往復変換、二因子での作用を確認する。2 Site の XYZ・Heisenberg・BLBQ を関連関数から構築し、`J=1,K=0.25` の BLBQ では `|+,->` の期待値が0.5となることも確認する。bilinear 結合の正符号は強磁性的で、2 Site の周期辺は一度だけ数える。既定の normalize feature では正規化と正規化期待値、export feature では局所3×3と BLBQ 9×9の Vec・sprs 行列の一致を確認する。`--no-default-features` でも実行できる。

`geometry_relations` は `SiteMapper` と `Relation<G>` を使い、同じ Chain から最近接と連続3 Site を列挙する。サンプル内の正方格子では、独自の Geometry・Relation と、中心・方向を保持する Group の実装方法を示す。複数葉の Site に対する `leaves_cartesian()`・`leaves_zip()`・`combine_leaves()` も扱う。正方格子はこの使用例内の実装である。

`pauli_strings` は1 Site が1 spin-1/2 の場合に、恒等項を含む高々2 Site の Pauli string を列挙する。3 Site で37項を構築し、`Sx/Sy/Sz` を Pauli 行列としてそのまま使う。支持集合の全組み合わせは通常の Iterator で列挙し、geometry の Relation とは分ける。構築に係数の評価能力や Index の Clone を要求しない例にもなっている。

この二つも `--no-default-features` で実行できる。サンプル内のテストは次のコマンドで実行する。Pauli string の行列比較テスト2件は `export` feature を使う。

```sh
cargo test -p examples --no-default-features --bin geometry_relations --bin pauli_strings
cargo test -p examples --features export --bin geometry_relations --bin pauli_strings
```

`normalize` はデフォルトで有効、`export` は任意の feature である。次のように正規化を無効にした構成も実行できる。

```sh
cargo run -p examples --bin heisenberg --no-default-features
cargo run -p examples --bin heisenberg --no-default-features --features export
```

`basis_expansion` は `type Spin = SpinHalf<Complex<f64>>` で係数型を一度選び、`ExpressionFactory` の型エイリアスから `Spin::state`・`Spin::operator` で局所式を作る。`heisenberg` は Index を利用側で発行して Site 順の局所空間列を作り、`SpinHalfModels<Spin>` の `heisenberg`・`neel` から構築した式を正規化・密疎 export まで使う。単精度へ切り替える場合は Spin の定義を `Complex<f32>` に変更する。

戻り値型には `qmb_core::StateOf<Spin,I>`・`OperatorOf<Spin,I>`、局所式には `LocalStateOf<Spin>`・`LocalOperatorOf<Spin>` を使える。行列 export の `MatrixArgs::<SpinHalfState,_>` は作用する状態 primitive 型の指定であり、係数精度の再指定ではない。出力先コンテナは引き続き `Exported<Vec<_>,_>` などで選ぶ。

移行前の `Chain<I>`・`ToSpace`・`qmb_spin_half::utility` の自由関数は撤去した。現在の Chain は `Chain::new(n, boundary)` として構築し、配置は外部の `&[HilbertSpace<I>]` を渡す。

式の構築・集約は `Spin::state_sum`・`Spin::tensor_in`・`Spin::operator_sum`・`Spin::operator_product` に揃えられる。状態は既存空間を明示し、`product_state(&space, primitives)` は葉の順序に primitive を配置、`tensor_in(&space, states)` は配置済み状態を合成する。演算子の `_at` は `(&factor, primitive)` の列を受け取る。集約の `try_` 版は Result の列を受け取り、最初の入力 ExpressionError を優先し、その後の要素を取り出さない。全成功時に空間・束縛を検証する。

`scaled` は `impl Into<C>` を受け取るため、選択済みの複素係数型に `.scaled(0.5)` と実数を直接渡せる。旧 scaled の期待型に依存した `.scaled(value.into())` は `.scaled(value)` とし、Default の生成型や constructor の C を推論できない場合は明示する。`*` は従来の係数型を直接受け取る契約を維持する。

`tensor_products` は `TensorKind` を空間の構築時に指定し、共通の作用・trace API を使う例である。Spin-half は全 primitive が even、P=I なので、ordinary / graded で同じ数値結果になる。共通機構の契約は[ordinary / graded の共通評価](../docs/design/composite-expression.md#12-ordinary--graded-の共通評価と移行)を参照する。
