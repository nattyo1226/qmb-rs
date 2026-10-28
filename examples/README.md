# qmb-rs examples

qmb-rs の利用例を、workspace に属する Cargo パッケージの bin としてまとめる。各利用例は `qmb-core`、`qmb-basic`、`qmb-geometry`、`qmb-spin-half` などを直接利用し、依存設定は workspace から継承する。

リポジトリ直下から実行する。

```sh
cargo run -p examples --bin two_spins
cargo run -p examples --bin basis_expansion
cargo run -p examples --bin heisenberg
cargo run -p examples --bin heisenberg --features export
```

`two_spins` は `Spin::state_tensor_at`・`Spin::operator_product_at`・`Spin::try_operator_sum` から2スピンの一重項・三重項と交換 Hamiltonian `H = S_a · S_b` を手で組み立てる。norm squared、エネルギー期待値（一重項 −0.75、三重項 0.25）、固有状態の残差、状態間の内積を計算・検証する。追加 feature は不要で、`--no-default-features` でも実行できる。

`normalize` はデフォルトで有効、`export` は任意の feature である。次のように正規化を無効にした構成も実行できる。

```sh
cargo run -p examples --bin heisenberg --no-default-features
cargo run -p examples --bin heisenberg --no-default-features --features export
```

`basis_expansion` は `type Spin = SpinHalf<Complex<f64>>` で係数型を一度選び、`QuantumSystem` の `Spin::state`・`Spin::operator` で局所式を作る。`heisenberg` は Index を利用側で発行して Site 順の局所空間列を作り、`Spin::heisenberg`・`Spin::neel` から構築した式を正規化・密疎 export まで使う。単精度へ切り替える場合は Spin の定義を `Complex<f32>` に変更する。

戻り値型には `qmb_core::StateOf<Spin,I>`・`OperatorOf<Spin,I>`、局所式には `LocalStateOf<Spin>`・`LocalOperatorOf<Spin>` を使える。行列 export の `MatrixArgs::<SpinHalfState,_>` は作用する状態 primitive 型の指定であり、係数精度の再指定ではない。出力先コンテナは引き続き `Exported<Vec<_>,_>` などで選ぶ。

移行前の `Chain<I>`・`ToSpace`・`qmb_spin_half::utility` の自由関数は撤去した。現在の Chain は `Chain::new(n, boundary)` として構築し、配置は外部の `&[HilbertSpace<I>]` を渡す。

式の構築・集約は `Spin::state_sum`・`Spin::state_tensor`・`Spin::operator_sum`・`Spin::operator_product` に揃えられる。`_at` は `(&factor, primitive)` の列を受け取って配置も行い、裸の基本名は構築済みの式を受け取る。集約の `try_` 版は Result の列を受け取り、最初の入力 ExpressionError を優先し、その後の要素を取り出さない。全成功時に空間・束縛を検証する。

`scaled` は `impl Into<C>` を受け取るため、選択済みの複素係数型に `.scaled(0.5)` と実数を直接渡せる。旧 scaled の期待型に依存した `.scaled(value.into())` は `.scaled(value)` とし、Default の生成型や constructor の C を推論できない場合は明示する。`*` は従来の係数型を直接受け取る契約を維持する。
