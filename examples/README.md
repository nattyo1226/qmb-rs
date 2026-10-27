# qmb-rs examples

qmb-rs の利用例を、workspace に属する Cargo パッケージの bin としてまとめる。各利用例は `qmb-core`、`qmb-basic`、`qmb-geometry`、`qmb-spin-half` などを直接利用し、依存設定は workspace から継承する。

リポジトリ直下から実行する。

```sh
cargo run -p examples --bin basis_expansion
cargo run -p examples --bin heisenberg
cargo run -p examples --bin heisenberg --features export
```

`normalize` はデフォルトで有効、`export` は任意の feature である。次のように正規化を無効にした構成も実行できる。

```sh
cargo run -p examples --bin heisenberg --no-default-features
cargo run -p examples --bin heisenberg --no-default-features --features export
```

`examples` ディレクトリ内からは `-p examples` を省略できる。必要に応じて、リポジトリの direnv 環境も利用できる。

```sh
direnv exec . cargo run -p examples --bin heisenberg
```
