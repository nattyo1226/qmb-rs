# 数値 export の設計

更新日: 2026-10-27

## 1. 状態と責任分担

**設計採用済み・未実装（2026-10-27）**。[実装計画](../implementation-plan.md)の Step 6 に向けた推奨案をすべて採用し、一時メモから本書へ統合した。本書の型・API・feature・依存指定は実装する契約であり、製品コードや依存への反映、受入確認の完了を示すものではない。

数学的契約は[要件第6節](../requirements.md#6-基底export)、既存の基底列挙・座標取得 API は[基底の設計](basis.md)、初期版の完成条件は[初期版の設計](initial-release.md)に従う。

代数的な基底列挙・作用・座標取得は `qmb-core`、数値データの格納と backend への接続は新設する `qmb-export` が担当する。core の非公開の式木へアクセスせず、評価規則を再実装しない。core に backend の依存を加えず、既存の `no_std + alloc` を維持する。export は std を利用できる一つの crate にまとめる。

## 2. 出力形式と結果型

### 2.1. 数値コンテナ

| 対象 | backend に依存しない密出力 | 外部 crate の数値コンテナ |
| --- | --- | --- |
| 状態ベクトル | `Vec<C>` とメタデータ | `ndarray::Array1<C>`、`sprs::CsVec<C>` |
| 演算子行列 | 平坦な `Vec<C>` とメタデータ | `ndarray::Array2<C>`、`sprs::CsMat<C>` |

密ベクトル・密行列の数値データはどちらも `Vec<C>` とする。密行列の shape（行数・列数）と格納順序はメタデータに持たせ、密行列専用の数値コンテナ型を新設しない。利用者は Vec や slice とメタデータを使い、別の外部 crate へ接続できる。

Vec とメタデータへの出力、入力から ndarray への直接出力、取得済みの Vec 出力結果から ndarray への変換をすべて公開する。ndarray への直接出力も内部では Vec の出力・変換処理を共有する。疎表現は Coordinates から sprs を構築し、密出力を中間生成しない。

### 2.2. Exported とメタデータ

数値データ `D` と解釈情報 `M` を保持する共通の `Exported<D, M>` を提供し、具体的な組み合わせには型別名を提供する。これは `qmb-export` が新設する結果型である。

- `data()`・`metadata()` でそれぞれを借用する。
- `into_data()` で数値コンテナを、`into_parts()` で数値データとメタデータを所有値として取り出す。
- Vec の slice は `data().as_slice()` で取得する。
- データ長・shape の対応を維持するため、コンテナ全体を自由に伸縮・reshape できる可変参照は初期 API に含めない。加工する場合は所有値を取り出す。

メタデータはベクトル用の `VectorMetadata` と行列用の `MatrixMetadata` に分ける。ベクトルは長さ、行列は行数・列数・格納順序を持ち、`shape()` と、密行列では `layout()` を提供する。密表現では Vec の長さと shape の要素数が一致することを構築・変換時に確認する。

メタデータから、使用した始域・終域、基底、因子順序を確認できるようにする。基底は既存の Coordinates と同じく借用し、複合式の空間は `HilbertSpace` ハンドルを共有する。局所式では因子数・局所次元を記録する。結果は元の式を借用せず、基底の Clone や Index の Clone・Hash・Ord を追加要求しない。精度は係数型 `C` で表し、精度ラベルを重複保持しない。

行列は同じ状態基底を入力・出力に使い、選択した同一の全作用空間を始域・終域とする。状態型 `S` は式からの export と行列のメタデータに型として対応付け、保持のために S の実体・Clone・数値能力を要求しない。

Coordinates 自体は対象の primitive 型を持たない。その直接 export は基底と空間に対応する座標ベクトルとして扱い、元の状態型や演算子型を Coordinates だけから復元したとみなさない。演算子基底の座標ベクトルを状態へ作用する行列と解釈しない。

## 3. 公開 API と所有権

### 3.1. Export trait

「型 T へ export できる」という能力を、次の trait として提供する。

```rust
pub trait Export<T, Args = ()> {
    type Error;

    fn export(self, args: Args) -> Result<T, Self::Error>;
}
```

`T` はメタデータを含む結果型を基本とする。変換先は返却型から推論させるか、`Export::<Target, _>::export(source, args)` で明示する。引数も trait の型引数とし、同じ変換先に対して異なる基底型を扱えるようにする。能力の制約は各実装で要求する。

| 変換元 | 結果が持つ数値データ | Args |
| --- | --- | --- |
| `&LocalState`・`&State` | Vec、ndarray の Array1、sprs の疎ベクトル | `&basis` |
| 所有する `Coordinates`・`&Coordinates` | Vec、ndarray の Array1、sprs の疎ベクトル | `()` |
| `&LocalOperator`・`&Operator` | Vec、ndarray の Array2、sprs の疎行列 | `MatrixArgs` |
| 所有する Vec 出力結果（ベクトル用メタデータ） | ndarray の Array1 | `()` |
| 所有する Vec 出力結果（行列用メタデータ） | ndarray の Array2 | `()` |

式への実装は借用型に置き、入力を保持する。所有する Coordinates からは係数を移動し、借用する Coordinates からは必要な係数だけを Clone する。既存の `into_entries()` と変換前に取得した基底・空間情報を利用する。

Vec 出力結果から ndarray への変換は所有値を受け取り、同じ係数型のまま数値データを移動する。基底・空間・shape の解釈も引き継ぐ。ベクトル用と行列用のメタデータによって Array1／Array2 への対応を定める。

以下は採用した API の利用イメージであり、実装・コンパイル確認済みの例ではない。

```rust
// Vec とメタデータへの出力。
let dense: Exported<Vec<_>, _> = (&psi).export(&basis)?;
let shape = dense.metadata().shape();
let values = dense.data().as_slice();

// 取得済みの Vec 出力結果を消費して ndarray へ変換する。
let array: Exported<ndarray::Array1<_>, _> = dense.export(())?;

// 内部で同じ Vec 出力・変換処理を使う直接出力。
let direct: Exported<ndarray::Array1<_>, _> = (&psi).export(&basis)?;
```

### 3.2. 演算子の引数と空間選択

`Operator<I, O, C>` 自体が持たない状態型 S を、`MatrixArgs::<S, _>::new(&basis)` で明示する。明示空間が必要な場合は `.in_space(&ambient)` を付ける。

設定型は `MatrixArgs<'b, S, B, A = ()>` とし、`&'b B`、空間指定 A、S を表す `PhantomData<fn() -> S>` を持つ。`.in_space` は `A = &'h HilbertSpace<I>` の引数型を返す。空間指定がない場合は A を `()` とし、不要な Index 型の指定・推論を要求しない。局所演算子の実装は `A = ()` の引数を受け付ける。

演算子の返却型は `Exported<Vec<_>, _>`、`Exported<ndarray::Array2<_>, _>`、または疎行列を持つ結果型として選択する。複合演算子は明示空間を優先し、省略時は既存の束縛を使う。未束縛かつ明示空間なしはエラーとし、局所配置から全作用空間を推論しない。局所演算子は局所状態基底を使い、大域空間の指定を要求しない。

空間の選択と参照の検証は出力領域の確保・列生成より先に行い、公開 API の `bind`・`act`・`act_in` へ委譲する。zero を理由に不正な参照を見逃さない。基底列挙では個数・番号の整合性を検査するが、完全性・線型独立性を一般の利用者実装について自動証明することは要求しない。

## 4. 変換経路と数値的契約

### 4.1. 密表現

状態は `coordinates_in` で得た疎な座標を Vec に配置する。行列は入力基底状態を一つずつ列挙し、作用結果を状態基底で座標化して列を得る。この列生成を密・疎で共有し、全基底状態・全列の Coordinates を一括保持しない。演算子基底への全展開や、密な一列の中間生成は必須にしない。

密行列は行優先に固定し、成分 `(row, col)` の位置は `row * ncols + col` とする。格納順序は仕様とメタデータで明示する。因子の基底番号は既存の右端が最速の mixed-radix とし、行列のメモリ配置とは区別する。

ndarray への変換では、ベクトルは `Array1::from_vec`、行列はメタデータの shape と Vec を `Array2::from_shape_vec` に渡す。追加の数値データのコピーを避け、結果には同じ解釈のメタデータを保持する。参照: [ndarray 0.17.2 の配置・構築 API](https://docs.rs/ndarray/0.17.2/ndarray/struct.ArrayBase.html)。

### 4.2. 疎表現と sprs の構築 API

疎ベクトルは Coordinates の成分から index・値の配列を作り、`CsVec::try_new(dimension, indices, data)` で構築する。疎行列の初期 export は CSC・index は `usize` に固定する。列ごとに得た座標の平坦番号順を使い、重複のない成分から `indptr`・`indices`・`data` を作って `CsMat::try_new_csc` へ渡す。検査付き API を利用し、unsafe や非公開の検査省略経路を使わない。

sprs に qmb の Coordinates をそのまま渡す既存 API があるという意味ではない。`qmb-export` の Export 実装が、`flat_entries()`・`into_entries()` などと sprs の constructor を接続する。密表現用の Vec 中間層や、新しい独自の疎コンテナを必須経路にしない。

sprs には座標形式で行列を作る `TriMat::from_triplets`・`add_triplet` と、`to_csc`・`to_csr` もある。採用した列生成では整列・集約済み座標から CSC を構築し、TriMat の経由を必須にしない。CSR が必要な利用者は取り出した sprs 行列を変換する。export 時の形式選択や別の index 整数型は後続の必要性に応じて追加する。

参照: [sprs の疎ベクトル構築 API](https://docs.rs/sprs/0.11.4/sprs/struct.CsVecBase.html)、[CSC 構築・形式変換 API](https://docs.rs/sprs/0.11.4/sprs/struct.CsMatBase.html)、[座標形式の構築 API](https://docs.rs/sprs/0.11.4/sprs/struct.TriMatBase.html)。

疎行列には各軸の次元・列数に応じたポインタ列・非零成分数の表現可能性を要求し、密行列の要素数 `nrows * ncols` の整数化を必須にしない。CSC の列ポインタと入力基底の走査は必要であり、メモリ量・計算量が常に非零成分数だけで決まるとは保証しない。

### 4.3. 係数・精度・非正規直交基底

係数は C のまま扱い、標準の `Complex<f32>`・`Complex<f64>` を検証する。能力は必要な実装にだけ要求し、Float・内積・演算子 Product などを一律に要求しない。

- 所有する Coordinates → 疎ベクトルでは、係数の Clone・Zero を追加要求しない。
- Coordinates → 密ベクトルでは、未格納成分を埋めるため Zero を要求する。借用する Coordinates では必要な係数の Clone も要求する。
- 式 → 座標取得は既存 API の制約を引き継ぐ。
- Vec 出力結果 → ndarray の所有値変換には係数の算術を追加要求しない。

厳密な零のみを省略する既存の座標の契約を維持し、近似誤差による切り捨てや非有限値の独自処理は追加しない。各 API は同じ C を維持し、精度変更は数値データを取り出した後の明示的な成分変換で扱う。

非正規直交基底でも、行列の列は作用結果の展開座標から得る。座標ベクトルのユークリッド内積を Hilbert 空間の内積と取り違えない。数値期待値の検証には、座標 c、Gram 行列 G、作用行列 A に対して `c† G A c / (c† G c)` を用いるか、独立な正規直交基底に再構成する。Gram 行列生成の新しい製品 API は実装の前提にしない。

## 5. エラーとメモリ確保

`ExportError` を提供し、既存の基底・空間・次元エラーを source として保持する。追加で、要素数・バイト数・圧縮配列長の overflow、export 自身が確保する領域の確保失敗、backend の shape・圧縮構造の不正を区別する。backend の不正は検査付き constructor の結果から伝える。

密領域では `nrows * ncols` と格納可能サイズを確保前に確認し、疎領域では各配列に必要な長さを確認する。主要な Vec 確保には `try_reserve` 系を用いる。既存の core 評価器・backend 内部を含むすべてのメモリ不足から回復する保証は設けない。

空間構築・代数的評価の契約は変更せず、数値コンテナが必要とするサイズ検証を export の境界で行う。総次元が overflow しても、既存の疎な座標保持・基底展開を利用できる契約は維持する。

## 6. feature・再公開・依存

### 6.1. 有効化と利用窓口

`qmb-export` の既定 feature は空とし、Export・Vec 出力・メタデータを常時提供する。`ndarray`・`sprs` の feature で各 API と任意依存を有効にする。ndarray feature は式・Coordinates からの直接出力と、所有する Vec 出力結果からの変換の両方を有効にする。

`qmb` の任意の `export` feature は両 backend を有効にし、`qmb::export` として crate を再公開する。Export は root と prelude からも提供する。Vec のみ、または片方の API を使う場合は `qmb-export` の feature を直接選択する。backend の型は対応 feature の下で `qmb_export::ndarray`・`qmb_export::sprs` として再公開する。

sprs は通常依存に ndarray を含むため、sprs feature だけでも ndarray 自体は依存に入る。この分離は公開する変換 API の選択である。sprs の既定の alga・multi_thread は無効とし、ndarray は std の範囲で使う。BLAS・rayon・保存用 feature は初期導入に含めない。

### 6.2. 採用する依存指定と検証状態

| 依存 | 採用する指定 |
| --- | --- |
| sprs | `version = "=0.11.4", default-features = false` |
| ndarray | `version = "0.17.2"`（std のみ） |

MSRV 1.85 を維持する。sprs 0.11.4 の公開 manifest は Rust 1.85.0 と ndarray `>=0.15.0, <0.18`、ndarray 0.17.2 は Rust 1.64 を指定する。sprs 0.11.5 は Rust 1.88.0 を指定し、README の Rust 1.64 という説明と一致しないため、manifest を基準にして 0.11.4 を完全一致で指定する。

参照: [sprs 0.11.4 manifest](https://docs.rs/crate/sprs/0.11.4/source/Cargo.toml)、[sprs 0.11.5 manifest](https://docs.rs/crate/sprs/0.11.5/source/Cargo.toml.orig)、[ndarray 0.17.2 manifest](https://github.com/rust-ndarray/ndarray/blob/0.17.2/Cargo.toml)、[sprs の feature](https://docs.rs/crate/sprs/0.11.4/features)。

依存指定の採用と実ビルドの検証は区別する。Cargo.toml・Cargo.lock への追加、sprs と直接依存での ndarray のバージョン統一、推移依存を含む MSRV 1.85 のビルド・テストは未実施であり、Step 6 の実装時に確認する。

## 7. 受入確認

以下は実施する検証であり、成功済みの結果ではない。

- 独立に定義した spin の局所行列・二因子の 4×4 行列を用い、複素係数を含む状態・積・作用・随伴・trace を代数計算・座標・密・疎の各経路で比較する。演算子基底で展開した式からの経路とも比較する。
- 状態・演算子それぞれで、ndarray への直接出力と Vec 出力結果からの変換について、成分・shape・格納順序・基底・空間・因子順序・係数精度の一致を確認する。非対称かつ虚部のある行列で flatten・転置の取り違えを検出する。
- Vec の slice とメタデータだけで独立な参照計算を行えること、Vec → ndarray の所有値変換で数値領域のポインタが維持されることを確認する。
- 三次元代数・非正規直交基底の座標と再構成、Gram 行列による内積・数値期待値を確認する。spin 固有の分岐を共通変換へ追加しない。
- 空空間・単因子・入れ子・zero・identity・spectator、明示空間による束縛の上書き、未束縛・不正な参照・基底・overflow を確認する。0因子は既存の単位空間の規約に従う。
- Heisenberg／Neel は長さ0〜4の open・periodic を基本とし、解析的期待値、独立な行列、密・疎による計算を比較する。実行可能な縦断例を完成させる。
- 比較式は要件第5.6節に従い、f32 の atol・rtol は初期値を各 `2e-5`、f64 は各 `1e-12` とする。調整する場合は計算規模と誤差の根拠を記録する。
- 暗黙の精度変更の拒否と、取り出した数値データの明示的な精度変換を検証する。保持・移動に不要な係数・Index・primitive の制約を追加していないことも確認する。
- 疎ベクトルは、非零成分の少ない大きな空間で密領域を確保しないことを確認する。疎行列は列ごとに処理し、密変換・全列の一括保持を行わないことを対象テストと実装レビューで確認する。
- backend 無効・ndarray のみ・sprs のみ・両方、および qmb の export 有効・無効を検証する。関連テスト、workspace の test・format・Clippy・警告をエラーとするドキュメント生成、core 単独チェック、MSRV 1.85 での検証を行う。

## 8. 実装時に具体化する内部事項と延期範囲

本書の公開契約・依存指定・受入条件は採用済みとする。メタデータ型の完全な型引数一覧、具体的な型別名、エラー variant の名称、非公開の補助型・module 分割は、この契約を維持して実装時に具体化する。

Partial trace、汎用的な写像の表現行列・表現係数、異なる入出力基底を扱う汎用機構、Graded・異種局所代数・追加 backend は既存の初期版の範囲に従って延期する。通常の状態・演算子の数値 export は、これらを前提条件にしない。
