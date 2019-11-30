# frasyr
  [![Travis build status](https://travis-ci.com/ichimomo/frasyr.svg?branch=master)](https://travis-ci.com/ichimomo/frasyr)
- Fisheries Research Agency (FRA) provides the method for calculating sustainable yield (SY) with R
- 水研機構で開発した，MSYを基礎とした目標管理基準値を計算するためのRのパッケージです．開発途中のものであること，ご承知おきください．

# 使い方

```
# devtoolsをインストールしていない人はインストールする
install.packages("devtools") 

# マスター版（最新・安定版）をインストールする場合
devtools::install_github("ichimomo/frasyr") 

# 開発中の最新版をインストールする場合（バグ可能性あり！）
# ref=""で開発中のブランチを指定します。だいたい、"dev"ブランチに開発中のものがあります
devtools::install_github("ichimomo/frasyr", ref="dev") 

# 過去の安定版を指定してインストールする場合
# @以下にリリースバージョンを指定します
devtools::install_github("ichimomo/frasyr@v1.00")

# 以上の操作をしてfrasyrをインストールしてから、以下のコマンドで呼び出します
library(frasyr) 

```

# リリースバージョン
- v1.00 : future-rvpaから移動してきたほぼそのままのバージョン
- v1.10 : future.vpaにuse.MSEオプションを追加

# マニュアル
- VPAによる資源量推定　https://ichimomo.github.io/frasyr/doc/vpa.html
- 新ルールのもとでの将来予測計算 https://ichimomo.github.io/frasyr/doc/future.html
- 管理基準値の計算 https://ichimomo.github.io/frasyr/doc/estMSY.html
- 再生産関係決定ガイドライン https://ichimomo.github.io/frasyr/doc/SRR-guidline.html
- 今後追加予定



# 開発ワークフロー（開発者向け）
## ブランチ構成
- master: 公開用のブランチ．いつでも動くようにしておく
- dev：開発用ブランチ．新規変更はこちらのブランチにアップすること
- dev - 個人の名前（たとえばichimomo）：個人の作業用ブランチ．このブランチで作業して，ある程度たまったらdevに変更をアップ
- bug_fix: masterにすぐに反映させたい細かいバグの修正
- web_edit: githubのウェブ上でファイルを編集する用．たとえばreadmeなど．編集が終わったらmaster?にmergeする・
## ワークフロー
1. 本レポジトリを直接cloneするか，自分のところにfolkしてからcloneする
2. 自分の環境下でコードを修正
   - インポートするパッケージの追加などは，自分の関数の近くに @import パッケージ名，または@importFrom パッケージ名 関数名として定義しておくと，checkのときにroxygen2が自動的にNAMESPACEを置き換えてくれるので，NAMESPACEは直接いじらない
3. パッケージをビルド，テスト，インストール
```{r}
devtools::load_all() 
```
4. 2,3を繰り返してコードの修正を完了させる
5. テストコードを走らせて確認
```{r}
devtools::test()
```
6. 仕上げ（vignetteも作り直す）
```{r}
devtools::check()
```
7. 修正した一連の変更を 個人の名前のブランチichimomoに push する
8. Githubのウェブ上で，ichimomoからdevにpull request => merge
9. 安定版が完成した段階でmasterも変更

