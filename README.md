# frasyr
- Fisheries Research Agency (FRA) provides the method for calculating sustainable yield (SY) with R
- 水研機構で開発した，MSYを基礎とした目標管理基準値を計算するためのRのパッケージです．開発途中のものであること，ご承知おきください．

# 使い方

```
# install.pakcages("devtools") # <-- devtoolsをインストールしていない人はインストールする
devtools::install_github("ichimomo/frasyr") # frasyrのインストール
devtools::install_github("ichimomo/frasyr",rev="dev") # 開発中バージョンのインストール
library(frasyr) # frasyrの呼び出し
```

# マニュアル
- VPAによる資源量推定　https://ichimomo.github.io/frasyr/doc/vpa.html
- 新ルールのもとでの将来予測計算 https://ichimomo.github.io/frasyr/doc/future.html
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
devtoos::test()
```
6. 仕上げ（vignetteも作り直す）
```{r}
devtoos::check()
```
7. 修正した一連の変更を 個人の名前のブランチichimomoに push する
8. Githubのウェブ上で，ichimomoからdevにpull request => merge
9. 安定版が完成した段階でmasterも変更

