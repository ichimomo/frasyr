# frasyr
- Fisheries Research Agency (FRA) provides the method for calculating sustainable yield (SY) with R
- 水研機構で開発した，MSYを基礎とした目標管理基準値を計算するためのRのパッケージです．開発途中のものであること，ご承知おきください．

# 使い方

```
# install.pakcages("devtools") # <-- devtoolsをインストールしていない人
devtools::install_github("ichimomo/frasyr") # frasyrのインストール
library(frasyr) # frasyrの呼び出し
```

# マニュアル
- VPAによる資源量推定　https://ichimomo.github.io/frasyr/doc/vpa.html
- 新ルールのもとでの将来予測計算 https://ichimomo.github.io/frasyr/doc/future.html
- 今後追加予定



# 開発ワークフロー（開発者向け）
1. 本レポジトリを直接cloneするか，自分のところにfolkしてからcloneする
2. 自分の環境下でコードを修正
```{r}
3. devtools::load_all() # パッケージをビルド，テスト，インストール
```
4. 2,3を繰り返してコードの修正を完了する

5. テストコードを走らせる
```{r}
devtoos::test()
```

6. 仕上げ（vignetteも作り直す）
```{r}
devtoos::check()
```
7. 修正したファイルを git add ファイル名，git commit -m "コメント" 
8 folkした場合は，pull request

