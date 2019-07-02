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



# 開発ワークフロー（開発者向け）
- 本レポジトリをcloneするかfolkする
- 自分の環境下でコードを修正
```{r}
devtools::loac_all() # パッケージをビルド，テスト，インストール
```
- 修正したファイルを git add ファイル名，git commit -m "コメント" 
- folkした場合は，pull request
- などなど
