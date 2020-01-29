# frasyr
  [![Travis build status](https://travis-ci.com/ichimomo/frasyr.svg?branch=master)](https://travis-ci.com/ichimomo/frasyr)
- Fisheries Research Agency (FRA) provides the method for calculating sustainable yield (SY) with R
- 水研機構で開発した，MSYを基礎とした目標管理基準値を計算するためのRのパッケージです．開発途中のものであること，ご承知おきください．

# 使い方など
https://ichimomo.github.io/main/ に一括した情報へのリンクがあります

# インストール方法

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

--- 

frasyr_tool群の全体説明は[こちら](https://ichimomo.github.io/main/)



