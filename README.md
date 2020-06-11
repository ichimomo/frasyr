# frasyr
  [![Travis build status](https://travis-ci.com/ichimomo/frasyr.svg?branch=master)](https://travis-ci.com/ichimomo/frasyr)
  [![Codecov test coverage](https://codecov.io/gh/ichimomo/frasyr/branch/master/graph/badge.svg)](https://codecov.io/gh/ichimomo/frasyr?branch=master)
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
- v1.20 : 2019年アジ・イワシ事前検討会用
- v2.00 : 2020年アジ・イワシ研究機関会議用(事前配布版)
   - 将来予測のメイン関数をfuture.vpaからfuture_vpaに移行。基本的な使い方は以下のマニュアルのリンク先を参照のこと
- v2.01 : 2020年アジ・イワシ研究機関会議用(最終版)
   - fit.parでL1の場合のSDをRMSEに統一
- v2.10 : 6月上旬ごろアップ予定   

# マニュアル
- VPAによる資源量推定　[vignette](https://ichimomo.github.io/frasyr/doc/vpa.html)
- VPAモデル診断スクリプト　[vignette](https://ichimomo.github.io/frasyr/doc/Diagnostics-for-VPA.html)
- 再生産関係のモデル診断スクリプト [wiki](https://ichimomo.github.io/frasyr/doc/SRR-guidline.html)
- 将来予測関数の使い方：[wiki](https://github.com/ichimomo/frasyr/wiki/future_new)


<!--
- 新ルールのもとでの将来予測計算 https://ichimomo.github.io/frasyr/doc/future.html
- 管理基準値の計算 https://ichimomo.github.io/frasyr/doc/estMSY.html
-->



=======

frasyr_tool群の全体説明は[こちら](https://ichimomo.github.io/main/)
