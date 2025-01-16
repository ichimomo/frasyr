# frasyr
[![Build status](https://github.com/ichimomo/frasyr/actions/workflows/check-standard.yml/badge.svg)](https://github.com/ichimomo/frasyr/actions/workflows/check-standard.yml)
  [![Codecov test coverage](https://codecov.io/gh/ichimomo/frasyr/branch/dev/graph/badge.svg)](https://codecov.io/gh/ichimomo/frasyr?branch=dev)
- Fisheries Research Agency (FRA) provides the method for calculating sustainable yield (SY) with R
- VPAを用いた資源量推定と，その推定結果をもとにしてMSYを基礎とした目標管理基準値を計算するためのRのパッケージです．開発途中のものであること，ご承知おきください．

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
- v2.1.0.0 : 2020年度資源評価会議用プロトタイプバージョン
   - VPAを用いたモデル診断スクリプトをvignetteに追加
   - 途中でプラスグループが変わる（対馬マイワシ）VPAの計算を修正
   - ほか、関数のテストを充実。

# マニュアル
- VPAによる資源量推定　[vignette](https://ichimomo.github.io/frasyr/articles/vpa.html)
- VPAモデル診断スクリプト　[vignette](https://ichimomo.github.io/frasyr/articles/Diagnostics-for-VPA.html)
- fitSR関数による再生産関係推定　[vignette](https://ichimomo.github.io/frasyr/articles/fittingSR.html)

これらはRコマンドで以下のようにしても見れます。
```
devtools::install_github("ichimomo/frasyr", ref="dev", build_vignettes=TRUE) # インストールするときにvignetteを作る（時間かかります）
library(frasyr)
vignette(package="frasyr") # 利用可能なvignetteを調べる
vignette("vpa",package="frasyr") # VPAの実施のしかた
vignette("Diagnostics-for-VPA",package="frasyr") # VPAのモデル診断
```

- 再生産関係のモデル診断 [wiki](https://github.com/ichimomo/frasyr/wiki/Diagnostics-for-Stock-Recruitment-Relationships)
- 将来予測関数の使い方：[wiki](https://github.com/ichimomo/frasyr/wiki/future_new)


<!--
- 新ルールのもとでの将来予測計算 https://ichimomo.github.io/frasyr/doc/future.html
- 管理基準値の計算 https://ichimomo.github.io/frasyr/doc/estMSY.html
-->



=======

frasyr_tool群の全体説明は[こちら](https://ichimomo.github.io/main/)
