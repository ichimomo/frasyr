# frasyr
[![Build status](https://github.com/ichimomo/frasyr/actions/workflows/check-standard.yaml/badge.svg)](https://github.com/ichimomo/frasyr/actions/workflows/check-standard.yaml)
  [![Codecov test coverage](https://codecov.io/gh/ichimomo/frasyr/branch/dev/graph/badge.svg)](https://codecov.io/gh/ichimomo/frasyr?branch=dev)
- Fisheries Research Agency (FRA) provides the method for calculating sustainable yield (SY) with R
- VPAを用いた資源量推定と，その推定結果をもとにしてMSYを基礎とした目標管理基準値を計算するためのRのパッケージです．開発途中のものであること，ご承知おきください．

# 使い方など
https://ichimomo.github.io/main/ に一括した情報へのリンクがあります

# インストール方法

pak（推奨）または devtools でインストールできます。

```
# --- pak を使う場合（推奨） ---
# pakをインストールしていない人はインストールする
install.packages("pak")

# 開発中の最新版（devブランチ）をインストールする
pak::pkg_install("github::ichimomo/frasyr@dev")

# 特定の年・ブランチやリリースを指定する場合は @ 以下を変える
# pak::pkg_install("github::ichimomo/frasyr@dev2026")  # 2026年資源評価用ブランチ
# pak::pkg_install("github::ichimomo/frasyr@v1.00")    # 過去の安定版

# --- devtools を使う場合 ---
# install.packages("devtools")
# devtools::install_github("ichimomo/frasyr", ref="dev")
# devtools::install_github("ichimomo/frasyr@v1.00")    # 過去の安定版

# インストール後、以下のコマンドで呼び出します
library(frasyr)

```

# マニュアル
- パッケージ解説サイト（pkgdown）: https://ichimomo.github.io/frasyr/
- 関数のヘルプ一覧（リファレンス）: https://ichimomo.github.io/frasyr/reference/

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

- 将来予測関数の使い方：[wiki](https://github.com/ichimomo/frasyr/wiki/future_new)


<!--
- 新ルールのもとでの将来予測計算 https://ichimomo.github.io/frasyr/doc/future.html
- 管理基準値の計算 https://ichimomo.github.io/frasyr/doc/estMSY.html
-->



=======

frasyr_tool群の全体説明は[こちら](https://ichimomo.github.io/main/)
