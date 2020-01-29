library(frasyr)

context("new future_vpa option")

test_that("utility function check",{

    # check apply_year_colum function
    tmpres <- as.numeric(apply_year_colum(matrix(1:20,4,5,dimnames=list(1:4,1:5)),target_year=-1:-2))
    for(i in 1:4) expect_equal(tmpres[i],i+14)
    tmpres <- as.numeric(apply_year_colum(matrix(1:20,4,5,dimnames=list(1:4,1:5)),target_year=4:5))
    for(i in 1:4) expect_equal(tmpres[i],i+14)

    # check sample_backward function
    set.seed(1)
    res <- sample_backward(rep(1:5,each=5), 30, 5)
    try(expect_equal(apply(matrix(res,5,6),2,min),c(5,4,3,2,2,1)))
    
})

test_that("future_vpa function",{

data(res_vpa)
data(res_sr_HSL2)

data_future_test <- make_future_data(res_vpa, # VPAの結果
	         nsim = 1000, # シミュレーション回数
                 nyear = 50, # 将来予測の年数
                 future_initial_year_name = 2017, # 年齢別資源尾数を参照して将来予測をスタートする年
                 start_F_year_name = 2018, # この関数で指定したFに置き換える最初の年
                 start_biopar_year_name=2018, # この関数で指定した生物パラメータに置き換える最初の年
                 start_random_rec_year_name = 2018, # この関数で指定した再生産関係からの加入の予測値に置き換える最初の年
                 # biopar setting
                 waa_year=2015:2017, waa=NULL, # 将来の年齢別体重の設定。過去の年を指定し、その平均値を使うか、直接ベクトルで指定するか。以下も同じ。
                 waa_catch_year=2015:2017, waa_catch=NULL,
                 maa_year=2015:2017, maa=NULL,
                 M_year=2015:2017, M=NULL,
                 # faa setting
                 faa_year=2015:2017, # currentF, futureFが指定されない場合だけ有効になる。将来のFを指定の年の平均値とする
                 currentF=NULL,futureF=NULL, # 将来のABC.year以前のFとABC.year以降のFのベクトル 
                 # HCR setting (not work when using TMB)
                 start_ABC_year_name=2019, # HCRを適用する最初の年
                 HCR_beta=1, # HCRのbeta
                 HCR_Blimit=-1, # HCRのBlimit
                 HCR_Bban=-1, # HCRのBban
                 HCR_year_lag=0, # HCRで何年遅れにするか
                 # SR setting
                 res_SR=res_sr_HSL2, # 将来予測に使いたい再生産関係の推定結果が入っているfit.SRの返り値
                 seed_number=1, # シード番号
                 resid_type="lognormal", # 加入の誤差分布（"lognormal": 対数正規分布、"resample": 残差リサンプリング）
                 resample_year_range=0, # リサンプリングの場合、残差をリサンプリングする年の範囲
                 bias_correction=TRUE, # バイアス補正をするかどうか
                 recruit_intercept=0, # 移入や放流などで一定の加入がある場合に足す加入尾数
                 # Other
                 Pope=res_vpa$input$Pope,
		 fix_recruit=list(year=c(2020,2021),rec=c(1000,2000)),
		 fix_wcatch=list(year=c(2020,2021),wcatch=c(1000,2000))		 
                 ) 
# data_future_testには、将来予測に使うデータ(data)とdata_futureを作るときに使った引数一覧(input)が入っている
names(data_future_test)
# data_future_test$dataには年齢別資源尾数naa_matなど。naa_matの将来予測部分にはまだNAが入っており、次のfuture_vpa関数でこのNAを埋める
names(data_future_test$data)

# backward-resamplingの場合
data_future_backward <- make_future_data(res_vpa, # VPAの結果
	         nsim = 1000, # シミュレーション回数
                 nyear = 50, # 将来予測の年数
                 future_initial_year_name = 2017, # 年齢別資源尾数を参照して将来予測をスタートする年
                 start_F_year_name = 2018, # この関数で指定したFに置き換える最初の年
                 start_biopar_year_name=2018, # この関数で指定した生物パラメータに置き換える最初の年
                 start_random_rec_year_name = 2018, # この関数で指定した再生産関係からの加入の予測値に置き換える最初の年
                 # biopar setting
                 waa_year=2015:2017, waa=NULL, # 将来の年齢別体重の設定。過去の年を指定し、その平均値を使うか、直接ベクトルで指定するか。以下も同じ。
                 waa_catch_year=2015:2017, waa_catch=NULL,
                 maa_year=2015:2017, maa=NULL,
                 M_year=2015:2017, M=NULL,
                 # faa setting
                 faa_year=2015:2017, # currentF, futureFが指定されない場合だけ有効になる。将来のFを指定の年の平均値とする
                 currentF=NULL,futureF=NULL, # 将来のABC.year以前のFとABC.year以降のFのベクトル 
                 # HCR setting (not work when using TMB)
                 start_ABC_year_name=2019, # HCRを適用する最初の年
                 HCR_beta=1, # HCRのbeta
                 HCR_Blimit=-1, # HCRのBlimit
                 HCR_Bban=-1, # HCRのBban
                 HCR_year_lag=0, # HCRで何年遅れにするか
                 # SR setting
                 res_SR=res_sr_HSL2, # 将来予測に使いたい再生産関係の推定結果が入っているfit.SRの返り値
                 seed_number=1, # シード番号
                 resid_type="backward", # 加入の誤差分布（"lognormal": 対数正規分布、"resample": 残差リサンプリング）
                 resample_year_range=0, # リサンプリングの場合、残差をリサンプリングする年の範囲
                 backward_duration=5,
                 bias_correction=TRUE, # バイアス補正をするかどうか
                 recruit_intercept=0, # 移入や放流などで一定の加入がある場合に足す加入尾数
                 # Other
                 Pope=res_vpa$input$Pope,
		 fix_recruit=list(year=c(2020,2021),rec=c(1000,2000)),
		 fix_wcatch=list(year=c(2020,2021),wcatch=c(1000,2000))		 
                 )     

# 単なる将来予測の場合
res_future_test <- future_vpa(tmb_data=data_future_test$data,
		              optim_method="none", 
                    	      multi_init = 1) 

# 単なる将来予測の場合
res_future_backward <- future_vpa(tmb_data=data_future_backward$data, 
		              optim_method="none", 
                    	      multi_init = 1) 

# MSY計算の場合
res_future_test_R <- future_vpa(tmb_data=data_future_test$data, 
		              optim_method="R", 
                    	      multi_init  = 1,
			      multi_lower = 0.001, multi_upper = 5,
			      objective="MSY")
# [1] 0.5269326
expect_equal(round(res_future_test_R$multi,3),0.527)

if(sum(installed.packages()[,1]=="TMB")){
    res_future_test_tmb <- future_vpa(tmb_data=data_future_test$data,
                                      optim_method="tmb", 
                                      multi_init  = 1,
                                      multi_lower = 0.001, multi_upper = 5,
                                      objective="MSY")
    expect_equal(round(res_future_test_tmb$multi,3),0.527)
}


})

