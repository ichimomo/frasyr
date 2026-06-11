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
    # 1-5が1つずつ並んでいるベクトルを30年分5年ブロックでbackward-resampingする
    # 最初の5年は5のみ、次の5年は4と5、次の5年は3-5...とリサンプリングの範囲が広くなるはず
    resid_test <- rep(1:5,each=5)
    res <- purrr::map_dfr(1:100, function(x)
        sample_backward(resid_test, 30, 5) %>%
        matrix(5,6) %>% as_tibble(.name_repair=~stringr::str_c("V",1:6)))

    # 長さがdurationの倍数でない場合
    resid_test2 <- c(1,1,1,rep(1:5,each=5))
    res2 <- purrr::map_dfr(1:10000, function(x)
        sample_backward(resid_test2, 30, 5) %>%
        matrix(5,6) %>% as_tibble(.name_repair=~stringr::str_c("V",1:6)))
    
    # Rのバージョンによってsampleの内部が変わっているので、seedを同じにしても、バージョンの違うR間で異なる結果が得られるらしい　https://community.rstudio.com/t/getting-different-results-with-set-seed/31624/4
    # そのため、100回繰り返して（乱数の影響を減らすため）各ブロックの最小値をテストする
    try(expect_equivalent(apply(res,2,min),c(5,4,3,2,1,1)))
    try(expect_equivalent(apply(res2,2,min),c(5,4,3,2,1,1)))

    # 長さが異なる場合でも最後のブロックの残差の平均値は、リサンプリングする残差の平均値にだいたい一致するはず
    try(expect_equivalent(round(mean(res2$V6),1),round(mean(resid_test2),1)))
    
})

test_that("make_future_data: special_setting で指定配列を上書きできる (level 1)", {
  data(res_vpa_org)
  data(res_sr_HSL2)

  base_args <- list(
    res_vpa = res_vpa_org, nsim = 5, nyear = 5,
    future_initial_year_name = 2017,
    start_F_year_name = 2018, start_biopar_year_name = 2018,
    start_random_rec_year_name = 2018,
    waa_year = 2015:2017, waa_catch_year = 2015:2017,
    maa_year = 2015:2017, M_year = 2015:2017, faa_year = 2015:2017,
    start_ABC_year_name = 2019, res_SR = res_sr_HSL2, silent = TRUE
  )

  # 上書きなしの基準データ
  d0 <- do.call(make_future_data, base_args)

  # M_mat と同じ形・全要素を sentinel 値で埋めた置換配列
  M_replace <- d0$data$M_mat
  M_replace[] <- 0.123

  d1 <- do.call(make_future_data,
                c(base_args, list(special_setting = list(M_mat = M_replace))))

  # special_setting で M_mat が上書きされている
  expect_equal(as.numeric(d1$data$M_mat), rep(0.123, length(d1$data$M_mat)))
  # 他の配列（naa_mat）は基準と一致＝副作用がない
  expect_equal(d1$data$naa_mat, d0$data$naa_mat)
})

test_that("make_future_data: special_setting に存在しない名前を渡すとエラー (level 1)", {
  data(res_vpa_org); data(res_sr_HSL2)
  expect_error(
    make_future_data(
      res_vpa = res_vpa_org, nsim = 2, nyear = 3,
      future_initial_year_name = 2017,
      start_F_year_name = 2018, start_biopar_year_name = 2018,
      start_random_rec_year_name = 2018,
      waa_year = 2015:2017, waa_catch_year = 2015:2017,
      maa_year = 2015:2017, M_year = 2015:2017, faa_year = 2015:2017,
      start_ABC_year_name = 2019, res_SR = res_sr_HSL2, silent = TRUE,
      special_setting = list(no_such_array = 1)
    )
  )
})
