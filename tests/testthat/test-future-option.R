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



