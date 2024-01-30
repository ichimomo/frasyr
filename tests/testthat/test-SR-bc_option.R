load(system.file("extdata","res_vpa_pma.rda",package = "frasyr"))
SRdata_pma_check_weight <- get.SRdata(res_vpa_pma, weight.year=1990:2000)

SRdata_pma_check_weight$R[] <- rnorm(30,log(mean(SRdata_pma_check_weight$R)),0.2) %>% exp()
    
res_SR1 <- fit.SR(SRdata_pma_check_weight, SR="HS", AR=0, method="L2", bias_correct=FALSE)
res_SR2 <- fit.SR(SRdata_pma_check_weight, SR="HS", AR=0, method="L2", bias_correct=TRUE)
expect_equal(res_SR2$pars$a*res_SR2$pars$b * exp(-0.5*res_SR2$pars$sd^2),
             res_SR1$pars$a*res_SR1$pars$b, tol=0.001)

res_SR1 <- fit.SR(SRdata_pma_check_weight, SR="BH", AR=0, method="L2", bias_correct=FALSE)
res_SR2 <- fit.SR(SRdata_pma_check_weight, SR="BH", AR=0, method="L2", bias_correct=TRUE, p0=res_SR1$opt$par)

#plot(res_SR1$pred$SSB,xx1 <- res_SR1$pred$R)
#points(res_SR2$pred$SSB,xx2 <- res_SR2$pred$R * exp(-0.5*res_SR2$pars$sd^2),col=2)
xx1 <- res_SR1$pred$R
xx2 <- res_SR2$pred$R * exp(-0.5*res_SR2$pars$sd^2)
expect_equal(xx1,xx2,tol=0.1)

# test for BHS (BHSでうまく推定できない)
SRdata_pma_check_weight <- get.SRdata(res_vpa_pma, years=1990:2000)
res_SR1 <- fit.SR(SRdata_pma_check_weight, SR="HS",  AR=0, method="L2", bias_correct=FALSE)
res_SR2 <- fit.SR(SRdata_pma_check_weight, SR="BHS", AR=0, method="L2", bias_correct=FALSE, gamma=1, p0=res_SR1$opt$par)

plot_SR(res_SR1)

# 正しいパラメータを入れるとちゃんと描画されるので、関数が間違っていることはない
res_SR2$pars[2] <- res_SR1$pars[2]
plot_SR(res_SR2)

