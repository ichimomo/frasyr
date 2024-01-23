
res_SR1 <- fit.SR(SRdata_pma_check_weight, SR="HS", AR=0, method="L2", bias_correct=FALSE)
res_SR2 <- fit.SR(SRdata_pma_check_weight, SR="HS", AR=0, method="L2", bias_correct=TRUE)
expect_equal(res_SR2$pars$a*res_SR2$pars$b * exp(-0.5*res_SR2$pars$sd^2),
             res_SR1$pars$a*res_SR1$pars$b, tol=0.001)

res_SR1 <- fit.SR(SRdata_pma_check_weight, SR="BH", AR=0, method="L2", bias_correct=FALSE)
res_SR2 <- fit.SR(SRdata_pma_check_weight, SR="BH", AR=0, method="L2", bias_correct=TRUE)
bind_rows(res_SR1$par, res_SR2$par)

res_SR1 <- fit.SR(SRdata_pma_check_weight, SR="RI", AR=0, method="L2", bias_correct=FALSE)
res_SR2 <- fit.SR(SRdata_pma_check_weight, SR="RI", AR=0, method="L2", bias_correct=TRUE)
bind_rows(res_SR1$par, res_SR2$par)
