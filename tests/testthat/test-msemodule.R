test_that("VPA bootstrap procedure", {

  data("res_vpa_example2")
  res.boot <- boo.vpa(res_vpa_example2, B_ite=5, type="index", method="p")
  SR.boot1  <- do_bootSR(res.boot, res_vpa_example2, type="fix", SR="BH")
  SR.boot2  <- do_bootSR(res.boot, res_vpa_example2, type="auto")
  SR.boot3  <- do_bootSR(res.boot, res_vpa_example2, type="regime", regime.year=2000)    
  
  
  SR.boot1 %>% map_dfr(function(x) x$pars)
  SR.boot2 %>% map_dfr(function(x) x$pars)
  SR.boot3 %>% map_dfr(function(x) x$regime_pars, .id="id")

 
})
