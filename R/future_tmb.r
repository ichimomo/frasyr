#### --- age-structured
tmb_future <- function(res_vpa,
                       nsim = 1000, # number of simulation
                       nyear = 50, # number of future year
                       future_initial_year_name = 2017,
                       start_F_year_name = 2018,
                       start_random_rec_year_name = 2018,
                       compile=FALSE) ## ここはvpa_nyear以下、任意
{

    # define year
    nage  <- nrow(res_vpa$naa)
    vpa_nyear <- ncol(res_vpa$naa)
    future_initial_year <- which(colnames(res_vpa$naa)==future_initial_year_name)
    total_nyear <- future_initial_year + nyear
    allyear_name <- min(as.numeric(colnames(res_vpa$naa)))+c(0:(total_nyear-1))
    allyear_label <- c(rep("VPA",future_initial_year),rep("future",nyear))
    start_F_year <- which(allyear_name==start_F_year_name)
    start_random_rec_year <- which(allyear_name==start_random_rec_year_name)    
    age_name <- as.numeric(rownames(res_vpa$naa))
    print(tibble(allyear_name, allyear_label))

    recruit_age <- min(as.numeric(rownames(res_vpa$naa)))

    waa_mat <- M_mat <- maa_mat <- naa_mat <- faa_mat <- array(0, dim=c(nage, total_nyear, nsim), dimnames=list(age=age_name, year=allyear_name, nsim=1:nsim)) 
    resid_mat <- array(0, dim=c(total_nyear, nsim))

    # input VPA part parameter (VPA設定を使うか将来予測設定を使うか、２つのオプションあり。ここではデータがある年のものはそれを使うやりかた)
    waa_mat[,1:vpa_nyear,] <- as.matrix(res_vpa$input$dat$waa)
    maa_mat[,1:vpa_nyear,] <- as.matrix(res_vpa$input$dat$maa)
    M_mat  [,1:vpa_nyear,] <- as.matrix(res_vpa$input$dat$M)
    naa_mat[,1:vpa_nyear,] <- as.matrix(res_vpa$naa)
    faa_mat[,1:vpa_nyear,] <- as.matrix(res_vpa$faa)
    resid_mat[1:future_initial_year,] <- SRmodel.base$resid[1:future_initial_year]

    # assume future parameter（いろいろオプションがあるが、最終年のものを使う）
    random_rec_year_period <- (start_random_rec_year):total_nyear
    waa_mat[,(vpa_nyear+1):total_nyear,] <- waa_mat[,vpa_nyear,]
    maa_mat[,(vpa_nyear+1):total_nyear,] <- maa_mat[,vpa_nyear,]
    M_mat  [,(vpa_nyear+1):total_nyear,] <- M_mat  [,vpa_nyear,]
    faa_mat[,start_F_year:total_nyear,] <- faa_mat[,vpa_nyear,]
    resid_mat[random_rec_year_period,] <-
        rnorm(nsim*length(random_rec_year_period),
              mean=-0.5 * (SRmodel.base$pars$sd)^2, sd=SRmodel.base$pars$sd)

    tmb_data_msy <- list(naa_mat=naa_mat,
                         #                     deviance_init=0,
                         SR = 0,
                         rec_par_a = SRmodel.base$pars$a,
                         rec_par_b = SRmodel.base$pars$b,
                         rec_par_rho = 0,
                         bias_corrected_mean = -0.5 * (SRmodel.base$pars$sd)^2, # bias correction factorを入れる
                         rec_resid_mat = resid_mat,
                         waa_mat = waa_mat,
                         maa_mat = maa_mat,                     
                         M_mat = M_mat,
                         faa_mat = faa_mat,                     
                         Pope = ifelse(isTRUE(res_vpa$input$Pope),1,0),
                         total_nyear = total_nyear,                     
                         future_initial_year = future_initial_year,
                         start_F_year=start_F_year,
                         nsim = nsim,
                         nage = nage,
                         recruit_age = recruit_age,                                          
                         obj_catch = 0, # 0: mean, 1:geomean
                         objective = 0, # 0: MSY, 1: PGY, 2: percentB0 or Bempirical
                         objective_value = 12000,
                         num_to_mass_scale = 1
                         )
    x_init <- 0; x_lower <- -3; x_upper <- 4

    # comple & load cpp file
    use_rvpa_tmb(TmbFile = "est_MSY_tmb",
                 CppDir = system.file("executable",package="frasyr"),
                 RunDir = getwd(), overwrite=compile) 

    objAD <- TMB::MakeADFun(tmb_data_msy, list(x=x_init), DLL="est_MSY_tmb")

    msy_optim <- nlminb(objAD$par, objAD$fn, gr=objAD$gr,
                        lower=list(x=x_lower), upper=list(x=x_upper))#,contol=nlminb_control)

    multi_msy <- as.numeric(exp(msy_optim$par))
    ssb_msy <- 
    msy <- exp(-msy_optim$objective)

    #    dimnames(N_mat) <- dimnames(F_mat) <- list()

    ssb <- objAD$report()$spawner_mat
    naa <- objAD$report()$N_mat
    faa <- objAD$report()$F_mat

    dimnames(ssb) <- list(year=allyear_name, nsim=1:nsim)
    dimnames(naa) <- dimnames(faa) <-
        list(age=age_name, year=allyear_name, nsim=1:nsim)
    
    # F=0, sd=0の場合の平衡状態のnumbers at age= 1486.19  901.42  546.74  842.79 一致
    # F=0, sd=0の場合の平衡状態のssb= 491083.1 一致
    #objAD$report()$spawner_mat[100,1] 
    #objAD$report()$spawner_mat[100,1]
    # sd=0のときのFmultiplier　frasyr; F multiplier= 0.5268391 , tmb= 0.5268342 (6桁目でずれるがまあまあOK)
    list(multi_msy=multi_msy,
         msy=msy,         
         ssb=ssb,
         naa=naa,
         faa=faa)
}


