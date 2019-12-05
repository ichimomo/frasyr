#### --- age-structured
tmb_future <- function(res_vpa,
                       SRmodel,
                       nsim = 1000, # number of simulation
                       nyear = 50, # number of future year
                       optim_method = "tmb", # or "R" or "both"
                       future_initial_year_name = 2017,
                       start_F_year_name = 2018,
                       start_random_rec_year_name = 2018,
                       x_init = 0,
                       x_lower = -3,
                       x_upper = 4,
                       compile=FALSE,
                       skip_setting=FALSE,
                       tmb_data=NULL,
                       seed_number=1
                       ) ## ここはvpa_nyear以下、任意
{

    # define age
    nage        <- nrow(res_vpa$naa)
    age_name    <- as.numeric(rownames(res_vpa$naa))
    recruit_age <- min(as.numeric(rownames(res_vpa$naa)))
    
   
    if(!isTRUE(skip_setting)){
        # define year
        vpa_nyear           <- ncol(res_vpa$naa)
        future_initial_year <- which(colnames(res_vpa$naa)==future_initial_year_name)
        total_nyear         <- future_initial_year + nyear
        allyear_name        <- min(as.numeric(colnames(res_vpa$naa)))+c(0:(total_nyear-1))
        allyear_label       <- c(rep("VPA",future_initial_year),rep("future",nyear))
        print(tibble(allyear_name, allyear_label))

        # define empty array
        waa_mat <- M_mat <- maa_mat <- naa_mat <- faa_mat <-
            array(0, dim=c(nage, total_nyear, nsim),
                  dimnames=list(age=age_name, year=allyear_name, nsim=1:nsim)) 
        resid_mat <- array(0, dim=c(total_nyear, nsim))

        # define biologial parameter
        waa_mat[,1:vpa_nyear,] <- as.matrix(res_vpa$input$dat$waa)
        maa_mat[,1:vpa_nyear,] <- as.matrix(res_vpa$input$dat$maa)
        M_mat  [,1:vpa_nyear,] <- as.matrix(res_vpa$input$dat$M)
        
        # assume future parameter（いろいろオプションがあるが、最終年のものを使う）
        waa_mat[,(vpa_nyear+1):total_nyear,] <- waa_mat[,vpa_nyear,]
        maa_mat[,(vpa_nyear+1):total_nyear,] <- maa_mat[,vpa_nyear,]
        M_mat  [,(vpa_nyear+1):total_nyear,] <- M_mat  [,vpa_nyear,]

        # define F scenario
        start_F_year <- which(allyear_name==start_F_year_name)    
        faa_mat[,1:vpa_nyear,] <- as.matrix(res_vpa$faa)    
        faa_mat[,start_F_year:total_nyear,] <- faa_mat[,vpa_nyear,]

        # define recruitment deviation
        set.seed(seed_number)
        start_random_rec_year  <- which(allyear_name==start_random_rec_year_name)
        random_rec_year_period <- (start_random_rec_year):total_nyear
        resid_mat[1:future_initial_year,] <- SRmodel$resid[1:future_initial_year]    
        resid_mat[random_rec_year_period,] <-
            rnorm(nsim*length(random_rec_year_period),
                  mean=-0.5 * (SRmodel$pars$sd)^2, sd=SRmodel$pars$sd)

        # define naa
        naa_mat[,1:vpa_nyear,] <- as.matrix(res_vpa$naa)

        # set data and parameter
        tmb_data <- list(naa_mat=naa_mat,
                         #                     deviance_init=0,
                         SR = 0,
                         rec_par_a = SRmodel$pars$a,
                         rec_par_b = SRmodel$pars$b,
                         rec_par_rho = 0,
                         bias_corrected_mean = -0.5 * (SRmodel$pars$sd)^2, # bias correction factorを入れる
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
                         objective_value = 12000
                         )
    }

    if(optim_method=="tmb" | optim_method=="both"){
        
        # comple & load cpp file
        use_rvpa_tmb(TmbFile = "est_MSY_tmb",
                     CppDir = system.file("executable",package="frasyr"),
                     RunDir = getwd(), overwrite=compile) 

        objAD <- TMB::MakeADFun(tmb_data, list(x=x_init), DLL="est_MSY_tmb")
        
        msy_optim <- nlminb(objAD$par, objAD$fn, gr=objAD$gr,
                            lower=list(x=x_lower), upper=list(x=x_upper))#,contol=nlminb_control)

        multi_msy <- as.numeric(exp(msy_optim$par))
        msy <- exp(-msy_optim$objective)
        ad_report <- objAD$report()

        ssb <- ad_report$spawner_mat
        caa <- ad_report$catch_mat
        naa <- ad_report$N_mat
        faa <- ad_report$F_mat

        dimnames(ssb) <- list(year=allyear_name, nsim=1:nsim)
        dimnames(naa) <- dimnames(faa) <- dimnames(caa) <- 
            list(age=age_name, year=allyear_name, nsim=1:nsim)
        
        res_future_tmb <- tibble::lst(multi_msy,msy,ssb,naa,faa,tmb_data)
    }

    if(optim_method=="R" | optim_method=="both"){
        R_obj_fun <- function(x, tmb_data, what_return="obj"){
            tmb_data$x <- x
            tmb_data$what_return <- what_return
            obj <- do.call(est_MSY_R, tmb_data)
            return(obj)
        }
        msy_optim <- nlminb(start=0, objective=R_obj_fun, tmb_data=tmb_data,
                            lower=list(x=x_lower), upper=list(x=x_upper))

#        msy_optim <- optimize(f=R_obj_fun, tmb_data=tmb_data,
        #                            lower=x_lower, upper=x_upper)
        multi_msy <- exp(msy_optim$par)
        res_future_R <- R_obj_fun(msy_optim$par, tmb_data=tmb_data,
                                  what_return="stat")
        res_future_R$multi_msy <- multi_msy
    }

    if(optim_method=="tmb") return(res_future_tmb)
    if(optim_method=="R")   return(res_future_R)
    if(optim_method=="both") return(list(tmb=res_future_tmb, R=res_future_R))
}

est_MSY_R <- function(naa_mat,
                      SR,
                      rec_par_a,
                      rec_par_b,
                      rec_par_rho,
                      bias_corrected_mean,
                      rec_resid_mat,
                      waa_mat,
                      M_mat,
                      maa_mat,                      
                      faa_mat,
                      Pope,
                      total_nyear,
                      future_initial_year,
                      start_F_year,
                      nsim,
                      nage,
                      recruit_age,
                      obj_catch,
                      objective,
                      objective_value,
                      x,
                      what_return="obj"
                      ){

    F_mat <- N_mat <- catch_mat <- naa_mat
    F_mat[] <- N_mat[] <- catch_mat <- 0
    rec_deviance_mat <- array(0, dim=c(total_nyear, nsim))
    
    N_mat[,1:future_initial_year,] <- naa_mat[,1:future_initial_year,]
    spawner_mat <- apply(N_mat * waa_mat * maa_mat, c(2,3) , sum)

    F_mat[,1:(start_F_year-1),] <- faa_mat[,1:(start_F_year-1),]
    F_mat[,start_F_year:total_nyear,] <- faa_mat[,start_F_year:total_nyear,] * exp(x)

    for(i in 1:nsim){
        for(t in future_initial_year:total_nyear){
            spawner_mat[t,i] <- sum(N_mat[,t,i] * waa_mat[,t,i] * maa_mat[,t,i])

            if(t>future_initial_year){
                if(SR == 0) { # Hockey-stick
                    rec_pred1 <- spawner_mat[t-recruit_age,i]*rec_par_a
                    rec_pred2 <- rec_par_b*rec_par_a
                    N_mat[1,t,i] <- min(rec_pred1, rec_pred2)
                }
                if(SR == 1) { # Beverton-Holt
                    N_mat[1,t,i] <- rec_par_a*spawner_mat[t-recruit_age,i]/
                        (1+rec_par_b*spawner_mat[t-recruit_age,i]);
                }
                if(SR == 2) { # Ricker
                    N_mat[1,t,i] <- rec_par_a*spawner_mat[t-recruit_age,i]*exp(-rec_par_b*spawner_mat[t-recruit_age,i]);
                }
                N_mat[1,t,i] <- N_mat[1,t,i]*exp(rec_par_rho*rec_deviance_mat[t,i])*
                    exp(rec_resid_mat[t,i]);
            }

            # forward calculation 
            if(t<total_nyear){
                for(iage in 1:(nage-1)) {
                    N_mat[iage+1,t+1,i] <- N_mat[iage,t,i]*exp(-M_mat[iage,t,i]-F_mat[iage,t,i])
                }
                N_mat[nage,t+1,i] <- N_mat[nage,t+1,i] + N_mat[nage,t,i]*exp(-M_mat[nage,t,i]-F_mat[nage,t,i])
            }
        }
    }

    if(Pope){
        catch_mat <- N_mat*(1-exp(-F_mat))*exp(-M_mat/2) * waa_mat
    }
    else{
        catch_mat <- naa*(1-exp(-faa-M))*faa/(faa+M) * waa_mat
    }

    if(objective<2){
        if(obj_catch==0) obj <- mean(catch_mat[,total_nyear,])
        if(obj_catch==1) obj <- geomean(catch_mat[,total_nyear,])
    }
    else{
        if(obj_catch==0) obj <- mean(spawner_mat[,total_nyear,])
        if(obj_catch==1) obj <- geomean(spawner_mat[,total_nyear,])
    }

    if(objective==0) obj <- -log(obj)
    else{
        obj <- (log(obj/objective_value))^2
    }

    if(what_return=="obj")  return(obj)
    if(what_return=="stat")  return(list(naa=N_mat, faa=F_mat,ssb=spawner_mat))    
}

