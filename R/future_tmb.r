make_array <- function(d3_mat, pars, pars.year, year_replace_future){
    years <- dimnames(d3_mat)[[2]]
    if(is.null(pars)){
        pars.future <- rowMeans(d3_mat[,years%in%pars.year,1])
    }
    else{
        if(length(pars)==dim(d3_mat)[[1]]) pars.future <- pars
        else stop("length of parameter is bad..")
    }
    d3_mat[,which(year_replace_future==years):length(years),] <- pars.future
    d3_mat
}

set_SR_mat_lognormal <- function(res_vpa, res_SR, SR_mat, seed_number, start_random_rec_year_name, future_initial_year){
    
    allyear_name <- dimnames(SR_mat)[[1]]
    start_random_rec_year  <- which(allyear_name==start_random_rec_year_name)
    random_rec_year_period <- (start_random_rec_year):length(allyear_name)
    
  
    # define SR function
    if(res_SR$input$SR=="HS") SR_mat[random_rec_year_period,,"SR_type"] <- 1
    if(res_SR$input$SR=="BH") SR_mat[random_rec_year_period,,"SR_type"] <- 2
    if(res_SR$input$SR=="RI") SR_mat[random_rec_year_period,,"SR_type"] <- 3

    # define SR parameter
    SR_mat[random_rec_year_period,,"a"] <- res_SR$pars$a
    SR_mat[random_rec_year_period,,"b"] <- res_SR$pars$b
    SR_mat[random_rec_year_period,,"rho"] <- res_SR$pars$rho

    # define recruitment deviation
    # あとはバイアス補正をすること
    set.seed(seed_number)
    sd_with_AR <- res_SR$pars$sd/(1-res_SR$pars$rho^2)
    SR_mat[1:(start_random_rec_year-1),,"deviance"] <- res_SR$resid[1:(start_random_rec_year-1)]    
    SR_mat[random_rec_year_period,,"rand_resid"] <- rnorm(dim(SR_mat)[[2]]*length(random_rec_year_period), mean=0, sd=res_SR$pars$sd)
    for(t in random_rec_year_period){
        SR_mat[t, ,"deviance"] <- SR_mat[t-1, ,"deviance"]*SR_mat[t,,"rho"] + SR_mat[t, ,"rand_resid"] - - 0.5* sd_with_AR^2
    }
    SR_mat[random_rec_year_period,,"deviance"]   <- SR_mat[random_rec_year_period,,"deviance"] - 0.5* sd_with_AR^2
    ### バイアス補正のやりかたがよくわからない。。。じっくり考えること
    return(SR_mat)
}

#### make data for future projection
make_tmb_data <- function(res_vpa,
                          nsim = 1000, # number of simulation
                          nyear = 50, # number of future year
                          future_initial_year_name = 2017,                          
                          # biopar setting
                          waa_year, waa=NULL,
                          waa.catch_year, waa.catch=NULL,
                          maa_year, maa=NULL,
                          M_year, M=NULL,
                          start_biopar_year_name=2018,
                          # faa setting
                          faa_year, faa=NULL,
                          start_F_year_name = 2018,                          
                          
                          # SR setting
                          res_SR=NULL,                       
                          start_random_rec_year_name = 2018,
                          seed_number=1
                          ) 
{

    argname <- ls()
    input <- lapply(argname,function(x) eval(parse(text=x)))
    names(input) <- argname
    
    # define age and year
    nage        <- nrow(res_vpa$naa)
    age_name    <- as.numeric(rownames(res_vpa$naa))
    recruit_age <- min(as.numeric(rownames(res_vpa$naa)))
   
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
    SR_mat <- array(0, dim=c(total_nyear, nsim, 6),
                    dimnames=list(year=allyear_name, nsim=1:nsim, par=c("rand_resid","deviance","SR_type","a","b","rho")))

    # fill vpa data 
    waa_mat[,1:vpa_nyear,] <- as.matrix(res_vpa$input$dat$waa)
    maa_mat[,1:vpa_nyear,] <- as.matrix(res_vpa$input$dat$maa)
    M_mat  [,1:vpa_nyear,] <- as.matrix(res_vpa$input$dat$M)
    faa_mat[,1:vpa_nyear,] <- as.matrix(res_vpa$faa)
    naa_mat[,1:vpa_nyear,] <- as.matrix(res_vpa$naa)

    waa_mat <- make_array(waa_mat, waa, waa_year, start_biopar_year_name)
    maa_mat <- make_array(maa_mat, maa, maa_year, start_biopar_year_name)
    M_mat   <- make_array(M_mat  , M  , M_year  , start_biopar_year_name)
    faa_mat <- make_array(faa_mat, faa, faa_year, start_F_year_name)
    start_F_year <- which(allyear_name==start_F_year_name)
    #    waa.catch_mat <- make_array(waa.catch_mat, waa.catch, waa.catch_year)

    SR_mat <- set_SR_mat_lognormal(res_vpa, res_SR, SR_mat, seed_number, start_random_rec_year_name, future_initial_year)    

    sd_with_AR <- res_SR$pars$sd/(1-res_SR$pars$rho^2)
    
    # set data and parameter
    tmb_data <- list(naa_mat=naa_mat,
                     SR_mat = SR_mat[,,"SR_type"],
                     rec_par_a_mat =  SR_mat[,,"a"],
                     rec_par_b_mat =  SR_mat[,,"b"],
                     rec_par_rho_mat =  SR_mat[,,"rho"],
                     bias_corrected_mean = -0.5 * (sd_with_AR)^2, # bias correction factorを入れる
                     rec_resid_mat = SR_mat[,,"deviance"],
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

    argname <- ls()
    arglist <- lapply(argname,function(x) eval(parse(text=x)))
    names(arglist) <- argname    
    
    return(tibble::lst(data=tmb_data,input=input))
}

future_tmb <- function(tmb_data,optim_method="tmb",
                       x_init = 0,
                       x_lower = -3,
                       x_upper = 4,
                       compile=FALSE){
                           
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

        dimnames(ssb) <- dimnames(tmb_data$naa_mat)[2:3]
        dimnames(naa) <- dimnames(faa) <- dimnames(caa) <- 
            dimnames(tmb_data$naa_mat)
        
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

