
#' ある系群から別の系群への移入が想定される場合の再生産関係の推定（マアジ用）
#'
#' @param SRdata1 移出が生じている系群の\code{SRdata}のオブジェクト
#' @param SRdata2 移入が生じている系群の\code{SRdata}のオブジェクト
#' @param mig_type 0: 移出入なし, 1: 一定の割合で移入が生じる, 2: 平均的な割合からランダムに変動する, 3:移出する割合が自己回帰モデルで変化する, 4: 移出する割合がランダムウォーク
#' @param mig_year 移出入が生じている期間
#' @encoding UTF-8
#' @export
fit.SRmig = function(
  SRdata1,
  SRdata2,
  SR = c("HS","HS"),
  AR = c(FALSE,FALSE),
  mig_year = NULL,
  mig_type = 0,
  one_patch = FALSE,
  # shared_error = FALSE,
  sigma_omega_ratio = 0, # if 0, omega is freely estimated independently of sigma
  omega_fix = NA,
  phase = FALSE,
  p0_list = NULL
) {
  argname <- ls()
  arglist <- lapply(argname,function(xx) eval(parse(text=xx)))
  names(arglist) <- argname
  
  rec1 = SRdata1$R
  rec2 = SRdata2$R
  ssb1 = SRdata1$SSB
  ssb2 = SRdata2$SSB
  
  # if (shared_error) {
  #   if (one_patch) stop("shared_error is effective only when one_patch = FALSE")
  #   # if (mig_type > 1) stop("shared_error is effective only when mig_type = 1")
  # }
  
  if (one_patch) {
    if (mig_type==0) {
      stop("mig_type == 0 cannot be used when one_patch = TRUE")
    }
  }
  
  if (any(SRdata2$year %in% SRdata1$year) == FALSE) {
    stop("all years of SRdata2 must be included in SRdata1$year")
  }
  iy1 = SRdata1$year-min(SRdata1$year)
  iy = iy1[SRdata1$year %in% SRdata2$year]
  if (is.null(mig_year)) {
    mig_occurrence = rep(1,length(iy1))
  } else {
      mig_occurrence = as.numeric(SRdata1$year %in% mig_year)
    }
  if (one_patch == 0) {
    if (length(SR)==1) {
      message("The same stock-recruit function is applied to two patches")
      SR = rep(SR,2)
    }
    if (length(AR)==1) {
      message("The same assumption on AR is applied to two patches")
      AR = rep(AR,2)
    }
  } else {
    if (length(SR)!=1) {
      message("The first SR is used when one_patch = TRUE")
      SR = SR[1]
    }
    if (length(AR)!=1) {
      message("The first AR is used when one_patch = TRUE")
      AR = AR[1]
    }
  }
  SR_type = NULL
  imax = ifelse(one_patch,1,2)
  for (i in 1:imax) {
    if (SR[i]=="HS") SR_type = c(SR_type,0)
    if (SR[i]=="BH") SR_type = c(SR_type,1)
    if (SR[i]=="RI") SR_type = c(SR_type,2)
    if (SR[i]=="RPS") SR_type = c(SR_type,3)
  }
  if (sigma_omega_ratio>0 && !is.na(omega_fix)) {
    stop("sigma_omega_ratio and omega_fix cannot be set simultaneously")
  }
  data_list=list(rec1=rec1,ssb1=ssb1,rec2=rec2,ssb2=ssb2,iy=iy,
                 SR=SR_type,mig_type=mig_type,
                 one_patch = as.numeric(one_patch),sigma_omega_ratio = sigma_omega_ratio, mig_occurrence = mig_occurrence)
  if (is.null(p0_list)) {
    rec_loga = rec_logb = log_sigma = NULL
    for (i in 1:imax) {
      if (SR[i]!="RPS") {
        if (i==1) {
          res0 = fit.SR(SRdata1,SR=SR[i],AR=TRUE,method="L2",out.AR=FALSE)
        } else {
          res0 = fit.SR(SRdata2,SR=SR[i],AR=TRUE,method="L2",out.AR=FALSE)
        }
      } else {
        if (i==1) {
          res0 = fit.SR(SRdata1,SR="BH",AR=TRUE,method="L2",out.AR=FALSE)
        } else {
          res0 = fit.SR(SRdata2,SR="BH",AR=TRUE,method="L2",out.AR=FALSE)
        }
      }
      rec_loga = c(rec_loga,res0$opt$par[1])
      rec_logb = c(rec_logb,res0$opt$par[2])
      log_sigma = c(log_sigma,log(res0$pars$sd))
    }
    if (one_patch) {
      logit_d = rec2/(rec1+rec2)
      logit_d = log(logit_d)-log(1-logit_d)
    } else {
      logit_d = sapply(1:length(rec1), function(i) {
        if (iy1[i] %in% iy) {
          p_half = 0.5*rec2[which(iy==i-1)]/rec1[i]
        } else {
          p_half = median(0.5*rec2/rec1[SRdata1$year %in% SRdata2$year])
        }
        log(p_half/(1-p_half))
      })
    }
    if (!is.na(omega_fix)) {
      log_omega = log(omega_fix)
    } else {
      log_omega = log(sd(logit_d))
    }
    param_init=list(rec_loga=rec_loga,rec_logb=rec_logb,log_sigma=log_sigma,trans_rho=rep(0,imax),
                    logit_d=logit_d,intercept_d=mean(logit_d),trans_phi=0,log_omega=log_omega)
    # if (mig_type == 4) {
    #   param_init$log_omega = param_init$log_omega +2
    # }
  } else{
    param_init = p0_list
  }
  
  map = list()
  map$trans_rho = 0:(imax-1)
  map$trans_rho[AR==FALSE] <- NA
  map$trans_rho = factor(map$trans_rho)
  if (mig_type==0) { # no migration
    map$intercept_d = factor(NA)
    map$trans_phi = factor(NA)
    map$log_omega = factor(NA)
    map$logit_d = rep(factor(NA),length(param_init$logit_d))
  }
  if (mig_type==1) { # constant
    map$trans_phi = factor(NA)
    map$log_omega = factor(NA)
    map$logit_d = rep(factor(NA),length(param_init$logit_d))
  }
  if (mig_type==2) { # white noise
    map$trans_phi = factor(NA)
  }
  if (mig_type==4) { # random walk
    map$intercept_d = factor(NA)
    map$trans_phi = factor(NA)
  }
  # if (one_patch) {
  #   if (mig_type==5) { # estimated as fixed effect
  #     map$intercept_d = factor(NA)
  #     map$trans_phi = factor(NA)
  #     map$log_omega = factor(NA)
  #   } else {
  #     map$logit_d = rep(factor(NA),length(param_init$logit_d))
  #   }
  # }
  map$rec_logb = 0:(imax-1)
  map$rec_logb[SR == "RPS"] <- NA
  map$rec_logb = factor(map$rec_logb)
  
  if (sigma_omega_ratio>0) map$log_omega = factor(NA) 
  if (!is.na(omega_fix)) map$log_omega = factor(NA)
  
  if (mig_type>0 && phase){
    par_list = param_init
    # objective = 99999
    # for (i in 1:10) {
      # phase 1
      map1 = map
      map1$intercept_d <- factor(NA)
      map1$logit_d <- rep(factor(NA), length(param_init$logit_d))
      f = MakeADFun(data=data_list, parameters=par_list, random = c("logit_d"), 
                    map = map1, DLL="TwoPatchSR")
      fit = nlminb(f$par, f$fn, f$gr)
      par_list = f$env$parList(fit$par)
      
      # phase 2
      map2 = map
      map_name = unique(names(fit$par))
      for(i in 1:length(map_name)) {
        map2[[map_name[i]]] <- rep(factor(NA), length(param_init[[map_name[i]]]))
      }
      f = MakeADFun(data=data_list, parameters=par_list, random = c("logit_d"), 
                    map = map2, DLL="TwoPatchSR")
      fit = nlminb(f$par, f$fn, f$gr)
      par_list = f$env$parList(fit$par)
      # if (abs(fit0$objective-objective) < 0.001) {break} else {objective = fit$objective}
    # }
    # phase 3
    f = MakeADFun(data=data_list, parameters=par_list, random = c("logit_d"), 
                  map = map, DLL="TwoPatchSR")
    fit = nlminb(f$par, f$fn, f$gr)
  } else {
    f = MakeADFun(data=data_list, parameters=param_init, random = c("logit_d"), 
                  map = map, DLL="TwoPatchSR")
    fit = nlminb(f$par, f$fn, f$gr,control=list(eval.max=2000,iter.max=2000))
      }
  
  Res <- list()
  Res$input <- arglist
  Res$opt <- fit
  Res$MakeADFun <- f
  rep = sdreport(f)
  Res$rep = rep
  Res$par_list = f$env$parList(fit$par)
  Res$pars = cbind(rep$value[names(rep$value)=="rec_a"],rep$value[names(rep$value)=="rec_b"],rep$value[names(rep$value)=="sigma"],rep$value[names(rep$value)=="rho"])
  rownames(Res$pars) = NULL
  colnames(Res$pars) = c("a","b","sd","rho")
  Res$pars = data.frame(Res$pars)
  for(i in 1:imax) {
    if (SR[i]=="RPS") Res$pars$b[i] <- NA
  }
  
  if (mig_type<3) phi=0
  if (mig_type==3) phi=as.numeric(rep$value[names(rep$value)=="phi"])
  if (mig_type==4) phi=1
  
  if (mig_type<2) {
    omega=0
  } else {
    omega = as.numeric(rep$value[names(rep$value)=="omega"])
  }
  intercept_d = f$env$parList(fit$par)$intercept_d
  intercept_d = 1/(1+exp(-intercept_d))
  Res$mig_pars = c(intercept_d,phi,omega)
  names(Res$mig_pars) = c("intercept_d","phi","omega")
  Res$mig_pars = as.data.frame(t(Res$mig_pars))
  if (mig_type ==0) Res$mig_pars$intercept_d = 0
  if (mig_type ==4) Res$mig_pars$intercept_d = NA
  if (mig_type ==5) Res$mig_pars$intercept_d = NA
  
  d = rep$value[names(rep$value)=="d"]
  if (one_patch) {
    # d_obs=rep$value[names(rep$value)=="d_obs"]
    d_obs=rec2/(rec1+rec2)
    Res$migration = data.frame(year=SRdata1$year,dispersal_rate=d,dispersal_rate_obs = d_obs)
  } else {
    migrant = rep$value[names(rep$value)=="migrant"]
    Res$migration = data.frame(year=SRdata1$year,dispersal_rate=d,migrant=migrant)
  }
  
  repro1 = rep$value[names(rep$value)=="repro1"]
  pred_repro1 = rep$value[names(rep$value)=="pred_repro1"]
  deviance1 = rep$value[names(rep$value)=="delta1"]
  repro2 = rep$value[names(rep$value)=="repro2"]
  pred_repro2 = rep$value[names(rep$value)=="pred_repro2"]
  # if (shared_error) {
  #   repro2 = repro2-migrant
  #   pred_repro2 = pred_repro2-migrant
  # }
  deviance2 = rep$value[names(rep$value)=="delta2"]
  
  if (one_patch) {
    SRdata_both = data.frame(year = SRdata1$year,R1 = SRdata1$R, SSB1 = SRdata1$SSB,
                             R2 = SRdata2$R, SSB2 = SRdata2$SSB) %>%
      mutate(R_both = R1+R2, SSB_both = SSB1+SSB2,predR = pred_repro1,
             deviance = deviance1,dispersal_rate = d,dispersal_rate_obs = d_obs)
    
  } else {
    SRdata1_est = as.data.frame(SRdata1)
    SRdata1_est = SRdata1_est %>%
      mutate(dispersal_rate = d,migrant=migrant,
             reproduction=repro1,pred_reproduction=pred_repro1,
             deviance=deviance1)
    
    SRdata2_est = as.data.frame(SRdata2)
    SRdata2_est = SRdata2_est %>%
      mutate(reproduction=repro2,pred_reproduction=pred_repro2,
             deviance=deviance2) %>%
      mutate(migrant=R-repro2) %>%
      mutate(prop_migrant=migrant/R) %>%
      select(year,SSB,R,migrant,prop_migrant,everything())
    
    Res$SRdata1_est = SRdata1_est
    Res$SRdata2_est = SRdata2_est
  }
  
  NN = length(SRdata1$R)+length(SRdata2$R)
  Res$loglik <- loglik <- -fit$objective
  Res$k <- k <- length(fit$par)
  Res$AIC <- -2*loglik+2*k
  Res$AICc <- Res$AIC+2*k*(k+1)/(NN-k-1)
  Res$BIC <- -2*loglik+k*log(NN)
  
  return(Res)
}
