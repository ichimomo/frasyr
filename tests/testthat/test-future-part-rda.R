library(frasyr)

context("future ref.F")

test_that("oututput value check",{
  load(system.file("extdata","res_vpa_pma.rda",package = "frasyr"))

  res_ref_f_pma_check <- ref.F(res_vpa_pma,Fcurrent=NULL,waa=NULL,maa=NULL,M=NULL,waa.catch=NULL,M.year=NULL,
                               waa.year=NULL,maa.year=NULL,rps.year = NULL,max.age = Inf,min.age = 0,
                               d = 0.001,Fem.init = 0.5,Fmax.init = 1.5,F0.1.init = 0.7,pSPR = seq(10,90,by=10),
                               iterlim=1000,plot=TRUE,Pope=FALSE,F.range = seq(from=0,to=2,length=101) )

  res.ref.f_2times <- ref.F(res_vpa_pma,Fcurrent=res_vpa_pma$Fc.at.age*2,
                            waa=NULL,maa=NULL,M=NULL,waa.catch=NULL,M.year=NULL,
                            waa.year=NULL,maa.year=NULL,rps.year = NULL,max.age = Inf,min.age = 0,
                            d = 0.001,Fem.init = 0.5,Fmax.init = 1.5,F0.1.init = 0.7,
                            pSPR = c(seq(10,90,by=10),res_ref_f_pma_check$currentSPR$perSPR*100),
                            iterlim=1000,plot=TRUE,Pope=FALSE,F.range = seq(from=0,to=2,length=101) )


  times2_check <- as.numeric(table(unlist(res_ref_f_pma_check$summary/res.ref.f_2times$summary[,1:16])))

  expect_equal(times2_check[1],2)
  expect_equal(times2_check[2],31)
  expect_equal(times2_check[3],15)
  expect_equal(round(res.ref.f_2times$summary[3,17],2),0.5)


  #上記引数での計算結果を読み込み
  load(system.file("extdata","res_ref_f_pma.rda",package = "frasyr"))

  #結果の数値を照合
  for(i in 1:length(res_ref_f_pma_check$sel)){
    expect_equal(res_ref_f_pma$sel[i], res_ref_f_pma_check$sel[i])
  }
  expect_equal(res_ref_f_pma$max.age,res_ref_f_pma_check$max.age)
  expect_equal(res_ref_f_pma$min.age,res_ref_f_pma_check$min.age)
  for(i in 1:length(res_ref_f_pma_check$rps.q)){
    expect_equal(res_ref_f_pma$rps.q[i], res_ref_f_pma_check$rps.q[i])
  }
  for(i in 1:length(res_ref_f_pma_check$spr.q)){
    expect_equal(res_ref_f_pma$spr.q[i], res_ref_f_pma_check$spr.q[i])
  }
  for(i in 1:length(res_ref_f_pma_check$Fcurrent)){
    expect_equal(res_ref_f_pma$Fcurrent[i], res_ref_f_pma_check$Fcurrent[i])
  }
  for(i in 1:length(res_ref_f_pma_check$Fmed)){
    expect_equal(res_ref_f_pma$Fmed[i], res_ref_f_pma_check$Fmed[i])
  }
  for(i in 1:length(res_ref_f_pma_check$Flow)){
    expect_equal(res_ref_f_pma$Flow[i], res_ref_f_pma_check$Flow[i])
  }
  for(i in 1:length(res_ref_f_pma_check$Fhigh)){
    expect_equal(res_ref_f_pma$Fhigh[i], res_ref_f_pma_check$Fhigh[i])
  }
  for(i in 1:length(res_ref_f_pma_check$Fmax)){
    expect_equal(res_ref_f_pma$Fmax[i], res_ref_f_pma_check$Fmax[i])
  }
  for(i in 1:length(res_ref_f_pma_check$F0.1)){
    expect_equal(res_ref_f_pma$F0.1[i], res_ref_f_pma_check$F0.1[i])
  }
  for(i in 1:length(res_ref_f_pma_check$Fmean)){
    expect_equal(res_ref_f_pma$Fmean[i], res_ref_f_pma_check$Fmean[i])
  }
  expect_equal(res_ref_f_pma$rps.data,res_ref_f_pma_check$rps.data)
  for(i in 1:nrow(res_ref_f_pma_check$FpSPR)){
    for(j in 1:ncol(res_ref_f_pma_check$FpSPR)){
      expect_equal(res_ref_f_pma$FpSPR[i,j],res_ref_f_pma_check$FpSPR[i,j])
    }
  }
  expect_equal(res_ref_f_pma$summary,res_ref_f_pma_check$summary)

  expect_equal(res_ref_f_pma$ypr.spr,res_ref_f_pma_check$ypr.spr)
  for(i in 1:length(res_ref_f_pma_check$waa)){
    expect_equal(res_ref_f_pma$waa[i], res_ref_f_pma_check$waa[i])
  }
  for(i in 1:length(res_ref_f_pma_check$waa.catch)){
    expect_equal(res_ref_f_pma$waa.catch[i], res_ref_f_pma_check$waa.catch[i])
  }
  for(i in 1:length(res_ref_f_pma_check$maa)){
    expect_equal(res_ref_f_pma$maa[i], res_ref_f_pma_check$maa[i])
  }
  expect_equal(res_ref_f_pma$spr0,res_ref_f_pma_check$spr0)
})

context("future SRdata")

test_that("oututput value check",{
  load(system.file("extdata","res_vpa_pma.rda",package = "frasyr"))

  SRdata0 <- get.SRdata(R.dat=exp(rnorm(10)),SSB.dat=exp(rnorm(10)))
  SRdata0usingPeriodFrom1990To2000 <- get.SRdata(res_vpa_pma,years=1990:2000)

  #上記引数での計算結果を読み込み
  SRdata0_pma_check <- read.csv(system.file("extdata","future_SRdata0_pma_check.csv",package="frasyr"),row.names=1)
  SRdata0usingPeriodFrom1990To2000_pma_check <- read.csv(system.file("extdata","future_SRdata0usingPeriodFrom1990To2000_pma_check.csv",package="frasyr"),row.names=1)

  #結果の数値を照合

  expect_equal(SRdata0$year, SRdata0_pma_check$year)
  #expect_equal(SRdata0$SSB, SRdata0_pma_check$SSB)
  #expect_equal(SRdata0$R, SRdata0_pma_check$R)

  expect_equal(SRdata0usingPeriodFrom1990To2000$year, SRdata0usingPeriodFrom1990To2000_pma_check$year)
  expect_equal(SRdata0usingPeriodFrom1990To2000$SSB, SRdata0usingPeriodFrom1990To2000_pma_check$SSB)
  expect_equal(SRdata0usingPeriodFrom1990To2000$R, SRdata0usingPeriodFrom1990To2000_pma_check$R)

})

context("future future.vpa")

test_that("oututput value check (iteration for future sim is fixed as 2) ",{

  load(system.file("extdata","res_vpa_pma.rda",package = "frasyr"))
  load(system.file("extdata","SRdata_pma.rda",package = "frasyr"))

  SRmodel.list <- expand.grid(SR.rel = c("HS","BH","RI"), AR.type = c(0, 1), out.AR=c(TRUE,FALSE), L.type = c("L1", "L2"))
  SR.list <- list()

  for (i in 1:nrow(SRmodel.list)) {
    SR.list[[i]] <- fit.SR(SRdata_pma, SR = SRmodel.list$SR.rel[i], method = SRmodel.list$L.type[i],
                           AR = SRmodel.list$AR.type[i], out.AR=SRmodel.list$out.AR[i], hessian = FALSE)
  }

  for (i in 1:nrow(SRmodel.list)) {
    infile.name <- sprintf("SRpma_%s_%s_AR%d_outAR%d.rda",SRmodel.list$SR.rel[i],SRmodel.list$L.type[i], SRmodel.list$AR.type[i],SRmodel.list$out.AR[i])
    resfitSR <- load(system.file("extdata",infile.name,package = "frasyr"))

    fittedSR <- eval(parse(text=resfitSR))
    fres_pma_recarg_list <-list(a=fittedSR$pars$a,b=fittedSR$pars$b,
                                rho=fittedSR$pars$rho, # ここではrho=0なので指定しなくてもOK
                                sd=fittedSR$pars$sd,resid=fittedSR$resid)

    selectedrecSR <- sprintf("%s.recAR",fittedSR$input$SR[1])

    fres_pma_check <-future.vpa(res0=res_vpa_pma,
                                multi=1, # res.pma$Fc.at.ageに掛ける乗数
                                nyear=50, # 将来予測の年数
                                start.year=2012, # 将来予測の開始年
                                N=2, # 確率的計算の繰り返し回数
                                ABC.year=2013, # ABCを計算する年
                                waa.year=2009:2011, # 生物パラメータの参照年
                                maa.year=2009:2011,
                                M.year=2009:2011,
                                is.plot=TRUE, # 結果をプロットするかどうか
                                seed=1,
                                silent=TRUE,
                                recfunc=eval(parse(text=selectedrecSR)), # 再生産関係の関数
                                # recfuncに対する引数
                                rec.arg=fres_pma_recarg_list)

    assign(sprintf("fres_pma_%s_%s_AR%d_outAR%d_check",SRmodel.list$SR.rel[i],SRmodel.list$L.type[i], SRmodel.list$AR.type[i],SRmodel.list$out.AR[i]),fres_pma_check)

  }

  # HS L1 AR0 ----
  load(system.file("extdata","fres_pma_HS_L1_AR0_outAR0.rda",package = "frasyr"))
  for( k in 1:length(fres_pma_HS_L1_AR0_outAR0_check$faa[1,1,])){
    for( i in 1:nrow(fres_pma_HS_L1_AR0_outAR0_check$faa[,,k])){
      for( j in 1:ncol(fres_pma_HS_L1_AR0_outAR0_check$faa[,,k])){
        expect_equal(fres_pma_HS_L1_AR0_outAR0$faa[i,j,k], fres_pma_HS_L1_AR0_outAR0_check$faa[i,j,k])

      }
    }
  }

  for( k in 1:length(fres_pma_HS_L1_AR0_outAR0_check$naa[1,1,])){
    for( i in 1:nrow(fres_pma_HS_L1_AR0_outAR0_check$naa[,,k])){
      for( j in 1:ncol(fres_pma_HS_L1_AR0_outAR0_check$naa[,,k])){
        expect_equal(fres_pma_HS_L1_AR0_outAR0$naa[i,j,k], fres_pma_HS_L1_AR0_outAR0_check$naa[i,j,k])

      }
    }
  }

  for( k in 1:length(fres_pma_HS_L1_AR0_outAR0_check$biom[1,1,])){
    for( i in 1:nrow(fres_pma_HS_L1_AR0_outAR0_check$biom[,,k])){
      for( j in 1:ncol(fres_pma_HS_L1_AR0_outAR0_check$biom[,,k])){
        expect_equal(fres_pma_HS_L1_AR0_outAR0$biom[i,j,k], fres_pma_HS_L1_AR0_outAR0_check$biom[i,j,k])

      }
    }
  }

  for( k in 1:length(fres_pma_HS_L1_AR0_outAR0_check$naa[1,1,])){
    for( i in 1:nrow(fres_pma_HS_L1_AR0_outAR0_check$naa[,,k])){
      for( j in 1:ncol(fres_pma_HS_L1_AR0_outAR0_check$baa[,,k])){
        expect_equal(fres_pma_HS_L1_AR0_outAR0$baa[i,j,k], fres_pma_HS_L1_AR0_outAR0_check$baa[i,j,k])

      }
    }
  }

  for( k in 1:length(fres_pma_HS_L1_AR0_outAR0_check$ssb[1,1,])){
    for( i in 1:nrow(fres_pma_HS_L1_AR0_outAR0_check$ssb[,,k])){
      for( j in 1:ncol(fres_pma_HS_L1_AR0_outAR0_check$ssb[,,k])){
        expect_equal(fres_pma_HS_L1_AR0_outAR0$ssb[i,j,k], fres_pma_HS_L1_AR0_outAR0_check$ssb[i,j,k])

      }
    }
  }

  for( k in 1:length(fres_pma_HS_L1_AR0_outAR0_check$wcaa[1,1,])){
    for( i in 1:nrow(fres_pma_HS_L1_AR0_outAR0_check$wcaa[,,k])){
      for( j in 1:ncol(fres_pma_HS_L1_AR0_outAR0_check$wcaa[,,k])){
        expect_equal(fres_pma_HS_L1_AR0_outAR0$wcaa[i,j,k], fres_pma_HS_L1_AR0_outAR0_check$wcaa[i,j,k])

      }
    }
  }

  for( k in 1:length(fres_pma_HS_L1_AR0_outAR0_check$M[1,1,])){
    for( i in 1:nrow(fres_pma_HS_L1_AR0_outAR0_check$M[,,k])){
      for( j in 1:ncol(fres_pma_HS_L1_AR0_outAR0_check$M[,,k])){
        expect_equal(fres_pma_HS_L1_AR0_outAR0$M[i,j,k], fres_pma_HS_L1_AR0_outAR0_check$M[i,j,k])

      }
    }
  }

  for( i in 1:nrow(fres_pma_HS_L1_AR0_outAR0_check$rps)){
    for( j in 1:ncol(fres_pma_HS_L1_AR0_outAR0_check$rps)){
      expect_equal(fres_pma_HS_L1_AR0_outAR0$rps[i,j], fres_pma_HS_L1_AR0_outAR0_check$rps[i,j])

    }
  }

  for( k in 1:length(fres_pma_HS_L1_AR0_outAR0_check$maa[1,1,])){
    for( i in 1:nrow(fres_pma_HS_L1_AR0_outAR0_check$maa[,,k])){
      for( j in 1:ncol(fres_pma_HS_L1_AR0_outAR0_check$maa[,,k])){
        expect_equal(fres_pma_HS_L1_AR0_outAR0$maa[i,j,k], fres_pma_HS_L1_AR0_outAR0_check$maa[i,j,k])

      }
    }
  }

  for( i in 1:nrow(fres_pma_HS_L1_AR0_outAR0_check$vbiom)){
    for( j in 1:ncol(fres_pma_HS_L1_AR0_outAR0_check$vbiom)){
      expect_equal(fres_pma_HS_L1_AR0_outAR0$vbiom[i,j], fres_pma_HS_L1_AR0_outAR0_check$vbiom[i,j])

    }
  }

  for( i in 1:nrow(fres_pma_HS_L1_AR0_outAR0_check$recruit)){
    for( j in 1:ncol(fres_pma_HS_L1_AR0_outAR0_check$recruit)){
      expect_equal(fres_pma_HS_L1_AR0_outAR0$recruit[i,j], fres_pma_HS_L1_AR0_outAR0_check$recruit[i,j])

    }
  }

  for( i in 1:nrow(fres_pma_HS_L1_AR0_outAR0_check$eaa)){
    for( j in 1:ncol(fres_pma_HS_L1_AR0_outAR0_check$eaa)){
      expect_equal(fres_pma_HS_L1_AR0_outAR0$eaa[i,j], fres_pma_HS_L1_AR0_outAR0_check$eaa[i,j])

    }
  }

  for( i in 1:nrow(fres_pma_HS_L1_AR0_outAR0_check$alpha)){
    for( j in 1:ncol(fres_pma_HS_L1_AR0_outAR0_check$alpha)){
      expect_equal(fres_pma_HS_L1_AR0_outAR0$alpha[i,j], fres_pma_HS_L1_AR0_outAR0_check$alpha[i,j])

    }
  }

  for( i in 1:nrow(fres_pma_HS_L1_AR0_outAR0_check$thisyear.ssb)){
    for( j in 1:ncol(fres_pma_HS_L1_AR0_outAR0_check$thisyear.ssb)){
      expect_equal(fres_pma_HS_L1_AR0_outAR0$thisyear.ssb[i,j], fres_pma_HS_L1_AR0_outAR0_check$thisyear.ssb[i,j])

    }
  }

  for( k in 1:length(fres_pma_HS_L1_AR0_outAR0_check$waa[1,1,])){
    for( i in 1:nrow(fres_pma_HS_L1_AR0_outAR0_check$waa[,,k])){
      for( j in 1:ncol(fres_pma_HS_L1_AR0_outAR0_check$waa[,,k])){
        expect_equal(fres_pma_HS_L1_AR0_outAR0$waa[i,j,k], fres_pma_HS_L1_AR0_outAR0_check$waa[i,j,k])

      }
    }
  }

  for( k in 1:length(fres_pma_HS_L1_AR0_outAR0_check$waa.catch[1,1,])){
    for( i in 1:nrow(fres_pma_HS_L1_AR0_outAR0_check$waa.catch[,,k])){
      for( j in 1:ncol(fres_pma_HS_L1_AR0_outAR0_check$waa.catch[,,k])){
        expect_equal(fres_pma_HS_L1_AR0_outAR0$waa.catch[i,j,k], fres_pma_HS_L1_AR0_outAR0_check$waa.catch[i,j,k])

      }
    }
  }

  expect_equal(fres_pma_HS_L1_AR0_outAR0$currentF, fres_pma_HS_L1_AR0_outAR0_check$currentF)

  for( i in 1:nrow(fres_pma_HS_L1_AR0_outAR0_check$vssb)){
    for( j in 1:ncol(fres_pma_HS_L1_AR0_outAR0_check$vssb)){
      expect_equal(fres_pma_HS_L1_AR0_outAR0$vssb[i,j], fres_pma_HS_L1_AR0_outAR0_check$vssb[i,j])

    }
  }

  for( i in 1:nrow(fres_pma_HS_L1_AR0_outAR0_check$vwcaa)){
    for( j in 1:ncol(fres_pma_HS_L1_AR0_outAR0_check$vwcaa)){
      expect_equal(fres_pma_HS_L1_AR0_outAR0$vwcaa[i,j], fres_pma_HS_L1_AR0_outAR0_check$vwcaa[i,j])

    }
  }

  for( k in 1:length(fres_pma_HS_L1_AR0_outAR0_check$naa_all[1,1,])){
    for( i in 1:nrow(fres_pma_HS_L1_AR0_outAR0_check$naa_all[,,k])){
      for( j in 1:ncol(fres_pma_HS_L1_AR0_outAR0_check$naa_all[,,k])){
        expect_equal(fres_pma_HS_L1_AR0_outAR0$naa_all[i,j,k], fres_pma_HS_L1_AR0_outAR0_check$naa_all[i,j,k])

      }
    }
  }

  expect_equal(fres_pma_HS_L1_AR0_outAR0$years, fres_pma_HS_L1_AR0_outAR0_check$years)

  expect_equal(fres_pma_HS_L1_AR0_outAR0$fyear.year, fres_pma_HS_L1_AR0_outAR0_check$fyear.year)

  expect_equal(fres_pma_HS_L1_AR0_outAR0$ABC, fres_pma_HS_L1_AR0_outAR0_check$ABC)

  expect_equal(fres_pma_HS_L1_AR0_outAR0$waa.year, fres_pma_HS_L1_AR0_outAR0_check$waa.year)

  expect_equal(fres_pma_HS_L1_AR0_outAR0$maa.year, fres_pma_HS_L1_AR0_outAR0_check$maa.year)

  expect_equal(fres_pma_HS_L1_AR0_outAR0$multi, fres_pma_HS_L1_AR0_outAR0_check$multi)

  expect_equal(fres_pma_HS_L1_AR0_outAR0$multi.year, fres_pma_HS_L1_AR0_outAR0_check$multi.year)


  # HS L1 AR1 ----
  load(system.file("extdata","fres_pma_HS_L1_AR1_outAR0.rda",package = "frasyr"))
  for( k in 1:length(fres_pma_HS_L1_AR1_outAR0_check$faa[1,1,])){
    for( i in 1:nrow(fres_pma_HS_L1_AR1_outAR0_check$faa[,,k])){
      for( j in 1:ncol(fres_pma_HS_L1_AR1_outAR0_check$faa[,,k])){
        expect_equal(fres_pma_HS_L1_AR1_outAR0$faa[i,j,k], fres_pma_HS_L1_AR1_outAR0_check$faa[i,j,k])

      }
    }
  }

  for( k in 1:length(fres_pma_HS_L1_AR1_outAR0_check$naa[1,1,])){
    for( i in 1:nrow(fres_pma_HS_L1_AR1_outAR0_check$naa[,,k])){
      for( j in 1:ncol(fres_pma_HS_L1_AR1_outAR0_check$naa[,,k])){
        expect_equal(fres_pma_HS_L1_AR1_outAR0$naa[i,j,k], fres_pma_HS_L1_AR1_outAR0_check$naa[i,j,k])

      }
    }
  }

  for( k in 1:length(fres_pma_HS_L1_AR1_outAR0_check$biom[1,1,])){
    for( i in 1:nrow(fres_pma_HS_L1_AR1_outAR0_check$biom[,,k])){
      for( j in 1:ncol(fres_pma_HS_L1_AR1_outAR0_check$biom[,,k])){
        expect_equal(fres_pma_HS_L1_AR1_outAR0$biom[i,j,k], fres_pma_HS_L1_AR1_outAR0_check$biom[i,j,k])

      }
    }
  }

  for( k in 1:length(fres_pma_HS_L1_AR1_outAR0_check$naa[1,1,])){
    for( i in 1:nrow(fres_pma_HS_L1_AR1_outAR0_check$naa[,,k])){
      for( j in 1:ncol(fres_pma_HS_L1_AR1_outAR0_check$baa[,,k])){
        expect_equal(fres_pma_HS_L1_AR1_outAR0$baa[i,j,k], fres_pma_HS_L1_AR1_outAR0_check$baa[i,j,k])

      }
    }
  }

  for( k in 1:length(fres_pma_HS_L1_AR1_outAR0_check$ssb[1,1,])){
    for( i in 1:nrow(fres_pma_HS_L1_AR1_outAR0_check$ssb[,,k])){
      for( j in 1:ncol(fres_pma_HS_L1_AR1_outAR0_check$ssb[,,k])){
        expect_equal(fres_pma_HS_L1_AR1_outAR0$ssb[i,j,k], fres_pma_HS_L1_AR1_outAR0_check$ssb[i,j,k])

      }
    }
  }

  for( k in 1:length(fres_pma_HS_L1_AR1_outAR0_check$wcaa[1,1,])){
    for( i in 1:nrow(fres_pma_HS_L1_AR1_outAR0_check$wcaa[,,k])){
      for( j in 1:ncol(fres_pma_HS_L1_AR1_outAR0_check$wcaa[,,k])){
        expect_equal(fres_pma_HS_L1_AR1_outAR0$wcaa[i,j,k], fres_pma_HS_L1_AR1_outAR0_check$wcaa[i,j,k])

      }
    }
  }

  for( k in 1:length(fres_pma_HS_L1_AR1_outAR0_check$M[1,1,])){
    for( i in 1:nrow(fres_pma_HS_L1_AR1_outAR0_check$M[,,k])){
      for( j in 1:ncol(fres_pma_HS_L1_AR1_outAR0_check$M[,,k])){
        expect_equal(fres_pma_HS_L1_AR1_outAR0$M[i,j,k], fres_pma_HS_L1_AR1_outAR0_check$M[i,j,k])

      }
    }
  }

  for( i in 1:nrow(fres_pma_HS_L1_AR1_outAR0_check$rps)){
    for( j in 1:ncol(fres_pma_HS_L1_AR1_outAR0_check$rps)){
      expect_equal(fres_pma_HS_L1_AR1_outAR0$rps[i,j], fres_pma_HS_L1_AR1_outAR0_check$rps[i,j])

    }
  }

  for( k in 1:length(fres_pma_HS_L1_AR1_outAR0_check$maa[1,1,])){
    for( i in 1:nrow(fres_pma_HS_L1_AR1_outAR0_check$maa[,,k])){
      for( j in 1:ncol(fres_pma_HS_L1_AR1_outAR0_check$maa[,,k])){
        expect_equal(fres_pma_HS_L1_AR1_outAR0$maa[i,j,k], fres_pma_HS_L1_AR1_outAR0_check$maa[i,j,k])

      }
    }
  }

  for( i in 1:nrow(fres_pma_HS_L1_AR1_outAR0_check$vbiom)){
    for( j in 1:ncol(fres_pma_HS_L1_AR1_outAR0_check$vbiom)){
      expect_equal(fres_pma_HS_L1_AR1_outAR0$vbiom[i,j], fres_pma_HS_L1_AR1_outAR0_check$vbiom[i,j])

    }
  }

  for( i in 1:nrow(fres_pma_HS_L1_AR1_outAR0_check$recruit)){
    for( j in 1:ncol(fres_pma_HS_L1_AR1_outAR0_check$recruit)){
      expect_equal(fres_pma_HS_L1_AR1_outAR0$recruit[i,j], fres_pma_HS_L1_AR1_outAR0_check$recruit[i,j])

    }
  }

  for( i in 1:nrow(fres_pma_HS_L1_AR1_outAR0_check$eaa)){
    for( j in 1:ncol(fres_pma_HS_L1_AR1_outAR0_check$eaa)){
      expect_equal(fres_pma_HS_L1_AR1_outAR0$eaa[i,j], fres_pma_HS_L1_AR1_outAR0_check$eaa[i,j])

    }
  }

  for( i in 1:nrow(fres_pma_HS_L1_AR1_outAR0_check$alpha)){
    for( j in 1:ncol(fres_pma_HS_L1_AR1_outAR0_check$alpha)){
      expect_equal(fres_pma_HS_L1_AR1_outAR0$alpha[i,j], fres_pma_HS_L1_AR1_outAR0_check$alpha[i,j])

    }
  }

  for( i in 1:nrow(fres_pma_HS_L1_AR1_outAR0_check$thisyear.ssb)){
    for( j in 1:ncol(fres_pma_HS_L1_AR1_outAR0_check$thisyear.ssb)){
      expect_equal(fres_pma_HS_L1_AR1_outAR0$thisyear.ssb[i,j], fres_pma_HS_L1_AR1_outAR0_check$thisyear.ssb[i,j])

    }
  }

  for( k in 1:length(fres_pma_HS_L1_AR1_outAR0_check$waa[1,1,])){
    for( i in 1:nrow(fres_pma_HS_L1_AR1_outAR0_check$waa[,,k])){
      for( j in 1:ncol(fres_pma_HS_L1_AR1_outAR0_check$waa[,,k])){
        expect_equal(fres_pma_HS_L1_AR1_outAR0$waa[i,j,k], fres_pma_HS_L1_AR1_outAR0_check$waa[i,j,k])

      }
    }
  }

  for( k in 1:length(fres_pma_HS_L1_AR1_outAR0_check$waa.catch[1,1,])){
    for( i in 1:nrow(fres_pma_HS_L1_AR1_outAR0_check$waa.catch[,,k])){
      for( j in 1:ncol(fres_pma_HS_L1_AR1_outAR0_check$waa.catch[,,k])){
        expect_equal(fres_pma_HS_L1_AR1_outAR0$waa.catch[i,j,k], fres_pma_HS_L1_AR1_outAR0_check$waa.catch[i,j,k])

      }
    }
  }

  expect_equal(fres_pma_HS_L1_AR1_outAR0$currentF, fres_pma_HS_L1_AR1_outAR0_check$currentF)

  for( i in 1:nrow(fres_pma_HS_L1_AR1_outAR0_check$vssb)){
    for( j in 1:ncol(fres_pma_HS_L1_AR1_outAR0_check$vssb)){
      expect_equal(fres_pma_HS_L1_AR1_outAR0$vssb[i,j], fres_pma_HS_L1_AR1_outAR0_check$vssb[i,j])

    }
  }

  for( i in 1:nrow(fres_pma_HS_L1_AR1_outAR0_check$vwcaa)){
    for( j in 1:ncol(fres_pma_HS_L1_AR1_outAR0_check$vwcaa)){
      expect_equal(fres_pma_HS_L1_AR1_outAR0$vwcaa[i,j], fres_pma_HS_L1_AR1_outAR0_check$vwcaa[i,j])

    }
  }

  for( k in 1:length(fres_pma_HS_L1_AR1_outAR0_check$naa_all[1,1,])){
    for( i in 1:nrow(fres_pma_HS_L1_AR1_outAR0_check$naa_all[,,k])){
      for( j in 1:ncol(fres_pma_HS_L1_AR1_outAR0_check$naa_all[,,k])){
        expect_equal(fres_pma_HS_L1_AR1_outAR0$naa_all[i,j,k], fres_pma_HS_L1_AR1_outAR0_check$naa_all[i,j,k])

      }
    }
  }

  expect_equal(fres_pma_HS_L1_AR1_outAR0$years, fres_pma_HS_L1_AR1_outAR0_check$years)

  expect_equal(fres_pma_HS_L1_AR1_outAR0$fyear.year, fres_pma_HS_L1_AR1_outAR0_check$fyear.year)

  expect_equal(fres_pma_HS_L1_AR1_outAR0$ABC, fres_pma_HS_L1_AR1_outAR0_check$ABC)

  expect_equal(fres_pma_HS_L1_AR1_outAR0$waa.year, fres_pma_HS_L1_AR1_outAR0_check$waa.year)

  expect_equal(fres_pma_HS_L1_AR1_outAR0$maa.year, fres_pma_HS_L1_AR1_outAR0_check$maa.year)

  expect_equal(fres_pma_HS_L1_AR1_outAR0$multi, fres_pma_HS_L1_AR1_outAR0_check$multi)

  expect_equal(fres_pma_HS_L1_AR1_outAR0$multi.year, fres_pma_HS_L1_AR1_outAR0_check$multi.year)

  # HS L1 AR1 outAR1 ----
  load(system.file("extdata","fres_pma_HS_L1_AR1_outAR1.rda",package = "frasyr"))
  for( k in 1:length(fres_pma_HS_L1_AR1_outAR1_check$faa[1,1,])){
    for( i in 1:nrow(fres_pma_HS_L1_AR1_outAR1_check$faa[,,k])){
      for( j in 1:ncol(fres_pma_HS_L1_AR1_outAR1_check$faa[,,k])){
        expect_equal(fres_pma_HS_L1_AR1_outAR1$faa[i,j,k], fres_pma_HS_L1_AR1_outAR1_check$faa[i,j,k])

      }
    }
  }

  for( k in 1:length(fres_pma_HS_L1_AR1_outAR1_check$naa[1,1,])){
    for( i in 1:nrow(fres_pma_HS_L1_AR1_outAR1_check$naa[,,k])){
      for( j in 1:ncol(fres_pma_HS_L1_AR1_outAR1_check$naa[,,k])){
        expect_equal(fres_pma_HS_L1_AR1_outAR1$naa[i,j,k], fres_pma_HS_L1_AR1_outAR1_check$naa[i,j,k])

      }
    }
  }

  for( k in 1:length(fres_pma_HS_L1_AR1_outAR1_check$biom[1,1,])){
    for( i in 1:nrow(fres_pma_HS_L1_AR1_outAR1_check$biom[,,k])){
      for( j in 1:ncol(fres_pma_HS_L1_AR1_outAR1_check$biom[,,k])){
        expect_equal(fres_pma_HS_L1_AR1_outAR1$biom[i,j,k], fres_pma_HS_L1_AR1_outAR1_check$biom[i,j,k])

      }
    }
  }

  for( k in 1:length(fres_pma_HS_L1_AR1_outAR1_check$naa[1,1,])){
    for( i in 1:nrow(fres_pma_HS_L1_AR1_outAR1_check$naa[,,k])){
      for( j in 1:ncol(fres_pma_HS_L1_AR1_outAR1_check$baa[,,k])){
        expect_equal(fres_pma_HS_L1_AR1_outAR1$baa[i,j,k], fres_pma_HS_L1_AR1_outAR1_check$baa[i,j,k])

      }
    }
  }

  for( k in 1:length(fres_pma_HS_L1_AR1_outAR1_check$ssb[1,1,])){
    for( i in 1:nrow(fres_pma_HS_L1_AR1_outAR1_check$ssb[,,k])){
      for( j in 1:ncol(fres_pma_HS_L1_AR1_outAR1_check$ssb[,,k])){
        expect_equal(fres_pma_HS_L1_AR1_outAR1$ssb[i,j,k], fres_pma_HS_L1_AR1_outAR1_check$ssb[i,j,k])

      }
    }
  }

  for( k in 1:length(fres_pma_HS_L1_AR1_outAR1_check$wcaa[1,1,])){
    for( i in 1:nrow(fres_pma_HS_L1_AR1_outAR1_check$wcaa[,,k])){
      for( j in 1:ncol(fres_pma_HS_L1_AR1_outAR1_check$wcaa[,,k])){
        expect_equal(fres_pma_HS_L1_AR1_outAR1$wcaa[i,j,k], fres_pma_HS_L1_AR1_outAR1_check$wcaa[i,j,k])

      }
    }
  }

  for( k in 1:length(fres_pma_HS_L1_AR1_outAR1_check$M[1,1,])){
    for( i in 1:nrow(fres_pma_HS_L1_AR1_outAR1_check$M[,,k])){
      for( j in 1:ncol(fres_pma_HS_L1_AR1_outAR1_check$M[,,k])){
        expect_equal(fres_pma_HS_L1_AR1_outAR1$M[i,j,k], fres_pma_HS_L1_AR1_outAR1_check$M[i,j,k])

      }
    }
  }

  for( i in 1:nrow(fres_pma_HS_L1_AR1_outAR1_check$rps)){
    for( j in 1:ncol(fres_pma_HS_L1_AR1_outAR1_check$rps)){
      expect_equal(fres_pma_HS_L1_AR1_outAR1$rps[i,j], fres_pma_HS_L1_AR1_outAR1_check$rps[i,j])

    }
  }

  for( k in 1:length(fres_pma_HS_L1_AR1_outAR1_check$maa[1,1,])){
    for( i in 1:nrow(fres_pma_HS_L1_AR1_outAR1_check$maa[,,k])){
      for( j in 1:ncol(fres_pma_HS_L1_AR1_outAR1_check$maa[,,k])){
        expect_equal(fres_pma_HS_L1_AR1_outAR1$maa[i,j,k], fres_pma_HS_L1_AR1_outAR1_check$maa[i,j,k])

      }
    }
  }

  for( i in 1:nrow(fres_pma_HS_L1_AR1_outAR1_check$vbiom)){
    for( j in 1:ncol(fres_pma_HS_L1_AR1_outAR1_check$vbiom)){
      expect_equal(fres_pma_HS_L1_AR1_outAR1$vbiom[i,j], fres_pma_HS_L1_AR1_outAR1_check$vbiom[i,j])

    }
  }

  for( i in 1:nrow(fres_pma_HS_L1_AR1_outAR1_check$recruit)){
    for( j in 1:ncol(fres_pma_HS_L1_AR1_outAR1_check$recruit)){
      expect_equal(fres_pma_HS_L1_AR1_outAR1$recruit[i,j], fres_pma_HS_L1_AR1_outAR1_check$recruit[i,j])

    }
  }

  for( i in 1:nrow(fres_pma_HS_L1_AR1_outAR1_check$eaa)){
    for( j in 1:ncol(fres_pma_HS_L1_AR1_outAR1_check$eaa)){
      expect_equal(fres_pma_HS_L1_AR1_outAR1$eaa[i,j], fres_pma_HS_L1_AR1_outAR1_check$eaa[i,j])

    }
  }

  for( i in 1:nrow(fres_pma_HS_L1_AR1_outAR1_check$alpha)){
    for( j in 1:ncol(fres_pma_HS_L1_AR1_outAR1_check$alpha)){
      expect_equal(fres_pma_HS_L1_AR1_outAR1$alpha[i,j], fres_pma_HS_L1_AR1_outAR1_check$alpha[i,j])

    }
  }

  for( i in 1:nrow(fres_pma_HS_L1_AR1_outAR1_check$thisyear.ssb)){
    for( j in 1:ncol(fres_pma_HS_L1_AR1_outAR1_check$thisyear.ssb)){
      expect_equal(fres_pma_HS_L1_AR1_outAR1$thisyear.ssb[i,j], fres_pma_HS_L1_AR1_outAR1_check$thisyear.ssb[i,j])

    }
  }

  for( k in 1:length(fres_pma_HS_L1_AR1_outAR1_check$waa[1,1,])){
    for( i in 1:nrow(fres_pma_HS_L1_AR1_outAR1_check$waa[,,k])){
      for( j in 1:ncol(fres_pma_HS_L1_AR1_outAR1_check$waa[,,k])){
        expect_equal(fres_pma_HS_L1_AR1_outAR1$waa[i,j,k], fres_pma_HS_L1_AR1_outAR1_check$waa[i,j,k])

      }
    }
  }

  for( k in 1:length(fres_pma_HS_L1_AR1_outAR1_check$waa.catch[1,1,])){
    for( i in 1:nrow(fres_pma_HS_L1_AR1_outAR1_check$waa.catch[,,k])){
      for( j in 1:ncol(fres_pma_HS_L1_AR1_outAR1_check$waa.catch[,,k])){
        expect_equal(fres_pma_HS_L1_AR1_outAR1$waa.catch[i,j,k], fres_pma_HS_L1_AR1_outAR1_check$waa.catch[i,j,k])

      }
    }
  }

  expect_equal(fres_pma_HS_L1_AR1_outAR1$currentF, fres_pma_HS_L1_AR1_outAR1_check$currentF)

  for( i in 1:nrow(fres_pma_HS_L1_AR1_outAR1_check$vssb)){
    for( j in 1:ncol(fres_pma_HS_L1_AR1_outAR1_check$vssb)){
      expect_equal(fres_pma_HS_L1_AR1_outAR1$vssb[i,j], fres_pma_HS_L1_AR1_outAR1_check$vssb[i,j])

    }
  }

  for( i in 1:nrow(fres_pma_HS_L1_AR1_outAR1_check$vwcaa)){
    for( j in 1:ncol(fres_pma_HS_L1_AR1_outAR1_check$vwcaa)){
      expect_equal(fres_pma_HS_L1_AR1_outAR1$vwcaa[i,j], fres_pma_HS_L1_AR1_outAR1_check$vwcaa[i,j])

    }
  }

  for( k in 1:length(fres_pma_HS_L1_AR1_outAR1_check$naa_all[1,1,])){
    for( i in 1:nrow(fres_pma_HS_L1_AR1_outAR1_check$naa_all[,,k])){
      for( j in 1:ncol(fres_pma_HS_L1_AR1_outAR1_check$naa_all[,,k])){
        expect_equal(fres_pma_HS_L1_AR1_outAR1$naa_all[i,j,k], fres_pma_HS_L1_AR1_outAR1_check$naa_all[i,j,k])

      }
    }
  }

  expect_equal(fres_pma_HS_L1_AR1_outAR1$years, fres_pma_HS_L1_AR1_outAR1_check$years)

  expect_equal(fres_pma_HS_L1_AR1_outAR1$fyear.year, fres_pma_HS_L1_AR1_outAR1_check$fyear.year)

  expect_equal(fres_pma_HS_L1_AR1_outAR1$ABC, fres_pma_HS_L1_AR1_outAR1_check$ABC)

  expect_equal(fres_pma_HS_L1_AR1_outAR1$waa.year, fres_pma_HS_L1_AR1_outAR1_check$waa.year)

  expect_equal(fres_pma_HS_L1_AR1_outAR1$maa.year, fres_pma_HS_L1_AR1_outAR1_check$maa.year)

  expect_equal(fres_pma_HS_L1_AR1_outAR1$multi, fres_pma_HS_L1_AR1_outAR1_check$multi)

  expect_equal(fres_pma_HS_L1_AR1_outAR1$multi.year, fres_pma_HS_L1_AR1_outAR1_check$multi.year)



  # HS L2 AR0 ----
  load(system.file("extdata","fres_pma_HS_L2_AR0_outAR0.rda",package = "frasyr"))
  for( k in 1:length(fres_pma_HS_L2_AR0_outAR0_check$faa[1,1,])){
    for( i in 1:nrow(fres_pma_HS_L2_AR0_outAR0_check$faa[,,k])){
      for( j in 1:ncol(fres_pma_HS_L2_AR0_outAR0_check$faa[,,k])){
        expect_equal(fres_pma_HS_L2_AR0_outAR0$faa[i,j,k], fres_pma_HS_L2_AR0_outAR0_check$faa[i,j,k])

      }
    }
  }

  for( k in 1:length(fres_pma_HS_L2_AR0_outAR0_check$naa[1,1,])){
    for( i in 1:nrow(fres_pma_HS_L2_AR0_outAR0_check$naa[,,k])){
      for( j in 1:ncol(fres_pma_HS_L2_AR0_outAR0_check$naa[,,k])){
        expect_equal(fres_pma_HS_L2_AR0_outAR0$naa[i,j,k], fres_pma_HS_L2_AR0_outAR0_check$naa[i,j,k])

      }
    }
  }

  for( k in 1:length(fres_pma_HS_L2_AR0_outAR0_check$biom[1,1,])){
    for( i in 1:nrow(fres_pma_HS_L2_AR0_outAR0_check$biom[,,k])){
      for( j in 1:ncol(fres_pma_HS_L2_AR0_outAR0_check$biom[,,k])){
        expect_equal(fres_pma_HS_L2_AR0_outAR0$biom[i,j,k], fres_pma_HS_L2_AR0_outAR0_check$biom[i,j,k])

      }
    }
  }

  for( k in 1:length(fres_pma_HS_L2_AR0_outAR0_check$naa[1,1,])){
    for( i in 1:nrow(fres_pma_HS_L2_AR0_outAR0_check$naa[,,k])){
      for( j in 1:ncol(fres_pma_HS_L2_AR0_outAR0_check$baa[,,k])){
        expect_equal(fres_pma_HS_L2_AR0_outAR0$baa[i,j,k], fres_pma_HS_L2_AR0_outAR0_check$baa[i,j,k])

      }
    }
  }

  for( k in 1:length(fres_pma_HS_L2_AR0_outAR0_check$ssb[1,1,])){
    for( i in 1:nrow(fres_pma_HS_L2_AR0_outAR0_check$ssb[,,k])){
      for( j in 1:ncol(fres_pma_HS_L2_AR0_outAR0_check$ssb[,,k])){
        expect_equal(fres_pma_HS_L2_AR0_outAR0$ssb[i,j,k], fres_pma_HS_L2_AR0_outAR0_check$ssb[i,j,k])

      }
    }
  }

  for( k in 1:length(fres_pma_HS_L2_AR0_outAR0_check$wcaa[1,1,])){
    for( i in 1:nrow(fres_pma_HS_L2_AR0_outAR0_check$wcaa[,,k])){
      for( j in 1:ncol(fres_pma_HS_L2_AR0_outAR0_check$wcaa[,,k])){
        expect_equal(fres_pma_HS_L2_AR0_outAR0$wcaa[i,j,k], fres_pma_HS_L2_AR0_outAR0_check$wcaa[i,j,k])

      }
    }
  }

  for( k in 1:length(fres_pma_HS_L2_AR0_outAR0_check$M[1,1,])){
    for( i in 1:nrow(fres_pma_HS_L2_AR0_outAR0_check$M[,,k])){
      for( j in 1:ncol(fres_pma_HS_L2_AR0_outAR0_check$M[,,k])){
        expect_equal(fres_pma_HS_L2_AR0_outAR0$M[i,j,k], fres_pma_HS_L2_AR0_outAR0_check$M[i,j,k])

      }
    }
  }

  for( i in 1:nrow(fres_pma_HS_L2_AR0_outAR0_check$rps)){
    for( j in 1:ncol(fres_pma_HS_L2_AR0_outAR0_check$rps)){
      expect_equal(fres_pma_HS_L2_AR0_outAR0$rps[i,j], fres_pma_HS_L2_AR0_outAR0_check$rps[i,j])

    }
  }

  for( k in 1:length(fres_pma_HS_L2_AR0_outAR0_check$maa[1,1,])){
    for( i in 1:nrow(fres_pma_HS_L2_AR0_outAR0_check$maa[,,k])){
      for( j in 1:ncol(fres_pma_HS_L2_AR0_outAR0_check$maa[,,k])){
        expect_equal(fres_pma_HS_L2_AR0_outAR0$maa[i,j,k], fres_pma_HS_L2_AR0_outAR0_check$maa[i,j,k])

      }
    }
  }

  for( i in 1:nrow(fres_pma_HS_L2_AR0_outAR0_check$vbiom)){
    for( j in 1:ncol(fres_pma_HS_L2_AR0_outAR0_check$vbiom)){
      expect_equal(fres_pma_HS_L2_AR0_outAR0$vbiom[i,j], fres_pma_HS_L2_AR0_outAR0_check$vbiom[i,j])

    }
  }

  for( i in 1:nrow(fres_pma_HS_L2_AR0_outAR0_check$recruit)){
    for( j in 1:ncol(fres_pma_HS_L2_AR0_outAR0_check$recruit)){
      expect_equal(fres_pma_HS_L2_AR0_outAR0$recruit[i,j], fres_pma_HS_L2_AR0_outAR0_check$recruit[i,j])

    }
  }

  for( i in 1:nrow(fres_pma_HS_L2_AR0_outAR0_check$eaa)){
    for( j in 1:ncol(fres_pma_HS_L2_AR0_outAR0_check$eaa)){
      expect_equal(fres_pma_HS_L2_AR0_outAR0$eaa[i,j], fres_pma_HS_L2_AR0_outAR0_check$eaa[i,j])

    }
  }

  for( i in 1:nrow(fres_pma_HS_L2_AR0_outAR0_check$alpha)){
    for( j in 1:ncol(fres_pma_HS_L2_AR0_outAR0_check$alpha)){
      expect_equal(fres_pma_HS_L2_AR0_outAR0$alpha[i,j], fres_pma_HS_L2_AR0_outAR0_check$alpha[i,j])

    }
  }

  for( i in 1:nrow(fres_pma_HS_L2_AR0_outAR0_check$thisyear.ssb)){
    for( j in 1:ncol(fres_pma_HS_L2_AR0_outAR0_check$thisyear.ssb)){
      expect_equal(fres_pma_HS_L2_AR0_outAR0$thisyear.ssb[i,j], fres_pma_HS_L2_AR0_outAR0_check$thisyear.ssb[i,j])

    }
  }

  for( k in 1:length(fres_pma_HS_L2_AR0_outAR0_check$waa[1,1,])){
    for( i in 1:nrow(fres_pma_HS_L2_AR0_outAR0_check$waa[,,k])){
      for( j in 1:ncol(fres_pma_HS_L2_AR0_outAR0_check$waa[,,k])){
        expect_equal(fres_pma_HS_L2_AR0_outAR0$waa[i,j,k], fres_pma_HS_L2_AR0_outAR0_check$waa[i,j,k])

      }
    }
  }

  for( k in 1:length(fres_pma_HS_L2_AR0_outAR0_check$waa.catch[1,1,])){
    for( i in 1:nrow(fres_pma_HS_L2_AR0_outAR0_check$waa.catch[,,k])){
      for( j in 1:ncol(fres_pma_HS_L2_AR0_outAR0_check$waa.catch[,,k])){
        expect_equal(fres_pma_HS_L2_AR0_outAR0$waa.catch[i,j,k], fres_pma_HS_L2_AR0_outAR0_check$waa.catch[i,j,k])

      }
    }
  }

  expect_equal(fres_pma_HS_L2_AR0_outAR0$currentF, fres_pma_HS_L2_AR0_outAR0_check$currentF)

  for( i in 1:nrow(fres_pma_HS_L2_AR0_outAR0_check$vssb)){
    for( j in 1:ncol(fres_pma_HS_L2_AR0_outAR0_check$vssb)){
      expect_equal(fres_pma_HS_L2_AR0_outAR0$vssb[i,j], fres_pma_HS_L2_AR0_outAR0_check$vssb[i,j])

    }
  }

  for( i in 1:nrow(fres_pma_HS_L2_AR0_outAR0_check$vwcaa)){
    for( j in 1:ncol(fres_pma_HS_L2_AR0_outAR0_check$vwcaa)){
      expect_equal(fres_pma_HS_L2_AR0_outAR0$vwcaa[i,j], fres_pma_HS_L2_AR0_outAR0_check$vwcaa[i,j])

    }
  }

  for( k in 1:length(fres_pma_HS_L2_AR0_outAR0_check$naa_all[1,1,])){
    for( i in 1:nrow(fres_pma_HS_L2_AR0_outAR0_check$naa_all[,,k])){
      for( j in 1:ncol(fres_pma_HS_L2_AR0_outAR0_check$naa_all[,,k])){
        expect_equal(fres_pma_HS_L2_AR0_outAR0$naa_all[i,j,k], fres_pma_HS_L2_AR0_outAR0_check$naa_all[i,j,k])

      }
    }
  }

  expect_equal(fres_pma_HS_L2_AR0_outAR0$years, fres_pma_HS_L2_AR0_outAR0_check$years)

  expect_equal(fres_pma_HS_L2_AR0_outAR0$fyear.year, fres_pma_HS_L2_AR0_outAR0_check$fyear.year)

  expect_equal(fres_pma_HS_L2_AR0_outAR0$ABC, fres_pma_HS_L2_AR0_outAR0_check$ABC)

  expect_equal(fres_pma_HS_L2_AR0_outAR0$waa.year, fres_pma_HS_L2_AR0_outAR0_check$waa.year)

  expect_equal(fres_pma_HS_L2_AR0_outAR0$maa.year, fres_pma_HS_L2_AR0_outAR0_check$maa.year)

  expect_equal(fres_pma_HS_L2_AR0_outAR0$multi, fres_pma_HS_L2_AR0_outAR0_check$multi)

  expect_equal(fres_pma_HS_L2_AR0_outAR0$multi.year, fres_pma_HS_L2_AR0_outAR0_check$multi.year)


  # HS L2 AR1 ----
  load(system.file("extdata","fres_pma_HS_L2_AR1_outAR0.rda",package = "frasyr"))
  for( k in 1:length(fres_pma_HS_L2_AR1_outAR0_check$faa[1,1,])){
    for( i in 1:nrow(fres_pma_HS_L2_AR1_outAR0_check$faa[,,k])){
      for( j in 1:ncol(fres_pma_HS_L2_AR1_outAR0_check$faa[,,k])){
        expect_equal(fres_pma_HS_L2_AR1_outAR0$faa[i,j,k], fres_pma_HS_L2_AR1_outAR0_check$faa[i,j,k])

      }
    }
  }

  for( k in 1:length(fres_pma_HS_L2_AR1_outAR0_check$naa[1,1,])){
    for( i in 1:nrow(fres_pma_HS_L2_AR1_outAR0_check$naa[,,k])){
      for( j in 1:ncol(fres_pma_HS_L2_AR1_outAR0_check$naa[,,k])){
        expect_equal(fres_pma_HS_L2_AR1_outAR0$naa[i,j,k], fres_pma_HS_L2_AR1_outAR0_check$naa[i,j,k])

      }
    }
  }

  for( k in 1:length(fres_pma_HS_L2_AR1_outAR0_check$biom[1,1,])){
    for( i in 1:nrow(fres_pma_HS_L2_AR1_outAR0_check$biom[,,k])){
      for( j in 1:ncol(fres_pma_HS_L2_AR1_outAR0_check$biom[,,k])){
        expect_equal(fres_pma_HS_L2_AR1_outAR0$biom[i,j,k], fres_pma_HS_L2_AR1_outAR0_check$biom[i,j,k])

      }
    }
  }

  for( k in 1:length(fres_pma_HS_L2_AR1_outAR0_check$naa[1,1,])){
    for( i in 1:nrow(fres_pma_HS_L2_AR1_outAR0_check$naa[,,k])){
      for( j in 1:ncol(fres_pma_HS_L2_AR1_outAR0_check$baa[,,k])){
        expect_equal(fres_pma_HS_L2_AR1_outAR0$baa[i,j,k], fres_pma_HS_L2_AR1_outAR0_check$baa[i,j,k])

      }
    }
  }

  for( k in 1:length(fres_pma_HS_L2_AR1_outAR0_check$ssb[1,1,])){
    for( i in 1:nrow(fres_pma_HS_L2_AR1_outAR0_check$ssb[,,k])){
      for( j in 1:ncol(fres_pma_HS_L2_AR1_outAR0_check$ssb[,,k])){
        expect_equal(fres_pma_HS_L2_AR1_outAR0$ssb[i,j,k], fres_pma_HS_L2_AR1_outAR0_check$ssb[i,j,k])

      }
    }
  }

  for( k in 1:length(fres_pma_HS_L2_AR1_outAR0_check$wcaa[1,1,])){
    for( i in 1:nrow(fres_pma_HS_L2_AR1_outAR0_check$wcaa[,,k])){
      for( j in 1:ncol(fres_pma_HS_L2_AR1_outAR0_check$wcaa[,,k])){
        expect_equal(fres_pma_HS_L2_AR1_outAR0$wcaa[i,j,k], fres_pma_HS_L2_AR1_outAR0_check$wcaa[i,j,k])

      }
    }
  }

  for( k in 1:length(fres_pma_HS_L2_AR1_outAR0_check$M[1,1,])){
    for( i in 1:nrow(fres_pma_HS_L2_AR1_outAR0_check$M[,,k])){
      for( j in 1:ncol(fres_pma_HS_L2_AR1_outAR0_check$M[,,k])){
        expect_equal(fres_pma_HS_L2_AR1_outAR0$M[i,j,k], fres_pma_HS_L2_AR1_outAR0_check$M[i,j,k])

      }
    }
  }

  for( i in 1:nrow(fres_pma_HS_L2_AR1_outAR0_check$rps)){
    for( j in 1:ncol(fres_pma_HS_L2_AR1_outAR0_check$rps)){
      expect_equal(fres_pma_HS_L2_AR1_outAR0$rps[i,j], fres_pma_HS_L2_AR1_outAR0_check$rps[i,j])

    }
  }

  for( k in 1:length(fres_pma_HS_L2_AR1_outAR0_check$maa[1,1,])){
    for( i in 1:nrow(fres_pma_HS_L2_AR1_outAR0_check$maa[,,k])){
      for( j in 1:ncol(fres_pma_HS_L2_AR1_outAR0_check$maa[,,k])){
        expect_equal(fres_pma_HS_L2_AR1_outAR0$maa[i,j,k], fres_pma_HS_L2_AR1_outAR0_check$maa[i,j,k])

      }
    }
  }

  for( i in 1:nrow(fres_pma_HS_L2_AR1_outAR0_check$vbiom)){
    for( j in 1:ncol(fres_pma_HS_L2_AR1_outAR0_check$vbiom)){
      expect_equal(fres_pma_HS_L2_AR1_outAR0$vbiom[i,j], fres_pma_HS_L2_AR1_outAR0_check$vbiom[i,j])

    }
  }

  for( i in 1:nrow(fres_pma_HS_L2_AR1_outAR0_check$recruit)){
    for( j in 1:ncol(fres_pma_HS_L2_AR1_outAR0_check$recruit)){
      expect_equal(fres_pma_HS_L2_AR1_outAR0$recruit[i,j], fres_pma_HS_L2_AR1_outAR0_check$recruit[i,j])

    }
  }

  for( i in 1:nrow(fres_pma_HS_L2_AR1_outAR0_check$eaa)){
    for( j in 1:ncol(fres_pma_HS_L2_AR1_outAR0_check$eaa)){
      expect_equal(fres_pma_HS_L2_AR1_outAR0$eaa[i,j], fres_pma_HS_L2_AR1_outAR0_check$eaa[i,j])

    }
  }

  for( i in 1:nrow(fres_pma_HS_L2_AR1_outAR0_check$alpha)){
    for( j in 1:ncol(fres_pma_HS_L2_AR1_outAR0_check$alpha)){
      expect_equal(fres_pma_HS_L2_AR1_outAR0$alpha[i,j], fres_pma_HS_L2_AR1_outAR0_check$alpha[i,j])

    }
  }

  for( i in 1:nrow(fres_pma_HS_L2_AR1_outAR0_check$thisyear.ssb)){
    for( j in 1:ncol(fres_pma_HS_L2_AR1_outAR0_check$thisyear.ssb)){
      expect_equal(fres_pma_HS_L2_AR1_outAR0$thisyear.ssb[i,j], fres_pma_HS_L2_AR1_outAR0_check$thisyear.ssb[i,j])

    }
  }

  for( k in 1:length(fres_pma_HS_L2_AR1_outAR0_check$waa[1,1,])){
    for( i in 1:nrow(fres_pma_HS_L2_AR1_outAR0_check$waa[,,k])){
      for( j in 1:ncol(fres_pma_HS_L2_AR1_outAR0_check$waa[,,k])){
        expect_equal(fres_pma_HS_L2_AR1_outAR0$waa[i,j,k], fres_pma_HS_L2_AR1_outAR0_check$waa[i,j,k])

      }
    }
  }

  for( k in 1:length(fres_pma_HS_L2_AR1_outAR0_check$waa.catch[1,1,])){
    for( i in 1:nrow(fres_pma_HS_L2_AR1_outAR0_check$waa.catch[,,k])){
      for( j in 1:ncol(fres_pma_HS_L2_AR1_outAR0_check$waa.catch[,,k])){
        expect_equal(fres_pma_HS_L2_AR1_outAR0$waa.catch[i,j,k], fres_pma_HS_L2_AR1_outAR0_check$waa.catch[i,j,k])

      }
    }
  }

  expect_equal(fres_pma_HS_L2_AR1_outAR0$currentF, fres_pma_HS_L2_AR1_outAR0_check$currentF)

  for( i in 1:nrow(fres_pma_HS_L2_AR1_outAR0_check$vssb)){
    for( j in 1:ncol(fres_pma_HS_L2_AR1_outAR0_check$vssb)){
      expect_equal(fres_pma_HS_L2_AR1_outAR0$vssb[i,j], fres_pma_HS_L2_AR1_outAR0_check$vssb[i,j])

    }
  }

  for( i in 1:nrow(fres_pma_HS_L2_AR1_outAR0_check$vwcaa)){
    for( j in 1:ncol(fres_pma_HS_L2_AR1_outAR0_check$vwcaa)){
      expect_equal(fres_pma_HS_L2_AR1_outAR0$vwcaa[i,j], fres_pma_HS_L2_AR1_outAR0_check$vwcaa[i,j])

    }
  }

  for( k in 1:length(fres_pma_HS_L2_AR1_outAR0_check$naa_all[1,1,])){
    for( i in 1:nrow(fres_pma_HS_L2_AR1_outAR0_check$naa_all[,,k])){
      for( j in 1:ncol(fres_pma_HS_L2_AR1_outAR0_check$naa_all[,,k])){
        expect_equal(fres_pma_HS_L2_AR1_outAR0$naa_all[i,j,k], fres_pma_HS_L2_AR1_outAR0_check$naa_all[i,j,k])

      }
    }
  }

  expect_equal(fres_pma_HS_L2_AR1_outAR0$years, fres_pma_HS_L2_AR1_outAR0_check$years)

  expect_equal(fres_pma_HS_L2_AR1_outAR0$fyear.year, fres_pma_HS_L2_AR1_outAR0_check$fyear.year)

  expect_equal(fres_pma_HS_L2_AR1_outAR0$ABC, fres_pma_HS_L2_AR1_outAR0_check$ABC)

  expect_equal(fres_pma_HS_L2_AR1_outAR0$waa.year, fres_pma_HS_L2_AR1_outAR0_check$waa.year)

  expect_equal(fres_pma_HS_L2_AR1_outAR0$maa.year, fres_pma_HS_L2_AR1_outAR0_check$maa.year)

  expect_equal(fres_pma_HS_L2_AR1_outAR0$multi, fres_pma_HS_L2_AR1_outAR0_check$multi)

  expect_equal(fres_pma_HS_L2_AR1_outAR0$multi.year, fres_pma_HS_L2_AR1_outAR0_check$multi.year)

  # HS L2 AR1 outAR1 ----
  load(system.file("extdata","fres_pma_HS_L2_AR1_outAR1.rda",package = "frasyr"))
  for( k in 1:length(fres_pma_HS_L2_AR1_outAR1_check$faa[1,1,])){
    for( i in 1:nrow(fres_pma_HS_L2_AR1_outAR1_check$faa[,,k])){
      for( j in 1:ncol(fres_pma_HS_L2_AR1_outAR1_check$faa[,,k])){
        expect_equal(fres_pma_HS_L2_AR1_outAR1$faa[i,j,k], fres_pma_HS_L2_AR1_outAR1_check$faa[i,j,k])

      }
    }
  }

  for( k in 1:length(fres_pma_HS_L2_AR1_outAR1_check$naa[1,1,])){
    for( i in 1:nrow(fres_pma_HS_L2_AR1_outAR1_check$naa[,,k])){
      for( j in 1:ncol(fres_pma_HS_L2_AR1_outAR1_check$naa[,,k])){
        expect_equal(fres_pma_HS_L2_AR1_outAR1$naa[i,j,k], fres_pma_HS_L2_AR1_outAR1_check$naa[i,j,k])

      }
    }
  }

  for( k in 1:length(fres_pma_HS_L2_AR1_outAR1_check$biom[1,1,])){
    for( i in 1:nrow(fres_pma_HS_L2_AR1_outAR1_check$biom[,,k])){
      for( j in 1:ncol(fres_pma_HS_L2_AR1_outAR1_check$biom[,,k])){
        expect_equal(fres_pma_HS_L2_AR1_outAR1$biom[i,j,k], fres_pma_HS_L2_AR1_outAR1_check$biom[i,j,k])

      }
    }
  }

  for( k in 1:length(fres_pma_HS_L2_AR1_outAR1_check$naa[1,1,])){
    for( i in 1:nrow(fres_pma_HS_L2_AR1_outAR1_check$naa[,,k])){
      for( j in 1:ncol(fres_pma_HS_L2_AR1_outAR1_check$baa[,,k])){
        expect_equal(fres_pma_HS_L2_AR1_outAR1$baa[i,j,k], fres_pma_HS_L2_AR1_outAR1_check$baa[i,j,k])

      }
    }
  }

  for( k in 1:length(fres_pma_HS_L2_AR1_outAR1_check$ssb[1,1,])){
    for( i in 1:nrow(fres_pma_HS_L2_AR1_outAR1_check$ssb[,,k])){
      for( j in 1:ncol(fres_pma_HS_L2_AR1_outAR1_check$ssb[,,k])){
        expect_equal(fres_pma_HS_L2_AR1_outAR1$ssb[i,j,k], fres_pma_HS_L2_AR1_outAR1_check$ssb[i,j,k])

      }
    }
  }

  for( k in 1:length(fres_pma_HS_L2_AR1_outAR1_check$wcaa[1,1,])){
    for( i in 1:nrow(fres_pma_HS_L2_AR1_outAR1_check$wcaa[,,k])){
      for( j in 1:ncol(fres_pma_HS_L2_AR1_outAR1_check$wcaa[,,k])){
        expect_equal(fres_pma_HS_L2_AR1_outAR1$wcaa[i,j,k], fres_pma_HS_L2_AR1_outAR1_check$wcaa[i,j,k])

      }
    }
  }

  for( k in 1:length(fres_pma_HS_L2_AR1_outAR1_check$M[1,1,])){
    for( i in 1:nrow(fres_pma_HS_L2_AR1_outAR1_check$M[,,k])){
      for( j in 1:ncol(fres_pma_HS_L2_AR1_outAR1_check$M[,,k])){
        expect_equal(fres_pma_HS_L2_AR1_outAR1$M[i,j,k], fres_pma_HS_L2_AR1_outAR1_check$M[i,j,k])

      }
    }
  }

  for( i in 1:nrow(fres_pma_HS_L2_AR1_outAR1_check$rps)){
    for( j in 1:ncol(fres_pma_HS_L2_AR1_outAR1_check$rps)){
      expect_equal(fres_pma_HS_L2_AR1_outAR1$rps[i,j], fres_pma_HS_L2_AR1_outAR1_check$rps[i,j])

    }
  }

  for( k in 1:length(fres_pma_HS_L2_AR1_outAR1_check$maa[1,1,])){
    for( i in 1:nrow(fres_pma_HS_L2_AR1_outAR1_check$maa[,,k])){
      for( j in 1:ncol(fres_pma_HS_L2_AR1_outAR1_check$maa[,,k])){
        expect_equal(fres_pma_HS_L2_AR1_outAR1$maa[i,j,k], fres_pma_HS_L2_AR1_outAR1_check$maa[i,j,k])

      }
    }
  }

  for( i in 1:nrow(fres_pma_HS_L2_AR1_outAR1_check$vbiom)){
    for( j in 1:ncol(fres_pma_HS_L2_AR1_outAR1_check$vbiom)){
      expect_equal(fres_pma_HS_L2_AR1_outAR1$vbiom[i,j], fres_pma_HS_L2_AR1_outAR1_check$vbiom[i,j])

    }
  }

  for( i in 1:nrow(fres_pma_HS_L2_AR1_outAR1_check$recruit)){
    for( j in 1:ncol(fres_pma_HS_L2_AR1_outAR1_check$recruit)){
      expect_equal(fres_pma_HS_L2_AR1_outAR1$recruit[i,j], fres_pma_HS_L2_AR1_outAR1_check$recruit[i,j])

    }
  }

  for( i in 1:nrow(fres_pma_HS_L2_AR1_outAR1_check$eaa)){
    for( j in 1:ncol(fres_pma_HS_L2_AR1_outAR1_check$eaa)){
      expect_equal(fres_pma_HS_L2_AR1_outAR1$eaa[i,j], fres_pma_HS_L2_AR1_outAR1_check$eaa[i,j])

    }
  }

  for( i in 1:nrow(fres_pma_HS_L2_AR1_outAR1_check$alpha)){
    for( j in 1:ncol(fres_pma_HS_L2_AR1_outAR1_check$alpha)){
      expect_equal(fres_pma_HS_L2_AR1_outAR1$alpha[i,j], fres_pma_HS_L2_AR1_outAR1_check$alpha[i,j])

    }
  }

  for( i in 1:nrow(fres_pma_HS_L2_AR1_outAR1_check$thisyear.ssb)){
    for( j in 1:ncol(fres_pma_HS_L2_AR1_outAR1_check$thisyear.ssb)){
      expect_equal(fres_pma_HS_L2_AR1_outAR1$thisyear.ssb[i,j], fres_pma_HS_L2_AR1_outAR1_check$thisyear.ssb[i,j])

    }
  }

  for( k in 1:length(fres_pma_HS_L2_AR1_outAR1_check$waa[1,1,])){
    for( i in 1:nrow(fres_pma_HS_L2_AR1_outAR1_check$waa[,,k])){
      for( j in 1:ncol(fres_pma_HS_L2_AR1_outAR1_check$waa[,,k])){
        expect_equal(fres_pma_HS_L2_AR1_outAR1$waa[i,j,k], fres_pma_HS_L2_AR1_outAR1_check$waa[i,j,k])

      }
    }
  }

  for( k in 1:length(fres_pma_HS_L2_AR1_outAR1_check$waa.catch[1,1,])){
    for( i in 1:nrow(fres_pma_HS_L2_AR1_outAR1_check$waa.catch[,,k])){
      for( j in 1:ncol(fres_pma_HS_L2_AR1_outAR1_check$waa.catch[,,k])){
        expect_equal(fres_pma_HS_L2_AR1_outAR1$waa.catch[i,j,k], fres_pma_HS_L2_AR1_outAR1_check$waa.catch[i,j,k])

      }
    }
  }

  expect_equal(fres_pma_HS_L2_AR1_outAR1$currentF, fres_pma_HS_L2_AR1_outAR1_check$currentF)

  for( i in 1:nrow(fres_pma_HS_L2_AR1_outAR1_check$vssb)){
    for( j in 1:ncol(fres_pma_HS_L2_AR1_outAR1_check$vssb)){
      expect_equal(fres_pma_HS_L2_AR1_outAR1$vssb[i,j], fres_pma_HS_L2_AR1_outAR1_check$vssb[i,j])

    }
  }

  for( i in 1:nrow(fres_pma_HS_L2_AR1_outAR1_check$vwcaa)){
    for( j in 1:ncol(fres_pma_HS_L2_AR1_outAR1_check$vwcaa)){
      expect_equal(fres_pma_HS_L2_AR1_outAR1$vwcaa[i,j], fres_pma_HS_L2_AR1_outAR1_check$vwcaa[i,j])

    }
  }

  for( k in 1:length(fres_pma_HS_L2_AR1_outAR1_check$naa_all[1,1,])){
    for( i in 1:nrow(fres_pma_HS_L2_AR1_outAR1_check$naa_all[,,k])){
      for( j in 1:ncol(fres_pma_HS_L2_AR1_outAR1_check$naa_all[,,k])){
        expect_equal(fres_pma_HS_L2_AR1_outAR1$naa_all[i,j,k], fres_pma_HS_L2_AR1_outAR1_check$naa_all[i,j,k])

      }
    }
  }

  expect_equal(fres_pma_HS_L2_AR1_outAR1$years, fres_pma_HS_L2_AR1_outAR1_check$years)

  expect_equal(fres_pma_HS_L2_AR1_outAR1$fyear.year, fres_pma_HS_L2_AR1_outAR1_check$fyear.year)

  expect_equal(fres_pma_HS_L2_AR1_outAR1$ABC, fres_pma_HS_L2_AR1_outAR1_check$ABC)

  expect_equal(fres_pma_HS_L2_AR1_outAR1$waa.year, fres_pma_HS_L2_AR1_outAR1_check$waa.year)

  expect_equal(fres_pma_HS_L2_AR1_outAR1$maa.year, fres_pma_HS_L2_AR1_outAR1_check$maa.year)

  expect_equal(fres_pma_HS_L2_AR1_outAR1$multi, fres_pma_HS_L2_AR1_outAR1_check$multi)

  expect_equal(fres_pma_HS_L2_AR1_outAR1$multi.year, fres_pma_HS_L2_AR1_outAR1_check$multi.year)

   })

context("future est MSY")

test_that("oututput value check",{
  load(system.file("extdata","res_vpa_pma.rda",package = "frasyr"))
  load(system.file("extdata","SRdata_pma.rda",package = "frasyr"))

  SRmodel.list <- expand.grid(SR.rel = c("HS","BH","RI"), AR.type = c(0, 1), out.AR=c(TRUE,FALSE), L.type = c("L1", "L2"))
  SR.list <- list()

  for (i in 1:nrow(SRmodel.list)) {
    SR.list[[i]] <- fit.SR(SRdata_pma, SR = SRmodel.list$SR.rel[i], method = SRmodel.list$L.type[i],
                           AR = SRmodel.list$AR.type[i], out.AR =SRmodel.list$out.AR[i], hessian = FALSE)
  }

  SRmodel.list$AICc <- sapply(SR.list, function(x) x$AICc)
  SRmodel.list$delta.AIC <- SRmodel.list$AICc - min(SRmodel.list$AICc)
  SR.list <- SR.list[order(SRmodel.list$AICc)]  # AICの小さい順に並べたもの
  (SRmodel.list <- SRmodel.list[order(SRmodel.list$AICc), ]) # 結果

  SRmodel.base <- SR.list[[1]] # AIC最小モデルを今後使っていく

  selectedSR <- sprintf("%s.recAR",SRmodel.base$input$SR[1])
 # future fcurrent ----
  res_future_Fcurrent_pma <- future.vpa(res_vpa_pma,
                                        multi=1,
                                        nyear=50, # 将来予測の年数
                                        start.year=2012, # 将来予測の開始年
                                        N=100, # 確率的計算の繰り返し回数
                                        ABC.year=2013, # ABCを計算する年
                                        waa.year=2009:2011, # 生物パラメータの参照年
                                        maa.year=2009:2011,
                                        M.year=2009:2011,
                                        is.plot=TRUE, # 結果をプロットするかどうか
                                        seed=1,
                                        silent=TRUE,
                                        recfunc=eval(parse(text=selectedSR))
                                        , # 再生産関係の関数
                                        # recfuncに対する引数
                                        rec.arg=list(a=SRmodel.base$pars$a,b=SRmodel.base$pars$b,
                                                     rho=SRmodel.base$pars$rho, # ここではrho=0なので指定しなくてもOK
                                                     sd=SRmodel.base$pars$sd,resid=SRmodel.base$resid))

  # est MSY ----
  res_MSY_pma_check <- est.MSY(res_vpa_pma, # VPAの計算結果
                         res_future_Fcurrent_pma$input, # 将来予測で使用した引数
                         seed=res_future_Fcurrent_pma$input$seed,
                         resid.year=0, # ARありの場合、最近何年分の残差を平均するかをここで指定する。ARありの設定を反映させたい場合必ずここを１以上とすること（とりあえず１としておいてください）。
                         N=100, # 確率的計算の繰り返し回数=>実際の計算では1000~5000回くらいやってください
                         calc.yieldcurve=TRUE,
                         PGY=c(0.95,0.9,0.6,0.1), # 計算したいPGYレベル。上限と下限の両方が計算される
                         onlylower.pgy=FALSE, # TRUEにするとPGYレベルの上限は計算しない（計算時間の節約になる）
                         B0percent=c(0.2,0.3,0.4),
                         Bempirical=c(round(tail(colSums(res_vpa_pma$ssb),n=1)),
                                      round(max(colSums(res_vpa_pma$ssb))),
                                      24000, # 現行Blimit
                                      SRmodel.base$pars$b) # HSの折れ点
  )

  # 上記設定の結果を読み込み ----
  load(system.file("extdata","res_MSY_pma.rda",package = "frasyr"))
  # 結果の比較----
  expect_equal(res_MSY_pma$summary,res_MSY_pma_check$summary)
  expect_equal(res_MSY_pma$summaryAR,res_MSY_pma_check$summaryAR)
  expect_equal(res_MSY_pma$F.msy,res_MSY_pma_check$F.msy)
  expect_equal(res_MSY_pma$all.stat,res_MSY_pma_check$all.stat)
  expect_equal(res_MSY_pma$all.statAR,res_MSY_pma_check$all.statAR)
  expect_equal(res_MSY_pma$all.stat,res_MSY_pma_check$all.stat)
  expect_equal(res_MSY_pma$ssb.ar.mean,res_MSY_pma_check$ssb.ar.mean)
  expect_equal(res_MSY_pma$SPR.msy,res_MSY_pma_check$SPR.msy)

 })

context("future HCR")

test_that("oututput value check",{
  load(system.file("extdata","res_MSY_pma.rda",package = "frasyr"))
  refs_all_pma <- res_MSY_pma$summary

  refs_all_pma$RP.definition[refs_all_pma$RP_name=="B0-20%" & refs_all_pma$AR==FALSE] <- "Btarget1"  # たとえばBtargetの代替値としてB020%も候補に残しておきたい場合
  refs_all_pma$RP.definition[refs_all_pma$RP_name=="PGY_0.95_lower" & refs_all_pma$AR==FALSE] <- "Btarget2"
  refs_all_pma$RP.definition[refs_all_pma$RP_name=="Ben-19431" & refs_all_pma$AR==FALSE] <- "Bcurrent"
  refs_all_pma$RP.definition[refs_all_pma$RP_name=="Ben-63967" & refs_all_pma$AR==FALSE] <- "Bmax"
  refs_all_pma$RP.definition[refs_all_pma$RP_name=="Ben-24000" & refs_all_pma$AR==FALSE] <- "Blimit1"
  refs_all_pma$RP.definition[refs_all_pma$RP_name=="Ben-51882" & refs_all_pma$AR==FALSE] <- "B_HS"

  refs_all_pma %>% dplyr::select(RP_name,RP.definition)

  refs_base_pma_check <- refs_all_pma %>%
    dplyr::filter(!is.na(RP.definition)) %>% # RP.definitionがNAでないものを抽出
    arrange(desc(SSB)) %>% # SSBを大きい順に並び替え
    dplyr::select(RP.definition,RP_name,SSB,SSB2SSB0,Catch,Catch.CV,U,Fref2Fcurrent)

  load(system.file("extdata","refs_base_pma.rda",package = "frasyr"))

  expect_equal(refs_base_pma,refs_base_pma_check)

})


context("future future.vpa (option of futureF)")

test_that("oututput value check (iteration for future sim is fixed as 2) ",{
  caa <- read.csv(system.file("extdata","caa_pma.csv",package="frasyr"),row.names=1)
  waa <- read.csv(system.file("extdata","waa_pma.csv",package="frasyr"),row.names=1)
  maa <- read.csv(system.file("extdata","maa_pma.csv",package="frasyr"),row.names=1)
  dat <- data.handler(caa=caa, waa=waa, maa=maa, M=0.5)
  res.pma <- vpa(dat,fc.year=2009:2011,rec=585,rec.year=2011,tf.year = 2008:2010,
                 term.F="max",stat.tf="mean",Pope=TRUE,tune=FALSE,p.init=1.0)
  SRdata <- get.SRdata(res.pma)

  HS.par0 <- fit.SR(SRdata,SR="HS",method="L2",AR=0,hessian=FALSE)
  HS.par1 <- fit.SR(SRdata,SR="HS",method="L2",AR=1,hessian=FALSE)

  currentF.test <- 1:4/10
  futureF.test <- 5:8/10
  fres.HS.check <- future.vpa(res.pma,
                              multi=2,
                              currentF=currentF.test,
                              futureF=5:8/10,
                              nyear=50,
                              start.year=2012,
                              N=2, ABC.year=2013,
                        waa.year=2009:2011,
                        maa.year=2009:2011,
                        M.year=2009:2011,
                        is.plot=TRUE,
                        seed=1,
                        silent=TRUE,
                        recfunc=HS.recAR,
                        rec.arg=list(a=HS.par0$pars$a,b=HS.par0$pars$b,
                                     rho=HS.par0$pars$rho,
                                     sd=HS.par0$pars$sd,resid=HS.par0$resid)
                        )


  # faaが想定通りに入っていればOKくらいのテストです
  for(i in 1:4){
      expect_true(mean(fres.HS.check$faa[i,dimnames(fres.HS.check$faa)[[2]]==2012,])==currentF.test[i])
      expect_true(mean(fres.HS.check$faa[i,dimnames(fres.HS.check$faa)[[2]]>2012,])==futureF.test[i]*2)
  }
})


