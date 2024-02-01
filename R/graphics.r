# 結果のプロットのための関数 ----

# 色や線の太さの設定・会議用の図の設定 ----

col.SBtarget    <- "#00533E"
col.SBlim       <- "#edb918"
col.SBlimit     <- "#edb918"
col.SBban       <- "#C73C2E"
col.Ftarget     <- "#714C99"
col.betaFtarget <- "#505596"
pt1             <- 0.3528

#' 会議用の図のフォーマット
#'
#' @export
#' @encoding UTF-8

theme_SH <- function(legend.position="none",base_size=12){

  ## if(isTRUE(stringr::str_detect(version$os, pattern="darwin"))){
  ##   font_MAC <- "HiraginoSans-W3"
  ##   theme_bw(base_size=base_size) +
  ##     theme(panel.grid = element_blank(),
  ##           axis.text.x=element_text(size=11,color="black"),
  ##           axis.text.y=element_text(size=11,color="black"),
  ##           axis.line.x=element_line(size= 0.3528),
  ##           axis.line.y=element_line(size= 0.3528),
  ##           legend.position=legend.position, text =element_text(family = font_MAC) )
  ## }
  ## else{
    theme_bw(base_size=base_size) +
      theme(panel.grid = element_blank(),
            axis.text.x=element_text(size=11,color="black"),
            axis.text.y=element_text(size=11,color="black"),
            axis.line.x=element_line(size= 0.3528),
            axis.line.y=element_line(size= 0.3528),
            legend.position=legend.position)

#  }
}

#' 会議用の図の出力関数（大きさ・サイズの指定済）：通常サイズ
#'
#' @export
#'

ggsave_SH <- function(...){
  ggsave(width=150,height=85,dpi=600,units="mm",...)
}

#' 会議用の図の出力関数（大きさ・サイズの指定済）：大きいサイズ
#'
#' @export
#' @encoding UTF-8
#'

ggsave_SH_large <- function(...){
  ggsave(width=150,height=120,dpi=600,units="mm",...)
}

# VPA関数用 ----
#' 複数のVPAの結果を重ね書きする
#'
#' @param vpalist vpaの返り値をリストにしたもの; 単独でも可
#' @param vpatibble tibble形式のVPA結果も可。この場合、convert_vpa_tibble関数の出力に準じる。複数のVPA結果がrbindされている場合は、列名"id"で区別する。
#' @param what.plot 何の値をプロットするか. NULLの場合、全て（SSB, biomass, U, catch, Recruitment, fish_number, fishing_mortality, weight, maturity, catch_number）をプロットする。what.plot=c("SSB","Recruitment")とすると一部のみプロット。ここで指定された順番にプロットされる。
#' @param legend.position 凡例の位置。"top" (上部), "bottom" (下部), "left" (左), "right" (右), "none" (なし)。
#' @param vpaname 凡例につけるVPAのケースの名前。vpalistと同じ長さにしないとエラーとなる
#' @param ncol 図の列数
#'
#' @examples
#' \dontrun{
#' data(res_vpa)
#' res_vpa2 <- res_vpa
#' res_vpa2$naa <- res_vpa2$naa*1.2
#'
#' plot_vpa(list(res_vpa, res_vpa2), vpaname=c("True","Dummy"))
#' plot_vpa(list(res_vpa, res_vpa2), vpaname=c("True","Dummy"),
#'                  what.plot=c("SSB","fish_number"))
#'
#' }
#'
#' @encoding UTF-8
#'
#' @export
#'

plot_vpa <- function(vpalist,
                     vpatibble = NULL,
                     what.plot = NULL,
                     plot_year = NULL,  # 浜辺加筆(2020/06/30)
                     legend.position = "top",
                     vpaname = NULL,
                     ncol = 2,
                     scale_value = NULL # 浜辺加筆(2020/07/29)
                     ){

  if(!is.null(vpaname)){
    if(length(vpaname)!=length(vpalist)) stop("Length of vpalist and vpaname is different")
    names(vpalist) <- vpaname
  }

  if(is.null(vpatibble)){
    if(isTRUE("naa" %in% names(vpalist))) vpalist <- list(vpalist)
    vpadata <- vpalist %>% purrr::map_dfr(convert_vpa_tibble ,.id="id") %>%
      mutate(age=factor(age))
  }
  else{
    vpadata <- vpatibble %>%
      mutate(age=factor(age))
    if("id" %in% !names(vpadata)) vpadata$id <- "vpa1"
  }
  if(!is.null(what.plot)) vpadata <- vpadata %>%  dplyr::filter(stat%in%what.plot)

  biomass_factor <- vpadata %>% dplyr::filter(is.na(age)) %>%
    select(stat) %>% unique() %>% unlist()
  age_factor <- vpadata %>% dplyr::filter(!is.na(age)) %>%
    select(stat) %>% unique() %>% unlist()

  if(!is.null(what.plot)){
    vpadata <- vpadata %>%
      mutate(stat=factor(stat,levels=what.plot))
  }
  else{
    vpadata <- vpadata %>%
      mutate(stat=factor(stat,levels=c(biomass_factor, age_factor)))
  }

  # scale_shape_manual関数のvalueを設定(浜辺'20/07/29)
  if(is.null(scale_value))scale_value <- 1:(length(unique(vpadata$id)))
  if(!(length(scale_value)==length(unique(vpadata$id)))){
    scale_value <- 1:(length(unique(vpadata$id)))
    print("Note! Length of scale value was different although plot was done...")
  }
  # シナリオの違いを線種+shapeにした(浜辺'20/06/30)
  g1 <- vpadata %>% ggplot()
  if(all(is.na(vpadata$age))){
    g1 <- g1+ geom_line(aes(x=year, y=value,lty=id))
    g1 <- try(g1 + geom_point(aes(x=year, y=value, shape=id)) +
                scale_shape_manual(values = scale_value))
  }
  else{
    g1 <- g1+ geom_line(aes(x=year, y=value, color=age, lty=id))
    g1 <- try(g1 + geom_point(aes(x=year, y=value, color=age, shape=id)) +
                scale_shape_manual(values = scale_value))
  }
  # 上のがダメな場合にオリジナルで対応
  if(class(g1[1])=="try-error"){ # g1の長さは2あって警告が頻発するので[1]を加えた
    g1 <- vpadata %>% ggplot()
    if(all(is.na(vpadata$age))){
      g1 <- g1+ geom_line(aes(x=year, y=value,lty=id))
    }
    else{
      g1 <- g1+ geom_line(aes(x=year, y=value,color=age,lty=id))
    }
  }

  g1 <- g1 +
    facet_wrap(~stat, scale="free_y", ncol=ncol) + ylim(0,NA) +
    theme_SH(legend.position=legend.position) +
    ylab("value") + xlab("Year")+
    guides(color=guide_legend(nrow=2))

  if(!(is.null(plot_year))){
    g2 <- g1 + xlim(plot_year[1], max(plot_year))
  } else {
    g2 <- g1
  } # もし動かない場合はこれまで通りに作図(浜辺'20/06/30)

  g2
}

#' F currentをプロットする
#'
#' @param vpares VPAの結果のオブジェクト
#' @encoding UTF-8
#'
#' @export

plot_Fcurrent <- function(vpares,
                          Fcurrent=NULL,
                          year.range=NULL){

  if(is.null(year.range)) year.range <- min(as.numeric(colnames(vpares$naa))):max(as.numeric(colnames(vpares$naa)))
  vpares_tb <- convert_vpa_tibble(vpares)

  faa_history <- vpares_tb %>%
    dplyr::filter(stat=="fishing_mortality", year%in%year.range) %>%
    mutate(F=value,year=as.character(year),type="History") %>%
    select(-stat,-sim,-value) %>%
    group_by(year) %>%
    dplyr::filter(!is.na(F)) %>%
    mutate(Year=as.numeric(year)) %>%
    arrange(age) %>%
    mutate(age_name=ifelse(max(age)==age,str_c(age,"+"),age)) %>%
    mutate(age_name=fct_inorder(age_name))

  if(is.null(Fcurrent)){
    fc_at_age_current <- vpares$Fc.at.age
  }
  else{
    fc_at_age_current <- Fcurrent
  }
  fc_at_age_current <- tibble(F=fc_at_age_current,age=as.numeric(rownames(vpares$faa)),
                              year="0",type="currentF") %>%
      dplyr::filter(!is.na(F)) %>%
      arrange(age) %>%
      mutate(age_name=ifelse(max(age)==age,str_c(age,"+"),age))%>%
      mutate(age_name=fct_inorder(age_name))

  pal <- c("#3B9AB2", "#56A6BA", "#71B3C2", "#9EBE91", "#D1C74C",
           "#E8C520", "#E4B80E", "#E29E00", "#EA5C00", "#F21A00")

    g <- faa_history  %>% ggplot() +
    geom_path(aes(x=age_name,y=F,color=Year,group=Year),lwd=1.5) +
    scale_color_gradientn(colors = pal)+
    geom_path(data=fc_at_age_current,
              mapping=aes(x=age_name,y=F,group=type),color="black",lwd=1.5,lty=1)+
    geom_point(data=fc_at_age_current,
               mapping=aes(x=age_name,y=F,shape=type),color="black",size=3)+
    coord_cartesian(ylim=c(0,max(faa_history$F,na.rm=T)))+
    ##                        xlim=range(as.numeric(faa_history$age_name,na.rm=T))+c(-0.5,0.5)
    #                        )+
    xlab("年齢")+ylab("漁獲係数(F)")+theme_SH(legend.position="right")+
    scale_shape_discrete(name="", labels=c("F current"))

  return(g)
}

# ref.F用 ----

#' ref.Fの出力をプロットするための関数
#'
#' @param rres ref.Fの出力結果
#' @encoding UTF-8
#'
#' @export

plot_Fref <- function(rres,xlabel="max", # or, "mean","Fref/Fcur"
                      vline.text=c("F0.1","Fmax","Fcurrent","Fmed") # and "FpSPR.20.SPR" etc..
){
  #old.par <- par()
  #par(mar=c(4,4,1,4))
  F.range <- rres$ypr.spr$F.range
  if(xlabel=="Fref/Fcur") F.range <- F.range/rres$summary$Fcurrent[1]*rres$summary$Fcurrent[3]
  if(xlabel=="mean") F.range <- F.range/rres$summary$Fcurrent[1]*rres$summary$Fcurrent[2]
  spr <- rres$ypr.spr$pspr
  ypr <- rres$ypr.spr$ypr
  plot(F.range,spr,xlab=xlabel,ylab="%SPR",type="l",ylim=c(0,max(spr)))
  par(new=T)
  plot(F.range,ypr,axes=F,xlab="",ylab="",lty=2,type="l",ylim=c(0,max(ypr)))
  axis(side=4)
  mtext("YPR",side=4,line=2)
  n.line <- which(rownames(rres$summary) %in% xlabel)
  abline(v=xx <- c(rres$summary[vline.text][n.line,]))
  text(xx,max(ypr)*seq(from=0.5,to=0.3,length=length(vline.text)),vline.text)
  legend("topright",lty=1:2,legend=c("SPR","YPR"))

  invisible(data.frame(F.range=F.range,spr=spr,ypr=ypr))
  #old.par[c("cin","cxy","csi","cra","csy","din","page")] <- NULL
  #par(old.par)
}

# 再生産関係推定 ----
#' SRdataをプロットする
#'
#' @param vpares VPAの結果のオブジェクト
#' @param frag   MSY計算時の将来予測で用いる引数のリスト
#' @param type "classic"(通常プロット) or "gg"(ggplot)
#'
#' @encoding UTF-8
#'
#'
#' @export
#'
#'

plot_SRdata <- function(SRdata, type=c("classic","gg")[1]){

  is_release_data <- "release" %in% str_sub(names(SRdata),1,7)
  if(is_release_data){
    if(!"release_alive" %in% names(SRdata)) SRdata$release_alive <- SRdata$release
    allR <- SRdata$allR <- SRdata$R + SRdata$release_alive
  }
  else{
    allR <- SRdata$R
  }
  if(type=="classic"){
    plot(SRdata$SSB,SRdata$R,xlab="SSB",ylab="R",xlim=c(0,max(SRdata$SSB)),ylim=c(0,max(allR)))
    if(is_release_data)  points(SRdata$SSB,allR,col=2,pch=3)
  }
  if(type=="gg"){
    if(!"id"%in%names(SRdata)){
      g0 <- ggplot2::qplot(y=R,x=SSB,data=as_tibble(SRdata),
                           xlab="SSB",ylab="R",xlim=c(0,NA),ylim=c(0,NA)) + theme_SH()
      if(is_release_data){
        g0 <- g0 + geom_point(aes(y=allR,x=SSB), color=2, shape=3)
      }
    }
    else{
      g0 <- ggplot(SRdata)+
        geom_point(aes(y=R,x=SSB,color=id),
                    xlab="SSB",ylab="R",
                    xlim=c(0,NA),ylim=c(0,NA)) +
        theme_SH()+
        theme(legend.position="top")
      if(is_release_data){
        g0 <- g0 + geom_point(aes(y=allR,x=SSB, color=id),shape=3)
      }
    }
    return(g0)
  }
}


#' 再生産関係をプロットする関数
#'
#' @param SR_result fit.SRの結果のオブジェクト
#' @param refs 管理基準値 (list(Blimit=0, Bmsy=10, Bban=0))
#' @param
#' @param
#' @param
#' @param
#'
#' @encoding UTF-8
#'
#' @examples
#' \dontrun{
#'   data(res_sr_HSL1)
#'   plot_SR(res_sr_HSL1)
#'   plot_SR(res_sr_HSL1, refs=list(Blimit=20000, Bmsy=60000, Bban=0),
#'           recruit_intercept=100, plot_CI=TRUE)
#'
#'   SRdata2 <- res_sr_HSL1$input$SRdata %>% as.data.frame()
#'   SRdata2$weight <- 1
#'   SRdata2$release <- 100 # 放流データがある場合
#'   SRdata2$weight[length(SRdata2$weight)] <- 0 # 最後の年をフィットに使わない設定
#'   res_SR2 <- fit.SR(SRdata=SRdata2, method="L1", AR=0)
#'   plot_SR(res_SR2)
#' }
#'
#' @export
#'

plot_SR <- function(SR_result,refs=NULL,xscale=1000,xlabel="千トン",yscale=1,ylabel="尾",
                    labeling.year=NULL,add.info=TRUE, recruit_intercept=0,
                    plot_CI=FALSE, CI=0.9, shape_custom=c(21,3),box.padding=0,
                    add_graph=NULL){
  font_MAC <- "HiraginoSans-W3"#"Japan1GothicBBB"#

  if(is.null(refs$Blimit) && !is.null(refs$Blim)) refs$Blimit <- refs$Blim

  if (SR_result$input$SR=="HS") SRF <- function(SSB,a,b,recruit_intercept=0) (ifelse(SSB*xscale>b,b*a,SSB*xscale*a)+recruit_intercept)/yscale
  if (SR_result$input$SR=="BH") SRF <- function(SSB,a,b,recruit_intercept=0) (a*SSB*xscale/(1+b*SSB*xscale)+recruit_intercept)/yscale
  if (SR_result$input$SR=="RI") SRF <- function(SSB,a,b,recruit_intercept=0) (a*SSB*xscale*exp(-b*SSB*xscale)+recruit_intercept)/yscale
  if (SR_result$input$SR=="Mesnil") SRF <- function(SSB,a,b,gamma) (0.5*a*(SSB*xscale+sqrt(b^2+gamma^2/4)-sqrt((SSB*xscale-b)^2+gamma^2/4))+recruit_intercept)/yscale
  if (SR_result$input$SR=="Shepherd") SRF <- function(SSB,a,b,gamma,recruit_intercept=0) (a*SSB*xscale/(1+(b*SSB*xscale)^gamma)+recruit_intercept)/yscale
  if (SR_result$input$SR=="Cushing") SRF <- function(SSB,a,b,recruit_intercept=0) (a*(SSB*xscale)^b+recruit_intercept)/yscale
  if (SR_result$input$SR=="BHS") SRF <- function(SSB,a,b,gamma){
    SSB <- SSB*xscale
    res <- ifelse(SSB < b, a*b*(SSB/b)^{1-(SSB/b)^gamma}, a*b)
    return(res/yscale)    
}


  SRF_CI <- function(CI,sigma,sign,...){
    exp(log(SRF(...))+qnorm(1-(1-CI)/2)*sigma*sign)
  }

  SRdata <- as_tibble(SR_result$input$SRdata) %>%
      mutate(type="obs")
  if(is.null(SRdata$weight)) SRdata$weight <- SR_result$input$w
  SRdata <- SRdata %>%
    mutate(weight=ifelse(weight>0, 1, 0)) %>%
    mutate(weight=factor(weight,levels=c("1","0")))
  SRdata.pred <- as_tibble(SR_result$pred) %>%
    mutate(type="pred", year=NA, R=R)

  is_release_data <- "release" %in% str_sub(names(SRdata),1,7)
  if(is_release_data){
    if(!"release_alive" %in% names(SRdata)){
      SR_result$input$SRdata$release_alive <- SRdata$release
      SR_result$input$SRdata$release       <- NULL
    }
    else{
      SR_result$input$SRdata$release_all <- NULL
      SR_result$input$SRdata$release_aliverate <- NULL
    }
    SRdata.release <- SR_result$input$SRdata %>% mutate(allR=release_alive+R) %>%
      mutate(weight=ifelse(weight>0, 1, 0)) %>%
      select(-R, -release_alive) %>% rename(R=allR) %>% mutate(type="release", weight=factor(weight, levels=c("0","1")))
  }
  else{
    SRdata.release <- NULL
  }
  alldata <- bind_rows(SRdata,SRdata.pred,SRdata.release) %>%
    mutate(R=R/yscale,SSB=SSB/xscale)
  ymax <- max(alldata$R)
  year.max <- max(alldata$year,na.rm=T)
  tmp <- 1950:2030
  if(is.null(labeling.year)) labeling.year <- c(tmp[tmp%%5==0],year.max)
    alldata <- alldata %>% mutate(pick.year=ifelse(year%in%labeling.year,year,""))

  use_gamma <- SR_result$input$SR %in% c("Mesnil", "Shepherd", "BHS")

  if(!use_gamma){
    g1 <- ggplot(data=alldata,mapping=aes(x=SSB,y=R)) +
      stat_function(fun=SRF,args=list(a=SR_result$pars$a,
                                      b=SR_result$pars$b),color="deepskyblue3",lwd=1.3,
                    n=5000)
  }else{
    g1 <- ggplot(data=alldata,mapping=aes(x=SSB,y=R)) +
      stat_function(fun=SRF,args=list(a=SR_result$pars$a,
                                      b=SR_result$pars$b,
                                      gamma=SR_result$input$gamma),color="deepskyblue3",lwd=1.3,
                    n=5000)
  }

  if(!is.null(add_graph)) g1 <- g1+add_graph

  if(isTRUE(plot_CI)){
    if(!use_gamma){
      g1 <- g1+
        stat_function(fun=SRF_CI,
                      args=list(a=SR_result$pars$a,
                                b=SR_result$pars$b,
                                sigma=SR_result$pars$sd,
                                sign=-1,
                                CI=CI),
                      color="deepskyblue3",lty=3,n=5000)+
        stat_function(fun=SRF_CI,
                      args=list(a=SR_result$pars$a,
                                b=SR_result$pars$b,
                                sigma=SR_result$pars$sd,
                                sign=1,
                                CI=CI),
                      color="deepskyblue3",lty=3,n=5000)
    }else{
      g1 <- g1+
        stat_function(fun=SRF_CI,
                      args=list(a=SR_result$pars$a,
                                b=SR_result$pars$b,
                                gamma=SR_result$input$gamma,
                                sigma=SR_result$pars$sd,
                                sign=-1,
                                CI=CI),
                      color="deepskyblue3",lty=3,n=5000)+
        stat_function(fun=SRF_CI,
                      args=list(a=SR_result$pars$a,
                                b=SR_result$pars$b,
                                gamma=SR_result$input$gamma,
                                sigma=SR_result$pars$sd,
                                sign=1,
                                CI=CI),
                      color="deepskyblue3",lty=3,n=5000)
    }
  }

    if(!isTRUE(stringr::str_detect(version$os, pattern="darwin"))){
  g1 <- g1+  geom_path(data=dplyr::filter(alldata,type=="obs"),
                       aes(y=R,x=SSB),color="black") +
    geom_point(data=dplyr::filter(alldata,type=="obs"),
               aes(y=R,x=SSB,shape=weight),fill="white") +
    scale_shape_manual(values = shape_custom) +
      #    ggrepel::geom_text_repel(data=dplyr::filter(alldata,type=="obs"),
      #                             segment.alpha=0.5,nudge_y=5,
      #                             aes(y=R,x=SSB,label=pick.year)) +
    ggrepel::geom_text_repel(data=dplyr::filter(alldata,type=="obs"),
                             box.padding=box.padding,segment.color="gray",nudge_y=5,
                             aes(y=R,x=SSB,label=pick.year)) +
    theme_bw(base_size=14)+
    theme(legend.position = 'none') +
    theme(panel.grid = element_blank()) +
    xlab(str_c("親魚量 (",xlabel,")"))+
    ylab(str_c("加入量 (",ylabel,")"))+
    coord_cartesian(ylim=c(0,ymax*1.05),expand=0)
    }else{
      g1 <- g1+  geom_path(data=dplyr::filter(alldata,type=="obs"),
                           aes(y=R,x=SSB),color="black") +
        geom_point(data=dplyr::filter(alldata,type=="obs"),
                   aes(y=R,x=SSB,shape=weight),fill="white") +
        scale_shape_manual(values = shape_custom) +
        #    ggrepel::geom_text_repel(data=dplyr::filter(alldata,type=="obs"),
        #                             segment.alpha=0.5,nudge_y=5,
        #                             aes(y=R,x=SSB,label=pick.year)) +
      ggrepel::geom_text_repel(data=dplyr::filter(alldata,type=="obs"),
                               box.padding=box.padding,segment.color="gray",nudge_y=5,
                               aes(y=R,x=SSB,label=pick.year)) +
        theme_bw(base_size=14)+
        theme(legend.position = 'none') +
        theme(panel.grid = element_blank()) +
        theme(text = element_text(family = font_MAC)) +
        xlab(str_c("親魚量 (",xlabel,")"))+
        ylab(str_c("加入量 (",ylabel,")"))+
        coord_cartesian(ylim=c(0,ymax*1.05),expand=0)
    }

  if(is_release_data){
    g1 <- g1 + geom_point(data=dplyr::filter(alldata,type=="release"),
                          aes(y=R,x=SSB,shape=weight),color=gray(0.5))
  }

  if(recruit_intercept>0){
    if(!use_gamma){
      g1 <- g1+stat_function(fun=SRF,
                             args=list(a=SR_result$pars$a,
                                       b=SR_result$pars$b,
                                       recruit_intercept=recruit_intercept),
                             color="deepskyblue3",lwd=1.3,lty=2)
    }else{
      g1 <- g1+stat_function(fun=SRF,
                             args=list(a=SR_result$pars$a,
                                       b=SR_result$pars$b,
                                       gamma=SR_result$input$gamma,
                                       recruit_intercept=recruit_intercept),
                             color="deepskyblue3",lwd=1.3,lty=2)
    }

  }

  if(add.info){
    cap1 <- str_c("関数形: ",SR_result$input$SR,", 自己相関: ",SR_result$input$AR,
                  ", 最適化法",SR_result$input$method,", AICc: ",round(SR_result$AICc,2))
    if(sum(SRdata$weight=="0")>0) cap1 <- str_c(cap1, "\n パラメータ推定に利用（丸）,利用していない（バツ） ")
    if(is_release_data) cap1 <- str_c(cap1, "\n 灰色：放流＋天然、黒：天然のみ")

    if(!isTRUE(stringr::str_detect(version$os, pattern="darwin"))){
      g1 <- g1+labs(caption=cap1)
    }else{
      g1 <- g1+labs(caption=cap1,family = font_MAC)
    }
  }

  if(!is.null(refs)){
    g1 <- g1+geom_vline(xintercept=c(refs$Bmsy/xscale,refs$Blimit/xscale,refs$Bban/xscale),
                        linetype=2,
                        col=c(col.SBtarget,col.SBlimit,col.SBban))
  }
  g1
}

#'
#' @export

SRplot_gg <- plot.SR <- function(...){
    plot_SR(...)
}

#' 複数の再生産関係を比較する関数
#'
#' @param SRlist 再生産関係の推定結果のリスト。
#' @param biomass.unit 資源量の単位
#' @param number.unit 尾数の単位
#'
#' @examples
#' \dontrun{
#' data(res_sr_HSL1)
#' data(res_sr_HSL2)
#'
#' (g1 <- compare_SRfit(list(HSL1=res_sr_HSL1, HSL2=res_sr_HSL2),
#'                      biomass.unit=1000, number.unit=1000))
#'
#' }
#'
#' @export
#' @encoding UTF-8
#'
#'

compare_SRfit <- function(SRlist, biomass.unit=1000, number.unit=1000, newplot=TRUE, output_folder=""){

  font_MAC <- "HiraginoSans-W3"#"Japan1GothicBBB"#

  if(newplot){
    if(!is.null(SRlist[[1]]$input)){
      SRdata <- purrr::map_dfr(SRlist[], function(x){
        x$input$SRdata %>%
          as_tibble() %>%
          mutate(SSB=SSB, R=R)
      },.id="id")
    }
    else{ # for model average
      SRdata <- purrr::map_dfr(SRlist, function(x){
        x[[1]]$input$SRdata %>%
          as_tibble() %>%
          mutate(SSB=SSB, R=R)
      },.id="id")
    }

    if(is.null(SRlist)) names(SRlist) <- 1:length(SRlist)

        SRpred <- purrr::map_dfr(SRlist,
                             function(x) x$pred, .id="SR_type")
        #SRpred$再生産関係 <- as.factor(SRpred$再生産関係)
    font_MAC <- "HiraginoSans-W3"#"Japan1GothicBBB"#

    g1 <- ggplot(data=SRpred)
    g1 <- g1 + geom_line(data=SRpred,
                         mapping=aes(x=SSB/biomass.unit,y=R/number.unit, linetype=SR_type, col=SR_type))
    g1 <- g1 + geom_point(data=SRdata, mapping=aes(x=SSB/biomass.unit, y=R/number.unit), color="black")
    if(!isTRUE(stringr::str_detect(version$os, pattern="darwin"))){
    g1 <- g1 + xlim(c(0,max(SRdata$SSB/biomass.unit))) + ylim(c(0,max(SRdata$R/number.unit))) +
      labs(x = "親魚量（千トン）", y = "加入尾数（百万尾)") + theme_SH(legend.position="top")
    }else{
      g1 <- g1 + xlim(c(0,max(SRdata$SSB/biomass.unit))) + ylim(c(0,max(SRdata$R/number.unit))) +
        labs(x = "親魚量（千トン）", y = "加入尾数（百万尾)") + theme_SH(legend.position="top")+theme(text = element_text(family = font_MAC))
    }
    #g1
    ggsave_SH(g1, file=paste("./",output_folder,"/resSRcomp.png",sep=""))
    g1
  }
  else{
    if(!is.null(SRlist[[1]]$input)){
      SRdata <- purrr::map_dfr(SRlist[], function(x){
        x$input$SRdata %>%
          as_tibble() %>%
          mutate(SSB=SSB/biomass.unit, R=R/number.unit)
      },.id="id")
    }
    else{ # for model average
      SRdata <- purrr::map_dfr(SRlist, function(x){
        x[[1]]$input$SRdata %>%
          as_tibble() %>%
          mutate(SSB=SSB/biomass.unit, R=R/number.unit)
      },.id="id")
    }

    if(is.null(SRlist)) names(SRlist) <- 1:length(SRlist)

  g1 <- plot_SRdata(SRdata,type="gg")

  SRpred <- purrr::map_dfr(SRlist,
                           function(x) x$pred, .id="SR_type")
  g1 <- g1+geom_line(data=SRpred,mapping=aes(x=SSB/biomass.unit,y=R/number.unit,col=SR_type)) +
    theme(legend.position="top") +
    xlab(str_c("SSB (x",biomass.unit,")")) +
    ylab(str_c("Number (x",number.unit,")"))
  g1
  }
}

#' fit.SRregimeの結果で得られた再生産関係をプロットするための関数
#'
#' レジームごとの観察値と予測線が描かれる
#' @param resSRregime \code{fit.SRregime}の結果
#' @param xscale X軸のスケーリングの値（親魚量をこの値で除す）
#' @param xlabel X軸のラベル
#' @param yscale Y軸のスケーリングの値（加入量をこの値で除す）
#' @param ylabel Y軸のラベル
#' @param labeling.year ラベルに使われる年
#' @param show.legend 凡例を描くかどうか
#' @param legend.title 凡例のタイトル（初期設定は\code{"Regime"}）
#' @param regime.name 凡例に描かれる各レジームの名前（レジームの数だけ必要）
#' @param base_size \code{ggplot}のベースサイズ
#' @param add.info \code{AICc}や\code{regime.year}, \code{regime.key}などの情報を加えるかどうか
#' @param use.fit.SR パラメータの初期値を決めるのに\code{frasyr::fit.SR}を使う（時間の短縮になる）
#' @examples
#' \dontrun{
#' data(res_vpa)
#' SRdata <- get.SRdata(res_vpa, weight.year=1988:2015)
#' resSRregime <- fit.SRregime(SRdata, SR="HS", method="L2",
#'                             regime.year=c(1994,2003), regime.key=c(0,1,0),
#'                             regime.par = c("a","b","sd")[2:3])
#' g1 <- SRregime_plot(resSRregime, regime.name=c("Low","High"))
#' g1
#' }
#' @encoding UTF-8
#' @export
#'

SRregime_plot <- plot_SRregime <- function (SRregime_result,xscale=1000,xlabel="SSB",yscale=1,ylabel="R",
                           labeling.year = NULL, show.legend = TRUE, legend.title = "Regime",regime.name = NULL,
                           base_size = 16, add.info = TRUE, themeSH = FALSE) {
  pred_data = SRregime_result$pred %>% mutate(Category = "Pred")
  obs_data = select(SRregime_result$pred_to_obs, -Pred, -resid) %>% mutate(Category = "Obs")
  if(!is.null(SRregime_result$input$SRdata$weight)){
    obs_data$weight <- factor(SRregime_result$input$SRdata$weight,levels=c("0","1"))
  }else{
    obs_data$weight <- factor(1,levels=c("0","1"))
  }

  combined_data = full_join(pred_data, obs_data) %>%
    mutate(Year = as.double(Year))
  combined_data = rename(combined_data,Regime.num=Regime)
  if (is.null(regime.name)) {
    regime.name = c(unique(as.character(combined_data$Regime.num)))
  }
  Regime <- as.character(regime.name[combined_data$Regime.num])
  combined_data <- combined_data %>% mutate(Regime)

  if(prod(as.numeric(as.character(combined_data$weight)),na.rm=T)==0){
    scaleshapeval <- c(3,20)
  }else{
    scaleshapeval <- c(20)
  }

  Weight <- as.character(combined_data$weight)
  Weight[which(Weight=="0")] <- "unweighted"
  Weight[which(Weight=="1")] <- "weighted"
  #Weight[is.na(as.numeric(combined_data$weight))] <- "NA"
  combined_data <- combined_data %>% mutate(Weight)

  if (is.null(labeling.year)) labeling.year <- c(min(obs_data$Year),obs_data$Year[obs_data$Year %% 5 == 0],max(obs_data$Year))
  combined_data = combined_data %>%
    mutate(label=if_else(is.na(Year),as.numeric(NA),if_else(Year %in% labeling.year, Year, as.numeric(NA)))) %>%
    mutate(SSB = SSB/xscale, R = R/yscale)

  g1 = ggplot(combined_data, aes(x=SSB,y=R,label=label)) +
    geom_path(data=dplyr::filter(combined_data, Category=="Pred"),aes(group=Regime,colour=Regime,linetype=Regime),size=2, show.legend = show.legend)+
    geom_point(data=dplyr::filter(combined_data, Category=="Obs"),aes(group=Regime,colour=Regime,shape=Weight),size=3, show.legend = show.legend) +
    scale_shape_manual(values = scaleshapeval) +
    scale_color_manual(values = seq(length(unique(regime.name)))) +
    geom_path(data=dplyr::filter(combined_data, Category=="Obs"),colour="darkgray",size=1)+
    xlab(xlabel)+ylab(ylabel)+
    ggrepel::geom_label_repel()+
    theme_bw(base_size=base_size)+
    coord_cartesian(ylim=c(0,max(combined_data$R)*1.05),expand=0)
  # if (show.legend) {
  #   if (is.null(regime.name)) {
  #     regime.name = unique(combined_data$Regime)
  #   }
  #   g1 = g1 + #scale_colour_hue(name=legend.title, labels = regime.name) +
  #     scale_linetype_discrete(name=legend.title, labels = regime.name)
  # }
  if(themeSH) g1 = g1 + theme_SH()
  if (add.info) {
    if (is.null(SRregime_result$input$regime.year)) {
      g1 = g1 +
        # labs(caption=str_c(SRregime_result$input$SR,"(",SRregime_result$input$method,
        #                    "), regime_year: ", paste0(SRregime_result$input$regime.year,collapse="&"),
        #                    ", regime_key: ",paste0(SRregime_result$input$regime.key,collapse="->"),", AICc: ",round(SRregime_result$AICc,2)))
        labs(caption=str_c(SRregime_result$input$SR,"(",SRregime_result$input$method,
                           "), No Regime", ", AICc: ",round(SRregime_result$AICc,2))
        )

    }  else {
      g1 = g1 +
        # labs(caption=str_c(SRregime_result$input$SR,"(",SRregime_result$input$method,
        #                    "), regime_year: ", paste0(SRregime_result$input$regime.year,collapse="&"),
        #                    ", regime_key: ",paste0(SRregime_result$input$regime.key,collapse="->"),", AICc: ",round(SRregime_result$AICc,2)))
        labs(caption=str_c(SRregime_result$input$SR,"(",SRregime_result$input$method,
                           "), ", "regime_par: ", paste0(SRregime_result$input$regime.par,collapse="&"),", ",
                           paste0(SRregime_result$input$regime.year,collapse="&"),
                           ", ",paste0(SRregime_result$input$regime.key,collapse="->"),
                           ", AICc: ",round(SRregime_result$AICc,2))
        )

    }
  }
  g1
}

#'
#' @export

plot_SRregime <- function(...){
  SRregime_plot(...)
}

# 将来予測用 ----

#' 将来予測の複数の結果をggplotで重ね書きする
#'
#' @param vpares VPAの結果のオブジェクト(NULLでもOK)
#' @param future.list 将来予測の結果をリストで並べたもの(こちらは必須)
#' @param future.name 将来予測のリストの名前。ない場合はfuture.listについている名前を使う
#' @param CI_range 予測区間の範囲。デフォルトは８０\%でc(0.1,0.9)
#' @param maxyear 表示する年の最大
#' @param minyear 表示する年の最小
#' @param is.plot.CIrange 予測区間を表示するかどうか
#' @param is.plot.CIline 予測区間の上限と下限の線を引くか
#' @param what.plot Recruitment,SSB,biomass,catch,beta_gamma,U,Fratioのうち何をプロットするか。これらの文字列のベクトルで指定する
#' @param n_example 個々のシミュレーションの例を示す数
#' @param width_example 個々のシミュレーションをプロットする場合の線の太さ (default=0.7)
#' @param future.replicate どのreplicateを選ぶかを選択する。この場合n_exampleによる指定は無効になる
#' @param number.unit 尾数（加入尾数とか）のときの単位
#' @param biomass.unit 量の単位
#' @param number.name 尾数の凡例をどのように表示するか（「億尾」とか）
#' @param RP_name 管理基準値をどのように名前つけるか
#' @param Btarget Btargetの値
#' @param Blimit Blimitの値
#' @param Bban Bbanの値
#' @param MSY MSYの値
#' @param Umsy  Umsyの値
#' @param SPRtarget MSYのときのSPRの値
#' @param exclude.japanese.font 日本語を図中に表示しない
#' @param font.size フォントの大きさ
#' @param ncol 図を並べるときの列数
#' @param legend.position 凡例の位置。top, right, left, bottomなど
#' @param average_lwd 将来予測の平均値の線の太さ. 基本は1.
#' @param remove.last.vpa.year VPAの最終年のデータのプロットを除くかどうか（last.catch.zero用オプション）
#'
#' @encoding UTF-8
#' @export

plot_futures <- function(vpares=NULL,
                         future.list=NULL,
                         future.name=names(future.list),
                         future_tibble=NULL,
                         CI_range=c(0.1,0.9),
                         maxyear=NULL,
                         minyear=NULL,
                         is.plot.CIrange=TRUE,
                         is.plot.CIline = TRUE,
                         what.plot=c("Recruitment","SSB","biomass","cbiomass","catch","beta_gamma","U","Fratio"),
                         biomass.unit=1,
                         number.unit=1,
                         number.name="",
                         RP_name=c("Btarget","Blimit","Bban"),
                         Btarget=0,Blimit=0,Bban=0,#Blow=0,
                         MSY=0,Umsy=0,
                         SPRtarget=NULL,
                         exclude.japanese.font=FALSE, # english version
                         average_lwd=1,
                         n_example=3, # number of examples
                         example_width=0.7, # line width of examples
                         future.replicate=NULL,
                         seed=1, # seed for selecting the above example
                         legend.position="top",
                         type="detail",
                         font.size=16,
                         ncol=3,
                         remove.last.vpa.year = FALSE
){

  col.SBtarget <- "#00533E"
  col.SBlim <- "#edb918"
  col.SBban <- "#C73C2E"
  col.MSY <- "black"
  col.Ftarget <- "#714C99"
  col.betaFtarget <- "#505596"

  for(i in 1:length(future.list)){
    if(class(future.list[[i]])=="future_new")
      future.list[[i]] <- format_to_old_future(future.list[[i]])
    det.run <- FALSE
  }

  if(!isTRUE(exclude.japanese.font)){
    junit <- c("","十","百","千","万")[log10(biomass.unit)+1]

    if(type=="detail"){
        rename_list <- tibble(stat=c("Recruitment","SSB","biomass","cbiomass","catch","beta_gamma","U","Fratio"),
                              jstat=c(str_c("加入尾数(",number.name,")"),
                                      str_c("親魚量 (",junit,"トン)"),
                                      str_c("資源量 (",junit,"トン)"),
                                      str_c("漁獲資源量 (",junit,"トン)"),
                                      str_c("漁獲量 (",junit,"トン)"),
                                      "beta_gamma(F/Fmsy)",
                                      "漁獲割合(%)",
                                      "漁獲圧の比(F/Fmsy)"))
    }
    if(type=="simple"){
        rename_list <- tibble(stat=c("Recruitment","SSB","biomass","cbiomass","catch","beta_gamma","U","Fratio"),
                              jstat=c(str_c("加入尾数(",number.name,")"),
                                      str_c("将来の親魚量 (",junit,"トン)"),
                                      str_c("資源量 (",junit,"トン)"),
                                      str_c("漁獲資源量 (",junit,"トン)"),
                                      str_c("将来の漁獲量 (",junit,"トン)"),
                                      "beta_gamma(F/Fmsy)",
                                      "漁獲割合(%)",
                                      "漁獲圧の比(F/Fmsy)"))
    }
  }
  else{
    junit <- c("","10","100","1000","10,000")[log10(biomass.unit)+1]
    #    require(tidyverse,quietly=TRUE)

    rename_list <- tibble(stat=c("Recruitment","SSB","biomass","cbiomass","catch","beta_gamma","U","Fratio"),
                          jstat=c(str_c("Recruits(",number.name,"fish)"),
                                  str_c("SB (",junit,"MT)"),
                                  str_c("Biomass (",junit,"MT)"),
                                  str_c("cBiomass (",junit,"MT)"),
                                  str_c("Catch (",junit,"MT)"),
                                  "multiplier to Fmsy",
                                  "Catch/Biomass (U)",
                                  "F ratio (F/Fmsy)"))
  }

  # define unit of value
  rename_list <- rename_list %>%
    mutate(unit=dplyr::case_when(stat%in%c("SSB","biomass","cbiomass","catch") ~ biomass.unit,
                                 stat%in%c("Recruitment")           ~ number.unit,
                                 stat%in%c("U")                     ~ 0.01,
                                 TRUE                               ~ 1))

  rename_list <- rename_list %>% dplyr::filter(stat%in%what.plot)

  if(!is.null(future.list)){
    if(is.null(future.name)) future.name <- str_c("s",1:length(future.list))
    names(future.list) <- future.name
  }
  else{
    if(is.null(future.name)) future.name <- str_c("s",1:length(unique(future_tibble$HCR_name)))
  }

  if(is.null(future_tibble)) future_tibble <- purrr::map_dfr(future.list,convert_future_table,.id="scenario")

  future_tibble <-
    future_tibble %>%
    dplyr::filter(stat%in%rename_list$stat) %>%
    mutate(stat=factor(stat,levels=rename_list$stat)) %>%
    left_join(rename_list) %>%
    mutate(value=value/unit)

  if(is.null(future.replicate)){
    set.seed(seed)
    future.replicate <- sample(2:max(future_tibble$sim),n_example)
  }
  future.example <- future_tibble %>%
    dplyr::filter(sim%in%future.replicate) %>%
    group_by(sim,scenario)

  if(is.null(maxyear)) maxyear <- max(future_tibble$year)
  if(is.null(minyear)) minyear <- min(future_tibble$year)

  #  min.age <- as.numeric(rownames(vpares$naa)[1])
  if(!is.null(vpares)){
    vpa_tb <- convert_vpa_tibble(vpares,SPRtarget=SPRtarget) %>%
      mutate(scenario=type,year=as.numeric(year),
             stat=factor(stat,levels=rename_list$stat),
             mean=value,sim=0)%>%
      dplyr::filter(stat%in%rename_list$stat) %>%
      left_join(rename_list) %>%
      mutate(value=value/unit,mean=mean/unit)

    if ( isTRUE(remove.last.vpa.year) ) vpa_tb = vpa_tb %>% filter(year<max(year))

    # 将来と過去をつなげるためのダミーデータ
    tmp <- vpa_tb %>% group_by(stat) %>%
      summarise(value=tail(value[!is.na(value)],n=1,na.rm=T),
                year=tail(year[!is.na(value)],n=1,na.rm=T),sim=0)
    future.dummy <- purrr::map_dfr(future.name,function(x) mutate(tmp,scenario=x))
  }
  else{
    future.dummy <- NULL
    vpa_tb <- NULL
  }

  org.warn <- options()$warn
  options(warn=-1)
  future_tibble <-
    bind_rows(future_tibble,vpa_tb,future.dummy) %>%
    mutate(stat=factor(stat,levels=rename_list$stat)) %>%
    mutate(scenario=factor(scenario,levels=c(future.name,"VPA"))) #%>%
  #        mutate(value=ifelse(stat%in%c("beta_gamma","U"),value,value/biomass.unit))

  future_tibble.qt <-
    future_tibble %>% group_by(scenario,year,stat) %>%
    summarise(low=quantile(value,CI_range[1],na.rm=T),
              high=quantile(value,CI_range[2],na.rm=T),
              median=median(value,na.rm=T),
              mean=mean(value))

  # make dummy for y range
  dummy <- future_tibble %>% group_by(stat) %>% summarise(max=max(value)) %>%
    mutate(value=0,year=min(future_tibble$year,na.rm=T)) %>%
    select(-max)

  dummy2 <- future_tibble %>% group_by(stat) %>%
    summarise(max=max(quantile(value,CI_range[2],na.rm=T))) %>%
    mutate(value=max*1.1,
           year=min(future_tibble$year,na.rm=T)) %>%
    select(-max)

  future_tibble.qt <- left_join(future_tibble.qt,rename_list) %>%
    mutate(jstat=factor(jstat,levels=rename_list$jstat))

  dummy     <- left_join(dummy,rename_list,by="stat") %>% dplyr::filter(!is.na(stat))
  dummy2    <- left_join(dummy2,rename_list,by="stat") %>% dplyr::filter(!is.na(stat))

  if("SSB" %in% what.plot){
    ssb_RP <- tibble(jstat = dplyr::filter(rename_list, stat == "SSB") %>%
                       dplyr::pull(jstat),
                     value = c(Btarget, Blimit, Bban) / biomass.unit,
                     RP_name = RP_name)
  }
  if("catch" %in% what.plot){
    catch_RP <- tibble(jstat=dplyr::filter(rename_list, stat == "catch") %>%
                         dplyr::pull(jstat),
                       value=MSY/biomass.unit,
                       RP_name="MSY")
  }
  if("U" %in% what.plot){
    U_RP <- tibble(jstat=dplyr::filter(rename_list, stat == "U") %>%
                     dplyr::pull(jstat),
                   value=Umsy,
                   RP_name="U_MSY")
  }

  options(warn=org.warn)

  if (average_lwd == example_width) {
      warning("average_lwd と example_width が同じ太さです. 図の脚注との整合を確認してください.")
  }

  g1 <- future_tibble.qt %>% dplyr::filter(!is.na(stat)) %>%
    ggplot()

  if(isTRUE(is.plot.CIrange)){
    if (isTRUE(is.plot.CIline)) {
      g1 <- g1+
        geom_line(data=dplyr::filter(future_tibble.qt,!is.na(stat) & scenario!="VPA" & year %in% minyear:maxyear),
                  mapping=aes(x=year,y=high,lty=scenario,color=scenario))+
        geom_line(data=dplyr::filter(future_tibble.qt,!is.na(stat) & scenario!="VPA" & year %in% minyear:maxyear),
                  mapping=aes(x=year,y=low,lty=scenario,color=scenario))
    }
    g1 <- g1 +
      geom_ribbon(data=dplyr::filter(future_tibble.qt,!is.na(stat) & scenario!="VPA" & year %in% minyear:maxyear),
                  mapping=aes(x=year,ymin=low,ymax=high,fill=scenario),alpha=0.4)+
      geom_line(data=dplyr::filter(future_tibble.qt,!is.na(stat) & scenario!="VPA" & year %in% minyear:maxyear),
                mapping=aes(x=year,y=mean,color=scenario),lwd=average_lwd)
  }
  #    else{
  #        g1 <- g1+
  #            geom_line(data=dplyr::filter(future_tibble.qt,!is.na(stat) & scenario=="VPA"),
  #                      mapping=aes(x=year,y=mean,color=scenario),lwd=1)#+
  #    }

  g1 <- g1+
    geom_blank(data=dummy,mapping=aes(y=value,x=year))+
    geom_blank(data=dummy2,mapping=aes(y=value,x=year))+
    #theme_bw(base_size=font.size) +
    #        coord_cartesian(expand=0)+
    scale_y_continuous(expand=expand_scale(mult=c(0,0.05)),labels = scales::comma)+
    facet_wrap(~factor(jstat,levels=rename_list$jstat),scales="free_y",ncol=ncol)+
    xlab("年")+ylab("")+ labs(fill = "",linetype="",color="")+
    xlim(minyear,maxyear)

  if("SSB" %in% what.plot && Btarget*Blimit*Bban>0){
    g1 <- g1 + geom_hline(data = ssb_RP,
                          aes(yintercept = value,linetype=RP_name),
						  color = c(col.SBtarget, col.SBlim, col.SBban))+
						  scale_linetype_manual(name="",values=c("solid","dashed",unlist(format_type()[1,3])[[1]],unlist(format_type()[1,3])[[1]],unlist(format_type()[3,3])[[1]],unlist(format_type()[2,3])[[1]],unlist(format_type()[1,3])[[1]]))
  }

  if("catch" %in% what.plot && MSY!=0 ){
    g1 <- g1 + geom_hline(data = catch_RP,
                          aes(yintercept = value, linetype = RP_name),
                          color = c(col.MSY))
  }

  if("U" %in% what.plot && Umsy!=0){
    g1 <- g1 + geom_hline(data = U_RP,
                          aes(yintercept = value, linetype = RP_name),
                          color = c(col.MSY))
  }

  if(n_example>0){
      if(n_example>1){
        tmpdata <- dplyr::filter(future.example,year <= maxyear)
        g1 <- g1 + geom_line(data=tmpdata,
                           mapping=aes(x=year,y=value,
                                       alpha=as.factor(sim),
                                       color=scenario),
                           lwd=example_width)
    }
    else{
      g1 <- g1 + geom_line(data=dplyr::filter(future.example,year %in% minyear:maxyear),
                           mapping=aes(x=year,y=value,
                                       color=scenario),
                           lwd=example_width)
    }
    g1 <- g1+scale_alpha_discrete(guide="none")
  }

  if(!exclude.japanese.font){
    caption_string <- str_c("(塗り:", CI_range[1]*100,"-",CI_range[2]*100,
                          "%予測区間, ",
                         ifelse(average_lwd < example_width, "細い", "太い"),
                         "実線: 平均値",
                         dplyr::case_when(n_example>0 & example_width < average_lwd ~ ", 細い実線: シミュレーションの1例)",
                                          n_example>0 & example_width > average_lwd ~ ", 太い実線: シミュレーションの1例)",
                                          n_example>0 & example_width == average_lwd ~ ", 薄い実線: シミュレーションの1例)",
                                          TRUE ~ ")"))
  }
  else{
    caption_string <- str_c("(Fill:", CI_range[1]*100,"-",CI_range[2]*100,
                          "%prediction interval, ",
                         ifelse(average_lwd < example_width, "thin", "thick"),
                         "sold line: average",
                         dplyr::case_when(n_example>0 & example_width < average_lwd ~ ", thin solid line: an example in the simulation)",
                                          n_example>0 & example_width > average_lwd ~ ", bold solid line: an example in the simulation)",
                                          n_example>0 & example_width == average_lwd ~ ", pale solid line: an example in the simulation)",
                                          TRUE ~ ")"))
  }

  g1 <- g1 + guides(lty=guide_legend(ncol=3),
                    fill=guide_legend(ncol=3),
                    col=guide_legend(ncol=3))+
    theme_SH(base_size=font.size,legend.position=legend.position)+
    scale_color_hue(l=40)+
    labs(caption = caption_string)

  g1 <- g1 +
    geom_line(data=dplyr::filter(future_tibble.qt,!is.na(stat) & scenario=="VPA"),
              mapping=aes(x=year,y=mean),lwd=1,color="black")# VPAのプロット
  return(g1)
}


#' 複数の将来予測の結果をプロットする（ggplotは使わず）
#'
#' @param fres.list future.vpaからの出力結果をリストで並べたもの
#' @encoding UTF-8
#'
#'
#'
#' @export

plot_futures_simple <- function(fres.list,conf=c(0.1,0.5,0.9),target="SSB",legend.text="",xlim.tmp=NULL,y.scale=1,det.run=TRUE){

  if(legend.text=="") legend.text <- names(fres.list)
  if(is.null(legend.text)) legend.text <- 1:length(fres.list)

  for(i in 1:length(fres.list)){
    if(class(fres.list[[i]])=="future_new")
      fres.list[[i]] <- format_to_old_future(fres.list[[i]])
    det.run <- FALSE
  }

  if(isTRUE(det.run)) select_col <- -1  else select_col <- TRUE

  if(target=="SSB")  aa <- lapply(fres.list,function(x) apply(x$vssb[,select_col],1,quantile,probs=conf))
  if(target=="Biomass") aa <- lapply(fres.list,function(x) apply(x$vbiom[,select_col],1,quantile,probs=conf))
  if(target=="Catch") aa <- lapply(fres.list,function(x) apply(x$vwcaa[,select_col],1,quantile,probs=conf))
  if(target=="Recruit"){
    if(is.null(x$recruit)) x$recruit <- x$naa
    aa <- lapply(fres.list,function(x) apply(x$recruit[,select_col],1,quantile,probs=conf))
  }

  if(is.null(xlim.tmp)) xlim.tmp <- as.numeric(range(unlist(sapply(aa,function(x) colnames(x)))))
  plot(0,max(unlist(aa)),type="n",xlim=xlim.tmp,
       ylim=y.scale*c(0,max(unlist(aa))),xlab="Year",ylab=target)
  lapply(1:length(aa),function(i) matpoints(colnames(aa[[i]]),t(aa[[i]]),col=i,type="l",lty=c(2,1,2)))
  legend("bottomright",col=1:length(fres.list),legend=legend.text,lty=1)
  invisible(aa)
}

#' @export
#'

plot.future <- function(...){
    plot_future_simple(...)
}

#' @export
#'

plot.futures <- function(...){
    plot_futures_simple(...)
}

#' 一つの将来予測の結果をプロットする（ggplotは使わず）
#'
#' @param fres0 future.vpaからの出力結果
#' @encoding UTF-8
#'
#'
#'
#' @export

plot_future_simple <- function(fres0,ylim.tmp=NULL,xlim.tmp=NULL,vpares=NULL,what=c(TRUE,TRUE,TRUE),conf=0.1,N.line=0,det.run=TRUE,
                        label=c("Biomass","SSB","Catch"),is.legend=TRUE,add=FALSE,col=NULL,...){
  ## 暗黙に、vssbなどのmatrixの1列目は決定論的なランの結果と仮定している

  if(class(fres0)=="future_new"){
    fres0 <- format_to_old_future(fres0)
    det.run <- FALSE
  }

  if(is.null(col)) col <- 1
  matplot2 <- function(x,add=FALSE,...){
    if(add==FALSE) matplot(rownames(x),x,type="l",lty=c(2,1,2),col=col,xlab="Year",...)
    if(add==TRUE) matpoints(rownames(x),x,type="l",lty=c(2,1,2),col=col,xlab="Year",...)
  }

  if(is.null(xlim.tmp)) xlim.tmp <- range(as.numeric(rownames(fres0$vssb)))

  if(isTRUE(det.run)) select_col <- -1  else select_col <- TRUE

  if(what[1]){
    matplot2(x <- t(apply(fres0$vbiom[,select_col],1,
                          quantile,probs=c(conf,0.5,1-conf))),
             ylim=c(0,ifelse(is.null(ylim.tmp),max(x),ylim.tmp[1])),
             xlim=xlim.tmp,
             ylab=label[1],main=label[1],add=add,...)
    points(rownames(fres0$vbiom),apply(fres0$vbiom[,select_col],1,mean),type="b",pch=1)
    if(isTRUE(det.run)) points(rownames(fres0$vbiom),as.numeric(fres0$vbiom[,1]),type="b",pch=3)
    if(!is.null(vpares)){
      points(colnames(vpares$baa),colSums(vpares$baa),type="o",pch=20)
    }
    if(N.line>0) matpoints(rownames(fres0$vbiom),fres0$vbiom[,2:(N.line+1)],col="gray",type="l",lty=1)
  }

  if(what[2]){
    matplot2(x <- t(apply(fres0$vssb[,select_col],1,quantile,
                          probs=c(conf,0.5,1-conf))),
             ylim=c(0,ifelse(is.null(ylim.tmp),max(x),ylim.tmp[2])),
             xlim=xlim.tmp,
             ylab=label[2],main=label[2],add=add,...)
    points(rownames(fres0$vssb),apply(fres0$vssb[,select_col],1,mean),type="b",pch=1)
    if(isTRUE(det.run)) points(rownames(fres0$vssb),as.numeric(fres0$vssb[,1]),type="b",pch=3)
    if(!is.null(fres0$input$Frec))
      if(!is.null(fres0$input$Frec$scenario))
        if(fres0$input$Frec$scenario!="catch.mean"){
          abline(h=fres0$input$Frec$Blimit,col=2)
          abline(v=fres0$input$Frec$future.year,col=2)
        }
    if(!is.null(vpares)){
      points(colnames(vpares$ssb),colSums(vpares$ssb),type="o",pch=20)
    }
    if(N.line>0) matpoints(rownames(fres0$vssb),fres0$vssb[,2:(N.line+1)],col="gray",type="l",lty=1)
  }

  if(what[3]){
    matplot2(x <- t(apply(fres0$vwcaa[,select_col],1,
                          quantile,probs=c(conf,0.5,1-conf))),
             ylim=c(0,ifelse(is.null(ylim.tmp),max(x),ylim.tmp[3])),
             xlim=xlim.tmp,
             ylab=label[3],main=label[3],add=add,...)
    points(rownames(fres0$vwcaa),apply(fres0$vwcaa[,select_col],1,mean),type="b",pch=1)
    if(isTRUE(det.run)) points(rownames(fres0$vwcaa),as.numeric(fres0$vwcaa[,1]),type="b",pch=3)
    if(!is.null(fres0$input$Frec))
      if(fres0$input$Frec$scenario=="catch.mean"){
        abline(h=fres0$input$Frec$Blimit,col=2)
        abline(v=fres0$input$Frec$future.year,col=2)
      }
    if(!is.null(vpares)){
      points(colnames(vpares$baa),colSums(vpares$input$dat$caa*vpares$input$dat$waa),type="o",pch=20)
    }
    if(N.line>0) matpoints(rownames(fres0$vwcaa),fres0$vwcaa[,2:(N.line+1)],col="gray",type="l",lty=1)
  }
  if(is.legend){
    if(sum(what)>1) plot(1:10,type = "n",ylab = "", xlab = "", axes = F)
    legend("topleft",lty=c(NA,NA,1,2),legend=c("Deterministic","Mean","Median",paste(100-(conf*2)*100,"%conf")),pch=c(3,1,NA,NA))
  }

}



#'
#' @export
#'

plot_waa <- function(vres){
  lm.list <- list()
  nage <- nrow(vres$naa)
  col.tmp <- rainbow(nage)
  logx <- log(unlist(vres$naa))
  logy <- log(unlist(vres$input$dat$waa))
  ages <- as.numeric(rep(rownames(vres$naa),ncol(vres$naa)))
  u.age <- unique(ages)
  plot(logx,logy,col=col.tmp[1+ages],xlab="log(N)",ylab="log(weight)")
  for(i in 1:length(u.age)){
    tmp <- ages==u.age[i] & logy>-Inf & logx>-Inf
    if(sum(tmp,na.rm=TRUE)>0){
      lm.list[[i]] <- lm(logy[tmp]~logx[tmp])
      l.type <- ifelse(summary(lm.list[[i]])$coeff[2,4]<0.05,1,2)
      if(!is.na(l.type)) abline(lm.list[[i]],col=col.tmp[1+ages[i]],lty=l.type)
    }
  }
  title(vres$stockid,line=0.2)
  legend("bottomleft",lty=c(1:2,rep(1,nage)),
         col=c(1,1,col.tmp),
         legend=c("p<0.05","p>0.05",paste("Age",u.age)))
  return(lm.list)
}

# 管理基準値・ステークホルダー向け ----

#' 漁獲量曲線(yield curve)を書く
#'
#' @param MSY_obj est.MSYの結果のオブジェクト
#' @param refs_base est.MSYから得られる管理基準値の表
#' @param future 将来予測結果のリスト。与えられると将来予測の結果を重ね書きする
#' @param future.replicat 将来予測結果から特定のreplicateのみを示す。futureで与えたリストの長さのベクトルを与える。
#' @param past  VPA結果。与えられると過去の推定値を重ね書きする
#' @encoding UTF-8
#'
#' @export

plot_yield <- function(MSY_obj,refs_base,
                       refs.label=NULL, # label for reference point
                       refs.color=c("#00533E","#edb918","#C73C2E"),
                       AR_select=FALSE,xlim.scale=1.1,
                       biomass.unit=1,labeling=TRUE,lining=TRUE,
                       age.label.ratio=0.6, # 年齢のラベルを入れる位置（xの最大値からの割合)
                       #                       family = "JP1",
                       ylim.scale=1.2,future=NULL,
                       future.replicate=NULL,
                       past=NULL,
                       past_year_range=NULL,
                       plus_group=TRUE,
                       future.name=NULL){

  junit <- c("","十","百","千","万")[log10(biomass.unit)+1]

  if ("trace" %in% names(MSY_obj)) {
    trace.msy <- MSY_obj$trace
  } else {
    trace.msy <- MSY_obj
  }

  #    require(tidyverse,quietly=TRUE)
  #    require(ggrepel)

  trace <- get.trace(trace.msy) %>%
    mutate("年齢"=age,ssb.mean=ssb.mean/biomass.unit,value=value/biomass.unit) %>%
    dplyr::filter(!is.na(value))

  refs_base <- refs_base %>%
    mutate(RP.definition=ifelse(is.na(RP.definition),"",RP.definition)) %>%
    mutate(SSB=SSB/biomass.unit)
  if("AR"%in%names(refs_base)) refs_base <- refs_base %>% dplyr::filter(AR==AR_select)

  ymax <- trace %>%
    group_by(ssb.mean) %>%
    summarise(catch.mean=sum(value))
  ymax <- max(ymax$catch.mean)


  g1 <- trace %>%
    ggplot2::ggplot()
 
  if(is.null(future.name)) future.name <- 1:length(future)

  if(is.null(refs.label)) {
    refs.label <- str_c(refs_base$RP_name,":",refs_base$RP.definition)
    refs.color <- 1:length(refs.label)
  }
    refs_base$refs.label <- refs.label

  plus.char <- ifelse(plus_group==TRUE,"+","")

  xmax <- max(trace$ssb.mean,na.rm=T)
  age.label.position <- trace$ssb.mean[which.min((trace$ssb.mean-xmax*xlim.scale*age.label.ratio)^2)]
  age.label <- trace %>% dplyr::filter(round(age.label.position,1)==round(ssb.mean,1))%>%
    mutate(cumcatch=cumsum(value)-value/2)%>%
    mutate(age=as.numeric(as.character(age)))
  age.label <- age.label %>%
    mutate(age_name=str_c(age,ifelse(age.label$age==max(age.label$age),plus.char,""),"歳"))
  
  legend.labels <- as.vector(age.label$age_name)
  
  nb.cols <- length(unique(trace$age)) # 年齢グループが多い場合に対応できるように変更
   mycolors <- grDevices::colorRampPalette(c("#F7FBFF", "#DEEBF7", "#C6DBEF", "#9ECAE1", "#6BAED6", "#4292C6",
"#2171B5", "#084594"))(nb.cols)
  
  g1 <- g1 + geom_area(aes(x=ssb.mean,y=value,fill=`年齢`),col="black",alpha=0.5,lwd=1*0.3528,stat="identity") +
    #    geom_line(aes(x=ssb.mean,y=catch.CV,fill=age)) +
    #    scale_y_continuous(sec.axis = sec_axis(~.*5, name = "CV catch"))+
    scale_fill_manual(values=mycolors,labels=rev(legend.labels)) +
    theme_bw() +
    #theme(legend.position = 'none') +
    #    geom_point(data=refs_base,aes(y=Catch,x=SSB,shape=refs.label,color=refs.label),size=4)+
    #形は塗りつぶしができる形にすること
    scale_shape_manual(values = c(21, 24,5,10)) +
    #塗りつぶし色を指定する
    scale_color_manual(values = c("darkgreen", "darkblue","orange","yellow"))+
    theme(panel.grid = element_blank(),axis.text=element_text(color="black")) +
    coord_cartesian(xlim=c(0,xmax*xlim.scale),
                    ylim=c(0,ymax*ylim.scale),expand=0) +
    #geom_text(data=age.label,
     #         mapping = aes(y = cumcatch, x = ssb.mean, label = age_name)#,
              #                            family = family
    #) +
    #    geom_text_repel(data=refs_base,
    #                     aes(y=Catch,x=SSB,label=refs.label),
    #                     size=4,box.padding=0.5,segment.color="gray",
    #                     hjust=0,nudge_y      = ymax*ylim.scale-refs_base$Catch,
    #    direction    = "y",
    #    angle        = 0,
    #    vjust        = 0,
    #        segment.size = 1)+
    xlab(str_c("平均親魚量 (",junit,"トン)")) + ylab(str_c("平均漁獲量 (",junit,"トン)"))

  if(!is.null(future)){
    futuredata <- NULL
    for(j in 1:length(future)){
      if(class(future[[j]])=="future_new"){
        future_init <- future[[j]]$input$tmb_data$future_initial_year
        future_init <- as.numeric(dimnames(future[[j]]$naa)[[2]][future_init])
        future[[j]] <- format_to_old_future(future[[j]])
      }
      else{
        future_init <- 0
      }
      if(is.null(future.replicate)){
        futuredata <- bind_rows(futuredata,
                                tibble(
                                  year        =as.numeric(rownames(future[[j]]$vssb)),
                                  ssb.future  =apply(future[[j]]$vssb[,-1],1,mean)/biomass.unit,
                                  catch.future=apply(future[[j]]$vwcaa[,-1],1,mean)/biomass.unit,
                                  scenario=future.name[j]))
      }
      else{
        futuredata <- bind_rows(futuredata,
                                tibble(
                                  year        =as.numeric(rownames(future[[j]]$vssb)),
                                  ssb.future  =future[[j]]$vssb[,future.replicate[j]]/biomass.unit,
                                  catch.future=future[[j]]$vwcaa[,future.replicate[j]]/biomass.unit,
                                  scenario=future.name[j]))
      }
      futuredata <- futuredata %>% group_by(scenario) %>%
        dplyr::filter(year > future_init)
      g1 <- g1 +
        geom_path(data=futuredata,
                  mapping=aes(x       =ssb.future,
                              y       = catch.future,
                              linetype=factor(scenario),
                              color   = factor(scenario)),
                  lwd=1)+
        geom_point(data=futuredata,
                   mapping=aes(x    =ssb.future,
                               y    =catch.future,
                               shape=factor(scenario),
                               color=factor(scenario)),
                   size=3)
    }
  }

  if(!is.null(past)){
    catch.past = unlist(colSums(past$input$dat$caa*past$input$dat$waa, na.rm=T)/biomass.unit)
    if (past$input$last.catch.zero && !is.null(future)) {
      catch.past[length(catch.past)] = apply(future[[1]]$vwcaa[,-1],1,mean)[1]/biomass.unit
    }
    pastdata <- tibble(
      year      =as.numeric(colnames(past$ssb)),
      ssb.past  =unlist(colSums(past$ssb, na.rm=T))/biomass.unit,
      catch.past=catch.past
    )

    if(past_year_range[1] > 0 && !is.null(past_year_range))
      pastdata <- pastdata %>% dplyr::filter(year%in%past_year_range)

    g1 <- g1 +
      geom_path(data=pastdata,
                mapping=aes(x=ssb.past,y=catch.past),
                color="darkred",lwd=1,alpha=0.9)
  }

  if(isTRUE(lining)){
    #        ylim.scale.factor <- rep(c(0.94,0.97),ceiling(length(refs.label)/2))[1:length(refs.label)]
    g1 <- g1 + geom_vline(xintercept=refs_base$SSB,lty=c(unlist(format_type()[1,3]),unlist(format_type()[2,3]),unlist(format_type()[3,3])),lwd=0.6,color=refs.color)+
      ggrepel::geom_label_repel(data=refs_base,
                                aes(y=ymax*ylim.scale*0.85,
                                    x=SSB,label=refs.label),
                                direction="x",size=11*0.282,nudge_y=ymax*ylim.scale*0.9)
  }

  if(isTRUE(labeling)){
    g1 <- g1 +
      geom_point(data=refs_base,
                 aes(y=Catch,x=SSB))+
      ggrepel::geom_label_repel(data=refs_base,
                                aes(y=Catch,x=SSB,label=refs.label),
                                #                            size=4,box.padding=0.5,segment.color="black",
                                hjust=0,#nudge_y      = ymax*ylim.scale-refs_base$Catch/2,
                                direction="y",angle=0,vjust        = 0,segment.size = 1)
    #             geom_label_repel(data=tibble(x=c(1,limit.ratio,ban.ratio),
    #                                          y=max.U,
    #                                          label=c("目標管理基準値","限界管理基準値","禁漁水準")),
    #                              aes(x=x,y=y,label=label),
    #                              direction="y",angle=0,nudge_y=max.U
  }


  return(g1)

}

#' Kobe plotを書く
#'
#' @param FBdata tibble(year=, Bratio=, Fratio=, Uratio=, DBratio= )の形式のデータ。UratioはU/Umsy、DBratioはB/(dynamic Bmsy)
#' @param xcol_name x軸に利用するデータの列名
#' @param ycol_name y軸に利用するデータの列名
#' @param refs_base est_MSYRPから得られる管理基準値の表
#' @param vpares 上記のdataを与えずにVPAの結果のオブジェクトなどからBratioを計算する場合のVPAの結果
#' @param Fratio 上記のdataを与えずにVPAの結果のオブジェクトなどからFratioを計算する場合、かつylab.type="F"の場合のFratioの値
#' @param Btarget est_MSYRPから得られる管理基準値の表の中のRP.definitionの列でtargetとする行のラベル
#' @param Blimit est_MSYRPから得られる管理基準値の表の中のRP.definitionの列でlimitとする行のラベル
#' @param Bban est_MSYRPから得られる管理基準値の表の中のRP.definitionの列でbanとする行のラベル
#'
#' @encoding UTF-8
#'
#' @export
#'

plot_kobe_gg <- plot_kobe <- function(FBdata=NULL,
                                      vpares=NULL,
                                      refs_base=NULL,
                                      roll_mean=1,
                                      ylab_name="Fratio",
                                      xlab_name="Bratio",
                                      Btarget=c("Btarget0"),
                                      Blimit=c("Blimit0"),
                                      Bban=c("Bban0"),
                                      write.vline=TRUE,
#                                      ylab.type="F", # or "U"
                                      labeling.year=NULL,
                                      RP.label=c("目標管理基準値","限界管理基準値","禁漁水準"),
                                      refs.color=c("#00533E","#edb918","#C73C2E"),
                                      Fratio=NULL,
                                      yscale=1.2,xscale=1.2,
                                      HCR.label.position=c(1,1),
                                      beta=NULL,
                                      plot.year="all"){

  if(!is.null(FBdata)) assertthat::assert_that(all(c(ylab_name, xlab_name) %in% colnames(FBdata)), TRUE)

  target.RP <- derive_RP_value(refs_base,Btarget)
  limit.RP <- derive_RP_value(refs_base,Blimit)
#  low.RP <- derive_RP_value(refs_base,Blow)
  ban.RP <- derive_RP_value(refs_base,Bban)

#  low.ratio <- low.RP$SSB/target.RP$SSB
  limit.ratio <- limit.RP$SSB/target.RP$SSB
  ban.ratio <- ban.RP$SSB/target.RP$SSB

  if(is.null(FBdata)){
    vpa_tb <- convert_vpa_tibble(vpares)
    FBdata <- vpa_tb %>% dplyr::filter(stat=="U" | stat=="SSB") %>%
      spread(key=stat,value=value) %>%
      mutate(Uratio=RcppRoll::roll_mean(U/target.RP$U,n=roll_mean,fill=NA,align="right"),
             Bratio=RcppRoll::roll_mean(SSB/target.RP$SSB,n=roll_mean,fill=NA,align="right")) %>%
      arrange(year)
    if(!is.null(Fratio)) FBdata <- FBdata %>% mutate(Fratio=Fratio)
  }

  FBdata <- FBdata %>% mutate(year_group=1)
  if(plot.year[1]!="all") {
    diff.year <- plot.year[which(diff(plot.year)>1)+1]
    FBdata <- FBdata %>% filter(year %in% plot.year)

    if (length(diff.year)>0) {
      for(i in 1:length(diff.year)){
        FBdata <- FBdata %>%
          mutate(year_group = ifelse(year >= diff.year[i], year_group+1, year_group))
      }
    }
  }

  if(is.null(labeling.year)){
    years <- unique(FBdata$year)
    labeling.year <- c(years[years%%5==0],max(years))
  }

  FBdata$Bratio <- FBdata[xlab_name][[1]] %>% as.numeric()
  FBdata$Fratio <- FBdata[ylab_name][[1]] %>% as.numeric()

  FBdata <- FBdata %>%
    mutate(year.label=ifelse(year%in%labeling.year,year,""))

  max.B <- max(c(FBdata$Bratio,xscale),na.rm=T)
  max.F <- max(c(FBdata$Fratio,yscale),na.rm=T)

  red.color <- "indianred1" # rgb(238/255,121/255,72/255)
  yellow.color <- "khaki1" # rgb(245/255,229/255,107/255)
  green.color <- "olivedrab2" # rgb(175/255,209/255,71/255) #"olivedrab2"#rgb(58/255,180/255,131/255)

  g4 <- ggplot(data=FBdata) +theme(legend.position="none")+
    geom_polygon(data=tibble(x=c(-1,1,1,-1),
                             y=c(-1,-1,1,1)),
                 aes(x=x,y=y),fill=yellow.color)+
    geom_polygon(data=tibble(x=c(1,20,20,1),
                             y=c(-1,-1,1,1)),
                 aes(x=x,y=y),fill=green.color)+
    geom_polygon(data=tibble(x=c(1,20,20,1),
                             y=c(1,1,20,20)),
                 aes(x=x,y=y),fill=yellow.color)+
    geom_polygon(data=tibble(x=c(-1,1,1,-1),
                             y=c(1,1,20,20)),
                 aes(x=x,y=y),fill=red.color) +
    geom_polygon(data=tibble(x=c(-1,1,1,-1),
                             y=c(-1,-1,1,1)),aes(x=x,y=y),fill=yellow.color)

  if(write.vline){
    g4 <- g4 + geom_vline(xintercept=c(1,limit.ratio,ban.ratio),color=refs.color,lty=c(unlist(format_type()[1,3]),unlist(format_type()[2,3]),unlist(format_type()[3,3])),lwd=0.7)+
      ggrepel::geom_label_repel(data=tibble(x=c(1,limit.ratio,ban.ratio),
                                            y=max.F*0.85,
                                            label=RP.label),
                                aes(x=x,y=y,label=label),
                                direction="x",nudge_y=max.F*0.9,size=11*0.282)
  }

  if(!is.null(beta)){
    ### HCRのプロット用の設定
    #Setting of the function to multiply current F for SSB
    multi2currF = function(x){
      if(x > limit.ratio) {multi2currF=beta}
      else if (x < ban.ratio) {multi2currF=0}
      else { multi2currF = beta*(x - ban.ratio)/(limit.ratio - ban.ratio) }
      return(multi2currF)
    }

    #Function setting for drawing.
    h=Vectorize(multi2currF)
    ####
    x.pos <- max.B*HCR.label.position[1]
    y.pos <- multi2currF(1.05)*HCR.label.position[2]
    g4 <- g4+stat_function(fun = h,lwd=1.5,color="black",n=5000)+
      annotate("text",x=x.pos,y=y.pos,
               label=str_c("漁獲管理規則\n(β=",beta,")"))
  }

  g4 <- g4 +
    geom_path(mapping=aes(x=Bratio,y=Fratio,group=year_group)) +
    geom_point(mapping=aes(x=Bratio,y=Fratio,group=year_group),shape=21,fill="white") +
    coord_cartesian(xlim=c(0,max.B*1.1),ylim=c(0,max.F*1.15),expand=0) +
    ylab("漁獲割合の比 (U/Umsy)") + xlab("親魚量の比 (SB/SBmsy)")  +
    ggrepel::geom_text_repel(#data=dplyr::filter(FBdata,year%in%labeling.year),
      aes(x=Bratio,y=Fratio,label=year.label),
      size=4,box.padding=0.5,segment.color="gray")

  if(ylab_name=="Fratio"){
    g4 <- g4 + ylab("漁獲圧の比 (F/Fmsy)")
  }

  g4 <- g4 + theme_SH()

  return(g4)
}


#' HCRを書く
#'
#' @param SBtarget 目標管理基準値
#' @param SBlim    限界管理基準値
#' @param SBlim    禁漁水準
#' @param Ftarget  Ftarget
#' @param is.text ラベルを記入するかどうか（FALSEにすると後で自分で書き換えられる）
#' @encoding UTF-8
#'
#' @export

plot_HCR <- function(SBtarget,SBlim,SBban,Ftarget,
                     Fcurrent=-1,
                     biomass.unit=1,
                     beta=0.8,col.multi2currf="black",col.SBtarget="#00533E",
                     col.SBlim="#edb918",col.SBban="#C73C2E",col.Ftarget="black",
                     col.betaFtarget="gray",is.text = TRUE,
                     RP.label=c("目標管理基準値","限界管理基準値","禁漁水準")){

  # Arguments; SBtarget,SBlim,SBban,Ftarget,beta,col.multi2currf,col.SBtarget,col.SBlim,col.SBban,col.Ftarget,col.betaFtarget.
  # col.xx means the line color for xx on figure.
  # beta and col.xx have default values.
  # Default setting for beta = 0.8, therefore define this as (beta <-0.8) outside this function if the beta-value changes frequently.
  # Default color setting for each parameter; Function(col.multi2currf="blue"), SBtarget(col.SBtarget = "green"), SBlimit(col.SBlim = "yellow"),SBban(col.SBban = "red"),Ftarget(col.Ftarget = "black"), β Ftarget(col.betaFtarget = "gray")

  junit <- c("","十","百","千","万")[log10(biomass.unit)+1]
  SBtarget <- SBtarget/biomass.unit
  SBlim <- SBlim/biomass.unit
  SBban <- SBban/biomass.unit

  #Setting of the function to multiply current F for SSB
  multi2currF = function(x){
    if(x > SBlim) {multi2currF=beta*Ftarget}
    else if (x < SBban) {multi2currF=0}
    else { multi2currF = (x - SBban)* beta*Ftarget/(SBlim - SBban) }
    return(multi2currF)
  }

  #Function setting for drawing.
  h=Vectorize(multi2currF)

  #Drawing of the funciton by ggplot2
  ggplct <- ggplot(data.frame(x = c(0,1.5*SBtarget),y= c(0,1.5*Ftarget)), aes(x=x)) +
    stat_function(fun = h,lwd=1.5,color=col.multi2currf, n=5000)
  g <- ggplct  + geom_vline(xintercept = SBtarget, size = 0.9, linetype = unlist(format_type()[1,3]), color = col.SBtarget) +
    geom_vline(xintercept = SBlim, size = 0.9, linetype = unlist(format_type()[2,3]), color = col.SBlim) +
    geom_vline(xintercept = SBban, size = 0.9, linetype = unlist(format_type()[3,3]), color = col.SBban) +
    geom_hline(yintercept = Ftarget, size = 0.9, linetype = "43", color = col.Ftarget) +
    geom_hline(yintercept = beta*Ftarget, size = 0.7, linetype = "43", color = col.betaFtarget) +
    labs(x = str_c("親魚量 (",junit,"トン)"),y = "漁獲圧の比(F/Fmsy)",color = "") +
    theme_bw(base_size=12)+
    theme(legend.position="none",panel.grid = element_blank())+
    stat_function(fun = h,lwd=1,color=col.multi2currf)

  if(Fcurrent>0){
    g <- g+geom_hline(yintercept = Fcurrent, size = 0.7, linetype = 1, color = "gray")+
      geom_label(label="Fcurrent", x=SBtarget*1.1, y=Fcurrent)

  }

  if(is.text) {
    RPdata <- tibble(RP.label=RP.label, value=c(SBtarget, SBlim, SBban), y=c(1.1,1.05,1.05))
    g <- g + ggrepel::geom_label_repel(data=RPdata,
                                       mapping=aes(x=value, y=y, label=RP.label),
                                       box.padding=0.5, nudge_y=0.05) +
      geom_label(label="Fmsy", x=SBtarget*1.3, y=Ftarget)+
      geom_label(label=str_c(beta,"Fmsy"), x=SBtarget*1.3, y=beta*Ftarget)+
      ylim(0,1.3)
  }

  return(g)

  #Drawing in a classical way
  # curve(h,
  #       xlim=c(0,2*SBtarget),  # range for x-axis is from 0 to 2*SBtarget
  #       ylim=c(0,1.2*Ftarget), # range for y-axis is from 0 to 1.2*Ftarget
  #       main="",
  #       xlab="SSB(×1000[t])",
  #       ylab="multipliyer to current F",
  #       lwd=2,
  #       col=col.multi2currf
  # )
  #Adding extention lines.
  # abline(v=SBtarget,lty=2,lwd=2,col=col.SBtarget)
  # abline(v=SBlim,lty=2,lwd=2,col=col.SBlim)
  # abline(v=SBban,lty=2,lwd=2,col=col.SBban)
  # abline(h=Ftarget,lty=2,col=col.Ftarget)
  # abline(h=beta*Ftarget,lty=3,col=col.betaFtarget)

  #Display legends at bottom right of the figure.
  # legend("bottomright",
  #        legend=c("SBtarget","SBlimit","SBban","Ftarget","β Ftarget"),
  #        lty=c(2,2,2,2,3),
  #        lwd=c(2,2,2,1,1),
  #        col=c(col.SBtarget, "yellow", "red","black","gray"),
  #        bty="n"
  # )

  #Setting each legend manually.
  # legend(SBtarget, 1.1*Ftarget,legend='SBtarget',bty="n")
  # legend(SBlim, 1.1*Ftarget, legend='SBlimit',bty="n")
  # legend(SBban, 1.1*Ftarget, legend='SBban',bty="n")
  # legend(0, Ftarget, legend='Ftarget',bty="n")
  # legend(0, beta*Ftarget, legend='β Ftarget',bty="n")

}

#' 縦軸が漁獲量のHCRを書く（traceの結果が必要）
#'
#' @param trace
#' @param fout 将来予測のアウトプット（finputがない場合)
#' @param Fvector Fのベクトル
#' @encoding UTF-8
#' @export

plot_HCR_by_catch <- function(trace,
                              fout0.8,
                              SBtarget,SBlim,SBban,Fmsy_vector,MSY,
                              M_vector,
                              biomass.unit=1,
                              beta=0.8,col.multi2currf="black",col.SBtarget="#00533E",
                              col.SBlim="#edb918",col.SBban="#C73C2E",col.Ftarget="black",
                              col.betaFtarget="gray",is.text = TRUE,
                              HCR_function_name="HCR_default",
                              Pope=TRUE,
                              RP.label=c("目標管理基準値","限界管理基準値","禁漁水準")){
  # 本当は途中までplot_HCRと統合させたい
  junit <- c("","十","百","千","万")[log10(biomass.unit)+1]
  biomass_comp <- trace %>% dplyr::select(starts_with("TB-mean-"))
  biomass_comp <- biomass_comp[,Fmsy_vector>0]
  M_vector <- M_vector[Fmsy_vector>0]
  Fmsy_vector <- Fmsy_vector[Fmsy_vector>0]

  calc_catch <- function(B, M, Fvec, Pope=TRUE){
    if(isTRUE(Pope)){
      total.catch <- B*(1-exp(-Fvec))*exp(-M/2)
    }
    else{
      total.catch <- B*(1-exp(-Fvec-M))*Fvec/(Fvec+M)
    }
    return(sum(total.catch))
  }

  n <- nrow(trace)
  HCR_function <- get(HCR_function_name)
  gamma <- HCR_function(as.numeric(trace$ssb.mean),
                       Blimit=rep(SBlim,n),Bban=rep(SBban,n),beta=rep(beta,n))
  F_matrix <- outer(gamma, Fmsy_vector)
  trace$catch_HCR <- purrr::map_dbl(1:nrow(trace), function(x)
    calc_catch(biomass_comp[x,],M_vector, F_matrix[x,], Pope=Pope))

  trace <- trace %>% dplyr::arrange(ssb.mean) %>%
    dplyr::filter(ssb.mean < SBtarget*1.5)

  g <- trace %>%
    ggplot()+
    geom_line(aes(x=ssb.mean/biomass.unit,y=catch_HCR/biomass.unit),lwd=1)+
    theme_SH()+
    geom_vline(xintercept = SBtarget/biomass.unit, size = 0.9, linetype = unlist(format_type()[1,3]), color = col.SBtarget) +
    geom_vline(xintercept = SBlim/biomass.unit, size = 0.9, linetype = unlist(format_type()[2,3]), color = col.SBlim) +
    geom_vline(xintercept = SBban/biomass.unit, size = 0.9, linetype = unlist(format_type()[3,3]), color = col.SBban) +
    #      geom_hline(yintercept = MSY/biomass.unit,color="gray")+
    xlab(str_c("親魚量 (",junit,"トン)"))+
    ylab(str_c("漁獲量 (",junit,"トン)"))

  if(is.text) {
    RPdata <- tibble(RP.label=RP.label, value=c(SBtarget, SBlim, SBban)/biomass.unit,
                     y=rep(max(trace$catch_HCR)*0.9,3)/biomass.unit)
    g <- g + ggrepel::geom_label_repel(data=RPdata,
                                       mapping=aes(x=value, y=y, label=RP.label),
                                       box.padding=0.5, nudge_y=1) #+
    # geom_label(label="MSY", x=SBtarget*1.4/biomass.unit, y=MSY/biomass.unit)
    #      geom_label(label=str_c(beta,"Fmsy"), x=SBtarget*1.3, y=beta*Ftarget)+
    #        ylim(0,1.3)
  }


}

#' F一定の場合で平衡状態になったときの統計量をx軸、y軸にプロットして比較する
#'
#'
#' 例えば、横軸を平均親魚量(ssb.mean)、縦軸を平均漁獲量(catch.mean)にすると漁獲量曲線が得られる。どの統計量がプロットできるかはest.MSYの返り値res_MSYの$trace以下の名前を参照(names(res_MSY$trace))。
#'
#' @param MSYlist est.MSYの返り値をリストにしたもの; 単独でも可
#' @param MSYname 凡例につけるMSYのケースの名前。MSYlistと同じ長さにしないとエラーとなる
#' @param x_axis_name x軸になにをとるか？("ssb.mean": 親魚の平均資源量, "fmulti": current Fに対する乗数、など)
#' @param y_axis_name y軸になにをとるか？("ssb.mean": 親魚の平均資源量, "catch.mean": 平均漁獲量、"rec.mean": 加入量の平均値など） get.statの返り値として出される値（またはMSYの推定結果のtrace内の表）に対応
#' @param plot_CI80 TRUE or FALSE, 平衡状態における信頼区間も追記する(現状では、縦軸が親魚量・漁獲量・加入尾数のときのみ対応)
#'
#' @examples
#' \dontrun{
#' data(res_MSY_HSL1)
#' data(res_MSY_HSL2)
#' MSY_list <- tibble::lst(res_MSY_HSL1, res_MSY_HSL2)
#' # 縦軸を漁獲量、横軸をFの大きさ
#' g1 <- compare_eq_stat(MSY_list,x_axis_name="fmulti",y_axis_name="catch.mean")
#' # 縦軸を親魚量にする
#' g2 <- compare_eq_stat(MSY_list,x_axis_name="fmulti",y_axis_name="ssb.mean")
#' # 縦軸を加入量
#' g3 <- compare_eq_stat(MSY_list,x_axis_name="fmulti",y_axis_name="rec.mean")
#' gridExtra::grid.arrange(g1,g2,g3,ncol=1)
#'
#' g3.withCI <- compare_eq_stat(MSY_list,x_axis_name="fmulti",y_axis_name="rec.mean",plot_CI80=TRUE)
#'
#' }
#'
#' @encoding UTF-8
#'
#' @export
#'

compare_eq_stat <- function(MSYlist,
                            x_axis_name="fmulti",
                            y_axis_name="catch.mean",
                            legend.position="top",
                            is_MSY_line=TRUE,
                            is.scale=FALSE,
                            MSYname=NULL,
                            plot_CI80=FALSE
){

  if(!is.null(MSYname)){
    if(length(MSYname)!=length(MSYlist)) stop("Length of MSYlist and MSYname is different")
    names(MSYlist) <- MSYname
  }
  if(isTRUE("summary" %in% names(MSYlist))) MSYlist <- list(MSYlist)

  data_yield <- purrr::map_dfr(MSYlist, function(x){
    x$trace %>% mutate(catch.order= rank(-catch.mean),
                       catch.max  = max(catch.mean)  ,
                       ssb.max    = max(ssb.mean))
  }
  ,.id="id")

  if(isTRUE(is.scale)) data_yield <- data_yield %>% mutate(catch.mean=catch.mean/catch.max,
                                                           ssb.mean=ssb.mean/ssb.max)

  g1 <- data_yield %>% ggplot()+
    geom_line(aes(x=get(x_axis_name), y=get(y_axis_name[1]), color=id))+
    theme_SH(legend.position=legend.position)+
    xlab(x_axis_name)+ylab(str_c(y_axis_name))+
    geom_vline(data=dplyr::filter(data_yield,catch.order==1),
               aes(xintercept=get(x_axis_name),color=id),lty=2)

  if(isTRUE(plot_CI80)){
    y_axis_name_L10 <- dplyr::case_when(
      y_axis_name == "catch.mean" ~ "catch.L10",
      y_axis_name == "ssb.mean"   ~ "ssb.L10",
      y_axis_name == "rec.mean"   ~ "rec.L10")
    y_axis_name_H10 <- dplyr::case_when(
      y_axis_name == "catch.mean" ~ "catch.H10",
      y_axis_name == "ssb.mean"   ~ "ssb.H10",
      y_axis_name == "rec.mean"   ~ "rec.H10")
    g1 <- g1 +
      geom_line(aes(x=get(x_axis_name), y=get(y_axis_name_L10), color=id),lty=2)+
      geom_line(aes(x=get(x_axis_name), y=get(y_axis_name_H10), color=id),lty=3)
  }

  return(g1)
}


#' 複数の管理基準値の推定結果を重ね書きする
#'
#' @param MSYlist est.MSYの返り値をリストにしたもの; 単独でも可
#' @param MSYname 凡例につけるMSYのケースの名前。MSYlistと同じ長さにしないとエラーとなる
#' @param legend.position 凡例の位置
#' @param yaxis
#'
#' @examples
#' \dontrun{
#' data(res_MSY)
#' MSY_list <- tibble::lst(res_MSY_HSL1, res_MSY_HSL2)
#' g1 <- compare_MSY(list(res_MSY, res_MSY))
#' }
#'
#' @encoding UTF-8
#'
#' @export
#'

compare_MSY <- function(MSYlist,
                        legend.position="top",
                        MSYname=NULL,
                        yaxis="Fref2Fcurrent"){

  if(!is.null(MSYname)){
    if(length(MSYname)!=length(MSYlist)) stop("Length of MSYlist and MSYname is different")
    names(MSYlist) <- MSYname
  }

  if(isTRUE("summary" %in% names(MSYlist))) MSYlist <- list(MSYlist)

  data_summary <- purrr::map_dfr(MSYlist, function(x) x$summary, .id="id")   %>%
    dplyr::filter(!is.na(RP.definition)) %>%
    mutate(label=stringr::str_c(id, RP.definition, sep="-")) %>%
    mutate(perSPR_rev=1-perSPR)

  g1 <- data_summary %>% ggplot()+
    geom_point(aes(x=SSB, y=get(yaxis), color=id))+
    ggrepel::geom_label_repel(aes(x=SSB, y=get(yaxis), color=id, label=label))+
    theme_SH(legend.position=legend.position)

  return(g1)
}

# test plot
#Fig_Fish_Manage_Rule(SBtarget,SBlim,SBban,Ftarget,col.multi2currf = "#093d86", col.SBtarget = "#00533E", col.SBlim = "#edb918",col.SBban = "#C73C2E",col.Ftarget = "#714C99", col.betaFtarget = "#505596")
# function;ruri-rio, sbtarget;moegi-iro, sblim;koki-ki; sbban;hi-iro, ftarget;sumire-iro, betaftarget;kikyou-iro

# MSE用 ----
#'
#' @export
#'

plot_bias_in_MSE <- function(fout, out="graph", error_scale="log", yrange=NULL){

  recruit_dat  <- convert_2d_future(fout$SR_MSE[,,"recruit"], name="Recruits", label="estimate") %>%
    rename(value_est=value)
  tmp <- convert_2d_future(fout$naa[1,,],            name="Recruits", label="true")
  recruit_dat$value_true <- tmp$value

  real_ABC_dat <- convert_2d_future(fout$HCR_realized[,,"wcatch"],    name="realABC", label="estimate") %>%
    rename(value_est=value)
  tmp  <- convert_2d_future(fout$SR_MSE[,,"real_true_catch"], name="realABC", label="true")
  real_ABC_dat$value_true <- tmp$value

  pseudo_ABC_dat <- convert_2d_future(fout$HCR_realized[,,"wcatch"],    name="pseudoABC", label="estimate") %>%
    rename(value_est=value)
  tmp  <- convert_2d_future(fout$SR_MSE[,,"pseudo_true_catch"], name="pseudoABC", label="true")
  pseudo_ABC_dat$value_true <- tmp$value

  alldat <- bind_rows(recruit_dat, real_ABC_dat, pseudo_ABC_dat) %>%
    dplyr::filter(value_est>0) %>%
    mutate(Relative_error_log=(log(value_est)-log(value_true))/log(value_true),
           Relative_error_normal=((value_est)-(value_true))/(value_true),
           Year=factor(year))

  g1 <- alldat %>% ggplot() +
    geom_boxplot(aes(x=Year,
                     #                         y=Relative_error_log,
                     y={if(error_scale=="log") Relative_error_log else Relative_error_normal},
                     fill=stat)) +
    facet_wrap(~stat,scale="free_y") +
    geom_hline(yintercept=0) +
    theme_SH() +
    ylab(str_c("Relative error (",error_scale,")"))

  if(!is.null(yrange)) g1 <- g1 + coord_cartesian(ylim=yrange)

  if(out=="graph"){
    return(g1)
  }
  else{
    return(alldat)
  }

}


#' 複数のKobe II tableの結果を重ね書きする
#'
#' @param kobeII_list betaをいろいろ変えたbeta.simulationの結果のオブジェクトのリスト
#' @param target_stat 目的とする統計量。デフォルトは"prob.over.ssbtarget"と"prob.over.ssblimit"。beta.simulationの返り値オブジェクトのリストの名前を必要に応じてとってくる
#' @param legend.position 凡例の位置
#' @param target_beta kobeII table から取り出すベータの値。デフォルトは0.8。
#'
#' @examples
#' \dontrun{
#' }
#'
#' @encoding UTF-8
#'
#' @export
#'
#'

compare_kobeII <- function(kobeII_list,
                           target_stat = c("prob.over.ssbtarget","prob.over.ssblimit",
                                           "prob.over.ssbban"),
                           legend_position = "top",
                           target_beta = 0.8){

  prob_data <- NULL

  for(i in 1:length(target_stat)){
    tmp_data <- purrr::map_dfr(kobeII_list[!is.na(kobeII_list)],
                                function(x){
                                  x[target_stat[i]][[1]] %>%
                                    dplyr::filter(beta==target_beta) %>%
                                    tidyr::pivot_longer(col=c(-HCR_name,-beta,-stat_name),
                                                        names_to="Year",values_to="Percent")},
                                .id="data_type") %>%
      mutate(Year=as.numeric(Year))
    prob_data <- rbind(prob_data,tmp_data)
  }

  g1 <- prob_data %>% ggplot() +
    geom_line(aes(x=Year,y=Percent,color=data_type))+
    facet_wrap(.~stat_name,scale="free_y")+
    theme_SH(legend="top")+ylim(0,NA)

  return(g1)
}

#' Plot SPR figure
#'
#' @inheritParams make_stock_table
#' @param years 全ての年をプロットしない場合、\code{years=1985:2011}のようにプロットする年を指定する
#' @export
plot_sprypr <- function(result_vpa, type, years=NULL) {
  df <- get.SPR(result_vpa, target.SPR = 30, Fmax = 10)$ysdata
  if (type == "perspr") {
    df$Year <- as.numeric(rownames(df))
    if (!is.null(years)) df = df %>% dplyr::filter(Year %in% years)
    plot <- df %>%
      ggplot(aes(Year, perSPR)) +
      geom_line() +
      geom_point() +
      scale_y_continuous(trans = "reverse",
                         breaks = seq(100, 0, -20),
                         limits = c(100, 0)) +
      xlab("年") +
      ylab("%SPR")
  }
  force(plot)
}



#' ARの効果を差し引いて再生産関係をプロットする
#'
#' @param res fit.SRから返されるオブジェクト
#' @examples
#' \dontrun{
#'
#' }
#'
#' @encoding UTF-8
#' @export


plot_SR_AReffect <- function(res){
    SRdata <- as_tibble(res$input$SRdata)

#    pred0 <- exp(log(SRdata$R)-res$resid)
#    pred1 <- exp(log(pred0)+res$pars$rho*c(0,res$resid[1:(length(res$resid)-1)]))


    SRdata2 <- bind_rows(SRdata %>% mutate(Data="Observed"),
                         SRdata %>% mutate(R=exp(log(SRdata$R)- # res$resid -
                                                   res$pars$rho*c(0,res$resid[1:(length(res$resid)-1)])),
                                           Data="Observed_without_AR"),
                         )


    SRdata2 <- SRdata2 %>%
      mutate(SB=SSB, Data=factor(Data,levels=c("Observed","Observed_without_AR")))

#    if(obs_change==FALSE){
#      g <- ggplot(SRdata)+
#        geom_point(aes(x=SSB, y=R))+
#        geom_line(aes(x=res$pred[,1],y=res$pred[,2]),data=res$pred)+
#        geom_path(aes(x=SSB, y=pred1), color="red",linetype=1)+
#        xlab("SB")+ylab("R")+theme_bw()
#    }
#    else{
      g <- ggplot(SRdata2)+
        geom_line(aes(x=res$pred[,1],y=res$pred[,2]),data=res$pred)+
        geom_point(aes(x=SB, y=R, color=Data))+
        xlab("SB")+ylab("Number of reruits")+
        theme_SH(legend.position="top") +
        scale_color_manual(values=c("black", "red"))
#    }
    return(g)
}

#' @encoding UTF-8
#' @export


plot_worm <- function(kobe_data){

#  HCR_selection <- read_csv("HCR_selection.csv") %>%
#    rename(HCR_category="HCR category",
#           num=`Serial number`,
#           stock=`Appied stock`) %>%
#    dplyr::filter(!is.na(HCR_category)) %>%
#    mutate(size=ifelse(Type=="B"|Type=="SS", 3, 1)) %>%
#    mutate(size=factor(size))

    mean_data <- bind_rows(kobe_data$catch.mean,
                           kobe_data$ssb.mean  ,
                           kobe_data$ssb.lower05percent) %>%
      pivot_longer(cols=c(-HCR_name,-beta,-stat_name)) %>%
      rename(Year=name) %>%
      mutate(Year=as.numeric(Year), MT=value/1000) #%>%
#      left_join(HCR_selection) %>%
#      left_join(tibble(stat_name =c("catch.mean","ssb.mean","ssb.ci10"),
#                       stat_name2=c("Catch (average)",     "SB (average)", "SB (L10%)"))) %>%
#      mutate(Type=factor(Type, levels=c("B","S","SS","A"))) %>%
#      dplyr::filter(use==1)

    g_worm <- mean_data %>%
      ggplot() +
      geom_line(aes(x=Year, y=MT, color=HCR_name, group=HCR_name),
                alpha=0.8) +
      ylim(0,NA) +
      facet_wrap(.~stat_name, scale="free_y") + theme_SH(base_size=14) +
      coord_cartesian(ylim=c(0,NA)) +
      ylab("1000 MT") +
#      scale_color_manual(values=c(1,gray(0.2),2,3)) +
#      scale_size_manual(values=c(1,0.5,1,0.5)) +
        theme(legend.position="bottom")

    g_worm

}
