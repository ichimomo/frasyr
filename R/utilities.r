col.SBtarget    <- "#00533E"
col.SBlim       <- "#edb918"
col.SBlimit     <- "#edb918"
col.SBban       <- "#C73C2E"
col.Ftarget     <- "#714C99"
col.betaFtarget <- "#505596"

pt1             <- 0.3528


convert_df <- function(df,name){
    df %>%
        as_tibble %>%  
        mutate(age = as.numeric(rownames(df))) %>% 
        gather(key=year, value=value, -age, convert=TRUE) %>%
        group_by(year) %>%
#        summarise(value=sum(value)) %>%
        mutate(type="VPA",sim="s0",stat=name)    
}

#' future.vpaの結果オブジェクトをtibble形式に変換する関数
#'
#' @param fout future.vpaの結果のオブジェクト
#' 
#' @encoding UTF-8
#' @export
#' 

convert_future_table <- function(fout,label="tmp"){
    ssb <- fout$vssb %>%
        as_tibble %>%
        mutate(year=rownames(fout$vssb)) %>%
        gather(key=sim, value=value, -year, convert=TRUE) %>%
        mutate(year=as.numeric(year),stat="SSB",label=label)

    catch <- fout$vwcaa %>%
        as_tibble %>%
        mutate(year=rownames(fout$vssb)) %>%
        gather(key=sim, value=value, -year, convert=TRUE) %>%
        mutate(year=as.numeric(year),stat="catch",label=label)

    biomass <- fout$vbiom %>%
        as_tibble %>%
        mutate(year=rownames(fout$vbiom)) %>%
        gather(key=sim, value=value, -year, convert=TRUE) %>%
        mutate(year=as.numeric(year),stat="biomass",label=label)

    U_table <- fout$vwcaa/fout$vbiom %>% as_tibble
    U_table <- U_table %>%
        mutate(year=rownames(fout$vbiom)) %>%
        gather(key=sim, value=value, -year, convert=TRUE) %>%
        mutate(year=as.numeric(year),stat="U",label=label)
    
    alpha_value <- fout$alpha %>%
        as_tibble %>%
        mutate(year=rownames(fout$alpha)) %>%
        gather(key=sim, value=value, -year, convert=TRUE) %>%
        mutate(year=as.numeric(year),stat="alpha",label=label)

    if(is.null(fout$Fsakugen)) fout$Fsakugen <- -(1-fout$faa[1,,]/fout$input$res0$Fc.at.age[1])
    Fsakugen <- fout$Fsakugen %>%
        as_tibble %>%
        mutate(year=rownames(fout$Fsakugen)) %>%
        gather(key=sim, value=value, -year, convert=TRUE) %>%
        mutate(year=as.numeric(year),stat="Fsakugen",label=label)

    Fsakugen_ratio <- Fsakugen %>%
        mutate(value=value+1)
    Fsakugen_ratio$stat <- "Fsakugen_ratio"

    if(is.null(fout$recruit)) fout$recruit <- fout$naa[1,,]
    Recruitment <- fout$recruit %>%                                    #追加
        as_tibble %>%                                                   #追加
        mutate(year=rownames(fout$recruit)) %>%                             #追加
        gather(key=sim, value=value, -year, convert=TRUE) %>%           #追加
        mutate(year=as.numeric(year),stat="Recruitment",label=label)
    
    bind_rows(ssb,catch,biomass,alpha_value,Fsakugen,Fsakugen_ratio,Recruitment, U_table)
}
        
    
convert_vector <- function(vector,name){
    vector %>%
        as_tibble %>%  
        mutate(year = as.integer(names(vector))) %>% 
        mutate(type="VPA",sim="s0",stat=name,age=NA) 
} 

#' VPAの結果オブジェクトをtibble形式に変換する関数
#'
#' @param vpares vpaの結果のオブジェクト
#' @encoding UTF-8
#' 
#'
#' @export

convert_vpa_tibble <- function(vpares){

  if (is.null(vpares$input$dat$waa.catch)) {
    total.catch <- colSums(vpares$input$dat$caa*vpares$input$dat$waa,na.rm=T)
  } else {
    total.catch <- colSums(vpares$input$dat$caa*vpares$input$dat$waa.catch,na.rm=T)
  }
    U <- total.catch/colSums(vpares$baa)

    SSB <- convert_vector(colSums(vpares$ssb,na.rm=T),"SSB") %>%
        dplyr::filter(value>0&!is.na(value))
    Biomass <- convert_vector(colSums(vpares$baa,na.rm=T),"biomass") %>%
        dplyr::filter(value>0&!is.na(value))
    FAA <- convert_df(vpares$faa,"fishing_mortality") %>%
        dplyr::filter(value>0&!is.na(value))
    Recruitment <- convert_vector(colSums(vpares$naa[1,,drop=F]),"Recruitment") %>%
        dplyr::filter(value>0&!is.na(value))    
    
    all_table <- bind_rows(SSB,
                           Biomass,
                           convert_vector(U[U>0],"U"),
                           convert_vector(total.catch[total.catch>0],"catch"),
                           convert_df(vpares$naa,"fish_number"),
                           FAA, 
                           convert_df(vpares$input$dat$waa,"weight"),
                           convert_df(vpares$input$dat$maa,"maturity"),
                           convert_df(vpares$input$dat$caa,"catch_number"),
                           Recruitment)
}


#' 再生産関係をプロットする関数
#'
#' @param SR_result fit.SRの結果のオブジェクト
#' @encoding UTF-8
#'
#' @export
#' 


SRplot_gg <- plot.SR <- function(SR_result,refs=NULL,xscale=1000,xlabel="千トン",yscale=1,ylabel="尾",
                      labeling.year=NULL,add.info=TRUE){
#    require(tidyverse,quietly=TRUE)    
    #    require(ggrepel)

    if(is.null(refs$Blimit) && !is.null(refs$Blim)) refs$Blimit <- refs$Blim

    if (SR_result$input$SR=="HS") SRF <- function(SSB,a,b) (ifelse(SSB*xscale>b,b*a,SSB*xscale*a))/yscale
    if (SR_result$input$SR=="BH") SRF <- function(SSB,a,b) (a*SSB*xscale/(1+b*SSB*xscale))/yscale
    if (SR_result$input$SR=="RI") SRF <- function(SSB,a,b) (a*SSB*xscale*exp(-b*SSB*xscale))/yscale
    
    SRdata <- as_tibble(SR_result$input$SRdata) %>%
        mutate(type="obs")
    SRdata.pred <- as_tibble(SR_result$pred) %>%
        mutate(type="pred",year=NA)    
    alldata <- bind_rows(SRdata,SRdata.pred) %>%
        mutate(R=R/yscale,SSB=SSB/xscale)
    ymax <- max(alldata$R)
    year.max <- max(alldata$year,na.rm=T)
    tmp <- 1950:2030
    if(is.null(labeling.year)) labeling.year <- c(tmp[tmp%%5==0],year.max)
    alldata <- alldata %>% mutate(pick.year=ifelse(year%in%labeling.year,year,""))

    g1 <- ggplot(data=alldata,mapping=aes(x=SSB,y=R)) +
#        geom_line(data=dplyr::filter(alldata,type=="pred"),
#                      aes(y=R,x=SSB),color="deepskyblue3",lwd=1.3) +
        stat_function(fun=SRF,args=list(a=SR_result$pars$a,
                                        b=SR_result$pars$b),color="deepskyblue3",lwd=1.3)+
    geom_path(data=dplyr::filter(alldata,type=="obs"),
                  aes(y=R,x=SSB),color=1) +
        geom_point(data=dplyr::filter(alldata,type=="obs"),
                   aes(y=R,x=SSB),shape=21,fill="white") +
#        scale_shape_discrete(solid=T)+        
#        geom_label_repel(data=dplyr::filter(alldata,type=="obs" & (year%%10==0|year==year.max)),
#                         aes(y=R,x=SSB,label=year),
    #                         size=3,box.padding=3,segment.color="black") +
    #        geom_text_repel(aes(y=R,x=SSB,label=pickyear)) +
    ggrepel::geom_text_repel(data=dplyr::filter(alldata,type=="obs"),
                    segment.alpha=0.5,nudge_y=5,
                    aes(y=R,x=SSB,label=pick.year)) +                
        theme_bw(base_size=14)+
    theme(legend.position = 'none') +
    theme(panel.grid = element_blank()) +
        xlab(str_c("親魚資源量 (",xlabel,")"))+
        ylab(str_c("加入尾数 (",ylabel,")"))+        
    coord_cartesian(ylim=c(0,ymax*1.05),expand=0)

    if(add.info){
        g1 <- g1+labs(caption=str_c("関数形: ",SR_result$input$SR,", 自己相関: ",SR_result$input$AR,
                           ", 最適化法",SR_result$input$method,", AICc: ",round(SR_result$AICc,2)))
    }

    if(!is.null(refs)){
        g1 <- g1+geom_vline(xintercept=c(refs$Bmsy/xscale,refs$Blimit/xscale,refs$Bban/xscale),
                            linetype=2,
                            col=c(col.SBtarget,col.SBlimit,col.SBban))
    }
    g1
}

get.trace <- function(trace){
    trace <- trace  %>% as_tibble() %>%
        select(starts_with("TC-mean"),ssb.mean,fmulti,catch.CV) %>%
        mutate(label=as.character(1:nrow(.)))

    trace <- trace %>% gather(value=value,key=age,-label,-fmulti,-ssb.mean,-catch.CV) %>%
        mutate(age=str_extract(age, "(\\d)+")) %>%
        mutate(age=factor(age)) %>%
        mutate(age=fct_reorder(age,length(age):1))
    return(trace)
}

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
                       age.label.ratio=0.9, # 年齢のラベルを入れる位置（xの最大値からの割合)
                       family = "JP1",
                       ylim.scale=1.2,future=NULL,
                       future.replicate=NULL,
                       past=NULL,future.name=NULL){
    
    junit <- c("","十","百","千","万")[log10(biomass.unit)+1]
   
    if ("trace" %in% names(MSY_obj)) {
      trace.msy <- MSY_obj$trace
    } else {
      trace.msy <- MSY_obj
    }
        
#    require(tidyverse,quietly=TRUE)
#    require(ggrepel)    

    trace <- get.trace(trace.msy) %>%
        mutate("年齢"=age,ssb.mean=ssb.mean/biomass.unit,value=value/biomass.unit)

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

    xmax <- max(trace$ssb.mean,na.rm=T)
    age.label.position <- trace$ssb.mean[which.min((trace$ssb.mean-xmax*xlim.scale*age.label.ratio)^2)]
    age.label <- trace %>% dplyr::filter(round(age.label.position,1)==round(ssb.mean,1))%>%
        mutate(cumcatch=cumsum(value)-value/2)%>%
        mutate(age=as.numeric(as.character(age)))
    age.label <- age.label %>%
        mutate(age_name=str_c(age,ifelse(age.label$age==max(age.label$age),"+",""),"歳"))
   
    g1 <- g1 + geom_area(aes(x=ssb.mean,y=value,fill=年齢),col="black",alpha=0.5,lwd=1*0.3528) +
#    geom_line(aes(x=ssb.mean,y=catch.CV,fill=age)) +
#    scale_y_continuous(sec.axis = sec_axis(~.*5, name = "CV catch"))+
    scale_fill_brewer() +
    theme_bw() +
    theme(legend.position = 'none') +        
#    geom_point(data=refs_base,aes(y=Catch,x=SSB,shape=refs.label,color=refs.label),size=4)+
    #形は塗りつぶしができる形にすること
    scale_shape_manual(values = c(21, 24,5,10)) +
    #塗りつぶし色を指定する
    scale_color_manual(values = c("green", "pink","orange","yellow"))+
    theme(panel.grid = element_blank(),axis.text=element_text(color="black")) +
    coord_cartesian(xlim=c(0,xmax*xlim.scale),
                    ylim=c(0,ymax*ylim.scale),expand=0) +
    geom_text(data=age.label,
              mapping = aes(y = cumcatch, x = ssb.mean, label = age_name),
                            family = family
              ) +
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
        tmpdata <- NULL
        for(j in 1:length(future)){
            if(is.null(future.replicate)){
               tmpdata <- bind_rows(tmpdata,
                   tibble(
                   year        =as.numeric(rownames(future[[j]]$vssb)),
                   ssb.future  =apply(future[[j]]$vssb[,-1],1,mean)/biomass.unit,
                   catch.future=apply(future[[j]]$vwcaa[,-1],1,mean)/biomass.unit,
                   scenario=future.name[j]))
            }
            else{
               tmpdata <- bind_rows(tmpdata,
                   tibble(
                   year        =as.numeric(rownames(future[[j]]$vssb)),
                   ssb.future  =future[[j]]$vssb[,future.replicate[j]]/biomass.unit,
                   catch.future=future[[j]]$vwcaa[,future.replicate[j]]/biomass.unit,
                   scenario=future.name[j]))                
            }
            tmpdata <- tmpdata %>% group_by(scenario)
            g1 <- g1 +
                geom_path(data=tmpdata,
                          mapping=aes(x       =ssb.future,
                                      y       = catch.future,
                                      linetype=factor(scenario),
                                      color   = factor(scenario)),
                          lwd=1)+
                geom_point(data=tmpdata,
                           mapping=aes(x    =ssb.future,
                                       y    =catch.future,
                                       shape=factor(scenario),
                                       color=factor(scenario)),
                           size=3)

            
        }
    }
    
    if(!is.null(past)){
      catch.past = unlist(colSums(past$input$dat$caa*past$input$dat$waa)/biomass.unit)
      if (past$input$last.catch.zero && !is.null(future)) {
        catch.past[length(catch.past)] = apply(future[[1]]$vwcaa[,-1],1,mean)[1]/biomass.unit
      }
        tmpdata <- tibble(
            year      =as.numeric(colnames(past$ssb)),
            ssb.past  =unlist(colSums(past$ssb))/biomass.unit,
            catch.past=catch.past
        )

        g1 <- g1 +
#            geom_point(data=tmpdata,mapping=aes(x=ssb.past,y=catch.past,
#                                                alpha=year),shape=2) +
    geom_path(data=tmpdata,
              mapping=aes(x=ssb.past,y=catch.past),
              color=col.SBban,lwd=1,alpha=0.9)
    }
    
    if(isTRUE(lining)){
#        ylim.scale.factor <- rep(c(0.94,0.97),ceiling(length(refs.label)/2))[1:length(refs.label)]
        g1 <- g1 + geom_vline(xintercept=refs_base$SSB,lty="41",lwd=0.6,color=refs.color)+
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
#                            size=4,box.padding=0.5,segment.color=1,
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

#' 管理基準値の表を作成する
#'
#' @param refs_base est.MSYから得られる管理基準値の表
#' @encoding UTF-8
#'
#' @export
#'

make_RP_table <- function(refs_base){
#    require(formattable)
#    require(tidyverse,quietly=TRUE)
    table_output <- refs_base %>%
        select(-RP_name) %>% # どの列を表示させるか選択する
        # 各列の有効数字を指定
        mutate(SSB=round(SSB,-floor(log10(min(SSB)))),
               SSB2SSB0=round(SSB2SSB0,2),                              
               Catch=round(Catch,-floor(log10(min(Catch)))),
               Catch.CV=round(Catch.CV,2),
               U=round(U,2),
               Fref2Fcurrent=round(Fref2Fcurrent,2)) %>%
        rename("管理基準値"=RP.definition,"親魚資源量"=SSB,"B0に対する比"=SSB2SSB0,
               "漁獲量"=Catch,"漁獲量の変動係数"=Catch.CV,"漁獲率"=U,"努力量の乗数"=Fref2Fcurrent)
    
   table_output  %>%    
        # 表をhtmlで出力
        formattable::formattable(list(親魚資源量=color_bar("olivedrab"),
                                  漁獲量=color_bar("steelblue"),
                              漁獲率=color_bar("orange"),
                              努力量の乗数=color_bar("tomato")))

#    return(table_output)
    
}

#' 管理基準値表から目的の管理基準値を取り出す関数
#'
#' @param refs_base est.MSYから得られる管理基準値の表
#' @param RP_name 取り出したい管理基準値の名前
#' @encoding UTF-8
#'
#' @export
#' 

derive_RP_value <- function(refs_base,RP_name){
#    refs_base %>% dplyr::filter(RP.definition%in%RP_name)
#    subset(refs_base,RP.definition%in%RP_name)
    refs_base[refs_base$RP.definition%in%RP_name,]    
}

#' Kobe II matrixを計算するための関数
#'
#' @param fres_base future.vpaの結果のオブジェクト
#' @param refs_base est.MSYから得られる管理基準値の表
#' @encoding UTF-8
#'
#' @export

calc_kobeII_matrix <- function(fres_base,
                              refs_base,
                              Btarget=c("Btarget0"),
                              Blimit=c("Blimit0"),
#                              Blow=c("Blow0"),
                              Bban=c("Bban0"),
                              year.lag=0,
                              beta=seq(from=0.5,to=1,by=0.1)){
#    require(tidyverse,quietly=TRUE)    
# HCRの候補を網羅的に設定
#    HCR_candidate1 <- expand.grid(
#        Btarget_name=refs_base$RP.definition[str_detect(refs_base$RP.definition,Btarget)],
#        Blow_name=refs_base$RP.definition[str_detect(refs_base$RP.definition,Blow)],    
#        Blimit_name=refs_base$RP.definition[str_detect(refs_base$RP.definition,Blimit)],
#        Bban_name=refs_base$RP.definition[str_detect(refs_base$RP.definition,Bban)],
    #        beta=beta)

    refs.unique <- unique(c(Btarget,Blimit,Bban))
    tmp <- !refs.unique%in%refs_base$RP.definition    
    if(sum(tmp)>0) stop(refs.unique[tmp]," does not appear in column of RP.definition\n")

    HCR_candidate1 <- expand.grid(
        Btarget_name=derive_RP_value(refs_base,Btarget)$RP.definition,
#        Blow_name=derive_RP_value(refs_base,Blow)$RP.definition,    
        Blimit_name=derive_RP_value(refs_base,Blimit)$RP.definition,
        Bban_name=derive_RP_value(refs_base,Bban)$RP.definition,
        beta=beta)    

    HCR_candidate2 <- expand.grid(
        Btarget=derive_RP_value(refs_base,Btarget)$SSB,
#        Blow=derive_RP_value(refs_base,Blow)$SSB,    
        Blimit=derive_RP_value(refs_base,Blimit)$SSB,    
        Bban=derive_RP_value(refs_base,Bban)$SSB,   
        beta=beta) %>% select(-beta)

    HCR_candidate <- bind_cols(HCR_candidate1,HCR_candidate2) %>% as_tibble()
    
    HCR_candidate <- refs_base %>% #dplyr::filter(str_detect(RP.definition,Btarget)) %>%
        dplyr::filter(RP.definition%in%Btarget) %>%
        mutate(Btarget_name=RP.definition,Fmsy=Fref2Fcurrent) %>%
        select(Btarget_name,Fmsy) %>%
        left_join(HCR_candidate) %>%
        arrange(Btarget_name,Blimit_name,Bban_name,desc(beta))
    
    HCR_candidate$HCR_name <- str_c(HCR_candidate$Btarget_name,
                                    HCR_candidate$Blimit_name,
                                    HCR_candidate$Bban_name,sep="-")
    fres_base$input$outtype <- "FULL"
    kobeII_table <- HCR.simulation(fres_base$input,HCR_candidate,year.lag=year.lag)

    cat(length(unique(HCR_candidate$HCR_name)), "HCR is calculated: ",
        unique(HCR_candidate$HCR_name),"\n")

    kobeII_data <- left_join(kobeII_table,HCR_candidate)
    return(kobeII_data)
}

#'
#' @export
#' 

make_kobeII_table0 <- function(kobeII_data,
                              res_vpa,
                              year.catch,
                              year.ssb,                              
                              year.Fsakugen,
                              year.ssbtarget,
                              year.ssblimit,
                              year.ssbban,
                              year.ssbmin,
                              year.ssbmax,                              
                              year.aav){
    # 平均漁獲量
    (catch.table <- kobeII.data %>%
         dplyr::filter(year%in%year.catch,stat=="catch") %>% # 取り出す年とラベル("catch")を選ぶ
         group_by(HCR_name,beta,year) %>%
         summarise(catch.mean=round(mean(value))) %>%  # 値の計算方法を指定（漁獲量の平均ならmean(value)）
         # "-3"とかの値で桁数を指定
         spread(key=year,value=catch.mean) %>% ungroup() %>%
         arrange(HCR_name,desc(beta)) %>% # HCR_nameとbetaの順に並び替え
         mutate(stat_name="catch.mean"))

    # 平均親魚
    (ssb.table <- kobeII.data %>%
         dplyr::filter(year%in%year.ssb,stat=="SSB") %>% 
         group_by(HCR_name,beta,year) %>%
         summarise(ssb.mean=round(mean(value))) %>%  
         spread(key=year,value=ssb.mean) %>% ungroup() %>%
         arrange(HCR_name,desc(beta)) %>% # HCR_nameとbetaの順に並び替え
         mutate(stat_name="ssb.mean"))    

    # 1-currentFに乗じる値=currentFからの努力量の削減率の平均値（実際には確率分布になっている）
    (Fsakugen.table <- kobeII.data %>%
         dplyr::filter(year%in%year.Fsakugen,stat=="Fsakugen") %>% # 取り出す年とラベル("catch")を選ぶ
         group_by(HCR_name,beta,year) %>%
         summarise(Fsakugen=round(mean(value),2)) %>%
         spread(key=year,value=Fsakugen) %>% ungroup() %>%
         arrange(HCR_name,desc(beta)) %>% # HCR_nameとbetaの順に並び替え
         mutate(stat_name="Fsakugen"))

    # SSB>SSBtargetとなる確率
    ssbtarget.table <- kobeII.data %>%
        dplyr::filter(year%in%year.ssbtarget,stat=="SSB") %>%
        group_by(HCR_name,beta,year) %>%
        summarise(ssb.over=round(100*mean(value>Btarget))) %>%
        spread(key=year,value=ssb.over) %>%
        ungroup() %>%
        arrange(HCR_name,desc(beta))%>%
        mutate(stat_name="Pr(SSB>SSBtarget)")

    # SSB>SSBlimとなる確率
    ssblimit.table <- kobeII.data %>%
        dplyr::filter(year%in%year.ssblimit,stat=="SSB") %>%
        group_by(HCR_name,beta,year) %>%
        summarise(ssb.over=round(100*mean(value>Blimit))) %>%
        spread(key=year,value=ssb.over)%>%
        ungroup() %>%
        arrange(HCR_name,desc(beta))%>%
        mutate(stat_name="Pr(SSB>SSBlim)")

    # SSB>SSBbanとなる確率
    ssbban.table <- kobeII.data %>%
        dplyr::filter(year%in%year.ssbban,stat=="SSB") %>%
        group_by(HCR_name,beta,year) %>%
        summarise(ssb.over=round(100*mean(value>Bban))) %>%
        spread(key=year,value=ssb.over)%>%
        ungroup() %>%
        arrange(HCR_name,desc(beta))%>%
        mutate(stat_name="Pr(SSB>SSBban)")    

    # SSB>SSBmin(過去最低親魚量を上回る確率)
    ssb.min <- min(unlist(colSums(res_vpa$ssb)))
    ssbmin.table <- kobeII.data %>%
        dplyr::filter(year%in%year.ssbmin,stat=="SSB") %>%
        group_by(HCR_name,beta,year) %>%
        summarise(ssb.over=round(100*mean(value>ssb.min))) %>%
        spread(key=year,value=ssb.over)%>%
        ungroup() %>%
        arrange(HCR_name,desc(beta))%>%
        mutate(stat_name="Pr(SSB>SSBmin)")

    # SSB>SSBmax(過去最低親魚量を上回る確率)
    ssb.max <- max(unlist(colSums(res_vpa$ssb)))
    ssbmax.table <- kobeII.data %>%
        dplyr::filter(year%in%year.ssbmax,stat=="SSB") %>%
        group_by(HCR_name,beta,year) %>%
        summarise(ssb.over=round(100*mean(value>ssb.max))) %>%
        spread(key=year,value=ssb.over)%>%
        ungroup() %>%
        arrange(HCR_name,desc(beta))%>%
        mutate(stat_name="Pr(SSB>SSBmax)")    

    # オプション: Catch AAV mean 
    calc.aav <- function(x)sum(abs(diff(x)))/sum(x[-1])
    catch.aav.table <- kobeII.data %>%
        dplyr::filter(year%in%year.aav,stat=="catch") %>%
        group_by(HCR_name,beta,sim) %>%
        dplyr::summarise(catch.aav=(calc.aav(value))) %>%
        group_by(HCR_name,beta) %>%
        summarise(catch.aav.mean=mean(catch.aav)) %>%
        arrange(HCR_name,desc(beta))%>%
        mutate(stat_name="catch.csv (recent 5 year)")

    res_list <- list(average.catch   = catch.table,
                     average.ssb     = ssb.table,
                     prob.ssbtarget  = ssbtarget.table,
                     prob.ssblimit   = ssblimit.table,
                     prob.ssbban     = ssbban.table,                     
                     prob.ssbmin     = ssbmin.table,
                     prob.ssbmax     = ssbmax.table,                     
                     catch.aav       = catch.aav.table)    
    return(res_list)
                
}

#'
#' @export
#' 

make_kobeII_table <- function(kobeII_data,
                              res_vpa,
                              year.catch,
                              year.ssb,                              
                              year.Fsakugen,
                              year.ssbtarget,
                              year.ssblimit,
                              year.ssbban,
                              year.ssbmin,
                              year.ssbmax,                              
                              year.aav,
                              Btarget=0,
                              Blimit=0,
                              Bban=0){
    # 平均漁獲量
    (catch.mean <- kobeII.data %>%
         dplyr::filter(year%in%year.catch,stat=="catch") %>% # 取り出す年とラベル("catch")を選ぶ
         group_by(HCR_name,beta,year) %>%
         summarise(catch.mean=mean(value)) %>%  # 値の計算方法を指定（漁獲量の平均ならmean(value)）
         # "-3"とかの値で桁数を指定
         spread(key=year,value=catch.mean) %>% ungroup() %>%
         arrange(HCR_name,desc(beta)) %>% # HCR_nameとbetaの順に並び替え
         mutate(stat_name="catch.mean"))

    # 平均親魚
    (ssb.mean <- kobeII.data %>%
         dplyr::filter(year%in%year.ssb,stat=="SSB") %>% 
         group_by(HCR_name,beta,year) %>%
         summarise(ssb.mean=mean(value)) %>%  
         spread(key=year,value=ssb.mean) %>% ungroup() %>%
         arrange(HCR_name,desc(beta)) %>% # HCR_nameとbetaの順に並び替え
         mutate(stat_name="ssb.mean"))

    # 親魚, 下10%
    (ssb.ci10 <- kobeII.data %>%
         dplyr::filter(year%in%year.ssb,stat=="SSB") %>% 
         group_by(HCR_name,beta,year) %>%
         summarise(ssb.ci10=quantile(value,probs=0.1)) %>%  
         spread(key=year,value=ssb.ci10) %>% ungroup() %>%
         arrange(HCR_name,desc(beta)) %>% # HCR_nameとbetaの順に並び替え
         mutate(stat_name="ssb.ci10"))

    # 親魚, 上10%
    (ssb.ci90 <- kobeII.data %>%
         dplyr::filter(year%in%year.ssb,stat=="SSB") %>% 
         group_by(HCR_name,beta,year) %>%
         summarise(ssb.ci90=quantile(value,probs=0.9)) %>%  
         spread(key=year,value=ssb.ci90) %>% ungroup() %>%
         arrange(HCR_name,desc(beta)) %>% # HCR_nameとbetaの順に並び替え
         mutate(stat_name="ssb.ci90"))            

    # 1-currentFに乗じる値=currentFからの努力量の削減率の平均値（実際には確率分布になっている）
    (Fsakugen.table <- kobeII.data %>%
         dplyr::filter(year%in%year.Fsakugen,stat=="Fsakugen") %>% # 取り出す年とラベル("catch")を選ぶ
         group_by(HCR_name,beta,year) %>%
         summarise(Fsakugen=round(mean(value),2)) %>%
         spread(key=year,value=Fsakugen) %>% ungroup() %>%
         arrange(HCR_name,desc(beta)) %>% # HCR_nameとbetaの順に並び替え
         mutate(stat_name="Fsakugen"))

    # SSB>SSBtargetとなる確率
    ssbtarget.table <- kobeII.data %>%
        dplyr::filter(year%in%year.ssbtarget,stat=="SSB") %>%
        group_by(HCR_name,beta,year) %>%
        summarise(ssb.over=round(100*mean(value>Btarget))) %>%
        spread(key=year,value=ssb.over) %>%
        ungroup() %>%
        arrange(HCR_name,desc(beta))%>%
        mutate(stat_name="Pr(SSB>SSBtarget)")

    # SSB>SSBlimとなる確率
    ssblimit.table <- kobeII.data %>%
        dplyr::filter(year%in%year.ssblimit,stat=="SSB") %>%
        group_by(HCR_name,beta,year) %>%
        summarise(ssb.over=round(100*mean(value>Blimit))) %>%
        spread(key=year,value=ssb.over)%>%
        ungroup() %>%
        arrange(HCR_name,desc(beta))%>%
        mutate(stat_name="Pr(SSB>SSBlim)")

    # SSB>SSBbanとなる確率
    ssbban.table <- kobeII.data %>%
        dplyr::filter(year%in%year.ssbban,stat=="SSB") %>%
        group_by(HCR_name,beta,year) %>%
        summarise(ssb.over=round(100*mean(value>Bban))) %>%
        spread(key=year,value=ssb.over)%>%
        ungroup() %>%
        arrange(HCR_name,desc(beta))%>%
        mutate(stat_name="Pr(SSB>SSBban)")    

    # SSB>SSBmin(過去最低親魚量を上回る確率)
    ssb.min <- min(unlist(colSums(res_vpa$ssb)))
    ssbmin.table <- kobeII.data %>%
        dplyr::filter(year%in%year.ssbmin,stat=="SSB") %>%
        group_by(HCR_name,beta,year) %>%
        summarise(ssb.over=round(100*mean(value>ssb.min))) %>%
        spread(key=year,value=ssb.over)%>%
        ungroup() %>%
        arrange(HCR_name,desc(beta))%>%
        mutate(stat_name="Pr(SSB>SSBmin)")

    # SSB>SSBmax(過去最低親魚量を上回る確率)
    ssb.max <- max(unlist(colSums(res_vpa$ssb)))
    ssbmax.table <- kobeII.data %>%
        dplyr::filter(year%in%year.ssbmax,stat=="SSB") %>%
        group_by(HCR_name,beta,year) %>%
        summarise(ssb.over=round(100*mean(value>ssb.max))) %>%
        spread(key=year,value=ssb.over)%>%
        ungroup() %>%
        arrange(HCR_name,desc(beta))%>%
        mutate(stat_name="Pr(SSB>SSBmax)")    

    # オプション: Catch AAV mean 
    calc.aav <- function(x)sum(abs(diff(x)))/sum(x[-1])
    catch.aav.table <- kobeII.data %>%
        dplyr::filter(year%in%year.aav,stat=="catch") %>%
        group_by(HCR_name,beta,sim) %>%
        dplyr::summarise(catch.aav=(calc.aav(value))) %>%
        group_by(HCR_name,beta) %>%
        summarise(catch.aav.mean=mean(catch.aav)) %>%
        arrange(HCR_name,desc(beta))%>%
        mutate(stat_name="catch.csv (recent 5 year)")

    res_list <- list(catch.mean   = catch.mean,
                     ssb.mean         = ssb.mean,
                     ssb.lower10percent            = ssb.ci10,
                     ssb.upper90percent            = ssb.ci90,
                     prob.over.ssbtarget  = ssbtarget.table,
                     prob.over.ssblimit   = ssblimit.table,
                     prob.over.ssbban     = ssbban.table,                     
                     prob.over.ssbmin     = ssbmin.table,
                     prob.over.ssbmax     = ssbmax.table,                     
                     catch.aav       = catch.aav.table)    
    return(res_list)
                
}


HCR.simulation <- function(finput,HCRtable,year.lag=year.lag){
    
    tb <- NULL
    
    for(i in 1:nrow(HCRtable)){
        HCR_base <- HCRtable[i,]
        finput$multi <- HCR_base$Fmsy
        finput$HCR <- list(Blim=HCR_base$Blimit,Bban=HCR_base$Bban,
                           beta=HCR_base$beta,year.lag=year.lag)
        finput$is.plot <- FALSE
        finput$silent <- TRUE
        fres_base <- do.call(future.vpa,finput) # デフォルトルールの結果→図示などに使う
        tmp <- convert_future_table(fres_base,label=HCRtable$HCR_name[i]) %>%
            rename(HCR_name=label) 
        tmp$beta <- HCR_base$beta
        tb <- bind_rows(tb,tmp)
    }
    tb <- tb %>% mutate(HCR_name=str_c("beta",beta)) %>%
        mutate(scenario=HCR_name)
    return(tb)
}

#' kobeII matrixの簡易版（Btarget, Blimitは決め打ちでβのみ変える)
#'
#' @encoding UTF-8
#' @export
#'
#' 

beta.simulation <- function(finput,beta_vector,year.lag=0){
    
    tb <- NULL
    
    for(i in 1:length(beta_vector)){
        finput$HCR$beta <- beta_vector[i]
        finput$is.plot <- FALSE
        finput$silent <- TRUE
        fres_base <- do.call(future.vpa,finput) # デフォルトルールの結果→図示などに使う
        tmp <- convert_future_table(fres_base,label=beta_vector[i]) %>%
            rename(HCR_name=label)  %>% mutate(beta=beta_vector[i])
        tb <- bind_rows(tb,tmp)
    }
    return(tb)
}


get.stat4 <- function(fout,Brefs,
                      refyear=c(2019:2023,2028,2038)){
    col.target <- ifelse(fout$input$N==0,1,-1)
    years <- as.numeric(rownames(fout$vwcaa))

    if(is.null(refyear)){
        refyear <- c(seq(from=min(years),to=min(years)+5),
                           c(min(years)+seq(from=10,to=20,by=5)))
    }

    catch.mean <- rowMeans(fout$vwcaa[years%in%refyear,col.target])
    names(catch.mean) <- str_c("Catch",names(catch.mean)) 
    catch.mean <- as_tibble(t(catch.mean))
    
    Btarget.prob <- rowMeans(fout$vssb[years%in%refyear,col.target]>Brefs$Btarget) %>%
        t() %>% as_tibble() 
    names(Btarget.prob) <- str_c("Btarget_prob",names(Btarget.prob))

#    Blow.prob <- rowMeans(fout$vssb[years%in%refyear,col.target]>Brefs$Blow) %>%
#        t() %>% as_tibble() 
#    names(Blow.prob) <- str_c("Blow_prob",names(Blow.prob))

    Blimit.prob <- rowMeans(fout$vssb[years%in%refyear,col.target]<Brefs$Blimit) %>%
        t() %>% as_tibble() 
    names(Blimit.prob) <- str_c("Blimit_prob",names(Blimit.prob))

    Bban.prob <- rowMeans(fout$vssb[years%in%refyear,col.target]<Brefs$Bban) %>%
        t() %>% as_tibble() 
    names(Bban.prob) <- str_c("Bban_prob",names(Bban.prob))             

    return(bind_cols(catch.mean,Btarget.prob,Blimit.prob,Bban.prob))
}

#' Kobe plotを書く
#'
#' @param vpares VPAの結果のオブジェクト
#' @param refs_base est.MSYから得られる管理基準値の表
#' @encoding UTF-8
#'
#' @export
#' 

plot_kobe_gg <- plot_kobe <- function(vpares,refs_base,roll_mean=1,
                         category=4,# 削除予定オプション
                         Btarget=c("Btarget0"),
                         Blimit=c("Blimit0"),
                         Blow=c("Blow0"), # 削除予定オプション
                         Bban=c("Bban0"), 
                         write.vline=TRUE,
                         ylab.type="U", # or "U"
                         labeling.year=NULL,
                         RP.label=c("目標管理基準値","限界管理基準値","禁漁水準"),
                         refs.color=c("#00533E","#edb918","#C73C2E"),                         
                         Fratio=NULL, 
                         yscale=1.2,xscale=1.2,
                         HCR.label.position=c(1,1),  
                         beta=NULL,
                         plot.year="all"){

    target.RP <- derive_RP_value(refs_base,Btarget)
    limit.RP <- derive_RP_value(refs_base,Blimit)
    low.RP <- derive_RP_value(refs_base,Blow) 
    ban.RP <- derive_RP_value(refs_base,Bban)

    low.ratio <- low.RP$SSB/target.RP$SSB
    limit.ratio <- limit.RP$SSB/target.RP$SSB
    ban.ratio <- ban.RP$SSB/target.RP$SSB

    vpa_tb <- convert_vpa_tibble(vpares)
    UBdata <- vpa_tb %>% dplyr::filter(stat=="U" | stat=="SSB") %>%
        spread(key=stat,value=value) %>%
        mutate(Uratio=RcppRoll::roll_mean(U/target.RP$U,n=roll_mean,fill=NA,align="right"),
               Bratio=RcppRoll::roll_mean(SSB/target.RP$SSB,n=roll_mean,fill=NA,align="right")) %>%
        arrange(year)
    if(ylab.type=="F") UBdata <- UBdata %>% mutate(Uratio=Fratio)
    
    if(is.null(labeling.year)){
        years <- unique(UBdata$year)
        labeling.year <- c(years[years%%5==0],max(years))
    }

    UBdata <- UBdata %>%
        mutate(year.label=ifelse(year%in%labeling.year,year,""))
    
    if (plot.year[1]!="all") {
      UBdata <- UBdata %>% filter(year %in% plot.year)
    }

    max.B <- max(c(UBdata$Bratio,xscale),na.rm=T)
    max.U <- max(c(UBdata$Uratio,yscale),na.rm=T)

    g4 <- ggplot(data=UBdata) +theme(legend.position="none")+
        geom_polygon(data=tibble(x=c(-1,1,1,-1),
                                 y=c(-1,-1,1,1)),
                     aes(x=x,y=y),fill="khaki1")+
        geom_polygon(data=tibble(x=c(1,20,20,1),
                                 y=c(-1,-1,1,1)),
                     aes(x=x,y=y),fill="olivedrab2")+
        geom_polygon(data=tibble(x=c(1,20,20,1),
                                 y=c(1,1,20,20)),
                     aes(x=x,y=y),fill="khaki1")+
        geom_polygon(data=tibble(x=c(-1,1,1,-1),
                                 y=c(1,1,20,20)),
                     aes(x=x,y=y),fill="indianred1") +
        geom_polygon(data=tibble(x=c(-1,1,1,-1),
                                 y=c(-1,-1,1,1)),aes(x=x,y=y),fill="khaki1")

    if(write.vline){
        g4 <- g4 + geom_vline(xintercept=c(1,limit.ratio,ban.ratio),color=refs.color,lty="41",lwd=0.7)+
         ggrepel::geom_label_repel(data=tibble(x=c(1,limit.ratio,ban.ratio),
                                          y=max.U*0.85,
                                          label=RP.label),
                              aes(x=x,y=y,label=label),
                              direction="x",nudge_y=max.U*0.9,size=11*0.282)
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
        g4 <- g4+stat_function(fun = h,lwd=1.5,color=1,n=1000)+
            annotate("text",x=x.pos,y=y.pos,
                     label=str_c("漁獲管理規則\n(β=",beta,")"))
    }
   
    g4 <- g4 +
        geom_path(mapping=aes(x=Bratio,y=Uratio)) +        
        geom_point(mapping=aes(x=Bratio,y=Uratio),shape=21,fill="white") +
        coord_cartesian(xlim=c(0,max.B*1.1),ylim=c(0,max.U*1.15),expand=0) +
        ylab("漁獲率の比 (U/Umsy)") + xlab("親魚量の比 (SB/SBmsy)")  +
        ggrepel::geom_text_repel(#data=dplyr::filter(UBdata,year%in%labeling.year),
                         aes(x=Bratio,y=Uratio,label=year.label),
                         size=4,box.padding=0.5,segment.color="gray")

    if(ylab.type=="F"){
        g4 <- g4 + ylab("漁獲圧の比 (F/Fmsy)")        
    }
    
    g4 <- g4 + theme_SH()
    
    return(g4) 
}

#' 将来予測の複数の結果をggplotで重ね書きする
#'
#' @param vpares VPAの結果のオブジェクト
#' @param future.list 将来予測の結果をリストで並べたもの
#' @param n_example 個々のシミュレーションの例を示す数
#' @param width_example 個々のシミュレーションをプロットする場合の線の太さ (default=0.7)
#' @param future.replicate どのreplicateを選ぶかを選択する。この場合n_exampleによる指定は無効になる
#' @encoding UTF-8
#' @export

plot_futures <- function(vpares,
                         future.list=NULL,
                         future.name=names(future.list),
                         future_tibble=NULL,
                         CI_range=c(0.1,0.9),
                         maxyear=NULL,
                         font.size=18,
                         ncol=3,
                         is.plot.CIrange=TRUE,
                         what.plot=c("Recruitment","SSB","biomass","catch","Fsakugen","Fsakugen_ratio"),
                         biomass.unit=1,
                         RP_name=c("Btarget","Blimit","Bban"),
                         Btarget=0,Blimit=0,Bban=0,#Blow=0,
                         MSY=0,
                         exclude.japanese.font=FALSE, # english version
                         n_example=3, # number of examples
                         example_width=0.7, # line width of examples
                         future.replicate=NULL, 
                         seed=1 # seed for selecting the above example
                         ){

    col.SBtarget <- "#00533E"
    col.SBlim <- "#edb918"
    col.SBban <- "#C73C2E"
    col.MSY <- "black"
    col.Ftarget <- "#714C99"
    col.betaFtarget <- "#505596"

    if(!isTRUE(exclude.japanese.font)){
        junit <- c("","十","百","千","万")[log10(biomass.unit)+1]
        #    require(tidyverse,quietly=TRUE)

        rename_list <- tibble(stat=c("Recruitment","SSB","biomass","catch","Fsakugen","Fsakugen_ratio","alpha"),
                              jstat=c(str_c("加入尾数"),
                                      str_c("親魚量 (",junit,"トン)"),
                                      str_c("資源量 (",junit,"トン)"),
                                      str_c("漁獲量 (",junit,"トン)"),
                                      "努力量の削減率",
                                      "Fcurrentに対する乗数",
                                      "alpha"))
    }
    else{
        junit <- c("","10","100","1000","10,000")[log10(biomass.unit)+1]
        #    require(tidyverse,quietly=TRUE)

        rename_list <- tibble(stat=c("Recruitment","SSB","biomass","catch","Fsakugen","Fsakugen_ratio","alpha"),
                              jstat=c(str_c("Recruits"),
                                      str_c("SB (",junit,"MT)"),
                                      str_c("Biomass (",junit,"MT)"),
                                      str_c("Catch (",junit,"MT)"),
                                      "Effort reduction",
                                      "multiplier to Fcurrent",
                                      "alpha"))        
        }

    rename_list <- rename_list %>% dplyr::filter(stat%in%what.plot)
    
    if(!is.null(future.list)){
        if(is.null(future.name)) future.name <- str_c("s",1:length(future.list))
        names(future.list) <- future.name
    }
    else{
        if(is.null(future.name)) future.name <- str_c("s",1:length(unique(future_tibble$HCR_name)))
    }

    if(is.null(future_tibble)) future_tibble <- purrr::map_dfr(future.list,convert_future_table,.id="scenario")

    future.table <-
        future_tibble %>%
        dplyr::filter(stat%in%rename_list$stat) %>%
        mutate(stat=factor(stat,levels=rename_list$stat))

    set.seed(seed)
    if(!is.null(future.replicate)){
        future.example <- future.table %>%
            dplyr::filter(sim%in%future.replicate)
    }
    else{
        future.example <- future.table %>%
            dplyr::filter(sim%in%sample(2:max(future.table$sim),n_example))
    }
    
    future.example <- future.example %>%
        mutate(stat = as.character(stat),
             value=ifelse((stat=="Fsakugen"|stat=="Fsakugen_ratio"),
                          value,value/biomass.unit)) %>%
      left_join(rename_list) %>%
      group_by(sim,scenario)
        

    if(is.null(maxyear)) maxyear <- max(future.table$year)

    min.age <- as.numeric(rownames(vpares$naa)[1])
    vpa_tb <- convert_vpa_tibble(vpares) %>%
        dplyr::filter(stat=="SSB"|stat=="biomass"|stat=="catch"|stat=="Recruitment") %>%
        mutate(scenario=type,year=as.numeric(year),
               stat=factor(stat,levels=rename_list$stat),
               mean=value,sim=0)
    # 将来と過去をつなげるためのダミーデータ
    tmp <- vpa_tb %>% group_by(stat) %>%
        summarise(value=tail(value[!is.na(value)],n=1,na.rm=T),year=tail(year[!is.na(value)],n=1,na.rm=T),sim=0) 
    future.dummy <- purrr::map_dfr(future.name,function(x) mutate(tmp,scenario=x))

    org.warn <- options()$warn
    options(warn=-1)
    future.table <-
        bind_rows(future.table,vpa_tb,future.dummy) %>%
        mutate(stat=factor(stat,levels=rename_list$stat)) %>%
        mutate(scenario=factor(scenario,levels=c("VPA",future.name))) %>%
        mutate(value=ifelse(stat%in%c("Fsakugen","Fsakugen_ratio","alpha"),value,value/biomass.unit))

    future.table.qt <- 
        future.table %>% group_by(scenario,year,stat) %>%
        summarise(low=quantile(value,CI_range[1],na.rm=T),
                  high=quantile(value,CI_range[2],na.rm=T),
                  median=median(value,na.rm=T),
                  mean=mean(value))

    # make dummy for y range
    dummy <- future.table %>% group_by(stat) %>% summarise(max=max(value)) %>%
        mutate(value=0,year=min(future.table$year,na.rm=T)) %>%
        select(-max)

    dummy2 <- future.table %>% group_by(stat) %>%
        summarise(max=max(quantile(value,CI_range[2],na.rm=T))) %>%
        mutate(value=max*1.1,
               year=min(future.table$year,na.rm=T)) %>%
        select(-max)

    future.table.qt <- left_join(future.table.qt,rename_list) %>%
        mutate(jstat=factor(jstat,levels=rename_list$jstat))


    dummy     <- left_join(dummy,rename_list) %>% dplyr::filter(!is.na(stat))
    dummy2    <- left_join(dummy2,rename_list) %>% dplyr::filter(!is.na(stat))
    ssb_table <- tibble(jstat = dplyr::filter(rename_list, stat == "SSB") %>%
                          dplyr::pull(jstat),
                        value = c(Btarget, Blimit, Bban) / biomass.unit,
                        RP_name = RP_name)
    if(!is.null(MSY)){
        ssb_table <- tibble(jstat=dplyr::filter(rename_list, stat == "catch") %>%
                                  dplyr::pull(jstat),
                              value=MSY/biomass.unit,
                              RP_name="MSY") %>% bind_rows(ssb_table)
    }
    
    options(warn=org.warn)
    
    g1 <- future.table.qt %>% dplyr::filter(!is.na(stat)) %>%
        ggplot()+
        geom_line(data=dplyr::filter(future.table.qt,!is.na(stat) & scenario=="VPA"),
                  mapping=aes(x=year,y=mean),lwd=1,color=1)# VPAのプロット                

    if(isTRUE(is.plot.CIrange)){
        g1 <- g1+
            geom_ribbon(data=dplyr::filter(future.table.qt,!is.na(stat) & scenario!="VPA" & year <= maxyear),
                        mapping=aes(x=year,ymin=low,ymax=high,fill=scenario),alpha=0.4)+
            geom_line(data=dplyr::filter(future.table.qt,!is.na(stat) & scenario!="VPA" & year <= maxyear),
                      mapping=aes(x=year,y=mean,color=scenario),lwd=1)
    }
#    else{
#        g1 <- g1+
#            geom_line(data=dplyr::filter(future.table.qt,!is.na(stat) & scenario=="VPA"),
#                      mapping=aes(x=year,y=mean,color=scenario),lwd=1)#+        
#    }
    
    g1 <- g1+
        geom_blank(data=dummy,mapping=aes(y=value,x=year))+
        geom_blank(data=dummy2,mapping=aes(y=value,x=year))+
        theme_bw(base_size=font.size) +
        #        coord_cartesian(expand=0)+
        scale_y_continuous(expand=expand_scale(mult=c(0,0.05)))+
        theme(legend.position="top",panel.grid = element_blank())+
        facet_wrap(~factor(jstat,levels=rename_list$jstat),scales="free_y",ncol=ncol)+        
        xlab("年")+ylab("")+ labs(fill = "",linetype="",color="")+
        xlim(min(future.table$year),maxyear)+
        geom_hline(data = ssb_table,
                   aes(yintercept = value, linetype = RP_name),
                   color = c(col.SBtarget, col.SBlim, col.SBban,col.MSY))


    if(n_example>0){
        if(n_example>1){
            g1 <- g1 + geom_line(data=dplyr::filter(future.example,year <= maxyear),
                                 mapping=aes(x=year,y=value,
                                             alpha=factor(sim),
                                             color=scenario),
                                 lwd=example_width) 
        }
        else{
            g1 <- g1 + geom_line(data=dplyr::filter(future.example,year <= maxyear),
                                 mapping=aes(x=year,y=value,
                                             color=scenario),
                                 lwd=example_width) 
        }
        g1 <- g1+scale_alpha_discrete(guide=FALSE)            
    }
    return(g1)
}

#' F currentをプロットする
#'
#' @param vpares VPAの結果のオブジェクト
#' @encoding UTF-8
#'
#' @export

plot_Fcurrent <- function(vpares,
                          year.range=NULL){

    if(is.null(year.range)) year.range <- min(as.numeric(colnames(vpares$naa))):max(as.numeric(colnames(vpares$naa)))
    vpares_tb <- convert_vpa_tibble(vpares)

    fc_at_age <- vpares_tb %>%
        dplyr::filter(stat=="fishing_mortality", year%in%year.range) %>%
        mutate(F=value,year=as.character(year)) %>%
        select(-stat,-sim,-type,-value)
    fc_at_age_current <- tibble(F=vpares$Fc.at.age,age=as.numeric(rownames(vpares$naa)),
                                year="currentF")
    fc_at_age <- bind_rows(fc_at_age,fc_at_age_current) %>%
        mutate(F_name=c("gray","tomato")[as.numeric(year=="currentF")+1]) %>%
        group_by(year)
    
    g <- fc_at_age %>% ggplot() +
        geom_line(aes(x=age,y=as.numeric(F),alpha=year,linetype=F_name,color=F_name),lwd=1.5) +
        #        geom_line(data=fc_at_age_current,mapping=aes(x=age,y=as.numeric(F)),color="tomato",lwd=1.5)+        
        #        geom_point(data=fc_at_age_current,mapping=aes(x=age,y=as.numeric(F)),color="tomato",size=2)+
        #        scale_color_gradient(low="gray",high="blue")+
        scale_colour_identity()+
        theme_bw()+theme(legend.position="none")+
        coord_cartesian(expand=0,ylim=c(0,max(fc_at_age$F)*1.1),xlim=range(fc_at_age$age)+c(-0.5,0.5))+
        labs(year="年",color="",labels=c(gray="",tomato="Current F"))+
        #        theme(#legend.position="bottom",
        #            panel.grid = element_blank())+
        xlab("年齢")+ylab("漁獲係数(F)")#+
#    scale_color_discrete(name="F type",breaks=c())
#    scale_colour_manual(
#        values = c(
#            col1  = "gray",
#            col2  = "tomato",
#            col3  = "blue3",
    #            col4  = "yellow3")    )
    return(g)
}


library(ggplot2)

#Setting parameter values.
#SBtarget <- 250
#SBban <- 0.1*SBtarget
#SBlim <- 0.4*SBtarget
#Ftarget <-1.5
#beta <- 0.8

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
                     biomass.unit=1,
                     beta=0.8,col.multi2currf="black",col.SBtarget="#00533E",
                     col.SBlim="#edb918",col.SBban="#C73C2E",col.Ftarget="black",
                     col.betaFtarget="gray",is.text = TRUE){
    
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
        stat_function(fun = h,lwd=2,color=col.multi2currf)
    g <- ggplct  + geom_vline(xintercept = SBtarget, size = 0.9, linetype = "41", color = col.SBtarget) +
        geom_vline(xintercept = SBlim, size = 0.9, linetype = "41", color = col.SBlim) +
        geom_vline(xintercept = SBban, size = 0.9, linetype = "41", color = col.SBban) +
        geom_hline(yintercept = Ftarget, size = 0.9, linetype = "43", color = col.Ftarget) +
        geom_hline(yintercept = beta*Ftarget, size = 0.7, linetype = "43", color = col.betaFtarget) +
        labs(title = "",subtitle = "", caption =  "", x = str_c("親魚量 (",junit,"トン)"),
             y = "努力量の乗数",color = "") +
        theme_bw(base_size=12)+
        theme(legend.position="none",panel.grid = element_blank())+
        stat_function(fun = h,lwd=1.5,color=col.multi2currf)
    if (is.text) {
      g <- g +
        annotate("text", label="目標水準", x=SBtarget, y=1.2*Ftarget) +
        annotate("text", label="限界水準", x=SBlim, y=1.1*Ftarget) +
        annotate("text", label="禁漁水準", x=SBban, y=1.2*Ftarget)+
        annotate("text", label="Ftarget", x=SBtarget/15, y=0.95*Ftarget)+
        annotate("text", label=str_c(beta,"Ftarget"), x=SBtarget/15, y=0.95*beta*Ftarget)
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

# test plot
#Fig_Fish_Manage_Rule(SBtarget,SBlim,SBban,Ftarget,col.multi2currf = "#093d86", col.SBtarget = "#00533E", col.SBlim = "#edb918",col.SBban = "#C73C2E",col.Ftarget = "#714C99", col.betaFtarget = "#505596")
# function;ruri-rio, sbtarget;moegi-iro, sblim;koki-ki; sbban;hi-iro, ftarget;sumire-iro, betaftarget;kikyou-iro

#' MSYを達成するときの\%SPRを計算する
#'
#' @param finput 将来予測インプット
#' @param Fvector Fのベクトル
#' @encoding UTF-8
#' @export
calc_perspr <- function(finput,
                        Fvector,
                        Fmax=10,
                        max.age=Inf,
                        target.col=NULL
                        ){
    res_vpa <- finput$res0
    # MSYにおける将来予測計算をやりなおし
    finput$outtype <- "FULL"
    fout.tmp <- do.call(future.vpa,finput)
    # 生物パラメータはその将来予測で使われているものを使う
    if(is.null(target.col)){
        waa.tmp           <- fout.tmp$waa[,dim(fout.tmp$waa)[[2]],1]
        waa.catch.tmp <- fout.tmp$waa.catch[,dim(fout.tmp$waa.catch)[[2]],1]    
        maa.tmp           <- fout.tmp$maa[,dim(fout.tmp$maa)[[2]],1]
        M.tmp                <- fout.tmp$M[,dim(fout.tmp$M)[[2]],1]
    }
    else{
        waa.tmp           <- fout.tmp$waa[,target.col,1]
        waa.catch.tmp <- fout.tmp$waa.catch[,target.col,1]    
        maa.tmp           <- fout.tmp$maa[,target.col,1]
        M.tmp               <- fout.tmp$M[,target.col,1]
    }

    # SPRを計算
    spr.current <- ref.F(res_vpa,Fcurrent=as.numeric(Fvector),
                         waa=waa.tmp,
                         waa.catch=waa.catch.tmp,pSPR=NULL,
                         maa=maa.tmp,M=M.tmp,rps.year=as.numeric(colnames(res_vpa$naa)),
                         F.range=c(seq(from=0,to=ceiling(max(res_vpa$Fc.at.age,na.rm=T)*Fmax),
                                       length=101),max(res_vpa$Fc.at.age,na.rm=T)),
                         plot=FALSE,max.age=max.age)$currentSPR$perSPR
    spr.current
}

#' kobeIItable から任意の表を指名して取り出す
#'
#' @param kobeII_table \code{make_kobeII_table}の出力
#' @param name \code{kobeII_table}の要素名
#'
#' @encoding UTF-8
pull_var_from_kobeII_table <- function(kobeII_table, name) {
  table <- kobeII.table[[name]]
  table %>%
    dplyr::arrange(desc(beta)) %>%
    dplyr::select(-HCR_name, -stat_name)
}

#' kobeIItableから取り出した表を整形
#'
#' - 報告書に不要な列を除去する
#' - 単位を千トンに変換
#' @param beta_table \code{pull_var_from_kobeII_table}で取得した表
#' @param divide_by 表の値をこの値で除する．トンを千トンにする場合には1000
#' @param round TRUEなら値を丸める．漁獲量は現状整数表示なのでデフォルトはTRUE
format_beta_table <- function(beta_table, divide_by = 1, round = TRUE) {
  beta   <- beta_table %>%
    dplyr::select(beta) %>%
    magrittr::set_colnames("\u03B2") # greek beta in unicode
  values <- beta_table %>%
    dplyr::select(-beta) / divide_by
  if (round == TRUE) return(cbind(beta, round(values)))
  cbind(beta, values)
}

#' 値の大小に応じて表の背景にグラデーションをつける
#' @param beta_table \code{format_beta_table}で整形したβの表
#' @param color 表の背景となる任意の色
colorize_table <- function(beta_table, color) {
  beta_table %>%
    formattable::formattable(list(formattable::area(col = -1) ~
                                    formattable::color_tile("white", color)))
}

#' 表を画像として保存
#'
#' @inheritParams \code{\link{formattable::as.htmlwidget}}
#' @inheritParams \code{\link{htmltools::html_print}}
#' @inheritParams \code{\link{webshot::webshot}}
#' @param table ファイルとして保存したい表
#' @examples
#' \dontrun{
#' your_table %>%
#'  export_formattable(file = "foo.png")
#' }
#' @export
export_formattable <- function(table, file, width = "100%", height = NULL,
                               background = "white", delay = 0.1) {
  widget <- formattable::as.htmlwidget(table, width = width, height = height)
  path   <- htmltools::html_print(widget, background = background, viewer = NULL)
  url    <- paste0("file:///", gsub("\\\\", "/", normalizePath(path)))
  webshot::webshot(url,
                   file = file,
                   selector = ".formattable_widget",
                   delay = delay)
}

#' kobeIItableから任意の表を取得し画像として保存
#'
#' @inheritParams \code{\link{pull_var_from_kobeII_table}}
#' @inheritParams \code{\link{format_beta_table}}
#' @inheritParams \code{\link{colorize_table}}
#' @inheritParams \code{\link{export_formattable}}
export_kobeII_table <- function(name, divide_by, color, fname, kobeII_table) {
  kobeII_table %>%
    pull_var_from_kobeII_table(name) %>%
    format_beta_table(divide_by = divide_by) %>%
    colorize_table(color) %>%
    export_formattable(fname)
}

#' β調整による管理効果を比較する表を画像として一括保存
#'
#' @inheritParams \code{\link{pull_var_from_kobeII_table}}
#' @param fname_ssb 「平均親魚量」の保存先ファイル名
#' @param fname_catch 「平均漁獲量」の保存先ファイル名
#' @param fname_ssb_above_target 「親魚量が目標管理基準値を上回る確率」の保存先ファイル名
#' @param fname_ssb_above_limit 「親魚量が限界管理基準値を上回る確率」の保存先ファイル名
#' @examples
#' \dontrun{
#' export_kobeII_tables(kobeII.table)
#' }
#' @export
export_kobeII_tables <- function(kobeII_table,
                                 fname_ssb = "tbl_ssb.png",
                                 fname_catch = "tbl_catch.png",
                                 fname_ssb_above_target = "tbl_ssb>target.png",
                                 fname_ssb_above_limit = "tbl_ssb>limit.png") {
  blue   <- "#96A9D8"
  green  <- "#B3CE94"
  yellow <- "#F1C040"

  purrr::pmap(list(name = c("ssb.mean", "catch.mean",
                            "prob.over.ssbtarget", "prob.over.ssblimit"),
                   divide_by = c(1000, 1000, 1, 1),
                   color = c(blue, green, yellow, yellow),
                   fname = c(fname_ssb, fname_catch,
                             fname_ssb_above_target, fname_ssb_above_limit)),
              .f = export_kobeII_table,
              kobeII_table = kobeII_table)
}

#' 会議用の図のフォーマット
#'
#' @export
#' 

theme_SH <- function(){
    theme_bw(base_size=12) +
    theme(panel.grid = element_blank(),
          axis.text.x=element_text(size=11,color="black"),
          axis.text.y=element_text(size=11,color="black"),
          axis.line.x=element_line(size= 0.3528),
          axis.line.y=element_line(size= 0.3528),
          legend.position="none")
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
#' 

ggsave_SH_large <- function(...){
    ggsave(width=150,height=120,dpi=600,units="mm",...)
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
#' SRdata <- get.SRdata(res_vpa)
#' resSRregime <- fit.SRregime(SRdata, SR="HS", method="L2", 
#'                             regime.year=c(1977,1989), regime.key=c(0,1,0),
#'                             regime.par = c("a","b","sd")[2:3])
#' g1 <- SRregime_plot(resSRregime, regime.name=c("Low","High"))
#' g1
#' }
#' @encoding UTF-8
#' @export
#' 

SRregime_plot <- function (SRregime_result,xscale=1000,xlabel="SSB",yscale=1,ylabel="R",
                           labeling.year = NULL, show.legend = TRUE, legend.title = "Regime",regime.name = NULL,
                           base_size = 16, add.info = TRUE) {
  pred_data = SRregime_result$pred %>% mutate(Category = "Pred")
  obs_data = select(SRregime_result$pred_to_obs, -Pred, -resid) %>% mutate(Category = "Obs")
  combined_data = full_join(pred_data, obs_data)
  if (is.null(labeling.year)) labeling.year <- c(min(obs_data$Year),obs_data$Year[obs_data$Year %% 5 == 0],max(obs_data$Year))
  combined_data = combined_data %>% 
    mutate(label=if_else(is.na(Year),as.numeric(NA),if_else(Year %in% labeling.year, Year, as.numeric(NA)))) %>%
    mutate(SSB = SSB/xscale, R = R/yscale)
  g1 = ggplot(combined_data, aes(x=SSB,y=R,label=label)) + 
    geom_path(data=dplyr::filter(combined_data, Category=="Pred"),aes(group=Regime,colour=Regime,linetype=Regime),size=2, show.legend = show.legend)+
    geom_point(data=dplyr::filter(combined_data, Category=="Obs"),aes(group=Regime,colour=Regime),size=3, show.legend = show.legend)+
    geom_path(data=dplyr::filter(combined_data, Category=="Obs"),colour="darkgray",size=1)+
    xlab(xlabel)+ylab(ylabel)+
    ggrepel::geom_label_repel()+
    theme_bw(base_size=base_size)+
    coord_cartesian(ylim=c(0,combined_data$R*1.05),expand=0)
  if (show.legend) {
    if (is.null(regime.name)) {
      regime.name = unique(combined_data$Regime)
    }
    g1 = g1 + scale_colour_hue(name=legend.title, labels = regime.name) +
      scale_linetype_discrete(name=legend.title, labels = regime.name)
  }
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
                       "), ", ", regime_par: ", paste0(SRregime_result$input$regime.par,collapse="&"),", ",
                       paste0(SRregime_result$input$regime.year,collapse="&"), 
                       ", ",paste0(SRregime_result$input$regime.key,collapse="->"),
                       ", AICc: ",round(SRregime_result$AICc,2))
         )
    
  }
  g1
}}