
col.SBtarget <- "#00533E"
col.SBlim <- "#edb918"
col.SBban <- "#C73C2E"
col.Ftarget <- "#714C99"
col.betaFtarget <- "#505596"

pt1 <- 0.3528


convert_df <- function(df,name){
    df %>%
        as_tibble %>%  
        mutate(age = as.numeric(rownames(df))) %>% 
        gather(key=year, value=value, -age, convert=TRUE) %>%
        group_by(year) %>%
#        summarise(value=sum(value)) %>%
        mutate(type="VPA",sim="s0",stat=name)    
}

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
    
    bind_rows(ssb,catch,biomass,alpha_value,Fsakugen,Fsakugen_ratio,Recruitment)
}
        
    
convert_vector <- function(vector,name){
    vector %>%
        as_tibble %>%  
        mutate(year = as.integer(names(vector))) %>% 
        mutate(type="VPA",sim="s0",stat=name,age=NA) 
} 

convert_vpa_tibble <- function(vpares){

    total.catch <- colSums(vpares$input$dat$caa*vpares$input$dat$waa,na.rm=T)
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

SRplot_gg <- function(SR_result,refs=NULL,xscale=1000,xlabel="千トン",yscale=1,ylabel="尾",
                      labeling.year=NULL,add.info=TRUE){
    require(tidyverse,quietly=TRUE)    
    require(ggrepel)
    
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


    g1 <- ggplot() +
        geom_line(data=dplyr::filter(alldata,type=="pred"),
                  aes(y=R,x=SSB),color="deepskyblue3",lwd=1.3) +
        geom_path(data=dplyr::filter(alldata,type=="obs"),
                  aes(y=R,x=SSB),color=1) +
        geom_point(data=dplyr::filter(alldata,type=="obs"),
                   aes(y=R,x=SSB),shape=21,fill="white") +
#        scale_shape_discrete(solid=T)+        
#        geom_label_repel(data=dplyr::filter(alldata,type=="obs" & (year%%10==0|year==year.max)),
#                         aes(y=R,x=SSB,label=year),
    #                         size=3,box.padding=3,segment.color="black") +
    #        geom_text_repel(aes(y=R,x=SSB,label=pickyear)) +
    geom_text_repel(data=dplyr::filter(alldata,type=="obs"),
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
        g1 <- g1+geom_vline(xintercept=c(refs$Bmsy,refs$Blim,refs$Bban),linetype=2)
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
    
plot_yield <- function(MSY_obj,refs_base,
                       refs.label=NULL, # label for reference point
                       refs.color=c("#00533E","#edb918","#C73C2E"),
                       AR_select=FALSE,xlim.scale=1.1,
                       biomass.unit=1,labeling=TRUE,lining=TRUE,
                       age.label.ratio=0.9, # 年齢のラベルを入れる位置（xの最大値からの割合)
                       family = "JP1",
                       ylim.scale=1.2,future=NULL,past=NULL,future.name=NULL){
    
    junit <- c("","十","百","千","万")[log10(biomass.unit)+1]
   
    if ("trace" %in% names(MSY_obj)) {
      trace.msy <- MSY_obj$trace
    } else {
      trace.msy <- MSY_obj
    }
        
    require(tidyverse,quietly=TRUE)
    require(ggrepel)    

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
            tmpdata <- bind_rows(tmpdata,
                tibble(
                year        =as.numeric(rownames(future[[j]]$vssb)),
                ssb.future  =apply(future[[j]]$vssb[,-1],1,mean)/biomass.unit,
                catch.future=apply(future[[j]]$vwcaa[,-1],1,mean)/biomass.unit,
                scenario=future.name[j]))
            }
        tmpdata <- tmpdata %>% group_by(scenario)
        g1 <- g1 +
            geom_path(data=tmpdata,
                      mapping=aes(x=ssb.future,y=catch.future,
                                  linetype=factor(scenario)),
                      lwd=1,col=col.SBtarget)+
            geom_point(data=tmpdata,
                       mapping=aes(x=ssb.future,y=catch.future,
                                   shape=factor(scenario)),color=col.SBtarget,size=3)


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
            geom_label_repel(data=refs_base,
                             aes(y=ymax*ylim.scale*0.85,
                                 x=SSB,label=refs.label),
                             direction="x",size=11*0.282,nudge_y=ymax*ylim.scale*0.9)  
    }

    if(isTRUE(labeling)){
        g1 <- g1 +
            geom_point(data=refs_base,
                        aes(y=Catch,x=SSB))+
            geom_label_repel(data=refs_base,
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

make_RP_table <- function(refs_base){
    require(formattable)
    require(tidyverse,quietly=TRUE)
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

derive_RP_value <- function(refs_base,RP_name){
#    refs_base %>% dplyr::filter(RP.definition%in%RP_name)
#    subset(refs_base,RP.definition%in%RP_name)
    refs_base[refs_base$RP.definition%in%RP_name,]    
}


calc_kobeII_matrix <- function(fres_base,
                              refs_base,
                              Btarget=c("Btarget0"),
                              Blimit=c("Blimit0"),
                              Blow=c("Blow0"),
                              Bban=c("Bban0"),
                              year.lag=0,
                              beta=seq(from=0.5,to=1,by=0.1)){
    require(tidyverse,quietly=TRUE)    
# HCRの候補を網羅的に設定
#    HCR_candidate1 <- expand.grid(
#        Btarget_name=refs_base$RP.definition[str_detect(refs_base$RP.definition,Btarget)],
#        Blow_name=refs_base$RP.definition[str_detect(refs_base$RP.definition,Blow)],    
#        Blimit_name=refs_base$RP.definition[str_detect(refs_base$RP.definition,Blimit)],
#        Bban_name=refs_base$RP.definition[str_detect(refs_base$RP.definition,Bban)],
    #        beta=beta)

    refs.unique <- unique(c(Btarget,Blimit,Blow,Bban))
    tmp <- !refs.unique%in%refs_base$RP.definition    
    if(sum(tmp)>0) stop(refs.unique[tmp]," does not appear in column of RP.definition\n")

    HCR_candidate1 <- expand.grid(
        Btarget_name=derive_RP_value(refs_base,Btarget)$RP.definition,
        Blow_name=derive_RP_value(refs_base,Blow)$RP.definition,    
        Blimit_name=derive_RP_value(refs_base,Blimit)$RP.definition,
        Bban_name=derive_RP_value(refs_base,Bban)$RP.definition,
        beta=beta)    

    HCR_candidate2 <- expand.grid(
        Btarget=derive_RP_value(refs_base,Btarget)$SSB,
        Blow=derive_RP_value(refs_base,Blow)$SSB,    
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

    kobeII_table <- left_join(kobeII_table,HCR_candidate)
    kobeII_table    
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
    tb <- tb %>% mutate(scenario=str_c(HCR_name,beta))
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

    Blow.prob <- rowMeans(fout$vssb[years%in%refyear,col.target]>Brefs$Blow) %>%
        t() %>% as_tibble() 
    names(Blow.prob) <- str_c("Blow_prob",names(Blow.prob))

    Blimit.prob <- rowMeans(fout$vssb[years%in%refyear,col.target]<Brefs$Blimit) %>%
        t() %>% as_tibble() 
    names(Blimit.prob) <- str_c("Blimit_prob",names(Blimit.prob))

    Bban.prob <- rowMeans(fout$vssb[years%in%refyear,col.target]<Brefs$Bban) %>%
        t() %>% as_tibble() 
    names(Bban.prob) <- str_c("Bban_prob",names(Bban.prob))             

    return(bind_cols(catch.mean,Btarget.prob,Blow.prob,Blimit.prob,Bban.prob))
}


plot_kobe_gg <- function(vpares,refs_base,roll_mean=1,
                         category=4,# 4区分か、6区分か
                         Btarget=c("Btarget0"),
                         Blimit=c("Blimit0"),
                         Blow=c("Blow0"),
                         Bban=c("Bban0"),write.vline=TRUE,
                         ylab.type="U", # or "U"
                         labeling.year=NULL,
                         Fratio=NULL, # ylab.type=="F"のとき
                         yscale=1.2,xscale=1.2,
                         HCR.label.position=c(1,1), # デフォルトはx軸方向が1, y軸方向が1の相対値です。様子を見ながら調整してください
                         refs.color=c("#00533E","#edb918","#C73C2E"),
                         beta=NULL){

   
    require(tidyverse,quietly=TRUE)
    require(ggrepel,quietly=TRUE)    

    target.RP <- derive_RP_value(refs_base,Btarget)
    limit.RP <- derive_RP_value(refs_base,Blimit)
    low.RP <- derive_RP_value(refs_base,Blow) 
    ban.RP <- derive_RP_value(refs_base,Bban)

    low.ratio <- low.RP$SSB/target.RP$SSB
    limit.ratio <- limit.RP$SSB/target.RP$SSB
    ban.ratio <- ban.RP$SSB/target.RP$SSB

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
    
    require(RcppRoll)
    vpa_tb <- convert_vpa_tibble(vpares)
    UBdata <- vpa_tb %>% dplyr::filter(stat=="U" | stat=="SSB") %>%
        spread(key=stat,value=value) %>%
        mutate(Uratio=roll_mean(U/target.RP$U,n=roll_mean,fill=NA,align="right"),
               Bratio=roll_mean(SSB/target.RP$SSB,n=roll_mean,fill=NA,align="right")) %>%
        arrange(year)
    if(ylab.type=="F") UBdata <- UBdata %>% mutate(Uratio=Fratio)
    
    if(is.null(labeling.year)){
        years <- unique(UBdata$year)
        labeling.year <- c(years[years%%5==0],max(years))
    }

    UBdata <- UBdata %>%
        mutate(year.label=ifelse(year%in%labeling.year,year,""))

    max.B <- max(c(UBdata$Bratio,xscale),na.rm=T)
    max.U <- max(c(UBdata$Uratio,yscale),na.rm=T)

    g6 <- ggplot(data=UBdata) + theme(legend.position="none")+
        geom_polygon(data=tibble(x=c(-1,low.ratio,low.ratio,-1),
                                 y=c(-1,-1,1,1)),
                     aes(x=x,y=y),fill="khaki1")+
        geom_polygon(data=tibble(x=c(low.ratio,10,10,low.ratio),
                                 y=c(-1,-1,1,1)),
                     aes(x=x,y=y),fill="olivedrab2")+
        geom_polygon(data=tibble(x=c(low.ratio,10,10,low.ratio),
                                 y=c(1,1,10,10)),
                     aes(x=x,y=y),fill="khaki1")+
        geom_polygon(data=tibble(x=c(-1,limit.ratio,limit.ratio,-1),
                                 y=c(1,1,10,10)),
                     aes(x=x,y=y),fill="indianred1") +
        geom_polygon(data=tibble(x=c(limit.ratio,low.ratio,low.ratio,limit.ratio),
                                 y=c(1,1,10,10)),aes(x=x,y=y),fill="tan1") +
        geom_polygon(data=tibble(x=c(-1,limit.ratio,limit.ratio,-1),
                                 y=c(-1,-1,1,1)),aes(x=x,y=y),fill="khaki2") +
        geom_polygon(data=tibble(x=c(limit.ratio,low.ratio,low.ratio,limit.ratio),
                                 y=c(-1,-1,1,1)),aes(x=x,y=y),fill="khaki1")+
        geom_vline(xintercept=c(1,ban.ratio),linetype=2)   

    g4 <- ggplot(data=UBdata) +theme(legend.position="none")+
        geom_polygon(data=tibble(x=c(-1,low.ratio,low.ratio,-1),
                                 y=c(-1,-1,1,1)),
                     aes(x=x,y=y),fill="khaki1")+
        geom_polygon(data=tibble(x=c(1,10,10,1),
                                 y=c(-1,-1,1,1)),
                     aes(x=x,y=y),fill="olivedrab2")+
        geom_polygon(data=tibble(x=c(1,10,10,1),
                                 y=c(1,1,10,10)),
                     aes(x=x,y=y),fill="khaki1")+
        geom_polygon(data=tibble(x=c(-1,1,1,-1),
                                 y=c(1,1,10,10)),
                     aes(x=x,y=y),fill="indianred1") +
        geom_polygon(data=tibble(x=c(-1,1,1,-1),
                                 y=c(-1,-1,1,1)),aes(x=x,y=y),fill="khaki1")

    if(write.vline){
     if(low.ratio<1){
        g6 <- g6 + geom_text(data=tibble(x=c(ban.ratio,limit.ratio,low.ratio,1),
                              y=rep(0.1,4),
                              label=c("Bban","Blimit","Blow","Btarget")),
                             aes(x=x,y=y,label=label))
        g4 <- g4 + geom_vline(xintercept=c(ban.ratio,limit.ratio,low.ratio,1),linetype=2)+
        geom_text(data=tibble(x=c(ban.ratio,limit.ratio,low.ratio,1),
                              y=rep(0.1,4),
                              label=c("Bban","Blimit","Blow","Btarget")),
                  aes(x=x,y=y,label=label))
     }else{
         
        g6 <- g6 + geom_text(data=tibble(x=c(ban.ratio,limit.ratio,1),
                                         y=max.U*c(1.05,1,1.05),
                                         label=c("禁漁水準","限界管理基準値","目標管理基準値")),
                             aes(x=x,y=y,label=label))
        g4 <- g4 + geom_vline(xintercept=c(1,limit.ratio,ban.ratio),color=refs.color,lty="41",lwd=0.7)+
#            geom_text(data=tibble(x=c(1,limit.ratio,ban.ratio),
#                                  y=max.U*c(1.05,1.1,1.05),
#                                  label=c("目標管理基準値","限界管理基準値","禁漁水準")),
    #                      aes(x=x,y=y,label=label),hjust=0)
             geom_label_repel(data=tibble(x=c(1,limit.ratio,ban.ratio),
                                          y=max.U*0.85,
                                          label=c("目標管理基準値","限界管理基準値","禁漁水準")),
                              aes(x=x,y=y,label=label),
                              direction="x",nudge_y=max.U*0.9,size=11*0.282)
    }}    

    if(!is.null(beta)){
        x.pos <- max.B*HCR.label.position[1]
        y.pos <- multi2currF(1.05)*HCR.label.position[2]
        g6 <- g6+stat_function(fun = h,lwd=1.5,color=1,n=1000)+
            annotate("text",x=x.pos,y=y.pos,            
                     label=str_c("漁獲管理規則\n(β=",beta,")"))            
        g4 <- g4+stat_function(fun = h,lwd=1.5,color=1,n=1000)+
            annotate("text",x=x.pos,y=y.pos,
                     label=str_c("漁獲管理規則\n(β=",beta,")"))
#        if(abs(HCR.label.position[1]-1)+abs(HCR.label.position[2]-1)>0.3){
#            label.line <- tibble(x=c(max.B,x.pos),
#                                 y=c(multi2currF(1.05),y.pos))
#            g6 <- g6 + geom_path(data=label.line,mapping=aes(x=x,y=y),color="gray")
#            g4 <- g4 + geom_path(data=label.line,mapping=aes(x=x,y=y),color="gray")
#        }
    }
   
    g6 <- g6 +
        geom_path(mapping=aes(x=Bratio,y=Uratio)) +
        geom_point(mapping=aes(x=Bratio,y=Uratio),shape=21,fill="white") +        
        coord_cartesian(xlim=c(0,max.B*1.1),ylim=c(0,max.U*1.15),expand=0) +
        ylab("漁獲率の比 (U/Umsy)") + xlab("親魚量の比 (SB/SBmsy)")  +
        geom_text_repel(#data=dplyr::filter(UBdata,year%in%labeling.year),
                         aes(x=Bratio,y=Uratio,label=year.label),
                         size=4,box.padding=0.5,segment.color="gray")

    g4 <- g4 +
        geom_path(mapping=aes(x=Bratio,y=Uratio)) +        
        geom_point(mapping=aes(x=Bratio,y=Uratio),shape=21,fill="white") +
        coord_cartesian(xlim=c(0,max.B*1.1),ylim=c(0,max.U*1.15),expand=0) +
        ylab("漁獲率の比 (U/Umsy)") + xlab("親魚量の比 (SB/SBmsy)")  +
        geom_text_repel(#data=dplyr::filter(UBdata,year%in%labeling.year),
                         aes(x=Bratio,y=Uratio,label=year.label),
                         size=4,box.padding=0.5,segment.color="gray")

    if(ylab.type=="F"){
        g6 <- g6 + ylab("漁獲圧の比 (F/Fmsy)")
        g4 <- g4 + ylab("漁獲圧の比 (F/Fmsy)")        
    }
    
    if(category==4) return(g4) else return(g6)
}

plot_futures <- function(vpares,
                         future.list=NULL,
                         future.name=names(future.list),
                         future_tibble=NULL,
                         CI_range=c(0.1,0.9),
                         maxyear=NULL,font.size=18,
                         ncol=3,
                         what.plot=c("Recruitment","SSB","biomass","catch","Fsakugen","Fsakugen_ratio"),
                         biomass.unit=1,RP_name=c("Btarget","Blimit","Bban"),
                         Btarget=0,Blimit=0,Bban=0,#Blow=0,
                         n_example=3, # number of examples
                         seed=1 # seed for selecting the above example
                         ){

    junit <- c("","十","百","千","万")[log10(biomass.unit)+1]
    require(tidyverse,quietly=TRUE)
    rename_list <- tibble(stat=c("Recruitment","SSB","biomass","catch","Fsakugen","Fsakugen_ratio","alpha"),
                          jstat=c(str_c("加入尾数"),
                              str_c("親魚量 (",junit,"トン)"),
                              str_c("資源量 (",junit,"トン)"),
                              str_c("漁獲量 (",junit,"トン)"),
                              "努力量の削減率",
                              "Fcurrentに対する乗数",
                              "alpha"))

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
    future.example <- future.table %>%
      dplyr::filter(sim%in%sample(2:max(future.table$sim),n_example)) %>%
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
    tmp <- vpa_tb %>% group_by(stat) %>%
        summarise(value=tail(value[!is.na(value)],n=1,na.rm=T),year=tail(year[!is.na(value)],n=1,na.rm=T),sim=0) 
    future.dummy <- purrr::map_dfr(future.name,function(x) mutate(tmp,scenario=x))

    org.warn <- options()$warn
    options(warn=-1)
    future.table <- bind_rows(future.table,vpa_tb,future.dummy) %>%
        mutate(stat=factor(stat,levels=rename_list$stat)) %>%
        mutate(scenario=factor(scenario,levels=c("VPA",future.name))) %>%
        mutate(value=ifelse(stat%in%c("Fsakugen","Fsakugen_ratio","alpha"),value,value/biomass.unit))

    future.table.qt <- future.table %>% group_by(scenario,year,stat) %>%
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
    
    options(warn=org.warn)
    
    g1 <- future.table.qt %>% dplyr::filter(!is.na(stat)) %>%
        ggplot() +
        geom_ribbon(aes(x=year,ymin=low,ymax=high,fill=scenario),alpha=0.4)+        
        geom_line(aes(x=year,y=mean,color=scenario),lwd=1)+
        geom_line(aes(x=year,y=mean,color=scenario),linetype=2,lwd=1)+
        geom_blank(data=dummy,mapping=aes(y=value,x=year))+
        geom_blank(data=dummy2,mapping=aes(y=value,x=year))+
        theme_bw(base_size=font.size) +
        coord_cartesian(expand=0)+
        theme(legend.position="top",panel.grid = element_blank())+
        facet_wrap(~factor(jstat,levels=rename_list$jstat),scales="free_y",ncol=ncol)+        
        xlab("年")+ylab("")+ labs(fill = "",linetype="",color="")+
        xlim(min(future.table$year),maxyear)+
        geom_hline(data = ssb_table,
                   aes(yintercept = value, linetype = RP_name),
                   color = c(col.SBtarget, col.SBlim, col.SBban))


    if(n_example>0){
        g1 <- g1 + geom_line(data=future.example,
                             mapping=aes(x=year,y=value,
                                         alpha=factor(sim),
                                         color=scenario)) +
            scale_alpha_discrete(guide=FALSE)
            
    }
    return(g1)
}

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

plot_HCR <- function(SBtarget,SBlim,SBban,Ftarget,
                     biomass.unit=1,
                     beta=0.8,col.multi2currf="black",col.SBtarget="#00533E",
                     col.SBlim="#edb918",col.SBban="#C73C2E",col.Ftarget="black",
                     col.betaFtarget="gray"){
    
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
        annotate("text", label="目標水準", x=SBtarget, y=1.2*Ftarget) +
        annotate("text", label="限界水準", x=SBlim, y=1.1*Ftarget) +
        annotate("text", label="禁漁水準", x=SBban, y=1.2*Ftarget)+
        annotate("text", label="Ftarget", x=SBtarget/15, y=0.95*Ftarget)+
        annotate("text", label=str_c(beta,"Ftarget"), x=SBtarget/15, y=0.95*beta*Ftarget)+
        theme_bw(base_size=12)+
        theme(legend.position="none",panel.grid = element_blank())+
        stat_function(fun = h,lwd=1.5,color=col.multi2currf)        

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


# MSYを達成するときの%SPRを計算する
calc_MSY_spr <- function(MSYres,Fmax=10,max.age=Inf){
    dres <- MSYres$input$msy$res0
    # MSYにおける将来予測計算をやりなおし
    MSYres$input.list$msy$outtype <- "FULL"
    fout.msy <- do.call(future.vpa,MSYres$input.list$msy)
    # 生物パラメータはその将来予測で使われているものを使う
    waa.msy <- fout.msy$waa[,dim(fout.msy$waa)[[2]],1]
    maa.msy <- fout.msy$maa[,dim(fout.msy$maa)[[2]],1]
    M.msy <- fout.msy$M[,dim(fout.msy$M)[[2]],1]
    # F.msyの定義
    F.msy <- MSYres$input$msy$multi*MSYres$input$msy$res0$Fc.at.age

    # PPRを計算
    dres$Fc.at.age <- F.msy
    spr.msy <- ref.F(dres,waa=waa.msy,maa=maa.msy,M=M.msy,rps.year=as.numeric(colnames(dres$naa)),
                     F.range=c(seq(from=0,to=ceiling(max(dres$Fc.at.age,na.rm=T)*Fmax),
                                   length=101),max(dres$Fc.at.age,na.rm=T)),plot=FALSE,max.age=max.age)$ypr.spr
    target.SPR <- spr.msy[spr.msy$Frange2Fcurrent==1,]$spr[1]
    target.SPR
}
