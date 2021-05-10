#' make.plots: plot output from APAM
#'
#' Plots output from APAM.
#'
#' @importFrom ggplot2 ggplot aes geom_line geom_ribbon labs theme theme_bw element_text facet_wrap
#' @importFrom ggplot2 scale_x_continuous scale_y_continuous geom_point guides
#' @importFrom rlang .data
#'
#' @param mfits object from \code{\link{make.fit}}
#'
#' @return a list containing all ggplots
#'   \describe{
#'     \item{\code{$indexfit1}}{observed vs predicted survey indices at age for ages 1-7}
#'     \item{\code{$indexfit2}}{observed vs predicted survey indices at age for ages 8-15}
#'     \item{\code{$resid_index}}{survey Z residuals}
#'     \item{\code{$pop_process}}{estimated population processes (recruitment, average fishing mortality rates,
#'     spawning stock biomass and population abundance)}
#'     \item{\code{$index_bplot}}{bubble plot of survey index residuals}
#'     \item{\code{$crl_fit}}{observed vs predicted crl fits}
#'     \item{\code{$resid_crl}}{crl Z residuals}
#'     \item{\code{$crl_bplot}}{bubble plot of crl residuals}
#'     \item{\code{$landings}}{estimated model landings with upper/lower bounds}
#'     \item{\code{$landratio}}{estimated landings ratios}
#'     \item{\code{$F_at_age}}{model predicted fishing mortality rates at age}
#'     \item{\code{$log_F_dev}}{model predicted log fishing mortality rate deviations}
#'     \item{\code{$process_errors}}{model predicted process errors}
#'     \item{\code{$pe_bplot}}{bubble plot of model predicted process errors}
#'     \item{\code{$catchability}}{estimated survey cathcability}
#'     }
#'
#' @examples
#' \dontrun{
#' plots<- make.plots(mfits)
#' }
#' @export
make.plots = function(mfits){

  #define variables
  years=1960:2017
  endyear <- years[mfits$tmb.data$Y]
  nages <- length(1:mfits$tmb.data$A)

  #plot pop processes
  value.names <- names(mfits$sdrep$value)
  pop_process = c("log_Rec","log_biomass","log_ssb","log_aveF_914","log_N")

  Popest = matrix(NA, nrow = mfits$tmb.data$Y, ncol = (length(pop_process)*3))
  colnames(Popest) = c("recruit","recruit.L95","recruit.U95","biomass","biomass.L95","biomass.U95",
                       "ssb","ssb.L95","ssb.U95","aveF_914","aveF_914.L95","aveF_914.U95",
                       "abund","abund.L95","abund.U95")

  j=1
  for(i in 1:length(pop_process)){

    ind <- value.names == pop_process[i]

    if(pop_process[i]!="log_N"){
      Popest[,j]=exp(mfits$sdrep$value[ind])
      Popest[,j+1]=exp(mfits$sdrep$value[ind] - stats::qnorm(0.975)*mfits$sdrep$sd[ind])
      Popest[,j+2]=exp(mfits$sdrep$value[ind] + stats::qnorm(0.975)*mfits$sdrep$sd[ind])
    }

    if(pop_process[i]=="log_N"){
      N = matrix(exp(mfits$sdrep$value[ind]), nrow = mfits$tmb.data$Y, ncol=  mfits$tmb.data$A)
      NU95 = matrix(c(exp(mfits$sdrep$value[ind] - stats::qnorm(0.975)*mfits$sdrep$sd[ind])), nrow = mfits$tmb.data$Y, ncol= mfits$tmb.data$A)
      NL95 = matrix(c(exp(mfits$sdrep$value[ind] + stats::qnorm(0.975)*mfits$sdrep$sd[ind])) , nrow = mfits$tmb.data$Y, ncol= mfits$tmb.data$A)
      Popest[,j] = rowSums(N);
      Popest[,j+1] = rowSums(NU95)
      Popest[,j+2] = rowSums(NL95)}

    j= j + 3
  }

  Popest = data.frame(Popest)
  Popest$year=1960:endyear

  #population processes
  Abund = Popest %>% ggplot(aes(x = .data$year, y = .data$abund))+
    geom_line(col = "black",size=0.75)+
    geom_ribbon(aes(ymin = .data$abund.L95, ymax = .data$abund.U95), fill = "black", alpha = 0.10)+
    scale_x_continuous(breaks = seq(1960,2017,8))+
    labs(y="Abund (000's)", x = "")+
    theme_bw()+
    theme(text = element_text(size=15))

  Recruit = Popest %>% ggplot(aes(x = .data$year, y = .data$recruit))+
    geom_line(col = "black",size=0.75)+
    geom_ribbon(aes(ymin = .data$recruit.L95, ymax = .data$recruit.U95), fill = "black", alpha = 0.10)+
    scale_x_continuous(breaks = seq(1960,2018,8))+
    labs(y="Rec (000's tons)", x = "Year")+
    theme_bw()+
    theme(text = element_text(size=15))

  AveF = Popest %>% ggplot(aes(x = .data$year, y = .data$aveF_914))+
    geom_line(col = "black",size=0.75)+
    geom_ribbon(aes(ymin = .data$aveF_914.L95, ymax = .data$aveF_914.U95), fill = "black", alpha = 0.10)+
    scale_x_continuous(breaks = seq(1960,2018,8))+
    labs(y="Ave F (9-14) ", x = "Year")+
    theme_bw()+
    theme(text = element_text(size=15))

  SSB = Popest %>% ggplot(aes(x = .data$year, y = .data$ssb))+
    geom_line(col = "black",size=0.75)+
    geom_ribbon(aes(ymin = .data$ssb.L95, ymax = .data$ssb.U95), fill = "black", alpha = 0.10)+
    scale_x_continuous(breaks = seq(1960,2018,8))+
    labs(y="SSB (000's tons)",x="")+
    theme_bw()+
    theme(text = element_text(size=15))


  pp <- cowplot::plot_grid(Abund,SSB,Recruit,AveF)

  #to plot predicted indices
  indexfit = index %>% dplyr::filter(.data$Year<=endyear)
  indexfit = indexfit %>% dplyr::mutate(YC = .data$Year-.data$Age,
                                 Elog_index = mfits$rep$Elog_index,
                                 Eindex = mfits$rep$Eindex,
                                 std.Resid =  mfits$rep$std_resid_index,
                                 resid.index = mfits$rep$resid_index,
                                 log_index = log(indexfit$index),
                                 Zresid = mfits$rep$index_Zresid,
                                 colr ='deepskyblue')

  indexfit$colr[indexfit$Zresid>0]<-'firebrick2'

  index_fit1 <- indexfit %>% dplyr::filter(.data$Age<=7) %>%
    dplyr::select(.data$Age, .data$survey, .data$index, .data$Eindex, .data$Year) %>%
    reshape2::melt(id = c("Age","Year","survey")) %>%
    ggplot( aes(x = .data$Year, y = .data$value, col = .data$variable) )+
    geom_line(size = 1)+
    labs(y="Abundance (millions)")+
    ggplot2::scale_color_manual(values=c("black", "red"), labels=c("Observed","Predicted"), name = "", position = "top")+
    facet_wrap(.data$survey~.data$Age,nrow=3,scales = "free")+
    scale_x_continuous(breaks = seq(1960,2018,5))+
    theme_bw()+
    theme(legend.position = 'none',text = element_text(size=15),
          axis.text.x = element_text(angle = 90, hjust = 1), plot.margin = ggplot2::margin(10, 10, 10, 12))


  index_fit2  <-indexfit %>% dplyr::filter(.data$Age>7) %>%
    dplyr::select(.data$Age,.data$survey, .data$index, .data$Eindex, .data$Year) %>%
    reshape2::melt(id = c("Age","Year","survey")) %>%
    ggplot(aes(x = .data$Year, y = .data$value, col = .data$variable) )+
    geom_line(size = 1)+
    labs(y="Abundance (millions)")+
    ggplot2::scale_color_manual(values=c("black", "red"),labels=c("Observed","Predicted"), name = "", position = "top")+
    facet_wrap(.data$survey~.data$Age,nrow=3, scales = "free")+
    scale_x_continuous(breaks = seq(1960,2018,5))+
    theme_bw()+
    theme(legend.position = 'none',text = element_text(size=15),axis.text.x = element_text(angle = 90, hjust = 1),
          plot.margin = ggplot2::margin(10, 10, 10, 12))

  #plot adjusted survey residuals
  xx<-NULL
  yy<-NULL
  zz<-NULL
  j=1
  for(j in 1:3){

    SS<- unique(indexfit$survey)[j]

    pp1<-  indexfit %>% dplyr::filter(.data$survey==SS) %>%
      ggplot(aes(x=.data$Year, y = .data$Zresid))+
      ggplot2::geom_text(aes(label = .data$Age), size = 1.2)+
      ggplot2::geom_hline(yintercept = 0)+
      geom_line(data = indexfit  %>% dplyr::filter(.data$survey==SS) %>%
                  dplyr::group_by(.data$Year) %>%
                  dplyr::summarise(yy=mean(.data$Zresid),.groups="drop"),
                aes(x=.data$Year, y=.data$yy), color="red", size=1)+
      scale_x_continuous(breaks = seq(1960,2018,12))+
      labs(title=paste(SS,sep=" "),y = "")+
      theme_bw()

    pp2<-indexfit %>% dplyr::filter(.data$survey==SS) %>%
      ggplot(aes(x=.data$YC, y = .data$Zresid))+
      geom_point(shape=3)+
      ggplot2::geom_hline(yintercept = 0)+
      geom_line(data = indexfit  %>% dplyr::filter(.data$survey==SS) %>%
                  dplyr::group_by(.data$YC) %>%
                  dplyr::summarise(yy=mean(.data$Zresid),.groups="drop"),
                aes(x=.data$YC, y=.data$yy), color="red", size=1)+
      scale_x_continuous(breaks = seq(1960,2018,12))+
      labs(y="")+
      theme_bw()

    pp3<-indexfit %>% dplyr::filter(.data$survey==SS) %>%
      ggplot(aes(x=.data$Age, y = .data$Zresid))+
      geom_point(shape=3)+
      ggplot2::geom_hline(yintercept = 0)+
      geom_line(data = indexfit  %>% dplyr::filter(.data$survey==SS) %>%
                  dplyr::group_by(.data$Age) %>%
                  dplyr::summarise(yy=mean(.data$Zresid),.groups="drop"),
                aes(x=.data$Age, y=.data$yy), color="black", size=1)+
      labs(x = "Age", y="")+
      theme_bw()

    pp4<-indexfit %>% dplyr::filter(.data$survey==SS) %>%
      ggplot(aes(x=.data$Eindex, y = .data$Zresid))+
      geom_point(shape=3)+
      ggplot2::geom_hline(yintercept = 0)+
      labs(y="")+
      theme_bw()

    if(j==1){xx<- cowplot::plot_grid(pp1,pp2,pp3,pp4, ncol = 1)}
    if(j==2){yy<- cowplot::plot_grid(pp1,pp2,pp3,pp4, ncol = 1)}
    # if(j==3){zz<- (pp1/pp2/pp3/pp4)}
  }
  resid_index<- cowplot::plot_grid(xx,yy)

  ###index bubble plots
  index_bplot <-  indexfit %>% ggplot( aes(x = .data$Year, y = .data$Age,size = abs(.data$Zresid), color = .data$colr))+
    geom_point(shape = 20, alpha = 0.7)+
    ggplot2::scale_size(range = c(0,15))+
    ggplot2::facet_grid(cols = ggplot2::vars(.data$survey))+
    scale_x_continuous(breaks = seq(1960,2018,5))+
    scale_y_continuous(breaks = seq(1,15,1))+
    ggplot2::scale_color_manual(values = c(  "deepskyblue", "firebrick2"))+
    guides(size = FALSE,color = FALSE)+
    labs(title="red bubbles are positive, blue are negative")+
    theme_bw()

  #continuation ratio logit plots
  temp<-as.data.frame(mfits$tmb.data$crl,row.names = 1960:endyear-1) %>% stats::setNames(5:14)
  CRL.vec<- reshape2::melt(as.matrix(temp))

  temp<-as.data.frame(mfits$rep$pred_crl,row.names = 1960:endyear-1) %>% stats::setNames(5:14)
  CRL.vec$pred <-reshape2::melt(as.matrix(temp))$value

  temp<-as.data.frame(mfits$rep$std_resid_crl,row.names = 1960:endyear-1) %>% stats::setNames(5:14)
  CRL.vec$std_resid_crl<-reshape2::melt(as.matrix(temp))$value

  temp<-as.data.frame(mfits$rep$crl_Zresid,row.names = 1960:endyear-1) %>% stats::setNames(5:14)
  CRL.vec$Zresid<-reshape2::melt(as.matrix(temp))$value

  CRL.vec$colr = 'deepskyblue'
  CRL.vec$colr[CRL.vec$Zresid>0]='firebrick2'
  CRL.vec$cohort = CRL.vec$Var1-CRL.vec$Var2

  crl_fit <- CRL.vec %>% ggplot(aes(x = .data$Var1, y = .data$value) )+
    geom_line(size = 1,aes(col = "red"))+
    geom_line(aes(x=.data$Var1,y=.data$pred,col = "black"), size = 1)+
    ggplot2::scale_color_manual(values=c("red", "black"), labels=c("Predicted","Observed"), name = "")+
    facet_wrap(~as.factor(.data$Var2), nrow = 3, scales = "free")+
    scale_x_continuous(breaks = seq(1960,2018,11))+
    labs(x="Year", y="CRLs Ages 5-14")+
    theme_bw()+
    theme(legend.position = "none",text = element_text(size=18))

  #######################

  #crl residuals
  cc1<-  CRL.vec %>%
    ggplot(aes(x=.data$Var1, y = .data$Zresid))+
    ggplot2::geom_text(aes(label = .data$Var2), size = 1.2)+
    ggplot2::geom_hline(yintercept = 0)+
    geom_line(data =  CRL.vec   %>%
                dplyr::group_by(.data$Var1) %>%
                dplyr::summarise(yy=mean(.data$Zresid),.groups="drop"),
              aes(x=.data$Var1, y=.data$yy), color="red", size=1)+
    scale_x_continuous(breaks = seq(1960,2018,12))+
    labs(y="",x="Year")+
    theme_bw()

  cc2<-CRL.vec %>%
    ggplot(aes(x=.data$cohort, y = .data$Zresid))+
    geom_point(shape=3)+
    ggplot2:: geom_hline(yintercept = 0)+
    geom_line(data = CRL.vec %>%
                dplyr::group_by(.data$cohort) %>%
                dplyr::summarise(yy=mean(.data$Zresid),.groups="drop"),
              aes(x=.data$cohort, y=.data$yy), color="red", size=1)+
    scale_x_continuous(breaks = seq(1960,2018,12))+
    labs(y="",x="Cohort")+
    theme_bw()

  cc3<-CRL.vec %>%
    ggplot(aes(x=.data$Var2, y = .data$Zresid))+
    geom_point(shape=3)+
    ggplot2::geom_hline(yintercept = 0)+
    geom_line(data = CRL.vec %>%
                dplyr::group_by(.data$Var2) %>%
                dplyr::summarise(yy=mean(.data$Zresid),.groups="drop"),
              aes(x=.data$Var2, y=.data$yy), color="black", size=1)+
    labs(y="",x="Age")+
    theme_bw()

  cc4<-CRL.vec %>%
    ggplot(aes(x=.data$pred, y =.data$Zresid))+
    geom_point(shape=3)+
    ggplot2::geom_hline(yintercept = 0)+
    labs(y="",x="Pred")+
    theme_bw()

  resid_crl <- cowplot::plot_grid(cc1,cc2,cc3,cc4,ncol=1)

  #crl resid bubble plot
  crl_bplot= CRL.vec %>%   ggplot(aes(x = .data$Var1, y = .data$Var2,size = abs(.data$Zresid), color = .data$colr))+
    geom_point(shape = 20, alpha = 0.7)+
    ggplot2::scale_size(range = c(0,15))+
    scale_x_continuous(breaks = seq(1960,2018,5))+
    scale_y_continuous(breaks = seq(1,15,1))+
    ggplot2::scale_color_manual(values = c("deepskyblue", "firebrick2"))+
    guides(size = FALSE,color = FALSE)+
    labs(title="red bubbles are positive, blue are negative",y="Age",x="Year")+
    theme_bw()

  #landings plots
  catch.tot <-  data.frame(year=1960:endyear,obs = mfits$tmb.data$olandings,
                           pred = exp(mfits$rep$log_landings_pred),std.res = exp(mfits$rep$std_landings_resid))

  catch.tot <- catch.tot %>% dplyr::mutate(pdiff = 100*(catch.tot$obs/catch.tot$pred-1),
                                    upper =  mfits$tmb.data$landings_U,
                                    lower = mfits$tmb.data$landings_L)

  catch.tot <- catch.tot %>% dplyr::mutate(ratio = (catch.tot$pred-catch.tot$obs)/catch.tot$obs,
                                    bound =(catch.tot$upper-catch.tot$obs)/catch.tot$obs)

  land_plot <- catch.tot %>% ggplot(aes(x = .data$year, y = log(.data$pred)))+
    ggplot2::geom_ribbon(aes(ymin = log(.data$obs), ymax = log(.data$upper)), fill = "grey", alpha = 0.50)+
    geom_line(aes(x = .data$year, y = log(.data$pred)), size = 0.75,linetype = "solid")+
    scale_x_continuous(breaks = seq(1960,2018,4))+
    scale_y_continuous(breaks = seq(-1,6,0.50))+
    labs(y="Log catch", x = "Year")+
    theme_bw()

  land_plot2 <- catch.tot %>%  ggplot(aes(x = .data$year, y = .data$ratio))+
    geom_line(aes(x = .data$year, y = .data$ratio), size = 0.75,linetype = "solid")+
    geom_line(aes(x = .data$year, y = .data$bound), size = 0.75,linetype = "dashed",color="red")+
    scale_x_continuous(breaks = seq(1960,2018,4))+
    scale_y_continuous(breaks = seq(-1,6,0.50))+
    labs(y="(pred-obs)/obs", x = "Year")+
    theme_bw()

  #fishing mortality rates
  temp<-as.data.frame(mfits$rep$F,row.names = 1960:endyear) %>% stats::setNames(paste("Age",5:15,sep=""))
  F.Res<- reshape2::melt(as.matrix(temp))

  F_plot <- F.Res %>% ggplot(aes(x = .data$Var1, y = .data$value) )+
    geom_line(size = 1)+
    facet_wrap(~.data$Var2, nrow = 3)+
    labs(x ="Year", y = "Fishing mortality rates")+
    theme_bw()+
    theme(text = element_text(size=18))+
    scale_x_continuous(breaks = seq(1960,2018,11))

  #log fishing mortality deviations
  temp<-as.data.frame(mfits$rep$log_F_dev,row.names = 1960:endyear) %>% stats::setNames(paste("Age",5:15,sep=""))
  F.dev<- reshape2::melt(as.matrix(temp))

  F_devplot <- F.dev %>%
    ggplot(aes(x = .data$Var1, y = .data$value) )+
    geom_line(size = 1)+
    facet_wrap(~.data$Var2, nrow = 3)+
    labs(x ="Year", y = "Log F deviations")+
    theme_bw()+
    theme(text = element_text(size=18))+
    scale_x_continuous(breaks = seq(1960,2018,10))

  #process error
  temp<-as.data.frame(mfits$rep$pe,row.names = 1960:(endyear-1)) %>% stats::setNames(1:14)
  pe.vec<- reshape2::melt(as.matrix(temp))
  pe.vec$colr = 'deepskyblue'
  pe.vec$colr[pe.vec$value>0]='firebrick2'

  pe_plot <- pe.vec %>% ggplot( aes(x = .data$Var1, y = .data$value) )+
    geom_line(size = 0.75)+
    ggplot2::scale_color_manual(values=c("red", "black"), labels=c("Predicted","Observed"), name = "", position = "top")+
    facet_wrap(~.data$Var2)+
    scale_x_continuous(breaks = seq(1960,2018,10))+
    ggplot2::geom_hline(yintercept = 0, linetype = 2, size=1)+
    labs(x ="Year", y = "Process errors")+
    theme_bw()+
    theme(text = element_text(size=18))

  pe_bubbleplot <- pe.vec %>%  ggplot( aes(x = .data$Var1, y = .data$Var2,size = abs(.data$value), color = .data$colr))+
    geom_point(shape = 20, alpha = 0.7)+
    ggplot2::scale_size(range = c(0,15))+
    scale_x_continuous(breaks = seq(1960,2018,5))+
    scale_y_continuous(breaks = seq(1,15,1))+
    ggplot2::scale_color_manual(values = c("deepskyblue", "firebrick2"))+
    guides(size = FALSE,color = FALSE)+
    labs(title="red bubbles are positive, blue are negative",x="Year", y="Age")

  #catchability estimates
  q.runs <- data.frame(pred =as.numeric(mfits$rep$q[,1:14]), name  = as.character(mfits$map$m_q),
                       m_q=  as.numeric(mfits$rep$m_q))
  vtnames = 'm_q'
  ind = names(mfits$sdrep$value) %in% vtnames
  q.runs$SE=  mfits$sdrep$sd[ind]

  q.runs = q.runs %>%
    dplyr::mutate(survey = stringr::str_extract(.data$name, "[^_][:lower:]+"),
                  name2 = stringr::str_extract(.data$name, "[^_]+$"),
                  gear =  stringr::str_extract(.data$name, "[:upper:]") )

  q.runs = q.runs %>%
    dplyr::mutate(gear = replace(.data$gear, .data$gear!= "E", "C"))

  qa1<-ggplot()+
    geom_line(data =q.runs %>% dplyr::filter(!is.na(.data$name2)) %>%
    dplyr::filter(.data$survey =="Fall") ,aes(x = .data$name2, y = exp(.data$pred),group=.data$gear, linetype=.data$gear))+
    labs(x = "Fall", y = "")+
    ggplot2::ylim(0,8)+
    theme_bw() +
    theme(legend.position = "none")

  qa2<-ggplot()+
    geom_line(data =q.runs %>% dplyr::filter(!is.na(.data$name2)) %>%
    dplyr::filter(.data$survey =="Spr") ,aes(x = .data$name2, y = exp(.data$pred),group=.data$gear, linetype=.data$gear))+
    labs(x = "Spring", y = "")+
    ggplot2::ylim(0,8)+
    theme_bw() +
    ggplot2::scale_linetype_discrete(labels = c("Campelen", "Engel"))+
    theme(legend.position = "bottom")

  qq<- cowplot::plot_grid(qa1,qa2,ncol=1)

  ret = list(index_fit1 = index_fit1,index_fit2=index_fit2, resid_index = resid_index, pop_process = pp,
             index_bplot=index_bplot, crl_fit = crl_fit,resid_crl=resid_crl, crl_bplot=crl_bplot,
             landings = land_plot,landratio = land_plot2,F_at_age = F_plot, log_F_dev  =   F_devplot,
             process_errors = pe_plot, pe_bplot=pe_bubbleplot,catchability = qq )
return(ret)
}
