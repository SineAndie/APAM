#' #' make.curv.plots: plot curvature
#' #'
#' #' Plot output from curvature to check suitability of local influence diagnostics.
#' #'
#' #' @param LocInf object from \code{\link{make.LI}}
#' #' @param curv object from \code{\link{make.curv}}
#' #'
#' #' @examples
#' #' \dontrun{
#' #' curv.plots <- make.curv.plots(mfits, curv)
#' #' }
#' #' @export
#' #'
#' #'
#' make.curv.plot = function(LocInf,curv){
#'
#'   LI_dat= reshape2::melt(as.matrix(LocInf$LI[,5]))
#'
#'   if(LocInf$type=="age"){
#'     Si_plot=    Si_dat %>%
#'       ggplot(aes(x= .data$Var1, y = .data$value)) +
#'       geom_line()+
#'       scale_x_continuous(breaks = seq(1,15,1))+
#'       facet_wrap(~.data$Var2)+
#'       labs(y="Total slope",x="Age")+
#'       ggplot2::geom_hline(yintercept = 0, linetype=2)+
#'       theme_bw()+theme(legend.position = "none",text = element_text(size=15))
#'
#'   }
#'
#'   if(LocInf$type=="year"){
#'     Si_plot = Si_dat %>%
#'       ggplot(aes(x= .data$Var1, y = .data$value)) +
#'       geom_line()+
#'       scale_x_continuous(breaks = seq(1960,2019,8))+
#'       facet_wrap(~.data$Var2)+
#'       labs("Total slope",x="Year")+
#'       ggplot2::geom_hline(yintercept = 0, linetype=2)+
#'       theme_bw()+theme(legend.position = "none",text = element_text(size=15))
#'   }
#'
#'     if(age==TRUE){
#'       Si_age = data.frame( Si = unlist(curv), Age = Aname[pert])
#'       Si_age$LI <- unlist(LocInf$results$Si)
#'
#'       Si_agep =   Si_age%>%
#'         ggplot(aes(x= Age, y = Si)) +
#'         geom_line()+
#'         scale_x_continuous(breaks = seq(1,15,1))+
#'         ylab("Curvature")+
#'         geom_hline(yintercept = 0, linetype=2)+
#'         theme_bw()+theme(legend.position = "none",text = element_text(size=15))
#'
#'       cur_vs_LI <-ggplot(Si_age, aes(x=Si, y=LI))+
#'         geom_point(alpha=0.50)+
#'         ylab("Influence slope")+
#'         xlab("Curvature")+
#'         geom_hline(yintercept = 0,lty=2)+
#'         geom_hline(yintercept = 1,lty=2,color="red")+
#'         geom_hline(yintercept = -1,lty=2,color="red")+
#'         theme_bw()+theme(legend.position = "bottom",text = element_text(size=15))
#'
#'       Curvature = list(results = Curvature,  Curvature_plot=Si_agep, Curv_vs_LI_plot =cur_vs_LI)
#'     }
#'
#'     if(year==TRUE){
#'       Si_year = data.frame( Si = unlist(curv),  Year = Yname[pert])
#'       Si_year$LI <- unlist(LocInf$results$Si)
#'
#'       Si_yearp =
#'         ggplot(data = Si_year, aes(x= Year, y = Si)) +
#'         geom_line()+
#'         scale_x_continuous(breaks = seq(1960,2019,8))+
#'         ylab("Curvature")+
#'         geom_hline(yintercept = 0, linetype=2)+
#'         theme_bw()+theme(legend.position = "none",text = element_text(size=15))
#'
#'       cur_vs_LI <-ggplot(Si_year, aes(x=Si, y=LI))+
#'         geom_point(alpha=0.50)+
#'         ylab("Influence slope")+
#'         xlab("Curvature")+
#'         geom_hline(yintercept = 0,lty=2)+
#'         geom_hline(yintercept = 1,lty=2,color="red")+
#'         geom_hline(yintercept = -1,lty=2,color="red")+
#'         theme_bw()+theme(legend.position = "bottom",text = element_text(size=15))
#'
#'       Curvature = list(results = Curvature,  Curvature_plot=Si_yearp, Curv_vs_LI_plot =cur_vs_LI)
#'     }}
#'   return(Curvature)
#'  }
