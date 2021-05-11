#' make.curv.plots: plot curvature
#'
#' Plot output from curvature.
#'
#' @param LocInf object from \code{\link{make.LI}}
#' @param curv object from \code{\link{make.curv}}
#'
#' @return a list containing  ggplot
#'   \describe{
#'     \item{\code{$Curv_vs_LI_plot}}{plots curvature vs influence slopes}
#'     }
#'
#' @examples
#' \dontrun{
#' curv.plots <- make.curv.plots(mfits, curv)
#' }
#' @export
make.curv.plots = function(LocInf,curv){

  LI_dat= reshape2::melt(as.matrix(LocInf$LI[,5]))
  LI_dat$curv = unlist(curv$curv)

  cur_vs_LI <-  LI_dat %>% ggplot( aes(x=.data$curv, y=.data$value))+
    geom_point(alpha=0.50)+
    labs(y="Influence slope",x="Curvature")+
    ggplot2::geom_hline(yintercept = 0,lty=2)+
    ggplot2::geom_hline(yintercept = 1,lty=2,color="red")+
    ggplot2::geom_hline(yintercept = -1,lty=2,color="red")+
    theme_bw()+theme(legend.position = "bottom",text = element_text(size=15))

  Curvature = list(Curv_vs_LI_plot =cur_vs_LI)

  return(Curvature)
}
