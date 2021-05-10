#' make.profile.plots: plot output from profile likelihoods
#'
#' Plot outputs from APAM profile likelihoods
#'
#' @param prof object from \code{\link{make.profile}}
#'
#' @return a list containing all ggplots
#'   \describe{
#'     \item{\code{$nll_plot}}{profile likelihoods vs M perturbations}
#'     \item{\code{$nll_dev_plot}}{profile likelihood deviations vs M perturbations}
#'     }
#'
#' @examples
#' \dontrun{
#' prof.plots<- make.profile.plots(prof)
#' }
#' @export
make.profile.plots = function(prof){

  #full profiles
  pdat1=stats::setNames(reshape2::melt(prof$conditional[1:4,]), c('component', 'del_m', 'nll'))
  pdat2=stats::setNames(reshape2::melt(prof$remainder[1:4,]), c('component', 'del_m', 'nll'))
  pdat3=stats::setNames(reshape2::melt(prof$marginal), c('component', 'del_m', 'nll'))

  pp= ggplot(data = pdat1 %>% dplyr::filter(.data$component!="Spanish index"),
             aes(x=.data$del_m, y=.data$nll)) +
    geom_line(aes(colour=.data$component),lwd=1)+
    geom_line(data = pdat2%>% dplyr::filter(.data$component!="Spanish index"),
              aes(x = .data$del_m, y = .data$nll, color=.data$component), linetype=2,lwd=1)+
    labs(x = "M deviation",y = "nll deviation")+
    scale_x_continuous(breaks = seq(-0.5,0.5,0.05))+
    ggplot2::coord_cartesian(xlim = c(-0.15,0.40))+
    theme_bw()+theme(legend.position = "bottom",text = element_text(size=15))

  pp2 = pp + geom_line(data = pdat3,aes(x=.data$del_m, y=.data$nll, color = .data$component),lwd=1 )+
    ggplot2::scale_color_manual(values= c("#FF7F00","#E41A1C","#984EA3","#377EB8","#4DAF4A","#FFFF33"))

  #plot nll deviations
  pdat4=stats::setNames(reshape2::melt(prof$conditional_dev[1:4,]), c('component', 'del_m', 'nll'))
  pdat5=stats::setNames(reshape2::melt(prof$remainder_dev[1:4,]), c('component', 'del_m', 'nll'))
  pdat6=stats::setNames(reshape2::melt(prof$marginal_dev), c('component', 'del_m', 'nll'))

  pp3= ggplot(data = pdat4 %>% dplyr::filter(.data$component!="Spanish index"),
              aes(x=.data$del_m, y=.data$nll)) +
    geom_line(aes(colour=.data$component),lwd=1)+
    geom_line(data = pdat5%>% dplyr::filter(.data$component!="Spanish index"),
              aes(x = .data$del_m, y = .data$nll, color=.data$component), linetype=2,lwd=1)+
    labs(x = "M deviation",y = "nll deviation")+scale_y_continuous(breaks=seq(0,6000,2.5))+
    scale_x_continuous(breaks = seq(-0.5,0.5,0.05))+
    ggplot2::coord_cartesian(xlim = c(-0.15,0.40), ylim = c(0,20))+
    theme_bw()+theme(legend.position = "bottom",text = element_text(size=15))

  pp4 = pp3 + geom_line(data = pdat6,aes(x=.data$del_m, y=.data$nll, color = .data$component),lwd=1 )+
    ggplot2::scale_color_manual(values= c("#FF7F00","#E41A1C","#984EA3","#377EB8","#4DAF4A","#FFFF33"))

  ret = list( nll_plot = pp2, nll_dev_plot = pp4)

  return(ret)
}
