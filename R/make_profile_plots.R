#' Plot output from \code{\link{make_profile}} object
#'
#' Outputs profile plots
#' @importFrom ggplot2 ggplot aes geom_line
#' @param prof object from \code{\link{make_profile}}
#'
make.profile.plots = function(profile){

pdat2=setNames(reshape2::melt(profile$conditional_dev[1:4,]), c('component', 'del_m', 'nll'))
pdat3=setNames(reshape2::melt(profile$remainder_dev[1:4,]), c('component', 'del_m', 'nll'))
pdat4=setNames(reshape2::melt(profile$marginal_dev), c('component', 'del_m', 'nll'))

pp= ggplot(data = pdat2, aes(x=del_m, y=nll)) +
  geom_line(aes(colour=component),lwd=1)+
  geom_line(data = pdat3, aes(x = del_m, y = nll, color=component), linetype=2,lwd=1)
#
# pp = pp + labs(x = "M deviation",y = "nll deviation")+scale_y_continuous(breaks=seq(0,6000,2.5))+scale_x_continuous(breaks = seq(-0.5,0.5,0.05))
#
# pp2=pp+ coord_cartesian(xlim = c(-0.05,0.35), ylim = c(0,20))+ theme_bw()+theme(legend.position = "bottom",text = element_text(size=15))
# pp2 = pp2 + geom_line(data = pdat4,aes(x=del_m, y=nll, color = component),lwd=1 )
# pp2 = pp2+  scale_color_manual(values= c("#FF7F00","#E41A1C","#984EA3","#377EB8","#4DAF4A","#FFFF33"))
#

ret = list( nll_dev_plot = pp)
# ret=NULL
return(ret)
}
