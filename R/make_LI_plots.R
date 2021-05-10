#' make.LI.plots: plot output from local influence diagnostics
#'
#' Plot outputs from APAM local influence diagnostics.
#'
#' @param LocInf object from \code{\link{make.LI}}
#'
#' @return a list containing  ggplot
#'   \describe{
#'     \item{\code{$LI_plot}}{local influence plot}
#'     }
#'
#' @examples
#' \dontrun{
#' LI.plots <- make.LI.plots(LocInf)
#' }
#' @export
 make.LI.plots = function(LocInf){

  if(LocInf$type=="age"){
    Si_dat= reshape2::melt(as.matrix(LocInf$LI))
    Si_plot=    Si_dat %>%
      ggplot(aes(x= .data$Var1, y = .data$value)) +
      geom_line()+
      scale_x_continuous(breaks = seq(1,15,1))+
      facet_wrap(~.data$Var2)+
      labs(y="Total slope",x="Age")+
      ggplot2::geom_hline(yintercept = 0, linetype=2)+
      theme_bw()+theme(legend.position = "none",text = element_text(size=15))

  }

  if(LocInf$type=="year"){
    rownames(LocInf$LI) = c(1960:2017)
    Si_dat= reshape2::melt(as.matrix(LocInf$LI))
    Si_plot = Si_dat %>%
      ggplot(aes(x= .data$Var1, y = .data$value)) +
      geom_line()+
      scale_x_continuous(breaks = seq(1960,2019,8))+
      facet_wrap(~.data$Var2)+
      labs(y="Total slope",x="Year")+
      ggplot2::geom_hline(yintercept = 0, linetype=2)+
      theme_bw()+theme(legend.position = "none",text = element_text(size=15))
  }

  if(LocInf$type=="all"){
    temp <- as.data.frame(matrix(unlist(LocInf$LI), nrow=58, ncol=15,byrow=F),
                            row.names = 1960:2017) %>% stats::setNames(1:15)
    Si_dat<- reshape2::melt(as.matrix(temp)) %>%
              dplyr::mutate(colr = 'deepskyblue')

    Si_dat$colr = 'deepskyblue'
    Si_dat$colr[Si_dat$value>0]='firebrick2'

    Si_plot = Si_dat%>%ggplot(aes(x= .data$Var1, y = .data$Var2,size=abs(.data$index),color=.data$colr)) +
      geom_point(shape = 20, alpha = 0.7)+
      ggplot2::scale_size(range = c(0,15))+
      scale_x_continuous(breaks = seq(1960,2018,5))+
      scale_y_continuous(breaks = seq(1,15,1))+
      ggplot2::scale_color_manual(values = c("deepskyblue", "firebrick2"))+
      guides(size = FALSE,color = FALSE)+
      labs(y="Age",x="Year")+
      theme_bw()
  }

  return(ret=list(LI_plot = Si_plot))
}
