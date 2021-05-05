#' #' #' Plot output from \code{\link{make_LI}} object
#' #' #'
#' #' #' Outputs plots of interest from APAM
#' #' #' @importFrom ggplot2 ggplot aes geom_line scale_x_continuous labs theme_bw
#' #' #' @importFrom ggplot2 theme geom_ribbon
#' #'
#' #' @param mfits object from make_fit
#' #' @param pert define which perturbations to run; default=all (for all ages/years)
#' #' @param all  to plot all influence slopes
#' #' @param age to plot age influence slopes
#' #' @param year to plot year influence slopes
#'
#'  make.plots = function(mfits,all=T, age=F,year=F){
#'
#'   Aname <- c(1:15)
#'   Yname <- c(1960:2017)
#'
#'   # if(all){
#'   #
#'   #   Si_matnll <- matrix(unlist(Si), nrow=Y, ncol=A,byrow=F)
#'   #   colnames(Si_matnll) <- colnames(tmb_data$weight)
#'   #   rownames(Si_matnll) <- rownames(tmb_data$weight)
#'   #   Si_matnll <- as.data.frame(Si_matnll)
#'   #   Si_matnll$Year <- as.numeric(rownames(Si_matnll))
#'   #
#'   #   Si_vecnll <- vec_func(Si_matnll)
#'   #   Si_vecnll$index <- Si_vecnll$index
#'   #   Si_vecnll$colr <- 'deepskyblue'
#'   #   Si_vecnll$colr[Si_vecnll$index>0] <- 'firebrick2'
#'   #
#'   #   Si_nll = Si_vecnll%>%ggplot(aes(x= Year, y = Age,size=abs(index))) +
#'   #     geom_point(aes(shape = colr),alpha=0.70)+
#'   #     scale_size( range = c(1,15))+
#'   #     scale_x_continuous(breaks = seq(1960,2018,5))+
#'   #     scale_y_continuous(breaks = seq(1,15,1))+
#'   #     scale_shape_manual(values=c("-","+"))+
#'   #     ylab("Age")+
#'   #     xlab("Year")+
#'   #     guides(size = FALSE,color = FALSE)+
#'   #     theme_bw()+theme(legend.position = "none",text = element_text(size=15))
#'   #
#'   #   #line plot of slopes
#'   #   Si_matnll <- matrix(unlist(Si), nrow=Y, ncol=A,byrow=F)
#'   #   Si_year <- apply(Si_matnll, 1,sum)
#'   #   Si_year <- data.frame( Si = Si_year,  Year = Yname)
#'   #
#'   #   Si_yearp =
#'   #     ggplot(data = Si_year, aes(x= Year, y = Si)) +
#'   #     geom_line()+
#'   #     scale_x_continuous(breaks = seq(1960,2019,8))+
#'   #     ylab("Total slope")+
#'   #     geom_hline(yintercept = 0, linetype=2)+
#'   #     theme_bw()+theme(legend.position = "none",text = element_text(size=15))
#'   #
#'   #   Si_age <- apply(Si_matnll, 2,sum)
#'   #   Si_age <- data.frame( Si = Si_age, Age = c(1:15))
#'   #
#'   #   Si_agep =   Si_age%>%
#'   #     ggplot(aes(x= Age, y = Si)) +
#'   #     geom_line()+
#'   #     scale_x_continuous(breaks = seq(1,15,1))+
#'   #     # scale_y_continuous(breaks = seq(-1000,1000,10))+
#'   #     ylab("Total slope")+
#'   #     geom_hline(yintercept = 0, linetype=2)+
#'   #     theme_bw()
#'   #
#'   #
#'   #   LI = list(results = LI,  Si_plot=Si_nll, Si_age_plot=Si_agep, Si_year_plot=Si_yearp)
#'   # }
#'
#'   if(age){
#'     Si_age = data.frame( Si = unlist(Si), Age = Aname[pert])
#'
#'     Si_agep =   Si_age%>%
#'       ggplot(aes(x= Age, y = Si)) +
#'       geom_line()+
#'       scale_x_continuous(breaks = seq(1,15,1))+
#'       ylab("Total slope")+
#'       geom_hline(yintercept = 0, linetype=2)+
#'       theme_bw()+theme(legend.position = "none",text = element_text(size=15))
#'
#'     LI = list(results = LI,  Si_age_plot=Si_agep)
#'   }
#'
#'   if(year){
#'     Si_year = data.frame( Si = unlist(Si),  Year = Yname[pert])
#'
#'     Si_yearp =
#'       ggplot(data = Si_year, aes(x= Year, y = Si)) +
#'       geom_line()+
#'       scale_x_continuous(breaks = seq(1960,2019,8))+
#'       ylab("Total slope")+
#'       geom_hline(yintercept = 0, linetype=2)+
#'       theme_bw()+theme(legend.position = "none",text = element_text(size=15))
#'
#'     LI = list(results = LI, Si_year_plot=Si_yearp)
#'   }
#'   return(Si_age_plot)
#' }
