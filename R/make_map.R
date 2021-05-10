#' make.map: prepares map for APAM
#'
#' Prepares map for use in APAM.
#'
#' @param data object returned from \code{\link{make.tmb.data}}
#' @param setmap (optional) list, can be used to change default map settings. See details.
#' @param crl.block (optional) T/F, turn on/off crl sd split? Default = \code{T}. See details.
#' @param no.pe (optional) T/F, turn on/off process errors? Default = \code{F}.
#' @param no.logits (optional) T/F, turn on/off all logit estimation? Default = \code{F}.
#' @param no.Flogits  (optional) T/F, turn on/off Flogit estimation? Default = \code{F}.
#'
#'@details
#'   \describe{
#'     \item{\code{setmap}}{list that contains parameter mapping. Default \code{setmap=}
#'     \itemize{
#'     \item \code{meanF= c("5","6"),}
#'     \item \code{stdF = c("5",rep("6+",length =   data$A-5)),}
#'     \item \code{ageFall = c("1",rep("2-11",10),rep("12-15",4)),}
#'     \item \code{ageSpring = c("1","2",rep("3-13",11),rep("14-15",2)),}
#'     \item \code{ageSpanish = NULL,}
#'     \item \code{stdcrl  = c(rep("5-6",2),rep("7-11",5),rep("12-14",3)),}
#'     \item \code{stdpe = rep("all", data$A-1),}
#'     \item \code{mapq = c(1:7,rep(NA,length = (data$A-1)-7))} }}
#'     \item{\code{crl.block}}{if \code{TRUE}, a separate crl sd is estimated pre/post 1993.}
#'     \item{\code{no.pe}}{if \code{TRUE}, sdpe is mapped out.}
#'     \item{\code{no.logits }}{if \code{TRUE}, all logit parameters are mapped out.}
#'     \item{\code{no.Flogits}}{if \code{TRUE}, all F logit parameters are mapped out. }}
#'
#'
#' @examples
#' \dontrun{
#' map <- make.map(data)
#' }
#' @export
make.map = function(data, setmap=NULL,crl.block=TRUE,no.pe=FALSE,no.logits=FALSE,no.Flogits=FALSE){


  if(is.null(setmap)){
    setmap <- list(
      meanF= c("5","6"),
      stdF = c("5",rep("6+",length =   data$A-5)),
      ageFall = c("1",rep("2-11",10),rep("12-15",4)),
      ageSpring = c("1","2",rep("3-13",11),rep("14-15",2)),
      ageSpanish = NULL,
      stdcrl  = c(rep("5-6",2),rep("7-11",5),rep("12-14",3)),
      stdpe = rep("all", data$A-1),
      mapq = c(1:7,rep(NA,length = (data$A-1)-7))
    )
  }


  #mean F split pre/post 1995
  mapF <- matrix(NA,nrow=data$Y,ncol=2,byrow=T)
  mapF[1:35,] <- matrix(paste('F_',setmap[[1]],'_pre',sep=''),nrow=35,ncol=2,byrow=T)
  mapF[36:data$Y,] <- matrix(paste('F_',setmap[[1]],'_post',sep=''),nrow=length(36:data$Y),ncol=2,byrow=T)

  #sd F split
  map_log_std_logF <- matrix(setmap[[2]],nrow=data$Y,ncol=data$A-4,byrow=T)

  #cv index split
  map_index_cv <- c(paste('Fall_',setmap[[3]],sep=""),
                    # paste('Span_',setmap[[5]],sep=""),
                    paste('Sprg_',setmap[[4]],sep=""))

  #crl sd split with separate year block
  map_crl_vec <- matrix(setmap[[6]],nrow=data$Y,ncol=data$A-5,byrow=T)
  if(crl.block){map_crl_vec[31:40,] <- "90-99"}

  #sd pe
  map_std_pe <- matrix(setmap[[7]],nrow=data$Y-1,ncol=data$A-1)
  if(no.pe){map_std_pe <- matrix(NA,nrow=data$Y-1,ncol=data$A-1)}

  #q split with age aggregation
  nq_NA <- length(which(is.na(setmap[[8]])))
  nq <-14 - nq_NA

  #identifier survey/year
  q1 <- matrix(paste('Fall_',1:nq,sep=""),nrow=data$NsF, ncol =nq,byrow=T)
  q3 <- matrix(paste('Spr_',1:nq,sep=""),nrow=(data$Ns - data$NsF - data$NsSpan), ncol = nq,byrow=T)

  n1_NA <- matrix(NA,nrow=data$Ns, ncol = nq_NA, byrow=T)

  map_q <- rbind(q1,q3)
  map_q <- cbind(map_q, n1_NA)

  #Engels/Campelen trawl split
  map_q[1:5,1:4] <- paste('E',map_q[1:5,1:4],sep="_")
  map_q[(data$NsF+data$NsSpan+1):(data$NsF+data$NsSpan+11),1:4] <- paste('E',map_q[(data$NsF+data$NsSpan+1):(data$NsF+data$NsSpan+11),1:4],sep="_")

  map = list(
    logit_ar_crl = factor(c(NA,NA)),
    logit_ar_pe = factor(c(NA,NA)),
    log_F_mean = factor(mapF),
    log_std_log_F = factor(map_log_std_logF),
    log_cv_index = factor(map_index_cv),
    log_std_crl = factor(map_crl_vec),
    log_std_pe = factor(map_std_pe),
    m_q = factor(map_q),
    h = factor(NA),
    resid_index_res  = factor(rep(NA, length(data$index))),
    resid_crl_res = factor(matrix(NA, nrow = 58, ncol = 10))

  )

  if(no.logits){
    map$logit_ar_index_age <- factor(rep(NA,3))
    map$logit_ar_logF <- factor(rep(NA,3))
    map$logit_ar_logRec <- factor(NA)
  }


  if(no.Flogits){
    map$logit_ar_logF <- factor(rep(NA,3))
  }

  return(map)

  }
