#' make.parm: prepares parameters for use in APAM
#'
#' Prepares parmameterss for use in APAM.
#'
#' @param data object returned from \code{\link{make.tmb.data}}
#' @param map object returned from \code{\link{make.map}}
#' @param parms (optional) list, can be used to change default starting parameters. See details.
#' @param parmsL (optional) list, can be used to change default parameter lower bounds. See details.
#' @param parmsU (optional) list, can be used to change default parameter upper bounds.. See details.
#' @param no.pe  (optional) if \code{no.pe=T} \code{\link{make.map}} set \code{TRUE} to fix \code{sdpe=0.01}
#' @param no.logits (optional) if \code{no.logits=T} \code{\link{make.map}} set \code{TRUE} to fix \code{logits}. See details.
#' @param no.Flogits (optional) if \code{no.Flogits=T} \code{\link{make.map}} set \code{TRUE} to fix \code{Flogits}. See details.
#'
#'@details
#'   \describe{
#'     \item{\code{parms}}{list that contains starting parameter values. Default \code{parm=}
#'     \itemize{
#'     \item \code{log_R=c(10,10),}
#'     \item \code{m_q=matrix(c(rep(log(0.1),t1),rep(-30,14-t1)),nrow=length(data$isurvey1), ncol = 14,byrow=T),}
#'     \item \code{log_cv_index=rep(0.5,length(unique(data$isd))),}
#'     \item \code{log_std_log_R=log(1),}
#'     \item \code{log_F_mean = matrix(log(0.1),nrow=Y,ncol=2,byrow=T),}
#'     \item \code{log_std_log_F = matrix(log(0.1),nrow=Y,ncol=A-4,byrow=T),}
#'     \item \code{log_std_pe = matrix(log(0.2),nrow=Y-1,ncol=A-1),}
#'     \item \code{log_std_crl = matrix(log(0.2),nrow=Y,ncol=A-5),}
#'     \item \code{logit_ar_index_age = rep(0,length(unique(data$isurvey1))),}
#'     \item \code{logit_ar_logF = c(0, 0, 0),}
#'     \item \code{logit_ar_crl = c(2.2,-10),}
#'     \item \code{logit_ar_pe = rep(-10,2),}
#'     \item \code{logit_ar_logRec = 0,}
#'     \item \code{log_F_devt=t(matrix(log(0.001),nrow=Y,ncol=A-4,byrow=T)}
#'     \item \code{log_Nt=t(matrix(log(10000),nrow=Y,ncol=A,byrow=T)),}
#'     \item \code{h=0}}}
#'     \item{\code{parmsL}}{list that contains starting parameter values. Default \code{parmL=}
#'     \itemize{
#'     \item \code{log_R=c(-10,-10),}
#'     \item \code{m_q=matrix(-30,nrow=length(data$isurvey1), ncol = A-1),}
#'     \item \code{log_cv_index=rep(log(0.001),length(unique(data$isd))),}
#'     \item \code{log_std_log_R=-10,}
#'     \item \code{log_F_mean = matrix(log(0.0001),nrow=Y,ncol=2,byrow=T),}
#'     \item \code{log_std_log_F =matrix(-Inf,nrow=Y,ncol=A-4,byrow=T),}
#'     \item \code{log_std_pe = matrix(log(0.001),nrow=Y-1,ncol=A-1),}
#'     \item \code{log_std_crl = matrix(log(0.001),nrow=Y,ncol=A-5),}
#'     \item \code{logit_ar_index_age = rep(-5,length(unique(data$isurvey1))),}
#'     \item \code{logit_ar_logF = rep(-5,3),}
#'     \item \code{logit_ar_crl = rep(-5,2),}
#'     \item \code{logit_ar_pe = rep(-5,2),}
#'     \item \code{logit_ar_logRec = -5}}}
#'     \item{\code{parmsu}}{list that contains starting parameter values. Default \code{parmU=}
#'     \itemize{
#'     \item \code{log_R=c(Inf,Inf),}
#'     \item \code{m_q=matrix(Inf,nrow=length(data$isurvey1), ncol = A-1),}
#'     \item \code{log_cv_index=rep(log(5),length(unique(data$isd))),}
#'     \item \code{log_std_log_R=Inf,}
#'     \item \code{log_F_mean = matrix(Inf,nrow=Y,ncol=2,byrow=T),}
#'     \item \code{log_std_log_F = matrix(Inf,nrow=Y,ncol=A-4,byrow=T),}
#'     \item \code{log_std_pe = matrix(10,nrow=Y-1,ncol=A-1),}
#'     \item \code{log_std_crl = matrix(5,nrow=Y,ncol=A-5),}
#'     \item \code{logit_ar_index_age = rep(5,length(unique(data$isurvey1))),}
#'     \item \code{logit_ar_logF = c(5,5,5),}
#'     \item \code{logit_ar_crl = c(5,5),}
#'     \item \code{logit_ar_pe = c(10,10),}
#'     \item \code{logit_ar_logRec = 5}}}
#'     \item{\code{no.logits }}{if \code{TRUE}, all survey logits = 0.73, all F logits = 0.95 and rec logit = 0.32 }
#'     \item{\code{no.Flogits}}{if \code{TRUE}, all F logits = 0.95}}
#'
#' @examples
#' \dontrun{
#' params <- make.parm(data,map)
#'}
#' @export
make.parm = function(data,map,parms=NULL, parmsL=NULL, parmsU=NULL, no.pe=F,no.logits=F,no.Flogits=F){

  #define variables
  Y<-data$Y
  A<-data$A
  t1<- (length(unique(levels(map$m_q)))-8)/2

  #set parameter starting values
  if(is.null(parms)){parms <- list(
    log_R=c(10,10),
    m_q=matrix(c(rep(log(0.1),t1),rep(-30,14-t1)),nrow=length(data$isurvey1), ncol = 14,byrow=T),
    log_cv_index=rep(0.5,length(unique(data$isd))),
    log_std_log_R=log(1),
    log_F_mean = matrix(log(0.1),nrow=Y,ncol=2,byrow=T),
    log_std_log_F = matrix(log(0.1),nrow=Y,ncol=A-4,byrow=T),
    log_std_pe = matrix(log(0.2),nrow=Y-1,ncol=A-1),
    log_std_crl = matrix(log(0.2),nrow=Y,ncol=A-5),
    logit_ar_index_age = rep(0,length(unique(data$isurvey1))),
    logit_ar_logF = c(0, 0, 0),
    logit_ar_crl = c(2.2,-10),
    logit_ar_pe = rep(-10,2),
    logit_ar_logRec = 0,
    log_F_devt=t(matrix(log(0.001),nrow=Y,ncol=A-4,byrow=T)),
    log_Nt=t(matrix(log(10000),nrow=Y,ncol=A,byrow=T)),
    h=0
  )}

  #formulation with fixed logits
  if(no.logits){
    parms$logit_ar_index_age <- rep(1,length(unique(data$isurvey1)))
    parms$logit_ar_logF <- c(3, 3, 3)
    parms$logit_ar_logRec <- -0.75
  }

  #formulation with fixed Flogits
  if(no.Flogits){
    parms$logit_ar_logF <- c(3, 3, 3)
  }

  #formulation with process error turned off
  if(no.pe){parms$log_std_pe <- matrix(log(0.01),nrow=Y-1,ncol=A-1)}

  #set lower bounds
  if(is.null(parmsL)){parmsL <- list(
    log_R=c(-10,-10),
    m_q=matrix(-30,nrow=length(data$isurvey1), ncol = A-1),
    log_cv_index=rep(log(0.001),length(unique(data$isd))),
    log_std_log_R=-10,
    log_F_mean = matrix(log(0.0001),nrow=Y,ncol=2,byrow=T),
    log_std_log_F =matrix(-Inf,nrow=Y,ncol=A-4,byrow=T),
    log_std_pe = matrix(log(0.001),nrow=Y-1,ncol=A-1),
    log_std_crl = matrix(log(0.001),nrow=Y,ncol=A-5),
    logit_ar_index_age = rep(-5,length(unique(data$isurvey1))),
    logit_ar_logF = rep(-5,3),
    logit_ar_crl = rep(-5,2),
    logit_ar_pe = rep(-5,2),
    logit_ar_logRec = -5
  )}

  #set upper bounds
  if(is.null(parmsU)){parmsU <- list(
    log_R=c(Inf,Inf),
    m_q=matrix(Inf,nrow=length(data$isurvey1), ncol = A-1),
    log_cv_index=rep(log(5),length(unique(data$isd))),
    log_std_log_R=Inf,
    log_F_mean = matrix(Inf,nrow=Y,ncol=2,byrow=T),
    log_std_log_F = matrix(Inf,nrow=Y,ncol=A-4,byrow=T),
    log_std_pe = matrix(10,nrow=Y-1,ncol=A-1),
    log_std_crl = matrix(5,nrow=Y,ncol=A-5),
    logit_ar_index_age = rep(5,length(unique(data$isurvey1))),
    logit_ar_logF = c(5,5,5),
    logit_ar_crl = c(5,5),
    logit_ar_pe = c(10,10),
    logit_ar_logRec = 5
  )}

  #to set lower bounds for maps
  tp=parmsL;
  tp$logit_ar_crl = NULL;
  tp$logit_ar_pe = NULL;
  tp$log_F_mean=rep(log(0.000001),length(unique(as.vector(map$log_F_mean))));
  tp$log_std_log_F=rep(-Inf,length(unique(map$log_std_log_F)));
  tp$log_cv_index=rep(-Inf,length(unique(map$log_cv_index)));
  tp$m_q=rep(-30,length(unique(map$m_q))-1);
  tp$log_std_crl=rep(log(0.001),length(unique(map$log_std_crl)));
  tp$log_std_pe=rep(log(0.001),length(unique(map$log_std_pe)));

  if(no.logits){tp$logit_ar_index_age<- NULL;
  tp$logit_ar_logF  <- NULL;
  tp$logit_ar_logRec  <- NULL;}

  if(no.Flogits){tp$logit_ar_logF  <- NULL;}


  if(no.pe){tp$log_std_pe<-NULL;}

  lower = unlist(tp);

  #to set upper bounds for maps
  tp=parmsU;
  tp$logit_ar_crl <- NULL;
  tp$logit_ar_pe <- NULL;
  tp$log_F_mean=rep(log(10),length(unique(as.vector(map$log_F_mean))));
  tp$log_std_log_F=rep(log(5),length(unique(map$log_std_log_F)));
  tp$log_cv_index=rep(log(5),length(unique(map$log_cv_index)));
  tp$m_q=rep(Inf,length(unique(map$m_q))-1);
  tp$log_std_crl=rep(log(10),length(unique(map$log_std_crl)));
  tp$log_std_pe=rep(log(10),length(unique(map$log_std_pe)));



  if(no.logits){tp$logit_ar_index_age <- NULL;
  tp$logit_ar_logF <- NULL;
  tp$logit_ar_logRec  <- NULL;}

  if(no.Flogits){tp$logit_ar_logF <- NULL;}

  if(no.pe){tp$log_std_pe<-NULL;}

  upper = unlist(tp);

  parms$resid_index_res <- rep(0, length(data$index))
  parms$resid_crl_res  <- matrix(0, nrow = 58, ncol = 10)

  ret = list(parms=parms,parmsL=parmsL,parmsU=parmsU,lower=lower,upper=upper)

  return(ret)


}

