#' Prepare parameters for APAM
#'
#' Prepares parameters for use in APAM
#'
#' @importFrom rlang .data
#' @param data input from make_tmb_data
#' @param parmap input from make_map
#' @param parameters can redefine starting parameter values
#' @param parmL can redefinie lower parameter bound
#' @param parmU can redefinie upper parameter bound
#' @param no.pe  turn on/off pe, i.e fix pe sd to 0.01
#' @param no.logits turn on/off logits
#' @param no.Flogits turn on/off Flogit estimation


make.parm = function(data,parmap,parameters=NULL, parmL=NULL, parmU=NULL, no.pe=F,no.logits=F,no.Flogits=F){

  #define variables
  Y<-data$Y
  A<-data$A
  t1<- (length(unique(levels(parmap$m_q)))-8)/2

  #set parameter starting values
  if(is.null(parameters)){parameters <- list(
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
    parameters$logit_ar_index_age <- rep(1,length(unique(data$isurvey1)))
    parameters$logit_ar_logF <- c(3, 3, 3)
    parameters$logit_ar_logRec <- -0.75
  }

  #formulation with fixed Flogits
  if(no.Flogits){
    parameters$logit_ar_logF <- c(3, 3, 3)
  }

  #formulation with process error turned off
  if(no.pe){parameters$log_std_pe <- matrix(log(0.01),nrow=Y-1,ncol=A-1)}

  #set lower bounds
  if(is.null(parmL)){parmL <- list(
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
  if(is.null(parmU)){parmU <- list(
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
  tp=parmL;
  tp$logit_ar_crl = NULL;
  tp$logit_ar_pe = NULL;
  tp$log_F_mean=rep(log(0.000001),length(unique(as.vector(parmap$log_F_mean))));
  tp$log_std_log_F=rep(-Inf,length(unique(parmap$log_std_log_F)));
  tp$log_cv_index=rep(-Inf,length(unique(parmap$log_cv_index)));
  tp$m_q=rep(-30,length(unique(parmap$m_q))-1);
  tp$log_std_crl=rep(log(0.001),length(unique(parmap$log_std_crl)));
  tp$log_std_pe=rep(log(0.001),length(unique(parmap$log_std_pe)));

  if(no.logits){tp$logit_ar_index_age<- NULL;
  tp$logit_ar_logF  <- NULL;
  tp$logit_ar_logRec  <- NULL;}

  if(no.Flogits){tp$logit_ar_logF  <- NULL;}


  if(no.pe){tp$log_std_pe<-NULL;}

  lower = unlist(tp);

  #to set upper bounds for maps
  tp=parmU;
  tp$logit_ar_crl <- NULL;
  tp$logit_ar_pe <- NULL;
  tp$log_F_mean=rep(log(10),length(unique(as.vector(parmap$log_F_mean))));
  tp$log_std_log_F=rep(log(5),length(unique(parmap$log_std_log_F)));
  tp$log_cv_index=rep(log(5),length(unique(parmap$log_cv_index)));
  tp$m_q=rep(Inf,length(unique(parmap$m_q))-1);
  tp$log_std_crl=rep(log(10),length(unique(parmap$log_std_crl)));
  tp$log_std_pe=rep(log(10),length(unique(parmap$log_std_pe)));



  if(no.logits){tp$logit_ar_index_age <- NULL;
  tp$logit_ar_logF <- NULL;
  tp$logit_ar_logRec  <- NULL;}

  if(no.Flogits){tp$logit_ar_logF <- NULL;}

  if(no.pe){tp$log_std_pe<-NULL;}

  upper = unlist(tp);

  ret = list(parameters=parameters,parmL=parmL,parmU=parmU,lower=lower,upper=upper)

  return(ret)


}

