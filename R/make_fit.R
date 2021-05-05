#' Prepare parameters for APAM
#'
#' Prepares parameters for use in APAM
#'@useDynLib APAM
#'
#' @importFrom rlang .data
#' @param data input from make_tmb_data
#' @param map input from make_map
#' @param parameters input from make_parameters
#' @param do.sd get sdreport
#' @param do.Zresid get Z residuals
#' @param do.NS  run extra Newton steps (default 4)
#' @param do.hess calculate hessian
#' @param NS default 4 extra Newton step
#' @export
#'
make.fit = function(data,map,parameters,do.sd=TRUE,do.Zresid=TRUE,do.NS=TRUE,do.hess=FALSE,NS=4){

  #fit model
  obj <- TMB::MakeADFun(data,parameters$parms, map=map,
                   random=c("log_F_devt","log_Nt"),
                   DLL = "APAM",
                   control = list(trace=10,eval.max=2000,iter.max=1000), silent = TRUE)

  opt1 <- stats::nlminb(obj$par,obj$fn,obj$gr,
                 upper=parameters$upper,lower=parameters$lower,
                 control = list(trace=0,eval.max=2000,iter.max=1000))


  opt <- stats::nlminb(opt1$par,obj$fn,obj$gr,
                upper=parameters$upper,lower=parameters$lower,
                control = list(trace=0,eval.max=2000,iter.max=1000))

  #take a few extra newton steps
  if(do.NS){ newtonsteps <- NS
  for(i in seq_len(newtonsteps)) {
    g <- as.numeric( obj$gr(opt$par) )
    h <- stats::optimHess(opt$par, obj$fn, obj$gr)
    opt$par <- opt$par - solve(h, g)
    opt$objective <- obj$fn(opt$par)
  }}

  #get model output
  rep <-  obj$report()

  sd.rep <- NULL
  sd.rep1 <- NULL

  if(do.sd){sd.rep <-  TMB::sdreport(obj)}
  hess <- NULL

  #return hessian
  if(do.hess){hess <- numDeriv::jacobian(obj$gr,opt$par)}

  #calculate z residuals
  if(do.Zresid){

    sp <- obj$env$parList(opt$par)

    map1 <- list(
      log_R=factor(c(NA,NA)),
      m_q=factor(matrix(NA,nrow=length(data$isurvey1), ncol = data$A-1)),
      log_cv_index=factor(rep(NA,length(unique(data$isd)))),
      log_std_log_R=factor(NA),
      log_F_mean=factor(rep(matrix(NA,nrow=data$Y,ncol=2,byrow=T))),
      log_std_log_F = factor(matrix(NA,nrow=data$Y,ncol=data$A-4,byrow=T)),
      log_std_pe = factor(matrix(NA,nrow=data$Y-1,ncol=data$A-1)),
      log_std_crl = factor(matrix(NA,nrow=data$Y,ncol=data$A-5)),
      logit_ar_index_age = factor(rep(NA,length(unique(data$isurvey1)))),
      logit_ar_logF = factor(rep(NA,3)),
      logit_ar_crl = factor(c(NA,NA)),
      logit_ar_pe = factor(c(NA,NA)),
      logit_ar_logRec = factor(NA),
      log_F_devt=factor(t(matrix(NA,nrow=data$Y,ncol=data$A-4,byrow=T))),
      log_Nt=factor(t(matrix(NA,nrow=data$Y,ncol=data$A,byrow=T))),
      h=factor(NA)
    )

    data$resid=1
    obj1 <- TMB::MakeADFun(data,sp,map=map1,DLL = "APAM", silent = TRUE)

    rep1 <- obj1$report();
    sd.rep1 <- TMB::sdreport(obj1)

    ind <- names(sd.rep1$value)=="resid_index"

    resid.cov.index <- sd.rep1$cov[ind,ind]
    ch.cov <- chol(resid.cov.index)

    rep$index_Zresid <- qr.solve(t(ch.cov), rep$resid_index)
    rep$resid.cov.index <- resid.cov.index

    ind <- names(sd.rep1$value)=="resid_crl"

    resid.cov.crl <- sd.rep1$cov[ind,ind]
    ch.cov <-chol(resid.cov.crl)

    crl_Zresid <- qr.solve(t(ch.cov), as.vector(rep$resid_crl))

    rep$crl_Zresid <-matrix(crl_Zresid,nrow=data$Y,ncol=data$A-5)
    rep$resid.cov.crl <- resid.cov.crl

  }

  data$resid=0

  ret = list(obj=obj,opt=opt,rep=rep,sdrep=sd.rep,sd.rep1=sd.rep1,
             tmb.data=data,map=map,parameters=parameters, hess=hess)

  return(ret)

}
