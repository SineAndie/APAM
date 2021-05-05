#' Helper functions for local influence diagnostics
#'
#' #function to compute Del in Eq. 22 of Perreault & Cadigan
#'#'@useDynLib APAM
#'
#' @param w weight for perturbation (set to 0 for null perturbation)
#' @param i inidcator for observation weight
#' @param mfits object from make_fit
#'
dFgrad_dtheta = function(w,i,mfits){
  mfits$tmb.data$d[i] <- 1
  param <- mfits$obj$env$parList(mfits$opt$par)
  param$h <- w

  obj <- TMB::MakeADFun(mfits$tmb.data,param,mfits$map,
                        random=mfits$obj$env$random, DLL = mfits$obj$env$DLL,
                        inner.control=list(maxit=100,trace=F),silent=T)
  return(obj$gr(obj$par))
}


#Del for slopes across ages
#'
#' @param w weight for perturbation (set to 0 for null perturbation)
#' @param i inidcator for observation weight
#' @param mfits object from make_fit
#'
dFgrad_dtheta_age = function(w,i,mfits){
  param <- mfits$obj$env$parList(mfits$opt$par)
  tmat <- matrix(0, nrow=mfits$tmb.data$Y, ncol=mfits$tmb.data$A,byrow=T)
  tmat[,i] <- 1

  mfits$tmb.data$d <- matrix(tmat, nrow=mfits$tmb.data$Y, ncol=mfits$tmb.data$A,byrow=F)
  param$h <- w

  obj <- TMB::MakeADFun(mfits$tmb.data,param,mfits$map,
                        random=mfits$obj$env$random, DLL = mfits$obj$env$DLL,
                        inner.control=list(maxit=100,trace=F),silent=T)
  return(obj$gr(obj$par))
}

#Del for slopes across years
#'
#' @param w weight for perturbation (set to 0 for null perturbation)
#' @param i inidcator for observation weight
#' @param mfits object from make_fit
#'
dFgrad_dtheta_year = function(w,i,mfits){
  param <-  mfits$obj$env$parList(mfits$opt$par)
  tmat <- matrix(0, nrow=mfits$tmb.data$Y, ncol=mfits$tmb.data$A,byrow=T)
  tmat[i,] <- 1

  mfits$tmb.data$d <- matrix(tmat, nrow=mfits$tmb.data$Y, ncol=mfits$tmb.data$A,byrow=F)
  param$h <- w

  obj <- TMB::MakeADFun(mfits$tmb.data,param,mfits$map,
                        random=mfits$obj$env$random, DLL = mfits$obj$env$DLL,
                        inner.control=list(maxit=100,trace=F),silent=T)
  return(obj$gr(obj$par))
}

#'To calculate influence function
#'
#' @param w weight for perturbation (set to 0 for null perturbation)
#' @param theta vector of fixed effects
#' @param mfits object from make_fit
#'
#influence function
gfunc = function(w, theta,mfits){
  param <- mfits$obj$env$parList(theta)
  param$h <- w

  obj <- TMB::MakeADFun(mfits$tmb.data,param,mfits$parmap,
                        random=mfits$obj$env$random, DLL = mfits$obj$env$DLL,
                        inner.control=list(maxit=100,trace=F),silent=T)
  return(obj$fn(obj$par))
}
#' To run local influence diagnostics for individual data components
#'
#' @param mfits object from make_fit
#' @param pert define which perturbations to run; default=all (for all ages/years)
#' @param full if running for data source, need to specify full model fit for Del calculation
#' @param all  to run individual perturbations
#' @param age to run age group perturbations
#' @param year tp run year group perturbations
#' @export

make.LI = function(mfits,pert=NULL,full=NULL,all=TRUE,age=FALSE,year=FALSE){

  full<-mfits
  tmb_data <- mfits$tmb.data

  n <- length(tmb_data$d)
  n.survey.sd <- length(unique(tmb_data$isd))
  Y<-tmb_data$Y
  A <- tmb_data$A
  unpert <- gfunc(0,mfits$opt$par,mfits)
  t1 <- (length(unique(levels(mfits$parmap$m_q)))-8)/2

  if(age){if(is.null(pert)){pert <- c(1:A)}
    all=FALSE}
  if(year){if(is.null(pert)){pert <- c(1:Y)}
    all=FALSE}
  if(all){if(is.null(pert)){pert <-c(1:n)} else{LI.plot=F}}

  LI_temp<-matrix(NA, nrow = length(pert), ncol = (length(tmb_data$nll_wt)+1))
  colnames(LI_temp) <- c("noFall", "noSpan", "noSpring","noland","agecomp","full")

  LI<-matrix(NA, nrow = length(pert), ncol = (length(tmb_data$nll_wt)+1))
  colnames(LI) <- c("Fall", "Span", "Spring","land","agecomp","full")

  for(i in 1:(length(mfits$tmb.data$nll_wt)+1)){

    if(i<5){  tmb_data$nll_wt[i]=0}
    if(i==5){ tmb_data$nll_wt = rep(0,5)
    tmb_data$nll_wt[i]=1}
    if(i==6){ tmb_data$nll_wt = rep(1,5)}

  #containers
  Si = vector("list", length(pert))
  for(j in 1:length(pert)){

    if(all){
      tmb_data$d <- as.vector(matrix(0, nrow=Y, ncol=A,byrow=T))
      tmb_data$d[pert[j]] <- 1

      if(is.null(full)){full<-mfits}
      dFtimeFn = function(num){
        return(numDeriv::jacobian(dFgrad_dtheta,0,,,,num,full))
      }}

    if(age){
      tmat <- matrix(0, nrow=Y, ncol=A,byrow=T)
      tmb_data$d <- as.vector(matrix(0, nrow=Y, ncol=A,byrow=T))
      tmat[,pert[j]] <- 1

      tmb_data$d <- matrix(tmat, nrow=Y, ncol=A,byrow=F)

      if(is.null(full)){full<-mfits}
      dFtimeFn = function(num){
        return(numDeriv::jacobian(dFgrad_dtheta_age,0,,,,num,full))
      }}

    if(year){
      tmat <- matrix(0, nrow=Y, ncol=A,byrow=T)
      tmb_data$d <- as.vector(matrix(0, nrow=Y, ncol=A,byrow=T))
      tmat[pert[j],] <- 1

      tmb_data$d <- matrix(tmat, nrow=Y, ncol=A,byrow=F)

      if(is.null(full)){full<-mfits}
      dFtimeFn = function(num){
        return(numDeriv::jacobian(dFgrad_dtheta_year,0,,,,num,full))
      }}

    Del <- t(sapply(pert[j],dFtimeFn))

    #dg_dw
    map2 <- list(
      log_R=factor(c(NA,NA)),
      m_q=factor(matrix(NA,nrow=length(tmb_data$isurvey1), ncol = 14,byrow=TRUE)),
      log_cv_index=factor(rep(NA,n.survey.sd)),
      log_std_log_R=factor(NA),
      log_F_mean = factor(matrix(NA,nrow=Y,ncol=2,byrow=T)),
      log_std_log_F = factor(matrix(NA,nrow=Y,ncol=A-4,byrow=T)),
      log_std_pe = factor(matrix(NA,nrow=Y-1,ncol=A-1)),
      log_std_crl = factor(matrix(NA,nrow=Y,ncol=A-5)),
      logit_ar_index_age = factor(rep(NA,length(unique(tmb_data$isurvey)))),
      logit_ar_logF = factor(c(NA, NA, NA)),
      logit_ar_crl = factor(c(NA,NA)),
      logit_ar_pe = factor(rep(NA,2)),
      logit_ar_logRec = factor(NA),
      resid_index_res = factor(rep(NA, length(tmb_data$index))),
      resid_crl_res = factor(matrix(NA, nrow = 58, ncol = 10))
    )

    obj2 <- TMB::MakeADFun(tmb_data, mfits$obj$env$parList(mfits$opt$par),map2,
                      random=c("log_F_devt","log_Nt"), DLL = mfits$obj$env$DLL,silent=T)

    dg_dw <- as.matrix(obj2$gr())

    #dg_dtheta
    obj3 <- TMB::MakeADFun(tmb_data, mfits$obj$env$parList(mfits$opt$par),mfits$map,
                      random=c("log_F_devt","log_Nt"), DLL = mfits$obj$env$DLL,silent=T)

    dg_dtheta <- as.matrix(obj3$gr())

    #results
    temp_Si <- dg_dw - Del%*%solve(mfits$hess)%*%dg_dtheta
    Si[[j]] <- temp_Si

  }
  LI_temp[,i] <- unlist(Si)
}
  LI[,1:4] = LI_temp[,6] - LI_temp[,1:4]
  LI[,5:6] = LI_temp[,5:6]
  return(LI)

}


