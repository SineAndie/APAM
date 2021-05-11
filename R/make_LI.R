#' dFgrad_dtheta: helper function for local influence diagnostics
#'
#' Helper function function to compute Del in Eq. 21 of Perreault & Cadigan (2021).
#' @useDynLib APAM
#'
#' @param w set weight for perturbation (0 for null perturbation).
#' @param i indicator for observation weight \code{1=on, 0=off}.
#' @param mfits object returned from \code{\link{make.fit}}
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


#' dFgrad_dtheta_age: helper function for local influence diagnostics
#'
#' Helper function function to compute Del in Eq. 21 of Perreault & Cadigan (2021).
#' @param w set weight for perturbation (0 for null perturbation).
#' @param i indicator for observation weight \code{1=on, 0=off}.
#' @param mfits object returned from \code{\link{make.fit}}
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

#' dFgrad_dtheta_year: helper function for local influence diagnostics
#'
#' Helper function function to compute Del in Eq. 21 of Perreault & Cadigan (2021)
#' @param w set weight for perturbation (0 for null perturbation).
#' @param i indicator for observation weight \code{1=on, 0=off}.
#' @param mfits object returned from \code{\link{make.fit}}
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

#'gfunc: helper function to return model fit
#'
#' Influence function for local influence diagnostics.
#' @param w set weight for perturbation (0 for null perturbation).
#' @param theta vector of fixed effect parameters
#' @param mfits object returned from \code{\link{make.fit}}
#'
gfunc = function(w, theta,mfits){
  param <- mfits$obj$env$parList(theta)
  param$h <- w

  obj <- TMB::MakeADFun(mfits$tmb.data,param,mfits$map,
                        random=mfits$obj$env$random, DLL = mfits$obj$env$DLL,
                        inner.control=list(maxit=100,trace=F),silent=T)
  return(obj$fn(obj$par))
}

#'make.LI: run local influence diagnostics for APAM.
#'
#' To run local influence diagnostics for full model and individual data components.
#' @param mfits object returned from \code{\link{make.fit}}
#' @param pert (optional) vector, to manually set which perturbations to run. Runs all by default. See details.
#' @param all  T/F, to run individual (i.e. age and year) perturbations. Default = \code{TRUE}.
#' @param age T/F, to run age group perturbations. Default = \code{FALSE}.
#' @param year T/F, to run year group perturbations. Default = \code{FALSE}.
#'
#'
#'@details
#'   \describe{
#'     \item{\code{pert}}{For all age and year group perturbations, pert is a vector of length (nyears*nages). Default \code{pert = c(1:870)}.
#'     For age group perturbations, pert has length = nages, and for year group perturbations, pert has length = nyears.}}
#'
#' @export
#' @examples
#' \dontrun{
#' LI <- make.LI(mfits)
#'
#' #to run age group perturbations
#' LI_age <- make.LI(mfits,age=T)
#' }
#' @export

make.LI = function(mfits,pert=NULL,all=TRUE,age=FALSE,year=FALSE){

  tmb_data <- mfits$tmb.data

  n.survey.sd <- length(unique(tmb_data$isd))
  Y<-tmb_data$Y
  A <- tmb_data$A
  unpert <- gfunc(0,mfits$opt$par,mfits)

  if(age){if(is.null(pert)){pert <- c(1:A)}
    all=FALSE
    type="age"}
  if(year){if(is.null(pert)){pert <- c(1:Y)}
    all=FALSE
    type='year'}
  if(all){if(is.null(pert)){pert <-c(1:length(tmb_data$d))}
    type="all"}

  LI_temp<-matrix(NA, nrow = length(pert), ncol = (length(tmb_data$nll_wt)))

  LI<-matrix(NA, nrow = length(pert), ncol = (length(tmb_data$nll_wt)))
  colnames(LI) <- c("Fall", "Spring","Landings","Age Comps","Full")
  rownames(LI) <- pert

  if(all){

    dFtimeFn = function(num){
      return(numDeriv::jacobian(dFgrad_dtheta,0,,,,num,mfits))}

    #containers
    Si = vector("list", length(pert))

    for(j in 1:length(pert)){

      tmb_data$d <- as.vector(matrix(0, nrow=Y, ncol=A,byrow=T))
      tmb_data$d[pert[j]] <- 1

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
        resid_crl_res = factor(matrix(NA, nrow = Y, ncol = A-5))
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

    return(results = list(LI= Si, type = type))
  }else if(age){

    dFtimeFn = function(num){
      return(numDeriv::jacobian(dFgrad_dtheta_age,0,,,,num,mfits))
    }

    for(i in 1:length(tmb_data$nll_wt)){

      tmb_data$nll_wt = rep(1,5)

      if(i<3){tmb_data$nll_wt[i]=0}

      if(i==3){tmb_data$nll_wt[i+1]=0}

      if(i==4){tmb_data$nll_wt = rep(0,5)
      tmb_data$nll_wt[i+1]=1}

      #containers
      Si = vector("list", length(pert))

      for(j in 1:length(pert)){

        tmat <- matrix(0, nrow=Y, ncol=A,byrow=T)
        tmb_data$d <- as.vector(matrix(0, nrow=Y, ncol=A,byrow=T))
        tmat[,pert[j]] <- 1

        tmb_data$d <- matrix(tmat, nrow=Y, ncol=A,byrow=F)

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
          resid_crl_res = factor(matrix(NA, nrow = Y, ncol = A-5))
        )

        obj2 <- TMB::MakeADFun(tmb_data, mfits$obj$env$parList(mfits$opt$par),map2,
                               random=c("log_F_devt","log_Nt"), DLL = mfits$obj$env$DLL,silent=T)

        dg_dw <- as.matrix(obj2$gr())

        #dg_dtheta
        obj3 <- TMB::MakeADFun(tmb_data, mfits$obj$env$parList(mfits$opt$par),mfits$map,
                               random=c("log_F_devt","log_Nt"), DLL = mfits$obj$env$DLL,silent=T)

        dg_dtheta <- as.matrix(obj3$gr())

        #results
        Si[[j]] <- dg_dw - Del%*%solve(mfits$hess)%*%dg_dtheta

      }
      LI_temp[,i] <- unlist(Si)
    }
    LI[,1:3] = LI_temp[,5] - LI_temp[,1:3]
    LI[,4:5] = LI_temp[,4:5]
    return(results = list(LI=LI, type = type))}else{
      for(i in 1:length(tmb_data$nll_wt)){

        dFtimeFn = function(num){
          return(numDeriv::jacobian(dFgrad_dtheta_year,0,,,,num,mfits))
        }

        tmb_data$nll_wt = rep(1,5)

        if(i<3){tmb_data$nll_wt[i]=0}

        if(i==3){ tmb_data$nll_wt[i+1]=0}

        if(i==4){tmb_data$nll_wt = rep(0,5)
        tmb_data$nll_wt[i+1]=1}

        #containers
        Si = vector("list", length(pert))

        for(j in 1:length(pert)){

          tmat <- matrix(0, nrow=Y, ncol=A,byrow=T)
          tmb_data$d <- as.vector(matrix(0, nrow=Y, ncol=A,byrow=T))
          tmat[pert[j],] <- 1

          tmb_data$d <- matrix(tmat, nrow=Y, ncol=A,byrow=F)

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
            resid_crl_res = factor(matrix(NA, nrow = Y, ncol = A-5))
          )

          obj2 <- TMB::MakeADFun(tmb_data, mfits$obj$env$parList(mfits$opt$par),map2,
                                 random=c("log_F_devt","log_Nt"), DLL = mfits$obj$env$DLL,silent=T)

          dg_dw <- as.matrix(obj2$gr())

          #dg_dtheta
          obj3 <- TMB::MakeADFun(tmb_data, mfits$obj$env$parList(mfits$opt$par),mfits$map,
                                 random=c("log_F_devt","log_Nt"), DLL = mfits$obj$env$DLL,silent=T)

          dg_dtheta <- as.matrix(obj3$gr())

          #results
          Si[[j]] <- dg_dw - Del%*%solve(mfits$hess)%*%dg_dtheta

        }
        LI_temp[,i] <- unlist(Si)
      }
      LI[,1:3] = LI_temp[,5] - LI_temp[,1:3]
      LI[,4:5] = LI_temp[,4:5]
      return(results = list(LI=LI, type = type))}
}
