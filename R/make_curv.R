#' make.curv: calculate curvature
#'
#' Calculates curvature for APAM for individual, age group and year group perturbations for full model fit.
#' Note that this functions currently does not calculate curvature for individual
#' data components.
#'
#' @param mfits object returned from \code{\link{make.fit}}
#' @param LocInf object returned from \code{\link{make.LI}}
#' @param pert  (optional) vector, to manually set which perturbations to run. Runs all by default.
#' @param tol  to set curvature tolerance. Default \code{h=0.10}.
#'
#' @examples
#' \dontrun{
#' curv<- make.curv(mfits,LocInf)
#' }
#' @export

make.curv = function(mfits,LocInf,pert=NULL,tol = NULL){

  tmb_data <- mfits$tmb.data
  par <- mfits$parameters
  par$parms <- mfits$obj$env$parList(mfits$opt$par)
  full<- gfunc(0,mfits$opt$par,mfits)

  Aname = c(1:15)
  Yname = c(1960:2017)
  Y=tmb_data$Y
  A=tmb_data$A
  n<- length(tmb_data$d)

  if(LocInf$type=="age"){if(is.null(pert)){pert<- c(1:A)}
    all=FALSE}
  if(LocInf$type=="year"){if(is.null(pert)){pert<- c(1:Y)}
    all=FALSE}
  if(LocInf$type=="all"){if(is.null(pert)){pert<-c(1:n)}}

  if(is.null(tol)){tol<-0.10}

  curv = vector("list", length(pert))

  for(j in 1:length(pert)){

    if(LocInf$type=="all"){
      tmb_data$d = as.vector(matrix(0, nrow=Y, ncol=A,byrow=T))
      tmb_data$d[pert[j]] = 1}

    if(LocInf$type=="age"){
      tmat <- matrix(0, nrow=Y, ncol=A,byrow=T)
      tmb_data$d = as.vector(matrix(0, nrow=Y, ncol=A,byrow=T))
      tmat[,pert[j]] = 1

      tmb_data$d <- matrix(tmat, nrow=Y, ncol=A,byrow=F)}

    if(LocInf$type=="year"){
      tmat <- matrix(0, nrow=Y, ncol=A,byrow=T)
      tmb_data$d = as.vector(matrix(0, nrow=Y, ncol=A,byrow=T))
      tmat[pert[j],] = 1

      tmb_data$d <- matrix(tmat, nrow=Y, ncol=A,byrow=F)}

    par$parms$h <- 2*tol

    gmax1<-make.fit(tmb_data,mfits$map,par,do.hess=F, do.sd = F,do.Zresid = F, do.NS = F)

    par$parms$h <- tol

    gmax2<-make.fit(tmb_data,mfits$map,par,do.hess=F, do.sd = F,do.Zresid = F, do.NS = F)

    par$parms$h <- -tol

    gmin1 <- make.fit(tmb_data,mfits$map,par,do.hess=F, do.sd = F,do.Zresid = F, do.NS = F)

    par$parms$h <- -2*tol

    gmin2 <- make.fit(tmb_data,mfits$map,par,do.hess=F, do.sd = F,do.Zresid = F, do.NS = F)

    temp1<- -gmax1$obj$fn(gmax1$opt$par) + 16*gmax2$obj$fn(gmax2$opt$par)  - 30*full + 16*gmin1$obj$fn(gmin1$opt$par) - gmin2$obj$fn(gmin2$opt$par)

    dg2_dh <- temp1/(12*tol^2)

    temp2<- (1+(as.numeric(LocInf$LI[j,5]))^2)^(3/2)

    tempc <- dg2_dh/temp2
    curv[j]=tempc
  }
  Curvature<- list(curv=curv)

  return(Curvature)
}
