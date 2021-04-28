#' Prepare parameters for APAM
#'
#' Prepares parameters for use in APAM
#'@useDynLib APAM

#' @importFrom rlang .data
#' @param data input from make_tmb_data
#' @param parmap input from make_map
#' @param params input from make_parameters
#' @param do.sd get sdreport
#' @param do.Zresid get Z residuals
#' @param do.NS  run extra Newton steps (default 4)
#' @param do.hess calculate hessian
#' @param NS defualt 4 extra Newton step

make.fit = function(data,parmap,params,do.sd=TRUE,do.Zresid=TRUE,do.NS=TRUE,do.hess=FALSE,NS=4){

  #fit model
  obj <- TMB::MakeADFun(data,params$parameters, map=parmap,
                   random=c("log_F_devt","log_Nt"),
                   DLL = "APAM",
                   control = list(trace=10,eval.max=2000,iter.max=1000), silent = TRUE)

  opt1 <- stats::nlminb(obj$par,obj$fn,obj$gr,
                 upper=params$upper,lower=params$lower,
                 control = list(trace=0,eval.max=2000,iter.max=1000))


  opt <- stats::nlminb(opt1$par,obj$fn,obj$gr,
                upper=params$upper,lower=params$lower,
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
  if(do.hess){hess <- numderiv::jacobian(obj$gr,opt$par)}

  ret = list(obj=obj,opt=opt,rep=rep,sdrep=sd.rep,sd.rep1=sd.rep1,
             tmb.data=data,parmap=parmap,params=params, hess=hess)

  return(ret)

}
