#' make.profile: calculate profile likelihoods
#'
#' Calculates profile likelihoods for APAM for individual data components.
#'@useDynLib APAM
#'
#' @param mfits object returned from \code{\link{make.fit}}
#' @param delta_M (optional) vector, to manually set range of M perturbations. Default \code{delta_M <- seq(-0.15,0.40, by = 0.01)}.
#'
#'@return Returns marginal profile for age composition data and remainder/conditional for
#' surveys and landings data. Additionally, components from the joint negative log-likelihood
#' are also returned.
#'   \describe{
#'     \item{\code{$conditional}}{conditional profile likelihoods for surveys and landings data}
#'     \item{\code{$remainder}}{remainder profile likelihoods for surveys and landings data}
#'     \item{\code{$marginal}}{marginal profile likelihoods for age composition data}
#'     \item{\code{$conditional_dev}}{conditional profile likelihood deviations for surveys and landings data}
#'     \item{\code{$remainder_dev}}{remainder profile likelihood deviations for surveys and landings data}
#'     \item{\code{$marginal_dev}}{marginal profile likelihood deviationsfor age composition data}
#'     \item{\code{$joint}}{joint profile likelihoods for all data components}
#'     \item{\code{$joint_dev}}{joint profile likelihood deviations for all data components}
#'     \item{\code{$conv_fit}}{returns max inner and outer gradients for M perturbed model fit and convergence message from optimization}
#'     \item{\code{$no.gr2}}{returns perturbation number and error message if error is returned during netwon steps}
#'     }
#'
#' @examples
#' \dontrun{
#' profile <- make.profile(mfits)
#' }
#' @export

make.profile = function(mfits, delta_M = NULL){

  tmb_data <- mfits$tmb.data
  tmb_data.orig <- tmb_data
  params.orig <- mfits$parameters

  if(is.null(delta_M)){delta_M <- seq(-0.15,0.40, by = 0.01)}

  ndel = length(delta_M)
  no.jnll = length(mfits$rep$jnll)

  prof_fit <- vector("list", ndel)
  conv_fit <- matrix(NA, nrow=3, ncol=ndel)
  rownames(conv_fit) <- c("Inner max grad","Outer max grad", "Fit conv")
  colnames(conv_fit) <- as.character(delta_M)

  no.gr1=NULL
  no.gr2=NULL
  no.gr3=NULL

  joint_res = matrix(NA, nrow = 8, ncol=ndel)

  for(n in 1:ndel){

    NS=0
    parm_F = 0.1
    count=1

    tmb_data$M <- tmb_data.orig$M + delta_M[n]
    params <- params.orig
    params$parms$log_F_mean  <- matrix(log(0.1-(parm_F*delta_M[n])),nrow=tmb_data$Y,ncol=2,byrow=T)

    temp_obj <- TMB::MakeADFun(tmb_data,params$parms, map=mfits$map, random=mfits$obj$env$random,
                               DLL = mfits$obj$env$DLL,control = list(trace=10,eval.max=2000,iter.max=1000),inner.control = list(maxit = 1000),silent = T)

    temp1<-temp_obj$fn(temp_obj$par)
    gr1<-max(abs(temp_obj$env$f(temp_obj$env$last.par, order=1)[temp_obj$env$random]))

    #grad inner optimization
    while(gr1>1e-8){

      temp_obj <- TMB::MakeADFun(tmb_data,parameters=temp_obj$env$parList(), map=mfits$map,
                                 random=mfits$obj$env$random, DLL = mfits$obj$env$DLL,
                                 control = list(trace=10,eval.max=2000,iter.max=1000), inner.control = list(maxit = 1000),silent = T)

      temp<-temp_obj$fn(temp_obj$par)
      gr1<-max(abs(temp_obj$env$f(temp_obj$env$last.par, order=1)[temp_obj$env$random]))
      count=count+1

      if(count==300){
        no.gr1=rbind(no.gr1,paste(n))
      }}

    temp_opt <-stats::nlminb(temp_obj$par,temp_obj$fn,temp_obj$gr,
                             upper=params$upper,lower=params$lower,
                             control = list(trace=0,eval.max=2000,iter.max=1000))

    gr2<-abs(max(temp_obj$gr(temp_opt$par)))

    #grad outer optimization
    # Take a few extra newton steps
    tryCatch(for(i in 1:4) {
      g <- as.numeric(temp_obj$gr(temp_opt$par))
      h <- stats::optimHess(temp_opt$par, temp_obj$fn, temp_obj$gr)
      temp_opt$par <- temp_opt$par - solve(h, g)
      temp_opt$objective <- temp_obj$fn(temp_opt$par)
    },error=function(e){no.gr2 <<- rbind(paste(n,conditionMessage(e),sep=":"),no.gr2)})

    gr2<-abs(max(temp_obj$gr(temp_opt$par)))

    conv_fit[1,n] <- gr1
    conv_fit[2,n] <- gr2
    conv_fit[3,n] <- temp_opt$message
    prof_fit[[n]] <- list(obj=temp_obj,
                          opt=temp_opt
    )
    rep=temp_obj$report()
    joint_res[,n] = rep$jnll
  }

  remain_prof = matrix(NA,nrow=5,ncol=ndel)
  cond_prof = matrix(NA,nrow=5,ncol=ndel)
  marg_prof = matrix(NA,nrow=1,ncol=ndel)

  for(j in 1:ndel){

    NS=1
    temp <- prof_fit[[j]]
    dat <- temp$obj$env$data
    nllr <- matrix(NA,1,5)
    nllc <- matrix(NA,1,5)
    nllm <-matrix(NA,1,1)

    for(k in 5:1){

      count=1
      dat$nll_wt[1:5] = 0

      if(k==5){
        dat$nll_wt[k] = 1
      }else if(k==4){
        dat$nll_wt[c(1:3,5)] = 1
      }else if(k==3){
        dat$nll_wt[c(1:2,4:5)] = 1
      }else if(k==2){
        dat$nll_wt[c(1,3:5)] = 1
      }else{dat$nll_wt[c(2:5)] = 1}

      obj1 <- TMB::MakeADFun(dat, parameters=temp$obj$env$parList(temp$opt$par),map=temp$obj$env$map,random=temp$obj$env$random,
                             DLL = temp$obj$env$DLL, silent = TRUE, inner.control = list(maxit = 1000))

      t2<-obj1$fn(obj1$par)
      gr3<-max(abs(obj1$env$f(obj1$env$last.par, order=1)[obj1$env$random]))

      while(gr3>1e-8){

        obj1 <- TMB::MakeADFun(dat,parameters=obj1$env$parList(), map=temp$obj$env$map,random=temp$obj$env$random,
                               DLL = temp$obj$env$DLL, silent = TRUE, inner.control = list(maxit = 1000))

        t2<- obj1$fn(obj1$par)
        gr3<-max(abs(obj1$env$f(obj1$env$last.par, order=1)[obj1$env$random]))
        count=count+1

        if(count==300){
          no.gr3=rbind(no.gr3,paste(j,k))
        }}

      opt1 <- stats::nlminb(obj1$par,obj1$fn,obj1$gr,
                            upper=params$upper,lower=params$lower,
                            control = list(trace=0,eval.max=2000,iter.max=1000))

      if(k<5){nllc[1,k] = temp$opt$objective - obj1$fn(obj1$par)
      nllr[1,k] =  temp$opt$objective - obj1$fn(opt1$par)
      }else{nllm<-obj1$fn(obj1$par)}


    }

    nllr[1,5] =  temp$opt$objective
    nllc[1,5] = temp$opt$objective

    remain_prof[,j] = nllr
    cond_prof[,j] = nllc
    marg_prof[,j] = nllm

  }

  rownames(remain_prof) = c("Fall index","Spring index","Spanish index","Landings","Total")
  colnames(remain_prof) = as.character(delta_M)

  rownames(cond_prof) = c("Fall index","Spring index","Spanish index","Landings","Total")
  colnames(cond_prof) = as.character(delta_M)

  rownames(marg_prof) = "Age comps"
  colnames(marg_prof) = as.character(delta_M)

  rownames(joint_res) = c("Fall index","Spring index","Spanish index","Landings","Age comps", "Recs","F", "PE")
  colnames(joint_res) = as.character(delta_M)

  remain_prof_dev = t(apply(remain_prof,1,function(x){x-min(x,na.rm=T)}))
  cond_prof_dev = t(apply(cond_prof,1,function(x){x-min(x,na.rm=T)}))
  marg_prof_dev = marg_prof-min(marg_prof)
  joint_prof_dev =  t(apply(joint_res,1,function(x){x-min(x,na.rm=T)}))

  results = list(conditional = cond_prof, remainder = remain_prof,marginal = marg_prof,
                 conditional_dev = cond_prof_dev,remainder_dev = remain_prof_dev, marginal_dev = marg_prof_dev,
                 joint = joint_res, joint_dev = joint_prof_dev,
                 conv_fit = t(conv_fit), no.gr2=no.gr2)
  return(results)
}
