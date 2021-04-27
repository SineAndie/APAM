#' Prepare data inputs for APAM
#'
#' Here is a nice explanation of what this does
#'
#'
#' @param do.retro T/F, turns on/off option to do retros
#' @param retro.year end year for retro.

make.tmb.data = function(do.retro=FALSE,retro.year=NULL,M.split=TRUE,C.bounds=NULL,M.matrix=NULL,sdL=NULL,d=NULL,SS=T){

  end.year<-max(landings$Year)
  if(do.retro){end.year<-retro.year}
  year<- 1960:end.year
  age<- 1:15

  ind.in <- indices %>% dplyr::filter(Year<=end.year)
  crl.in <- crl[1:length(year),]
  landings.in <- landings %>% dplyr::filter(Year<=end.year)
  catch.ind <- catch.wt[1:length(year),]
  mat.in <- mat.all[1:length(year),]
  sw.ind <- mean.sw.mat[1:length(year),]

  ######################
  #set M assumption
  if(is.null(M.matrix)){
    M.matrix <- matrix(0.2,nrow=length(year),ncol=length(age),byrow=T)
    M.matrix[,1:3] <- rep(0.50)
    M.matrix[,4] <- rep(0.30)}

  #M = 0.53 for all ages for years 1989-1996
  if(M.split){M.matrix[30:37,] <- M.matrix[30:37,] + 0.33}

  temp <- unique(ind.in$surv_year)

  #make tmb list
  tmb_data = list(
    M = M.matrix,
    weight = as.matrix(sw.ind),
    mat = as.matrix(mat.in),
    midy_weight = as.matrix(catch.ind),
    index = ind.in$index,
    olandings = landings.in$Landings/1000,
    iyear = ind.in$iyear,
    iage = ind.in$iage,
    isurvey = ind.in$isurvey,
    isd = ind.in$isd ,
    is_year = ind.in$is_year,
    fs = ind.in$fs,
    A = length(age),
    Y = length(year),
    Ns = length(unique(ind.in$is_year)),
    NsF = length(temp[substr(temp ,1,4) == 'Fall']),
    NsSpan = length(temp[substr(temp ,1,4) == 'Span']),
    isurvey1 = unlist(tapply(ind.in$isurvey,ind.in$is_year,unique)),
    crl = crl.in
  )

  #give names for surveys
  names(tmb_data$isd) <- indices$surv_age

  #to set landings bounds
  if(is.null(C.bounds)){C.bounds=c(2,1.2,1.5)}
  tmb_data$landings_L <- tmb_data$olandings

  # high uncertainty
  tmb_data$landings_U <- C.bounds[1]*tmb_data$olandings

  # low uncertainty
  ind = (landings.in$Year>=1977)&(landings.in$Year<=1982) | (landings.in$Year>=1994)&(landings.in$Year<=2010)
  tmb_data$landings_U[ind] <- C.bounds[2]*tmb_data$olandings[ind]

  # moderate uncertainty
  ind = (landings.in$Year>=1983)&(landings.in$Year<=1993) | (landings.in$Year>=2011)
  tmb_data$landings_U[ind] <- C.bounds[3]*tmb_data$olandings[ind]

  #format for model
  tmb_data$log_landings <- rep(0,length(log(tmb_data$olandings)))
  tmb_data$log_landings_L <- log(tmb_data$landings_L)
  tmb_data$log_landings_U <- log(tmb_data$landings_U)
  tmb_data$log_lowerM <- log(tmb_data$landings_L)
  tmb_data$log_upperM <- log(tmb_data$landings_U)

  #sd for landings
  if(is.null(sdL)){sdL=0.05}
  tmb_data$std_log_landings <- sdL

  #for local influence diagonstics (0 for regular model)
  if(is.null(d)){d <- as.vector(matrix(0, nrow=tmb_data$Y, ncol=tmb_data$A,byrow=T))}
  tmb_data$d <- d

  #for profile likelihoods (full model if all = 1)
  tmb_data$nll_wt = rep(1,5)

  return(tmb_data)

}

