#' Prepare data inputs for APAM
#'
#' Prepares data for use in APAM
#'
#' @importFrom rlang .data
#' @param do.retro T/F, turns on/off option to do retros
#' @param retro.year end year for retro.
#' @param M.split T/F, turns on/off M increase for years 1989-1996
#' @param M.matrix can redefine M input matrix;
#' @param C.bounds defines landings upper bounds; default uncertainty(high=2,low = 1.2, moderate=1.5)
#' @param sdL  default landings sd fixed at 0.05;
#' @param d matrix of 0s; can be changed for influence diagnostics

make.tmb.data = function(do.retro=FALSE,retro.year=NULL,M.split=TRUE,M.matrix=NULL,
                         C.bounds=NULL,sdL=NULL,d=NULL){

  end.year<-max(land$Year)
  if(do.retro){end.year<-retro.year}
  year<- 1960:end.year
  age<- 1:15

  indices <- index %>% dplyr::filter(.data$Year<=end.year)
  crl <- crl.mat[1:length(year),]
  landings <- land %>% dplyr::filter(.data$Year<=end.year)
  cw <- cw.mat[1:length(year),]
  mat <- matur.mat[1:length(year),]
  sw <- sw.mat[1:length(year),]

  #set M assumption
  if(is.null(M.matrix)){
    M.matrix <- matrix(0.2,nrow=length(year),ncol=length(age),byrow=T)
    M.matrix[,1:3] <- rep(0.50)
    M.matrix[,4] <- rep(0.30)}

  #M = 0.53 for all ages for years 1989-1996
  if(M.split){M.matrix[30:37,] <- M.matrix[30:37,] + 0.33}

  temp <- unique(indices$surv_year)

  #make tmb list
  tmb_data = list(
    M = M.matrix,
    weight = as.matrix(sw),
    mat = as.matrix(mat),
    midy_weight = as.matrix(cw),
    index = indices$index,
    olandings = landings$Landings/1000,
    iyear = indices$iyear,
    iage = indices$iage,
    isurvey = indices$isurvey,
    isd = indices$isd ,
    is_year = indices$is_year,
    fs = indices$fs,
    A = length(age),
    Y = length(year),
    Ns = length(unique(indices$is_year)),
    NsF = length(temp[substr(temp ,1,4) == 'Fall']),
    NsSpan = length(temp[substr(temp ,1,4) == 'Span']),
    isurvey1 = unlist(tapply(indices$isurvey,indices$is_year,unique)),
    crl = crl
  )

  #give names for surveys
  names(tmb_data$isd) <- indices$surv_age

  #to set landings bounds
  if(is.null(C.bounds)){C.bounds=c(2,1.2,1.5)}
  tmb_data$landings_L <- tmb_data$olandings

  # high uncertainty
  tmb_data$landings_U <- C.bounds[1]*tmb_data$olandings

  # low uncertainty
  ind = (landings$Year>=1977)&(landings$Year<=1982) | (landings$Year>=1994)&(landings$Year<=2010)
  tmb_data$landings_U[ind] <- C.bounds[2]*tmb_data$olandings[ind]

  # moderate uncertainty
  ind = (landings$Year>=1983)&(landings$Year<=1993) | (landings$Year>=2011)
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

  #for local influence diagnostics (0 for regular model)
  if(is.null(d)){d <- as.vector(matrix(0, nrow=tmb_data$Y, ncol=tmb_data$A,byrow=T))}
  tmb_data$d <- d

  #for profile likelihoods (full model if all = 1)
  tmb_data$nll_wt = rep(1,5)

  return(tmb_data)

}

