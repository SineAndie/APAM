#' make.tmb.data: prepare data inputs for APAM
#'
#' Prepares data for use in APAM. No inputs need to be changed to run the default model (see details).
#'
#' @importFrom rlang .data
#' @param do.retro (optional) T/F, turn on/off option to change retro year? Default = \code{F}.
#' @param retro.year  (optional) numeric, defines end year for retros if needed. Default = \code{NULL}.
#' @param M.split  (optional) T/F, turn on/off M increase for years 1989-1996? Default = \code{T}.
#' @param M.matrix (optional) a matrix to manually define M assumption. Default = \code{NULL} See details. .
#' @param C.bounds (optional) a vector to manually define landings upper bounds; Default = \code{NULL}. See details.
#' @param sdL  (optional) a positive scalar to  manually defined landings sd ; Default = \code{NULL}.
#' @param d (optional) a vector to turn on/off d for use in \code{make.LI()}. Default = \code{NULL}.
#'
#'
#'@details
#'   \describe{
#'     \item{\code{M.matrix}}{M.matrix has dimensions nrow=nyears and ncol=nages; M = 0.5 for ages 1-3, 0.3 for age 4 and 0.2 for ages 5+. }
#'     \item{\code{M.split}}{if \code{TRUE}, \code{M.matrix = M.matrix + 0.33} for years 1989-1996.}
#'     \item{\code{C.bounds}}{vector of length 3; Default \code{C.bounds = (2,1.2,1.5)} for high, low, and moderate uncertainty.}
#'     \item{\code{sdL}}{default sd landings = 0.05}
#'     \item{\code{d}}{Numeric vector of length(nyears*nages) of 0's. }}
#'
#'
#' @examples
#' \dontrun{
#' tmb.data <- make.tmb.data()
#'
#' #model with no M split
#' tmb.data.ns = make.tmb.data(M.split=FALSE)
#' }
#' @export
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
  tmb_data$resid=0

  return(tmb_data)

}

