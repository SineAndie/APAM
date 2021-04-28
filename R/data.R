#' 3LNO American plaice survey indices
#'
#' A dataset containing the survey data for 3LNO American plaice
#'
#'
#' @format A data frame with 840 observations and 12 variables
#' \describe{
#'   \item{Year}{survey year}
#'   \item{Age}{survey age in years}
#'   \item{index}{survey abundance index in millions}
#'   \item{survey}{Type of survey; Fall/Spring}
#'   \item{fs}{Fraction of year that survey takes places}
#'   \item{iyear}{year indicator}
#'   \item{iage}{age indicator}
#'   \item{isurvey}{survey indicator}
#'   \item{surv_age}{identifier survey and age combination}
#'   \item{isd}{survey/age identifier}
#'   \item{surv_year}{identifier survey and year combination}
#'   \item{is_year}{survey/year identifier}
#' }
"index"

#' Stock weights
#'
#' Stock weights estimated from from a spatiotemporal biphasic vonBertalanffy
#' growth model. accounted for the length-stratified age sampling design.
#' The 3LNO stock weights were combined for each division
#' by weighting the values for each division by the average
#' abundance index at age during 1975–2017. Stock weights
#' prior to 1975 were fixed at the mean values for 1975–77.
#'
#' @format A data frame with 58 observations and 15 variables; rows represent
#' years, columns represent ages and values are weight in kg.
"sw.mat"

#' Continuation ratio logit transformation
#'
#' The age composition data are fit using the continuation ratio
#' logit (crl) transformation (see Perreault et al., 2020) for details.
#'
#' @format A data frame with 58 observations and 11 variables; rows represent
#' years, columns represent ages and values are crl at that age and year
"crl.mat"

#' Commercial landings
#'
#' Nominal landings (t) of American plaice for NAFO divisions 3LNO
#'
#' @format A data frame with 58 observations and 2 variables
#'  \describe{
#'   \item{Landings}{landings estimates (t)}
#'   \item{Year}{landed year}
#'
#'   }
"land"

#' Maturity at age
#'
#' Estimated proportion mature for 3LNO American plaice
#'
#' @format A data frame with 59 observations and 12 variables; rows represent
#' years, columns represent ages and values are proportion mature
"matur.mat"

#' Maturity at age
#'
#' Estimated commercial catch weights for 3LNO American plaice
#'
#' @format A data frame with 58 observations and 12 variables; rows represent
#' years, columns represent ages and values are weight in kg.
"cw.mat"


