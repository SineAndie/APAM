% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/make_LI.R
\name{dFgrad_dtheta}
\alias{dFgrad_dtheta}
\title{dFgrad_dtheta: helper function for local influence diagnostics}
\usage{
dFgrad_dtheta(w, i, mfits)
}
\arguments{
\item{w}{set weight for perturbation (0 for null perturbation).}

\item{i}{indicator for observation weight \code{1=on, 0=off}.}

\item{mfits}{object returned from \code{\link{make.fit}}}
}
\description{
Helper function function to compute Del in Eq. 21 of Perreault & Cadigan (2021).
}
