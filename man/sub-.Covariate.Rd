% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/chromunity.R
\name{[.Covariate}
\alias{[.Covariate}
\title{title}
\usage{
\method{[}{Covariate}(obj, range)
}
\arguments{
\item{obj}{Covariate This is the Covariate to be subset}

\item{range}{vector This is the range of Covariates to return, like subsetting a vector. e.g. c(1,2,3,4,5)[3:4] == c(3,4)}
}
\value{
A new Covariate object that contains only the Covs within the given range
}
\description{
Overrides the subset operator x[] for use with Covariate to allow for vector like subsetting
}
\author{
Zoran Z. Gajic
}
