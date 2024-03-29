% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/chromunity.R
\name{c.Covariate}
\alias{c.Covariate}
\title{title}
\usage{
\method{c}{Covariate}(...)
}
\arguments{
\item{...}{A series of Covariates, note all objects must be of type Covariate}
}
\value{
Covariate object that can be passed directly into the Chromunity object constructor that contains all of the Covariate covariates
Passed in the ... param
}
\description{
Override the c operator for covariates so that you can merge them like a vector
}
\author{
Zoran Z. Gajic
}
