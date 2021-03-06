% Generated by roxygen2 (4.1.0): do not edit by hand
% Please edit documentation in R/generalFunctions.R
\name{bsp}
\alias{bsp}
\title{Builder for B-splines basis and penalty matrices typically used in DoubleGam and DoubleRobGam}
\usage{
bsp(x, nknots = 15, p = 3, center = TRUE, sm.p = 1, order = 2)
}
\arguments{
\item{x}{data vector}

\item{nknots}{number of knots used in building the B-splines}

\item{p}{polynomial degree of the basis}

\item{center}{logical indicator on whether the design matrix should be centered (and it should be the case when more than one covariate enter the model)}

\item{sm.p}{smoothing parameter used in the penalty matrix build}

\item{order}{penalty matrix order}
}
\description{
Builder for B-splines basis and penalty matrices typically used in DoubleGam and DoubleRobGam
}

