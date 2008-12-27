\name{CADFpvalues}
\alias{CADFpvalues}
\title{p-values of the CADF test.}
\description{
  The asymptotic p-values of the Hansen's (1995) Covariate-Augmented Dickey Fuller (CADF) test 
  for a unit root are computed using the approach outlined in Costantini, Lupi, and Popp (2007).
}
\usage{
CADFpvalues(t0, rho2 = 0.5, trend=c("c", "nc", "ct", "none", "drift", "trend"))
}
\arguments{
  \item{t0}{the value of the test statistic;}
  \item{rho2}{the estimated value of the long-run correlation;}
  \item{trend}{defines the deterministic kernel used in the test. It accepts the values used either in 
               \code{fUnitRoots} or \code{urca} packages. It specifies if the underlying model must be with 
               constant ("c" or "drift", the default), without constant ("nc" or "none"), or with constant and 
               trend ("ct" or "trend").}
}
\value{
  \code{p.value}, a scalar containing the estimated asymptotic p-value of the test.
}
\references{
    Hansen, BE (1995): "Rethinking the Univariate Approach to Unit Root Testing: Using
	Covariates to Increase Power", \emph{Econometric Theory}, \bold{11} (5), 1148--1171.

    Costantini M, Lupi C, Popp S (2007), "A Panel-CADF Test for Unit Roots", University of Molise,
    \emph{Economics \& Statistics Discussion Paper} 39/07, 
    URL \url{http://econpapers.repec.org/paper/molecsdps/esdp07039.htm}.
}
\author{Claudio Lupi}
\examples{
  CADFpvalues(t0=-1.7, rho2=0.20, trend="trend")
}
