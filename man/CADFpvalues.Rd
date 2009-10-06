\name{CADFpvalues}
\alias{CADFpvalues}
\title{p-values of the CADF test.}
\description{
  The asymptotic p-values of the Hansen's (1995) Covariate-Augmented Dickey Fuller (CADF) test 
  for a unit root are computed using the approach outlined in Costantini et al. (2007).
}
\usage{
CADFpvalues(t0, rho2 = 0.5, type=c("trend", "drift", "none")) 
}
\arguments{
  \item{t0}{the value of the test statistic;}
  \item{rho2}{the estimated value of the long-run correlation;}
  \item{type}{defines the deterministic kernel used in the test. It accepts the values used in 
               package \code{urca}. It specifies if the underlying model must be with 
               linear trend ("trend", the default), with constant ("drift") or without constant 
			   ("none").}
}
\value{
  \code{p.value}, a scalar containing the estimated asymptotic p-value of the test.
}
\references{
    Hansen BE (1995). Rethinking the Univariate Approach to Unit Root Testing: Using
	Covariates to Increase Power, \emph{Econometric Theory}, \bold{11}(5), 1148--1171.

    Costantini M, Lupi C, Popp S (2007), A Panel-CADF Test for Unit Roots, University of Molise,
    \emph{Economics \& Statistics Discussion Paper} 39/07.
    \url{http://econpapers.repec.org/paper/molecsdps/esdp07039.htm}
}
\author{Claudio Lupi}
\examples{
  CADFpvalues(t0=-1.7, rho2=0.20, type="trend")
}
