\name{update.CADFtest}
\alias{update.CADFtest}
\title{Function to update the formula of CADFtest objects}
\description{
  This function updates the formula and/or the other arguments of CADFtest object and re-run the test 
  using the updated arguments.
  It can be useful if one wants to see the effect of adding/removing stationary covariates or the effect
  of changing lags, kernel, etc. If covariates have to be added/removed, \code{update()} works only if
  \code{model} is passed as a formula.
}
\usage{
  \method{update}{CADFtest}(object, change, ...)
}
\arguments{
  \item{object}{an object belonging to the class \code{CADFtest}.}
  \item{change}{list of charater describing the changes to be applied to the existing model.}
  \item{...}{currently not used.}
}
\value{
	The function re-run the test and returns an object of class \code{CADFtest}. See \code{CADFtest()}.
}
\author{Claudio Lupi}
\examples{
  data(npext, package="urca")
  npext$unemrate <- exp(npext$unemploy)      # compute unemployment rate
  L <- ts(npext, start=1860)                 # time series of levels
  D <- diff(L)                               # time series of diffs
  S <- window(ts.intersect(L,D), start=1909) # select same sample as Hansen's
  CADFt <- CADFtest(L.gnpperca~D.unemrate, data=S, max.lag.y=3,
    kernel="Parzen", prewhite=FALSE)
  CADFt.2 <- update(CADFt, change=list("+ D.indprod", "max.lag.X=3", 
	"criterion='BIC'"))
}
