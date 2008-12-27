\name{summary.CADFtest}
\alias{summary.CADFtest}
\title{Function to print a summary of CADFtest objects}
\description{
  This function conveniently prints the detailed results of the Covariate-Augmented Dickey Fuller 
  test carried out in \code{CADFtest}.
}
\usage{
  \method{summary}{CADFtest}(object, ...)
}
\arguments{
  \item{object}{an object belonging to the class \code{CADFtest}.}
  \item{...}{other.}
}
\author{Claudio Lupi}
\examples{
  data(npext, package="urca")
  ADFt <- CADFtest(npext$realgnp, trend="ct")
  summary(ADFt)
}
