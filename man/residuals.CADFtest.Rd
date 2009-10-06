\name{residuals.CADFtest}
\alias{residuals.CADFtest}
\title{Function to extract the residuals from CADFtest objects}
\description{
  This function applies the \code{residuals()} method to an object of class
  \code{CADFtest}.
}
\usage{
  \method{residuals}{CADFtest}(object, ...)
}
\arguments{
  \item{object}{an object belonging to the class \code{CADFtest}.}
  \item{...}{currently not used.}
}
\author{Claudio Lupi}
\examples{
  data(npext, package="urca")
  ADFt <- CADFtest(npext$realgnp, type="trend")
  residuals(ADFt)
}
