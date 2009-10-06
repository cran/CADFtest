\name{plot.CADFtest}
\alias{plot.CADFtest}
\title{Function to plot CADFtest objects}
\description{
  This function conveniently plots the residuals of the Covariate Augmented Dickey-Fuller
  or the standard Augmented Dickey-Fuller regression carried out in \code{CADFtest}.
}
\usage{
  \method{plot}{CADFtest}(x, plots=(1:4), ...)
}
\arguments{
  \item{x}{an object belonging to the class \code{CADFtest}.}
  \item{plots}{the plots to be produced (all the four plots by default): 
		1: standardized residuals plot;
		2: density of the residuals, with an indication of the p-value of the Jarque-Bera test
			for normality;
		3: ACF of the residuals;
		4: partial ACF of the residuals.}
  \item{...}{currently not used.}
}
\author{Claudio Lupi}
\examples{
  data(npext, package="urca")
  ADFt <- CADFtest(npext$realgnp, type="trend")
  plot(ADFt, plots=c(3,4)) 
}
