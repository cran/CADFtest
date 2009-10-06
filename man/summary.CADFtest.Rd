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
  \item{...}{currently not used.}
}
\value{
	The function returns an object of class \code{CADFtestsummary} containing the main
	results of the test.
	\item{test.summary}{a matrix, containing the t-test statistic, the estimated value of \eqn{\rho^2},
		the p-value of the test, the max lag used for the differenced dependent variable, 
		the max lag of the stationary covariate, the max lead of the stationary covariate.
		When a standard ADF test is performed, only the t-test statistic, the p-value  and
		the max lag of the differenced dependent variable are reported.}
	\item{model.summary}{the summary of the test model, in the usual form. However, note that the 
		p-value of the lagged dependent is computed under the null of a unit root. Furthermore,
		differently from the common practice, the F-statistic refers to the joint significance of
		the stationary regressors. If no stationary regressors are used (no lagged differences 
		of the dependent, no stationary covariates) then the F-statistic is not computed and a
		\code{NA} value is returned.}
}
\author{Claudio Lupi}
\examples{
  data(npext, package="urca")
  ADFt <- CADFtest(npext$realgnp, type="trend")
  summary(ADFt)
}
