\name{CADFtest}
\alias{CADFtest}
\alias{CADFtest.formula}
\alias{CADFtest.default}
\title{Hansen's Covariate-Augmented Dickey Fuller (CADF) test}
\description{
  This function is an interface to \code{CADFtest.default} that computes the CADF unit root test
  proposed in Hansen (1995). The asymptotic p-values of the test are also computed. Automatic 
  model selection is allowed.
}
\usage{
CADFtest(model, X=NULL, trend=c("c", "nc", "ct", "none", "drift", "trend"), 
                             data=list(), max.lag.y=1, min.lag.X=0, max.lag.X=0, dname="",
                             Auto=FALSE, criterion=c("BIC", "AIC"), prewhite=FALSE,
                             kernel = c("Parzen", "Quadratic Spectral", "Truncated", "Bartlett", "Tukey-Hanning"))
}

\arguments{
  \item{model}{a formula of the kind \code{y ~ x1 + x2} containing the variable \code{y} to be tested and 
               the stationary covariate(s) to be used in the test. If \code{model=y} is specified where \code{y}
               is a vector or a time series, then an ordinary ADF test is performed. It is not the actual model, 
               it is used to simplify variable specification. The covariates are assumed to be stationary.}
  \item{X}{if \code{model=y}, a matrix or a vector time series of stationary covariates \code{X} can be passed 
           directly, instead of using the \code{formula} expression.} 
  \item{trend}{defines the deterministic kernel used in the test. It accepts the values used either in 
               \code{fUnitRoots} or \code{urca} packages. It specifies if the underlying model must be with 
               constant ("c" or "drift", the default), without constant ("nc" or "none"), or with constant and 
               trend ("ct" or "trend").}
  \item{data}{data to be used (optional).}
  \item{max.lag.y}{maximum number of lags allowed for the lagged differences of the variable to be tested.}
  \item{min.lag.X}{if negative it is maximum lead allowed for the covariates. If zero, it is the minimum lag 
                   allowed for the covariates.}
  \item{max.lag.X}{maximum lag allowed for the covariates.}
  \item{dname}{character. Data name. In general there is no need to modify the default value. The correct data 
               name is computed on the basis of the \code{model} passed to the function.} 
  \item{Auto}{logical. If \code{Auto==FALSE} then the test is performed using the given orders of lags and leads.
              If \code{Auto==TRUE} then the test is performed using the model that minimizes the selection 
              citerion defined in \code{criterion}. In this case, the max e min orders serve as upper and 
              lower bounds in the model selection.}
  \item{criterion}{it can be either \code{"BIC"} or \code{"AIC"}. It is effective only when \code{Auto==T}.}
  \item{prewhite}{logical or integer. Should the estimating functions be prewhitened? If TRUE or greater than 0 
             a VAR model of order \code{as.integer(prewhite)} is fitted via \code{ar} with method \code{"ols"} 
             and \code{demean = FALSE}. The default is to use no prewhitening (\code{prewhite = FALSE}).}
  \item{kernel}{a character specifying the kernel used. All kernels used are described in Andrews (1991). }
}
%\details{
%   ~~ If necessary, more details than the description above ~~
%}
\value{
  The function returns an object of class \code{c("CADFtest", "htest")} containing:
  \item{statistic}{the t test statistic.}
  \item{parameter}{the nuisance rho2 parameter.}
  \item{method}{the test performed: it can be either \code{ADF} or \code{CADF}.}
  \item{p.value}{the p-value of the test.}
  \item{data.name}{the data name.}
  \item{max.lag.y}{the maximum lag of the differences of the dependent variable.}
  \item{min.lag.X}{the maximum lead of the stationary covariate(s).}
  \item{max.lag.X}{the maximum lag of the stationary covariate(s).}
  \item{AIC}{the value of the AIC for the selected model.}
  \item{BIC}{the value of the BIC for the selected model.}
  \item{est.model}{the estimated model.}
  \item{estimate}{the estimated value of the parameter of the lagged dependent variable.}
  \item{null.value}{the value of the parameter of the lagged dependent variable under the null.}
  \item{alternative}{the alternative hypothesis.}
}
\references{
    Hansen, BE (1995): "Rethinking the Univariate Approach to Unit Root Testing: Using
	Covariates to Increase Power", \emph{Econometric Theory}, \bold{11} (5), 1148--1171.

    Costantini M, Lupi C, Popp S (2007), "A Panel-CADF Test for Unit Roots", University of Molise,
    \emph{Economics \& Statistics Discussion Paper} 39/07, 
    URL \url{http://econpapers.repec.org/paper/molecsdps/esdp07039.htm}.
}
\author{Claudio Lupi}
\note{If p-values are used, please cite Costantini, Lupi and Popp (2007).} 

\seealso{\code{fUnitRoots}, \code{urca}}
\examples{
##---- ADF test on extended Nelson-Plosser data ----
##--   Data taken from package urca
  data(npext, package="urca")
  ADFt   <- CADFtest(npext$gnpperca, max.lag.y=3, trend="trend")

##---- CADF test on extended Nelson-Plosser data ----
  data(npext, package="urca")
  npext$unemrate <- exp(npext$unemploy)      # compute unemployment rate
  L <- ts(npext, start=1860)                 # time series of levels
  D <- diff(L)                               # time series of diffs
  S <- window(ts.intersect(L,D), start=1909) # select same sample as Hansen's
  CADFt <- CADFtest(L.gnpperca~D.unemrate, data=S, trend="ct", max.lag.y=3)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ts}
\keyword{htest}
\keyword{univar}
