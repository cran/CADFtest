\name{CADFtest}
\alias{CADFtest}
\alias{CADFtest.formula}
\alias{CADFtest.default}
\title{Hansen's Covariate-Augmented Dickey Fuller (CADF) test for unit roots}
\description{
  This function is an interface to \code{CADFtest.default} that computes the CADF unit root test
  proposed in Hansen (1995). The asymptotic p-values of the test are also computed along the lines
  proposed in Costantini et al. (2007). Automatic model selection is allowed. A full description 
  and some applications can be found in Lupi (2009).
}
\usage{
CADFtest(model, X=NULL, type=c("trend", "drift", "none"), 
     data=list(), max.lag.y=1, min.lag.X=0, max.lag.X=0, 
     dname=NULL, criterion=c("none", "BIC", "AIC", "HQC", 
     "MAIC"), ...)
}

\arguments{
  \item{model}{a formula of the kind \code{y ~ x1 + x2} containing the variable \code{y} to be tested and 
               the stationary covariate(s) to be used in the test. If the model is specified as
			   \code{y ~ 1}, then an ordinary ADF is carried out. Note that the specification
			   \code{y ~ .} here does not imply a model with all the disposable regressors, 
			   but rather a model with no stationary covariate (which correspons to an ADF test). 
			   This is because the stationary covariates have to be explicitly indicated (they 
			   are usually one or two). An ordinary ADF is performed also if \code{model=y} 
			   is specified, where \code{y} is a vector or a time series. It should 
			   be noted that \code{model} is not the actual model, but rather a representation that is 
			   used to simplify variable specification. The covariates are assumed to be stationary.}
  \item{X}{if \code{model=y}, a matrix or a vector time series of stationary covariates \code{X} can be passed 
			directly, instead of using the \code{formula} expression. However, the \code{formula}
			expression should in general be preferred.} 
  \item{type}{defines the deterministic kernel used in the test. It accepts the values used in package 
               \code{urca}. It specifies if the underlying model must be with linear trend ("trend", the
			   default), with constant ("drift") or without constant ("none").}
  \item{data}{data to be used (optional). This argument is effective only when \code{model} is passed as a
              formula.}
  \item{max.lag.y}{maximum number of lags allowed for the lagged differences of the variable to be tested.}
  \item{min.lag.X}{if negative it is maximum lead allowed for the covariates. If zero, it is the minimum lag 
                   allowed for the covariates.}
  \item{max.lag.X}{maximum lag allowed for the covariates.}
  \item{dname}{NULL or character. It can be used to give a special name to the model. 
  If the NULL default is accepted and the model is specified using a formula notation, then
  dname is computed according to the used formula.}
  \item{criterion}{it can be either \code{"none"} (the default), \code{"BIC"}, \code{"AIC"},
			\code{"HQC"} or \code{"MAIC"}. If \code{criterion="none"}, no automatic model selection 
			is performed. Otherwise, automatic model selection is performed using the specified 
			criterion. In this case, the max and min orders serve as upper and lower bounds in the 
			model selection.}          
  \item{...}{Extra arguments that can be set to use special kernels, prewhitening, etc. in the estimation of 
             \eqn{\rho^2}. A Quadratic kernel with a VAR(1) prewhitening is the default choice. To set 
             these extra arguments to different values, see \code{kernHAC} in package \code{sandwich} 
             (Zeileis, 2004, 2006). If Hansen's results have to be duplicated, then 
             \code{kernel="Parzen"} and \code{prewhite=FALSE} must be specified.}
}
%\details{
%   ~~ If necessary, more details than the description above ~~
%}
\value{
  The function returns an object of class \code{c("CADFtest", "htest")} containing:
  \item{statistic}{the t test statistic.}
  \item{parameter}{the estimated nuisance parameter \eqn{\rho^2} (see Hansen, 1995, p. 1150).}
  \item{method}{the test performed: it can be either \code{ADF} or \code{CADF}.}
  \item{p.value}{the p-value of the test.}
  \item{data.name}{the data name.}
  \item{max.lag.y}{the maximum lag of the differences of the dependent variable.}
  \item{min.lag.X}{the maximum lead of the stationary covariate(s).}
  \item{max.lag.X}{the maximum lag of the stationary covariate(s).}
  \item{AIC}{the value of the AIC for the selected model.}
  \item{BIC}{the value of the BIC for the selected model.}
  \item{HQC}{the value of the HQC for the selected model.}
  \item{MAIC}{the value of the MAIC for the selected model.}
  \item{est.model}{the estimated model.}
  \item{estimate}{the estimated value of the parameter of the lagged dependent variable.}
  \item{null.value}{the value of the parameter of the lagged dependent variable under the null.}
  \item{alternative}{the alternative hypothesis.}
  \item{call}{the call to the function.}
  \item{type}{the deterministic kernel used.}
}
\references{
    Costantini M, Lupi C, Popp S (2007). A Panel-CADF Test for Unit Roots, University of Molise,
    \emph{Economics & Statistics Discussion Paper} 39/07.
    \url{http://econpapers.repec.org/paper/molecsdps/esdp07039.htm}

    Hansen BE (1995). Rethinking the Univariate Approach to Unit Root Testing: Using
	Covariates to Increase Power, \emph{Econometric Theory}, \bold{11}(5), 1148--1171.

    Lupi C (2009). Unit Root CADF Testing with R, \emph{Journal of Statistical Software},
    \bold{32}(2), 1--19. \url{http://www.jstatsoft.org/v32/i02/}

    Zeileis A (2004). Econometric Computing with HC and HAC
    Covariance Matrix Estimators, \emph{Journal of Statistical Software}, \bold{11}(10),
    1--17. \url{http://www.jstatsoft.org/v11/i10/}

    Zeileis A (2006). Object-Oriented Computation of Sandwich
    Estimators, \emph{Journal of Statistical Software}, \bold{16}(9), 1--16.
    \url{http://www.jstatsoft.org/v16/i09/}.
}
\author{Claudio Lupi}

\seealso{\code{fUnitRoots}, \code{urca}}
\examples{
##---- ADF test on extended Nelson-Plosser data ----
##--   Data taken from package urca
  data(npext, package="urca")
  ADFt <- CADFtest(npext$gnpperca, max.lag.y=3, type="trend")

##---- CADF test on extended Nelson-Plosser data ----
  data(npext, package="urca")
  npext$unemrate <- exp(npext$unemploy)      # compute unemployment rate
  L <- ts(npext, start=1860)                 # time series of levels
  D <- diff(L)                               # time series of diffs
  S <- window(ts.intersect(L,D), start=1909) # select same sample as Hansen's
  CADFt <- CADFtest(L.gnpperca~D.unemrate, data=S, max.lag.y=3,
    kernel="Parzen", prewhite=FALSE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ts}
\keyword{htest}
\keyword{univar}
