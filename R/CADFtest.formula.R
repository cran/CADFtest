CADFtest.formula <- function(model, X=NULL, trend=c("c", "nc", "ct", "none", "drift", "trend"), 
                             data=list(), max.lag.y=1, min.lag.X=0, max.lag.X=0, dname="",
                             Auto=FALSE, criterion=c("BIC", "AIC"), prewhite=FALSE,
                             kernel = c("Parzen", "Quadratic Spectral", "Truncated", "Bartlett", "Tukey-Hanning"))
{
# Author:       Claudio Lupi
# This version: December 12, 2008
# This function is an interface to function CADFtest.default that computes Hansen's (1995) CADF test.
# Reference:
# @ARTICLE{,
#   author = {Hansen, Bruce E.},
#   title = {Rethinking the Univariate Approach to Unit Root Testing: {U}sing Covariates to Increase Power},
#   journal = {Econometric Theory},
#   year = {1995},
#   volume = {11},
#   pages = {1148--1171},
#   number = {5},
# }
#
# Arguments:
# model:     a formula. If model is a vector or a time series then the 
#            standard ADF test is performed on the series described by model. If a CADF test is desired, then
#            model should specified in the form z0 ~ z1 + z2 + ... where z0 is the variable to be tested, 
#            while z1 and z2 are the stationary covariates to be used in the test. Note that the model is stylized
#            and all the variables are in levels. It is NOT the model equation on which the test is based.
#            However, the covariates must be STATIONARY. 
# trend:     it specifies if the underlying model must be with constant ("c", the default), without constant ("nc"),
#            or with constant and trend ("ct").
# max.lag.y: it specifies the number of lags of the dependent (\Delta y_t).
# min.lag.X: it specifies the maximum lead of the covariates (it must be negative or zero).
# max.lag.X: it specifies the maximum lag of the covariates (it must be positive or zero).
# Auto:      logical. If Auto==FALSE then the test is performed using the given orders of lags and leads. If Auto==TRUE
#            then the test is performed using the model that minimizes the selection citerion defined in
#            'criterion'. In this case, the max e min orders serve as upper and lower bounds in the model
#            selection.
# criterion: it can be either "BIC" or "AIC". It is effective only when Auto==TRUE.
# prewhite:  logical or integer. Should the estimating functions be prewhitened? If TRUE or greater than 0 
#            a VAR model of order as.integer(prewhite) is fitted via ar with method "ols" and demean = FALSE. 
#            The default is to use no prewhitening.
# kernel:    a character specifying the kernel used. All kernels used are described in Andrews (1991).
#
# The procedure to compute the CADF test p-value is proposed in Costantini et al. (2007). Please cite the paper
# when you use the present function.
# Reference:
# @TECHREPORT{,
#   author = {Costantini, Mauro and Lupi, Claudio and Popp, Stephan},
#   title = {A Panel-{CADF} Test for Unit Roots},
#   institution = {University of Molise},
#   year = {2007},
#   type = {Economics \& Statistics Discussion Paper},
#   number = {39/07},
#   url = {http://econpapers.repec.org/paper/molecsdps/esdp07039.htm}
# }

dname <- deparse(substitute(model))
mf    <- model.frame(model, data=data)
y     <- model.response(mf)
X     <- model.matrix(model, data=data)
X     <- X[,2:dim(X)[2]]

test.results <- CADFtest.default(model=y, X=X, trend=trend, max.lag.y=max.lag.y, 
                                 min.lag.X=min.lag.X, max.lag.X=max.lag.X, dname=dname,
                                 Auto=Auto, criterion=criterion, kernel=kernel, prewhite=prewhite)

return(test.results)
}

