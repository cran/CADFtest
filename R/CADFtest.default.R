CADFtest.default <- function(model, X=NULL, trend=c("c", "nc", "ct", "none", "drift", "trend"), 
                             data=list(), max.lag.y=1, min.lag.X=0, max.lag.X=0, dname="",
                             Auto=FALSE, criterion=c("BIC", "AIC"), prewhite=FALSE,
                             kernel = c("Parzen", "Quadratic Spectral", "Truncated", "Bartlett", "Tukey-Hanning"))
{
# Author:       Claudio Lupi
# This version: December 12, 2008
# This function computes Hansen's (1995) Covariate-Augmented Dickey-Fuller (CADF) test.
# The only required argument is y, the Tx1 time series to be tested (y can be a vector).
# If no time series of stationary covariates X is passed to the procedure, then an ordinary ADF test is performed.
# The test.types are no-constant (nc), constant (c, the default value), constant plus trend (ct).
#
# max.lag.y >= 0
# min.lag.X <= 0
# max.lag.X >= 0

if (dname=="") dname <- deparse(substitute(model))

method <- "CADF test"; if (is.null(X)) method <- "ADF test"

y <- model
trend     <- match.arg(trend)
criterion <- match.arg(criterion)
kernel    <- match.arg(kernel)

switch(trend,
       "none" = trend <- "nc",
       "drift" = trend <- "c",
       "trend" = trend <- "ct")

rho2 <- NULL # default value for rho^2
nX   <- 0    # default number of covariates. The exact number is computed below

if (is.ts(y)==FALSE) y <- ts(y)
t <- ts(1:length(y), start=start(y), frequency=frequency(y))

#############################################################################################################
if (Auto==FALSE)  # no automatic model selection
{
  model <- "d(y) ~ "
  if (trend=="ct") model <- paste(model, "t +", sep="")
  model <- paste(model, " L(y, 1)", sep="")

  if (max.lag.y > 0)
  {
    for (i in 1:max.lag.y) model <- paste(model, " + L(d(y), ",i,")", sep="")
  }

  if (is.null(X)==FALSE)
  {
    if (is.ts(X)==FALSE) X <- ts(X, start=start(y), frequency=frequency(y))
    nX <- 1; if (is.null(dim(X))==FALSE) nX <- dim(X)[2]  # number of covariates
    nX <- (max.lag.X - min.lag.X + 1)*nX              # number of X's (including the lags)
    if ((min.lag.X==0) & (max.lag.X==0)) model <- paste(model, " + L(X, 0)", sep="")
    if ((min.lag.X!=0) | (max.lag.X!=0))
    {
      for (i in min.lag.X:max.lag.X) model <- paste(model, " + L(X, ",i,")", sep="")
    }
  }

  if (trend=="nc") model <- paste(model, " -1", sep="")

  est.model      <- dynlm(formula=formula(model))
  summ.est.model <- summary(est.model)
  q              <- summ.est.model$df[1]
  TT             <- q + summ.est.model$df[2]

  model.AIC <- AIC(est.model)
  model.BIC <- model.AIC + (log(TT)-2)*q

  t.value <- summ.est.model$coefficients[(2 - as.numeric(trend=="nc") + as.numeric(trend=="ct")),3]

  if (is.null(X)) p.value <- punitroot(t.value, trend = trend, statistic = "t") # MacKinnon p-values

  if (is.null(X)==FALSE)
  {
    # Compute Hansen's p-value
    k <- length(est.model$coefficients)
    series  <- as.matrix(est.model$model)          # explanatory variables considered in the model (excluding constant)
    nseries <- dim(series)[2]                      # number of variables
    Xseries <- series[,(nseries-nX+1):nseries]     # the X's are the last nX columns of series
    if (nX==1) Xseries <- Xseries - mean(Xseries)  # demean the X's
    if (nX>1)  Xseries <- Xseries - apply(Xseries,2,mean)
    e <- as.matrix(est.model$residuals)
    if (nX==1) v <- Xseries * est.model$coefficients[k] + e
    if (nX>1)  v <- Xseries%*%est.model$coefficients[(k-nX+1):k] + e
    V <- cbind(e,v)
    mod <- lm(V~1)
    LRCM <- (kernHAC(mod, kernel=kernel, prewhite=prewhite))*nrow(V)
    rho2 <- LRCM[1,2]^2/(LRCM[1,1]*LRCM[2,2])

    p.value <- CADFpvalues(t.value, rho2, trend)
  }

  test.results <- list(statistic=t.value,
                       parameter=c("rho2" = rho2),
                       method=method,
                       p.value=as.vector(p.value),
                       data.name=dname,
                       max.lag.y=max.lag.y,
                       min.lag.X=min.lag.X,
                       max.lag.X=max.lag.X,
                       AIC=model.AIC,
                       BIC=model.BIC,
                       est.model=est.model) 
}
#############################################################################################################

#############################################################################################################
if (Auto==TRUE)  # automatic model selection
{
  all.models <- expand.grid(0:max.lag.y, min.lag.X:0, 0:max.lag.X)  # all possible models
  models.num <- dim(all.models)[1]                                  # number of models to be estimated 
  AICBICmatrix <- matrix(NA, models.num, 5)                         # matrix to store lag orders, AIC and BIC

  for (modeln in 1:models.num)
  {
    max.lag.y <- all.models[modeln, 1]
    min.lag.X <- all.models[modeln, 2]
    max.lag.X <- all.models[modeln, 3]

    model <- "d(y) ~ "
    if (trend=="ct") model <- paste(model, "t +", sep="")
    model <- paste(model, " L(y, 1)", sep="")

    if (max.lag.y > 0)
    {
      for (i in 1:max.lag.y) model <- paste(model, " + L(d(y), ",i,")", sep="")
    }

    if (is.null(X)==FALSE)
    {
      if (is.ts(X)==FALSE) X <- ts(X, start=start(y), frequency=frequency(y))
      nX <- 1; if (is.null(dim(X))==FALSE) nX <- dim(X)[2]  # number of covariates
      nX <- (max.lag.X - min.lag.X + 1)*nX              # number of X's (including the lags)
      if ((min.lag.X==0) & (max.lag.X==0)) model <- paste(model, " + L(X, 0)", sep="")
      if ((min.lag.X!=0) | (max.lag.X!=0))
      {
        for (i in min.lag.X:max.lag.X) model <- paste(model, " + L(X, ",i,")", sep="")
      }
    }

    if (trend=="nc") model <- paste(model, " -1", sep="")

    est.model      <- dynlm(formula=formula(model))

    summ.est.model <- summary(est.model)
    q              <- summ.est.model$df[1]
    TT             <- q + summ.est.model$df[2]

    model.AIC <- AIC(est.model)
    model.BIC <- model.AIC + (log(TT)-2)*q

    AICBICmatrix[modeln, ] <- c(max.lag.y, min.lag.X, max.lag.X, model.AIC, model.BIC)
  }
  
  if (criterion=="AIC") selected.model <- which(AICBICmatrix[,4]==min(AICBICmatrix[,4]))
  if (criterion=="BIC") selected.model <- which(AICBICmatrix[,5]==min(AICBICmatrix[,5]))
  
  if (length(selected.model) > 1) selected.model <- selected.model[1]

  max.lag.y <- AICBICmatrix[selected.model, 1]
  min.lag.X <- AICBICmatrix[selected.model, 2]
  max.lag.X <- AICBICmatrix[selected.model, 3]

  ################################## ESTIMATION & TEST WITH THE SELECTED MODEL ##############################
  model <- "d(y) ~ "
  if (trend=="ct") model <- paste(model, "t", sep="")
  model <- paste(model, " + L(y, 1)", sep="")

  if (max.lag.y > 0)
  {
    for (i in 1:max.lag.y) model <- paste(model, " + L(d(y), ",i,")", sep="")
  }

  if (is.null(X)==FALSE)
  {
    if (is.ts(X)==FALSE) X <- ts(X, start=start(y), frequency=frequency(y))
    nX <- 1; if (is.null(dim(X))==FALSE) nX <- dim(X)[2]  # number of covariates
    nX <- (max.lag.X - min.lag.X + 1)*nX              # number of X's (including the lags)
    if ((min.lag.X==0) & (max.lag.X==0)) model <- paste(model, " + L(X, 0)", sep="")
    if ((min.lag.X!=0) | (max.lag.X!=0))
    {
      for (i in min.lag.X:max.lag.X) model <- paste(model, " + L(X, ",i,")", sep="")
    }
  }

  if (trend=="nc") model <- paste(model, " -1", sep="")

  est.model      <- dynlm(formula=formula(model))
  summ.est.model <- summary(est.model)
  q              <- summ.est.model$df[1]
  TT             <- q + summ.est.model$df[2]

  model.AIC <- AIC(est.model)
  model.BIC <- model.AIC + (log(TT)-2)*q

  t.value <- summ.est.model$coefficients[(2 - as.numeric(trend=="nc") + as.numeric(trend=="ct")),3]

  if (is.null(X)) p.value <- punitroot(t.value, trend = trend, statistic = "t") # MacKinnon p-values

  if (is.null(X)==FALSE)
  {
    # Compute Hansen's p-value
    k <- length(est.model$coefficients)
    series  <- as.matrix(est.model$model)          # explanatory variables considered in the model (excluding constant)
    nseries <- dim(series)[2]                      # number of variables
    Xseries <- series[,(nseries-nX+1):nseries]     # the X's are the last nX columns of series
    if (nX==1) Xseries <- Xseries - mean(Xseries)  # demean the X's
    if (nX>1)  Xseries <- Xseries - apply(Xseries,2,mean)
    e <- as.matrix(est.model$residuals)
    if (nX==1) v <- Xseries * est.model$coefficients[k] + e
    if (nX>1)  v <- Xseries%*%est.model$coefficients[(k-nX+1):k] + e
    V <- cbind(e,v)
    mod <- lm(V~1)
    LRCM <- (kernHAC(mod, kernel=kernel, prewhite=prewhite))*nrow(V)
    rho2 <- LRCM[1,2]^2/(LRCM[1,1]*LRCM[2,2])

    p.value <- CADFpvalues(t.value, rho2, trend)
  }

  test.results <- list(statistic=t.value,
                       parameter=c("rho2" = rho2),
                       method=method,
                       p.value=as.vector(p.value),
                       data.name=dname,
                       max.lag.y=max.lag.y,
                       min.lag.X=min.lag.X,
                       max.lag.X=max.lag.X,
                       AIC=model.AIC,
                       BIC=model.BIC,
                       est.model=est.model) 
}

class(test.results) <- c("CADFtest", "htest")
if (is.null(X)){
  names(test.results$statistic) <- paste("ADF(",max.lag.y,")",sep="")}
else{
  names(test.results$statistic) <- paste("CADF(",max.lag.y,",",max.lag.X,",",min.lag.X,")",sep="")}
test.results$estimate <- c("delta" = as.vector(test.results$est.model$coefficients[(2 - as.numeric(trend=="nc") + as.numeric(trend=="ct"))])) 
test.results$null.value <- c("delta" = 0) 
test.results$alternative <- "less" 

return(test.results)
}
