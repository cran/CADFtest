summary.CADFtest <- function(object, ...)
{
  # object is an object of class CADFtest
  cat("Covariate-Augmented Dickey Fuller (CADF) test. \n")
  if (is.null(object$parameter)) 
  {
    cat("No covariate used in the analysis: standard ADF is performed. \n")
  }
  rval <- sprintf("Test statistic t:                       %s", round(object$statistic,4))
  cat(rval, "\n")
  if (!is.null(object$parameter))
  {
    rval <- sprintf("Estimated rho^2:                         %s", round(object$parameter,4))
    cat(rval, "\n")
  }
  rval <- sprintf("Test p-value:                            %s", round(object$p.value,4))
  cat(rval, "\n")
  rval <- sprintf("Max lag of the diff. dependent variable: %s", object$max.lag.y)
  cat(rval, "\n")
  if (!is.null(object$parameter))
  {
    rval <- sprintf("Max lag of the stationary covariate(s):  %s", object$max.lag.X)
    cat(rval, "\n")
    rval <- sprintf("Max lead of the stationary covariate(s): %s", object$min.lag.X)
    cat(rval, "\n")
  }
  
#   if (object$p.value > 0.10) cat("The I(1) null is not rejected. \n")
#   if ((object$p.value < 0.10)&(object$p.value > 0.05)) cat("The I(1) null is rejected at the 10% significance level. \n")
#   if ((object$p.value < 0.05)&(object$p.value > 0.01)) cat("The I(1) null is rejected at the 5% significance level. \n")
#   if (object$p.value < 0.01) cat("The I(1) null is rejected at the 1% significance level. \n")
    
  {
    s <- summary.lm(object$est.model)
    print(s, signif.stars=FALSE)
  }
}
