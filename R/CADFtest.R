CADFtest <- function(model, X=NULL, trend=c("c", "nc", "ct", "none", "drift", "trend"), 
                     data=list(), max.lag.y=1, min.lag.X=0, max.lag.X=0, dname="",
                     Auto=FALSE, criterion=c("BIC", "AIC"), prewhite=FALSE,
                     kernel = c("Parzen", "Quadratic Spectral", "Truncated", "Bartlett", "Tukey-Hanning"))
UseMethod("CADFtest")
  
