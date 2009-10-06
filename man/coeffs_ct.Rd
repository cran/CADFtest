\name{coeffs_ct}
\alias{coeffs_ct}
\alias{coeffs_nc}
\alias{coeffs_c}
\docType{data}
\title{Tables of coefficients to compute p-values as estimated by Costantini et al. (2007)}
\description{
  These tables contain the response surface coefficients needed to compute the p-value of Hansen's
  CADF test (see Costantini et al., 2007; Lupi, 2009). \code{coeffs_ct}, \code{coeffs_c}, 
  \code{coeffs_nc} are the relevant tables for the constant plus trend, constant, and 
  no constant case, respectively.  
}
\usage{
data("coeffs_ct")
}

\format{
  The tables are saved as binary data \code{.rda} objects. They are \eqn{(1005 \times 5)} matrices, 
  where the first column represents probabilities and the following four columns are 
  \eqn{\beta_0, \ldots, \beta_3} in eqn. (13) in Costantini et al. (2007) (see also Lupi, 2009).
}
\source{
    Costantini et al. (2007).
}
\references{
    Costantini M, Lupi C, Popp S (2007), A Panel-CADF Test for Unit Roots, University of Molise,
    \emph{Economics \& Statistics Discussion Paper} 39/07.
    \url{http://econpapers.repec.org/paper/molecsdps/esdp07039.htm}.

    Lupi C (2009). Unit Root CADF Testing with R, \emph{Journal of Statistical Software},
    \bold{32}(2), 1--19. \url{http://www.jstatsoft.org/v32/i02/}
}
\author{Claudio Lupi}
\note{If these tables are used, please cite Costantini et al. (2007).} 
\keyword{datasets}
