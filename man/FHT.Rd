\name{FHT}
\alias{FHT}
\docType{data}
\title{FHT data introduced in Simon et al. (2011).}
\description{
The \code{FHT} data set has n = 50 observations and p = 100 predictors. The covariance between predictors Xj and Xj' has the same correlation 0.5. See details in Simon et al. (2011).
}
\usage{data(FHT)}
\format{
  This \code{list} object named "FHT" contains the following data:
  \describe{
    \item{x}{a covariate matrix with 50 rows and 100 columns}
    \item{y}{the distinct failure times}
    \item{status}{the censoring indicator (status = 1 indicates no censoring and status = 0 indicates right censoring)}
  }
}
\references{Friedman, J., Hastie, T. and Tibshirani, R. (2008)
  \emph{Regularization Paths for Generalized Linear Models via Coordinate
    Descent},   \url{http://www.stanford.edu/~hastie/Papers/glmnet.pdf}\cr
  \emph{Journal of Statistical Software, Vol. 33(1), 1-22 Feb 2010}\cr
  \url{http://www.jstatsoft.org/v33/i01/}\cr\cr
  Simon, N., Friedman, J., Hastie, T., Tibshirani, R. (2011)
  \emph{Regularization Paths for Cox's Proportional Hazards Model via
    Coordinate Descent, Journal of Statistical Software, Vol. 39(5)
    1-13}\cr
  \url{http://www.jstatsoft.org/v39/i05/}
  }
\examples{
data(FHT)
}
\keyword{datasets}
