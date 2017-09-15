\name{print.cocktail}
\alias{print.cocktail}
\title{print a cocktail object}
\description{
Print a summary of the cocktail path at each step along the path. This function is modified based on the \code{print} function from the \code{glmnet} package.
}
\usage{
\method{print}{cocktail}(x, digits = max(3, getOption("digits") - 3), ...)
}
\arguments{
  \item{x}{fitted \code{\link{cocktail}} object}
  \item{digits}{significant digits in printout}
  \item{\dots}{additional print arguments}
}
\details{
The call that produced the \code{\link{cocktail}} object is printed, followed by a two-column matrix with columns \code{Df} and \code{Lambda}. The \code{Df} column is the number of nonzero coefficients.
}
\value{
a two-column matrix, the first columns is the number of nonzero coefficients and the second column is \code{Lambda}.}

\author{Yi Yang and Hui Zou\cr
Maintainer: Yi Yang  <yi.yang6@mcgill.ca>}

\references{
Yang, Y. and Zou, H. (2013), 
"A Cocktail Algorithm for Solving The Elastic Net Penalized Cox's Regression in High Dimensions", \emph{Statistics and Its Interface,} 6:2, 167-173.\cr
\url{https://github.com/emeryyi/fastcox}\cr

Friedman, J., Hastie, T. and Tibshirani, R. (2008)
  "Regularization Paths for Generalized Linear Models via Coordinate
    Descent",   \url{http://www.stanford.edu/~hastie/Papers/glmnet.pdf}\cr
  \emph{Journal of Statistical Software}, Vol. 33(1), 1-22 Feb 2010\cr
  \url{http://www.jstatsoft.org/v33/i01/}\cr\cr
  Simon, N., Friedman, J., Hastie, T., Tibshirani, R. (2011)
  "Regularization Paths for Cox's Proportional Hazards Model via
    Coordinate Descent", \emph{Journal of Statistical Software}, Vol. 39(5)
    1-13 \cr
  \url{http://www.jstatsoft.org/v39/i05/}
}

\examples{
data(FHT)
m1<-cocktail(x=FHT$x,y=FHT$y,d=FHT$status,alpha=0.5)
print(m1)
}
\keyword{models}
\keyword{regression}

