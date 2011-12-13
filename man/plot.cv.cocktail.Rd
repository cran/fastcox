\name{plot.cv.cocktail}
\alias{plot.cv.cocktail}
\title{plot the cross-validation curve produced by cv.cocktail}
\description{
Plots the cross-validation curve, and upper and lower standard deviation
curves, as a function of the \code{lambda} values used. This function is modified based on the \code{plot.cv} function from the \code{glmnet} package.}
\usage{
\method{plot}{cv.cocktail}(x, sign.lambda, ...)
}
\arguments{
	\item{x}{fitted \code{\link{cv.cocktail}} object}
	\item{sign.lambda}{either plot against \code{log(lambda)} (default) or
	its negative if \code{sign.lambda=-1}.}
	\item{\dots}{other graphical parameters to plot}
}
\details{A plot is produced.}

\author{Yi Yang and Hui Zou\cr
Maintainer: Yi Yang  <yiyang@umn.edu>}
\references{
Yang, Y. and Zou, H. (2012), 
"A Cocktail Algorithm for Solving The Elastic Net Penalized Cox's Regression in High Dimensions", \emph{Statistics and Its Interface.}\cr
\url{http://code.google.com/p/fastcox/}\cr


Friedman, J., Hastie, T., and Tibshirani, R. (2010), "Regularization paths for generalized
linear models via coordinate descent," \emph{Journal of Statistical Software, 33, 1.}\cr
\url{http://www.jstatsoft.org/v33/i01/}}

\seealso{\code{\link{cv.cocktail}}.}
\keyword{models}
\keyword{regression}
