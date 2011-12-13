\name{cv.cocktail}
\alias{cv.cocktail}
\alias{cv.survpath}
\title{Cross-validation for cocktail}
\description{Does k-fold cross-validation for cocktail, produces a plot,
and returns a value for \code{lambda}. This function is modified based on the \code{cv} function from the \code{glmnet} package.}
\usage{
cv.cocktail(x,y,d,lambda=NULL,...,nfolds=10,foldid)
}
\arguments{
\item{x}{matrix of predictors, of dimension \eqn{N \times p}{N*p}; each row is an observation vector.}

\item{y}{a survival time for Cox models. Currently tied failure times are not supported.}

\item{d}{censor status with 1 if died and 0 if right censored.}

		\item{lambda}{optional user-supplied lambda sequence; default is
		\code{NULL}, and \code{\link{cocktail}} chooses its own sequence.}
		\item{nfolds}{number of folds - default is 5. Although \code{nfolds}
		can be as large as the sample size (leave-one-out CV), it is not
		recommended for large datasets. Smallest value allowable is \code{nfolds=3}.}
		\item{foldid}{an optional vector of values between 1 and \code{nfold}
		identifying what fold each observation is in. If supplied,
		\code{nfold} can be missing.}
		\item{\dots}{other arguments that can be passed to cocktail.}

}
\details{The function runs \code{\link{cocktail}} \code{nfolds}+1 times; the
first to get the \code{lambda} sequence, and then the remainder to
compute the fit with each of the folds omitted. The average error and standard deviation over the
folds are computed. 
}
\value{an object of class \code{\link{cv.cocktail}} is returned, which is a
list with the ingredients of the cross-validation fit.
		\item{lambda}{the values of \code{lambda} used in the fits.}
		\item{cvm}{the mean cross-validated error - a vector of length
		\code{length(lambda)}.}
		\item{cvsd}{estimate of standard error of \code{cvm}.}
		\item{cvup}{upper curve = \code{cvm+cvsd}.}
		\item{cvlo}{lower curve = \code{cvm-cvsd}.}
		\item{nzero}{number of non-zero coefficients at each \code{lambda}.}
		\item{name}{a text string indicating partial likelihood (for plotting purposes).}
		\item{cocktail.fit}{a fitted \code{\link{cocktail}} object for the full data.}
		\item{lambda.min}{The optimal value of \code{lambda} that gives minimum
		cross validation error \code{cvm}.}
		\item{lambda.1se}{The largest value of \code{lambda} such that error is
		within 1 standard error of the minimum.}
}


\author{Yi Yang and Hui Zou\cr
Maintainer: Yi Yang  <yiyang@umn.edu>}
\references{
Yang, Y. and Zou, H. (2012), "A Cocktail Algorithm for Solving The Elastic Net Penalized Cox's Regression in High Dimensions", \emph{Statistics and Its Interface.}\cr
\url{http://code.google.com/p/fastcox/}\cr

Friedman, J., Hastie, T., and Tibshirani, R. (2010), "Regularization paths for generalized
linear models via coordinate descent," \emph{Journal of Statistical Software, 33, 1.}\cr
\url{http://www.jstatsoft.org/v33/i01/}
}


\seealso{\code{\link{cocktail}}, \code{\link{plot.cv.cocktail}}.}
\examples{
data(FHT)
cv1<-cv.cocktail(x=FHT$x[,1:10],y=FHT$y,d=FHT$status,alpha=0.5,nfolds=3)
cv1
plot(cv1)
}
\keyword{models}
\keyword{regression}

