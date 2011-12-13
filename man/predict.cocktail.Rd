\name{predict.cocktail}
\alias{predict.cocktail}
\alias{predict.survpath}
\title{make predictions from a "cocktail" object.}

\description{
Similar to other predict methods, this functions predicts fitted values, link function and more from a fitted \code{\link{cocktail}} object. This function is modified based on the \code{predict} function from the \code{glmnet} package.
}

\usage{
\method{predict}{cocktail}(object,newx,s=NULL,type=c("link","response","coefficients","nonzero"),...)
}

\arguments{
		\item{object}{fitted \code{\link{cocktail}} model object.}
		\item{newx}{matrix of new values for \code{x} at which predictions are
		to be made. Must be a matrix. This argument is not used for \code{type=c("coefficients","nonzero")}}
		\item{s}{value(s) of the penalty parameter \code{lambda} at which
		predictions are required. Default is the entire sequence used to
		create the model. }
		\item{type}{
		type of prediction required. 
		\itemize{
				\item Type \code{"link"} gives the linear predictors for Cox's model. 
				\item Type \code{"response"} gives the fitted relative-risk for Cox's model. 
                \item Type \code{"coefficients"} computes the coefficients at the requested values for s.
                \item Type \code{"nonzero"} returns a list of the indices of the nonzero coefficients for each value of s. }}
		
		\item{\dots}{Not used. Other arguments to predict. } }
\details{\code{s} is the new vector at which predictions are requested. If \code{s} is not in the lambda sequence used for fitting the model, the \code{predict} function will use linear interpolation to make predictions. The new values are interpolated using a fraction of predicted values from both left and right \code{lambda} indices. }
\value{The object returned depends on type.}


\author{Yi Yang and Hui Zou\cr
Maintainer: Yi Yang  <yiyang@umn.edu>}
\references{
Yang, Y. and Zou, H. (2012) 
\emph{A Cocktail Algorithm for Solving The Elastic Net Penalized Cox's Regression in High Dimensions Statistics and Its Interface.}\cr
\url{http://code.google.com/p/fastcox/}\cr

Friedman, J., Hastie, T. and Tibshirani, R. (2008)
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


\seealso{\code{\link{coef}} method}
\examples{
data(FHT)
m1<-cocktail(x=FHT$x,y=FHT$y,d=FHT$status,alpha=0.5)
predict(m1,type="nonzero")
predict(m1,s=c(0.01,0.005),type="coefficients")
predict(m1,newx=FHT$x[1:5,],type="response")
predict(m1,newx=FHT$x[1:3,],s=0.01,type="link")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
\keyword{regression}
