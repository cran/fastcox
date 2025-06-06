\name{cocktail}
\alias{cocktail}
\title{Fits the regularization paths for the elastic net penalized Cox's model}
\description{Fits a regularization path for the elastic net penalized Cox's model at a sequence of regularization parameters lambda.}
\usage{
cocktail(x,y,d,
	nlambda=100,
	lambda.min=ifelse(nobs<nvars,1e-2,1e-4),
	lambda=NULL, 
	alpha=1,
	pf=rep(1,nvars),
	exclude,
	dfmax=nvars+1,
	pmax=min(dfmax*1.2,nvars),
	standardize=TRUE,
	eps=1e-6,
	maxit=3e4)
}
\arguments{
		\item{x}{matrix of predictors, of dimension \eqn{N \times p}{N*p}; each row is an observation vector.}

		\item{y}{a survival time for Cox models. Currently tied failure times are not supported.}

		\item{d}{censor status with 1 if died and 0 if right censored.}

		\item{nlambda}{the number of \code{lambda} values - default is 100.}

		\item{lambda.min}{given as a fraction of \code{lambda.max} - the smallest value of \code{lambda} for which all coefficients are zero. The default depends on the relationship between \eqn{N} (the number of rows in the matrix of predictors) and \eqn{p} (the number of predictors). If \eqn{N > p}, the default is \code{0.0001},
		close to zero.  If \eqn{N<p}, the default is \code{0.01}.
		A very small value of \code{lambda.min} will lead to a saturated fit. It takes no effect if there is user-defined \code{lambda} sequence.} 

		\item{lambda}{a user supplied \code{lambda} sequence. Typically, by leaving this option unspecified users can have 
		the program compute its own \code{lambda} sequence based on
		\code{nlambda} and \code{lambda.min}. Supplying a value of
		\code{lambda} overrides this. It is better to supply
		a decreasing sequence of \code{lambda} values than a single (small) value, if not, the program will sort user-defined \code{lambda} sequence in decreasing order automatically.}

		\item{alpha}{The elasticnet mixing parameter, with \eqn{0 < \alpha <= 1}. See details.}

		\item{pf}{separate penalty weights can be applied to each coefficient of \eqn{\beta}{beta} to allow
		differential shrinkage. Can be 0 for some variables, which implies
		no shrinkage, and results in that variable always being included in the
		model. Default is 1 for all variables (and implicitly infinity for
		variables listed in \code{exclude}). See details.}

		\item{exclude}{indices of variables to be excluded from the
		model. Default is none. Equivalent to an infinite penalty factor.}

		\item{dfmax}{limit the maximum number of variables in the
		model. Useful for very large \eqn{p}, if a partial path is desired. Default is \eqn{p+1}.}

		\item{pmax}{limit the maximum number of variables ever to be nonzero. For example once \eqn{\beta} enters the model, no matter how many times it exits or re-enters model through the path, it will be counted only once. Default is \code{min(dfmax*1.2,p)}.}

		\item{standardize}{logical flag for variable standardization, prior to
		fitting the model sequence. If \code{TRUE}, x matrix is normalized such that sum squares of each column \eqn{\sum^N_{i=1}x_{ij}^2/N=1}{<Xj,Xj>/N=1}. Note that x is always centered (i.e. \eqn{\sum^N_{i=1}x_{ij}=0}{sum(Xj)=0}) no matter \code{standardize} is \code{TRUE} or \code{FALSE}. The coefficients are always returned on
		the original scale. Default is is \code{TRUE}.}

		\item{eps}{convergence threshold for coordinate majorization descent. Each inner
		coordinate majorization descent loop continues until the relative change in any
		coefficient (i.e. \eqn{\max_j|\beta_j^{new}-\beta_j^{old}|^2}{max(j)|beta_new[j]-beta_old[j]|^2}) is less than \code{eps}. Defaults value is \code{1e-6}.}

		\item{maxit}{maximum number of outer-loop iterations allowed at fixed lambda value. Default is 1e4. If models do not converge, consider increasing \code{maxit}.}
}

\details{
The algorithm estimates \eqn{\beta} based on observed data, through elastic net penalized log partial likelihood of Cox's model.
\deqn{\arg\min(-loglik(Data,\beta)+\lambda*P(\beta))}
It can compute estimates at a fine grid of values of \eqn{\lambda}{lambda}s in order to pick up a data-driven optimal \eqn{\lambda}{lambda} for fitting a 'best' final model. The penalty is a combination of l1 and l2 penalty:
\deqn{P(\beta)=(1-\alpha)/2||\beta||_2^2+\alpha||\beta||_1.} 
\code{alpha=1} is the lasso penalty.
For computing speed reason, if models are not converging or running slow, consider increasing \code{eps}, decreasing
\code{nlambda}, or increasing \code{lambda.min} before increasing
\code{maxit}.



\strong{FAQ:}

\bold{Question: }\dQuote{\emph{I am not sure how are we optimizing alpha. I can get optimal lambda for each value of alpha. But how do I select optimum alpha?}} 

\bold{Answer: } \code{cv.cocktail} only finds the optimal lambda given alpha fixed. So to
chose a good alpha you need to fit CV on a grid of alpha, say (0.1,0.3, 0.6, 0.9, 1) and let cv.cocktail choose the optimal lambda for each alpha, then you choose the (alpha, lambda) pair that corresponds to the lowest predicted deviance.

\bold{Question: }\dQuote{\emph{I understand your are referring to minimizing the quantity \code{cv.cocktail\$cvm}, the mean 'cross-validated error' to optimize alpha and lambda as you did in your implementation. However, I don't know what the equation of this error is and this error is not referred to in your paper either. Do you mind explaining what this is?
}}

\bold{Answer: } We first define the log partial-likelihood for the Cox model. Assume
\eqn{\hat{\beta}^{[k]}} is the estimate fitted on \eqn{k}-th fold, define the log partial likelihood function as 
\deqn{
L(Data,\hat{\beta}[k])=\sum_{s=1}^{S} x_{i_{s}}^{T}\hat{\beta}[k]-\log(\sum_{i\in R_{s}}\exp(x_{i}^{T}\hat{\beta}[k])).
}
Then the log partial-likelihood deviance of the \eqn{k}-th fold is defined
as
\deqn{
D[Data,k]=-2(L(Data,\hat{\beta}[k])).
}
We now define the measurement we actually use for cross validation:
it is the difference between the log partial-likelihood deviance evaluated
on the full dataset and that evaluated on the on the dataset with
\eqn{k}-th fold excluded. The cross-validated error is defined as
\deqn{
CV-ERR[k]=D(Data[full],k)-D(Data[k^{th}\,\,fold\,\,excluded],k).
}

}


\value{
An object with S3 class \code{\link{cocktail}}.
		\item{call}{the call that produced this object}
		\item{beta}{a \eqn{p*length(lambda)} matrix of coefficients, stored as a sparse matrix (\code{dgCMatrix} class, the standard class for sparse numeric matrices in the \code{Matrix} package.). To convert it into normal type matrix use \code{as.matrix()}.}
		\item{lambda}{the actual sequence of \code{lambda} values used}
		\item{df}{the number of nonzero coefficients for each value of
		\code{lambda}.}
		\item{dim}{dimension of coefficient matrix (ices)}
		\item{npasses}{total number of iterations (the most inner loop) summed over all lambda values}
		\item{jerr}{error flag, for warnings and errors, 0 if no error.}
}

\author{Yi Yang and Hui Zou\cr
Maintainer: Yi Yang  <yi.yang6@mcgill.ca>}
\references{
Yang, Y. and Zou, H. (2013), 
"A Cocktail Algorithm for Solving The Elastic Net Penalized Cox's Regression in High Dimensions", \emph{Statistics and Its Interface,} 6:2, 167-173.\cr
\url{https://github.com/archer-yang-lab/fastcox}\cr
}
\seealso{\code{plot.cocktail}}
\examples{
data(FHT)
m1<-cocktail(x=FHT$x,y=FHT$y,d=FHT$status,alpha=0.5)
predict(m1,type="nonzero")
plot(m1)
}
\keyword{models}
\keyword{regression}
