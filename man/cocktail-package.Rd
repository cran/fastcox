\name{fastcox-package}
\alias{fastcox-package}
\docType{package}
\title{Lasso and elastic-net penalized Cox's regression in high dimensions models using the cocktail algorithm}
\description{We introduce a cocktail algorithm, a good mixture of coordinate decent, the majorization-minimization principle and the strong rule, for computing the solution paths of the elastic net penalized Cox's proportional hazards model.}
\details{
\tabular{ll}{
Package:	\tab fastcox\cr
Type: 		\tab Package\cr
Version: 	\tab 1.0.0\cr
Date: 		\tab 2012-03-26\cr
Depends: 	\tab Matrix\cr
License: 	\tab GPL (version 2)\cr
URL: 		\tab http://code.google.com/p/fastcox/\cr
}
}

\author{Yi Yang and Hui Zou\cr
Maintainer: Yi Yang  <yiyang@umn.edu>}
\references{
Yang, Y. and Zou, H. (2012) 
\emph{A Cocktail Algorithm for Solving The Elastic Net Penalized Cox's Regression in High Dimensions Statistics and Its Interface.}\cr
\url{http://code.google.com/p/fastcox/}}

\examples{
data(FHT)
m1<-cocktail(x=FHT$x,y=FHT$y,d=FHT$status,alpha=0.5)
predict(m1,newx=FHT$x[1:5,],s=c(0.01,0.005))
predict(m1,type="nonzero")
plot(m1)
}

\keyword{package}

