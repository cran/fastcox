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
URL: 		\tab https://github.com/emeryyi/fastcox\cr
}
}

\author{Yi Yang and Hui Zou\cr
Maintainer: Yi Yang  <yi.yang6@mcgill.ca>}
\references{

Yang, Y. and Zou, H. (2013), 
"A Cocktail Algorithm for Solving The Elastic Net Penalized Cox's Regression in High Dimensions", \emph{Statistics and Its Interface,} 6:2, 167-173.\cr
\url{https://github.com/emeryyi/fastcox}\cr
}
	
\examples{
data(FHT)
m1<-cocktail(x=FHT$x,y=FHT$y,d=FHT$status,alpha=0.5)
predict(m1,type="nonzero")
plot(m1)
}

\keyword{package}

