######################################################################
## These functions are adapted/modified based on the functions from 
## the glmnet package:
##        Jerome Friedman, Trevor Hastie, Robert Tibshirani (2010). 
##        Regularization Paths for Generalized Linear Models via Coordinate Descent. 
##        Journal of Statistical Software, 33(1), 1-22. 
##        URL http://www.jstatsoft.org/v33/i01/.



plot.cocktail <- function(x, xvar = c("norm", "lambda"), color = FALSE, 
    label = FALSE, ...) {
    beta <- x$beta
    lambda <- x$lambda
    xvar <- match.arg(xvar)
    ##beta should be in 'dgCMatrix' format
    which <- nonzeroCoef(beta)
    beta <- as.matrix(beta[which, ])
    xvar <- match.arg(xvar)
    switch(xvar, norm = {
        index <- apply(abs(beta), 2, sum)
        iname <- "L1 Norm"
    }, lambda = {
        index <- log(lambda)
        iname <- "Log Lambda"
    })
    xlab <- iname
    ylab <- "Coefficients"
    dotlist <- list(...)
    type <- dotlist$type
    if (is.null(type)) {
        if (color == FALSE) 
            matplot(index, t(beta), lty = 1, xlab = xlab, ylab = ylab, type = "l", 
                pch = 200, col = gray.colors(12, start = 0.05, end = 0.7, gamma = 2.2), 
                ...) else matplot(index, t(beta), lty = 1, xlab = xlab, ylab = ylab, type = "l", 
            pch = 500, ...)
    } else matplot(index, t(beta), lty = 1, xlab = xlab, ylab = ylab, ...)
    if (label) {
        nnz <- length(which)
        xpos <- max(index)
        pos <- 4
        if (xvar == "lambda") {
            xpos <- min(index)
            pos <- 2
        }
        xpos <- rep(xpos, nnz)
        ypos <- beta[, ncol(beta)]
        text(xpos, ypos, paste(which), cex = 0.5, pos = pos)
    }
}




predict.cocktail <- function(object, newx, s = NULL, type = c("link", 
    "response", "coefficients", "nonzero"), ...) {
    type <- match.arg(type)
    if (missing(newx)) {
        if (!match(type, c("coefficients", "nonzero"), FALSE)) 
            stop("You need to supply a value for 'newx'")
    }
    nbeta <- object$beta
    if (!is.null(s)) {
        vnames <- dimnames(nbeta)[[1]]
        dimnames(nbeta) <- list(NULL, NULL)
        lambda <- object$lambda
        lamlist <- lambda.interp(lambda, s)
        nbeta <- nbeta[, lamlist$left, drop = FALSE] * lamlist$frac + nbeta[, 
            lamlist$right, drop = FALSE] * (1 - lamlist$frac)
        dimnames(nbeta) <- list(vnames, paste(seq(along = s)))
    }
    if (type == "coefficients") 
        return(nbeta)
    if (type == "nonzero") 
        return(nonzeroCoef(nbeta, bystep = TRUE))
    nfit <- as.matrix(newx %*% nbeta)
    switch(type, response = exp(nfit), link = nfit)
}




print.cocktail <- function(x, digits = max(3, getOption("digits") - 
    3), ...) {
    cat("\nCall: ", deparse(x$call), "\n\n")
    print(cbind(Df = x$df, Lambda = signif(x$lambda, digits)))
} 
