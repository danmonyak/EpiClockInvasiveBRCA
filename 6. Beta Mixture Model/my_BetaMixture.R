## betamix.R
## Fit a mixture of K beta distributions

#' Compute the (negative) log-likelihood
#'
#' @param param A vector of parameters of length 2K
#' @param vec A data vector of length N
#' @param z A matrix of size NxK of nonegative real numbers, with the sum of each row equal to 1.
NegBetaLogLike <- function(param, vec, z) {
  N <- length(vec)
  K <- ncol(z)
  W <- options("warn")
  options(warn = -1)
  on.exit(options(warn = W$warn))
  temp <- sapply(1:K, function(I) {
    a <- param[2*I - 1]
    b <- param[2*I]
    L <- lbeta(a, b)
    ifelse(is.infinite(L),
           .Machine$double.xmax,
           sum(z[, I] * (-L + (a - 1) * log(vec) + (b - 1) * log(1-vec))))
  })
  - sum(temp)
}


setClass("BetaMixture",
         slots = c(datavec = "numeric", # input data
                   phi = "numeric",     # mixture weights
                   mle = "matrix",      # fitted beta parameters
                   Z = "matrix",        # fitted latent mixture parameters
                   loglike = "numeric", # log likelihood for c(Z, mle)
                   converged = "logical"
         ))

setValidity("BetaMixture", function(object) {
  N <- length(object@datavec)
  K = length(object@phi)
  all(dim(object@mle) == c(2, K)) &
    all(dim(object@Z) == c(N, K)) &
    isTRUE(all.equal(sum(object@phi), 1.0)) &
    all(abs(1 -apply(object@Z, 1, sum)) < 1E-10) &
    length(object@converged) == 1
})

setMethod("summary", "BetaMixture", function(object, ...) {
  cat("An object of the 'BetaMixture' class with ",
      length(object@phi), " components using ",
      length(object@datavec),  " observations.\n", sep = "")
  cat("Mixing parameters (weights):\n")
  print(object@phi)
  cat("Beta component parameters:\n")
  print(object@mle)
  cat("The model did", ifelse(object@converged, " ", " not "),
      "converge and has a log-likelihood of ",
      object@loglike, "\n", sep = "")
  #  invisible(object)
})

setMethod("hist", "BetaMixture", function(x, mixcols = 1:7, ...) {
  hist(x@datavec, prob = TRUE, ...)
  while(length(mixcols) < ncol(x@mle)) mixcols <- c(mixcols, mixcols)
  xv <- seq(0, 1, length=502)[1:501]
  for (J in 1:ncol(x@mle)) {
    y  <- x@phi[J] * dbeta(xv, x@mle[1, J], x@mle[2, J])
    lines(xv, y, lwd=2, col = mixcols[J])
  }
  invisible(x)
})


#' Run en expectation-maximization algorithm to fit a mixture of beta distributions
#'
#' @param datavec The observed vector of data points (between 0 and 1)
#' @param K The number of mixture components to fit
#' @param forever The maximum number of iterations of the algorithm
#' @param epsilon The minimum change of parameters that must be observed in order to continue iterating
#' @param debug A local valuecontrolling whethe to print out intermediate steps
#'
#' @importFrom stats dbeta nlm
#' 
#' @export
BetaMixture_my <- function(datavec, K = 2, forever = 100,
                        epsilon = 0.001, relative = 0.001, debug = FALSE) {
  ### intitalize cluster assignments
  n <- length(datavec)
  starter <- cut(datavec, breaks = K, labels = FALSE, include.lowest = TRUE)

  #  starter <- sample(K, n, replace = TRUE)
  Z <- matrix(0, nrow = n, ncol = K)
  for (I in 1:length(starter)) Z[I, starter[I]]  <-  1
  ### Set up loop control
  lastlike <- -10^10
  currlike <- 0
  delta <- abs(currlike - lastlike)
  reld <- 2*relative
  count <- 0
  mle <- rep(1, 2*K)
  while ((count < forever) & (delta > epsilon) & (reld > relative)) {
    count <- count + 1
    if(debug) print(c(count, abs(lastlike - currlike), mle))
    ### M-step
    maxlike <- nlm(NegBetaLogLike, mle, vec = datavec, z = Z,
                   stepmax = 10000, print.level = 0)
    mle <- maxlike$estimate
    lastlike <- currlike
    currlike <- maxlike$minimum
    if (any(mle < 0)) stop("One of the mle's is negative!")
    ## E-step
    phi <- apply(Z, 2, mean)
    mix <- sapply(1:K, function(I) {
      phi[I] * dbeta(datavec, mle[2*I -1], mle[2*I])
    }) # should be a matrix with n rows and K columns

    Z <-  sweep(mix, 1, apply(mix, 1, sum), "/")
    delta <- abs(lastlike - currlike)
    reld <- abs(delta/currlike)
  }

  new("BetaMixture",
      datavec = datavec,
      phi = phi,
      mle = matrix(mle, nrow = 2),
      Z = Z,
      loglike = -currlike,
      converged = count < forever)
}



if(FALSE) {
  set.seed(73892)
  datavec <- c(rbeta(130, 1, 4),
               rbeta(170, 7, 4))
  model <- BetaMixture(datavec)
  summary(model)
  hist(model, breaks=35)
  
  datavec <- c(rbeta(130, 1, 4),
               rbeta(170, 7, 4),
               rbeta(200, 8, 8))
  
  model2 <- BetaMixture(datavec, K=2)
  summary(model2)
  hist(model2, breaks=35)
  
  model3 <- BetaMixture(datavec, K=3)
  summary(model3)
  hist(model3, breaks=35)
  
  model4 <- BetaMixture(datavec, K=4)
  summary(model4)
  hist(model4, breaks=35)
  
  model5 <- BetaMixture(datavec, K=5)
  summary(model5)
  hist(model5, breaks=35)
  
  
  fubar <- function(Z) {
    dex <- apply(Z, 1, which.max)
    W <- 0*Z
    for (I in 1:length(dex)) W[I, dex[I]]  <-  1
    W
  }
  mz2 <- fubar(model2$Z)
  
  nll <- c(NegBetaLogLike(model2$mle, datavec, model2$Z),
           NegBetaLogLike(model3$mle, datavec, model3$Z),
           NegBetaLogLike(model4$mle, datavec, model4$Z),
           NegBetaLogLike(model5$mle, datavec, model5$Z))
  
  f  <-  function(M) M*(2:5) + nll
  barf <- t(sapply(X <- seq(70, 120, 1), f))
  plot(X, barf[, 1], ylim=range(barf), type = "n")
  for (I in 1:4) lines(X, barf[,I], col = I, lwd=2, pch = 16, type = "b")
  
  
  g <- function(M) M*(2:5)*log(500) + nll
  barf <- t(sapply(X <- seq(10, 20, 1), g))
  plot(X, barf[, 1], ylim=range(barf), type = "n")
  for (I in 1:4) lines(X, barf[,I], col = I, lwd=2, pch=16, type = "b")
  
}