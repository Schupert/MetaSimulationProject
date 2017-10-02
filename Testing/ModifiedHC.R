anyNA <- function(x) {
  i <- 1
  repeat {
    if (is.na(x[i])) return(TRUE)
    i <- i + 1
    if (i > length(x)) return(FALSE)
  }
}

.psort <- function(x,y) {
  
  ### t(apply(xy, 1, sort)) would be okay, but problematic if there are NAs;
  ### either they are removed completely (na.last=NA) or they are always put
  ### first/last (na.last=FALSE/TRUE); but we just want to leave the NAs in
  ### their position!
  
  if (is.null(x) || length(x) == 0) ### need to catch this
    return(NULL)
  
  if (missing(y)) {
    if (is.matrix(x)) {
      xy <- x
    } else {
      xy <- rbind(x) ### in case x is just a vector
    }
  } else {
    xy <- cbind(x,y)
  }
  
  n <- nrow(xy)
  
  for (i in seq_len(n)) {
    if (anyNA(xy[i,]))
      next
    xy[i,] <- sort(xy[i,])
  }
  
  colnames(xy) <- NULL
  
  return(xy)
  
}

mod.hc <- function(object, digits, transf, targs, control, tau2est, ...) {
  
  if (!inherits(object, "rma.uni"))
    stop("Argument 'object' must be an object of class \"rma.uni\".")
  
  if (inherits(object, "rma.ls"))
    stop("Method not yet implemented for objects of class \"rma.ls\". Sorry!")
  
  x <- object
  
  if (!x$int.only)
    stop("Method only applicable for models without moderators.")
  
  if (missing(digits))
    digits <- x$digits
  
  if (missing(transf))
    transf <- FALSE
  
  if (missing(targs))
    targs <- NULL
  
  yi <- x$yi
  vi <- x$vi
  k  <- length(yi)
  
  if (k == 1)
    stop("Stopped because k = 1.")
  
  if (!x$allvipos)
    stop("Cannot use method when one or more sampling variances are non-positive.")
  
  level <- ifelse(x$level > 1, (100-x$level)/100, ifelse(x$level > .5, 1-x$level, x$level))
  
  if (missing(control))
    control <- list()
  
  #########################################################################
  
  ### set control parameters for uniroot() and possibly replace with user-defined values
  con <- list(tol=.Machine$double.eps^0.25, maxiter=1000, verbose=FALSE)
  con[pmatch(names(control), names(con))] <- control
  
  #########################################################################
  
  ### original code by Henmi & Copas (2012), modified by Michael Dewey, small adjustments
  ### for consistency with other functions in the metafor package by Wolfgang Viechtbauer
  
  wi <- 1/vi ### fixed effects weights
  
  W1 <- sum(wi)
  W2 <- sum(wi^2) / W1
  W3 <- sum(wi^3) / W1
  W4 <- sum(wi^4) / W1
  
  ### fixed-effects estimate of theta
  beta <- sum(wi*yi) / W1
  
  ### Q statistic
  Q <- sum(wi * ((yi - beta)^2))
  
  ### DL estimate of tau^2
  ###### Modified here to take REML tau2
  tau2 <- max(0, tau2est)
  
  vb  <- (tau2 * W2 + 1) / W1 ### estimated Var of b
  se  <- sqrt(vb)             ### estimated SE of b
  VR  <- 1 + tau2 * W2        ### estimated Var of R
  SDR <- sqrt(VR)             ### estimated SD of R
  
  ### conditional mean of Q given R=r
  EQ <- function(r)
    (k - 1) + tau2 * (W1 - W2) + (tau2^2)*((1/VR^2) * (r^2) - 1/VR) * (W3 - W2^2)
  
  ### conditional variance of Q given R=r
  VQ <- function(r) {
    rsq <- r^2
    recipvr2 <- 1 / VR^2
    2 * (k - 1) + 4 * tau2 * (W1 - W2) +
      2 * tau2^2 * (W1*W2 - 2*W3 + W2^2) +
      4 * tau2^2 * (recipvr2 * rsq - 1/VR) * (W3 - W2^2) +
      4 * tau2^3 * (recipvr2 * rsq - 1/VR) * (W4 - 2*W2*W3 + W2^3) +
      2 * tau2^4 * (recipvr2 - 2 * (1/VR^3) * rsq) * (W3 - W2^2)^2
  }
  
  scale <- function(r){VQ(r)/EQ(r)}   ### scale parameter of the gamma distribution
  shape <- function(r){EQ(r)^2/VQ(r)} ### shape parameter of the gamma distribution
  
  ### inverse of f
  finv <- function(f)
    (W1/W2 - 1) * ((f^2) - 1) + (k - 1)
  
  ### equation to be solved
  eqn <- function(x) {
    integrand <- function(r) {
      pgamma(finv(r/x), scale=scale(SDR*r), shape=shape(SDR*r))*dnorm(r)
    }
    integral <- integrate(integrand, lower=x, upper=Inf)$value
    val <- integral - level / 2
    #cat(val, "\n")
    val
  }
  
  t0 <- try(uniroot(eqn, lower=0, upper=2, tol=con$tol, maxiter=con$maxiter))
  
  if (inherits(t0, "try-error"))
    stop("Error in uniroot().")
  
  t0 <- t0$root
  u0 <- SDR * t0 ### (approximate) percentage point for the distribution of U
  
  #########################################################################
  
  ci.lb <- beta - u0 * se ### lower CI bound
  ci.ub <- beta + u0 * se ### upper CI bound
  
  beta.rma  <- x$beta
  se.rma    <- x$se
  ci.lb.rma <- x$ci.lb
  ci.ub.rma <- x$ci.ub
  
  ### if requested, apply transformation to yi's and CI bounds
  
  if (is.function(transf)) {
    if (is.null(targs)) {
      beta      <- sapply(beta, transf)
      beta.rma  <- sapply(beta.rma, transf)
      se        <- NA
      se.rma    <- NA
      ci.lb     <- sapply(ci.lb, transf)
      ci.ub     <- sapply(ci.ub, transf)
      ci.lb.rma <- sapply(ci.lb.rma, transf)
      ci.ub.rma <- sapply(ci.ub.rma, transf)
    } else {
      beta      <- sapply(beta, transf, targs)
      beta.rma  <- sapply(beta.rma, transf, targs)
      se        <- NA
      se.rma    <- NA
      ci.lb     <- sapply(ci.lb, transf, targs)
      ci.ub     <- sapply(ci.ub, transf, targs)
      ci.lb.rma <- sapply(ci.lb.rma, transf, targs)
      ci.ub.rma <- sapply(ci.ub.rma, transf, targs)
    }
  }
  
  ### make sure order of intervals is always increasing
  
  tmp <- .psort(ci.lb, ci.ub)
  ci.lb <- tmp[,1]
  ci.ub <- tmp[,2]
  
  tmp <- .psort(ci.lb.rma, ci.ub.rma)
  ci.lb.rma <- tmp[,1]
  ci.ub.rma <- tmp[,2]
  
  #########################################################################
  
  res <- list(beta=beta, se=se, ci.lb=ci.lb, ci.ub=ci.ub,
              beta.rma=beta.rma, se.rma=se.rma, ci.lb.rma=ci.lb.rma, ci.ub.rma=ci.ub.rma,
              method="DL", method.rma=x$method, tau2=tau2, tau2.rma=x$tau2, digits=digits)
  
  class(res) <- "hc.rma.uni"
  return(res)
  
}

temp.data <- Normal.Simulation[J(m, i, k, l, n)]
ma.reml <- rma.uni(temp.data$Study_estimate, temp.data$Study_sd^2  , method = "REML", control = list(stepadj = 0.5))
ma.DL <- rma.uni(temp.data$Study_estimate, temp.data$Study_sd^2  , method = "DL")
hc.DL <- hc(ma.DL)
hc.reml <- mod.hc(ma.reml, tau2est = ma.reml$tau2)
ma.fe <- rma.uni(temp.data$Study_estimate, temp.data$Study_sd^2  , method = "FE")
doi <- sum( ( as.vector(weights(ma.fe, type = "diagonal")/100)^2 ) * (temp.data$Study_sd^2 + ma.DL$tau2) )
mawd.lm <- lm(temp.data$Study_estimate ~ 1, weights = 1/(temp.data$Study_sd^2))
sm.mawd.lm <- summary(mawd.lm)
ifelse(mean(sm.mawd.lm$residuals^2) < 1, phi.est  <- 1, phi.est <- mean(sm.mawd.lm$residuals^2))
rma.uni(temp.data$Study_estimate, temp.data$Study_sd * sqrt(phi.est) , method = "FE")
phi.est
ma.DL$H2
