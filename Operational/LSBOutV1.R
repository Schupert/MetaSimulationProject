### Remove previous variables
rm(list = ls())

StartTimeTotal <- proc.time()

#### Libraries, set seed, set cores ----
library(data.table)
library(doParallel)
library(foreach)
library(doRNG)
library(copula)
library(compiler)
enableJIT(3)

set.seed(123)

# Number of cores for parallel
num.Cores <- detectCores() - 1
c1 <- makeCluster(num.Cores)

#### Declare variables ----

# Reps = number of repetitions of experiment
Reps = 30

# k = number of studies in series
Studies = c(3,5,10,30,50,100)

# subj = number of subjects in study
Subj <- list(as.integer(c(100,100)), as.integer(c(20,100)), as.integer(c(250, 1000)), as.numeric(c(4.7, 1.2)))

# sd = study level standard deviation
True.sd = sqrt(2)

# theta = population level log(OR) - this should be considered more purely on the log scale
theta = c(log(0.25), log(0.75), log(1), log(1.5), log(4))

# tau.sq = between studies variance (can be squared due to sqrt() in normal draw)
tau.sq = c(0, 0.01777778, 0.04, 3.04)

# Frequency of event averaged across 2 arms (before applying change due to theta) = EvFreq
EvFreq = c(0.1, 0.3, 0.5)

# controlProp = proportion of total sample in control arm
controlProp = 0.5

## Boundary of step function on p value, causing severity of publication bias
Severity.boundary <- c(0.05, 0.2)

# Set up strength of publication bias selection
Begg_a <- 0.5
Begg_b <- 3
Begg_c <- -0.3
Begg_sided <- 1

# Set up within study reporting bias
Tested.outcomes <- 5
Sd.split <- 0.8

# Size of per unit bias increase
Bias.multiple <- 0.9

#### Functions ----

### Log_Odds_Ratio function

Log_Odds_Ratio <- function(StudySize, Log_O_R, Heterogeneity, Control_Prop, mu){
  StudyLogOR <- rnorm(1, Log_O_R, sqrt(Heterogeneity))
  Group1Size <- as.integer(Control_Prop*StudySize)
  Group2Size <- Group1Size
  EventFreq <- log(mu / (1 - mu))
  Pic <- exp(EventFreq - 0.5*StudyLogOR) / (1 + exp(EventFreq - 0.5*StudyLogOR))
  Group1Out1 <- as.integer(rbinom(1, Group1Size, Pic))
  Group1Out2 <- as.integer(Group1Size - Group1Out1)
  Pit <- exp(EventFreq + 0.5*StudyLogOR)  / (1 + exp(EventFreq + 0.5*StudyLogOR))
  Group2Out1 <- as.integer(rbinom(1, Group2Size, Pit))
  Group2Out2 <- as.integer(Group2Size - Group2Out1)
  if (Group1Out1 == 0 | Group2Out1 == 0 | Group1Out2 == 0 | Group2Out2 == 0){
    Group1Out1 <- Group1Out1 + 0.5
    Group2Out1 <- Group2Out1 + 0.5
    Group1Out2 <- Group1Out2 + 0.5
    Group2Out2 <- Group2Out2 + 0.5
  }
  return(c(Group1Out1, Group1Out2, Group2Out1, Group2Out2, Group1Size, Group2Size))
}

### Function for multiple outcomes

LogOR_mult_out <- function(StudySize, Theta, Heterogeneity, Control_Prop, mu, frac, num.times){
  StudyLogOR <- rnorm(1, Theta, sqrt(Heterogeneity))
  Group1Size <- as.integer(Control_Prop*StudySize)
  Group2Size <- Group1Size
  EventFreq <- log(mu / (1 - mu))
  Pic <- exp(EventFreq - 0.5*StudyLogOR) / (1 + exp(EventFreq - 0.5*StudyLogOR))
  Pit <- exp(EventFreq + 0.5*StudyLogOR)  / (1 + exp(EventFreq + 0.5*StudyLogOR))
  z <- normalCopula(param = frac, dim = num.times)
  Z <- rCopula(Group1Size, z)
  ControlGroup <- qbinom(Z, size=1, prob=Pic)
  CGO1 <- apply(ControlGroup, 2, sum)
  CGO2 <- Group1Size - CGO1
  y <- normalCopula(param = frac, dim = num.times)
  Y <- rCopula(Group1Size, y)
  TreatmentGroup <- qbinom(Y,size=1,prob=Pit)
  TGO1 <- apply(TreatmentGroup, 2, sum)
  TGO2 <- Group2Size - TGO1
  o <- data.table(CGO1, CGO2, TGO1, TGO2)
  o[CGO1 == 0] <- o[CGO1 == 0] + 0.5
  o[CGO2 == 0] <- o[CGO2 == 0] + 0.5
  o[TGO1 == 0] <- o[TGO1 == 0] + 0.5
  o[TGO2 == 0] <- o[TGO2 == 0] + 0.5
  Study.est <- log((o$TGO1/o$TGO2) / (o$CGO1/o$CGO2))
  Study.se <- sqrt(1/o$TGO1 + 1/o$TGO2 + 1/o$CGO1 + 1/o$CGO2)
  Study.p.val <- pnorm(Study.est/Study.se)
  return(c(o[order(Study.p.val)], Group1Size))
}

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
  
  ###
  
  ### set control parameters for uniroot() and possibly replace with user-defined values
  con <- list(tol=.Machine$double.eps^0.25, maxiter=1000, verbose=FALSE)
  con[pmatch(names(control), names(con))] <- control
  
  ###
  
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
  
  ###
  
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
  
  ###
  
  res <- list(beta=beta, se=se, ci.lb=ci.lb, ci.ub=ci.ub,
              beta.rma=beta.rma, se.rma=se.rma, ci.lb.rma=ci.lb.rma, ci.ub.rma=ci.ub.rma,
              method="DL", method.rma=x$method, tau2=tau2, tau2.rma=x$tau2, digits=digits)
  
  class(res) <- "hc.rma.uni"
  return(res)
  
}

##### Parallel Loop Simulation ----

rng <- RNGseq(length(Subj)*length(theta), 4567)

StartTime <- proc.time()

registerDoParallel(c1)
set.seed(123)
r <- foreach (i = Subj, .combine=rbind, .packages = c("data.table", "copula"), 
              .export = c("Studies", "Subj", "True.sd",
                          "theta", "tau.sq", "controlProp", "Severity.boundary", "Begg_a", 
                          "Begg_b", "Begg_sided", "Tested.outcomes", "Sd.split",
                          "Bias.multiple", "Log_Odds_Ratio", "EvFreq", "LogOR_mult_out", "Begg_c")
) %:% foreach (k = theta, res = rng[(apply(sapply(Subj, function(vec) {i %in% vec}), 1, which.max)[1]-1)*length(theta) + 1:length(theta)], .combine=rbind, .packages = c("data.table", "copula"), 
               .export = c("Studies", "Subj", "True.sd",
                           "theta", "tau.sq", "controlProp", "Severity.boundary", "Begg_a", 
                           "Begg_b", "Begg_sided", "Tested.outcomes", "Sd.split",
                           "Bias.multiple", "Log_Odds_Ratio", "EvFreq", "LogOR_mult_out", "Begg_c")
) %dopar% {
  
  #set rng seed
  rngtools::setRNG(res)
  
  ID = length(tau.sq) * length(EvFreq) * Reps * sum(Studies) 
  
  LogOR.Simulation <- data.table(
    Unique_ID = integer(length = ID),
    Study_G1O1 = numeric(length = ID),
    Study_G2O1 = numeric(length = ID),
    Study_n = integer(length = ID),
        Study_rejectedMeans = list(length = ID),
        Study_rejectedSDs = list(length = ID)
    #     Study_Number.of.biases = integer(length = ID)
  )
  
  dummy.counter <- 1
  
  for (l in tau.sq){
    
    for (j in EvFreq){
      
      for (n in Studies){
        
        for (o in 1:n){
          
          for (m in 1:Reps){
            
            counter <- as.integer((apply(sapply(Subj, function(vec) {i %in% vec}), 1, which.max)[1]-1) * length(EvFreq) * length(theta) * length(tau.sq) * sum(Studies) * Reps + 
                                    (match(k, theta)-1) * length(tau.sq) * length(EvFreq) * sum(Studies) * Reps +
                                    (match(l, tau.sq)-1) * length(EvFreq) * sum(Studies) * Reps +
                                    (match(j, EvFreq)-1) * sum(Studies) * Reps +
                                    (sum(Studies[0:(match(n, Studies)-1)]) + o -1) * Reps +
                                    m
            )
            
            #Select sample size
            if (is.integer(i[1]) == TRUE){
              Study_patientnumber <- as.integer( (runif(1, sqrt(i[1]), sqrt(i[2])))^2 )
            } else {
              Study_patientnumber <- round(rlnorm(1, meanlog = i[1], sdlog = i[2]) + 2)
            }
            
            # Different function for simulation creates multiple correlated outcomes, ordered by p value
            x <- LogOR_mult_out(Study_patientnumber, k, l, controlProp, j, Sd.split, Tested.outcomes)
            x.N <- 2*x[[5]]
            
            LogOR.Simulation[dummy.counter, `:=` (Unique_ID = counter,
              Study_G1O1 = x[[1]][1], 
              Study_G2O1 = x[[3]][1], 
              Study_n = x.N,
              Study_rej_G1O1 = list(x[[1]][-1]),
              Study_rej_G2O1 = list(x[[3]][-1])
            )]
            
            dummy.counter <- dummy.counter + 1
            
          }
          
        }
      }
    }
  }
  LogOR.Simulation
}

#### Need to then sort final table and add values for rep number, rep subj, rep theta, rep tau2, rep numstudies
LogOR.Simulation <- r[order(Unique_ID)]

# ID = total number of data points required, also used as an ID number. 
ID =  length(Subj) * length(theta) * length(tau.sq) * length(EvFreq) * Reps * sum(Studies)

LogOR.Simulation$Rep_Number =  rep(1:Reps, times = ID/Reps)
intermediate <- integer()
for (i in Studies){intermediate <- append(intermediate, rep(i, times = i*Reps))}
LogOR.Simulation$Rep_NumStudies = rep(intermediate, times = ID/(Reps*sum(Studies)))
rm(intermediate)
LogOR.Simulation$Rep_ev_freq = rep(rep(EvFreq, each = Reps * sum(Studies)), times = ID/(Reps*sum(Studies)*length(EvFreq)))
LogOR.Simulation$Rep_tau.sq = rep(rep(tau.sq, each = Reps * sum(Studies)*length(EvFreq)), times = ID/(Reps*sum(Studies)*length(tau.sq)*length(EvFreq)))
LogOR.Simulation$Rep_theta = rep( rep(theta, each = Reps * sum(Studies) * length(tau.sq)*length(EvFreq)), times = length(Subj))

### Create keyable vector for Subj
Subj2 <- c(100, 20, 250, 4.7)
LogOR.Simulation$Rep_Subj = rep(Subj2, each = ID / length(Subj))

### Mean and sd

LogOR.Simulation[, Study_estimate := ifelse(LogOR.Simulation$Study_G1O1 %% 1 == 0.5, 
                                            log( (Study_G2O1 / ((Study_n/2 + 1) - Study_G2O1) ) / (Study_G1O1 / ((Study_n/2 + 1) - Study_G1O1) ) ),
                                            log( (Study_G2O1 / ((Study_n/2) - Study_G2O1) ) / (Study_G1O1 / ((Study_n/2) - Study_G1O1) ) ) )
                 ]

LogOR.Simulation[, Study_sd := ifelse(LogOR.Simulation$Study_G1O1 %% 1 == 0.5, 
                                      sqrt(1/Study_G1O1 + 1/((Study_n/2 + 1) - Study_G1O1) + 1/Study_G2O1 + 1/((Study_n/2 + 1) - Study_G2O1)),
                                      sqrt(1/Study_G1O1 + 1/(Study_n/2 - Study_G1O1) + 1/Study_G2O1 + 1/(Study_n/2 - Study_G2O1)) )
                 ]
TimeTakenSim <- proc.time() - StartTime

stopCluster(c1)

#write.csv(LogOR.Simulation, file = "LSBOutV1.csv")

df.LogOR.Simulation <- as.data.frame(LogOR.Simulation)
df.LogOR.Simulation$Study_rejectedMeans <- vapply(df.LogOR.Simulation$Study_rejectedMeans, paste, collapse = ", ", character(1L))
df.LogOR.Simulation$Study_rejectedSDs <- vapply(df.LogOR.Simulation$Study_rejectedSDs, paste, collapse = ", ", character(1L))

write.csv(df.LogOR.Simulation, file = "LSBOutV1.csv")
rm(r)
rm(df.LogOR.Simulation)

##### Parallel Loop Analysis ----

library(metafor)

Subj <- c(100, 20, 250, 4.7)

#LogOR.Simulation <- read.csv("LSBStepV1.csv")
#LogOR.Simulation <- data.table(LogOR.Simulation)
setkey(LogOR.Simulation, "Rep_Number", "Rep_Subj", "Rep_theta", "Rep_tau.sq", "Rep_NumStudies", "Rep_ev_freq")


StartTime <- proc.time()
c1 <- makeCluster(num.Cores)
registerDoParallel(c1)

r <- foreach (i = Subj, .combine=rbind, .packages = c("data.table", "metafor"), 
              .export = c("Studies", "Subj", "True.sd", "Reps",
                          "theta", "tau.sq", "EvFreq", "LogOR.Simulation", "anyNA", ".psort", "mod.hc")
) %:% foreach (k = theta, .combine=rbind, .packages = c("data.table", "metafor"), 
               .export = c("Studies", "Subj", "True.sd", "Reps",
                           "theta", "tau.sq", "EvFreq", "LogOR.Simulation", "anyNA", ".psort", "mod.hc")
) %dopar% {
  # ID different for analysis
  ID = length(tau.sq) * Reps * length(Studies) * length(EvFreq)
  
  LogOR.Sim.Results <- data.table(
    Unique_ID = integer(length = ID),
    FE_Estimate = numeric(length = ID),
    FE_se = numeric(length = ID),
    REML_Estimate = numeric(length = ID),
    REML_se = numeric(length = ID),
    REML_tau2 = numeric(length = ID),
    DL_Estimate = numeric(length = ID),
    DL_se = numeric(length = ID),
    DL_tau2 = numeric(length = ID),
    DL_I2 = numeric(length = ID),
    HC_DL_se = numeric(length = ID),
    HC_DL_CIlb = numeric(length = ID),
    HC_DL_CIub = numeric(length = ID),
    HC_REML_se = numeric(length = ID),
    HC_REML_CIlb = numeric(length = ID),
    HC_REML_CIub = numeric(length = ID),
    KH_REML_CIlb = numeric(length = ID),
    KH_REML_CIub = numeric(length = ID),
    KH_REML_se = numeric(length = ID),
    KH_DL_CIlb = numeric(length = ID),
    KH_DL_CIub = numeric(length = ID),
    KH_DL_se = numeric(length = ID),
    Moreno_Estimate = numeric(length = ID),
    Moreno_se = numeric(length = ID),
    Mult_se = numeric(length = ID),
    Num_exc = integer(length = ID)
  )
  dummy.counter <- 1
  
  for (l in tau.sq){
    
    for (j in EvFreq){
      
      for (n in Studies){
        
        for (m in 1:Reps){
          
          
          counter <- as.integer((match(i, Subj)-1) * length(EvFreq) * length(theta) * length(tau.sq) * length(Studies) * Reps + 
                                  (match(k, theta)-1) * length(tau.sq) * length(EvFreq) * length(Studies) * Reps +
                                  (match(l, tau.sq)-1) * length(EvFreq) * length(Studies) * Reps +
                                  (match(j, EvFreq)-1) * length(Studies) * Reps +
                                  (match(n, Studies)-1) * Reps +
                                  m
          )
          
          
          ### Temporary data.table
          temp.data <- LogOR.Simulation[J(m, i, k, l, n, j)]
          
          temp.data <- temp.data[ !((Study_G1O1 == 0.5 & Study_G2O1 == 0.5) | ( ((Study_n/2 + 1) - Study_G1O1) == 0.5 & ((Study_n/2 + 1) - Study_G2O1) == 0.5 ))]
          Excluded.studies <- n - dim(temp.data)[1]
          
          if (dim(temp.data)[1] < 2){
            LogOR.Sim.Results[dummy.counter, `:=` (Unique_ID = counter,
                                                   FE_Estimate = NA,
                                                   FE_se = NA,
                                                   REML_Estimate = NA,
                                                   REML_se = NA,
                                                   REML_tau2 = NA,
                                                   DL_Estimate = NA,
                                                   DL_se = NA,
                                                   DL_tau2 = NA,
                                                   DL_I2 = NA,
                                                   HC_DL_se = NA,
                                                   HC_DL_CIlb = NA,
                                                   HC_DL_CIub = NA,
                                                   HC_REML_se = NA,
                                                   HC_REML_CIlb = NA,
                                                   HC_REML_CIub = NA,
                                                   KH_REML_CIlb = NA,
                                                   KH_REML_CIub = NA,
                                                   KH_REML_se = NA,
                                                   KH_DL_CIlb = NA,
                                                   KH_DL_CIub = NA,
                                                   KH_DL_se = NA,
                                                   Moreno_Estimate = NA,
                                                   Moreno_se = NA,
                                                   Mult_se = NA,
                                                   Num_exc = Excluded.studies
            )]
            dummy.counter <- dummy.counter + 1
          } else {
            #Fixed and random effects
            ma.fe <- tryCatch({
              rma.uni(temp.data$Study_estimate, temp.data$Study_sd^2 , method = "FE")
            },
            error = function(e){
              return(list(b = NA,  se = NA))
            },
            warning = function(w){
              return(list(list(b = NA,  se = NA)))
            }
            )
            
            ma.reml <- tryCatch({
              rma.uni(temp.data$Study_estimate, temp.data$Study_sd^2  , method = "REML", control = list(stepadj = 0.5))
            },
            error = function(e){
              return(list(b = NA, tau2 = NA, se = NA))
            },
            warning = function(w){
              return(list(list(b = NA, tau2 = NA, se = NA)))
            }
            )
            
            ma.DL <- tryCatch({
              rma.uni(temp.data$Study_estimate, temp.data$Study_sd^2  , method = "DL")
            },error = function(e){
              return(list(b = NA, tau2 = NA, se = NA, I2 = NA))
            },warning = function(w){
              return(list(list(b = NA, tau2 = NA, se = NA, I2 = NA)))
            })
            
            # Henmi & Copas
            
            ma.hc.DL <- tryCatch({
              hc(ma.DL)
            },error = function(e){
              return(list(se = NA, ci.lb = NA, ci.ub = NA))
            },warning = function(w){
              return(list(list(se = NA, ci.lb = NA, ci.ub = NA)))
            })
            
            ma.hc.REML <- tryCatch({
              mod.hc(ma.reml, tau2est = ma.reml$tau2)
            },error = function(e){
              return(list(se = NA, ci.lb = NA, ci.ub = NA))
            },warning = function(w){
              return(list(list(se = NA, ci.lb = NA, ci.ub = NA)))
            })
            
            # Knapp Hartung
            
            ma.reml.kh <- tryCatch({
              rma.uni(temp.data$Study_estimate, temp.data$Study_sd^2  , method = "REML", knha = TRUE, control = list(stepadj = 0.5))
            },
            error = function(e){
              return(list(b = NA, tau2 = NA, se = NA, ci.lb = NA, ci.ub = NA))
            },
            warning = function(w){
              return(list(b = NA, tau2 = NA, se = NA, ci.lb = NA, ci.ub = NA))
            }
            )
            
            ma.DL.kh <- tryCatch({rma.uni(temp.data$Study_estimate, temp.data$Study_sd^2  , method = "DL", knha = TRUE)}
                                 , error = function(e){return(list(se = NA, ci.lb = NA, ci.ub = NA))
                                 })
            
            ## Doi
            # estimate is equal to fixed effect, as are weights
            doi.var.DL <- tryCatch(sum( ( as.vector(weights(ma.fe, type = "diagonal")/100)^2 ) * (temp.data$Study_sd^2 + ma.DL$tau2) ), error=function(err) NA)
            doi.var.REML <- tryCatch(sum( ( as.vector(weights(ma.fe, type = "diagonal")/100)^2 ) * (temp.data$Study_sd^2 + ma.reml$tau2) ), error=function(err) NA)
            
            ## Moreno (?D-var) - not exactly clear which implementation is being used is likely equation 2a
            moreno.est <- tryCatch({ma.moren <- regtest(ma.fe , predictor = "vi", model = "lm")
            ma.moren$fit[[5]][1]
            }, error=function(err) NA)
            
            ## Mawdesley
            
            ma.mult <- tryCatch({
              mawd.lm <- lm(temp.data$Study_estimate ~ 1, weights = 1/(temp.data$Study_sd^2))
              sm.mawd.lm <- summary(mawd.lm)
              ifelse(mean(sm.mawd.lm$residuals^2) < 1, phi.est  <- 1, phi.est <- mean(sm.mawd.lm$residuals^2))
              rma.uni(temp.data$Study_estimate, temp.data$Study_sd * sqrt(phi.est) , method = "FE")}
              , error = function(e){return(list(se = NA))
              })
            
            
            LogOR.Sim.Results[dummy.counter, `:=` (Unique_ID = counter,
                                                   FE_Estimate = ma.fe[[1]][1],
                                                   FE_se = ma.fe$se,
                                                   REML_Estimate = ma.reml$b[1],
                                                   REML_se = ma.reml$se,
                                                   REML_tau2 = ma.reml$tau2,
                                                   DL_Estimate = ma.DL[[1]][1],
                                                   DL_se = ma.DL$se,
                                                   DL_tau2 = ma.DL$tau2,
                                                   DL_I2 = ma.DL$I2,
                                                   HC_DL_se = ma.hc.DL$se,
                                                   HC_DL_CIlb = ma.hc.DL$ci.lb,
                                                   HC_DL_CIub = ma.hc.DL$ci.ub,
                                                   HC_REML_se = ma.hc.REML$se,
                                                   HC_REML_CIlb = ma.hc.REML$ci.lb,
                                                   HC_REML_CIub = ma.hc.REML$ci.ub,
                                                   KH_REML_CIlb = ma.reml.kh$ci.lb,
                                                   KH_REML_CIub = ma.reml.kh$ci.ub,
                                                   KH_REML_se = ma.reml.kh$se,
                                                   KH_DL_CIlb = ma.DL.kh$ci.lb,
                                                   KH_DL_CIub = ma.DL.kh$ci.ub,
                                                   KH_DL_se = ma.DL.kh$se,
                                                   IVHet_DL_var = doi.var.DL,
                                                   IVHet_REML_var = doi.var.REML,
                                                   Moreno_Estimate = moreno.est[1],
                                                   Moreno_se = moreno.est[2],
                                                   Mult_se = ma.mult$se,
                                                   Num_exc = Excluded.studies
            )]
            
            dummy.counter <- dummy.counter + 1
            
          }
        }
      }
    }
  }
  LogOR.Sim.Results
}

LogOR.Sim.Results <- r[order(Unique_ID)]

##### Need to re append values - specific to analysis
ID =  length(Subj) * length(theta) * length(tau.sq) * length(EvFreq) * Reps * length(Studies)

LogOR.Sim.Results$Rep_Number =  rep(1:Reps, times = ID/Reps)
LogOR.Sim.Results$Rep_NumStudies = rep(rep(Studies, each = Reps), times = ID/(Reps*length(Studies)))
LogOR.Sim.Results$Rep_ev_freq = rep(rep(EvFreq, each = Reps * length(Studies)), times = ID/(Reps*length(Studies)*length(EvFreq)))
LogOR.Sim.Results$Rep_tau.sq = rep(rep(tau.sq, each = Reps * length(Studies)*length(EvFreq)), times = ID/(Reps*length(Studies)*length(tau.sq)*length(EvFreq)))
LogOR.Sim.Results$Rep_theta = rep( rep(theta, each = Reps * length(Studies) * length(tau.sq)*length(EvFreq)), times = length(Subj))
LogOR.Sim.Results$Rep_Subj = rep(Subj, each = ID / length(Subj))

TimeTakenAn <- proc.time() - StartTime

write.csv(LogOR.Sim.Results, file = "LSBOutV1An.csv")

sum(is.na(LogOR.Sim.Results))
head(LogOR.Sim.Results[is.na(REML_Estimate)])

#setkey(LogOR.Simulation, "Rep_Number", "Rep_Subj", "Rep_theta", "Rep_tau.sq", "Rep_NumStudies", "Rep_ev_freq")

#asdf1 <- LogOR.Simulation[J(5, 20, log(0.75), 0.04, 3, 0.1)]

mean(LogOR.Sim.Results[Rep_theta == 0 & Rep_NumStudies == 100 & Rep_Subj == 100 & Rep_tau.sq == 3.04]$DL_Estimate)
mean(LogOR.Sim.Results[Rep_theta == 0 & Rep_NumStudies == 100 & Rep_Subj == 100 & Rep_tau.sq == 0.04]$DL_Estimate)
mean(LogOR.Sim.Results[Rep_theta == 0 & Rep_NumStudies == 100 & Rep_Subj == 100 & Rep_tau.sq == 0.01777778]$DL_Estimate)

#### Timings ----

TimeTakenTotal <- proc.time() - StartTimeTotal

TimeTakenSim
TimeTakenAn
TimeTakenTotal