### Set seed
set.seed(1)

### Import required libraries
require(boot)
require(data.table)

### Get data and set keyframes
#r <- read.csv("NormalSimulation1Analysis.csv")
#r <- data.frame(r)
#setkey(r, "Rep_Number", "Rep_Subj", "Rep_theta", "Rep_tau.sq", "Rep_NumStudies")

###### Select specific set of results
Working_copy1 <- r[Rep_Subj == 1 & Rep_theta == 1 & Rep_tau.sq ==3 & Rep_NumStudies == 3]$FE_Estimate

## Can increase multiplier when have more reps to estimate with up to available data/3
Target_MCE <- 0.01
Reps <- 1000
Mult <- 33


bootstrapMCE <- function(data, indices){
  d <- data[indices]
  return(mean(d))
}

### Function to bootstrap an estimate for MCE by finding multiple sds

estimate_MCE <- function(p, Boot_Repeats, Multiplier, data, indices){
  
  d <- data[indices]
  
  sd_values <- vector(length = p)
  predictor1 <- vector(length = p)
  
  for (j in 1:p){
    Rj <- j*Multiplier
    predictor1[j] <- 1/sqrt(Rj)
    X_Rj <- sample(d , size = Rj, replace = FALSE, prob = NULL)
    results <- boot(data = X_Rj, R = Boot_Repeats, statistic=bootstrapMCE)
    sd_values[j] <- sd(results$t)
  }
  
  model1 <- lm(sd_values~ 0 + predictor1)
  return(model1$coefficients[[1]])
}

### Actually calculates an sample size estimate to get fixed MCE
Total_Results1 <- boot(data = Working_copy1, R = Reps, statistic = estimate_MCE, p = 3, Boot_Repeats = 100, Multiplier = Mult)
R1 <- (mean(Total_Results1$t)/Target_MCE)^2
print(R1)