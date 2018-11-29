### Set seed
set.seed(1)

### Import required libraries
require(boot)
require(data.table)

#### UMD Declare variables ----

# Reps = number of repetitions of experiment
Reps = 10000

# k = number of studies in series
Studies = c(3,5,10,30,50,100)

# subj = number of subjects in study, likely to be distributed
Subj <- list(as.integer(c(60,60)), as.integer(c(20,100)), as.integer(c(250, 1000)), as.numeric(c(4.2, 1.1)))

# sd = study level standard deviation
True.sd = sqrt(1)

# theta = population level mean - need good sense of range for SMD
theta = c( -0.76,  -0.12,  0, 0.12, 0.76)

# tau.sq = between studies variance (can be squared due to sqrt() in normal draw), ?to be distributed
tau.sq = c(0, 0.005, 0.022, 1.676)

# controlProp = proportion of total sample in control arm
controlProp = 0.5

## Boundary of step function on p value, causing severity of publication bias
Severity.boundary <- c(0.05, 0.2)

# Set up strength of publication bias selection IF STILL USING
Begg_a <- 0.5
Begg_b <- 3
Begg_c <- -0.3
Begg_sided <- 1

# Set up within study reporting bias - this is now one sided
Tested.outcomes <- 5
Sd.split <- 0.6

# Size of per unit bias increase
Bias.multiple <- c(0, log(0.85)/(-1.81) * 2, log(0.7225)/(-1.81) * 2)

#####

system.time(Normal.Sim.Results <- readRDS(file = "NSTotalV2RDS"))

Normal.Sim.Results <- data.table(Normal.Sim.Results)
r <- Normal.Sim.Results[Rep_Subj == 4.2 & Rep_theta == theta[3] & Rep_tau.sq == tau.sq[1] & Rep_NumStudies == Studies[1]]
Working_copy1 <- r[1:50,]$FE_Estimate

### Get data and set keyframes
#r <- read.csv("NormalSimulation1Analysis.csv")
#r <- data.frame(r)
#setkey(r, "Rep_Number", "Rep_Subj", "Rep_theta", "Rep_tau.sq", "Rep_NumStudies")

###### Select specific set of results
#Working_copy1 <- r[Rep_Subj == 1 & Rep_theta == 1 & Rep_tau.sq ==3 & Rep_NumStudies == 3]$FE_Estimate
#Working_copy1 <- LogOR.Sim.Results[Rep_Subj == 250 & Rep_theta == log(1) & Rep_tau.sq == 1 & Rep_NumStudies == 100 & Rep_ev_freq == 0.5]$REML_Estimate

## Can increase multiplier when have more reps to estimate with up to available data/3
Target_MCE <- 0.001
Reps <- 1000
Mult <- 16


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