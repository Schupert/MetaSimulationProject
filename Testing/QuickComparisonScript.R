#### Remove previous variables ----
rm(list = ls())

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
theta = c(-0.76,  -0.12,  0, 0.12, 0.76)

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

sig.level <- (1 - 0.05/2)

#### Getting variables UMD -----
setwd("D:/Stats/AFP/Results/2019/UMD")

Sum_results2 <- data.table(read.csv("UMDResults_Looptest.csv"))

Sum_results2[Bias_type == "Method" & Rep_theta == theta[3] & Rep_tau.sq == tau.sq[3] & Rep_Subj == 4.2 & Rep_NumStudies == 5]

dim(Sum_results2[MSE.Moreno < MSE_FE & Bias_type == "Method"])
a <- Sum_results2[MSE.Moreno < MSE_FE & Bias_type == "Method"]
summary(as.factor(a$Rep_Subj))
summary(as.factor(a[Rep_Subj == 20]$Rep_NumStudies))


#### Remove previous variables ----
rm(list = ls())
#### LOR Declare variables ----

# Reps = number of repetitions of experiment
Reps = 10000

# k = number of studies in series
Studies = c(3,5,10,30,50,100)

# subj = number of subjects in study
Subj <- list(as.integer(c(100,100)), as.integer(c(20,100)), as.integer(c(250, 1000)), as.numeric(c(4.7, 1.2)))

# sd = study level standard deviation
True.sd = sqrt(1)

# theta = population level log(OR) - this should be considered more purely on the log scale
theta = c(log(0.25), log(0.8), log(1), log(1.25), log(4))

# tau.sq = between studies variance (can be squared due to sqrt() in normal draw)
tau.sq = c(0, 0.008, 0.04, 3.04)

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
Bias.multiple <- 0.85

#### Getting variables LOR -----

setwd("D:/Stats/AFP/Results/2019/LOR")

Sum_results2 <- data.table(read.csv("LogORResults_Looped.csv"))

Sum_results2[Bias_type == "Method" & Rep_theta == theta[3] & Rep_tau.sq == tau.sq[3] & Rep_Subj == 4.7 & Rep_NumStudies == 5 & Rep_ev_freq == 0.3]

dim(Sum_results2[MSE.Moreno < MSE_FE & Bias_type == "Step"])
a <- Sum_results2[MSE.Moreno < MSE_FE & Bias_type == "Step"]
summary(as.factor(a$Rep_ev_freq))
summary(as.factor(a[Rep_Subj == 20]$Rep_NumStudies))
