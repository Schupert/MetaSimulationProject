### Remove previous variables
rm(list = ls())

#### Libraries
library(data.table)
library(doParallel)
library(foreach)
library(doRNG)
library(metafor)

library(compiler)
enableJIT(3)

## Set working directory

#LogOR.Simulation <- read.csv("LogORSimulation1.csv")
#LogOR.Simulation <- data.table(LogOR.Simulation)
setkey(LogOR.Simulation, "Rep_Number", "Rep_Subj", "Rep_theta", "Rep_tau.sq", "Rep_NumStudies", "Rep_ev_freq")



# Number of cores for parallel
num.Cores <- detectCores() - 1
c1 <- makeCluster(num.Cores)

#### Declare variables

# Reps = number of repetitions of experiment
Reps = 10

# k = number of studies in series
Studies = c(3,5,10,30,50,100)
#Studies = c(3,5,10,30)

# subj = number of subjects in study, likely to be distributed
#Subj <- list(as.integer(c(100,100)), as.integer(c(30,40)), as.integer(c(250, 1000)), as.numeric(c(4.7, 1.2)))
Subj <- c(100, 30, 250, 4.7)

# sd = study level standard deviation
True.sd = 2

# theta = population level log(OR) - this should be considered more purely on the log scale
theta = c(log(0.5), log(1), log(1.5), log(3))

# tau.sq = between studies variance (can be squared due to sqrt() in normal draw), ?to be distributed
tau.sq = c(1,2,3)

# Frequency of event averaged across 2 arms (before applying change due to theta) = EvFreq
EvFreq = c(0.1, 0.3, 0.5)

# controlProp = proportion of total sample in control arm
controlProp = 0.5

## Boundary of step function on p value, causing severity of publication bias
Severity.boundary <- c(0.05, 0.2)

# Set up strength of publication bias selection IF STILL USING
Begg_a <- 1.5
Begg_b <- 4
Begg_sided <- 1

# Set up within study reporting bias - this is now one sided
Tested.outcomes <- 10
Chosen.outcomes <- 1
Sd.split <- 0.5

# Size of per unit bias increase
Bias.multiple <- 1/0.9


StartTime <- proc.time()

registerDoParallel(c1)

r <- foreach (i = Subj, .combine=rbind, .packages = c("data.table", "metafor"), 
              .export = c("Studies", "Subj", "True.sd",
                          "theta", "tau.sq", "EvFreq", "LogOR.Simulation")
) %:% foreach (k = theta, .combine=rbind, .packages = c("data.table", "metafor"), 
               .export = c("Studies", "Subj", "True.sd",
                           "theta", "tau.sq", "EvFreq", "LogOR.Simulation")
) %dopar% {
  # ID different for analysis
  ID = length(tau.sq) * Reps * length(Studies) * length(EvFreq)
  
  LogOR.Sim.Results <- data.table(
    Unique_ID = integer(length = ID),
    Rep_Number = integer(length = ID),
    Rep_Subj = numeric(length = ID),
    Rep_theta = numeric(length = ID),
    Rep_tau.sq = numeric(length = ID),
    Rep_NumStudies = numeric(length = ID),
    Rep_ev_freq = numeric(length = ID),
    FE_Estimate = numeric(length = ID),
    FE_Est_Low_CI = numeric(length = ID),
    FE_Est_Up_CI = numeric(length = ID),
    REML_Estimate = numeric(length = ID),
    REML_Est_Low_CI = numeric(length = ID),
    REML_Est_Up_CI = numeric(length = ID),
    REML_tau2 = numeric(length = ID),
    REML_I2 = numeric(length = ID),
    DL_Estimate = numeric(length = ID),
    DL_Est_Low_CI = numeric(length = ID),
    DL_Est_Up_CI = numeric(length = ID),
    DL_tau2 = numeric(length = ID),
    DL_I2 = numeric(length = ID)
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
          
          
          ma.fe <- rma.uni(temp.data$Study_estimate, temp.data$Study_sd , method = "FE")
          ma.reml <- rma.uni(temp.data$Study_estimate, temp.data$Study_sd  , method = "REML")
          ma.DL <- rma.uni(temp.data$Study_estimate, temp.data$Study_sd  , method = "DL")
          
          #         LogOR.Sim.Results[dummy.counter, `:=` (Rep_Number= m, Rep_Subj = i, 
          #                                                 Rep_theta = k, Rep_tau.sq = l, Rep_NumStudies = n,
          #                                                 FE_Estimate = ma.fe[[1]],
          #                                                 FE_Est_Low_CI = ma.fe[[5]],
          #                                                 FE_Est_Up_CI = ma.fe[[6]],
          #                                                 REML_Estimate = ma.reml[[1]],
          #                                                 REML_Est_Low_CI = ma.reml[[5]],
          #                                                 REML_Est_Up_CI = ma.reml[[6]],
          #                                                 REML_tau2 = ma.reml[[8]],
          #                                                 REML_I2 = ma.reml[[22]],
          #                                                 DL_Estimate = ma.DL[[1]],
          #                                                 DL_Est_Low_CI = ma.DL[[5]],
          #                                                 DL_Est_Up_CI = ma.DL[[6]],
          #                                                 DL_tau2 = ma.DL[[8]],
          #                                                 DL_I2 = ma.DL[[22]])]
          
          LogOR.Sim.Results[dummy.counter, `:=` (Unique_ID = counter,
                                                 FE_Estimate = ma.fe[[1]],
                                                 FE_Est_Low_CI = ma.fe[[5]],
                                                 FE_Est_Up_CI = ma.fe[[6]],
                                                 REML_Estimate = ma.reml[[1]],
                                                 REML_Est_Low_CI = ma.reml[[5]],
                                                 REML_Est_Up_CI = ma.reml[[6]],
                                                 REML_tau2 = ma.reml[[8]],
                                                 REML_I2 = ma.reml[[22]],
                                                 DL_Estimate = ma.DL[[1]],
                                                 DL_Est_Low_CI = ma.DL[[5]],
                                                 DL_Est_Up_CI = ma.DL[[6]],
                                                 DL_tau2 = ma.DL[[8]],
                                                 DL_I2 = ma.DL[[22]])]
          
          
          dummy.counter <- dummy.counter + 1
          
        }
      }
    }
  }
  LogOR.Sim.Results
}

LogOR.Sim.Results <- r[order(Unique_ID)]

##### Need to re append values - specific to analysis
ID =  length(Subj) * length(controlProp) * length(theta) * length(tau.sq) * length(EvFreq) * Reps * length(Studies)

LogOR.Sim.Results$Rep_Number =  rep(1:Reps, times = ID/Reps)
LogOR.Sim.Results$Rep_NumStudies = rep(Studies, times = ID/(Reps*length(Studies)))
LogOR.Sim.Results$Rep_ev_freq = rep(rep(EvFreq, each = Reps * length(Studies)), times = ID/(Reps*length(Studies)*length(EvFreq)))
LogOR.Sim.Results$Rep_tau.sq = rep(rep(tau.sq, each = Reps * length(Studies)*length(EvFreq)), times = ID/(Reps*length(Studies)*length(tau.sq)*length(EvFreq)))
LogOR.Sim.Results$Rep_theta = rep( rep(theta, each = Reps * length(Studies) * length(tau.sq)*length(EvFreq)), times = length(Subj))

# ### Create keyable vector for Subj
# Subj2 <- c(100, 30, 250, 4.7)
LogOR.Sim.Results$Rep_Subj = rep(Subj, each = ID / length(Subj))

TimeTaken <- proc.time() - StartTime

write.csv(LogOR.Sim.Results, file = "LogORSimAnalysis1.csv")


rm(r)

