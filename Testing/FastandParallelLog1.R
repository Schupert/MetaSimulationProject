### Remove previous variables
rm(list = ls())

#### Libraries
library(data.table)
library(doParallel)
library(foreach)
library(doRNG)

library(compiler)
enableJIT(3)

set.seed(123)

# Number of cores for parallel
num.Cores <- detectCores() - 1
c1 <- makeCluster(num.Cores)

#### Declare variables

# Reps = number of repetitions of experiment
Reps = 50

# k = number of studies in series
Studies = c(3,5,10,30,50,100)
#Studies = c(3,5,10,30)

# subj = number of subjects in study, likely to be distributed
Subj <- list(as.integer(c(100,100)), as.integer(c(20,100)), as.integer(c(250, 1000)), as.numeric(c(4.7, 1.2)))

# sd = study level standard deviation
True.sd = sqrt(2)

# theta = population level log(OR) - this should be considered more purely on the log scale
theta = c(log(0.25), log(0.75), log(1), log(1.5), log(4))

# tau.sq = between studies variance (can be squared due to sqrt() in normal draw), ?to be distributed
tau.sq = c(0, 0.01777778, 0.04, 3.04)

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
Tested.outcomes <- 5
Chosen.outcomes <- 1
Sd.split <- 0.8

# Size of per unit bias increase
Bias.multiple <- 0.9

### Log_Odds_Ratio changed to keep sample sizes exactly equal

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


##### Can't use doRNG on nested loops, see work around in vignette later

StartTime <- proc.time()

registerDoParallel(c1)
set.seed(123)
r <- foreach (i = Subj, .combine=rbind, .packages = c("data.table"), 
              .export = c("Studies", "Subj", "True.sd",
                          "theta", "tau.sq", "controlProp", "Severity.boundary", "Begg_a", 
                          "Begg_b", "Begg_sided", "Tested.outcomes", "Chosen.outcomes", "Sd.split",
                          "Bias.multiple", "Log_Odds_Ratio", "EvFreq")
) %:% foreach (k = theta, .combine=rbind, .packages = c("data.table"), 
               .export = c("Studies", "Subj", "True.sd",
                           "theta", "tau.sq", "controlProp", "Severity.boundary", "Begg_a", 
                           "Begg_b", "Begg_sided", "Tested.outcomes", "Chosen.outcomes", "Sd.split",
                           "Bias.multiple", "Log_Odds_Ratio", "EvFreq")
) %dopar% {
  
  ID = length(tau.sq) * Reps * sum(Studies) * length(EvFreq)
  
#   Normal.Simulation <- data.table(
#     Unique_ID = integer(length = ID),
#     # Rep_Number = integer(length = ID),
#     # Rep_Subj = list(length = ID),
#     # Rep_theta = numeric(length = ID),
#     # Rep_tau.sq = numeric(length = ID),
#     # Rep_NumStudies = numeric(length = ID),
#     # Study_ID = integer(length = ID),
#     Study_estimate = numeric(length = ID),
#     Study_sd = numeric(length = ID),
#     Study_n = integer(length = ID),
#     Study_rejectedMeans = list(length = ID),
#     Study_rejectedSDs = list(length = ID),
#     Study_Number.of.biases = integer(length = ID)
#   )
  
  LogOR.Simulation <- data.table(
    Unique_ID = integer(length = ID),
    #Rep_Number = integer(length = ID),
    #Rep_Subj = integer(length = ID),
    #Rep_ev_freq = numeric(length = ID),
    #Rep_theta = numeric(length = ID),
    #Rep_tau.sq = numeric(length = ID),
    #Rep_NumStudies = numeric(length = ID),
    #Study_ID = integer(length = ID),
    #Study_estimate = numeric(length = ID),
    #Study_sd = numeric(length = ID),
    Study_G1O1 = numeric(length = ID),
    Study_G2O1 = numeric(length = ID),
    Study_n = integer(length = ID),
#     Group1Outcome1 = numeric(length = ID),
#     Group1Outcome2 = numeric(length = ID),
#     Group2Outcome1 = numeric(length = ID),
#     Group2Outcome2 = numeric(length = ID),
#     Group1Size = integer(length = ID),
#     Group2Size = integer(length = ID),
    Study_rejectedMeans = list(length = ID),
    Study_rejectedSDs = list(length = ID),
    Study_Number.of.biases = integer(length = ID)
  )
  
  dummy.counter <- 1
  
  for (l in tau.sq){
    
    for (n in Studies){
      
      for (j in EvFreq){
      
      for (o in 1:n){
        
        for (m in 1:Reps){
          
          
#           counter <- as.integer((apply(sapply(Subj, function(vec) {i %in% vec}), 1, which.max)[1]-1) * length(theta) * length(tau.sq) * sum(Studies) * Reps +
#                                   (match(k, theta)-1) * length(tau.sq) * sum(Studies) * Reps +
#                                   (match(l, tau.sq)-1) * sum(Studies) * Reps +
#                                   (sum(Studies[0:(match(n, Studies)-1)]) + o -1) * Reps +
#                                   m
#           )
          
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
          
          x <- Log_Odds_Ratio(Study_patientnumber, k, l, controlProp, j)
          x.N <- 2*x[5]
          
          LogOR.Simulation[dummy.counter, `:=` (Unique_ID = counter,
                                                Study_G1O1 = x[1], 
                                                Study_G2O1 = x[3], 
                                                Study_n = x.N)]
          
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
ID =  length(Subj) * length(controlProp) * length(theta) * length(tau.sq) * length(EvFreq) * Reps * sum(Studies)

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
##################### Need to work out what to do about 0.5 cells

LogOR.Simulation[, Study_estimate := ifelse(LogOR.Simulation$Study_G1O1 %% 1 == 0.5, 
                                            log( (Study_G2O1 / ((Study_n/2 + 1) - Study_G2O1) ) / (Study_G1O1 / ((Study_n/2 + 1) - Study_G1O1) ) ),
                                            log( (Study_G2O1 / ((Study_n/2) - Study_G2O1) ) / (Study_G1O1 / ((Study_n/2) - Study_G1O1) ) ) )
                 ]

LogOR.Simulation[, Study_sd := ifelse(LogOR.Simulation$Study_G1O1 %% 1 == 0.5, 
                                      sqrt(1/Study_G1O1 + 1/((Study_n/2 + 1) - Study_G1O1) + 1/Study_G2O1 + 1/((Study_n/2 + 1) - Study_G2O1)),
                                      sqrt(1/Study_G1O1 + 1/(Study_n/2 - Study_G1O1) + 1/Study_G2O1 + 1/(Study_n/2 - Study_G2O1)) )
                 ]
TimeTaken <- proc.time() - StartTime

stopCluster(c1)

df.LogOR.Simulation <- as.data.frame(LogOR.Simulation)
df.LogOR.Simulation$Study_rejectedMeans <- vapply(df.LogOR.Simulation$Study_rejectedMeans, paste, collapse = ", ", character(1L))
df.LogOR.Simulation$Study_rejectedSDs <- vapply(df.LogOR.Simulation$Study_rejectedSDs, paste, collapse = ", ", character(1L))

write.csv(df.LogOR.Simulation, file = "LogORSimulation1.csv")
rm(r)