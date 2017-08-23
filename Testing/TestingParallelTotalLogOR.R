### Remove previous variables
rm(list = ls())

#### Libraries
library(data.table)
library(doParallel)
library(foreach)

# Number of cores for parallel
num.Cores <- detectCores() - 2
c1 <- makeCluster(num.Cores)

#### Set working directory - this is a problem with transferring code around computers. 
####    Current system works when opened using R project - sets wd to project wd/Results
#setwd("D:\\Stats\\AFP\\R Code")
if (length(regmatches(getwd(), gregexpr("/Results", getwd()))) == 0) 
{workingDirectory <- paste(getwd(), "/Results", sep = "")}
if (length(regmatches(workingDirectory, gregexpr("/Results", workingDirectory)))) {setwd(workingDirectory)}

#### set seed for reproduceability
set.seed(1234)

#### Declare variables

# Reps = number of repetitions of experiment
Reps = 100

# k = number of studies in series
Studies = c(2,4,6,8,10)

# subj = number of subjects in study, likely to be distributed
Subj = 100

# controlProp = proportion of total sample in control arm
controlProp = 0.5

# theta = population level log(OR) - this should be considered more purely on the log scale
theta = c(log(2), log(1))

# tau.sq = between studies variance (can be squared due to sqrt() in normal draw), ?to be distributed
tau.sq = c(1,2,3)

# Frequency of event averaged across 2 arms (before applying change due to theta) = EvFreq
EvFreq = c(0.05, 0.2, 0.5)

# ID = total number of data points required, also used as an ID number. WILL NEED UPDATING
ID =  length(Subj) * length(controlProp) * length(theta) * length(tau.sq) * length(EvFreq) * Reps * sum(Studies) 

Log_Odds_Ratio <- function(StudySize, Log_O_R, Heterogeneity, Control_Prop, mu){
  StudyLogOR <- rnorm(1, Log_O_R, sqrt(Heterogeneity))
  Group1Size <- as.integer(Control_Prop*StudySize)
  Group2Size <- as.integer(StudySize - Group1Size)
  Pic <- exp(mu - 0.5*StudyLogOR) / (1 + exp(mu - 0.5*StudyLogOR))
  Group1Out1 <- as.integer(rbinom(1, Group1Size, Pic))
  Group1Out2 <- as.integer(Group1Size - Group1Out1)
  Pit <- exp(mu + 0.5*StudyLogOR)  / (1 + exp(mu + 0.5*StudyLogOR))
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


StartTime <- proc.time()

## First parallel loop
registerDoParallel(c1)

r <- foreach (m = 1:Reps, 
              .combine=rbind, .packages = c("data.table"), 
              .export = c("Studies", "Subj", "controlProp", "theta", "tau.sq", "EvFreq", "Log_Odds_Ratio")
) %dopar% {
  
  ID.this.loop <- integer()
  
  for (i in Subj){
    
    for (j in EvFreq){
      
      for (k in theta){
        
        for (l in tau.sq){
          
          for (n in Studies){
            
            for (o in 1:n){
              
              new.this.loop <- as.integer((match(i, Subj)-1) * length(EvFreq) * length(theta) * length(tau.sq) * sum(Studies) * Reps + 
                                            (match(j, EvFreq)-1) * length(theta) * length(tau.sq) * sum(Studies) * Reps +
                                            (match(k, theta)-1) * length(tau.sq) * sum(Studies) * Reps +
                                            (match(l, tau.sq)-1) * sum(Studies) * Reps +
                                            (sum(Studies[0:(match(n, Studies)-1)]) + o -1) * Reps +
                                            m
              )
              
              ID.this.loop <- append(ID.this.loop, new.this.loop)
              
            }
          }
        }
      }
    }
  }
              
  ID <- length(ID.this.loop)
  
  LogOR.Simulation <- data.table(
    Unique_ID = c(ID.this.loop),
    Rep_Number = integer(length = ID),
    Rep_Subj = integer(length = ID),
    Rep_ev_freq = numeric(length = ID),
    Rep_theta = numeric(length = ID),
    Rep_tau.sq = numeric(length = ID),
    Rep_NumStudies = numeric(length = ID),
    Study_ID = integer(length = ID),
    Study_estimate = numeric(length = ID),
    Study_sd = numeric(length = ID),
    Study_n = integer(length = ID),
    Group1Outcome1 = numeric(length = ID),
    Group1Outcome2 = numeric(length = ID),
    Group2Outcome1 = numeric(length = ID),
    Group2Outcome2 = numeric(length = ID),
    Group1Size = integer(length = ID),
    Group2Size = integer(length = ID)
  )
  
  for (i in Subj){
    
    for (j in EvFreq){
      
      for (k in theta){
        
        for (l in tau.sq){
          
          for (n in Studies){
            
            for (o in 1:n){
              
              counter <- as.integer((match(i, Subj)-1) * length(EvFreq) * length(theta) * length(tau.sq) * sum(Studies) * Reps + 
                                      (match(j, EvFreq)-1) * length(theta) * length(tau.sq) * sum(Studies) * Reps +
                                      (match(k, theta)-1) * length(tau.sq) * sum(Studies) * Reps +
                                      (match(l, tau.sq)-1) * sum(Studies) * Reps +
                                      (sum(Studies[0:(match(n, Studies)-1)]) + o -1) * Reps +
                                      m
              )
              
              #Statement left in case of varying number of subjects later
              Study_patientnumber <- round(rlnorm(1, meanlog = 4.6, sdlog = 1))
              
              x <- Log_Odds_Ratio(Study_patientnumber, k, l, controlProp, j)
              
              LogOR.Simulation[Unique_ID == counter, `:=` (Rep_Number= m, Rep_Subj = i, Rep_ev_freq = j,
                                                           Rep_theta = k, Rep_tau.sq = l, Rep_NumStudies = n,
                                                           Study_ID = o, 
                                                           Study_estimate = log((x[3]/x[4])/(x[1]/x[2])), 
                                                           Study_sd = sqrt(1/x[1] + 1/x[2] + 1/x[3] + 1/x[4]), 
                                                           Study_n = Study_patientnumber,
                                                           Group1Outcome1 = x[1], Group1Outcome2 = x[2],
                                                           Group2Outcome1 = x[3], Group2Outcome2 = x[4],
                                                           Group1Size = x[5], Group2Size = x[6]
              )]
              
            }
          }
        }
      }
    }
  }
  LogOR.Simulation
}

write.csv(r, file = "LogORSimulation1.csv")

## Begg
registerDoParallel(c1)

r <- foreach (m = 1:Reps, 
              .combine=rbind, .packages = c("data.table"), 
              .export = c("Studies", "Subj", "controlProp", "theta", "tau.sq", "EvFreq", "Log_Odds_Ratio")
) %dopar% {
  
  # Set up strength of publication bias selection
  Begg_a <- 1.5
  Begg_b <- 4
  Begg_sided <- 1
  
  ID.this.loop <- integer()
  
  for (i in Subj){
    
    for (j in EvFreq){
      
      for (k in theta){
        
        for (l in tau.sq){
          
          for (n in Studies){
            
            for (o in 1:n){
              
              new.this.loop <- as.integer((match(i, Subj)-1) * length(EvFreq) * length(theta) * length(tau.sq) * sum(Studies) * Reps + 
                                            (match(j, EvFreq)-1) * length(theta) * length(tau.sq) * sum(Studies) * Reps +
                                            (match(k, theta)-1) * length(tau.sq) * sum(Studies) * Reps +
                                            (match(l, tau.sq)-1) * sum(Studies) * Reps +
                                            (sum(Studies[0:(match(n, Studies)-1)]) + o -1) * Reps +
                                            m
              )
              
              ID.this.loop <- append(ID.this.loop, new.this.loop)
              
            }
          }
        }
      }
    }
  }
  
  ID <- length(ID.this.loop)
  
  LogOR.Simulation <- data.table(
    Unique_ID = c(ID.this.loop),
    Rep_Number = integer(length = ID),
    Rep_Subj = integer(length = ID),
    Rep_ev_freq = numeric(length = ID),
    Rep_theta = numeric(length = ID),
    Rep_tau.sq = numeric(length = ID),
    Rep_NumStudies = numeric(length = ID),
    Study_ID = integer(length = ID),
    Study_estimate = numeric(length = ID),
    Study_sd = numeric(length = ID),
    Study_n = integer(length = ID),
    Group1Outcome1 = numeric(length = ID),
    Group1Outcome2 = numeric(length = ID),
    Group2Outcome1 = numeric(length = ID),
    Group2Outcome2 = numeric(length = ID),
    Group1Size = integer(length = ID),
    Group2Size = integer(length = ID)
  )
  
  for (i in Subj){
    
    for (j in EvFreq){
      
      for (k in theta){
        
        for (l in tau.sq){
          
          for (n in Studies){
            
            for (o in 1:n){
              
              counter <- as.integer((match(i, Subj)-1) * length(EvFreq) * length(theta) * length(tau.sq) * sum(Studies) * Reps + 
                                      (match(j, EvFreq)-1) * length(theta) * length(tau.sq) * sum(Studies) * Reps +
                                      (match(k, theta)-1) * length(tau.sq) * sum(Studies) * Reps +
                                      (match(l, tau.sq)-1) * sum(Studies) * Reps +
                                      (sum(Studies[0:(match(n, Studies)-1)]) + o -1) * Reps +
                                      m
              )
              
              #Statement left in case of varying number of subjects later
              Study_patientnumber <- round(rlnorm(1, meanlog = 4.6, sdlog = 1))
              
              ### Implement Begg publication bias
              
              repeat{
                
                x <- Log_Odds_Ratio(Study_patientnumber, k, l, controlProp, j)
                
                Study_mean <- log((x[3]/x[4])/(x[1]/x[2]))
                Study_StanDev <- sqrt(1/x[1] + 1/x[2] + 1/x[3] + 1/x[4])
                
                # Study mean changed to absolute - not described in paper
                Begg_weight <-exp(
                  -Begg_b * (
                    (Begg_sided * pnorm(- Study_mean/(Study_StanDev))) 
                    ^Begg_a ) 
                ) 
                
                if(rbinom(1,1, Begg_weight) == 1 ){break}
                
              }
              
              LogOR.Simulation[Unique_ID == counter, `:=` (Rep_Number= m, Rep_Subj = i, Rep_ev_freq = j,
                                                           Rep_theta = k, Rep_tau.sq = l, Rep_NumStudies = n,
                                                           Study_ID = o, 
                                                           Study_estimate = Study_mean, 
                                                           Study_sd = Study_StanDev, 
                                                           Study_n = Study_patientnumber,
                                                           Group1Outcome1 = x[1], Group1Outcome2 = x[2],
                                                           Group2Outcome1 = x[3], Group2Outcome2 = x[4],
                                                           Group1Size = x[5], Group2Size = x[6]
              )]
              
            }
          }
        }
      }
    }
  }
  LogOR.Simulation
}

write.csv(r, file = "LogORSimulationBegg1.csv")

## Methods

registerDoParallel(c1)

r <- foreach (m = 1:Reps, 
              .combine=rbind, .packages = c("data.table"), 
              .export = c("Studies", "Subj", "controlProp", "theta", "tau.sq", "EvFreq", "Log_Odds_Ratio")
) %dopar% {
  
  # Size of per unit bias increase
  Bias.multiple <- -0.2
  
  ID.this.loop <- integer()
  
  for (i in Subj){
    
    for (j in EvFreq){
      
      for (k in theta){
        
        for (l in tau.sq){
          
          for (n in Studies){
            
            for (o in 1:n){
              
              new.this.loop <- as.integer((match(i, Subj)-1) * length(EvFreq) * length(theta) * length(tau.sq) * sum(Studies) * Reps + 
                                            (match(j, EvFreq)-1) * length(theta) * length(tau.sq) * sum(Studies) * Reps +
                                            (match(k, theta)-1) * length(tau.sq) * sum(Studies) * Reps +
                                            (match(l, tau.sq)-1) * sum(Studies) * Reps +
                                            (sum(Studies[0:(match(n, Studies)-1)]) + o -1) * Reps +
                                            m
              )
              
              ID.this.loop <- append(ID.this.loop, new.this.loop)
              
            }
          }
        }
      }
    }
  }
  
  ID <- length(ID.this.loop)
  
  LogOR.Simulation <- data.table(
    Unique_ID = c(ID.this.loop),
    Rep_Number = integer(length = ID),
    Rep_Subj = integer(length = ID),
    Rep_ev_freq = numeric(length = ID),
    Rep_theta = numeric(length = ID),
    Rep_tau.sq = numeric(length = ID),
    Rep_NumStudies = numeric(length = ID),
    Study_ID = integer(length = ID),
    Study_estimate = numeric(length = ID),
    Study_sd = numeric(length = ID),
    Study_n = integer(length = ID),
    Group1Outcome1 = numeric(length = ID),
    Group1Outcome2 = numeric(length = ID),
    Group2Outcome1 = numeric(length = ID),
    Group2Outcome2 = numeric(length = ID),
    Group1Size = integer(length = ID),
    Group2Size = integer(length = ID),
    Study_Number.of.biases = integer(length = ID)
  )
  
  for (i in Subj){
    
    for (j in EvFreq){
      
      for (k in theta){
        
        for (l in tau.sq){
          
          for (n in Studies){
            
            for (o in 1:n){
              
              counter <- as.integer((match(i, Subj)-1) * length(EvFreq) * length(theta) * length(tau.sq) * sum(Studies) * Reps + 
                                      (match(j, EvFreq)-1) * length(theta) * length(tau.sq) * sum(Studies) * Reps +
                                      (match(k, theta)-1) * length(tau.sq) * sum(Studies) * Reps +
                                      (match(l, tau.sq)-1) * sum(Studies) * Reps +
                                      (sum(Studies[0:(match(n, Studies)-1)]) + o -1) * Reps +
                                      m
              )
              
              #Statement left in case of varying number of subjects later
              Study_patientnumber <- round(rlnorm(1, meanlog = 4.6, sdlog = 1))
              
              ## Draw from binomial how many methodological concerns study has
              #for (a in seq(50, 1000, 100)){ print(1/exp(a^0.15))}
              Number.of.biases <- rbinom(1, 3, 1/(exp(Study_patientnumber^0.15)))
              
              # Currently using alternative formulation with extra mu
              x <- Log_Odds_Ratio(Study_patientnumber, k + Number.of.biases*Bias.multiple, l, controlProp, j)
              
              LogOR.Simulation[Unique_ID == counter, `:=` (Rep_Number= m, Rep_Subj = i, Rep_ev_freq = j,
                                                           Rep_theta = k, Rep_tau.sq = l, Rep_NumStudies = n,
                                                           Study_ID = o, 
                                                           Study_estimate = log((x[3]/x[4])/(x[1]/x[2])), 
                                                           Study_sd = sqrt(1/x[1] + 1/x[2] + 1/x[3] + 1/x[4]), 
                                                           Study_n = Study_patientnumber,
                                                           Group1Outcome1 = x[1], Group1Outcome2 = x[2],
                                                           Group2Outcome1 = x[3], Group2Outcome2 = x[4],
                                                           Group1Size = x[5], Group2Size = x[6],
                                                           Study_Number.of.biases = Number.of.biases
              )]
              
            }
          }
        }
      }
    }
  }
  LogOR.Simulation
}

write.csv(r, file = "LogORSimulationMethods1.csv")

## Outcome bais
registerDoParallel(c1)

r <- foreach (m = 1:Reps, 
              .combine=rbind, .packages = c("data.table"), 
              .export = c("Studies", "Subj", "controlProp", "theta", "tau.sq", "EvFreq", "Log_Odds_Ratio")
) %dopar% {
  
  ID.this.loop <- integer()
  
  # Set up within study reporting bias - this is now one sided
  Tested.outcomes <- 10
  Chosen.outcomes <- 1
  
  for (i in Subj){
    
    for (j in EvFreq){
      
      for (k in theta){
        
        for (l in tau.sq){
          
          for (n in Studies){
            
            for (o in 1:n){
              
              new.this.loop <- as.integer((match(i, Subj)-1) * length(EvFreq) * length(theta) * length(tau.sq) * sum(Studies) * Reps + 
                                            (match(j, EvFreq)-1) * length(theta) * length(tau.sq) * sum(Studies) * Reps +
                                            (match(k, theta)-1) * length(tau.sq) * sum(Studies) * Reps +
                                            (match(l, tau.sq)-1) * sum(Studies) * Reps +
                                            (sum(Studies[0:(match(n, Studies)-1)]) + o -1) * Reps +
                                            m
              )
              
              ID.this.loop <- append(ID.this.loop, new.this.loop)
              
            }
          }
        }
      }
    }
  }
  
  ID <- length(ID.this.loop)
  
  LogOR.Simulation <- data.table(
    Unique_ID = c(ID.this.loop),
    Rep_Number = integer(length = ID),
    Rep_Subj = integer(length = ID),
    Rep_ev_freq = numeric(length = ID),
    Rep_theta = numeric(length = ID),
    Rep_tau.sq = numeric(length = ID),
    Rep_NumStudies = numeric(length = ID),
    Study_ID = integer(length = ID),
    Study_estimate = numeric(length = ID),
    Study_sd = numeric(length = ID),
    Study_n = integer(length = ID),
    Group1Outcome1 = numeric(length = ID),
    Group1Outcome2 = numeric(length = ID),
    Group2Outcome1 = numeric(length = ID),
    Group2Outcome2 = numeric(length = ID),
    Group1Size = integer(length = ID),
    Group2Size = integer(length = ID),
    Study_rejectedMeans = list(length = ID),
    Study_rejectedSDs = list(length = ID)
  )
  
  for (i in Subj){
    
    for (j in EvFreq){
      
      for (k in theta){
        
        for (l in tau.sq){
          
          for (n in Studies){
            
            for (o in 1:n){
              
              counter <- as.integer((match(i, Subj)-1) * length(EvFreq) * length(theta) * length(tau.sq) * sum(Studies) * Reps + 
                                      (match(j, EvFreq)-1) * length(theta) * length(tau.sq) * sum(Studies) * Reps +
                                      (match(k, theta)-1) * length(tau.sq) * sum(Studies) * Reps +
                                      (match(l, tau.sq)-1) * sum(Studies) * Reps +
                                      (sum(Studies[0:(match(n, Studies)-1)]) + o -1) * Reps +
                                      m
              )
              
              #Statement left in case of varying number of subjects later
              Study_patientnumber <- round(rlnorm(1, meanlog = 4.6, sdlog = 1))
              
              # In the log(OR) condition tau2 variance is split between study level and person level
              Study_mu <- rnorm(1, mean = k, sd = sqrt(0.5 * l))
              #Person_values <- rnorm(Study_patientnumber, mean = Study_mu, sd = sqrt(0.5*l))
              
              # Sample multiple outcome measures from same set of patients, using Person_values as mean
              Study_values <- replicate(Tested.outcomes, Log_Odds_Ratio(Study_patientnumber, Study_mu, 0.5 * l, controlProp, j) )
              
              Study_mean <- numeric(length = Tested.outcomes)
              Study_StanDev <- numeric(length = Tested.outcomes)
              Begg_p <- numeric(length = Tested.outcomes)
              
              for (z in 1:Tested.outcomes) {
                Study_mean[z] <- log((Study_values[3,z]/Study_values[4,z])/(Study_values[1,z]/Study_values[2,z]))
                Study_StanDev[z] <- sqrt(1/Study_values[1,z] + 1/Study_values[2,z] + 1/Study_values[3,z] + 1/Study_values[4,z])
                Begg_p[z] <- pnorm(-Study_mean[z]/(Study_StanDev[z]))
              }
              
              lv <- which.min(Begg_p)
              
              LogOR.Simulation[Unique_ID == counter, `:=` (Rep_Number= m, Rep_Subj = i, Rep_ev_freq = j,
                                                           Rep_theta = k, Rep_tau.sq = l, Rep_NumStudies = n,
                                                           Study_ID = o, 
                                                           Study_estimate = Study_mean[lv], 
                                                           Study_sd = Study_StanDev[lv],
                                                           Study_n = Study_patientnumber,
                                                           Group1Outcome1 = Study_values[1,lv], Group1Outcome2 = Study_values[2,lv],
                                                           Group2Outcome1 = Study_values[3,lv], Group2Outcome2 = Study_values[4,lv],
                                                           Group1Size = Study_values[5,lv], Group2Size = Study_values[6,lv],
                                                           Study_rejectedMeans = list(Study_mean[-lv]),
                                                           Study_rejectedSDs = list(Study_StanDev[-lv])
                                                           )
                               ]
              
            }
          }
        }
      }
    }
  }
  LogOR.Simulation
}

#write.csv(r, file = "LogORSimulationOutcome1.csv")

df.LogOR.Simulation <- as.data.frame(r)
df.LogOR.Simulation$Study_rejectedMeans <- vapply(df.LogOR.Simulation$Study_rejectedMeans, paste, collapse = ", ", character(1L))
df.LogOR.Simulation$Study_rejectedSDs <- vapply(df.LogOR.Simulation$Study_rejectedSDs, paste, collapse = ", ", character(1L))
write.csv(df.LogOR.Simulation, file = "NormalSimulationOutcomeBias1.csv")

## Step
registerDoParallel(c1)

r <- foreach (m = 1:Reps, 
              .combine=rbind, .packages = c("data.table"), 
              .export = c("Studies", "Subj", "controlProp", "theta", "tau.sq", "EvFreq", "Log_Odds_Ratio")
) %dopar% {
  
  # Boundary of step function on p value, causing severity of publication bias
  Severity.boundary <- c(0.05, 0.2)
  
  ID.this.loop <- integer()
  
  for (i in Subj){
    
    for (j in EvFreq){
      
      for (k in theta){
        
        for (l in tau.sq){
          
          for (n in Studies){
            
            for (o in 1:n){
              
              new.this.loop <- as.integer((match(i, Subj)-1) * length(EvFreq) * length(theta) * length(tau.sq) * sum(Studies) * Reps + 
                                            (match(j, EvFreq)-1) * length(theta) * length(tau.sq) * sum(Studies) * Reps +
                                            (match(k, theta)-1) * length(tau.sq) * sum(Studies) * Reps +
                                            (match(l, tau.sq)-1) * sum(Studies) * Reps +
                                            (sum(Studies[0:(match(n, Studies)-1)]) + o -1) * Reps +
                                            m
              )
              
              ID.this.loop <- append(ID.this.loop, new.this.loop)
              
            }
          }
        }
      }
    }
  }
  
  ID <- length(ID.this.loop)
  
  LogOR.Simulation <- data.table(
    Unique_ID = c(ID.this.loop),
    Rep_Number = integer(length = ID),
    Rep_Subj = integer(length = ID),
    Rep_ev_freq = numeric(length = ID),
    Rep_theta = numeric(length = ID),
    Rep_tau.sq = numeric(length = ID),
    Rep_NumStudies = numeric(length = ID),
    Study_ID = integer(length = ID),
    Study_estimate = numeric(length = ID),
    Study_sd = numeric(length = ID),
    Study_n = integer(length = ID),
    Group1Outcome1 = numeric(length = ID),
    Group1Outcome2 = numeric(length = ID),
    Group2Outcome1 = numeric(length = ID),
    Group2Outcome2 = numeric(length = ID),
    Group1Size = integer(length = ID),
    Group2Size = integer(length = ID)
  )
  
  for (i in Subj){
    
    for (j in EvFreq){
      
      for (k in theta){
        
        for (l in tau.sq){
          
          for (n in Studies){
            
            for (o in 1:n){
              
              counter <- as.integer((match(i, Subj)-1) * length(EvFreq) * length(theta) * length(tau.sq) * sum(Studies) * Reps + 
                                      (match(j, EvFreq)-1) * length(theta) * length(tau.sq) * sum(Studies) * Reps +
                                      (match(k, theta)-1) * length(tau.sq) * sum(Studies) * Reps +
                                      (match(l, tau.sq)-1) * sum(Studies) * Reps +
                                      (sum(Studies[0:(match(n, Studies)-1)]) + o -1) * Reps +
                                      m
              )
              
              #Statement left in case of varying number of subjects later
              Study_patientnumber <- round(rlnorm(1, meanlog = 4.6, sdlog = 1))
              
              ### Publication bias by step function - currently one sided as per Hedges
              
              repeat{
                
                x <- Log_Odds_Ratio(Study_patientnumber, k, l, controlProp, j)
                
                Study_mean <- log((x[3]/x[4])/(x[1]/x[2]))
                Study_StanDev <- sqrt(1/x[1] + 1/x[2] + 1/x[3] + 1/x[4])
                
                Begg_p <- pnorm(- Study_mean/(Study_StanDev))
                
                Step_weight <- ifelse(Begg_p < Severity.boundary[1], 1, ifelse(Begg_p < Severity.boundary[2], 0.75, 0.25))
                
                if(rbinom(1,1, Step_weight) == 1 ){break}
                
              }
              
              LogOR.Simulation[Unique_ID == counter, `:=` (Rep_Number= m, Rep_Subj = i, Rep_ev_freq = j,
                                                           Rep_theta = k, Rep_tau.sq = l, Rep_NumStudies = n,
                                                           Study_ID = o, 
                                                           Study_estimate = Study_mean, 
                                                           Study_sd = Study_StanDev, 
                                                           Study_n = Study_patientnumber,
                                                           Group1Outcome1 = x[1], Group1Outcome2 = x[2],
                                                           Group2Outcome1 = x[3], Group2Outcome2 = x[4],
                                                           Group1Size = x[5], Group2Size = x[6]
              )]
              
            }
          }
        }
      }
    }
  }
  LogOR.Simulation
}

write.csv(r, file = "LogORSimulationStep1.csv")

TimeTaken <- proc.time() - StartTime

stopCluster(c1)