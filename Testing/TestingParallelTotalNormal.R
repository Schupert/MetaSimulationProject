### Remove previous variables
rm(list = ls())

#### Libraries
library(data.table)
library(doParallel)
library(foreach)

# Number of cores for parallel
num.Cores <- detectCores() - 2
c1 <- makeCluster(num.Cores)

#### Declare variables

# Reps = number of repetitions of experiment
Reps = 1000

# k = number of studies in series
Studies = c(2,4,6,8,10)

# subj = number of subjects in study, likely to be distributed
Subj = c(3.5, 4)

# sd = study level standard deviation
True.sd = c(2,3)

# theta = population level mean
theta = c(0,5)

# tau.sq = between studies variance (can be squared due to sqrt() in normal draw), ?to be distributed
tau.sq = c(1,2,3)


#### Set working directory - this is a problem with transferring code around computers. 
####    Current system works when opened using R project - sets wd to project wd/Results
#setwd("D:\\Stats\\AFP\\R Code")
if (length(regmatches(getwd(), gregexpr("/Results", getwd()))) == 0) 
{workingDirectory <- paste(getwd(), "/Results", sep = "")}
if (length(regmatches(workingDirectory, gregexpr("/Results", workingDirectory)))) {setwd(workingDirectory)}

StartTime <- proc.time()

## First parallel loop
registerDoParallel(c1)

r <- foreach (m = 1:Reps, .combine=rbind, .packages = c("data.table"), .export = c("Studies", "Subj", "True.sd",
                                                                                   "theta", "tau.sq")) %dopar% {
  
  
  ID.this.loop <- integer()
  for (i in Subj){
    
    for (j in True.sd){
      
      for (k in theta){
        
        for (l in tau.sq){
          
          for (n in Studies){
            
            for (o in 1:n){
              new.this.loop <- as.integer((match(i, Subj)-1) * length(True.sd) * length(theta) * length(tau.sq) * sum(Studies) * Reps + 
                                            (match(j, True.sd)-1) * length(theta) * length(tau.sq) * sum(Studies) * Reps +
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
  
  Normal.Simulation <- data.table(
    Unique_ID = c(ID.this.loop),
    Rep_Number = integer(length = ID),
    Rep_Subj = numeric(length = ID),
    Rep_sd = numeric(length = ID),
    Rep_theta = numeric(length = ID),
    Rep_tau.sq = numeric(length = ID),
    Rep_NumStudies = numeric(length = ID),
    Study_ID = integer(length = ID),
    Study_estimate = numeric(length = ID),
    Study_sd = numeric(length = ID),
    Study_n = integer(length = ID)
  )
  
  
  for (i in Subj){
    
    for (j in True.sd){
      
      for (k in theta){
        
        for (l in tau.sq){
          
          
          
          for (n in Studies){
            
            for (o in 1:n){
              
              counter <- as.integer((match(i, Subj)-1) * length(True.sd) * length(theta) * length(tau.sq) * sum(Studies) * Reps + 
                                      (match(j, True.sd)-1) * length(theta) * length(tau.sq) * sum(Studies) * Reps +
                                      (match(k, theta)-1) * length(tau.sq) * sum(Studies) * Reps +
                                      (match(l, tau.sq)-1) * sum(Studies) * Reps +
                                      (sum(Studies[0:(match(n, Studies)-1)]) + o -1) * Reps +
                                      m
              )
            
              #Statement left in case of varying number of subjects later
              Study_patientnumber <- round(rlnorm(1, meanlog = 4.2, sdlog = 1.1)+2)
              
              Study_mu <- rnorm(1, mean = k, sd = sqrt(l))
              Study_values <- rnorm(Study_patientnumber, mean = Study_mu, sd = j)
              Study_mean <- mean(Study_values)
              Study_StanDev <- sd(Study_values)
              
              Normal.Simulation[Unique_ID == counter, `:=` (Rep_Number= m, Rep_Subj = i, Rep_sd = j,
                                                            Rep_theta = k, Rep_tau.sq = l, Rep_NumStudies = n,
                                                            Study_ID = o, 
                                                            Study_estimate = Study_mean, Study_sd = Study_StanDev, 
                                                            Study_n = Study_patientnumber)]
              
              
              
            }
          }
          
        }
      }
    }
  }
  
  Normal.Simulation
}

write.csv(r, file = "NormalSimulation1.csv")

## Begg parallel loop

registerDoParallel(c1)

r <- foreach (m = 1:Reps, .combine=rbind, .packages = c("data.table"), .export = c("Studies", "Subj", "True.sd",
                                                                                   "theta", "tau.sq")) %dopar% {
  
  
  # Set up strength of publication bias selection
  Begg_a <- 1.5
  Begg_b <- 4
  Begg_sided <- 1
  
  
  ID.this.loop <- integer()
  for (i in Subj){
    
    for (j in True.sd){
      
      for (k in theta){
        
        for (l in tau.sq){
          
          for (n in Studies){
            
            for (o in 1:n){
              new.this.loop <- as.integer((match(i, Subj)-1) * length(True.sd) * length(theta) * length(tau.sq) * sum(Studies) * Reps + 
                                            (match(j, True.sd)-1) * length(theta) * length(tau.sq) * sum(Studies) * Reps +
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
  
  Normal.Simulation <- data.table(
    Unique_ID = c(ID.this.loop),
    Rep_Number = integer(length = ID),
    Rep_Subj = numeric(length = ID),
    Rep_sd = numeric(length = ID),
    Rep_theta = numeric(length = ID),
    Rep_tau.sq = numeric(length = ID),
    Rep_NumStudies = numeric(length = ID),
    Study_ID = integer(length = ID),
    Study_estimate = numeric(length = ID),
    Study_sd = numeric(length = ID),
    Study_n = integer(length = ID)
  )
  
  
  for (i in Subj){
    
    for (j in True.sd){
      
      for (k in theta){
        
        for (l in tau.sq){
          
          
          
          for (n in Studies){
            
            for (o in 1:n){
              
              counter <- as.integer((match(i, Subj)-1) * length(True.sd) * length(theta) * length(tau.sq) * sum(Studies) * Reps + 
                                      (match(j, True.sd)-1) * length(theta) * length(tau.sq) * sum(Studies) * Reps +
                                      (match(k, theta)-1) * length(tau.sq) * sum(Studies) * Reps +
                                      (match(l, tau.sq)-1) * sum(Studies) * Reps +
                                      (sum(Studies[0:(match(n, Studies)-1)]) + o -1) * Reps +
                                      m
              )
              
              #Statement left in case of varying number of subjects later
              Study_patientnumber <- round(rlnorm(1, meanlog = 4.2, sdlog = 1.1)+2)
              
              ### Implement Begg and Mazumdar publication bias
              repeat{
                
                Study_mu <- rnorm(1, mean = k, sd = sqrt(l))
                Study_values <- rnorm(Study_patientnumber, mean = Study_mu, sd = j)
                Study_mean <- mean(Study_values)
                Study_StanDev <- sd(Study_values)
                
                # Study mean changed to absolute - not described in paper
                Begg_weight <-exp(
                  -Begg_b * (
                    (Begg_sided * pnorm(- Study_mean/(Study_StanDev))) 
                    ^Begg_a ) 
                ) 
                
                if(rbinom(1,1, Begg_weight) == 1 ){break}
                
              }
              
              Normal.Simulation[Unique_ID == counter, `:=` (Rep_Number= m, Rep_Subj = i, Rep_sd = j,
                                                            Rep_theta = k, Rep_tau.sq = l, Rep_NumStudies = n,
                                                            Study_ID = o, 
                                                            Study_estimate = Study_mean, Study_sd = Study_StanDev, 
                                                            Study_n = Study_patientnumber)]
              
              
              
            }
          }
          
        }
      }
    }
  }
  
  Normal.Simulation
}

write.csv(r, file = "NormalSimulationBegg1.csv")

## Third parallel loop - Step publication bias

registerDoParallel(c1)

r <- foreach (m = 1:Reps, .combine=rbind, .packages = c("data.table"), .export = c("Studies", "Subj", "True.sd",
                                                                                   "theta", "tau.sq")) %dopar% {
  
  ## Boundary of step function on p value, causing severity of publication bias
  #Severity.boundary <- c(0.05, 0.2)
  
  
  ID.this.loop <- integer()
  for (i in Subj){
    
    for (j in True.sd){
      
      for (k in theta){
        
        for (l in tau.sq){
          
          for (n in Studies){
            
            for (o in 1:n){
              new.this.loop <- as.integer((match(i, Subj)-1) * length(True.sd) * length(theta) * length(tau.sq) * sum(Studies) * Reps + 
                                            (match(j, True.sd)-1) * length(theta) * length(tau.sq) * sum(Studies) * Reps +
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
  
  Normal.Simulation <- data.table(
    Unique_ID = c(ID.this.loop),
    Rep_Number = integer(length = ID),
    Rep_Subj = numeric(length = ID),
    Rep_sd = numeric(length = ID),
    Rep_theta = numeric(length = ID),
    Rep_tau.sq = numeric(length = ID),
    Rep_NumStudies = numeric(length = ID),
    Study_ID = integer(length = ID),
    Study_estimate = numeric(length = ID),
    Study_sd = numeric(length = ID),
    Study_n = integer(length = ID)
  )
  
  
  for (i in Subj){
    
    for (j in True.sd){
      
      for (k in theta){
        
        for (l in tau.sq){
          
          
          
          for (n in Studies){
            
            for (o in 1:n){
              
              counter <- as.integer((match(i, Subj)-1) * length(True.sd) * length(theta) * length(tau.sq) * sum(Studies) * Reps + 
                                      (match(j, True.sd)-1) * length(theta) * length(tau.sq) * sum(Studies) * Reps +
                                      (match(k, theta)-1) * length(tau.sq) * sum(Studies) * Reps +
                                      (match(l, tau.sq)-1) * sum(Studies) * Reps +
                                      (sum(Studies[0:(match(n, Studies)-1)]) + o -1) * Reps +
                                      m
              )
              
              #Statement left in case of varying number of subjects later
              Study_patientnumber <- round(rlnorm(1, meanlog = 4.2, sdlog = 1.1)+2)
              
              ### Publication bias by step function - currently one sided as per Hedges
              
              repeat{
                
                Study_mu <- rnorm(1, mean = k, sd = sqrt(l))
                Study_values <- rnorm(Study_patientnumber, mean = Study_mu, sd = j)
                Study_mean <- mean(Study_values)
                Study_StanDev <- sd(Study_values)
                
                Begg_p <- pnorm(- Study_mean/(Study_StanDev))
                
                Step_weight <- ifelse(Begg_p < Severity.boundary[1], 1, ifelse(Begg_p < Severity.boundary[2], 0.75, 0.25))
                
                if(rbinom(1,1, Step_weight) == 1 ){break}
                
              }
              
              Normal.Simulation[Unique_ID == counter, `:=` (Rep_Number= m, Rep_Subj = i, Rep_sd = j,
                                                            Rep_theta = k, Rep_tau.sq = l, Rep_NumStudies = n,
                                                            Study_ID = o, 
                                                            Study_estimate = Study_mean, Study_sd = Study_StanDev, 
                                                            Study_n = Study_patientnumber)]
              
              
            }
          }
          
        }
      }
    }
  }
  
  Normal.Simulation
}

write.csv(r, file = "NormalSimulationStepBias1.csv")

############# Outcome bias

registerDoParallel(c1)

r <- foreach (m = 1:Reps, .combine=rbind, .packages = c("data.table"), .export = c("Studies", "Subj", "True.sd",
                                                                                   "theta", "tau.sq")) %dopar% {
  
  
  # Set up within study reporting bias - this is now one sided
  Tested.outcomes <- 10
  Chosen.outcomes <- 1
  
  
  ID.this.loop <- integer()
  for (i in Subj){
    
    for (j in True.sd){
      
      for (k in theta){
        
        for (l in tau.sq){
          
          for (n in Studies){
            
            for (o in 1:n){
              new.this.loop <- as.integer((match(i, Subj)-1) * length(True.sd) * length(theta) * length(tau.sq) * sum(Studies) * Reps + 
                                            (match(j, True.sd)-1) * length(theta) * length(tau.sq) * sum(Studies) * Reps +
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
  
  Normal.Simulation <- data.table(
    Unique_ID = c(ID.this.loop),
    Rep_Number = integer(length = ID),
    Rep_Subj = numeric(length = ID),
    Rep_sd = numeric(length = ID),
    Rep_theta = numeric(length = ID),
    Rep_tau.sq = numeric(length = ID),
    Rep_NumStudies = numeric(length = ID),
    Study_ID = integer(length = ID),
    Study_estimate = numeric(length = ID),
    Study_sd = numeric(length = ID),
    Study_n = integer(length = ID),
    Study_rejectedMeans = list(length = ID),
    Study_rejectedSDs = list(length = ID)
  )
  
  
  for (i in Subj){
    
    for (j in True.sd){
      
      for (k in theta){
        
        for (l in tau.sq){
          
          
          
          for (n in Studies){
            
            for (o in 1:n){
              
              counter <- as.integer((match(i, Subj)-1) * length(True.sd) * length(theta) * length(tau.sq) * sum(Studies) * Reps + 
                                      (match(j, True.sd)-1) * length(theta) * length(tau.sq) * sum(Studies) * Reps +
                                      (match(k, theta)-1) * length(tau.sq) * sum(Studies) * Reps +
                                      (match(l, tau.sq)-1) * sum(Studies) * Reps +
                                      (sum(Studies[0:(match(n, Studies)-1)]) + o -1) * Reps +
                                      m
              )
              
              #Statement left in case of varying number of subjects later
              Study_patientnumber <- round(rlnorm(1, meanlog = 4.2, sdlog = 1.1)+2)
              
              ### Implement Within study multiple outcomes bias - split variance by simulating values
              # for each individual which represent the between person variability, then sample from these
              # with sd to simulate testing variability
              
              
              Study_mu <- rnorm(1, mean = k, sd = sqrt(l))
              
              # Person level 'true' values with half the total sd
              Person_values <- rnorm(Study_patientnumber, mean = Study_mu, sd = (1/sqrt(2))*j)
              
              # Sample multiple outcome measures from same set of patients, using Person_values as mean
              Study_values <- as.vector(replicate(Tested.outcomes, rnorm(Study_patientnumber, mean = Person_values, sd = (1/sqrt(2))*j)))
              
              Study_mean <- numeric(length = Tested.outcomes)
              Study_StanDev <- numeric(length = Tested.outcomes)
              Begg_p <- numeric(length = Tested.outcomes)
              
              for (z in 1:Tested.outcomes) {
                Study_mean[z] <- mean( Study_values[((z-1)*Study_patientnumber + 1): (z*Study_patientnumber)] )
                Study_StanDev[z] <- sd( Study_values[((z-1)*Study_patientnumber + 1): (z*Study_patientnumber)] )
                Begg_p[z] <- pnorm(-Study_mean[z]/(Study_StanDev[z]))
              }
              
              lv <- which.min(Begg_p)
              
              Normal.Simulation[Unique_ID == counter, `:=` (Rep_Number= m, Rep_Subj = i, Rep_sd = j,
                                                            Rep_theta = k, Rep_tau.sq = l, Rep_NumStudies = n,
                                                            Study_ID = o, 
                                                            Study_estimate = Study_mean[lv], 
                                                            Study_sd = Study_StanDev[lv], 
                                                            Study_n = Study_patientnumber,
                                                            Study_rejectedMeans = list(Study_mean[-lv]),
                                                            Study_rejectedSDs = list(Study_StanDev[-lv])
              )]
            }
          }
          
        }
      }
    }
  }
  
  Normal.Simulation
}

df.Normal.Simulation <- as.data.frame(r)
df.Normal.Simulation$Study_rejectedMeans <- vapply(df.Normal.Simulation$Study_rejectedMeans, paste, collapse = ", ", character(1L))
df.Normal.Simulation$Study_rejectedSDs <- vapply(df.Normal.Simulation$Study_rejectedSDs, paste, collapse = ", ", character(1L))
write.csv(df.Normal.Simulation, file = "NormalSimulationOutcomeBias1.csv")

############# Methods bias

registerDoParallel(c1)

r <- foreach (m = 1:Reps, .combine=rbind, .packages = c("data.table"), .export = c("Studies", "Subj", "True.sd",
                                                                                   "theta", "tau.sq")) %dopar% {
  
  # Size of per unit bias increase
  Bias.multiple <- -0.2
  
  
  ID.this.loop <- integer()
  for (i in Subj){
    
    for (j in True.sd){
      
      for (k in theta){
        
        for (l in tau.sq){
          
          for (n in Studies){
            
            for (o in 1:n){
              new.this.loop <- as.integer((match(i, Subj)-1) * length(True.sd) * length(theta) * length(tau.sq) * sum(Studies) * Reps + 
                                            (match(j, True.sd)-1) * length(theta) * length(tau.sq) * sum(Studies) * Reps +
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
  
  Normal.Simulation <- data.table(
    Unique_ID = c(ID.this.loop),
    Rep_Number = integer(length = ID),
    Rep_Subj = numeric(length = ID),
    Rep_sd = numeric(length = ID),
    Rep_theta = numeric(length = ID),
    Rep_tau.sq = numeric(length = ID),
    Rep_NumStudies = numeric(length = ID),
    Study_ID = integer(length = ID),
    Study_estimate = numeric(length = ID),
    Study_sd = numeric(length = ID),
    Study_n = integer(length = ID),
    Study_Number.of.biases = integer(length = ID)
  )
  
  
  for (i in Subj){
    
    for (j in True.sd){
      
      for (k in theta){
        
        for (l in tau.sq){
          
          
          
          for (n in Studies){
            
            for (o in 1:n){
              
              counter <- as.integer((match(i, Subj)-1) * length(True.sd) * length(theta) * length(tau.sq) * sum(Studies) * Reps + 
                                      (match(j, True.sd)-1) * length(theta) * length(tau.sq) * sum(Studies) * Reps +
                                      (match(k, theta)-1) * length(tau.sq) * sum(Studies) * Reps +
                                      (match(l, tau.sq)-1) * sum(Studies) * Reps +
                                      (sum(Studies[0:(match(n, Studies)-1)]) + o -1) * Reps +
                                      m
              )
              
              #Statement left in case of varying number of subjects later
              Study_patientnumber <- round(rlnorm(1, meanlog = 4.2, sdlog = 1.1)+2)
              
              ## Draw from binomial how many methodological concerns study has
              Number.of.biases <- rbinom(1, 3, 1/(exp(Study_patientnumber^0.15)))
              
              Study_mu <- rnorm(1, mean = k + Number.of.biases*Bias.multiple, sd = sqrt(l))
              Study_values <- rnorm(Study_patientnumber, mean = Study_mu, sd = j)
              Study_mean <- mean(Study_values)
              Study_StanDev <- sd(Study_values)
              
              Normal.Simulation[Unique_ID == counter, `:=` (Rep_Number= m, Rep_Subj = i, Rep_sd = j,
                                                            Rep_theta = k, Rep_tau.sq = l, Rep_NumStudies = n,
                                                            Study_ID = o, 
                                                            Study_estimate = Study_mean, Study_sd = Study_StanDev, 
                                                            Study_n = Study_patientnumber,
                                                            Study_Number.of.biases = Number.of.biases)]
              
            }
          }
          
        }
      }
    }
  }
  
  Normal.Simulation
}

write.csv(r, file = "NormalSimulationMethods1.csv")

TimeTaken <- proc.time() - StartTime

stopCluster(c1)
