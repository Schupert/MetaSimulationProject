### Remove previous variables
rm(list = ls())

#### Libraries
library(data.table)

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
Subj = c(3.5, 4)

# sd = study level standard deviation
True.sd = c(2,3)

# theta = population level mean
theta = c(0,5)

# tau.sq = between studies variance (can be squared due to sqrt() in normal draw), ?to be distributed
tau.sq = c(1,2,3)

# ?need to state I.sq in advance?


# ID = total number of data points required, also used as an ID number. WILL NEED UPDATING
ID =  length(Subj) * length(True.sd) * length(theta) * length(tau.sq) * Reps * sum(Studies)

StartTime <- proc.time()

############# Basic normal simulation

### Set up data.table to store results

Normal.Simulation <- data.table(
  Unique_ID = c(1:ID),
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

### Simulate
# Order of loops is based on calculation for ID above
# If number in each study is being drawn from distribution will need extra coding

counter <- integer(1)

for (i in Subj){
  
  for (j in True.sd){
    
    for (k in theta){
      
      for (l in tau.sq){
        
        for (m in 1:Reps){
          
          for (n in Studies){
            
            for (o in 1:n){
              
              #Statement left in case of varying number of subjects later
              Study_patientnumber <- round(rlnorm(1, meanlog = 4.2, sdlog = 1.1)+0.6)
              
              Study_mu <- rnorm(1, mean = k, sd = sqrt(l))
              Study_values <- rnorm(Study_patientnumber, mean = Study_mu, sd = j)
              Study_mean <- mean(Study_values)
              Study_StanDev <- sd(Study_values)
              
              Normal.Simulation[Unique_ID == counter, `:=` (Rep_Number= m, Rep_Subj = i, Rep_sd = j,
                                                            Rep_theta = k, Rep_tau.sq = l, Rep_NumStudies = n,
                                                            Study_ID = o, 
                                                            Study_estimate = Study_mean, Study_sd = Study_StanDev, 
                                                            Study_n = Study_patientnumber)]
              
              counter <- as.integer(counter + 1)
            }
          }
        }
      }
    }
  }
}

write.csv(Normal.Simulation, file = "NormalSimulation1.csv")


############# Begg publication bias

# Set up strength of publication bias selection
Begg_a <- 1.5
Begg_b <- 4
Begg_sided <- 1

### Set up data.table to store results

Normal.Simulation <- data.table(
  Unique_ID = c(1:ID),
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

### Simulate
# Order of loops is based on calculation for ID above
# If number in each study is being drawn from distribution will need extra coding

counter <- 1

for (i in Subj){
  
  for (j in True.sd){
    
    for (k in theta){
      
      for (l in tau.sq){
        
        for (m in 1:Reps){
          
          for (n in Studies){
            
            for (o in 1:n){
              
              #Statement left in case of varying number of subjects later
              Study_patientnumber <- round(rlnorm(1, meanlog = 4.2, sdlog = 1.1)+0.6)
              
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
              
              counter <- counter + 1
            }
          }
        }
      }
    }
  }
}

write.csv(Normal.Simulation, file = "NormalSimulationBegg1.csv")

############# Step publication bias

# Boundary of step function on p value, causing severity of publication bias
Severity.boundary <- c(0.05, 0.2)

### Set up data.table to store results

Normal.Simulation <- data.table(
  Unique_ID = c(1:ID),
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

### Simulate
# Order of loops is based on calculation for ID above
# If number in each study is being drawn from distribution will need extra coding

counter <- 1

for (i in Subj){
  
  for (j in True.sd){
    
    for (k in theta){
      
      for (l in tau.sq){
        
        for (m in 1:Reps){
          
          for (n in Studies){
            
            for (o in 1:n){
              
              #Statement left in case of varying number of subjects later
              Study_patientnumber <- round(rlnorm(1, meanlog = 4.2, sdlog = 1.1)+0.6)
              
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
              
              counter <- counter + 1
            }
          }
        }
      }
    }
  }
}

write.csv(Normal.Simulation, file = "NormalSimulationStepBias1.csv")

############# Outcome bias

# Set up within study reporting bias - this is now one sided
Tested.outcomes <- 10
Chosen.outcomes <- 1

### Set up data.table to store results

Normal.Simulation <- data.table(
  Unique_ID = c(1:ID),
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

### Simulate
# Order of loops is based on calculation for ID above
# If number in each study is being drawn from distribution will need extra coding


counter <- 1

for (i in Subj){
  
  for (j in True.sd){
    
    for (k in theta){
      
      for (l in tau.sq){
        
        for (m in 1:Reps){
          
          for (n in Studies){
            
            for (o in 1:n){
              
              #Statement left in case of varying number of subjects later
              Study_patientnumber <- round(rlnorm(1, meanlog = 4.2, sdlog = 1.1)+0.6)
              
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
              
              counter <- counter + 1
            }
          }
        }
      }
    }
  }
}

#write.csv(Normal.Simulation, file = "NormalSimulationOutcomeBias1.csv")
df.Normal.Simulation <- as.data.frame(Normal.Simulation)
df.Normal.Simulation$Study_rejectedMeans <- vapply(df.Normal.Simulation$Study_rejectedMeans, paste, collapse = ", ", character(1L))
df.Normal.Simulation$Study_rejectedSDs <- vapply(df.Normal.Simulation$Study_rejectedSDs, paste, collapse = ", ", character(1L))
write.csv(df.Normal.Simulation, file = "NormalSimulationOutcomeBias1.csv")

############# Methods bias

# Size of per unit bias increase
Bias.multiple <- -0.2

### Set up data.table to store results

Normal.Simulation <- data.table(
  Unique_ID = c(1:ID),
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

### Simulate
# Order of loops is based on calculation for ID above
# If number in each study is being drawn from distribution will need extra coding

counter <- 1

for (i in Subj){
  
  for (j in True.sd){
    
    for (k in theta){
      
      for (l in tau.sq){
        
        for (m in 1:Reps){
          
          for (n in Studies){
            
            for (o in 1:n){
              
              #Statement left in case of varying number of subjects later
              Study_patientnumber <- round(rlnorm(1, meanlog = 4.2, sdlog = 1.1)+0.6)
              
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
              
              counter <- counter + 1
            }
          }
        }
      }
    }
  }
}

write.csv(Normal.Simulation, file = "NormalSimulationMethods1.csv")

TimeTaken <- proc.time() - StartTime