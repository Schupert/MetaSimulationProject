rm(list=ls(all=TRUE))

#### Libraries
library(data.table)

#### Set working directory
setwd("D:\\Stats\\AFP\\R Code")

#### set seed for reproduceability
#set.seed(1234)

#### Declare variables

# Reps = number of repetitions of experiment
Reps = 10000

# k = number of studies in series
Studies = c(1)

# subj = number of subjects in study, likely to be distributed
Subj = 100

# controlProp = proportion of total sample in control arm
controlProp = 0.5

# theta = population level log(OR) - this should be considered more purely on the log scale
theta = c(log(1))

# tau.sq = between studies variance (can be squared due to sqrt() in normal draw), ?to be distributed
tau.sq = c(2)

# ?need to state I.sq in advance?

# ID = total number of data points required, also used as an ID number. WILL NEED UPDATING
ID =  length(Subj) * length(controlProp) * length(theta) * length(tau.sq) * Reps * sum(Studies) 


### Function to simulate a single observational study. Control_prop represents the proportion of 
### the total sample that are controls. Pic is the proportion of controls that are expected to die (or alternative poor outcome)
## Draws Pic from uniform distribution [0.05,0.65] as per Sidik Jonkman
## Group1 = Controls, Group2 = Exposed
## Outcome1 = Death, Outcome2 = Life
# 
# Log_Odds_Ratio <- function(StudySize, Log_O_R, Heterogeneity, Control_Prop){
#   StudyLogOR <- rnorm(1, Log_O_R, sqrt(Heterogeneity))
#   Group1Size <- as.integer(Control_Prop*StudySize)
#   Group2Size <- as.integer(StudySize - Group1Size)
#   #Pic <- runif(1, 0.05, 0.65)
#   Pic <- 0.5
#   Group1Out1 <- as.integer(rbinom(1, Group1Size, Pic))
#   Group1Out2 <- as.integer(Group1Size - Group1Out1)
#   Pit <- (Pic * exp(StudyLogOR))/(1 - Pic + (Pic * exp(StudyLogOR)))
#   Group2Out1 <- as.integer(rbinom(1, Group2Size, Pit))
#   Group2Out2 <- as.integer(Group2Size - Group2Out1)
#   if (Group1Out1 == 0){Group1Out1 <- 0.5}
#   if (Group2Out1 == 0){Group2Out1 <- 0.5}
#   if (Group1Out2 == 0){Group1Out2 <- 0.5}
#   if (Group2Out2 == 0){Group2Out2 <- 0.5}
#   return(c(Group1Out1, Group1Out2, Group2Out1, Group2Out2, Group1Size, Group2Size))
# }

### Alternative formulation drawing Pic from a normal distribution with a stated mean and sd, rather than 
### from uniform distribution. ?Unclear if this is a good assumption

# Log_Odds_Ratio_alt <- function(StudySize, Log_O_R, Heterogeneity, Control_Prop, Av_Pic, Sd_Pic){
#   StudyLogOR <- rnorm(1, Log_O_R, sqrt(Heterogeneity))
#   Group1Size <- as.integer(Control_Prop*StudySize)
#   Group2Size <- as.integer(StudySize - Group1Size)
#   Pic <- rnorm(Av_Pic, Sd_Pic)
#   Group1Out1 <- as.integer(rbinom(1, Group1Size, Pic))
#   Group1Out2 <- as.integer(Group1Size - Group1Out1)
#   Pit <- (Pic * exp(StudyLogOR))/(1 - Pic + (Pic * exp(StudyLogOR)))
#   Group2Out1 <- as.integer(rbinom(1, Group2Size, Pit))
#   Group2Out2 <- as.integer(Group2Size - Group2Out1)
#   if (Group1Out1 == 0){Group1Out1 <- 0.5}
#   if (Group2Out1 == 0){Group2Out1 <- 0.5}
#   if (Group1Out2 == 0){Group1Out2 <- 0.5}
#   if (Group2Out2 == 0){Group2Out2 <- 0.5}
#   return(c(Group1Out1, Group1Out2, Group2Out1, Group2Out2, Group1Size, Group2Size))
# }

### Alternative formulation where rather than a set Pic, theta is split between increasing PT and reducing PC and 
# allocation changed so if there is any 0 cells, 0.5 is added to all cells

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

Log_Odds_Ratio_var <- function(StudySize, Log_O_R, Heterogeneity, Control_Prop, mu){
  dummy1 <- rnorm(1, Log_O_R, sqrt(Heterogeneity))
  Group1Size <- as.integer(Control_Prop*StudySize)
  Group2Size <- as.integer(StudySize - Group1Size)
  ### Add individual level variance
  dummy2 <- rnorm(Group1Size, dummy1, sqrt(0.5)*0.1)
  dummy3 <- rnorm(Group2Size, dummy1, sqrt(0.5)*0.1)
  StudyLogOR1 <- rnorm(Group1Size, dummy2, sqrt(0.5)*0.1)
  StudyLogOR2 <- rnorm(Group2Size, dummy3, sqrt(0.5)*0.1)
  
  Pic <- exp(mu - 0.5*StudyLogOR1) / (1 + exp(mu - 0.5*StudyLogOR1))
  Group1Out1 <- as.integer(sum(rbinom(Group1Size, 1, Pic)))
  Group1Out2 <- as.integer(Group1Size - Group1Out1)
  Pit <- exp(mu + 0.5*StudyLogOR2)  / (1 + exp(mu + 0.5*StudyLogOR2))
  Group2Out1 <- as.integer(sum(rbinom(Group2Size, 1, Pit)))
  Group2Out2 <- as.integer(Group2Size - Group2Out1)
  if (Group1Out1 == 0 | Group2Out1 == 0 | Group1Out2 == 0 | Group2Out2 == 0){
    Group1Out1 <- Group1Out1 + 0.5
    Group2Out1 <- Group2Out1 + 0.5
    Group1Out2 <- Group1Out2 + 0.5
    Group2Out2 <- Group2Out2 + 0.5
  }
  return(c(Group1Out1, Group1Out2, Group2Out1, Group2Out2, Group1Size, Group2Size))
}


### Set up data.table to store results

LogOR.Simulation <- data.table(
  Unique_ID = c(1:ID),
  Rep_Number = integer(length = ID),
  Rep_Subj = integer(length = ID),
  Rep_control_proportion = numeric(length = ID),
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

### Simulate
# Order of loops is based on calculation for ID above
# If number in each study is being drawn from distribution will need extra coding

StartTime <- proc.time()

counter <- 1

for (i in Subj){
  
  for (j in controlProp){
    
    for (k in theta){
      
      for (l in tau.sq){
        
        for (m in 1:Reps){
          
          for (n in Studies){
            
            for (o in 1:n){
              
              #Statement left in case of varying number of subjects later
              Study_patientnumber <- i
              
              # Currently using alternative formulation with extra mu
              # x <- Log_Odds_Ratio(Study_patientnumber, k, l, j, mu= 0.5)
              
              # Checking change in variance with individual level risk
              x <- Log_Odds_Ratio_var(Study_patientnumber, k, l, j, mu= 0.5)
              
              LogOR.Simulation[Unique_ID == counter, `:=` (Rep_Number= m, Rep_Subj = i, Rep_control_proportion = j,
                                                           Rep_theta = k, Rep_tau.sq = l, Rep_NumStudies = n,
                                                           Study_ID = o, 
                                                           Study_estimate = log((x[3]/x[4])/(x[1]/x[2])), 
                                                           Study_sd = sqrt(1/x[1] + 1/x[2] + 1/x[3] + 1/x[4]), 
                                                           Study_n = Study_patientnumber,
                                                           Group1Outcome1 = x[1], Group1Outcome2 = x[2],
                                                           Group2Outcome1 = x[3], Group2Outcome2 = x[4],
                                                           Group1Size = x[5], Group2Size = x[6]
              )]
              
              counter <- counter + 1
            }
          }
        }
      }
    }
  }
}

#write.csv(LogOR.Simulation, file = "LogORSimulation1.csv")

TimeTaken <- proc.time() - StartTime

hist(LogOR.Simulation$Study_estimate)
mean(LogOR.Simulation$Study_estimate)
sd(LogOR.Simulation$Study_estimate)
mean(LogOR.Simulation$Study_sd)
sd(LogOR.Simulation$Study_sd)
hist(LogOR.Simulation$Study_sd)
