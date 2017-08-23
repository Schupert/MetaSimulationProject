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



library(doParallel)
library(foreach)
num.Cores <- detectCores() - 4
c1 <- makeCluster(num.Cores)
registerDoParallel(c1)

# Reps = number of repetitions of experiment
Reps = 100

StartTime <- proc.time()
r <- foreach (m = 1:Reps, .combine=rbind, .packages = c("data.table")) %dopar% {
  
  
  
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
  
  # ID = total number of data points required, also used as an ID number. WILL NEED UPDATING
  #ID =  length(Subj) * length(True.sd) * length(theta) * length(tau.sq) * sum(Studies) * Reps
  
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
  
#   Normal.Simulation <- data.table(
#     Unique_ID = c(1:ID),
#     Rep_Number = integer(length = ID),
#     Rep_Subj = numeric(length = ID),
#     Rep_sd = numeric(length = ID),
#     Rep_theta = numeric(length = ID),
#     Rep_tau.sq = numeric(length = ID),
#     Rep_NumStudies = numeric(length = ID),
#     Study_ID = integer(length = ID),
#     Study_estimate = numeric(length = ID),
#     Study_sd = numeric(length = ID),
#     Study_n = integer(length = ID)
#   )
  

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
              # ID =  length(Subj) * length(True.sd) * length(theta) * length(tau.sq) * sum(Studies) * Reps
              
              
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
#Normal.Simulation <- r[Study_ID != 0 & Rep_sd != 0]

TimeTaken <- proc.time() - StartTime

write.csv(Normal.Simulation, file = "NormalSimulation1.csv")



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

TimeTaken <- proc.time() - StartTime




















#### Early playing

library(doParallel)
library(foreach)
c1 <- makeCluster(4)
registerDoParallel(c1)

x <- foreach(i=1:Reps, .combine='rbind') %do% rnorm(4)
x


##### Further parallel checks

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
Reps = 10

# k = number of studies in series
Studies = c(2,3)

# subj = number of subjects in study, likely to be distributed
Subj = c(3.5, 4)

# sd = study level standard deviation
True.sd = c(2,3)

# theta = population level mean
theta = c(0,5)

# tau.sq = between studies variance (can be squared due to sqrt() in normal draw), ?to be distributed
tau.sq = c(1,2)

# ?need to state I.sq in advance?


# ID = total number of data points required, also used as an ID number. WILL NEED UPDATING
ID =  length(Subj) * length(True.sd) * length(theta) * length(tau.sq) * Reps * sum(Studies)


counter <- integer()

for (i in Subj){
  
  for (j in True.sd){
    
    for (k in theta){
      
      for (l in tau.sq){
        
        for (m in 1:Reps){
          
          for (n in Studies){
            
            for (o in 1:n){
              a <- as.integer((match(i, Subj)-1) * length(True.sd) * length(theta) * length(tau.sq) * sum(Studies) * Reps + 
                                (match(j, True.sd)-1) * length(theta) * length(tau.sq) * sum(Studies) * Reps +
                                (match(k, theta)-1) * length(tau.sq) * sum(Studies) * Reps +
                                (match(l, tau.sq)-1) * sum(Studies) * Reps +
                                (sum(Studies[0:(match(n, Studies)-1)]) + o -1) * Reps +
                                m
              )
              counter<- append(counter, a)
            }
          }
        }
      }
    }
  }
}

sort(counter)
length(counter)
hist(counter)



library(doParallel)
library(foreach)
c1 <- makeCluster(4)
registerDoParallel(c1)

# Reps = number of repetitions of experiment
Reps = 10

StartTime <- proc.time()
r <- foreach (m = 1:Reps, .combine=c, .packages = c("data.table")) %dopar% {
  
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
  
  # ID = total number of data points required, also used as an ID number. WILL NEED UPDATING
  ID =  length(Subj) * length(True.sd) * length(theta) * length(tau.sq) * sum(Studies) * Reps
  
  
  counter <- integer()
  
  for (i in Subj){
    
    for (j in True.sd){
      
      for (k in theta){
        
        for (l in tau.sq){
          
          
          
          for (n in Studies){
            
            #for (o in 1:n){
              a <- as.integer((match(i, Subj)-1) * length(True.sd) * length(theta) * length(tau.sq) * length(Studies) * Reps + 
                                (match(j, True.sd)-1) * length(theta) * length(tau.sq) * length(Studies) * Reps +
                                (match(k, theta)-1) * length(tau.sq) * length(Studies) * Reps +
                                (match(l, tau.sq)-1) * length(Studies) * Reps +
                                (match(n, Studies)-1) * Reps +
                                m
              )
              counter<- append(counter, a)
            #}
          }
          
        }
      }
    }
  }
  counter 
}

#sort(r)
length(r)
hist(r, breaks = c(100))
