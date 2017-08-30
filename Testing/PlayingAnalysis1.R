### Remove previous variables
rm(list = ls())

#### Libraries
library(data.table)
library(doParallel)
library(foreach)
library(metafor)

## Set working directory

Normal.Simulation <- read.csv("NormalSimulation1.csv")
Normal.Simulation <- data.table(Normal.Simulation)
setkey(Normal.Simulation, "Rep_Number", "Rep_Subj", "Rep_sd", "Rep_theta", "Rep_tau.sq", "Rep_NumStudies")

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

####### Normal simulation analysis with no bias

registerDoParallel(c1)

r <- foreach (m = 1:Reps, 
              .combine=rbind, .packages = c("data.table"), 
              .export = c("Studies", "Subj", "controlProp", "theta", "tau.sq", "True.sd", "Normal.Simulation")
) %dopar% {
  
  ID.this.loop <- integer()
  
  for (i in Subj){
    
    for (j in True.sd){
      
      for (k in theta){
        
        for (l in tau.sq){
          
          for (n in Studies){
            
            new.this.loop <- as.integer((match(i, Subj)-1) * length(True.sd) * length(theta) * length(tau.sq) * length(Studies) * Reps + 
                                          (match(j, True.sd)-1) * length(theta) * length(tau.sq) * length(Studies) * Reps +
                                          (match(k, theta)-1) * length(tau.sq) * length(Studies) * Reps +
                                          (match(l, tau.sq)-1) * length(Studies) * Reps +
                                          (match(n, Studies)-1) * Reps +
                                          m
            )
            
            ID.this.loop <- append(ID.this.loop, new.this.loop)
            
          }
        }
      }
    }
  }
  
  ID <- length(ID.this.loop)
  
  Normal.Sim.Results <- data.table(
    Unique_ID = c(ID.this.loop),
    Rep_Number = integer(length = ID),
    Rep_Subj = numeric(length = ID),
    Rep_sd = numeric(length = ID),
    Rep_theta = numeric(length = ID),
    Rep_tau.sq = numeric(length = ID),
    Rep_NumStudies = numeric(length = ID),
    FE_Estimate = numeric(length = ID),
    FE_Est_Low_CI = numeric(length = ID),
    FE_Est_Up_CI = numeric(length = ID),
    REML_Estimate = numeric(length = ID),
    REML_Est_Low_CI = numeric(length = ID),
    REML_Est_Up_CI = numeric(length = ID),
    REML_tau2 = numeric(length = ID),
    REML_I2 = numeric(length = ID)
  )
  
  for (i in Subj){
    
    for (j in True.sd){
      
      for (k in theta){
        
        for (l in tau.sq){
          
          for (n in Studies){
            
            ### Temporary data.table
            temp.data <- Normal.Simulation[J(m, i, j, k, l, n)]
            
            ma.fe <- rma.uni(yi = temp.data$Study_estimate, vi = temp.data$Rep_sd^2, method = "FE")
            
            ma.reml <- rma.uni(yi = temp.data$Study_estimate, vi = temp.data$Rep_sd^2, method = "REML")
            
            Normal.Sim.Results[Unique_ID == counter, `:=` (Rep_Number= m, Rep_Subj = i, Rep_sd = j,
                                                          Rep_theta = k, Rep_tau.sq = l, Rep_NumStudies = n,
                                                          FE_Estimate = ma.fe[1],
                                                          FE_Est_Low_CI = ma.fe[5],
                                                          FE_Est_Up_CI = ma.fe[6],
                                                          REML_tau2 = ma.reml[8],
                                                          REML_I2 = ma.reml[22])]
            
          }
        }
      }
    }
  }
  Normal.Sim.Results
}
            