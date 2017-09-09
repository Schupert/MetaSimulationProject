### Remove previous variables
rm(list = ls())

#### Libraries
library(data.table)
library(doParallel)
library(foreach)
library(doRNG)

# Number of cores for parallel
num.Cores <- detectCores() - 1
c1 <- makeCluster(num.Cores)

#### Declare variables

# Reps = number of repetitions of experiment
Reps = 1

# k = number of studies in series
Studies = c(3,5,10,30,50,100)

# subj = number of subjects in study, likely to be distributed
Subj <- list(as.integer(c(60,60)), as.integer(c(30,40)), as.integer(c(250, 1000)), as.numeric(c(4.2, 1.1)))
#Subj <- list(1,2,3,4) for testing

# sd = study level standard deviation
True.sd = 2

# theta = population level mean - need good sense of range for SMD
theta = c(-0.5, 0, 0.5, 1)

# tau.sq = between studies variance (can be squared due to sqrt() in normal draw), ?to be distributed
tau.sq = c(1,2,3)

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

### Unstandardised mean differrence function

UMD <- function(StudySize, Theta, Heterogeneity, Control_Prop, sd){
  StudyUMD <- rnorm(1, Theta, sqrt(Heterogeneity))
  Group1Size <- as.integer(Control_Prop*StudySize)
  Group2Size <- as.integer(StudySize - Group1Size)
  ControlGroup <- rnorm(Group1Size, -StudyUMD/2, sd)
  TreatmentGroup <- rnorm(Group2Size, StudyUMD/2, sd)
  Studymean <- mean(TreatmentGroup) - mean(ControlGroup)
  Studysd <- sqrt( var(ControlGroup)/Group1Size + var(TreatmentGroup)/Group2Size )
  return(c(Studymean, Studysd))
}

### UMD function with multiple outcome bias with frac being sd in first level, num.times = number of outcomes simulated
# outputs vectors ordered by p-val

UMD.mult.out <- function(StudySize, Theta, Heterogeneity, Control_Prop, total.sd, frac, num.times){
  StudyUMD <- rnorm(1, Theta, sqrt(Heterogeneity))
  Group1Size <- as.integer(Control_Prop*StudySize)
  Group2Size <- as.integer(StudySize - Group1Size)
  ControlGroup1 <- rnorm(Group1Size, -StudyUMD/2, sqrt(frac) * total.sd)
  TreatmentGroup1 <- rnorm(Group2Size, mean = StudyUMD/2, sqrt(frac) * total.sd)
  ControlGroupAll <- replicate(num.times, rnorm(Group1Size, ControlGroup1, sqrt(1-frac) * total.sd), simplify = FALSE)
  TreatmentGroupAll <- replicate(num.times, rnorm(Group2Size, TreatmentGroup1, sqrt(1-frac) * total.sd), simplify = FALSE)
  Studymean <- sapply(TreatmentGroupAll, mean) - sapply(ControlGroupAll, mean)
  Studysd <- sqrt( sapply(ControlGroupAll, var)/Group1Size + sapply(TreatmentGroupAll, var)/Group2Size )
  Begg_p <- pnorm(-Studymean/Studysd)
  return(list(Studymean[order(Begg_p)], Studysd[order(Begg_p)]))
}

#### Set working directory - this is a problem with transferring code around computers. 
#
#
#
#

StartTime <- proc.time()

## First parallel loop
#registerDoParallel(c1)
set.seed(123)
#Rprof("test.out")
# r <- foreach (m = 1:Reps, .combine=rbind, .packages = c("data.table"), 
#               .export = c("Studies", "Subj", "True.sd",
#                           "theta", "tau.sq", "controlProp", "UMD", "Severity.boundary", "Begg_a", 
#                           "Begg_b", "Begg_sided", "Tested.outcomes", "Chosen.outcomes", "Sd.split",
#                           "Bias.multiple", "UMD.mult.out")
# ) %dorng% {
for (m in 1:Reps) { 
  
  ID.this.loop <- integer()
  for (i in Subj){
    
    for (k in theta){
      
      for (l in tau.sq){
        
        for (n in Studies){
          
          for (o in 1:n){
            new.this.loop <- as.integer((apply(sapply(Subj, function(vec) {i %in% vec}), 1, which.max)[1]-1) * length(theta) * length(tau.sq) * sum(Studies) * Reps + 
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
  
  ID.per.bias <- length(ID.this.loop)
  Types.of.bias.vector <- c("None", "Step", "PandN", "Outcome", "Method")
  ID <- ID.per.bias * length(Types.of.bias.vector)
  
  Normal.Simulation <- data.table(
    Bias_type = rep(Types.of.bias.vector, each = ID.per.bias),
    Unique_ID = rep(ID.this.loop, times = length(Types.of.bias.vector)),
    Rep_Number = integer(length = ID),
    Rep_Subj = list(length = ID),
    Rep_theta = numeric(length = ID),
    Rep_tau.sq = numeric(length = ID),
    Rep_NumStudies = numeric(length = ID),
    Study_ID = integer(length = ID),
    Study_estimate = numeric(length = ID),
    Study_sd = numeric(length = ID),
    Study_n = integer(length = ID),
    Study_rejectedMeans = list(length = ID),
    Study_rejectedSDs = list(length = ID),
    Study_Number.of.biases = integer(length = ID)
  )
  
  
  for (i in Subj){
    
    for (k in theta){
      
      for (l in tau.sq){
        
        for (n in Studies){
          
          for (o in 1:n){
            
            counter <- as.integer((apply(sapply(Subj, function(vec) {i %in% vec}), 1, which.max)[1]-1) * length(theta) * length(tau.sq) * sum(Studies) * Reps + 
                                    (match(k, theta)-1) * length(tau.sq) * sum(Studies) * Reps +
                                    (match(l, tau.sq)-1) * sum(Studies) * Reps +
                                    (sum(Studies[0:(match(n, Studies)-1)]) + o -1) * Reps +
                                    m
            )
            
            ##### No bias
            
            #Select sample size
            if (is.integer(i[1]) == TRUE){
              Study_patientnumber <- round(runif(1, i[1], i[2]))
            } else {
              Study_patientnumber <- round(rlnorm(1, meanlog = 4.2, sdlog = 1.1) + 4)
            }
            
            Study_summary <- UMD(Study_patientnumber, k, l, controlProp, True.sd)
            Study_mean <- Study_summary[1]
            Study_StanDev <- Study_summary[2]
            
            Normal.Simulation[Unique_ID == counter & Bias_type == "None", `:=` (Rep_Number= m, Rep_Subj = list(i),
                                                                                Rep_theta = k, Rep_tau.sq = l, Rep_NumStudies = n,
                                                                                Study_ID = o, 
                                                                                Study_estimate = Study_mean, Study_sd = Study_StanDev, 
                                                                                Study_n = Study_patientnumber)]
            
            
          }
        }
      }
    }
  }
  
  Normal.Simulation
}

# Need to write list values
#write.csv(r, file = "NormalSimulationAll1.csv")

TimeTaken <- proc.time() - StartTime

#Rprof(NULL)
#summaryRprof("test.out")

stopCluster(c1)