######## Half way through coding test to see how much quicker this is than the TestProfiling
#
#Remaining Jobs =
# 1. Set up data.table outside loop fully
# 2. Vectorise if possible to complete all repeats in 1 loop
# 3. Even if not possible, do simple assign at end, rather than == asign
# 4. Find fast way to merge outside of reps loop

# counter.func <- function(i, k, l, n, m){
#   as.integer((apply(sapply(Subj, function(vec) {i %in% vec}), 1, which.max)[1]-1) * length(theta) * length(tau.sq) * sum(Studies) * Reps +
#                (match(k, theta)-1) * length(tau.sq) * sum(Studies) * Reps +
#                (match(l, tau.sq)-1) * sum(Studies) * Reps +
#                (sum(Studies[0:(match(n, Studies)-1)]) + o -1) * Reps +
#                m
#   )
#   
# }
# 
# library(compiler)
# library(microbenchmark)
# 
# g <- cmpfun(counter.func)
# 
# 
# 
# mbm = microbenchmark(
#   rcode = counter.func(c(30,40), 0, 1, 3, 1),
#   ccode = g(c(30,40), 0, 1, 3, 1)
# )
# 
# library(ggplot2)
# 
# autoplot(mbm)

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
Reps = 10

# k = number of studies in series
Studies = c(3,5,10,30,50,100)
#Studies = c(3,5,10,30)

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


# ID = total number of data points required, also used as an ID number. WILL NEED UPDATING
 ID =  length(Subj) * length(theta) * length(tau.sq) * Reps * sum(Studies)

 vector1 <- 1:ID

Normal.Simulation <- data.table(
   Unique_ID = vector1,
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
 
 
 # base.table <- data.table(
 #   Unique_ID = vector1,
 #   Rep_Number = integer(length = ID),
 #   Rep_Subj = list(length = ID),
 #   Rep_theta = numeric(length = ID),
 #   Rep_tau.sq = numeric(length = ID),
 #   Rep_NumStudies = numeric(length = ID),
 #   Study_ID = integer(length = ID),
 #   Study_estimate = numeric(length = ID),
 #   Study_sd = numeric(length = ID),
 #   Study_n = integer(length = ID),
 #   Study_rejectedMeans = list(length = ID),
 #   Study_rejectedSDs = list(length = ID),
 #   Study_Number.of.biases = integer(length = ID)
 # )

# ID <- 0
# 
# base.table <- data.table(
#   Unique_ID = integer(),
#   Rep_Number = integer(),
#   Rep_Subj = list(),
#   Rep_theta = numeric(),
#   Rep_tau.sq = numeric(),
#   Rep_NumStudies = numeric(),
#   Study_ID = integer(),
#   Study_estimate = numeric(),
#   Study_sd = numeric(),
#   Study_n = integer(),
#   Study_rejectedMeans = list(),
#   Study_rejectedSDs = list(),
#   Study_Number.of.biases = integer()
# )

list.of.sample.sizes <- c(round(runif(ID/4, 30, 40)), rep(60, times = ID/4), round(runif(ID/4, 250, 1000)), round(rlnorm(ID/4, meanlog = 4.2, sdlog = 1.1) + 4) )

StartTime <- proc.time()



for (i in Subj){
  
  for (k in theta){
    
    for (l in tau.sq){
      
      for (n in Studies){
        
        for (o in 1:n){
          
          # counter <- as.integer((apply(sapply(Subj, function(vec) {i %in% vec}), 1, which.max)[1]-1) * length(theta) * length(tau.sq) * sum(Studies) * Reps +
          #                         (match(k, theta)-1) * length(tau.sq) * sum(Studies) * Reps +
          #                         (match(l, tau.sq)-1) * sum(Studies) * Reps +
          #                         (sum(Studies[0:(match(n, Studies)-1)]) + o -1) * Reps
          # )
          # 
          # vector1 <- c((counter+1):(counter+Reps))
          # 
          # ID <- length(vector1)
          # 
          # Normal.Simulation <- data.table(
          #   Unique_ID = vector1,
          #   Rep_Number = integer(length = ID),
          #   Rep_Subj = list(length = ID),
          #   Rep_theta = numeric(length = ID),
          #   Rep_tau.sq = numeric(length = ID),
          #   Rep_NumStudies = numeric(length = ID),
          #   Study_ID = integer(length = ID),
          #   Study_estimate = numeric(length = ID),
          #   Study_sd = numeric(length = ID),
          #   Study_n = integer(length = ID),
          #   Study_rejectedMeans = list(length = ID),
          #   Study_rejectedSDs = list(length = ID),
          #   Study_Number.of.biases = integer(length = ID)
          # )
          
          for (m in 1:Reps){
            
          
          counter <- as.integer((apply(sapply(Subj, function(vec) {i %in% vec}), 1, which.max)[1]-1) * length(theta) * length(tau.sq) * sum(Studies) * Reps +
                                  (match(k, theta)-1) * length(tau.sq) * sum(Studies) * Reps +
                                  (match(l, tau.sq)-1) * sum(Studies) * Reps +
                                  (sum(Studies[0:(match(n, Studies)-1)]) + o -1) * Reps +
                                  m
          )
            
#counter <- g(i, k, l, n, m)
          
          ##### No bias
          
          #Select sample size
          if (is.integer(i[1]) == TRUE){
            Study_patientnumber <- round(runif(1, i[1], i[2]))
          } else {
            Study_patientnumber <- round(rlnorm(1, meanlog = 4.2, sdlog = 1.1) + 4)
          }
          
          # Study_patientnumber <- list.of.sample.sizes[counter]
          
          # ifelse (is.integer(i[1]) == TRUE, 
          #         Study_patientnumber <- round(runif(1, i[1], i[2])), 
          #         round(rlnorm(1, meanlog = 4.2, sdlog = 1.1) + 4) )
          
          Study_summary <- UMD(Study_patientnumber, k, l, controlProp, True.sd)
          Study_mean <- Study_summary[1]
          Study_StanDev <- Study_summary[2]
          
          # Normal.Simulation[m , `:=` (Unique_ID = (counter + m), Rep_Number= m, Rep_Subj = list(i),
          #                                                                     Rep_theta = k, Rep_tau.sq = l, Rep_NumStudies = n,
          #                                                                     Study_ID = o,
          #                                                                     Study_estimate = Study_mean, Study_sd = Study_StanDev,
          #                                                                     Study_n = Study_patientnumber)]
          
          # Normal.Simulation[Unique_ID == counter, `:=` (Rep_Number= m, Rep_Subj = list(i),
          #                                                                     Rep_theta = k, Rep_tau.sq = l, Rep_NumStudies = n,
          #                                                                     Study_ID = o,
          #                                                                     Study_estimate = Study_mean, Study_sd = Study_StanDev,
          #                                                                     Study_n = Study_patientnumber)]
          
          # Normal.Simulation[counter, `:=` (Rep_Number= m, Rep_Subj = list(i),
          #                                               Rep_theta = k, Rep_tau.sq = l, Rep_NumStudies = n,
          #                                               Study_ID = o,
          #                                               Study_estimate = Study_mean, Study_sd = Study_StanDev,
          #                                               Study_n = Study_patientnumber)]
          
          Normal.Simulation[counter, `:=` (Study_estimate = Study_mean, Study_sd = Study_StanDev,
                                           Study_n = Study_patientnumber)]

          
          }
        #base.table <- rbindlist(list(base.table, Normal.Simulation), use.names = TRUE, fill = TRUE)
          #base.table[Unique_ID == ]
      }
    }
  }
}

}

TimeTaken <- proc.time() - StartTime
