# 
#
#
#
#
#
# CURRENTLY ALTERED TO HAVE UNORDERED OUTCOME IN UMD to allow estimate of covariance

## Remove variables
rm(list = ls())

#### Libraries
library(data.table)
library(metafor)
library(copula)

#### Set working directory - this is a problem with transferring code around computers. 
####    Current system works when opened using R project - sets wd to project wd/Results
#setwd("D:\\Stats\\AFP\\R Code")
if (length(regmatches(getwd(), gregexpr("/Results", getwd()))) == 0) 
{workingDirectory <- paste(getwd(), "/Results", sep = "")}
if (length(regmatches(workingDirectory, gregexpr("/Results", workingDirectory)))) {setwd(workingDirectory)}

#### set seed for reproduceability
#set.seed(1234)

#### Declare variables

# Reps = number of repetitions of experiment
Reps = 10

# k = number of studies in series
Studies = c(1)

# subj = number of subjects in study, likely to be distributed
Subj = 100

# sd = study level standard deviation
True.sd = 1

# theta = population level mean
theta = 0

# tau.sq = between studies variance (can be squared due to sqrt() in normal draw), ?to be distributed
tau.sq = c(0)

# ?need to state I.sq in advance?

Tested.outcomes <- 100
Chosen.outcomes <- 1
Sd.split <- 0.9

# controlProp = proportion of total sample in control arm
controlProp = 0.5

### UMD function with multiple outcome bias with frac being sd in first level, num.times = number of outcomes simulated
# outputs vectors ordered by p-val

UMD.mult.out <- function(StudySize, Theta, Heterogeneity, Control_Prop, total.sd, frac, num.times){
  StudyUMD <- rnorm(1, Theta, sqrt(Heterogeneity))
  Group1Size <- as.integer(Control_Prop*StudySize)
  Group2Size <- Group1Size
  z <- normalCopula(param = frac, dim = Group1Size)
  Z <- rCopula(num.times, z)
  ControlGroup <- qnorm(Z, mean = -StudyUMD/2, sd = total.sd)
  y <- normalCopula(param = frac, dim = Group1Size)
  Y <- rCopula(num.times, y)
  TreatmentGroup <- qnorm(Y, mean = StudyUMD/2, sd = total.sd)
  Studymean <- apply(TreatmentGroup,1,mean) - apply(ControlGroup, 1, mean)
  Studysd <- sqrt( (apply(TreatmentGroup, 1, var) * (Group1Size - 1) + apply(TreatmentGroup, 1, var) * (Group2Size-1))/ (Group1Size + Group2Size -2) )
  Begg_p <- pnorm(-Studymean/Studysd)
  return(list(Studymean[order(Begg_p)], Studysd[order(Begg_p)]))
}



# ID = total number of data points required, also used as an ID number. WILL NEED UPDATING
ID =  length(Subj) * length(True.sd) * length(theta) * length(tau.sq) * Reps * sum(Studies) 

### Set up data.table to store results

Normal.Simulation <- data.table(
  Unique_ID = c(1:ID),
  Rep_Number = integer(length = ID),
  Rep_Subj = integer(length = ID),
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

StartTime <- proc.time()

counter <- 1

for (i in Subj){
  
  for (j in True.sd){
    
    for (k in theta){
      
      for (l in tau.sq){
        
        for (m in 1:Reps){
          
          for (n in Studies){
            
            for (o in 1:n){
              
              #Statement left in case of varying number of subjects later
              Study_patientnumber <- round(rlnorm(1, meanlog = 4.2, sdlog = 1.1)+4)
              
              ### Implement Within study multiple outcomes bias - split variance by simulating values
              # for each individual which represent the between person variability, then sample from these
              # with sd to simulate testing variability
              
              
              # Sample multiple outcome measures from same set of patients, using Person_values as mean
              Study_values <- UMD.mult.out(Study_patientnumber, k, l, controlProp, True.sd, Sd.split, Tested.outcomes)
              
              Study_mean <- Study_values[[1]]
              Study_StanDev <- Study_values[[2]]
              
              ### Altered to give full list to calculate covariance
              
              Normal.Simulation[Unique_ID == counter, `:=` (Rep_Number= m, Rep_Subj = i, Rep_sd = j,
                                                            Rep_theta = k, Rep_tau.sq = l, Rep_NumStudies = n,
                                                            Study_ID = o, 
                                                            Study_estimate = Study_mean[1], 
                                                            Study_sd = Study_StanDev[1], 
                                                            Study_n = Study_patientnumber,
                                                            #### should be [-1]
                                                            Study_rejectedMeans = list(Study_mean),
                                                            Study_rejectedSDs = list(Study_StanDev)
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

TimeTaken <- proc.time() - StartTime

hist(Normal.Simulation$Study_estimate)
mean(Normal.Simulation$Study_estimate)
sd(Normal.Simulation$Study_estimate)
mean(Normal.Simulation$Study_sd)
hist(Normal.Simulation$Study_sd)

plot(Normal.Simulation$Study_estimate ~ Normal.Simulation$Study_n)
plot(Normal.Simulation$Study_sd ~ Normal.Simulation$Study_estimate)

testRes <- rma.uni(Study_estimate, Study_sd^2, data=Normal.Simulation, method="FE")
funnel(testRes)


#### Problems with list assignment
Normal.Simulation[Unique_ID == 1000, Study_rejectedSDs := Study_StanDev[-lv]]
typeof(Normal.Simulation[1000]$Study_rejectedSDs)

DT = data.table(a=1:2, b=list(1:5,1:10))
DT


sapply(DT$b, length)
typeof(DT$b)
typeof(DT[1]$b)

#### Testing correlation
Theta <- 0
Heterogeneity <- 2
Control_Prop <- 0.5
StudySize <- 100000
frac <- 0.5
total.sd <- 3
num.times <- 5


StudyUMD <- rnorm(1, Theta, sqrt(Heterogeneity))
Group1Size <- as.integer(Control_Prop*StudySize)
Group2Size <- as.integer(StudySize - Group1Size)
ControlGroup1 <- rnorm(Group1Size, 0, sqrt(frac) * total.sd)
TreatmentGroup1 <- rnorm(Group2Size, mean = StudyUMD, sqrt(frac) * total.sd)
ControlGroupAll <- replicate(num.times, rnorm(Group1Size, ControlGroup1, sqrt(1-frac) * total.sd), simplify = FALSE)
TreatmentGroupAll <- replicate(num.times, rnorm(Group2Size, TreatmentGroup1, sqrt(1-frac) * total.sd), simplify = FALSE)
Studymean <- sapply(TreatmentGroupAll, mean) - sapply(ControlGroupAll, mean)
Studysd <- sqrt( sapply(ControlGroupAll, var)/Group1Size + sapply(TreatmentGroupAll, var)/Group2Size )
Begg_p <- pnorm(-Studymean/Studysd)

output <- matrix(unlist(TreatmentGroupAll), ncol = 5)
cor(output)
output2 <- matrix(unlist(ControlGroupAll), ncol = 5)
cor(output2)

### Correlation of total simulation outputs ? wrong set up
# output <- matrix(unlist(Normal.Simulation$Study_rejectedMeans), ncol = 1000)
# asdf <- cor(output)
# mean(asdf[lower.tri(asdf)])

#### Likely correct correlation set up, estimates show equal to first level sd plit = sd.split
output <- matrix(unlist(Normal.Simulation$Study_rejectedMeans), ncol = Tested.outcomes, byrow = TRUE)
dim(output)
cor(output[1,], output[10,])
asdf <-cor(output)
mean(asdf[lower.tri(asdf)])
hist(asdf[lower.tri(asdf)])
plot(output[,8], output[,6])

### ? Is this correct - likely not
output2 <- matrix(unlist(Normal.Simulation$Study_rejectedMeans), ncol = Reps, byrow = FALSE)
asdf2 <- cor(output2)
mean(asdf2[lower.tri(asdf2)])
plot(output2[,8], output2[,6])
