## Remove variables
rm(list = ls())

#### Libraries
library(data.table)
library(metafor)
library(copula)

#### Set working directory - this is a problem with transferring code around computers. 
####    Current system works when opened using R project - sets wd to project wd/Results
#setwd("D:\\Stats\\AFP\\R Code")
if (lengths(regmatches(getwd(), gregexpr("/Results", getwd()))) == 0) 
{workingDirectory <- paste(getwd(), "/Results", sep = "")}
if (lengths(regmatches(workingDirectory, gregexpr("/Results", workingDirectory)))) {setwd(workingDirectory)}

#### set seed for reproduceability
set.seed(1234)

#### Declare variables

# Reps = number of repetitions of experiment
Reps = 1000

# k = number of studies in series
Studies = c(1)

# subj = number of subjects in study, likely to be distributed
#Subj = list(c(100,100))
Subj = 100

# controlProp = proportion of total sample in control arm
controlProp = 0.5

# theta = population level log(OR) - this should be considered more purely on the log scale
theta = c(log(1))

# tau.sq = between studies variance (can be squared due to sqrt() in normal draw), ?to be distributed
tau.sq = c(0)

# Frequency of event averaged across 2 arms (before applying change due to theta) = EvFreq
EvFreq = c(0.5)

# ?need to state I.sq in advance?

# Set up within study reporting bias - this is now one sided
Tested.outcomes <- 100
Chosen.outcomes <- 1
Sd.split <- 0.8

# ID = total number of data points required, also used as an ID number. WILL NEED UPDATING
ID =  length(Subj) * length(controlProp) * length(theta) * length(tau.sq) * length(EvFreq) * Reps * sum(Studies) 


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
#   Pic <- runif(1, 0.05, 0.65)
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

LogOR_mult_out <- function(StudySize, Theta, Heterogeneity, Control_Prop, mu, frac, num.times){
  StudyLogOR <- rnorm(1, Theta, sqrt(Heterogeneity))
  Group1Size <- as.integer(Control_Prop*StudySize)
  Group2Size <- Group1Size
  EventFreq <- log(mu / (1 - mu))
  Pic <- exp(EventFreq - 0.5*StudyLogOR) / (1 + exp(EventFreq - 0.5*StudyLogOR))
  Pit <- exp(EventFreq + 0.5*StudyLogOR)  / (1 + exp(EventFreq + 0.5*StudyLogOR))
  ### This doesn't create correlated outcomes
  z <- normalCopula(param = frac, dim = num.times)
  Z <- rCopula(Group1Size, z)
  ControlGroup <- qbinom(Z, size=1, prob=Pic)
  CGO1 <- apply(ControlGroup, 2, sum)
  CGO2 <- Group1Size - CGO1
  
  y <- normalCopula(param = frac, dim = num.times)
  Y <- rCopula(Group1Size, y)
  TreatmentGroup <- qbinom(Y,size=1,prob=Pit)
  TGO1 <- apply(TreatmentGroup, 2, sum)
  TGO2 <- Group2Size - TGO1
  
  o <- data.table(CGO1, CGO2, TGO1, TGO2)
  
#   if (CGO1 == 0 | TGO1 == 0 | CGO2 == 0 | TGO2 == 0){
#     CGO1 <- CGO1 + 0.5
#     TGO1 <- TGO1 + 0.5
#     CGO2 <- CGO2 + 0.5
#     TGO2 <- TGO2 + 0.5
#   }
  
  o[CGO1 == 0] <- o[CGO1 == 0] + 0.5
  o[CGO2 == 0] <- o[CGO2 == 0] + 0.5
  o[TGO1 == 0] <- o[TGO1 == 0] + 0.5
  o[TGO2 == 0] <- o[TGO2 == 0] + 0.5
  
  ############### Can use same selection mechanism but need to solve problem of what happens with 0 celled results
  
  Study.est <- log((o$TGO1/o$TGO2) / (o$CGO1/o$CGO2))
  Study.se <- sqrt(1/o$TGO1 + 1/o$TGO2 + 1/o$CGO1 + 1/o$CGO2)
  Study.p.val <- pnorm(Study.est/Study.se)
  
  #return(c(o$CGO1[order(Study.p.val)], o$CGO2[order(Study.p.val)], o$TGO1[order(Study.p.val)], o$TGO2[order(Study.p.val)], Group1Size, Group2Size))
#return(c(o[order(Study.p.val)], Group1Size))
  ## Unordered to check corr
  return(c(o, Group1Size, list(Study.est) ))
  }


### Set up data.table to store results

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
  Study_rej_G1O1 = list(length = ID),
  Study_rej_G2O1 = list(length = ID),
  Study_Number.of.biases = integer(length = ID),
  Study_rej_means = list(length = ID)
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
              
              for (p in EvFreq){
                
                #Statement left in case of varying number of subjects later
                #Study_patientnumber <- round(rlnorm(1, meanlog = 4.6, sdlog = 1))
                
                
                ### Implement Within study multiple outcomes bias - split variance by simulating values
                # for each individual which represent the between person variability, then sample from these
                # with sd to simulate testing variability
                
                #                   x <- Log_Odds_Ratio(Study_patientnumber, k, l, j, p)
                #                   
                #                   Study_mean <- log((x[3]/x[4])/(x[1]/x[2]))
                #                   Study_StanDev <- 1/x[1] + 1/x[2] + 1/x[3] + 1/x[4]
                #                   
                #                   Begg_p <- pnorm(- Study_mean/(Study_StanDev))
                #                   
                #                   Step_weight <- ifelse(Begg_p < Severity.boundary[1], 1, ifelse(Begg_p < Severity.boundary[2], 0.75, 0.25))
                #                   
                
                
#                 # In the log(OR) condition tau2 variance is split between study level and person level
#                 Study_mu <- rnorm(1, mean = k, sd = sqrt(0.5 * l))
#                 #Person_values <- rnorm(Study_patientnumber, mean = Study_mu, sd = sqrt(0.5*l))
#                 
#                 # Sample multiple outcome measures from same set of patients, using Person_values as mean
#                 #Study_values <- replicate(Tested.outcomes, rnorm(Study_patientnumber, mean = Person_values, sd = (1/sqrt(2))*j) )
#                 Study_values <- replicate(Tested.outcomes, Log_Odds_Ratio(Study_patientnumber, Study_mu, 0.5 * l, j, p) )
#                 
#                 Study_mean <- numeric(length = Tested.outcomes)
#                 Study_StanDev <- numeric(length = Tested.outcomes)
#                 Begg_p <- numeric(length = Tested.outcomes)
#                 
#                 for (z in 1:Tested.outcomes) {
#                   Study_mean[z] <- log((Study_values[3,z]/Study_values[4,z])/(Study_values[1,z]/Study_values[2,z]))
#                   Study_StanDev[z] <- sqrt(1/Study_values[1,z] + 1/Study_values[2,z] + 1/Study_values[3,z] + 1/Study_values[4,z])
#                   Begg_p[z] <- pnorm(-Study_mean[z]/(Study_StanDev[z]))
#                 }
#                 
#                 lv <- which.min(Begg_p)
                
                #Select sample size
#                 if (is.integer(i[1]) == TRUE){
#                   Study_patientnumber <- round(runif(1, i[1], i[2]))
#                 } else {
#                   Study_patientnumber <- round(rlnorm(1, meanlog = 4.7, sdlog = 1.2) + 2)
#                 }
                Study_patientnumber = i
                
                x <- LogOR_mult_out(Study_patientnumber, k, l, controlProp, j, Sd.split, Tested.outcomes)
                x.N <- 2*x[[5]]
                
                LogOR.Simulation[counter, `:=` (#Unique_ID = counter,
                                                      Study_G1O1 = x[[1]][1], 
                                                      Study_G2O1 = x[[3]][1], 
                                                      Study_n = x.N,
                                                      Study_rej_G1O1 = list(x[[1]][-1]),
                                                      Study_rej_G2O1 = list(x[[3]][-1]),
                                                      #this is extra
                                                      Study_rej_means = list(x[[6]])
                                                      )]
                
                
                
#                 LogOR.Simulation[Unique_ID == counter, `:=` (Rep_Number= m, Rep_Subj = i, Rep_control_proportion = j,
#                                                              Rep_theta = k, Rep_tau.sq = l, Rep_NumStudies = n,
#                                                              Study_ID = o, 
#                                                              Study_G1O1 = x[1], 
#                                                              Study_G2O1 = x[3], 
#                                                              Study_n = x.N
#                                                              Study_rejectedMeans = list(Study_mean[-lv]),
#                                                              Study_rejectedSDs = list(Study_StanDev[-lv])
#                 )]
#                 
                counter <- counter + 1
              }
            }
          }
        }
      }
    }
  }
}

#write.csv(LogOR.Simulation, file = "LogORSimulationOutcome1.csv")

# if(LogOR.Simulation$Study_G1O1 %% 1 == 0.5){
#   LogOR.Simulation[, Study_estimate := log( (Study_G2O1 / ((Study_n/2 + 1) - Study_G2O1) ) / (Study_G1O1 / ((Study_n/2 + 1) - Study_G1O1) ) )]
#   LogOR.Simulation[, Study_sd := sqrt(1/Study_G1O1 + 1/((Study_n/2 + 1) - Study_G1O1) + 1/Study_G2O1 + 1/((Study_n/2 + 1) - Study_G2O1))]
#   
# } else {
#   LogOR.Simulation[, Study_estimate := log( (Study_G2O1 / ((Study_n/2) - Study_G2O1) ) / (Study_G1O1 / ((Study_n/2) - Study_G1O1) ) )]
#   LogOR.Simulation[, Study_sd := sqrt(1/Study_G1O1 + 1/(Study_n/2 - Study_G1O1) + 1/Study_G2O1 + 1/(Study_n/2 - Study_G2O1))]
#   
# }
# 
# LogOR.Simulation$Study_estimate <- ifelse(LogOR.Simulation$Study_G1O1 %% 1 == 0.5, log( (Study_G2O1 / ((Study_n/2 + 1) - Study_G2O1) ) / (Study_G1O1 / ((Study_n/2 + 1) - Study_G1O1) ) ),
#                                           log( (Study_G2O1 / ((Study_n/2) - Study_G2O1) ) / (Study_G1O1 / ((Study_n/2) - Study_G1O1) ) ))
#        


LogOR.Simulation[, Study_estimate := ifelse(LogOR.Simulation$Study_G1O1 %% 1 == 0.5, 
                                            log( (Study_G2O1 / ((Study_n/2 + 1) - Study_G2O1) ) / (Study_G1O1 / ((Study_n/2 + 1) - Study_G1O1) ) ),
                                            log( (Study_G2O1 / ((Study_n/2) - Study_G2O1) ) / (Study_G1O1 / ((Study_n/2) - Study_G1O1) ) ) )
                 ]

LogOR.Simulation[, Study_sd := ifelse(LogOR.Simulation$Study_G1O1 %% 1 == 0.5, 
                                      sqrt(1/Study_G1O1 + 1/((Study_n/2 + 1) - Study_G1O1) + 1/Study_G2O1 + 1/((Study_n/2 + 1) - Study_G2O1)),
                                      sqrt(1/Study_G1O1 + 1/(Study_n/2 - Study_G1O1) + 1/Study_G2O1 + 1/(Study_n/2 - Study_G2O1)) )
                 ]


TimeTaken <- proc.time() - StartTime

hist(LogOR.Simulation$Study_estimate)
mean(LogOR.Simulation$Study_estimate)
sd(LogOR.Simulation$Study_estimate)
mean(LogOR.Simulation$Study_sd)
hist(LogOR.Simulation$Study_sd)

plot(LogOR.Simulation$Study_estimate ~ LogOR.Simulation$Study_n)
plot(LogOR.Simulation$Study_sd ~ LogOR.Simulation$Study_estimate)

testRes <- rma.uni(Study_estimate, Study_sd^2, data=LogOR.Simulation, method="FE")
funnel(testRes)



## testing splitting variance
a <- rnorm(10000, 0, 2)
b <- rnorm(10000, 0, sqrt(3.5))
c <- rnorm(10000, b, sqrt(0.5))
mean(a)
mean(c)
sd(a)
sd(c)

d <- rnorm(10000, b, sqrt(0.5))
e <- data.frame(d,c)
cor(e)

Log_Odds_Ratio_mean <- function(StudySize, Log_O_R, Heterogeneity, Control_Prop, mu){
  #StudyLogOR <- rnorm(1, Log_O_R, sqrt(Heterogeneity))
  Group1Size <- as.integer(Control_Prop*StudySize)
  Group2Size <- as.integer(StudySize - Group1Size)
  Pic <- exp(mu - 0.5*Log_O_R) / (1 + exp(mu - 0.5*Log_O_R))
  Group1Out1 <- as.integer(rbinom(length(Log_O_R), Group1Size, Pic))
  Group1Out2 <- as.integer(Group1Size - Group1Out1)
  Pit <- exp(mu + 0.5*Log_O_R)  / (1 + exp(mu + 0.5*Log_O_R))
  Group2Out1 <- as.integer(rbinom(length(Log_O_R), Group2Size, Pit))
  Group2Out2 <- as.integer(Group2Size - Group2Out1)
  #   if (Group1Out1 == 0 | Group2Out1 == 0 | Group1Out2 == 0 | Group2Out2 == 0){
  #     Group1Out1 <- Group1Out1 + 0.5
  #     Group2Out1 <- Group2Out1 + 0.5
  #     Group1Out2 <- Group1Out2 + 0.5
  #     Group2Out2 <- Group2Out2 + 0.5
  #   }
  
  #Group1Out1, Group1Out2, Group2Out1, Group2Out2
  Output_logor <- log((Group2Out1/Group2Out2)/(Group1Out1/Group1Out2))
  Output_logsd <- sqrt(1/Group1Out1 + 1/Group1Out2 + 1/Group2Out1 + 1/Group2Out2)
  return(c(Output_logor, Output_logsd))
}

f <- Log_Odds_Ratio_mean(100, a, 0, 0.5, 0.5)
hist(f[1:1000])

g <- Log_Odds_Ratio_mean(100, d, 0, 0.5, 0.5)
hist(g[1:1000])

mean(f[1:1000 & f<1000 & f>-1000])
mean(g[1:1000 & g<1000 & g>-1000])
sd(f[1:1000 & f<1000 & f>-1000])
sd(g[1:1000 & g<1000 & g>-1000])

f <- replicate(3, Log_Odds_Ratio(100, a, 0, 0.5, 0.5))

#### Likely correct correlation set up, estimates show equal to first level sd plit = sd.split

output <- matrix(unlist(LogOR.Simulation$Study_rej_means), ncol = Tested.outcomes, byrow = TRUE)
dim(output)
cor(output[1,], output[10,])
asdf <-cor(output)
mean(asdf[lower.tri(asdf)])
hist(asdf[lower.tri(asdf)])
plot(output[,8], output[,6])
abline(a = 0, b = 1)
