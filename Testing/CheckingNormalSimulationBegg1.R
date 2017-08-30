#### Libraries
library(data.table)
library(metafor)

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
Reps = 1000

# k = number of studies in series
Studies = c(1)

# subj = number of subjects in study, likely to be distributed
Subj = 100

# sd = study level standard deviation
True.sd = 5

# theta = population level mean
theta = 0

# tau.sq = between studies variance (can be squared due to sqrt() in normal draw), ?to be distributed
tau.sq = c(2)

# ?need to state I.sq in advance?

# Set up strength of publication bias selection
Begg_a <- 1.5
Begg_b <- 4
Begg_sided <- 2


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
  Study_n = integer(length = ID)
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
              Study_patientnumber <- round(rlnorm(1, meanlog = 4.2, sdlog = 1.1)+0.5)
              
              ### Implement Begg and Mazumdar publication bias
              repeat{
                
                Study_mu <- rnorm(1, mean = k, sd = sqrt(l))
                Study_values <- rnorm(Study_patientnumber, mean = Study_mu, sd = j)
                Study_mean <- mean(Study_values)
                Study_StanDev <- sd(Study_values)
                
                Begg_weight <-exp(
                  -Begg_b * (
                    (Begg_sided * pnorm(-abs(Study_mean)/(Study_StanDev^0.5))) 
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

#write.csv(Normal.Simulation, file = "NormalSimulationBegg1.csv")

TimeTaken <- proc.time() - StartTime

## Check distribution of results

hist(Normal.Simulation$Study_estimate)
mean(Normal.Simulation$Study_estimate)
sd(Normal.Simulation$Study_estimate)
mean(Normal.Simulation$Study_sd)
hist(Normal.Simulation$Study_sd)

plot(Normal.Simulation$Study_estimate ~ Normal.Simulation$Study_n)
plot(Normal.Simulation$Study_sd ~ Normal.Simulation$Study_estimate)

testRes <- rma.uni(Study_estimate, Study_sd^2, data=Normal.Simulation, method="FE")
funnel(testRes)

### Checking Begg weights against sample size

dummy1 <- numeric()
dummy2 <- numeric()


for (a in 10:500){
  dummy3 <- numeric()
  for (b in 1:1000){
  Study_mu <- rnorm(1, mean = 0, sd = 1)
  Study_values <- rnorm(a, mean = Study_mu, sd = 2)
  Study_mean <- mean(Study_values)
  Study_StanDev <- sd(Study_values)
  
  Begg_weight <-exp(
    -Begg_b * (
      (Begg_sided * dnorm(-abs(Study_mean)/(Study_StanDev))) 
      ^Begg_a ) )
  dummy3 <- append(dummy3, Begg_weight)
  }
  
  dummy1 <- append(dummy1, a)
  dummy2 <- append(dummy2, mean(dummy3))
  
}

plot(dummy1, dummy2)

### Checking Preston weights against sample size

dummy1 <- numeric()
dummy2 <- numeric()


for (a in 10:500){
  dummy3 <- numeric()
  for (b in 1:1000){
    Study_mu <- rnorm(1, mean = 0, sd = 1)
    Study_values <- rnorm(a, mean = Study_mu, sd = 2)
    Study_mean <- mean(Study_values)
    Study_StanDev <- sd(Study_values)
    
    Begg_weight <-exp(
      - Begg_b * Study_StanDev/(sqrt(a)) * (
        (1- dnorm(-abs(Study_mean)/(Study_StanDev))) 
        ^2 ) )
    dummy3 <- append(dummy3, Begg_weight)
  }
  
  dummy1 <- append(dummy1, a)
  dummy2 <- append(dummy2, mean(dummy3))
  
}

plot(dummy1, dummy2)


#### Rewritten to plot for p-value and n

dummy1 <- numeric()
dummy2 <- numeric()
dummy3 <- numeric()

for (a in 10:1000){
  for(b in 1:10){
    Study_mu <- rnorm(1, mean = 0, sd = 1)
    Study_values <- rnorm(a, mean = Study_mu, sd = 1)
    Study_mean <- mean(Study_values)
    Study_StanDev <- sd(Study_values)
    
    Begg_weight <-exp(
      - 30 * Study_StanDev /sqrt(a)  * (
        (dnorm(-abs(Study_mean)/(Study_StanDev))) 
        ^1.5 ) )
    
  
  
  dummy1 <- append(dummy1, a)
  dummy2 <- append(dummy2, Begg_weight)
  dummy3 <- append(dummy3, dnorm(-abs(Study_mean)/(Study_StanDev)))
}
}

plot(dummy1, dummy2)
plot(dummy3, dummy2)

asdf <- data.frame(dummy1, dummy2,dummy3)
with(asdf[asdf$dummy1 < 100,], plot(dummy3, dummy2))
with(asdf[asdf$dummy3 > 0.1,], plot(dummy1, dummy2))


## How Begg weight changes with p-value separated by n

dummy1 <- numeric()
dummy2 <- numeric()
dummy3 <- numeric()
for (a in seq(0,1, 0.01) ){
  for (b in c(1000,600, 200, 100, 50 , 30)){
  #for (b in c(1, 0.6, 0.4, 0.3, 0.2, 0.1, 0.05, 0.03, 0.01)){  
  
  Begg_weight <-exp(
    - 3  * (b^(-0.3))* (
      (a)) 
      ^0.5 ) 
  
  
  dummy1 <- append(dummy1, a)
  dummy2 <- append(dummy2, Begg_weight)
  dummy3 <- append(dummy3, b)
  
  }
}

plot(dummy1, dummy2, ylim = c(0,1))
abline(v = 0.05)


### Checking integrals

Begg_weight_fn <- function(p){
  return(exp(-30*(0.01)*p^1.5))
}

### ?legitimacy of using bound from 0.05 - see Begg and Mazumdar
integrate(Begg_weight_fn, lower = 0.00, upper = 1)


