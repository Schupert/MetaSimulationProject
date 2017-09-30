#### Declare variables

# Reps = number of repetitions of experiment
Reps = 2

# k = number of studies in series
Studies = c(3,5)
#Studies = c(3,5,10,30)

# subj = number of subjects in study, likely to be distributed
Subj <- c(100,20)

# sd = study level standard deviation
True.sd = sqrt(2)

# theta = population level log(OR) - this should be considered more purely on the log scale
theta = c(log(0.25), log(0.75))

# tau.sq = between studies variance (can be squared due to sqrt() in normal draw), ?to be distributed
tau.sq = c(0, 0.01777778)

# Frequency of event averaged across 2 arms (before applying change due to theta) = EvFreq
EvFreq = c(0.1, 0.3)

# controlProp = proportion of total sample in control arm
controlProp = 0.5

## Boundary of step function on p value, causing severity of publication bias
Severity.boundary <- c(0.05, 0.2)

# Set up strength of publication bias selection IF STILL USING
Begg_a <- 1.5
Begg_b <- 4
Begg_sided <- 1

# Set up within study reporting bias - this is now one sided
Tested.outcomes <- 5
Chosen.outcomes <- 1
Sd.split <- 0.8

# Size of per unit bias increase
Bias.multiple <- 0.9

a <- numeric(length = 0)
b <- numeric(length = 0)
c <- numeric(length = 0)
d <- numeric(length = 0)
e <- numeric(length = 0)
f <- numeric(length = 0)
g <- numeric(length = 0)
counter1 <- numeric(length = 0)

for (i in Subj){
  for (k in theta){
    for (l in tau.sq){
      
      for (j in EvFreq){
        
        for (n in Studies){
          

          
          for (o in 1:n){
            
            for (m in 1:Reps){
              a <- append(a, i)
              b <- append(b,k)
              c <- append(c,l)
              d <- append(d,n)
              e <- append(e,j)
              f <- append(f,o)
              g <- append(g,m)
              
              counter <- as.integer((match(i, Subj)-1) * length(EvFreq) * length(theta) * length(tau.sq) * sum(Studies) * Reps + 
                                      (match(k, theta)-1) * length(tau.sq) * length(EvFreq) * sum(Studies) * Reps +
                                      (match(l, tau.sq)-1) * length(EvFreq) * sum(Studies) * Reps +
                                      (match(j, EvFreq)-1) * sum(Studies) * Reps +
                                      (sum(Studies[0:(match(n, Studies)-1)]) + o -1) * Reps +
                                      m
              )
              counter1 <- append(counter1, counter)
            }
          }
        }
      }
    }
  }
}

asdf1 <- data.table(counter1, a,b,c,d,e,f,g)
asdf1 <- asdf1[order(counter1)]

##### Need to re append values - specific to analysis
ID =  length(Subj) * length(theta) * length(tau.sq) * length(EvFreq) * Reps * sum(Studies)

Rep_Number =  rep(1:Reps, times = ID/Reps)
intermediate <- integer()
for (i in Studies){intermediate <- append(intermediate, rep(i, times = i*Reps))}
Rep_NumStudies = rep(intermediate, times = ID/(Reps*sum(Studies)))
rm(intermediate)
Rep_ev_freq = rep(rep(EvFreq, each = Reps * sum(Studies)), times = ID/(Reps*sum(Studies)*length(EvFreq)))
Rep_tau.sq = rep(rep(tau.sq, each = Reps * sum(Studies)*length(EvFreq)), times = ID/(Reps*sum(Studies)*length(tau.sq)*length(EvFreq)))
Rep_theta = rep( rep(theta, each = Reps * sum(Studies) * length(tau.sq)*length(EvFreq)), times = length(Subj))

### Create keyable vector for Subj
Subj2 <- c(100, 20)
Rep_Subj = rep(Subj2, each = ID / length(Subj))

identical(asdf1$a, Rep_Subj)
identical(asdf1$b, Rep_theta)
identical(asdf1$c, Rep_tau.sq)
identical(asdf1$e, Rep_ev_freq)
identical(asdf1$d, Rep_NumStudies)
identical(as.integer(asdf1$g), Rep_Number)
