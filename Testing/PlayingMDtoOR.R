x <- replicate(10000, {
  MD <- rnorm(1, 2, 2)
  cont <- rnorm(1000, 0 - 0.5*MD, 1)
  treatment <- rnorm(1000, 0 + 0.5*MD, 1)
  #cont <- rnorm(1000, 0, 1)
  #treatment <- rnorm(1000, MD, 1)
  cont.dic <- as.numeric(cont > 0)
  treat.dic <- as.numeric(treatment > 0)
  a <- sum(cont.dic)
  b <- length(cont.dic) - a
  c <- sum(treat.dic)
  d <- length(treat.dic) - c
  if (a == 0 | b == 0 | c == 0 | d == 0){
    a <- a + 0.5
    b <- b + 0.5
    c <- c + 0.5
    c <- c + 0.5
  }
 log((c/d)/(a/b))
  
})

hist(x, freq = FALSE, breaks = c(20))

cont <- rnorm(1000, 0, 1)
treatment <- rnorm(1000, 1, 1)
cont.dic <- as.numeric(cont > 0)
treat.dic <- as.numeric(treatment > 0)
a <- sum(cont.dic)
b <- length(cont.dic) - a
c <- sum(treat.dic)
d <- length(treat.dic) - c

Study_estimate = log((c/d)/(a/b))

mu = 0.5
individuals <- rnorm(100, 0, 0.5)
Pic <- exp(mu - 0.5*individuals) / (1 + exp(mu - 0.5*individuals))



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

dummy1 <- numeric()

for i in 1:1000{
  dummy2 <- rbinom(1, )
  
}


dummy2 <- rnorm(100,0,1)

dummy3 <- rbinom(100, 1, exp(dummy2)/(1 + exp(dummy2)))
sum(dummy3)



length(dummy3)
