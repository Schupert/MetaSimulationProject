########## Calculating tau2 and I2

UMDI2 <- c(0,0.05, 0.5, 0.95)
LORI2 <- c(0,0.05,0.2, 0.95)

tau2UMD <- UMDI2 * 4 / 60 / (1 - UMDI2)
tau2LOR <- LORI2 * 16 / 100 / (1 - LORI2)

tau2UMD / (tau2UMD + 4 /20)
tau2UMD / (tau2UMD + 4 /60)
tau2UMD / (tau2UMD + 4 /100)
tau2UMD / (tau2UMD + 4 /250)
tau2UMD / (tau2UMD + 4 /1000)

tau2LOR / (tau2LOR + 16 /20)
tau2LOR / (tau2LOR + 16 /100)
tau2LOR / (tau2LOR + 16 /250)
tau2LOR / (tau2LOR + 16 /1000)

####### Calculating theta

LORtheta <- c(log(1/4), log(0.8), log(1), log(1.25), log(4))
SMDtheta <- LORtheta * ( - sqrt(3) / pi)

#### Basing UMD tau2 on LOR

newtau2UMD <- tau2LOR * ( sqrt(3) / pi)

newtau2UMD / (newtau2UMD + 4 /20)
newtau2UMD / (newtau2UMD + 4 /60)
newtau2UMD / (newtau2UMD + 4 /100)
newtau2UMD / (newtau2UMD + 4 /250)
newtau2UMD / (newtau2UMD + 4 /1000)