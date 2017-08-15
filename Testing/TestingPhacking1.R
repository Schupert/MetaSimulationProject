#### Testing p-hacking

### Counfounder
asdf <- rnorm(100, 100, 10)
conf1 <- rnorm(100, 0, 1)
model1 <- lm(asdf~ conf1)
summary(model1)

### Random sample

cat1 <- sample(c(1,2,3), 100, replace= TRUE)
asdf[cat1 == TRUE]
