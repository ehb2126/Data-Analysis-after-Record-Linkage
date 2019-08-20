###############################################################################

# simulate some data
set.seed(1019)
n <- 1000
blockindex <- sample(1:200, n, replace = TRUE)
d <- 5
X <- matrix(rnorm(n*d), nrow = n, ncol = d)
betastar <- c(1,1,-1,-1,0)

# randomly 
alpha <- 0.15
pistar <- c(sample(n*0.15), (n*0.15 + 1):n)

Xpistar <- X[pistar,]
sigma <- 0.2
xi <- sigma * rnorm(n)
y <- Xpistar %*% betastar  + xi

###############################################################################

source("../code/robust_regression.R")
library(MASS)

rlm0<- rlm(X, y)

robreg <- robust_regression(X, y, lambda = 1.345 * rlm0$s)

###############################################################################

plot(robreg$xi, type = "h")
abline(v = n*alpha, col = "red")

sqrt(sum((robreg$beta - betastar)^2)) # estimation error of robust estimator
sqrt(sum((coef(lm(y ~ X - 1)) - betastar)^2)) # estimation erro of naive estimator

###############################################################################


