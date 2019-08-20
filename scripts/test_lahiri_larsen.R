###############################################################################

set.seed(1019)
n <- 5000
blockindex <- sample(1:500, n, replace = TRUE)
d <- 5
X <- matrix(rnorm(n*d), nrow = n, ncol = d)
betastar <- c(1,1,-1,-1,0)

# generate a blockwise random permutation
Lstarlist <- lapply(split(1:n, blockindex), function(z) cbind(z, z * (length(z)==1) + (length(z) > 1) * sample(x=z, size=length(z))))
Lstar <- cbind(unlist(lapply(Lstarlist, function(z) z[,1])), unlist(lapply(Lstarlist, function(z) z[,2]))) 
pistar <- Lstar[order(Lstar[,1]), 2]

Xpistar <- X[pistar,]
sigma <- 0.2
xi <- sigma * rnorm(n)
y <- Xpistar %*% betastar  + xi

###############################################################################

source("../code/optimal_matching.R")
source("../code/lahiri_larsen.R")

# approach I (allows for general Q, but if Q follows block structure this is slow)
Q <- generate_Q_block(blockindex)
betaQ <- lahiri_larsen(X, y, Q)
coef(lm(y ~ Xpistar - 1))

# approach II (tailored to block structure, faster)
betaQcheck <- coef(lahiri_larsen_block(X, y, blockindex))

###############################################################################

pihat <- optimal_matching(X %*% betaQ, y, blockindex) 
pihat <- optimal_matching(X %*% betaQ, y, blockindex=rep(1,n))

# before and after re-matching

mean(pistar != (1:n)) # 90%mismatch

# before
plot(Xpistar %*% betastar, X %*% betastar)
# after
plot(Xpistar %*% betastar, X[pihat,] %*% betastar)
