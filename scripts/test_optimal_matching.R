###############################################################################

set.seed(1019)
n <- 1000
blockindex <- sample(1:200, n, replace = TRUE) # random division into blocks
x <- rnorm(n)

# generate random permutations within each block
Lstarlist <- lapply(split(1:n, blockindex), function(z) cbind(z, z * (length(z)==1) + (length(z) > 1) * sample(x=z, size=length(z))))
Lstar <- cbind(unlist(lapply(Lstarlist, function(z) z[,1])), unlist(lapply(Lstarlist, function(z) z[,2]))) 
pistar <- Lstar[order(Lstar[,1]), 2]

xpistar <- x[pistar]
sigma <- 0.001
xi <- sigma * rnorm(n)
y <- 2*xpistar + xi

###############################################################################

source("../code/optimal_matching.R")

pihat <- optimal_matching(x, y, blockindex) 
pihat <- optimal_matching(x, y, blockindex=rep(1,n))
sum(pihat != pistar)
