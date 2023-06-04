###############################################################################

elnino <- read.csv("../data/elnino-small/elnino.csv", header = TRUE)

# ground truth about the underlying permutation
pistar <- unlist(read.csv("../data/elnino-small/Pistar_elnino.csv", header = FALSE))
pistarinv <- unlist(read.csv("../data/elnino-small/Pistarinv_elnino.csv", header = FALSE))
# sum(pistar != order(pistarinv)) -- check

# block index based on (latitude, longitude)-pairs 
blockix <- unlist(read.csv("../data/elnino-small/clusterindex_elnino.csv", header = FALSE))

# center data 
elnino_c <- elnino
sel <- c("zon.winds", "mer.winds", "humidity", "air.temp", "s.s.temp")
elnino_c[,sel] <- scale(elnino[,sel], scale = FALSE)

# response (matched) vs response (original)
plot(elnino_c$air.temp, elnino_c$air.temp[pistar])

# linear model with no linkage error
lm0 <- lm(air.temp ~ zon.winds + mer.winds + humidity + s.s.temp - 1, data = elnino_c)

# linear model with linkage error 
lmshuff <- lm(air.temp[pistar] ~ zon.winds + mer.winds + humidity + s.s.temp - 1, data = elnino_c)

###

source("../code/lahiri_larsen.R")
source("../code/optimal_matching.R")

X <- model.matrix(lmshuff)
y <- elnino_c$air.temp[pistar]

beta_Q <- coef(lahiri_larsen_block(X, y, blockix))
Q <- generate_Q_block(blockix) 
beta_Q_check <- lahiri_larsen(X, y, Q)

sqrt(sum((beta_Q - coef(lm0))^2))
sqrt(sum((coef(lmshuff) - coef(lm0))^2))

pihat <- optimal_matching(drop(X %*% beta_Q), y, blockix)

# response (after re-matching) vs. original response
plot(elnino_c$air.temp[pihat], elnino_c$air.temp[pistar])

###

source("../code/mixture_model.R")

tausq <- mean(y^2)
f0 <- function(z) dnorm(z, mean = 0, sd = sqrt(tausq))

res <- fit_mixture(X, y, f0, control = list(init = "robust"))

sqrt(sum((res$betahat - coef(lm0))^2))
