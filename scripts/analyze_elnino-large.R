###############################################################################

elnino <- read.csv("../data/elnino-large/elnino.csv", header = TRUE)
elnino_merged <- read.csv("../data/elnino-large/elnino-merged.csv", header = TRUE)

plot(elnino$air.temp, elnino_merged$air.temp)

lm0 <- lm(air.temp ~ s.s.temp + zon.winds + mer.winds + humidity, data = elnino)
lmmerged <- lm(air.temp ~ s.s.temp + zon.winds + mer.winds + humidity, data = elnino_merged)

# this time block index is not provided. Here, we generate it from scratch --- quite time consuming (quadratic runtime) 
uniqlatlong <- unique(x = cbind(elnino$latitude, elnino$longitude))


blockix <- apply(elnino_merged, 1, function(z) which(z["latitude"] == uniqlatlong[,1] & z["longitude"] ==  uniqlatlong[,2])) #~takes about 30s
# CHECK: max(rowSums((uniqlatlong[blockix,] - elnino_merged[,c("latitude", "longitude")])^2))

source("../code/lahiri_larsen.R")
source("../code/optimal_matching.R")

X <- model.matrix(lmmerged)
X <- X[,!(colnames(X) %in% "(Intercept)")]
y <- model.extract(lmmerged$model, "response")
Xc <- apply(X, 2, function(z) z - mean(z))
yc <- y - mean(y)

beta_Q <- coef(lahiri_larsen_block(Xc, yc, blockix))
interc <- mean(y - X %*% beta_Q)

pihat <- optimal_matching(drop(X %*% beta_Q), y, blockix)

plot(elnino$air.temp[pihat], elnino_merged$air.temp)

mean((elnino$air.temp[pihat] - elnino_merged$air.temp)^2)
mean((elnino$air.temp - elnino_merged$air.temp)^2)

###

source("../code/mixture_model.R")

tausq <- mean(yc^2)
f0 <- function(z) dnorm(z, mean = 0, sd = sqrt(tausq))

res <- fit_mixture(Xc, yc, f0, control = list(init = "robust"))

interc <- mean(y - X %*% res$betahat)
