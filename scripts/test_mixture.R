### read simulate data with 10% mismatches 

X <- as.matrix(read.csv("../data/simulated/Xsim.csv", header = FALSE))
y <- as.matrix(read.csv("../data/simulated/ysim.csv", header = FALSE))

source("../code/mixture_model.R")

res0 <- fit_mixture0(X, y, control = list(init = "robust"))
tausq <- mean(y^2)

f0 <- function(z) dnorm(z, mean = 0, sd = sqrt(tausq))
res <- fit_mixture(X, y, f0, control = list(init = "robust"))

### analysis of large merged el nino data set

merged_elnino <- read.csv("../data/elnino-large/elnino-merged.csv", header = TRUE)

# work with centered data (simulated data were generated to be centered)
merged_elnino_c <- apply(merged_elnino, 2, function(z) z -mean(z))
tausq <- mean(merged_elnino_c[,"air.temp"]^2)
f0 <- function(z) dnorm(z, mean = 0, sd = sqrt(tausq))

# fully parametric model
res <- fit_mixture(merged_elnino_c[,c("zon.winds", "mer.winds", "humidity", "s.s.temp")], as.numeric(merged_elnino_c[,"air.temp"]), f0, control = list(init = "robust"))
interc <- mean(merged_elnino[,"air.temp"] - as.matrix(merged_elnino[,c("zon.winds", "mer.winds", "humidity", "s.s.temp")]) %*% res$betahat)


# partially non-parametric model 
dens <- density(merged_elnino_c[,"air.temp"])
f0dens <- approxfun(dens$x, dens$y)

resdens <- fit_mixture(merged_elnino_c[,c("zon.winds", "mer.winds", "humidity", "s.s.temp")], as.numeric(merged_elnino_c[,"air.temp"]), f0dens, control = list(init = "robust"))




