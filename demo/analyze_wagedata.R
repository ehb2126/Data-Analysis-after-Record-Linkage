###############################################################################

# 'ground truth': all variables with correct correspondence in a single file 
wages_f <- read.csv("../data/fakedata/wages_df.csv", header = TRUE)
lm0 <- lm(log(WAGE) ~ SEX + EXPERIENCE + I(EXPERIENCE^2)  + 
            EDUCATION + as.factor(OCCUPATION) + UNION,
          data = wages_f)
summary(lm0)

###############################################################################
# separate files (response in file A)

wages_dA <- read.csv("../data/fakedata/wages_dA.csv", header = TRUE)
wages_dB <- read.csv("../data/fakedata/wages_dB.csv", header = TRUE)

##############################################################################

# merge files and A based on probabilistic record linkage based on common variable "ZIPCODE"# (there are in fact more common variables) 
library(fastLink)

set.seed(1427)

wages_link<-fastLink(wages_dA, wages_dB, varnames = c("ZIPCODE"))

matched_wages <- getMatches(wages_dA, wages_dB, wages_link, combine.dfs = FALSE)

dA<-matched_wages$dfA.match; dB<-matched_wages$dfB.match


# to ensure unique column names in the merged file

commonvars <- intersect(colnames(dA), colnames(dB))
colnames(dA)[colnames(dA) %in% commonvars] <- paste("A.", colnames(dA)[colnames(dA) %in% commonvars], sep="")
colnames(dB)[colnames(dB) %in% commonvars] <- paste("B.", colnames(dB)[colnames(dB) %in% commonvars], sep="")


# merging the data mached by fastLink

merged_wages <-cbind.data.frame(dA, dB)

mean(merged_wages[,"A.ID"] != merged_wages[,"B.ID"]) # fraction of mismatches ~13%

plot(wages_f$WAGE[merged_wages[,"A.ID"]], wages_f$WAGE[merged_wages[,"B.ID"]]) # response before and after matching 

lm_merged <- lm(log(WAGE) ~ B.SEX + EXPERIENCE + I(EXPERIENCE^2)  + 
            EDUCATION + as.factor(OCCUPATION) + UNION,
          data = merged_wages)
summary(lm_merged)

### robust estimation

library(MASS)
rlm_merged <- rlm(log(WAGE) ~ B.SEX + EXPERIENCE + I(EXPERIENCE^2)  + 
            EDUCATION + as.factor(OCCUPATION) + UNION,
          data = merged_wages)

# comparing estimation error of naive and robust estimation 
sqrt(sum((coef(lm0) - coef(lm_merged))^2))
sqrt(sum((coef(lm0) - coef(rlm_merged))^2))

# mixture model
source("../code/mixture_model.R")
source("../code/mixture_model_sandwich.R")

X <- model.matrix(lm_merged)
X <- X[,!(colnames(X) %in% "(Intercept)")]
y <- model.extract(lm_merged$model, "response")
Xc <- apply(X, 2, function(z) z - mean(z))
yc <- y - mean(y)
tausq <- mean(yc^2)
f0 <- function(z) dnorm(z, mean = 0, sd = sqrt(tausq))

res <- fit_mixture(Xc, yc, f0, control = list(init = "robust"))
interc <- mean(y - X %*% res$betahat)
coef_mixture <- c(interc, res$betahat)

sqrt(sum((coef(lm0) - coef_mixture)^2))
sqrt(sum((coef(lm0) - coef(rlm_merged))^2))
sqrt(sum((coef(lm0) - coef(lm_merged))^2))

### blocking --- use information about linkage process
source("../code/lahiri_larsen.R")

blockix <- merged_wages$A.ZIPCODE 
barplot(table(as.numeric(table(blockix))))

beta_Q <- coef(lahiri_larsen_block(Xc, yc, blockix))
Q <- generate_Q_block(blockix) 
beta_Q_check <- lahiri_larsen(Xc, yc, Q)
    
interc <- mean(y - X %*% beta_Q)
coef_Q <- c(interc, beta_Q)

sqrt(sum((coef(lm0) - coef_Q)^2))

### re-matching based on sorting
source("../code/optimal_matching.R")

pihat <- optimal_matching(drop(Xc %*% res$betahat), yc, blockix)
pihatinv <- order(pihat)
pistar <- match(merged_wages[,"B.ID"], merged_wages[,"A.ID"])
plot(wages_f$WAGE[merged_wages[,"A.ID"]], wages_f$WAGE[merged_wages[,"B.ID"]])
points(wages_f$WAGE[merged_wages[,"A.ID"]], wages_f$WAGE[merged_wages[pihat,"B.ID"]], col = "red")

###
n <- nrow(merged_wages)
which(pihat != (1:n))
names(pihat) <- rep("", n)

### based on fit from mixture model

yhat <- Xc %*% res$betahat
resid <- yc - yhat

plot(yhat, resid, cex = 2.5, pch = 1)

ind <- merged_wages[,"A.ID"] != merged_wages[,"B.ID"]
points(yhat[ind], resid[ind], col = "red", cex = 2.5, pch = 16)
