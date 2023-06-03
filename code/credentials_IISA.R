#######################################################################
###  Washington HealthCare Credential Data revisited               ####
###                  --------------------------------------------######
#######################################################################
ls()

### [1] READ-in generously linked data set

df_mm <- read.csv("credentials_merged.csv", sep = ",", header = TRUE)
colnames(df_mm)[1] <- "idx"

# variable of interest:
summary(df_mm$temp_time)
mean(df_mm$temp_time < 0) # at least 7.5% are mismatches

# filter
df_mm_cl <- df_mm[-which(df_mm$temp_time < 0),]

# generate x,y
xperm <- (df_mm_cl$temp_issue_integer - min(df_mm_cl$temp_issue_integer)) / (max(df_mm_cl$temp_issue_integer))

yperm <- df_mm_cl$temp_time
nperm <- length(yperm) # 72,226
ylimperm <- c(0,1000)/2

pdf("scatter_data_perm.pdf")
par(mar = c(4, 5.3, 1, 3.5))
plot(xperm, yperm, pch = 16, cex = 2/3, ylim = ylimperm, xlab = "time", cex.axis = 1.9, cex.lab = 2.3,
     ylab = "Duration (days)")

mtext("01/01/2009", side = 1, cex = 2, col = "blue", line = 3, at = 0)
mtext("12/31/2021", side = 1, cex = 2, col = "blue", line = 3, at = 1)
dev.off()

### [2] READ-in restrictively linked data set

df_em <- read.table("credentials_merged_exactmatches.csv", sep = ",", header = TRUE,
                    quote="")

# still roughly 1.4% times are problematic (negative or more than 1.5 years of waiting)

# further processed (filtered) data
df_em_cl <- df_em[-which(df_em$temp_time < 0),]


xoracle  <- (df_em_cl$temp_issue_integer - min(df_em_cl$temp_issue_integer)) / (max(df_em_cl$temp_issue_integer))

yoracle <- df_em_cl$temp_time
noracle <- length(yoracle)

pdf("scatter_data_oracle.pdf")
par(mar = c(4, 5.3, 1, 3.5))
plot(xoracle, yoracle, pch = 16, cex = 2/3, ylim = ylimperm, xlab = "time", cex.axis = 1.9, cex.lab = 2.3,
     ylab = "Duration (days)")

mtext("01/01/2009", side = 1, cex = 2, col = "blue", line = 3, at = 0)
mtext("12/31/2021", side = 1, cex = 2, col = "blue", line = 3, at = 1)
dev.off()

### [---] knot sequence and grid for plotting curves

knots <- seq(from = 0, to = 1, length = 1000)
xgr <-  seq(from = 0, to = 0.99, by = 0.01)
xgr_df <- data.frame(x = xgr)

### [1.1] perform plain spline smoothing for generously linked data set

library(mgcv)

y <- yperm
x <- xperm

gamperm <- gam(y ~ s(x, fx=FALSE, k = 32, bs="cr"), knots = data.frame(x0 = knots))

pred_perm <- predict(gamperm, newdata = xgr_df, se.fit = TRUE)

                                        #, col = "blue", lwd = 3, type = "l")


# produce plot w/ and w/o pointwise confidence bands

pdf("fit_naive.pdf")
par(mar = c(4, 5.3, 1, 3.5))
plot(xgr, pred_perm$fit, col = "darkgreen", lwd = 3, type = "l", xlab = "time",
     ylab = "Duration (days)", cex.axis = 1.9, cex.lab = 2.3)
mtext("01/01/2009", side = 1, cex = 2, col = "blue", line = 3, at = 0)
mtext("12/31/2021", side = 1, cex = 2, col = "blue", line = 3, at = 1)
dev.off()

pdf("fit_naive_CI.pdf")

par(mar = c(4, 5.3, 1, 3.5))
plot(xgr, pred_perm$fit, col = "darkgreen", lwd = 3, type = "l", xlab = "time",
     ylab = "Duration (days)", cex.axis = 1.9, cex.lab = 2.3)

polygon(c(xgr,rev(xgr)), c(pred_perm$fit - 2 * pred_perm$se.fit, rev(pred_perm$fit + 2 * pred_perm$se.fit)), col = "grey", border = "grey")

lines(xgr, pred_perm$fit, col = "darkgreen", lwd = 3, type = "l")
dev.off()

# extract some relevant info from the GAM fit

gamobj_perm <- gam(y ~ s(x, fx=FALSE, k = 32, bs="cr"), knots = data.frame(x0 = knots), fit = FALSE)
Sgam_perm <- rbind(0, cbind(0, gamobj_perm$smooth[[1]]$S[[1]]))
Bgam_perm <- as.matrix(gamobj_perm$X)
lambda_perm <- gamperm$sp

### [1.2] attempt adjustment for generously linked data set

sigma <- sqrt(mean(residuals(gamperm)^2))
mucur <- fitted(gamperm)
fy <- dnorm(y, mean = mean(y), sd = sd(y))
fymu <- function(mu, sigma) dnorm(y, mu, sigma)

# preliminary checks

nloglik_marginal <- function(beta, sigma, alpha, sigma0) sum(-log((1-alpha) * (fymu(drop(Bgam_perm %*% beta), sigma) ) + alpha * fy)) + 0.5*(lambda_perm/sigma0^2)*sum(beta * (Sgam_perm %*% beta))

betacur <- coef(gamperm)
mucur <- fitted(gamperm)
alpha <- 0.2


sigma0 <- sqrt(mean(residuals(gamperm)^2))
sigma <- 0.95*sigma0

tol <- 1E-4
maxiter <- 1000
objs <- numeric(maxiter)
maxiter <- 1000
iter <- 1
objs[iter] <- nloglik_marginal(betacur, sigma, alpha, sigma0)


while(iter < maxiter){

    num <- (1-alpha) * fymu(mucur, sigma)
    denom <- num + alpha * fy
    pcur <- num/denom
    alpha <- 1 - mean(pcur)

     wglmfit <- gam(y ~ s(x, fx=FALSE, k = 32, bs="cr"), knots = data.frame(x0 = knots), weights = pcur, sp = (sigma/sigma0)^2 * lambda_perm)


    sigma <- sqrt(sum(residuals(wglmfit)^2) / sum(pcur))
    betacur <- coef(wglmfit)

    mucur <- Bgam_perm %*% betacur
    iter <- iter + 1
    objs[iter] <- nloglik_marginal(betacur, sigma, alpha, sigma0)
    if(objs[iter] + tol > objs[iter-1])
        break
}

wglmfit_naive_unrestricted <- wglmfit

#    objs[iter]

#}

# without any restriction on (sigma / sigma0), that ratio is 0.08,
# and alpha is .161 (which is larger than the suspected mismatch rate)
# the objective at the end equals 420676.5

pdf("fit_naive_sigma_unrestricted.pdf")
par(mar = c(4, 5.3, 1, 3.5))
plot(xgr, predict(gamperm, xgr_df), col = "darkgreen", lwd = 3, type = "l", ylim = ylimperm * 0.8,
     xlab = "time", ylab = "Duration (days)", cex.axis = 1.9, cex.lab = 2.)
lines(xgr, predict(wglmfit, xgr_df), col = "orange", lwd = 3, type = "l")
mtext("01/01/2009", side = 1, cex = 2, col = "blue", line = 3, at = 0)
mtext("12/31/2021", side = 1, cex = 2, col = "blue", line = 3, at = 1)
dev.off()

# now repeat the plot with pointwise confidence bands

pred_wglmfit <- predict(wglmfit, newdata = xgr_df, se.fit = TRUE)

pdf("fit_naive_sigma_unrestricted_CI.pdf")
par(mar = c(4, 5.3, 1, 3.5))
plot(xgr, pred_perm$fit, col = "darkgreen", lwd = 3, type = "l", xlab = "time",
     ylab = "Duration (days)", cex.axis = 1.9, cex.lab = 2.3)

polygon(c(xgr,rev(xgr)), c(pred_perm$fit - 2 * pred_perm$se.fit, rev(pred_perm$fit + 2 * pred_perm$se.fit)), col = "grey", border = "grey")

lines(xgr, pred_perm$fit, col = "darkgreen", lwd = 3, type = "l")

polygon(c(xgr,rev(xgr)), c(pred_wglmfit$fit - 2 * pred_wglmfit$se.fit, rev(pred_wglmfit$fit + 2 * pred_wglmfit$se.fit)), col = "grey", border = "grey")
lines(xgr, predict(wglmfit, xgr_df), col = "orange", lwd = 3, type = "l")
mtext("01/01/2009", side = 1, cex = 2, col = "blue", line = 3, at = 0)
mtext("12/31/2021", side = 1, cex = 2, col = "blue", line = 3, at = 1)
dev.off()


#  grid for (sigma / sigma0) ratio
sigmagr <- c(0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1)

final_obj <- numeric(length(sigmagr))
final_alpha <- numeric(length(sigmagr))


for(j in 1:length(sigmagr)){

betacur <- coef(gamperm)
mucur <- fitted(gamperm)
alpha <- 0.2


sigma0 <- sqrt(mean(residuals(gamperm)^2))
sigma <- sigma0 * sigmagr[j]

tol <- 1E-4
maxiter <- 1000
objs <- numeric(maxiter)
maxiter <- 1000
iter <- 1
objs[iter] <- nloglik_marginal(betacur, sigma, alpha, sigma0)


while(iter < maxiter){

    num <- (1-alpha) * fymu(mucur, sigma)
    denom <- num + alpha * fy
    pcur <- num/denom
    alpha <- 1 - mean(pcur)

    wglmfit <- gam(y ~ s(x, fx=FALSE, k = 32, bs="cr"), knots = data.frame(x0 = knots), weights = pcur, sp = (sigma/sigma0)^2 * lambda_perm)


    #sigma <- sqrt(sum(residuals(wglmfit)^2) / sum(pcur))
    betacur <- coef(wglmfit)

    mucur <- Bgam_perm %*% betacur
    iter <- iter + 1
    objs[iter] <- nloglik_marginal(betacur, sigma, alpha, sigma0)
    if(objs[iter] + tol > objs[iter-1])
        break
}

    final_obj[j] <- objs[iter]
    final_alpha[j] <- alpha

}

### after checking all fitted models, we have a winner:

jstar <- which(sigmagr == 0.4)

# refit that specific model and store the result


betacur <- coef(gamperm)
mucur <- fitted(gamperm)
alpha <- 0.2


sigma0 <- sqrt(mean(residuals(gamperm)^2))
sigma <- sigma0 * sigmagr[jstar]

tol <- 1E-4
maxiter <- 1000
objs <- numeric(maxiter)
maxiter <- 1000
iter <- 1
objs[iter] <- nloglik_marginal(betacur, sigma, alpha, sigma0)


while(iter < maxiter){

    num <- (1-alpha) * fymu(mucur, sigma)
    denom <- num + alpha * fy
    pcur <- num/denom
    alpha <- 1 - mean(pcur)

     wglmfit <- gam(y ~ s(x, fx=FALSE, k = 32, bs="cr"), knots = data.frame(x0 = knots), weights = pcur, sp = (sigma/sigma0)^2 * lambda_perm)


    betacur <- coef(wglmfit)

    mucur <- Bgam_perm %*% betacur
    iter <- iter + 1
    objs[iter] <- nloglik_marginal(betacur, sigma, alpha, sigma0)
    if(objs[iter] + tol > objs[iter-1])
        break
}

wglmfit_winner <- wglmfit
alphahat_winner <- alpha
sigmahat_winner <- sigma
betacur_winner <- betacur


pdf("fit_naive_sigma_restricted_winner.pdf")
par(mar = c(4, 5.3, 1, 3.5))
plot(xgr, predict(gamperm, xgr_df), col = "darkgreen", lwd = 3, type = "l", ylim = ylimperm * 0.8,
     xlab = "time", ylab = "Duration (days)", cex.axis = 1.9, cex.lab = 2.)
lines(xgr, predict(wglmfit_winner, xgr_df), col = "orange", lwd = 3, type = "l")
mtext("01/01/2009", side = 1, cex = 2, col = "blue", line = 3, at = 0)
mtext("12/31/2021", side = 1, cex = 2, col = "blue", line = 3, at = 1)
dev.off()

# now repeat the plot with pointwise confidence bands

pred_wglmfit_winner <- predict(wglmfit_winner, newdata = xgr_df, se.fit = TRUE)


pdf("fit_naive_sigma_restricted_winner_CI.pdf")
par(mar = c(4, 5.3, 1, 3.5))
plot(xgr, pred_perm$fit, col = "darkgreen", lwd = 3, type = "l", xlab = "time",
     ylab = "Duration (days)", cex.axis = 1.9, cex.lab = 2.3)

polygon(c(xgr,rev(xgr)), c(pred_perm$fit - 2 * pred_perm$se.fit, rev(pred_perm$fit + 2 * pred_perm$se.fit)), col = "grey", border = "grey")

lines(xgr, pred_perm$fit, col = "darkgreen", lwd = 3, type = "l")

polygon(c(xgr,rev(xgr)), c(pred_wglmfit_winner$fit - 2 * pred_wglmfit_winner$se.fit, rev(pred_wglmfit_winner$fit + 2 * pred_wglmfit_winner$se.fit)), col = "grey", border = "grey")
lines(xgr, predict(wglmfit_winner, xgr_df), col = "orange", lwd = 3, type = "l")
mtext("01/01/2009", side = 1, cex = 2, col = "blue", line = 3, at = 0)
mtext("12/31/2021", side = 1, cex = 2, col = "blue", line = 3, at = 1)
dev.off()

save.image("credentials_revisited_generously_linked.Rdata")

###*********************************************************************###
###*********************************************************************###
### >>>>  First part of analysis is complete                        <<<<###
###*********************************************************************###
###*********************************************************************###

### [1.3] perform plain spline smoothing for generously linked data set

y <- yoracle
x <- xoracle

gamoracle <- gam(y ~ s(x, fx=FALSE, k = 32, bs="cr"), knots = data.frame(x0 = knots))

pred_oracle <- predict(gamoracle, newdata = xgr_df, se.fit = TRUE)

# just oracle

pdf("fit_oracle.pdf")
par(mar = c(4, 5.3, 1, 3.5))
plot(xgr, pred_oracle$fit, col = "red", lwd = 3, type = "l", xlab = "time",
     ylab = "Duration (days)", cex.axis = 1.9, cex.lab = 2.3, ylim = c(0,375))
mtext("01/01/2009", side = 1, cex = 2, col = "blue", line = 3, at = 0)
mtext("12/31/2021", side = 1, cex = 2, col = "blue", line = 3, at = 1)
dev.off()

# oracle vs. naive

pdf("fit_naive_vs_oracle.pdf")
par(mar = c(4, 5.3, 1, 3.5))
plot(xgr, pred_perm$fit, col = "darkgreen", lwd = 3, type = "l", xlab = "time",
     ylab = "Duration (days)", cex.axis = 1.9, cex.lab = 2.3)
lines(xgr, pred_oracle$fit, col = "red", lwd = 3, type = "l", xlab = "time",
     ylab = "Duration (days)", cex.axis = 1.9, cex.lab = 2.3)
mtext("01/01/2009", side = 1, cex = 2, col = "blue", line = 3, at = 0)
mtext("12/31/2021", side = 1, cex = 2, col = "blue", line = 3, at = 1)
dev.off()

# oracle with CI

pdf("fit_oracle_CI.pdf")

par(mar = c(4, 5.3, 1, 3.5))
plot(xgr, pred_oracle$fit, col = "red", lwd = 3, type = "l", xlab = "time",
     ylab = "Duration (days)", cex.axis = 1.9, cex.lab = 2.3, ylim = c(0,375))

polygon(c(xgr,rev(xgr)), c(pred_oracle$fit - 2 * pred_oracle$se.fit, rev(pred_oracle$fit + 2 * pred_oracle$se.fit)), col = "grey", border = "grey")

lines(xgr, pred_oracle$fit, col = "red", lwd = 3, type = "l")
dev.off()

# oracle vs. naive, both with CI

pdf("fit_naive_vs_oracle_CI.pdf")

par(mar = c(4, 5.3, 1, 3.5))

plot(xgr, pred_perm$fit, col = "darkgreen", lwd = 3, type = "l", xlab = "time",
     ylab = "Duration (days)", cex.axis = 1.9, cex.lab = 2.3)

polygon(c(xgr,rev(xgr)), c(pred_perm$fit - 2 * pred_perm$se.fit, rev(pred_perm$fit + 2 * pred_perm$se.fit)), col = "grey", border = "grey")

lines(xgr, pred_perm$fit, col = "darkgreen", lwd = 3, type = "l")


polygon(c(xgr,rev(xgr)), c(pred_oracle$fit - 2 * pred_oracle$se.fit, rev(pred_oracle$fit + 2 * pred_oracle$se.fit)), col = "grey", border = "grey")

lines(xgr, pred_oracle$fit, col = "red", lwd = 3, type = "l")
dev.off()

# extract some relevant info from the GAM fit

gamobj_oracle <- gam(y ~ s(x, fx=FALSE, k = 32, bs="cr"), knots = data.frame(x0 = knots), fit = FALSE)
Sgam_oracle <- rbind(0, cbind(0, gamobj_oracle$smooth[[1]]$S[[1]]))
Bgam_oracle <- as.matrix(gamobj_oracle$X)
lambda_oracle <- gamoracle$sp

### [1.4] attempt adjustment for restrictively linked data set


sigma <- sqrt(mean(residuals(gamoracle)^2))
mucur <- fitted(gamoracle)
fy <- dnorm(y, mean = mean(y), sd = sd(y))
fymu <- function(mu, sigma) dnorm(y, mu, sigma)

# preliminary checks


nloglik_marginal <- function(beta, sigma, alpha, sigma0) sum(-log((1-alpha) * (fymu(drop(Bgam_oracle %*% beta), sigma) ) + alpha * fy)) + 0.5*(lambda_oracle/sigma0^2)*sum(beta * (Sgam_oracle %*% beta))

betacur <- coef(gamoracle)
mucur <- fitted(gamoracle)
alpha <- 0.2



sigma0 <- sqrt(mean(residuals(gamoracle)^2))
sigma <- 0.95*sigma0

tol <- 1E-4
maxiter <- 1000
objs <- numeric(maxiter)
maxiter <- 1000
iter <- 1
objs[iter] <- nloglik_marginal(betacur, sigma, alpha, sigma0)


while(iter < maxiter){

    num <- (1-alpha) * fymu(mucur, sigma)
    denom <- num + alpha * fy
    pcur <- num/denom
    alpha <- 1 - mean(pcur)

     wglmfit <- gam(y ~ s(x, fx=FALSE, k = 32, bs="cr"), knots = data.frame(x0 = knots), weights = pcur, sp = (sigma/sigma0)^2 * lambda_oracle)


    sigma <- sqrt(sum(residuals(wglmfit)^2) / sum(pcur))
    betacur <- coef(wglmfit)

    mucur <- Bgam_oracle %*% betacur
    iter <- iter + 1
    objs[iter] <- nloglik_marginal(betacur, sigma, alpha, sigma0)
    if(objs[iter] + tol > objs[iter-1])
        break
}

wglmfit_oracle_unrestricted <- wglmfit

# result when sigma is not further restricted

# results before and after adjustment
pdf("fit_oracle_sigma_unrestricted.pdf")
par(mar = c(4, 5.3, 1, 3.5))
plot(xgr, predict(gamoracle, xgr_df), col = "red", lwd = 3, type = "l", ylim = ylimperm * 0.8,
     xlab = "time", ylab = "Duration (days)", cex.axis = 1.9, cex.lab = 2.)
lines(xgr, predict(wglmfit, xgr_df), col = "purple", lwd = 3, type = "l") # wglmfit -> wglmfit_oracle_unrestricted
mtext("01/01/2009", side = 1, cex = 2, col = "blue", line = 3, at = 0)
mtext("12/31/2021", side = 1, cex = 2, col = "blue", line = 3, at = 1)
dev.off()

# results before and after adjustment (w/ confidence intervals)

pred_wglmfit <- predict(wglmfit, newdata = xgr_df, se.fit = TRUE)

#

pdf("fit_oracle_sigma_unrestricted_CI.pdf")
par(mar = c(4, 5.3, 1, 3.5))
plot(xgr, pred_oracle$fit, col = "red", lwd = 3, type = "l", xlab = "time",
     ylab = "Duration (days)", cex.axis = 1.9, cex.lab = 2.3, ylim = c(0, 375))

polygon(c(xgr,rev(xgr)), c(pred_oracle$fit - 2 * pred_oracle$se.fit, rev(pred_oracle$fit + 2 * pred_oracle$se.fit)), col = "grey", border = "grey")

lines(xgr, pred_oracle$fit, col = "red", lwd = 3, type = "l")

polygon(c(xgr,rev(xgr)), c(pred_wglmfit$fit - 2 * pred_wglmfit$se.fit, rev(pred_wglmfit$fit + 2 * pred_wglmfit$se.fit)), col = "grey", border = "grey")
lines(xgr, predict(wglmfit, xgr_df), col = "purple", lwd = 3, type = "l")
mtext("01/01/2009", side = 1, cex = 2, col = "blue", line = 3, at = 0)
mtext("12/31/2021", side = 1, cex = 2, col = "blue", line = 3, at = 1)
dev.off()

#

pdf("fit_oracle_vs_naive_sigma_unrestricted.pdf")
par(mar = c(4, 5.3, 1, 3.5))
plot(xgr, predict(gamperm, xgr_df), col = "darkgreen", lwd = 3, type = "l",
     xlab = "time", ylab = "Duration (days)", cex.axis = 1.9, cex.lab = 2.3)
lines(xgr, predict(gamoracle, xgr_df), col = "red", lwd = 3, type = "l",
      xlab = "time", ylab = "Duration (days)", cex.axis = 1.9, cex.lab = 2.3)
lines(xgr, predict(wglmfit_naive_unrestricted, xgr_df), col = "orange", lwd = 3, type = "l")
lines(xgr, predict(wglmfit_oracle_unrestricted, xgr_df), col = "purple", lwd = 3, type = "l")
mtext("01/01/2009", side = 1, cex = 2, col = "blue", line = 3, at = 0)
mtext("12/31/2021", side = 1, cex = 2, col = "blue", line = 3, at = 1)
dev.off()

# now re-create this plot with CIs

pred_wglmfit_naive_unrestricted <- predict(wglmfit_naive_unrestricted, newdata = xgr_df, se.fit = TRUE)
pred_wglmfit_oracle_unrestricted <- predict(wglmfit_oracle_unrestricted, newdata = xgr_df, se.fit = TRUE)

pdf("fit_oracle_vs_naive_sigma_unrestricted_CI.pdf")
par(mar = c(4, 5.3, 1, 3.5))
plot(xgr, predict(gamperm, xgr_df), col = "darkgreen", lwd = 3, type = "l",
     xlab = "time", ylab = "Duration (days)", cex.axis = 1.9, cex.lab = 2.3)

polygon(c(xgr,rev(xgr)), c(pred_oracle$fit - 2 * pred_oracle$se.fit, rev(pred_oracle$fit + 2 * pred_oracle$se.fit)), col = "grey", border = "grey")

polygon(c(xgr,rev(xgr)), c(pred_perm$fit - 2 * pred_perm$se.fit, rev(pred_perm$fit + 2 * pred_perm$se.fit)), col = "grey", border = "grey")

polygon(c(xgr,rev(xgr)), c(pred_wglmfit_naive_unrestricted$fit - 2 * pred_wglmfit_naive_unrestricted$se.fit, rev(pred_wglmfit_naive_unrestricted$fit + 2 * pred_wglmfit_naive_unrestricted$se.fit)), col = "grey", border = "grey")


polygon(c(xgr,rev(xgr)), c(pred_wglmfit_oracle_unrestricted$fit - 2 * pred_wglmfit_oracle_unrestricted$se.fit, rev(pred_wglmfit_oracle_unrestricted$fit + 2 * pred_wglmfit_oracle_unrestricted$se.fit)), col = "grey", border = "grey")


lines(xgr, predict(gamperm, xgr_df), col = "darkgreen", lwd = 3, type = "l",
      xlab = "time", ylab = "Duration (days)", cex.axis = 1.9, cex.lab = 2.)
lines(xgr, predict(gamoracle, xgr_df), col = "red", lwd = 3, type = "l",
      xlab = "time", ylab = "Duration (days)", cex.axis = 1.9, cex.lab = 2.)
lines(xgr, predict(wglmfit_naive_unrestricted, xgr_df), col = "orange", lwd = 3, type = "l")
lines(xgr, predict(wglmfit_oracle_unrestricted, xgr_df), col = "purple", lwd = 3, type = "l")
mtext("01/01/2009", side = 1, cex = 2, col = "blue", line = 3, at = 0)
mtext("12/31/2021", side = 1, cex = 2, col = "blue", line = 3, at = 1)

dev.off()



# now sigmagrid adjustment
sigmagr <- c(0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1)

final_obj <- numeric(length(sigmagr))
final_alpha <- numeric(length(sigmagr))


for(j in 1:length(sigmagr)){

betacur <- coef(gamoracle)
mucur <- fitted(gamoracle)
alpha <- 0.2


sigma0 <- sqrt(mean(residuals(gamoracle)^2))
sigma <- sigma0 * sigmagr[j]

tol <- 1E-4
maxiter <- 1000
objs <- numeric(maxiter)
maxiter <- 1000
iter <- 1
objs[iter] <- nloglik_marginal(betacur, sigma, alpha, sigma0)


while(iter < maxiter){

    num <- (1-alpha) * fymu(mucur, sigma)
    denom <- num + alpha * fy
    pcur <- num/denom
    alpha <- 1 - mean(pcur)

    wglmfit <- gam(y ~ s(x, fx=FALSE, k = 32, bs="cr"), knots = data.frame(x0 = knots), weights = pcur, sp = (sigma/sigma0)^2 * lambda_oracle)

    betacur <- coef(wglmfit)

    mucur <- Bgam_oracle %*% betacur
    iter <- iter + 1
    objs[iter] <- nloglik_marginal(betacur, sigma, alpha, sigma0)
    if(objs[iter] + tol > objs[iter-1])
        break
}

    final_obj[j] <- objs[iter]
    final_alpha[j] <- alpha

}


### after checking all fitted models, we have a winner:

jstar_oracle <- which(sigmagr == 0.7)

betacur <- coef(gamoracle)
mucur <- fitted(gamoracle)
alpha <- 0.2


sigma0 <- sqrt(mean(residuals(gamoracle)^2))
sigma <- sigma0 * sigmagr[jstar_oracle]

tol <- 1E-4
maxiter <- 1000
objs <- numeric(maxiter)
maxiter <- 1000
iter <- 1
objs[iter] <- nloglik_marginal(betacur, sigma, alpha, sigma0)


while(iter < maxiter){

    num <- (1-alpha) * fymu(mucur, sigma)
    denom <- num + alpha * fy
    pcur <- num/denom
    alpha <- 1 - mean(pcur)

   wglmfit <- gam(y ~ s(x, fx=FALSE, k = 32, bs="cr"), knots = data.frame(x0 = knots), weights = pcur, sp = (sigma/sigma0)^2 * lambda_oracle)

    betacur <- coef(wglmfit)

    mucur <- Bgam_oracle %*% betacur
    iter <- iter + 1
    objs[iter] <- nloglik_marginal(betacur, sigma, alpha, sigma0)
    if(objs[iter] + tol > objs[iter-1])
        break
}

wglmfit_winner_oracle  <- wglmfit
alphahat_winner_oracle <- alpha
sigmahat_winner_oracle <- sigma
betacur_winner_oracle  <- betacur

###

save.image("credentials_revisited_restrictively_linked.Rdata")

###

pdf("fit_oracle_sigma_restricted_winner.pdf")
par(mar = c(4, 5.3, 1, 3.5))
plot(xgr, predict(gamoracle, xgr_df), col = "red", lwd = 3, type = "l", ylim = ylimperm * 0.8,
     xlab = "time", ylab = "Duration (days)", cex.axis = 1.9, cex.lab = 2.)
lines(xgr, predict(wglmfit_winner_oracle, xgr_df), col = "purple", lwd = 3, type = "l")
mtext("01/01/2009", side = 1, cex = 2, col = "blue", line = 3, at = 0)
mtext("12/31/2021", side = 1, cex = 2, col = "blue", line = 3, at = 1)
dev.off()

# repeat with CI

pred_wglmfit_winner_oracle <- predict(wglmfit_winner_oracle, newdata = xgr_df, se.fit = TRUE)

pdf("fit_oracle_sigma_restricted_winner_CI.pdf")
par(mar = c(4, 5.3, 1, 3.5))
plot(xgr, pred_oracle$fit, col = "red", lwd = 3, type = "l", xlab = "time",
     ylab = "Duration (days)", cex.axis = 1.9, cex.lab = 2.3, ylim = ylimperm * 0.8)

polygon(c(xgr,rev(xgr)), c(pred_oracle$fit - 2 * pred_oracle$se.fit, rev(pred_oracle$fit + 2 * pred_oracle$se.fit)), col = "grey", border = "grey")

lines(xgr, pred_oracle$fit, col = "red", lwd = 3, type = "l")

polygon(c(xgr,rev(xgr)), c(pred_wglmfit_winner_oracle$fit - 2 * pred_wglmfit_winner_oracle$se.fit, rev(pred_wglmfit_winner_oracle$fit + 2 * pred_wglmfit_winner_oracle$se.fit)), col = "grey", border = "grey")
lines(xgr, predict(wglmfit_winner_oracle, xgr_df), col = "purple", lwd = 3, type = "l")
mtext("01/01/2009", side = 1, cex = 2, col = "blue", line = 3, at = 0)
mtext("12/31/2021", side = 1, cex = 2, col = "blue", line = 3, at = 1)
dev.off()

                                        #
pdf("fit_oracle_vs_naive_sigma_restricted.pdf")
par(mar = c(4, 5.3, 1, 3.5))
plot(xgr, predict(gamperm, xgr_df), col = "darkgreen", lwd = 3, type = "l",
     xlab = "time", ylab = "Duration (days)", cex.axis = 1.9, cex.lab = 2.3)
lines(xgr, predict(gamoracle, xgr_df), col = "red", lwd = 3, type = "l",
      xlab = "time", ylab = "Duration (days)", cex.axis = 1.9, cex.lab = 2.3)
lines(xgr, predict(wglmfit_winner, xgr_df), col = "orange", lwd = 3, type = "l")
lines(xgr, predict(wglmfit_winner_oracle, xgr_df), col = "purple", lwd = 3, type = "l")
mtext("01/01/2009", side = 1, cex = 2, col = "blue", line = 3, at = 0)
mtext("12/31/2021", side = 1, cex = 2, col = "blue", line = 3, at = 1)
dev.off()

# now re-create this plot with CIs


pdf("fit_oracle_vs_naive_sigma_restricted_CI.pdf")
par(mar = c(4, 5.3, 1, 3.5))
plot(xgr, predict(gamperm, xgr_df), col = "darkgreen", lwd = 3, type = "l",
     xlab = "time", ylab = "Duration (days)", cex.axis = 1.9, cex.lab = 2.3)

polygon(c(xgr,rev(xgr)), c(pred_oracle$fit - 2 * pred_oracle$se.fit, rev(pred_oracle$fit + 2 * pred_oracle$se.fit)), col = "grey", border = "grey")

polygon(c(xgr,rev(xgr)), c(pred_perm$fit - 2 * pred_perm$se.fit, rev(pred_perm$fit + 2 * pred_perm$se.fit)), col = "grey", border = "grey")

polygon(c(xgr,rev(xgr)), c(pred_wglmfit_winner$fit - 2 * pred_wglmfit_winner$se.fit, rev(pred_wglmfit_winner$fit + 2 * pred_wglmfit_winner$se.fit)), col = "grey", border = "grey")


polygon(c(xgr,rev(xgr)), c(pred_wglmfit_winner_oracle$fit - 2 * pred_wglmfit_winner_oracle$se.fit, rev(pred_wglmfit_winner_oracle$fit + 2 * pred_wglmfit_winner_oracle$se.fit)), col = "grey", border = "grey")


lines(xgr, predict(gamperm, xgr_df), col = "darkgreen", lwd = 3, type = "l",
      xlab = "time", ylab = "Duration (days)", cex.axis = 1.9, cex.lab = 2.)
lines(xgr, predict(gamoracle, xgr_df), col = "red", lwd = 3, type = "l",
      xlab = "time", ylab = "Duration (days)", cex.axis = 1.9, cex.lab = 2.)
lines(xgr, predict(wglmfit_winner, xgr_df), col = "orange", lwd = 3, type = "l")
lines(xgr, predict(wglmfit_winner_oracle, xgr_df), col = "purple", lwd = 3, type = "l")
mtext("01/01/2009", side = 1, cex = 2, col = "blue", line = 3, at = 0)
mtext("12/31/2021", side = 1, cex = 2, col = "blue", line = 3, at = 1)

dev.off()

### [1.5] BAYESIAN: attempt adjustment for generously linked data set

#rm(list = ls())
load("credentials_revisited_generously_linked.Rdata")
library(mgcv)
# now sigmagrid adjustment
sigmagr <- c(0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1)

for(jj in 1:length(sigmagr)){

set.seed(1104)

nsamples <- 1E3 # total number of samples
thinning <- 10  # take only every 10th sample
total_samples <- thinning * nsamples
burnin <- 100   # 100 "burn-in" iterations

d <- ncol(Bgam_perm)
dstar <-  d - sum(svd(Sgam_perm)$d < 1E-10)

# note superscript "mm" for mismatches
beta_samples_mm <- matrix(nrow = d, ncol = nsamples)
sigmasq_samples_mm <- numeric(nsamples)
tausq_samples_mm <- numeric(nsamples)
z_mm <- numeric(nperm) # z represents mismatch indicator
alpha_samples_mm <- numeric(nsamples)


beta_naive <- coef(gamperm)
RSS_it <- .5 * sum(residuals(gamperm)^2)
Qf_betait <- sum(beta_naive * (Sgam_perm %*% beta_naive))
z_it <- rep(0, nperm)
alpha_it <- 0.2
beta_it <- beta_naive
for(i in 1:burnin){


    # sigmasq-update
    sigmasq_it <- (sigmagr[jj] * sqrt(mean(residuals(gamperm)^2)))^2 #gamperm$sig2 #(1 / rgamma(1, shape = sum(z_it == 0)/2, scale = 2/RSS_it))

    # tausq-update
    tausq_it <- (1 / rgamma(1, shape = dstar/2, scale = 2/Qf_betait))

    # beta-update
    lambda_it <- sigmasq_it / tausq_it

    Omega_it <- crossprod(Bgam_perm[z_it == 0,]) + lambda_it * Sgam_perm
    svd_Omega_it <- svd(Omega_it)
    Gamma_it <- svd_Omega_it$u
    Theta_it <- svd_Omega_it$d
    rhs_it <- crossprod(Bgam_perm[z_it == 0,], yperm[z_it == 0])
    mu_it <- Gamma_it %*% (crossprod(Gamma_it, rhs_it)/Theta_it)
    covroot_it <- sqrt(sigmasq_it) * Gamma_it
    beta_it <- mu_it + (covroot_it %*% (rnorm(d)/sqrt(Theta_it)))

    Qf_betait <- sum(beta_it * (Sgam_perm %*% beta_it))


    # z-update
    res_it <- (yperm - Bgam_perm %*% beta_it)


    pis_num <- alpha_it * fy
    pis_denom <- pis_num + (1-alpha_it) * dnorm(res_it, sd = sqrt(sigmasq_it))
    pis <- pis_num / pis_denom
    z_it <- rbinom(nperm, size = 1, prob = pis)

    # alpha-update
    alpha_it <- rbeta(1, sum(z_it) + 1, nperm - sum(z_it) +1)

    #
    #res_it <- (yperm[z_it == 0] - Bgam_perm[z_it == 0,] %*% beta_it)
    #RSS_it <- sum(res_it^2)
}


for(i in 1:total_samples){

   # sigmasq-update
    #res_it <- (yperm[z_it == 0] - Bgam_perm[z_it == 0,] %*% beta_it)
    #RSS_it <- sum(res_it^2)
    sigmasq_it <- (sigmagr[jj] * sqrt(mean(residuals(gamperm)^2)))^2#(1 / rgamma(1, shape = sum(z_it == 0)/2, scale = 2/RSS_it))

    # tausq-update
    tausq_it <- (1 / rgamma(1, shape = dstar/2, scale = 2/Qf_betait))

    # beta-update
    lambda_it <- sigmasq_it / tausq_it

    Omega_it <- crossprod(Bgam_perm[z_it == 0,]) + lambda_it * Sgam_perm
    svd_Omega_it <- svd(Omega_it)
    Gamma_it <- svd_Omega_it$u
    Theta_it <- svd_Omega_it$d
    rhs_it <- crossprod(Bgam_perm[z_it == 0,], yperm[z_it == 0])
    mu_it <- Gamma_it %*% (crossprod(Gamma_it, rhs_it)/Theta_it)
    covroot_it <- sqrt(sigmasq_it) * Gamma_it
    beta_it <- mu_it + (covroot_it %*% (rnorm(d)/sqrt(Theta_it)))

    Qf_betait <- sum(beta_it * (Sgam_perm %*% beta_it))


    # z-update
    res_it <- (yperm - Bgam_perm %*% beta_it)
    pis_num <- alpha_it * fy
    pis_denom <- pis_num + (1-alpha_it) * dnorm(res_it, sd = sqrt(sigmasq_it))
    pis <- pis_num / pis_denom
    z_it <- rbinom(nperm, size = 1, prob = pis)

    # alpha-update
    alpha_it <- rbeta(1, sum(z_it) + 1, nperm - sum(z_it) +1)

   # store every 10th sample
   if((i %% 10) == 0){

       beta_samples_mm[,i/10] <- beta_it
       #sigmasq_samples_mm[i/10] <- sigmasq_it
       tausq_samples_mm[i/10] <- tausq_it
       alpha_samples_mm[i/10] <- alpha_it
       z_mm <- z_mm  + z_it / nsamples

   }
}

write.table(beta_samples_mm, file = paste("beta", jj, sep= "_"), row.names = TRUE, col.names = TRUE)
write.table(tausq_samples_mm, file = paste("tausq", jj, sep= "_"),row.names = TRUE, col.names = TRUE)
write.table(alpha_samples_mm, file = paste("alpha", jj, sep= "_"), row.names = TRUE, col.names = TRUE)
write.table(z_mm, file = paste("z", jj, sep= "_"), row.names = TRUE, col.names = TRUE)

}


### [1.6] BAYESIAN: attempt adjustment for restrictively linked data set

#rm(list = ls())
load("credentials_revisited_restrictively_linked.Rdata")
library(mgcv)
# now sigmagrid adjustment
sigmagr <- c(0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1)

for(jj in 1:length(sigmagr)){

set.seed(1104)

nsamples <- 1E3 # total number of samples
thinning <- 10  # take only every 10th sample
total_samples <- thinning * nsamples
burnin <- 100   # 100 "burn-in" iterations

d <- ncol(Bgam_oracle)
dstar <-  d - sum(svd(Sgam_oracle)$d < 1E-10)

# note superscript "mm" for mismatches
beta_samples_mm <- matrix(nrow = d, ncol = nsamples)
sigmasq_samples_mm <- numeric(nsamples)
tausq_samples_mm <- numeric(nsamples)
z_mm <- numeric(noracle) # z represents mismatch indicator
alpha_samples_mm <- numeric(nsamples)


beta_naive <- coef(gamoracle)
RSS_it <- .5 * sum(residuals(gamoracle)^2)
Qf_betait <- sum(beta_naive * (Sgam_oracle %*% beta_naive))
z_it <- rep(0, noracle)
alpha_it <- 0.2
beta_it <- beta_naive
for(i in 1:burnin){


    # sigmasq-update
    sigmasq_it <- (sigmagr[jj] * sqrt(mean(residuals(gamoracle)^2)))^2 #gamperm$sig2 #(1 / rgamma(1, shape = sum(z_it == 0)/2, scale = 2/RSS_it))

    # tausq-update
    tausq_it <- (1 / rgamma(1, shape = dstar/2, scale = 2/Qf_betait))

    # beta-update
    lambda_it <- sigmasq_it / tausq_it

    Omega_it <- crossprod(Bgam_oracle[z_it == 0,]) + lambda_it * Sgam_oracle
    svd_Omega_it <- svd(Omega_it)
    Gamma_it <- svd_Omega_it$u
    Theta_it <- svd_Omega_it$d
    rhs_it <- crossprod(Bgam_oracle[z_it == 0,], yoracle[z_it == 0])
    mu_it <- Gamma_it %*% (crossprod(Gamma_it, rhs_it)/Theta_it)
    covroot_it <- sqrt(sigmasq_it) * Gamma_it
    beta_it <- mu_it + (covroot_it %*% (rnorm(d)/sqrt(Theta_it)))

    Qf_betait <- sum(beta_it * (Sgam_oracle %*% beta_it))


    # z-update
    res_it <- (yoracle - Bgam_oracle %*% beta_it)


    pis_num <- alpha_it * fy
    pis_denom <- pis_num + (1-alpha_it) * dnorm(res_it, sd = sqrt(sigmasq_it))
    pis <- pis_num / pis_denom
    z_it <- rbinom(noracle, size = 1, prob = pis)

    # alpha-update
    alpha_it <- rbeta(1, sum(z_it) + 1, noracle - sum(z_it) +1)

    #
    #res_it <- (yperm[z_it == 0] - Bgam_perm[z_it == 0,] %*% beta_it)
    #RSS_it <- sum(res_it^2)
}


for(i in 1:total_samples){

   # sigmasq-update
    #res_it <- (yperm[z_it == 0] - Bgam_perm[z_it == 0,] %*% beta_it)
    #RSS_it <- sum(res_it^2)
    sigmasq_it <- (sigmagr[jj] * sqrt(mean(residuals(gamoracle)^2)))^2#(1 / rgamma(1, shape = sum(z_it == 0)/2, scale = 2/RSS_it))

    # tausq-update
    tausq_it <- (1 / rgamma(1, shape = dstar/2, scale = 2/Qf_betait))

    # beta-update
    lambda_it <- sigmasq_it / tausq_it

    Omega_it <- crossprod(Bgam_oracle[z_it == 0,]) + lambda_it * Sgam_oracle
    svd_Omega_it <- svd(Omega_it)
    Gamma_it <- svd_Omega_it$u
    Theta_it <- svd_Omega_it$d
    rhs_it <- crossprod(Bgam_oracle[z_it == 0,], yoracle[z_it == 0])
    mu_it <- Gamma_it %*% (crossprod(Gamma_it, rhs_it)/Theta_it)
    covroot_it <- sqrt(sigmasq_it) * Gamma_it
    beta_it <- mu_it + (covroot_it %*% (rnorm(d)/sqrt(Theta_it)))

    Qf_betait <- sum(beta_it * (Sgam_oracle %*% beta_it))


    # z-update
    res_it <- (yoracle - Bgam_oracle %*% beta_it)
    pis_num <- alpha_it * fy
    pis_denom <- pis_num + (1-alpha_it) * dnorm(res_it, sd = sqrt(sigmasq_it))
    pis <- pis_num / pis_denom
    z_it <- rbinom(noracle, size = 1, prob = pis)

    # alpha-update
    alpha_it <- rbeta(1, sum(z_it) + 1, noracle - sum(z_it) +1)

   # store every 10th sample
   if((i %% 10) == 0){

       beta_samples_mm[,i/10] <- beta_it
       #sigmasq_samples_mm[i/10] <- sigmasq_it
       tausq_samples_mm[i/10] <- tausq_it
       alpha_samples_mm[i/10] <- alpha_it
       z_mm <- z_mm  + z_it / nsamples

   }
}

write.table(beta_samples_mm, file = paste("beta-oracle", jj, sep= "_"), row.names = TRUE, col.names = TRUE)
write.table(tausq_samples_mm, file = paste("tausq-oracle", jj, sep= "_"), row.names = TRUE, col.names = TRUE)
write.table(alpha_samples_mm, file = paste("alpha-oracle", jj, sep= "_"), row.names = TRUE, col.names = TRUE)
write.table(z_mm, file = paste("z-oracle", jj, sep= "_"), row.names = TRUE, col.names = TRUE)

}





