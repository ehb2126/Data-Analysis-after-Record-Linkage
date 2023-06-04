###################
### LIFE-M DEMO ###
###################

# 1. Install and load "pldamixture" package
# download "pldamixture_0.0.0.9000.tar.gz" and set working directory to its location
install.packages("pldamixture_0.0.0.9000.tar.gz", repos = NULL, type = "source")
library(pldamixture)

# 2. Read in demo data (2 (HL): 1 (ML) SRS Ratio - 3,238 Individuals)
data("demodata")

# 3. Plain Linear Regression: Naïve (n = 3,238)
naive_fit <- lm(surv_age ~ poly(uyob,3,raw = TRUE), demodata)
summary(naive_fit)

# 4. Plain Linear Regression: Hand-Linked Only (n = 2,159)
hldemodata <- demodata[demodata$hndlnk == "Hand-Linked At Some Level",]
hl_fit <- lm(surv_age ~ poly(uyob,3,raw = TRUE), hldemodata)
summary(hl_fit)

# 5. Adjustment Method (n = 3,238)
?fit_mixture()
?summary_fit_mixture()
am_fit <- fit_mixture(surv_age ~ poly(uyob, 3, raw = TRUE), 	 		demodata, m.formula = ~ commf + comml, 
                      safematches = ifelse(demodata$hndlnk == 				              "Hand-Linked At Some Level", TRUE, FALSE), 
                      assumed_mrate = 0.05) 
summary_fit_mixture(am_fit)


#################################################
### PREDICTED AGES AT DEATH FOR BIRTH COHORTS ###
#################################################

# custom functions to obtain predictions + 95% point-wise CI
fitci <- function(dat, type, model){
  if (type == "lm"){
    x = dat$yob
    x = (x - min(x)) / (max(x) - min(x))
    y = dat$surv_age
    
    model <- lm(y ~ poly(x, 3, raw = TRUE))
    pred <- predict(model, dat, se.fit = TRUE, interval = "confidence", level = 0.95)$fit
  }
  
  if (type == "Adjustment"){
    x = dat$yob
    x = (x - min(x)) / (max(x) - min(x))
    
    # https://stackoverflow.com/questions/38109501/how-does-predict-lm-compute-confidence-interval-and-prediction-interval
    Xp <- cbind(intercept = rep(1,nrow(dat)), x, x2 = x^2, x3 = x^3) # design matrix
    b <- model$betahat
    yh <- c(Xp %*% b) # fitted values
    
    sig2 <- model$sigmahat^2
    V <- solve(t(Xp) %*% Xp)*sig2 #variance-covariance matrix (vcov())
    var.fit <- rowSums((Xp %*% V) * Xp)
    se.fit <- sqrt(var.fit)
    
    dof <- nrow(dat) - 3 - 1
    
    alpha <- 0.95
    qt <- c(-1,1) * qt((1-alpha) / 2, dof, lower.tail = FALSE)
    
    ci <- yh + outer(se.fit, qt)
    
    pred <- cbind(fit = yh, lwr = ci[,1], upr = ci[,2])
  }
  pred
}

# obtain predictions + 95% point-wise CI based on each analysis method
predam <- as.data.frame(cbind(yob = demodata$yob, fitci(dat = demodata, type = "Adjustment", model = am_fit)))
predn <- as.data.frame(cbind(yob = demodata$yob, fitci(dat = demodata, type = "lm")))
predh <- as.data.frame(cbind(yob = hldemodata$yob, fitci(dat = hldemodata, type = "lm")))

# obtain plot
library(ggplot2)
ggplot() +
  geom_line(data = predam, aes(x = yob, y = fit, color = "Adjustment"), linewidth = 1) +
  geom_ribbon(data = predam, aes(x = yob, y = fit, ymin = lwr, ymax = upr), fill = "#D55E00", alpha = 0.2) +
  geom_line(data = predn, aes(x = yob, y = fit, color = "Naive"), linewidth = 1, linetype = "longdash") +
  geom_line(data = predh, aes(x = yob, y = fit, color = "Hand-Linked Only"), linewidth = 1, linetype = "longdash") +
  geom_ribbon(data = predn, aes(x = yob, y = fit, ymin = lwr, ymax = upr), fill = "seagreen", alpha = 0.2) +
  geom_ribbon(data = predh, aes(x = yob, y = fit, ymin = lwr, ymax = upr), fill = "black", alpha = 0.2) +
  xlab("Year of Birth (x)") +
  ylab("Age at Death (y)") +
  theme(plot.background = element_rect(colour = "darkgray", linewidth = 2)) +
  scale_x_continuous(breaks = seq(1880,1910, 2)) +
  scale_y_continuous(breaks = seq(0,100, 5)) +
  theme_bw(base_size = 18) +
  labs(color = "Approach") +
  scale_color_manual(values = c("Naive" = "seagreen", 
                                "Hand-Linked Only" = "black",
                                "Adjustment" = "#D55E00")) +
  theme(legend.position = "bottom")