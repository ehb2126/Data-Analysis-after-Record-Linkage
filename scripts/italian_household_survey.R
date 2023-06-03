ls()
ihs <- read.csv("italian_household_survey.txt", header = FALSE)
colnames(ihs) <- c("x", "y_original", "y_mismatched")

# reproducing the results and figures on the slides

### w/o mismatches

plot(y_original ~ x, data = ihs, pch = 16, cex = .8, xlab = "Income 2008",
     ylab = "Income 2010", cex.lab = 1.8, cex.axis = 1.6)

lm_original <- lm(y_original ~ x, data = ihs)
abline(lm_original)
summary(lm_original)
summary(lm_original)$sigma^2

dev.off()
### w/ mismatches



plot(y_mismatched ~ x, data = ihs, pch = ifelse(ihs$y == ihs$y_mismatched, 16, 2),
     cex = .8, xlab = "Income 2008",
     ylab = "Income 2010", cex.lab = 1.8, cex.axis = 1.6)

lm_mismatched <- lm(y_mismatched ~ x, data = ihs)
abline(lm_mismatched,col = 'red')

summary(lm_mismatched)
summary(lm_mismatched)$sigma^2
