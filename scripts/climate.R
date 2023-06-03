setwd("/home/emanuelben-david/Desktop/OneDrive_1_6-2-2023")
X <- read.csv("climate_X.csv", header = TRUE)
Y <- read.csv("climate_Y.csv", header = TRUE)
Y_mismatched <- read.csv("climate_Y_mismatched.csv", header = TRUE)

### PCA based on correctly matched data
pca <- princomp(~., data = cbind(X,Y), cor = TRUE)
# proportion of explained variance
cumsum(pca$sdev^2) / sum(pca$sdev^2)

### PCA based on partially mismatched data
pca_mm <- princomp(~., data = cbind(X,Y_mismatched), cor = TRUE)
# proportion of explained variance
cumsum(pca_mm$sdev^2) / sum(pca_mm$sdev^2)
