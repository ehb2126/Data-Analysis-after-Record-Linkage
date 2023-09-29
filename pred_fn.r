
library(readr)
library(dplyr)


mse <- function(x1,x2){
  mean((x1 - x2)^2)
}

# Defining the bagging trees

library(rpart)

pred_bag_trees <- function(formula, data, newdata, ntrees=100, wghts){
  n <- nrow(data)
  m <- nrow(newdata)
  new_mu_yx <- numeric(m)
  for(i in 1:ntrees){
    samp <- sample(1:n, replace = TRUE)
    dsamp = data[samp,]
    tr <- rpart(formula,
                data = data, control = rpart.control(minsplit = 20),
                weights = wghts, method = "anova")
    yhat <- predict(tr, newdata = newdata)
    new_mu_yx <- new_mu_yx + yhat/ntrees
  }
  
  return(new_mu_yx)
}


# permute 20% of the data and see how the results change:
# Permutation function, permuting a specified fraction of y

y.perm <- function(y,fraction=.2){
  k = floor(fraction*length(y))
  yperm <- y
  yperm[1:k] <- y[sample(1:k)]
  return(yperm)
}






nloglik <- function(y, mu_yx,mu_y,sigma_yx, sigma_y, Alpha){
  -sum(log(diag(1 - Alpha)%*%dnorm(y,mu_yx, sd = sigma_yx)
           + diag(Alpha)%*%dnorm(y, mu_y, sd = sigma_y)))
}


# Defining adjusted bagging trees: using nloglik for assigning 
# better weights to bagging trees function 



pred_adj_bag_trees1 <-function(formula, data, newdata,  ntrees=100, wghts, Alpha){
  
  new_mu_yx <- pred_bag_trees(formula, data, newdata,  ntrees, wghts)
  mu_yxcur <- bag_trees(formula, data, ntrees, wghts = wghts)
  n <- nrow(data)
  m <- nrow(newdata)
  df <- model.frame(formula, data)
  y <- model.response(df)
  sigma_yx <- sigma(lm(formula, df))
  sigma_y <- sd(y)
  mu_y <- mean(y)
  tol <- 1E-4
  maxiter <- 1000
  iter <- 1
  objs <- numeric(maxiter)
  objs[1] <- nloglik(y, mu_yxcur, mu_y, sigma_yx, sigma_y, Alpha)
  while(iter < maxiter){
    
    num <- diag(1-Alpha)%*%dnorm(y,mu_yxcur,sigma_yx)
    denom <- num + diag(Alpha)%*%dnorm(y,mu_y,sigma_y)
    wghts <- num/denom
    Alpha <- 1 - wghts
    new_mu_yx <- pred_bag_trees(formula, data, newdata,  ntrees, wghts)
    mu_yxcur <- bag_trees(formula, data, ntrees, wghts = wghts) 
    iter <- iter + 1
    objs[iter] <- nloglik(y, mu_yxcur, mu_y, sigma_yx, sigma_y, Alpha)
    if(objs[iter] + tol > objs[iter-1])
      break
  }
  return(new_mu_yx)
}






# The adjusted rpart version 2 is defined as follows.

pred_adj_rpart2 <- function(formula, data, newdata,  wghts, Alpha){
  
  n <- nrow(data)
  m <- nrow(newdata)
  df <- model.frame(formula, data)
  y <- model.response(df)
  sigma_yx <- sigma(lm(formula, df))
  sigma_y <- sd(y)
  mu_y <- mean(y)
  # First the usual decision tree applied 
  tr <-  rpart(formula,
               data = df, control = rpart.control(minsplit = 20),
               weights = wghts, 
               method = "anova")
  new_mu_yx <- predict(tr, newdata = newdata)
  mu_yxcur <- predict(tr, newdata = data)
  
  tol <- 1E-4
  maxiter <- 1000
  objs <- numeric(maxiter)
  objs[1] <- nloglik(y, mu_yxcur, mu_y, sigma_yx, 
                     sigma_y, Alpha)
  iter <- 1
  while(iter < maxiter){
    
    num <- diag(1-Alpha)%*%dnorm(y,mu_yxcur, sd = sigma_yx)
    denom <- num + diag(Alpha)%*%dnorm(y,mu_y, sd = sigma_y)
    wghts <- num/denom
    Alpha <- 1 - wghts
    
    # need predicted values here
    tr <-  rpart(formula,
                 data = df, control = rpart.control(minsplit = 20), 
                 weights=wghts, method = "anova")
    mu_yxcur <- predict(tr, newdata = data)
    new_mu_yx <- predict(tr, newdata = newdata)
    iter <- iter + 1
    objs[iter] <- nloglik(y, mu_yxcur, mu_y, sigma_yx, sigma_y, Alpha)
    if(objs[iter] + tol > objs[iter-1])
      break
  }
  
  return(new_mu_yx)
}


# The version 2 of adjusted bagging trees:



pred_adj_bag_trees2 <- function(formula, data, newdata,  ntrees = 100, wghts, Alpha){
  n <- nrow(data)
  m <- nrow(newdata)
  new_mu_yx <- numeric(m)
  
  
  for(i in 1:ntrees){
    samp <- sample(1:n, replace = TRUE)
    dsamp = data[samp,]
    wghts = wghts[samp]
    yhat <- pred_adj_rpart2(formula, data =  dsamp, newdata = newdata, wghts, Alpha)
    new_mu_yx <- new_mu_yx + yhat/ntrees
    
  }
  
  new_mu_yx
}





# The 3rd version of adjusted bagging trees is based on 



adj_rpart_pred3 <- function(formula, data, newdata,  wghts, Alpha){
  
  # need predicted values
  n <- nrow(data)
  df <- model.frame(formula, data)
  y <- model.response(df)
  sigma_yx <- sigma(lm(formula, df))
  sigma_y <- sd(y)
  mu_y <- mean(y)
  
  samp <- sample(1:n, replace = TRUE)
  dsamp = df[samp,]
  wghts = wghts[samp]
  
  tr <-  rpart(formula,
               dsamp, control = rpart.control(minsplit = 20), wghts,  
               method = "anova")
  mu_yxcur <- predict(tr, newdata = df)
  Alpha_sum <- Alpha
  Alpha_hat <- numeric(n)
  
  for(i in 1:n){
    tol <- 1E-3
    maxiter <- 1000
    alpha <- Alpha[i]
    num <- (1-alpha) * dnorm(y[i], mu_yxcur[i], sd = sigma_yx)
    denom <- num + alpha * dnorm(y[i], mu_y, sd = sigma_y)
    objs <- numeric(maxiter)
    objs[1] <- Alpha[i]
    iter <- 1
    r <- 1/(1 - Alpha[i])
    while(iter < maxiter){
      mu_yxcur[i] <- r*mu_yxcur[i] + (1 - r)*mu_y
      num <- (1-Alpha[i]) * dnorm(y[i], mu_yxcur[i], sd = sigma_yx)
      denom <- num + Alpha[i] * dnorm(y[i], mu_y, sd = sigma_y)
      iter <- iter + 1
      objs[iter] <- 1 - 1/r
      r <- denom/num
      Alpha_sum[i] <- Alpha_sum[i] + (1 - 1/r)
      if(abs(objs[iter] - objs[iter-1]) < tol)
        break
      
    }
    Alpha_hat[i] <- Alpha_sum[i]/iter
  }
  mu_yxcur
}



# Defining adjusted bagging trees, replacing rpart with 
# *adjusted* rpart

pred_adj_bag_trees3 <- function(formula, data, newdata,  ntrees=100, wghts, Alpha){
  #browser()
  n <- nrow(data)
  mu_yxcur <- numeric(n)
  Alpha_hat <- numeric(n)
  for(i in 1:ntrees){
    samp <- sample(1:n, replace = TRUE)
    dsamp = data[samp,]
    wghts = wghts[samp]
    adj3 <- adj_rpart_pred3(formula, dsamp, wghts, Alpha)
    mu_yxcur <- mu_yxcur + adj3/ntrees
  }
  
  return(mu_yxcur)
}



pred_adj_rf1<-function(formula, data, newdata, Alpha, maxiter = 1000){
  
  tr <- ranger(formula, data)
  mu_yxcur <- tr$predictions
  new_mu_yx <- predict(tr, data = newdata)$predictions
  n <- nrow(data)
  df <- model.frame(formula, data)
  y <- model.response(df)
  sigma_yx <- sigma(lm(formula, df))
  sigma_y <- sd(y)
  mu_y <- mean(y)
  tol <- 1E-4
  iter <- 1
  objs <- numeric(maxiter)
  objs[1] <- nloglik(y, mu_yxcur, mu_y, sigma_yx, sigma_y, Alpha)
  while(iter < maxiter){
    
    num <- diag(1-Alpha)%*%dnorm(y,mu_yxcur,sigma_yx)
    denom <- num + diag(Alpha)%*%dnorm(y,mu_y,sigma_y)
    wghts <- num/denom
    Alpha <- 1 - wghts
    tr <- ranger(formula, data, case.weights = wghts)
    mu_yxcur <- tr$predictions
    new_mu_yx <- predict(tr, data = newdata)$predictions
    iter <- iter + 1
    objs[iter] <- nloglik(y, mu_yxcur, mu_y, sigma_yx, sigma_y, Alpha)
    if(objs[iter] + tol > objs[iter-1])
      break
  }
  return(new_mu_yx)
}




