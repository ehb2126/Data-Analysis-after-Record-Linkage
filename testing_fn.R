
library(readr)
library(dplyr)

mse <- function(x1,x2){
   mean((x1 - x2)^2)
}

# Defining the bagging trees

library(rpart)

bag_trees <- function(formula, data, ntrees=100, wghts){
  n <- nrow(data)
  muhat <- numeric(n)
  for(i in 1:ntrees){
    samp <- sample(1:n, replace = TRUE)
    dsamp = data[samp,]
    tr <- rpart(formula,
                data = dsamp, control = rpart.control(minsplit = 20),
                weights = wghts, method = "anova")
    yhat <- predict(tr, newdata = data)
    muhat <- muhat + yhat/ntrees
  }
  
  return(muhat)
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

adj_bag_trees1<-function(formula, data, ntrees=100, wghts, Alpha){
  
  mu_yxcur <- bag_trees(formula, data, ntrees)
  n <- nrow(data)
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
    mu_yxcur <- bag_trees(formula, data, ntrees, wghts)
    iter <- iter + 1
    objs[iter] <- nloglik(y, mu_yxcur, mu_y, sigma_yx, sigma_y, Alpha)
    if(objs[iter] + tol > objs[iter-1])
      break
  }
  return(mu_yxcur)
}




# The adjusted rpart version 2 is defined as follows.

adj_rpart2 <- function(formula, data, wghts, Alpha){
  
  n <- nrow(data)
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
  mu_yxcur <- predict(tr, newdata = df)
  
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
    mu_yxcur <- predict(tr, newdata = df)
    iter <- iter + 1
    objs[iter] <- nloglik(y, mu_yxcur, mu_y, sigma_yx, sigma_y, Alpha)
    if(objs[iter] + tol > objs[iter-1])
      break
  }
  
  return(tr)
}


# The version 2 of adjusted bagging trees:



adj_bag_trees2 <- function(formula, data, ntrees = 100, wghts, Alpha){
  
  n <- nrow(data)
  muhat <- numeric(n)
  
  
  for(i in 1:ntrees){
    samp <- sample(1:n, replace = TRUE)
    dsamp = data[samp,]
    wghts = wghts[samp]
    tr <- adj_rpart2(formula, dsamp, wghts, Alpha)
    yhat <- predict(tr, newdata = data)
    muhat <- muhat + yhat/ntrees
    
  }
  
  muhat
}





# The 3rd version of adjusted bagging trees is based on 



adj_rpart3 <- function(formula, data, wghts, Alpha){
  
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

adj_bag_trees3 <- function(formula, data, ntrees=100, wghts, Alpha){
  #browser()
  n <- nrow(data)
  mu_yxcur <- numeric(n)
  Alpha_hat <- numeric(n)
  for(i in 1:ntrees){
    samp <- sample(1:n, replace = TRUE)
    dsamp = data[samp,]
    wghts = wghts[samp]
    adj3 <- adj_rpart3(formula, dsamp, wghts, Alpha)
    mu_yxcur <- mu_yxcur + adj3/ntrees
  }
  
  return(mu_yxcur)
}



# The version 4th of rpart, in this version alpha is
# directly estimated based on maximizing the 
# likelihood function of the mixture density


optim_alpha_rpart <- function(formula, data){
  
  # need predicted values
  
  n <- nrow(data)
  df <- model.frame(formula, data)
  y <- model.response(df)
  sigma_yx <- sigma(lm(formula, df))
  sigma_y <- sd(y)
  mu_y <- mean(y)
  
  
  
  tr <-  rpart(formula,
               data = data, control = rpart.control(minsplit = 20),
               method = "anova")
  muhat <- predict(tr, newdata = df)
  
  mu_yxcur <- numeric(n)
  
  Alpha_hat <- numeric(n)
  for(i in 1:n){
    
    fn <- function(alpha){
      r = 1/(1 - alpha)
      mu_yxcur[i] <- r*muhat[i] + (1 - r)*mu_y
      num <- (1-alpha) * dnorm(y[i], mu_yxcur[i], sigma_yx)
      denom <- num + alpha * dnorm(y[i], mu_y, sigma_y)
      return(denom)
    }
   
    Alpha_hat[i] <- optimize(fn,c(.01,.9), 
                      maximum = TRUE, tol = 1E-4)$maximum
  }
  return(Alpha_hat)
}



# Adjusted version of bagging trees based on adj_rpart


optimal_Alpha <- function(formula, data, ntrees=100){
  n <- nrow(data)
  muhat <- numeric(n)
  Alpha_hat <- numeric(n)
  for(iter in 1:ntrees){
    samp <- sample(1:n, replace = TRUE)
    dsamp = data[samp,]
    Alph_optim <- optim_alpha_rpart(formula, dsamp)
    Alpha_hat <- Alpha_hat + Alph_optim/ntrees
  }
  return(Alpha_hat)
}

adj_alpha_method <- function(formula, data, Alpha, mu_ytildex, ntrees){
  df <- model.frame(formula, data)
  y <- model.response(df)
  mu_y <- mean(y)
  R <- 1/(1 - Alpha)
  mu_yx <- R*mu_ytildex + (1 - R)*mu_y
  return(mu_yx)
  }

# This is adjustment applied to random forest by reweighing case.weights in ranger 
# implementation of random forests

adj_rf1<-function(formula, data, Alpha, maxiter = 1000){
  
  mu_yxcur <- ranger(formula, data)$predictions
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
    iter <- iter + 1
    objs[iter] <- nloglik(y, mu_yxcur, mu_y, sigma_yx, sigma_y, Alpha)
    if(objs[iter] + tol > objs[iter-1])
      break
  }
  return(mu_yxcur)
}


# This is a less computationally intensive implementation of reweighing method.


adj_rf2 <- function(formula, data, Alpha, maxiter = 1000){
  
  n <- nrow(data)
  muhat_sum <- numeric(n)
  muhat <- numeric(n)
  df <- model.frame(formula, data)
  y <- model.response(df)
  sigma_yx <- sigma(lm(formula, df))
  sigma_y <- sd(y)
  mu_y <- mean(y)
  
  num_breaks = max(floor(n/200),5)
  
  folds <- cut(seq(1, nrow(sim_dat)), breaks = num_breaks, labels = FALSE)
  for(J in 1:num_breaks){
    trainIndexes <- which(folds == J, arr.ind = TRUE)
    trainData <- data[trainIndexes, ]
    
    
    tr <-  ranger(formula, data = trainData)
    
    mu_yxcur <- predict(tr, data = df)$predictions
    
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
      tr <-  ranger(formula, data = trainData, case.weights = wghts[trainIndexes])
      mu_yxcur <- predict(tr, data = df)$predictions
      muhat <- (muhat + mu_yxcur)/iter 
      iter <- iter + 1
      objs[iter] <- nloglik(y, mu_yxcur, mu_y, sigma_yx, sigma_y, Alpha)
      if(objs[iter] + tol > objs[iter-1])
        break
    }
    muhat_sum <- muhat_sum + muhat
  }
  return(muhat_sum/num_breaks)
}



optimal_mean_alpha_ranger <- function(formula, data){
  
  # need predicted values
  
  n <- nrow(data)
  df <- model.frame(formula, data)
  y <- model.response(df)
  sigma_yx <- sigma(lm(formula, df))
  sigma_y <- sd(y)
  mu_y <- mean(y)
  
  muhat <-  ranger(formula,data)$predictions
  mu_yxcur <- numeric(n)
  
  Alpha_hat <- numeric(n)
  
  for(i in 1:n){
    
    fn <- function(alpha){
      r = 1/(1 - alpha)
      mu_yxcur[i] <- r*muhat[i] + (1 - r)*mu_y
      num <- (1-alpha) * dnorm(y[i], mu_yxcur[i], sigma_yx)
      denom <- num + alpha * dnorm(y[i], mu_y, sigma_y)
      return(denom)
    }
    
    Alpha_hat[i] <- optimize(fn,c(.01,.9), 
                             maximum = TRUE, tol = 1E-4)$maximum
    
  }
  
  return(mean(Alpha_hat))
}

# Adjusted version of bagging trees based on adj_rpart


optimal_mean_alpha <- function(formula, data, ratio = .5, maxiter = 20){
  alpha_mean_sum <- 0
  for(iter in 1:maxiter){
    n <- nrow(data)
    k <- floor(ratio*n)
    samp <- sample(1:n, k, replace = TRUE)
    dsamp = data[samp, ]
    alpha_mean <- optimal_mean_alpha_ranger(formula, data = dsamp)
    alpha_mean_sum <- alpha_mean_sum + alpha_mean
  }
  alpha_mean <- alpha_mean_sum/maxiter
  return(alpha_mean)
}



