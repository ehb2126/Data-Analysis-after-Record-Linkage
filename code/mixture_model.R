### X:       --- n x p design matrix
### y:       --- n x 1 observed response
### f0:      --- density function of the marginal distribution of y
### control: --- additional parameters

fit_mixture <- function(X, y, f0, control = list(init = c("LS", "robust"), tol = 1E-4, maxit = 100, gamma = 0.95)){

    if(!is.null(control$init))
        init <- control$init[1]
    else
        init <- "LS"
    if(!is.null(control$tol))
        tol <- control$tol
    else
        tol <- 1E-4
    if(!is.null(control$maxit))
        maxit <- control$maxit
    else
        maxit <- 100
    if(!is.null(control$gamma))
        gamma <- control$gamma
    else
        gamma <- 0.95
   
    X <- as.matrix(X)
    n <- nrow(X)
    p <- ncol(X)
    
    y <- as.numeric(y)

    if(length(y) != n)
        stop("Error: Length of 'y' must be equal to the number of rows of 'X' \n")

    if(!is.function(f0))
        stop("Error: 'f0' must be a function (density) \n")

    
     if(!(init %in% c("LS", "robust")))
        stop("Error: 'init' must be either 'LS' or 'robust' \n")

    if(init == "LS"){
        betahat <- solve(crossprod(X), crossprod(X, y))
        sigmahat <- sqrt(sum((y - X %*% betahat)^2)/(n-p))
        alphahat <- 0.5    
    }

    if(init == "robust"){

        require(MASS)
        rlmsol <- rlm(x=X, y=y)
        betahat <- coef(rlmsol)
        sigmahat <- rlmsol$s
        betahat_LS <- solve(crossprod(X), crossprod(X, y))
        alphahat <- max(min(1, sqrt(sum(betahat_LS^2))/sqrt(sum(betahat^2))),0.01)     
    }


    nloglik <- function(beta, sigma, alpha){
               mean(-log(alpha * (1/sqrt(2 * pi * sigma^2)) * exp(-(y - (X %*% beta))^2/(2 * sigma^2)) 
                         + (1-alpha) * f0(y)))
    }


    it <- 0
    conv <- numeric(maxit)           
    nloglik_cur <- nloglik(betahat, sigmahat, alphahat)
    conv[it + 1] <- nloglik_cur
    
    ### iterations

    while(it < maxit){
        r <- X %*% betahat - y
        
        pnum <-  alphahat * (1/sqrt(2 * pi * sigmahat^2)) * exp(-r^2/(2*sigmahat^2))
        pdenom <- pnum + (1 - alphahat) * f0(y)
        pcur <- pnum / pdenom
        w <- pcur
       
        Xw <- X
        Xw <- sweep(X, 1, STATS = sqrt(w), FUN="*") 
        yw <- sqrt(w) * y

        betahatnew <- solve(crossprod(Xw), crossprod(Xw, yw))        
        sigmahatnew <- sqrt(sum((yw - Xw %*% betahatnew)^2) / sum(w))

        alphahatnew = mean(pcur);
        nloglik_new = nloglik(betahatnew, sigmahatnew, alphahatnew)

        if((nloglik_cur - nloglik_new) < -tol^2)
            break
        else{
            betahat = betahatnew;
            alphahat = alphahatnew;
            sigmahat = sigmahatnew;
            it <- it+1
            conv[it+1] = nloglik_new
            if((nloglik_cur - nloglik_new) < tol)
                break
            else{
                nloglik_cur <- nloglik_new
            }
       }

    }


     return(list(betahat = betahat, sigmahat = sigmahat, alphahat = alphahat, conv = conv[1:(it+1)]))

}




fit_mixture0 <- function(X, y, control = list(init = c("LS", "robust"), tol = 1E-4, maxit = 100, gamma = 0.95)){

    if(!is.null(control$init))
        init <- control$init[1]
    else
        init <- "LS"
    if(!is.null(control$tol))
        tol <- control$tol
    else
        tol <- 1E-4
    if(!is.null(control$maxit))
        maxit <- control$maxit
    else
        maxit <- 100
    if(!is.null(control$gamma))
        gamma <- control$gamma
    else
        gamma <- 0.95


    X <- as.matrix(X)
    n <- nrow(X)
    p <- ncol(X)
    
    y <- as.numeric(y)

    if(length(y) != n)
        stop("Error: Length of 'y' must be equal to the number of rows of 'X' \n")

    if(!(init %in% c("LS", "robust")))
        stop("Error: 'init' must be either 'LS' or 'robust' \n")

    if(init == "LS"){

        betahat <- solve(crossprod(X), crossprod(X, y))
        sigmahat <- sqrt(sum((y - X %*% betahat)^2)/(n-p))
        alphahat <- 0.5    
    }

    if(init == "robust"){

        require(MASS)
        rlmsol <- rlm(x=X, y=y)
        betahat <- coef(rlmsol)
        sigmahat <- rlmsol$s
        betahat_LS <- solve(crossprod(X), crossprod(X, y))
        alphahat <- max(min(1, sqrt(sum(betahat_LS^2))/sqrt(sum(betahat^2))),0.01)     
    }

    # functions for optimization

    normsq <- function(x) sum(x^2)
    nloglik <- function(beta, sigma, alpha){
               mean(-log(alpha * (1/sqrt(2 * pi * sigma^2)) * exp(-(y - (X %*% beta))^2/(2 * sigma^2)) 
                         + (1-alpha) * (1/sqrt(2 * pi * (normsq(beta) + sigma^2))) * exp(-y^2 / (2 * (normsq(beta) + sigma^2)))))
               }
    
    
    it <- 0
    conv <- numeric(maxit)           
    nloglik_cur <- nloglik(betahat, sigmahat, alphahat)
    conv[it + 1] <- nloglik_cur
    
    ### iterations

    while(it < maxit){
        r <- X %*% betahat - y
        
        pnum <-  alphahat * (1/sqrt(2 * pi * sigmahat^2)) * exp(-r^2/(2*sigmahat^2))
        pdenom <- pnum + (1 - alphahat) * (1/sqrt(2 * pi * (normsq(betahat) + sigmahat^2))) * exp(-y^2 / (2 * (normsq(betahat) + sigmahat^2)))
        pcur <- pnum / pdenom
        w <- pcur
       
        Xw <- X
        Xw <- sweep(X, 1, STATS = sqrt(w), FUN="*") 
        yw <- sqrt(w) * y

        rw <- Xw %*% betahat - yw
        
        gradbeta <- (1/sigmahat^2) * crossprod(Xw, rw)  - betahat * (sum((1-pcur)*y^2)/((sigmahat^2 + normsq(betahat))^2) - sum((1-pcur))/(sigmahat^2 + normsq(betahat)))
        gradsigma <- -normsq(rw)/(2*sigmahat^4) + sum(pcur)/(2*sigmahat^2) - sum((1-pcur)*y^2)/(2*(sigmahat^2 + normsq(betahat))^2) + 0.5*sum(1-pcur)/(sigmahat^2 + normsq(betahat))
        grad <- c(gradbeta, gradsigma)
 
        Ehess_bb <- crossprod(Xw)/sigmahat^2 + 2 * sum((1-pcur)/((sigmahat^2 + normsq(betahat))^2)) * outer(c(betahat), c(betahat))
        Ehess_bsigma <- sum((1-pcur)/((sigmahat^2 + normsq(betahat))^2)) * betahat
        Ehess_sigmasigma <- sum(0.5 * pcur / sigmahat^4) +  sum(0.5*(1-pcur)*(1/((sigmahat^2 + normsq(betahat))^2)))
        Ehess <- rbind(cbind(Ehess_bb, Ehess_bsigma), c(Ehess_bsigma, Ehess_sigmasigma))
        
        dir <- solve(Ehess, -grad)

        m <- 0
        
        while(TRUE){
            
            betahatnew <- betahat + gamma^m * dir[1:p]
            sigmahatnew <- sqrt(sigmahat^2 + (gamma^m) * dir[p + 1])
            alphahatnew <- sum(pcur)/n
            
            nloglik_new = nloglik(betahatnew, sigmahatnew, alphahatnew)
            if(is.na(nloglik_new))
                m <- m+1
            else{
                if(nloglik_cur - nloglik_new > tol){
                    nloglik_cur <- nloglik_new
                    betahat  <- betahatnew
                    sigmahat <- sigmahatnew
                    alphahat <- alphahatnew
                    break
                   }
                else{m <- m+1}
            }

              if(gamma^m < tol^2)
                break
        }

        if(gamma^m < tol^2)
                break
        
        
        it <- it+1
        conv[it+1] <- nloglik_cur;
        
    }


    return(list(betahat = c(betahat), sigmahat = sigmahat, alphahat = alphahat, conv = conv[1:(it+1)]))

}


### 
