### X:       --- n x p design matrix
### y:       --- n x 1 observed response
### lambda:  --- non-negative tuning parameter (proportional to noise level)
### control: --- additional parameters

robust_regression <- function(X, y, lambda, control = list(tol = 1E-4, maxit = 1E4)){

    tol <- control$tol
    maxit <- control$maxit

    
    X <- as.matrix(X)
    n <- nrow(X)
    p <- ncol(X)
    
    y <- as.numeric(y)

    if(length(y) != n)
        stop("Error: Length of 'y' must be equal to the number of rows of 'X' \n")

    lambda <- lambda[1]

    if(lambda[1] < 0)
        stop("Error: 'lambda' must be non-negative \n")


    ### intialization + compute QR decomposition of X
    
    qrX <- qr(X)

    yhat  <- qr.fitted(qrX, y)
    xi <- numeric(n)
    conv <- numeric(maxit)
    
    it <- 1
    conv[it] <- sum((y-yhat)^2)
   
    ### iterations

    while(TRUE){

        xi <- soft(y - yhat, lambda)
        z <- y - xi
        yhat <- qr.fitted(qrX, z)
        r <- y - yhat - xi
        it <- it+1

        conv[it] = sum(r^2) + 2*lambda*sum(abs(xi))
   
        if(conv[it-1] - conv[it] < tol)
            break
  
   
        if(it >= maxit)
            break 

    }

    conv <- conv[1:it]
    betahat <- qr.coef(qrX, z)

    return(list(betahat = betahat, xi = xi, conv = conv))

}


### soft thresholding function
soft <- function(z, t){

    zshrink <- (abs(z) > t) * (z - sign(z)*t)

    zshrink
    

}





