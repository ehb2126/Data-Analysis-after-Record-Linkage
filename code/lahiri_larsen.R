### X: --- n x p design matrix
### y: --- n x 1 observed response
### Q: --- n x n adjustment matrix as in Lahiri and Larsen (2005)

lahiri_larsen <- function(X, y, Q){

    X <- as.matrix(X)
    n <- nrow(X)
    p <- ncol(X)
    
    y <- as.numeric(y)
    Q <- as.matrix(Q)

    if(length(y) != n)
        stop("Error: Length of 'y' must be equal to the number of rows of 'X' \n")
    if(nrow(Q) != ncol(Q))
        stop("Error: 'Q' must be a square matrix \n")
    if(nrow(Q) != n)
        stop("Dimension of 'Q' must matcht the number of rows of 'X' \n")
    
    
    
    W <- Q %*% X
    
    betahat <- solve(crossprod(W), t(W) %*% y)
   

    return(betahat)

}

### generates the adjustment matrix Q based on blocking variables; each
### unique combination of blocking variables corresponds to an index "blockindex"

generate_Q_block <- function(blockindex){

    n <- length(blockindex)
    L <- split(1:n, blockindex)
    Q <- matrix(data = 0, nrow = n, ncol = n)
    lapply(L, function(z) Q[z,z] <<- 1/length(z))

    return(Q)

}

# fast versions Lahiri-Larsen in the blocking case
lahiri_larsen_block <- function(X, y, blockindex){

    listX <- lapply(split(data.frame(X = X), blockindex), colMeans)
    
    spY <- split(y, blockindex)
    listY <- lapply(spY, mean)
    ws <- unlist(lapply(spY, length))
    ybar <- unlist(listY)
    Xbar <- do.call(rbind, listX)

    return(lm(ybar ~ Xbar - 1, weights = ws))
    
}
