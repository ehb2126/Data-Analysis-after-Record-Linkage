### pred: numeric(n) --- predicted values
### resp: numeric(n) --- observed response
### blockindex: numeric(n) --- block index based on matching variables

optimal_matching <- function(pred, resp, blockindex){

    n <- length(pred)
    index <- 1:n
    if(length(resp) != n)
        stop("Error: length of 'pred' has to match length of 'resp' \n")

    if(length(blockindex) != n)
        stop("Error: length of 'blockindex' has to match length of 'resp' \n")

    xy <- data.frame(x=pred, y=resp, ix=index)
    list_xy <- split(xy, blockindex)
    Llist <- lapply(list_xy, function(z) {ox <- order(z$x); oy <- order(z$y); return(cbind(z$ix[oy], z$ix[ox]))})

    L <- cbind(unlist(lapply(Llist, function(z) z[,1])), unlist(lapply(Llist, function(z) z[,2]))) 

    pihat <- L[order(L[,1]),2]

    return(pihat)

}
