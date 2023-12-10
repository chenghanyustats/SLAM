#####
# matern kernel
#####
matern_ker <- function(x, y, nu = 1.5, tau = 1, l = 1) {
    if (!(nu %in% c(0.5, 1.5, 2.5))) {
        stop("p must be equal to 0.5, 1.5 or 2.5")
    }
    p <- nu - 0.5
    d <- abs(x - y)
    if (p == 0) {
        return(tau ^ 2 * exp(- d / l))
    } else if (p == 1) {
        b <- sqrt(3) * d / l
        return(tau ^ 2 * (1 + b) * exp(-b))
    } else {
        b <- sqrt(5) * d / l
        return(tau ^ 2 * (1 + b + b ^ 2 / 3) * exp(-b))
    } 
}