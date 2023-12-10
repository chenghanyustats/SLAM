#-------------------------------
# Functions used in the algorithm
#-------------------------------
get.tune <- function(tune, keep, k, target = target.accept.rate){
    # adaptive tuning
    a <- min(0.25, 1 / sqrt(k))
    exp(ifelse(keep < target, log(tune) - a, log(tune) + a))
}