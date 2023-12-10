################################################################################
# DGP multisubject and multigroup simulated data generation                    #                                                              #
# "./data-raw/data_replicates.R"                                               #
# Cheng-Han Yu                                                                 #
################################################################################
# library(here)
source("./data-raw/seed.R")


## size
n <- 200

## replicates
n_rep <- 100

## number of replicate (individuals) assume both groups have the same number
no_data <- 10

## standard deviation of error noise
sig <- 0.25
# sig <- 0.3

## input range
x_a <- 0
x_b <- 1


## true reg function

regfcn_cos <- function(x, k) {
    # reg <- sqrt(x * (1 - x)) * 12 * cos(2 * pi + 2 / (-0.035 * k + x + 1))
    reg <- cos(2*pi*x + k/10 + 1.2) - 3*x
    return(reg)
}

regfcn_sin <- function(x, k) {
    reg <- -2 * sin(2 * pi * x + k/15 - 0.3)
    return(reg)
}



## first derivative
fn_sin <- function(x, k) {
    -4 * pi * cos(2 * pi * x + k/15 - 0.3)
}

fn_cos <- function(x, k) {
    # 6 * (1 - 2 * x) * cos(2 / (x - 7*k/200 + 1)) / (sqrt(x * (1 - x))) +
    #     24 * sqrt(x * (1 - x)) * sin(2 / (x - 7*k/200 + 1)) / ((x - 7*k/200 + 1)^2)
    - 2 * pi * sin(2*pi*x + k/10+ 6/5) - 3
}

## stationary point
cri_pts_cos_lst <- vector("list", 10)
cri_pts_sin_lst <- vector("list", 10)

for (k in 1:10) {
    cri_pts_cos <- c(uniroot(fn_cos, c(x_a, 0.5), k = k)$root,
                     uniroot(fn_cos, c(0.5, x_b), k = k)$root)
    cri_pts_sin <- c(uniroot(fn_sin, c(x_a, 0.5), k = k)$root,
                     uniroot(fn_sin, c(0.5, x_b), k = k)$root)
    cri_pts_cos_lst[[k]] <- cri_pts_cos
    cri_pts_sin_lst[[k]] <- cri_pts_sin
}

extreme_val_cos_lst <- vector("list", 10)

for (k in 1:10) {
    cri_pts_cos <- c(uniroot(fn_cos, c(x_a, 0.5), k = k)$root,
                     uniroot(fn_cos, c(0.5, x_b), k = k)$root)

    extreme_val_cos_lst[[k]] <- regfcn_cos(cri_pts_cos, k)
}


# cri_pts_sin_lst <- vector("list", 10)


## produce data
x <- seq(x_a, x_b, length.out = n)
YYsin_lst <- vector("list", length = n_rep)
YYsin_lst <- lapply(YYsin_lst, function(x) {
    vector(mode = "list", length = no_data)
})

YYcos_lst <- vector("list", length = n_rep)
YYcos_lst <- lapply(YYcos_lst, function(x) {
    vector(mode = "list", length = no_data)
})

for (j in 1: n_rep) {
    for (k in 1:no_data) {
        set.seed(seeds[k])
        y <- regfcn_cos(x, k) + rnorm(n, 0, sig)
        YYcos_lst[[j]][[k]] <- list(x = x, y = y)
        set.seed(seeds[k])
        y <- regfcn_sin(x, k) + rnorm(n, 0, sig)
        YYsin_lst[[j]][[k]] <- list(x = x, y = y)
    }
}


# save(YYcos_lst, YYsin_lst, sig, x_a, x_b, cri_pts_cos_lst, cri_pts_sin_lst,
#      no_data, n_rep,
#      file = "./data/data_replicates.RData")
