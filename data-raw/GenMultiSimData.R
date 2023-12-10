################################################################################
# DGP multisubject and multigroup simulated data generation                    #                                                              #
# "./data-raw/GenMultiSimData.R"                                               #
# Cheng-Han Yu                                                                 #
################################################################################
# library(here)
# source("./data-raw/seed.R")

## data list
YYsin <- list()
YYcos <- list()


## size
n <- 100

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

for (k in 1:no_data) {
    set.seed(seeds[k])
    y <- regfcn_cos(x, k) + rnorm(n, 0, sig)
    YYcos[[k]] <- list(x = x, y = y)
    y <- regfcn_sin(x, k) + rnorm(n, 0, sig)
    YYsin[[k]] <- list(x = x, y = y)
}



## plotting
par(mfrow = c(1, 1), mar = c(4, 4, 2, 0))

## reg function

idx <- seq(x_a, x_b, length.out = 500)
truefcn_sin <- regfcn_sin(idx, k = 1)
par(mar = c(3, 3.5, 2, 0), mgp = c(2, 0.8, 0))
par(mfrow = c(1, 2))
plot(idx, truefcn_sin, col = "red", type = "l", xlim = c(x_a - 0.05, x_b + 0.05),
     lwd = 1, main = "Sine group 1", cex.main = 1, cex.lab = 1,
     ylim = c(-2.2, 2.2), xlab = "x", ylab = "f(x)")
# no_data <- 3
for (k in 2:no_data) {
    truefcn <- regfcn_sin(idx, k = k)
    # truefcn <- sqrt(l * idx * (1 - idx)) * 15 * cos(2 * pi + j / (m * k + idx + 1))
    lines(idx, truefcn, pch = 3, col = k)
}

truefcn_cos <- regfcn_cos(idx, k = 1)
plot(idx, truefcn_cos, col = "red", type = "l", cex.main = 1, cex.lab = 1,
     lwd = 1, main = "Cosine group 2", ylim = c(-3.2, 0.2),
     xlab = "x", ylab = "f(x)")

for (k in 2:no_data) {
    truefcn <- regfcn_cos(idx, k = k)
    lines(idx, truefcn, pch = 3, col = k)
}


## data
# par(mfrow = c(2, 5), mar = c(4, 4, 2, 0))
# for (k in 1:no_data) {
#     truefcn <- regfcn_cos(idx, k = k)
#     plot(idx, truefcn, col = "red", type = "l", xlim = c(x_a - 0.05, x_b + 0.05),
#          lwd = 2, main = paste("cos Subj", k),
#          ylim = c(-4, 3), las = 1, xlab = "x", ylab = "f(x)")
#     points(YYcos[[k]]$x, YYcos[[k]]$y, pch = 3, col = 1)
#     points(cri_pts_cos_lst[[k]], regfcn_cos(cri_pts_cos_lst[[k]], k = k),
#            pch = "_", lwd = .1,
#            col = "green", cex = 3)
# }
#
# for (k in 1:no_data) {
#     truefcn <- regfcn_sin(idx, k = k)
#     plot(idx, truefcn, col = "red", type = "l", xlim = c(x_a - 0.05, x_b + 0.05),
#          lwd = 2, main = paste("sin Subj", k),
#          ylim = c(-4, 3), las = 1, xlab = "x", ylab = "f(x)")
#     points(YYsin[[k]]$x, YYsin[[k]]$y, pch = 3, col = 1)
#     points(cri_pts_sin_lst[[k]], regfcn_sin(cri_pts_sin_lst[[k]], k = k),
#            pch = "_", lwd = .1,
#            col = "green", cex = 3)
# }


# save(YYcos, YYsin, sig, x_a, x_b, cri_pts_cos_lst, cri_pts_sin_lst, no_data,
#      file = "./data/multi_sim_data.RData")

# load("./data/multi_sim_data.RData")
