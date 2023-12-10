##############################################################################
# DGP multisubject and multigroup simulated data generation                    #                                                              #
# "./data-raw/data_from_model.R"                                               #
# Cheng-Han Yu                                                                 #
################################################################################

# library(here)
source("./data-raw/seed.R")
library(extraDistr)


## size
n <- 100

## replicates
n_rep <- 50

## parameters
sig <- 0.5
beta0_1 <- 0.3
beta1_1 <- -0.5
beta0_2 <- -0.3
beta1_2 <- 1

logistic <- function(x) {
    1 / (1 + exp(-x))
}

r1_1 <- logistic(beta0_1) ##young 1
r1_2 <- logistic(beta0_2) ##young 2
r2_1 <- logistic(beta0_1 + beta1_1) ##old 1
r2_2 <- logistic(beta0_2 + beta1_2) ##old 2


a_1 <- 0
b_1 <- 0.5
a_2 <- 0.5
b_2 <- 1


r_loc <- function(r, a = 0, b = 0.5) {
    (b - a) * r + a
}


r_loc(r1_1)
r_loc(r2_1)
r_loc(r1_2, a = a_2, b = b_2)
r_loc(r2_2, a = a_2, b = b_2)

eta <- 8


## data generating
no_data <- 10

t_sim_lst <- vector("list", length = n_rep)

t_sim_lst <- lapply(t_sim_lst, function(x) {
    y <- vector(mode = "list", length = 4)
    names(y) <- c("t1_1", "t1_2", "t2_1", "t2_2")
    return(y)
})


for (j in 1:n_rep) {
    set.seed(seeds[j])
    t_sim_lst[[j]]$t1_1 <- (b_1 - a_1) * rbeta(no_data, r1_1 * eta, (1 - r1_1) * eta) + a_1
    set.seed(seeds[j])
    t_sim_lst[[j]]$t1_2 <- (b_2 - a_2) * rbeta(no_data, r1_2 * eta, (1 - r1_2) * eta) + a_2
    set.seed(seeds[j])
    t_sim_lst[[j]]$t2_1 <- (b_1 - a_1) * rbeta(no_data, r2_1 * eta, (1 - r2_1) * eta) + a_1
    set.seed(seeds[j])
    t_sim_lst[[j]]$t2_2 <- (b_2 - a_2) * rbeta(no_data, r2_2 * eta, (1 - r2_2) * eta) + a_2
}

f1 <- function(x1, t1, amp = 20) {
    amp * (x1 - t1) ^ 2
}

f2 <- function(x2, t2, x1max, t1, amp = 20) {
    amp * (-(x2 - t2) ^ 2 + (x1max - t1) ^ 2 + (min(x2) - t2) ^ 2)
}



x1 <- seq(0, 0.5, by = 0.01)
x2 <- seq(0.5, 1, by = 0.01)

par(mfrow = c(1, 2))
par(mar = c(3.2, 3.2, 1.5, 0.5), mgp = c(2, 0.8, 0))
plot(c(x1,x2),
     c(f1(x1, t_sim_lst[[1]]$t1_1[1]),
       f2(x2, t_sim_lst[[1]]$t1_2[1], x1[length(x1)], t_sim_lst[[1]]$t1_1[1])),
     type = "l", col = 1, main = "group 1", ylim = c(-2, 6), las = 1,
     xlab = "x", ylab = "f(x)")
for (i in 2:no_data) {
    lines(c(x1,x2),
          c(f1(x1, t_sim_lst[[1]]$t1_1[i]),
            f2(x2, t_sim_lst[[1]]$t1_2[i],
               x1[length(x1)], t_sim_lst[[1]]$t1_1[i])),
          col = i)
}
abline(v = c(0.5, r_loc(r1_1), r_loc(r1_2, a = a_2, b = b_2)), lty = 2,
       col = c(1, 2, 2))


plot(c(x1,x2),
     c(f1(x1, t_sim_lst[[1]]$t2_1[1]),
       f2(x2, t_sim_lst[[1]]$t2_2[1], x1[length(x1)], t_sim_lst[[1]]$t2_1[1])),
     type = "l", col = 1, main = "group 2", ylim = c(-0.5, 6), las = 1,
     xlab = "x", ylab = "f(x)")
for (i in 2:no_data) {
    lines(c(x1,x2),
          c(f1(x1, t_sim_lst[[1]]$t2_1[i]),
            f2(x2, t_sim_lst[[1]]$t2_2[i], x1[length(x1)], t_sim_lst[[1]]$t2_1[i])),
          col = i)
}
abline(v = c(0.5, r_loc(r2_1), r_loc(r2_2, a = a_2, b = b_2)), lty = 2,
       col = c(1, 2, 2))


# produce data
# x <- seq(x_a, x_b, length.out = n)
x1 <- seq(0, 0.49, length = n/2)
x2 <- seq(0.51, 1, length = n/2)
xx <- c(x1, x2)


YY1_lst_la <- vector("list", length = n_rep)
YY1_lst_la <- lapply(YY1_lst_la, function(x) {
    vector(mode = "list", length = no_data)
})

YY2_lst_la <- vector("list", length = n_rep)
YY2_lst_la <- lapply(YY2_lst_la, function(x) {
    vector(mode = "list", length = no_data)
})


for (j in 1: n_rep) {
    for (k in 1:no_data) {
        set.seed(seeds[k])
        y <- c(f1(x1, t_sim_lst[[j]]$t1_1[k]),
               f2(x2, t_sim_lst[[j]]$t1_2[k], x1[length(x1)], t_sim_lst[[j]]$t1_1[k])) +
            extraDistr::rlaplace(n, mu = 0, sigma = sig)
        YY1_lst_la[[j]][[k]] <- list(x = xx, y = y)
        set.seed(seeds[k])
        y <- c(f1(x1, t_sim_lst[[j]]$t2_1[k]),
               f2(x2, t_sim_lst[[j]]$t2_2[k], x1[length(x1)], t_sim_lst[[j]]$t2_1[k])) +
            extraDistr::rlaplace(n, mu = 0, sigma = sig)
        YY2_lst_la[[j]][[k]] <- list(x = xx, y = y)
    }
}


par(mfrow = c(2, 2))
for (k in 1:no_data) {
    plot(c(x1,x2),
         c(f1(x1, t_sim_lst[[1]]$t1_1[k]),
           f2(x2, t_sim_lst[[1]]$t1_2[k], x1[length(x1)], t_sim_lst[[1]]$t1_1[k])),
         type = "l", col = 1, main = "group 1", ylim = c(-2, 6), las = 1,
         xlab = "x", ylab = "f(x)")
    points(YY1_lst_la[[1]][[k]]$x, YY1_lst_la[[1]][[k]]$y)
}


# save(YY1_lst_la, YY2_lst_la, xx, sig, beta0_1, beta1_1, beta0_2, beta1_2, no_data,
#      n_rep,
#      r1_1, r1_2, r2_1, r2_2, eta, a_1, b_1, a_2, b_2,
#      t1_1, t1_2, t2_1, t2_2,
#      f1, f2, t_sim_lst,
#      file = "./data/data_la.RData")
