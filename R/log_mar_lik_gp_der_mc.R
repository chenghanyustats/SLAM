## original marg_lik_gp_der_mc_new_h_1()


# log_mar_lik_gp_der_mc <- function(theta, y, x_vec, der_mc_mat, H0,
#                                        ga_shape = 6, ga_rate = 5,
#                                        is.sig.par = TRUE) {
#     nn <- ncol(der_mc_mat)
#     n <- length(y)
#     mu_vec <- rep(0L, n)
# 
#     Kff <- power_expo_h(H0, theta[2], theta[3])
#     
#     Kdf_mc <- apply(der_mc_mat, 1, computeCovDer1_h, idx2 = x_vec, 
#                      tau = theta[2], h = theta[3])
#     
#     den_value_mc <- apply(Kdf_mc, 2, function(x) {
#         Sigma <- Kff - tcrossprod(x) / ((theta[2] / theta[3]) ^ 2)
#         mvnfast::dmvn(y, mu = mu_vec, sigma = Sigma + diag(theta[1] ^ 2, n))
#     })
#     A <- sum(den_value_mc) / length(den_value_mc)
# 
#     if (is.sig.par) {
#         B <- invgamma::dinvgamma(theta[1] ^ 2, shape = ga_shape, 
#                                  rate = ga_rate, log = TRUE)
#         return(-log(A) - B)
#     } else {
#         return(-log(A))
#     }
# }

marg_lik_gp_der_mc_new_h_1 <- function(theta, y, x_vec, der_mc_mat, H0,
                                       ga_shape = 6, ga_rate = 5,
                                       is.sig.par = TRUE) {
    # dim_der <- dim(der_mc_mat)
    nn <- ncol(der_mc_mat)
    n <- length(y)
    # mu_vec <- rep(0, n)
    Kff <- se_ker(H0 = H0, tau = theta[2], h = theta[3])
    
    Kdf_lst <- as.list(data.frame(apply(der_mc_mat, 2, computeCovDer1_h, 
                                        idx2 = x_vec, 
                                        tau = theta[2], h = theta[3])))
    Kdf_lst_old <- lapply(Kdf_lst, function(x){
        matrix(x, nrow = nn, ncol = n)
    })
    
    # den_value_lst <- lapply(1:mm, function(m) {
    #     Sigma <- Kff - crossprod(Kdf_lst[[m]]) / ((theta[2] / theta[3]) ^ 2)
    #     diag(Sigma) <- diag(Sigma) + theta[1] ^ 2
    #     dmvn(y, mu = mu_vec, sigma = Sigma)
    # })
    den_value_lst_old <- lapply(Kdf_lst_old, function(m) {
        Sigma <- Kff_old - crossprod(m) / ((theta[2] / theta[3]) ^ 2)
        diag(Sigma) <- diag(Sigma) + theta[1] ^ 2
        dmvn(y, mu = rep(0, n), sigma = Sigma)
        # Sigma
    })
    # , log = TRUE
    A <- sum(unlist(den_value_lst)) / nrow(der_mc_mat)
    if (is.sig.par) {
        B <- dinvgamma(theta[1] ^ 2, shape = ga_shape, rate = ga_rate, log = TRUE)
        return(-log(A)-B)
    } else {
        return(-log(A))
    }
    # B <- dinvgamma(theta[1]^2, shape = 3, scale = 1, log = TRUE)
    # return(-log(A))
    # return(-A)
}


log_mar_lik_gp_der_mc <- function(theta, y, x_vec, der_mc_mat, H0,
                                  ga_shape = 6, ga_rate = 5,
                                  is.sig.par = TRUE) {
    nn <- nrow(der_mc_mat)
    n <- length(y)
    # mu_vec <- rep(0L, n)
    
    Kff <- se_ker(H0 = H0, tau = theta[2], h = theta[3])
    # Kff <- compute_cov_1d(idx1 = x_vec, tau = theta[2], h = theta[3])
    
    # Kdf_mc <- apply(der_mc_mat, 2, computeCovDer1, idx2 = x_vec, 
    #                 tau = theta[2], h = theta[3])
    
    
    Kdf_lst <- as.list(data.frame(apply(der_mc_mat, 2, computeCovDer1, 
                                        idx2 = x_vec, 
                                        tau = theta[2], h = theta[3])))
    
    Kdf_lst_new <- lapply(Kdf_lst, function(x){
        matrix(x, nrow = nn, ncol = n)
    })
    
    den_value_lst_new <- lapply(Kdf_lst_new, function(m) {
        # m <- matrix(m, nrow = nn, ncol = n)
        Sigma <- Kff - crossprod(m) / ((theta[2] / theta[3]) ^ 2)
        Psi <- Sigma + diag(theta[1] ^ 2, n)
        mvnfast::dmvn(y, mu = rep(0L, n), sigma = Psi)
        # Sigma <- Sigma + diag(theta[1] ^ 2, n)
        # Sigma
    })
    
    
    # den_value_mc <- apply(Kdf_mc, 2, function(x) {
    #     x <- matrix(x, nrow = nn, ncol = n)
    #     Sigma <- Kff - crossprod(x) / ((theta[2] / theta[3]) ^ 2)
    #     mvnfast::dmvn(y, mu = mu_vec, sigma = Sigma + diag(theta[1] ^ 2, n))
    # })
    
    # A <- sum(den_value_mc) / length(den_value_mc)
    
    # A <- sum(unlist(den_value_lst)) / ncol(der_mc_mat)
    if (is.sig.par) {
        B <- invgamma::dinvgamma(theta[1] ^ 2, shape = ga_shape, 
                                 rate = ga_rate, log = TRUE)
        return(-log(sum(unlist(den_value_lst_new)) / ncol(der_mc_mat)) - B)
    } else {
        return(-log(sum(unlist(den_value_lst_new)) / ncol(der_mc_mat)))
    }
    
    
    # if (is.sig.par) {
    #     # B <- invgamma::dinvgamma(theta[1] ^ 2, shape = ga_shape, 
    #     #                          rate = ga_rate, log = TRUE)
    #     return(-log(sum(den_value_mc) / length(den_value_mc)) - 
    #                invgamma::dinvgamma(theta[1] ^ 2, shape = ga_shape, 
    #                                    rate = ga_rate, log = TRUE))
    # } else {
    #     return(-log(sum(den_value_mc) / length(den_value_mc)))
    # }
}



log_mar_lik_gp_der_mc <- function(theta, y, x_vec, der_mc_mat, H0,
                                  ga_shape = 6, ga_rate = 5,
                                  a_h = 1, b_h = 1,
                                  is.sig.par = TRUE,
                                  is.h.par = TRUE) {
    nn <- nrow(der_mc_mat)
    n <- length(y)
    # mu_vec <- rep(0L, n)
    
    Kff <- se_ker(H0 = H0, tau = theta[2], h = theta[3])
    # Kff <- compute_cov_1d(idx1 = x_vec, tau = theta[2], h = theta[3])
    
    # Kdf_mc <- apply(der_mc_mat, 2, computeCovDer1, idx2 = x_vec, 
    #                 tau = theta[2], h = theta[3])
    
    
    Kdf_lst <- as.list(data.frame(apply(der_mc_mat, 2, computeCovDer1, 
                                        idx2 = x_vec, 
                                        tau = theta[2], h = theta[3])))
    
    Kdf_lst_new <- lapply(Kdf_lst, function(x){
        matrix(x, nrow = nn, ncol = n)
    })
    
    den_value_lst_new <- lapply(Kdf_lst_new, function(m) {
        # m <- matrix(m, nrow = nn, ncol = n)
        Sigma <- Kff - crossprod(m) / ((theta[2] / theta[3]) ^ 2)
        Psi <- Sigma + diag(theta[1] ^ 2, n)
        mvnfast::dmvn(y, mu = rep(0L, n), sigma = Psi)
        # Sigma <- Sigma + diag(theta[1] ^ 2, n)
        # Sigma
    })
    
    
    # den_value_mc <- apply(Kdf_mc, 2, function(x) {
    #     x <- matrix(x, nrow = nn, ncol = n)
    #     Sigma <- Kff - crossprod(x) / ((theta[2] / theta[3]) ^ 2)
    #     mvnfast::dmvn(y, mu = mu_vec, sigma = Sigma + diag(theta[1] ^ 2, n))
    # })
    
    # A <- sum(den_value_mc) / length(den_value_mc)
    
    # A <- sum(unlist(den_value_lst)) / ncol(der_mc_mat)
    if (is.sig.par) {
        B <- invgamma::dinvgamma(theta[1] ^ 2, shape = ga_shape, 
                                 rate = ga_rate, log = TRUE) + log(2 * theta[1])
        if (is.h.par) {
            C <- dgamma(theta[3], shape = a_h, rate = b_h, log = TRUE)
            return(-log(sum(unlist(den_value_lst_new)) / ncol(der_mc_mat)) - B - C)
        } else {
            return(-log(sum(unlist(den_value_lst_new)) / ncol(der_mc_mat)) - B)
        }
    } else {
        if (is.h.par) {
            C <- dgamma(theta[3], shape = a_h, rate = b_h, log = TRUE)
            return(-log(sum(unlist(den_value_lst_new)) / ncol(der_mc_mat)) - C)
        } else {
            return(-log(sum(unlist(den_value_lst_new)) / ncol(der_mc_mat)))
        }
        # return(-log(sum(unlist(den_value_lst_new)) / ncol(der_mc_mat)))
    }
}




# Reduce("+", your.list) 
# test_val1 <- 1.4 + rnorm(100, 0, 0.1)
# test_val2 <- 1.4 + rnorm(200, 0, 0.1)
# 
# test_mat1 <- matrix(test_val1, 1, 100)
# test_mat2 <- matrix(test_val2, 2, 100)
# 
# 
# 
# test <- apply(test_mat2, 2, computeCovDer1_h,
#       idx2 = YY[[1]]$x, tau = 1, h = 1)
# 
# # test_lst0 <- as.list(test)
# 
# 
# test_lst <- (as.list(data.frame(apply(test_mat2, 2, computeCovDer1_h,
#                                       idx2 = YY[[1]]$x, tau = 1, h = 1))))
# 
# # test_lst1 <- lapply(test_lst, function(x) matrix(x, nrow = 1, ncol = 50))
# 
# 
# microbenchmark(test <- apply(test_mat2, 2, computeCovDer1_h,
#                              idx2 = YY[[1]]$x, tau = 1, h = 1),
#                test_lst <- (as.list(data.frame(apply(t(test_mat2), 1, computeCovDer1_h,
#                                                      idx2 = YY[[1]]$x, tau = 1, h = 1)))), times = 100)
# # 
# # test_apply <- apply(test, 2, function(x) {
#     x <- matrix(x, nrow = 2, ncol = 50)
#     Sigma <- Kff - crossprod(x)
#     mu_vec <- rep(0L, 50)
#     Psi <- Sigma + diag(1, 50)
#     mvnfast::dmvn(y, mu = mu_vec, sigma = Psi)
# })
# # 
# microbenchmark(test_apply <- apply(test, 2, function(x) {
#     x <- matrix(x, nrow = 2, ncol = 50)
#     Sigma <- Kff - crossprod(x)
#     mu_vec <- rep(0L, 50)
#     Psi <- Sigma + diag(1, 50)
#     # mvnfast::dmvn(y, mu = mu_vec, sigma = Psi)
# }),
# den_value_lst <- lapply(test_lst, function(m) {
#     m <- matrix(m, nrow = 2, ncol = 50)
#     Sigma <- Kff - crossprod(m)
#     mu_vec <- rep(0L, 50)
#     Psi <- Sigma + diag(1, 50)
#     # mvnfast::dmvn(y, mu = mu_vec, sigma = Psi)
# }), times = 10)


# microbenchmark({Kdf_lst <- as.list(data.frame(apply(der_mc_mat, 1, computeCovDer1_h, idx2 = x_vec, 
#                                                    tau = theta[2], h = theta[3])))
#                Kdf_lst <- lapply(Kdf_lst, function(x){
#                    matrix(x, nrow = nn, ncol = n)
#                })
#                den_value_lst <- lapply(Kdf_lst, function(m) {
#                    Sigma <- Kff - crossprod(m) / ((theta[2] / theta[3]) ^ 2)
#                    diag(Sigma) <- diag(Sigma) + theta[1] ^ 2
#                    # dmvn(y, mu = mu_vec, sigma = Sigma)
#                })},
#                {Kdf_mc <- apply(der_mc_mat, 2, computeCovDer1, idx2 = x_vec, 
#                                 tau = theta[2], h = theta[3])
#                den_value_mc <- apply(Kdf_mc, 2, function(x) {
#                    x <- matrix(x, nrow = nn, ncol = n)
#                    Sigma <- Kff - crossprod(x) / ((theta[2] / theta[3]) ^ 2)
#                    Sigma <- Sigma + diag(theta[1] ^ 2, n)
#                    # mvnfast::dmvn(y, mu = mu_vec, sigma = Sigma + diag(theta[1] ^ 2, n))
#                })},
#                times = 100)



# Kdf_lst <- as.list(data.frame(apply(der_mc_mat, 1, computeCovDer1_h, idx2 = x_vec, 
#                                     tau = theta[2], h = theta[3])))
# Kdf_lst <- lapply(Kdf_lst, function(x){
#     matrix(x, nrow = nn, ncol = n)
# })
# den_value_lst <- lapply(Kdf_lst, function(m) {
#     Sigma <- Kff - crossprod(m) / ((theta[2] / theta[3]) ^ 2)
#     diag(Sigma) <- diag(Sigma) + theta[1] ^ 2
#     # dmvn(y, mu = mu_vec, sigma = Sigma)
# })




# microbenchmark(sum(unlist(den_value_lst)) / 100, 
#                sum(test_apply) / 100)
# 
# 
# 
# den_value_lst <- lapply(test_lst, function(m) {
#     m <- matrix(m, nrow = 1, ncol = 50)
#     Sigma <- Kff - crossprod(m)
#     mu_vec <- rep(0L, 50)
#     Psi <- Sigma + diag(1, 50)
#     mvnfast::dmvn(y, mu = mu_vec, sigma = Psi)
# })
# 
# 
# microbenchmark(test_apply <- apply(test, 2, function(x) {
#     x <- matrix(x, nrow = 1, ncol = 50)
#     Sigma <- Kff - crossprod(x)
#     mu_vec <- rep(0L, 50)
#     Psi <- Sigma + diag(1, 50)
#     mvnfast::dmvn(y, mu = mu_vec, sigma = Psi)
# }),
# test_apply <- apply(test, 2, function(x) {
#     # x <- matrix(x, nrow = 1, ncol = 50)
#     Sigma <- Kff - tcrossprod(x)
#     mu_vec <- rep(0L, 50)
#     Psi <- Sigma + diag(1, 50)
#     mvnfast::dmvn(y, mu = mu_vec, sigma = Psi)
# }), times = 1000)
# 
# 
# microbenchmark(mean(test_apply), sum(test_apply) / length(100))


