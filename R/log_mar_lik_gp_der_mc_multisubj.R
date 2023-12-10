log_mar_lik_gp_der_mc_multisubj <- function(theta, multi_y, x_vec, der_mc_mat, H0,
                                            ga_shape = 6, ga_rate = 5,
                                            a_h = 1, b_h = 1,
                                            is.sig.par = TRUE,
                                            is.h.par = FALSE) {
    
    ## der_mc_mat: matrix of dim n.mcmc by n_subj if assume 1 t
    
    # nn <- nrow(der_mc_mat)
    D <- nrow(der_mc_mat)
    nn <- 1
    n <- nrow(multi_y)
    n_subj <- ncol(multi_y)
    mu <- rep(0L, n)
    Kff <- se_ker(H0 = H0, tau = theta[2], h = theta[3])
    
    log_all <- sapply(1:n_subj, function(s) {
        Kdf_lst <- as.list(data.frame(apply(matrix(der_mc_mat[, s], nn, D), 2, 
                                            computeCovDer1, 
                                            idx2 = x_vec, 
                                            tau = theta[2], h = theta[3])))
        Kdf_lst_new <- lapply(Kdf_lst, function(x) {
            matrix(x, nrow = nn, ncol = n)
        })
        
        den_value_lst_new <- lapply(Kdf_lst_new, function(m) {
            Sigma <- Kff - crossprod(m) / ((theta[2] / theta[3]) ^ 2)
            Psi <- Sigma + diag(theta[1] ^ 2, n)
            mvnfast::dmvn(multi_y[, s], mu = mu, sigma = Psi)
        })
        return(log(sum(unlist(den_value_lst_new)) / D))
    })
    
    log_summ <- sum(log_all)
    
    if (is.sig.par) {
        B <- invgamma::dinvgamma(theta[1] ^ 2, shape = ga_shape, 
                                 rate = ga_rate, log = TRUE) + log(2 * theta[1])
        if (is.h.par) {
            C <- dgamma(theta[3], shape = a_h, rate = b_h, log = TRUE)
            return(-log_summ - B - C)
        } else {
            return(-log_summ - B)
        }
    } else {
        if (is.h.par) {
            C <- dgamma(theta[3], shape = a_h, rate = b_h, log = TRUE)
            return(-log_summ - C)
        } else {
            return(-log_summ)
        }
        # return(-log_summ)
    }
}
