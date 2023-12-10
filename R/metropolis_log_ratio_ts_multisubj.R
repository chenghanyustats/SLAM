library(truncnorm)
library(mvnfast)
metropolis_log_ratio_ts_multisubj <- function(ys, t_vec, w_vec = NULL, 
                                              ts_new, ws_new = NULL, 
                                              Sigma_s, Sigma_s_new,
                                              theta, r1, r2, s, eta, 
                                              a1, b1, b2) {
    # conditional of ts (Metropolis ratio)
    # ============================================
    n <- length(ys)
    # n_subj <- length(t_vec)
    sig_diag <- diag(theta[1] ^ 2, n)
    mu <- rep(0L, n)
    
    log_normal <- dmvn(ys, mu = mu, sigma = sig_diag + Sigma_s, log = TRUE)
    
    log_normal_new <- dmvn(ys, mu = mu, sigma = sig_diag + Sigma_s_new, 
                           log = TRUE)
    
    if (ts_new < b1) {
        r <- r1
    } else {
        r <- r2
    }
    # log_normal_t_new <- dnorm(ts_new, r, eta, log = TRUE)
    # log_normal_t <- dnorm(t_vec[s], r, eta, log = TRUE)
    
    
    log_trunc_normal_t_new <- log(dtruncnorm(ts_new, a = a1, b = b2, 
                                             mean = r, sd = eta))
    log_trunc_normal_t <- log(dtruncnorm(t_vec[s], a = a1, b = b2, 
                                         mean = r, sd = eta))
    
    # return(log_normal_new - log_normal + log_normal_t_new - log_normal_t)
    return(log_normal_new - log_normal + log_trunc_normal_t_new - log_trunc_normal_t)
}