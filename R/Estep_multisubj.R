## depend on metropolis_log_ratio_ts_fix_eta_r()

Estep_multisubj <- function(multi_y, x_vec, n_subj, theta_k, H0,
                            start.lst, eta = 1,
                            n.mcmc = 1000, burn = 100, thin = 1, 
                            name.par, adapt = TRUE, tune.lst, keep.lst,
                            tune.len = 50, target.accept.rate = 0.35, 
                            a = NULL, b = NULL, v1, v2, a1, b1, b2, xi) {
    
    
    #-------------------------------
    # Functions used in the algorithm
    #-------------------------------
    get.tune <- function(tune, keep, k, target = target.accept.rate){  
        # adaptive tuning
        a <- min(0.025, 1 / sqrt(k))
        exp(ifelse(keep < target, log(tune) - a, log(tune) + a))
    }
    
    #-------------------------------
    # Adaptive tuning
    #------------------------------- 
    Tb <- tune.len  # frequency of adaptive tuning
    # for t
    keep <- keep.lst
    keep.tmp <- keep  # track MH accpetance rate for adaptive tuning
    
    #-------------------------------
    # Here create some objects that are used later.
    #-------------------------------
    sig <- theta_k[1]
    tau <- theta_k[2]
    h <- theta_k[3]
    
    Kff <- se_ker(H0 = H0, tau = tau, h = h)
    
    
    #-------------------------------
    # Storage
    #-------------------------------
    sampleidx <- seq(from = (burn + thin), to = n.mcmc, by = thin)
    draws <- matrix(NA, nrow = length(sampleidx), ncol = length(name.par))
    colnames(draws) <- name.par
    
    
    #-------------------------------
    # Starting values
    #-------------------------------
    t_vec <- start.lst$t
    r1 <- start.lst$r1
    r2 <- start.lst$r2
    
    
    # if(!exists("Sigma_lst")) {
    #     Sigma_lst <- lapply(1:n_subj, function(s) {
    #         Kdf <- computeCovDer1(idx1 = t_vec[s], idx2 = x_vec, tau = tau, 
    #                               h = h) 
    #         
    #         Kdd <- computeCovDer2(idx1 = t_vec[s], tau = tau, h = h)
    #         Sigma <- Kff - quad.form.inv(Kdd, Kdf)
    #         return(Sigma)
    #     })
    # }
    
    
    #-------------------------------
    # M-H algorithm
    #------------------------------- 
    for (i in 1:n.mcmc) {
        if (i %% 5 == 0) cat("iter:", i, "\r")
        flush.console()
        
        # ====================================================================
        #------------------------------------------------------------------
        # Sampling distributuion of t given theta^(i)
        #------------------------------------------------------------------
        # ====================================================================
        
        Sigma_lst <- lapply(1:n_subj, function(s) {
            Kdf <- computeCovDer1(idx1 = t_vec[s], idx2 = x_vec, tau = tau, 
                                  h = h) 
            
            Kdd <- computeCovDer2(idx1 = t_vec[s], tau = tau, h = h)
            Sigma <- Kff - quad.form.inv(Kdd, Kdf)
            return(Sigma)
        })
        
        for (s in 1:n_subj) {
            # print(s)
            ts_star <- stats::runif(1, min = a, max = b)
            
            Kdf_new <- computeCovDer1(idx1 = ts_star, idx2 = x_vec, tau = tau, 
                                      h = h) 
            
            Kdd_new <- computeCovDer2(idx1 = ts_star, tau = tau, h = h)
            
            Sigma_new <- Kff - quad.form.inv(Kdd_new, Kdf_new)
            
            log_rho <- metropolis_log_ratio_ts_multisubj(ys = multi_y[, s],
                                                         t_vec = t_vec, 
                                                         ts_new = ts_star,
                                                         Sigma_s = Sigma_lst[[s]], 
                                                         Sigma_s_new = Sigma_new,
                                                         theta = theta_k,
                                                         r1 = r1, r2 = r2, 
                                                         s = s, eta = eta,
                                                         a1 = a1, b1 = b1,
                                                         b2 = b2)
            if (log(runif(1)) < log_rho) {
                t_vec[s] <- ts_star
                # w_vec[s] <- ws_star
                keep[[s]] <- keep[[s]] + 1
                keep.tmp[[s]] <- keep.tmp[[s]] + 1
                # Kdf <- computeCovDer1(idx1 = t_vec[s], idx2 = x_vec, tau = tau, 
                #                       h = h) 
                # Kdd <- computeCovDer2(idx1 = t_vec[s], tau = tau, h = h)
                # Sigma_lst[[s]] <- Kff - quad.form.inv(Kdd, Kdf)
            }
        }
        
        # ====================================================================
        # Sampling distributuion of r given t and theta^(i)
        # ====================================================================
        # r <- update_r(t_vec = t_vec, n_subj = n_subj, eta = eta)
        idx_1 <- t_vec < b1
        S1 <- sum(idx_1)
        # print(paste("S1 = ", S1))
        t_vec_1 <- t_vec[idx_1]
        t_vec_2 <- t_vec[!idx_1]
        S2 <- n_subj - S1
        # print(paste("S2 = ", S2))
        xi_1 <- 1 / xi ^ 2
        c1 <- 1 / ((S1 / eta ^ 2) + xi_1)
        # print(paste("c1 = ", c1))
        m1 <- c1 * ((sum(t_vec_1) / eta ^ 2) + (v1 / xi ^ 2))
        # print(paste("m1 = ", m1))
        c2 <- 1 / ((S2 / eta ^ 2) + xi_1)
        # print(paste("c2 = ", c2))
        m2 <- c2 * ((sum(t_vec_2) / eta ^ 2) + (v2 / xi ^ 2))
        # print(paste("m2 = ", m2))
        
        r1 <- rnorm(1, mean = m1, sd = sqrt(c1))
        r2 <- rnorm(1, mean = m2, sd = sqrt(c2))
        
        #  Save samples 
        # -----------------------------------------------------
        if (i > burn) {
            if (i %% thin == 0) {
                draws[(i - burn) %/% thin, ] <- c(t_vec, r1, r2)
            } 
        }
    }
    
    # Acceptance Probability
    #----------------------------
    keep <- lapply(keep, function(x) x / n.mcmc)
    # keep.r <- keep.r / n.mcmc
    # Write output
    #--------------
    return(list(draws = draws, 
                accept = keep,
                start = start.lst, 
                tune = tune.lst, 
                burn = burn, 
                thin = thin, 
                n.mcmc = n.mcmc, 
                sampleidx = sampleidx))
}