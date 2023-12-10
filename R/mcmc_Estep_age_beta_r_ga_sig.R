mcmc_Estep_age_beta_r_ga_sig <- function(multi_y_young, multi_y_old, x_vec, 
                                         n_subj_young, n_subj_old, theta_k, H0,
                                         start.lst,
                                         n.mcmc = 1000, burn = 100, thin = 1, 
                                         name.par, adapt = TRUE, 
                                         # tune.lst.young, keep.lst.young,
                                         # tune.lst.old, keep.lst.old,
                                         keep.young.eta.1, keep.young.eta.2,
                                         keep.old.eta.1, keep.old.eta.2,
                                         tune.young.eta.1, tune.young.eta.2,
                                         tune.old.eta.1, tune.old.eta.2,
                                         tune.len = 50, target.accept.rate = 0.35, 
                                         # a = NULL, b = NULL,
                                         a_young_1, b_young_1, a_young_2, b_young_2,
                                         a_old_1, b_old_1, a_old_2, b_old_2,
                                         a_1, b_1, a_2, b_2,
                                         alpha_r_young_1, beta_r_young_1,
                                         alpha_r_young_2, beta_r_young_2,
                                         alpha_r_old_1, beta_r_old_1,
                                         alpha_r_old_2, beta_r_old_2,
                                         a_eta_young_1, 
                                         b_eta_young_1 ,
                                         a_eta_young_2, 
                                         b_eta_young_2,
                                         a_eta_old_1, 
                                         b_eta_old_1,
                                         a_eta_old_2, 
                                         b_eta_old_2,
                                         ga_shape, ga_rate) {
    #-------------------------------
    # Functions used in the algorithm
    #-------------------------------
    get.tune <- function(tune, keep, k, target = target.accept.rate){  
        # adaptive tuning
        a <- min(0.25, 1 / sqrt(k))
        exp(ifelse(keep < target, log(tune) - a, log(tune) + a))
    }
    
    #-------------------------------
    # Adaptive tuning
    #------------------------------- 
    Tb <- tune.len  # frequency of adaptive tuning
    
    # for eta
    keep.tmp.young.eta.1 <- keep.young.eta.1  # track MH acceptance rate for adaptive tuning
    keep.tmp.young.eta.2 <- keep.young.eta.2
    keep.tmp.old.eta.1 <- keep.old.eta.1
    keep.tmp.old.eta.2 <- keep.old.eta.2
    
    #-------------------------------
    # Here create some objects that are used later.
    #-------------------------------
    tau <- theta_k[1]
    h <- theta_k[2]
    Kff <- se_ker(H0 = H0, tau = theta_k[1], h = theta_k[2])
    
    #-------------------------------
    # Storage
    #-------------------------------
    sampleidx <- seq(from = (burn + thin), to = n.mcmc, by = thin)
    draws <- matrix(NA, nrow = length(sampleidx), ncol = length(name.par))
    colnames(draws) <- name.par
    
    
    #-------------------------------
    # Starting values
    #-------------------------------
    t_mat_young <- start.lst$t_young
    t_mat_old <- start.lst$t_old
    r1_young <- start.lst$r1_young
    r2_young <- start.lst$r2_young
    r1_old <- start.lst$r1_old
    r2_old <- start.lst$r2_old
    eta1_young <- start.lst$eta1_young
    eta2_young <- start.lst$eta2_young
    eta1_old <- start.lst$eta1_old
    eta2_old <- start.lst$eta2_old
    xi1_young <- log(eta1_young)
    xi2_young <- log(eta2_young)
    xi1_old <- log(eta1_old)
    xi2_old <- log(eta2_old)
    sig <- start.lst$sig
    
    #-------------------------------
    # M-H algorithm
    #------------------------------- 
    for (i in 1:n.mcmc) {
        if (i %% 5 == 0) cat("iter:", i, "\r")
        flush.console()
        
        # ====================================================================
        #------------------------------------------------------------------
        # Sampling t given (r, theta^(i))
        #------------------------------------------------------------------
        # ====================================================================
        Sigma_lst_young <- lapply(1:n_subj_young, function(s) {
            Kdf <- computeCovDer1(idx1 = t_mat_young[, s], 
                                  idx2 = x_vec, tau = tau, h = h) 
            
            Kdd <- computeCovDer2(idx1 = t_mat_young[, s], tau = tau, h = h)
            Sigma <- Kff - quad.form.inv(Kdd, Kdf)
            Sigma <- (Sigma + t(Sigma)) / 2
            return(Sigma)
        })
        
        Sigma_lst_old <- lapply(1:n_subj_old, function(s) {
            Kdf <- computeCovDer1(idx1 = t_mat_old[, s],
                                  idx2 = x_vec, tau = tau, h = h) 
            
            Kdd <- computeCovDer2(idx1 = t_mat_old[, s], tau = tau, h = h)
            Sigma <- Kff - quad.form.inv(Kdd, Kdf)
            Sigma <- (Sigma + t(Sigma)) / 2
            return(Sigma)
        })

        # ---------------
        ## update t_young: M-H
        # ---------------
        # cat("\n", "update t_young")
        # if (i > burn) {
        #     tt_young <- draws[floor((i - burn) %/% thin), 1:2*n_subj_young]
        # } else {
        #     tt_young <- c(t_mat_young[1, ], t_mat_young[, 2])
        # }
        
        for (s in 1:n_subj_young) {
            # print(s)
            ts_star_young <- stats::runif(2,
                                          min = c(a_young_1, a_young_2),
                                          max = c(b_young_1, b_young_2))
            
            # ts_star_young <- stats::runif(
            #     2, 
            #     min = c(max(tt_young[s] - 0.2, a_young_1), 
            #             max(tt_young[s + n_subj_young] - 0.2, a_young_2)), 
            #     max = c(min(tt_young[s] + 0.2, b_young_1), 
            #             min(tt_young[s + n_subj_young] + 0.2, b_young_2)))
            
            # ts_star_young <- stats::runif(
            #     2,
            #     min = c(max(t_mat_young[1, s] - 0.2, a_young_1),
            #             max(t_mat_young[2, s] - 0.2, a_young_2)),
            #     max = c(min(t_mat_young[1, s] + 0.2, b_young_1),
            #             min(t_mat_young[2, s] + 0.2, b_young_2)))
            
            

            Kdf_new <- computeCovDer1(idx1 = ts_star_young, idx2 = x_vec, 
                                      tau = tau, h = h) 
            Kdd_new <- computeCovDer2(idx1 = ts_star_young, tau = tau, h = h)
            
            Sigma_new <- Kff - quad.form.inv(Kdd_new, Kdf_new)

            log_rho_young <- metropolis_log_ratio_ts_age_beta_multi_sig(
                ys = multi_y_young[, s],
                t_vec = t_mat_young[, s], 
                ts_new = ts_star_young,
                Sigma_s = Sigma_lst_young[[s]], 
                Sigma_s_new = Sigma_new,
                sig = sig,
                r1 = r1_young, 
                r2 = r2_young, 
                eta1 = eta1_young,
                eta2 = eta2_young,
                a1 = a_young_1, 
                a2 = a_young_2, 
                b1 = b_young_1, 
                b2 = b_young_2
                )

            if (log(runif(1)) < log_rho_young) {
                t_mat_young[, s] <- ts_star_young
                # print(dim(t_mat_young))
            }
        }
        
        # print(paste("t_mat_young", t_mat_young))
        # ---------------
        ## update t_old: M-H
        # ---------------
        # cat("\n", "update t_old")
        
        # if (i > burn) {
        #     tt_old <- draws[floor((i - burn) %/% thin), 1:2*n_subj_old]
        # } else {
        #     tt_old <- c(t_mat_old[1, ], t_mat_old[, 2])
        # }
        
        # tt_old <- draws[floor((i - burn) %/% thin), 
        #                 (2*n_subj_young + 1):(2*n_subj_young + 2*n_subj_old)]
        for (s in 1:n_subj_old) {
            ts_star_old <- stats::runif(2, min = c(a_old_1, a_old_2),
                                        max = c(b_old_1, b_old_2))
            
            # ts_star_old <- stats::runif(
            #     2, 
            #     min = c(max(tt_old[s] - 0.2, a_old_1), 
            #             max(tt_old[s + n_subj_old] - 0.2, a_old_2)), 
            #     max = c(min(tt_old[s] + 0.2, b_old_1), 
            #             min(tt_old[s + n_subj_old] + 0.2, b_old_2)))
            
            
            # ts_star_old <- stats::runif(
            #     2,
            #     min = c(max(t_mat_old[1, s] - 0.2, a_old_1),
            #             max(t_mat_old[2, s] - 0.2, a_old_2)),
            #     max = c(min(t_mat_old[1, s] + 0.2, b_old_1),
            #             min(t_mat_old[2, s] + 0.2, b_old_2)))
            
            
            Kdf_new <- computeCovDer1(idx1 = ts_star_old, idx2 = x_vec, 
                                      tau = tau, h = h) 
            Kdd_new <- computeCovDer2(idx1 = ts_star_old, tau = tau, h = h)
            Sigma_new <- Kff - emulator::quad.form.inv(Kdd_new, Kdf_new)
            
            log_rho_old <- metropolis_log_ratio_ts_age_beta_multi_sig(
                ys = multi_y_old[, s],
                t_vec = t_mat_old[, s], 
                ts_new = ts_star_old,
                Sigma_s = Sigma_lst_old[[s]], 
                Sigma_s_new = Sigma_new,
                sig = sig,
                r1 = r1_old, 
                r2 = r2_old, 
                eta1 = eta1_old,
                eta2 = eta2_old,
                a1 = a_old_1, 
                a2 = a_old_2, 
                b1 = b_old_1, 
                b2 = b_old_2
                )
            
            if (log(runif(1)) < log_rho_old) {
                t_mat_old[, s] <- ts_star_old
                # print(dim(t_mat_old))
            }
        }
        # print(paste("t_mat_old", t_mat_old))
        # ====================================================================
        # Sampling r given (t, theta^(i))
        # ====================================================================
        # cat("\n", "Sampling r")
        t_vec_1_young <- t_mat_young[1, ]
        t_vec_2_young <- t_mat_young[2, ]
        t_vec_1_old <- t_mat_old[1, ]
        t_vec_2_old <- t_mat_old[2, ]
        
        # ---------
        # r1_young: M-H
        # ---------
        r1_young_star <- stats::runif(1, min = 0, max = 1)
        
        log_rho_1_young <- metropolis_log_ratio_r_age_beta(
            t_vec = t_vec_1_young, 
            r_new = r1_young_star, 
            r = r1_young, 
            eta = eta1_young, 
            a_t = a_young_1, 
            b_t = b_young_1, 
            a_r = a_1,
            b_r = b_1,
            alpha = alpha_r_young_1, 
            beta = beta_r_young_1
            )
        
        if (log(runif(1)) < log_rho_1_young) {
            r1_young <- r1_young_star
        }

        # print(paste("r1_young", r1_young))
        # ---------
        # r2_young: M-H
        # ---------
        r2_young_star <- stats::runif(1, min = 0, max = 1)
        
        log_rho_2_young <- metropolis_log_ratio_r_age_beta(
            t_vec = t_vec_2_young, 
            r_new = r2_young_star, 
            r = r2_young, 
            eta = eta2_young, 
            a_t = a_young_2, 
            b_t = b_young_2, 
            a_r = a_2,
            b_r = b_2,
            alpha = alpha_r_young_2, 
            beta = beta_r_young_2)
        
        if (log(runif(1)) < log_rho_2_young) {
            r2_young <- r2_young_star
        }
        
        
        # print(paste("r2_young", r2_young))
        # ---------
        # r1_old: M-H
        # ---------
        r1_old_star <- stats::runif(1, min = 0, max = 1)
        
        log_rho_1_old <- metropolis_log_ratio_r_age_beta(
            t_vec = t_vec_1_old, 
            r_new = r1_old_star, 
            r = r1_old, 
            eta = eta1_old, 
            a_t = a_old_1, 
            b_t = b_old_1, 
            a_r = a_1,
            b_r = b_1,
            alpha = alpha_r_old_1, 
            beta = beta_r_old_1)
        
        if (log(runif(1)) < log_rho_1_old) {
            r1_old <- r1_old_star
        }
        # print(paste("r1_old", r1_old))
        # ---------
        # r2_old: M-H
        # ---------
        r2_old_star <-  stats::runif(1, min = 0, max = 1)
        
        log_rho_2_old <- metropolis_log_ratio_r_age_beta(
            t_vec = t_vec_2_old, 
            r_new = r2_old_star, 
            r = r2_old, 
            eta = eta2_old, 
            a_t = a_old_2, 
            b_t = b_old_2, 
            a_r = a_2,
            b_r = b_2,
            alpha = alpha_r_old_2, 
            beta = beta_r_old_2)
        
        if (log(runif(1)) < log_rho_2_old) {
            r2_old <- r2_old_star
        }
        # print(paste("r1_old", r1_old))
        # ====================================================================
        # Sampling eta given (t, r, theta^(i))
        # ====================================================================
        # Update tuning parameter
        #------------------------
        if (adapt == TRUE & i %% Tb == 0) { 
            # Adaptive tuning
            keep.tmp.young.eta.1 <- keep.tmp.young.eta.1 / Tb
            keep.tmp.young.eta.2 <- keep.tmp.young.eta.2 / Tb
            keep.tmp.old.eta.1 <- keep.tmp.old.eta.1 / Tb
            keep.tmp.old.eta.2 <- keep.tmp.old.eta.2 / Tb
            
            tune.young.eta.1 <- get.tune(tune.young.eta.1, keep.tmp.young.eta.1, 
                                         k = i)
            tune.young.eta.2 <- get.tune(tune.young.eta.2, keep.tmp.young.eta.2, 
                                         k = i)
            tune.old.eta.1 <- get.tune(tune.old.eta.1, keep.tmp.old.eta.1, 
                                       k = i)
            tune.old.eta.2 <- get.tune(tune.old.eta.2, keep.tmp.old.eta.2, 
                                       k = i)
            
            keep.tmp.young.eta.1 <- 0
            keep.tmp.young.eta.2 <- 0
            keep.tmp.old.eta.1 <- 0
            keep.tmp.old.eta.2 <- 0
        } 	

        
        # eta1_young
        xi1_young_star <- rnorm(1, xi1_young, tune.young.eta.1)
        log_rho_1_young_eta <- metropolis_log_ratio_eta_age_beta_ga(
            t_vec = t_vec_1_young, 
            r = r1_young, 
            xi_new = xi1_young_star, 
            xi = xi1_young, 
            a = a_young_1, 
            b = b_young_1,
            a_eta = a_eta_young_1, 
            b_eta = b_eta_young_1)
        
        if (log(runif(1)) < log_rho_1_young_eta) {
            xi1_young <- xi1_young_star
            eta1_young <- exp(xi1_young)
            keep.young.eta.1 <- keep.young.eta.1 + 1
            keep.tmp.young.eta.1 <- keep.tmp.young.eta.1 + 1
        }
        
        # print(paste("eta1_young", eta1_young))
        
        
        # eta2_young
        xi2_young_star <- rnorm(1, xi2_young, tune.young.eta.2)
        log_rho_2_young_eta <- metropolis_log_ratio_eta_age_beta_ga(
            t_vec = t_vec_2_young, 
            r = r2_young, 
            xi_new = xi2_young_star, 
            xi = xi2_young, 
            a = a_young_2, 
            b = b_young_2,
            a_eta = a_eta_young_2, 
            b_eta = b_eta_young_2)
        
        if (log(runif(1)) < log_rho_2_young_eta) {
            xi2_young <- xi2_young_star
            eta2_young <- exp(xi2_young)
            keep.young.eta.2 <- keep.young.eta.2 + 1
            keep.tmp.young.eta.2 <- keep.tmp.young.eta.2 + 1
        }
        
        # print(paste("eta2_young", eta2_young))
        
        
        # eta1_old
        xi1_old_star <- rnorm(1, xi1_old, tune.old.eta.1)
        log_rho_1_old_eta <- metropolis_log_ratio_eta_age_beta_ga(
            t_vec = t_vec_1_old, 
            r = r1_old, 
            xi_new = xi1_old_star, 
            xi = xi1_old, 
            a = a_old_1, 
            b = b_old_1,
            a_eta = a_eta_old_1,
            b_eta = b_eta_old_1)
        
        if (log(runif(1)) < log_rho_1_old_eta) {
            xi1_old <- xi1_old_star
            eta1_old <- exp(xi1_old)
            keep.old.eta.1 <- keep.old.eta.1 + 1
            keep.tmp.old.eta.1 <- keep.tmp.old.eta.1 + 1
        }
        
        # print(paste("eta1_old", eta1_old))
        
        # eta2_old
        xi2_old_star <- rnorm(1, xi2_old, tune.old.eta.2)
        log_rho_2_old_eta <- metropolis_log_ratio_eta_age_beta_ga(
            t_vec = t_vec_2_old, 
            r = r2_old, 
            xi_new = xi2_old_star, 
            xi = xi2_old, 
            a = a_old_2, 
            b = b_old_2,
            a_eta = a_eta_old_2, 
            b_eta = b_eta_old_2)
        
        if (log(runif(1)) < log_rho_2_old_eta) {
            xi2_old <- xi2_old_star
            eta2_old <- exp(xi2_old)
            keep.old.eta.2 <- keep.old.eta.2 + 1
            keep.tmp.old.eta.2 <- keep.tmp.old.eta.2 + 1
        }
        
        # print(paste("eta2_old", eta2_old))
        
        
        # ====================================================================
        # Sampling sigma^2 given (t, r, theta^(i))
        # ====================================================================
        
        sig <- sqrt(sample_sig2_subj_age(multi_y_young = multi_y_young, 
                                         multi_y_old = multi_y_old, 
                                         x_vec = x_vec, theta = theta_k, 
                                         Sigma_young = Sigma_lst_young, 
                                         Sigma_old = Sigma_lst_old, 
                                         ga_shape = ga_shape, 
                                         ga_rate = ga_rate))
        
        #  Save samples 
        # -----------------------------------------------------

        if (i > burn) {
            if (i %% thin == 0) {
                # t_vec_1_young <- t_mat_young[1, ]
                # t_vec_2_young <- t_mat_young[2, ]
                # t_vec_1_old <- t_mat_old[1, ]
                # t_vec_2_old <- t_mat_old[2, ]
                # print(length(c(t_vec_1_young,
                #                t_vec_2_young,
                #                t_vec_1_old,
                #                t_vec_2_old)))
                # 
                # print(length(c(r1_young, 
                #                r2_young,
                #                r1_old, 
                #                r2_old,
                #                eta1_young,
                #                eta2_young,
                #                eta1_old,
                #                eta2_old,
                #                sig)))
                # print(dim(draws))
                draws[(i - burn) %/% thin, ] <- c(t_mat_young[1, ], 
                                                  t_mat_young[2, ],
                                                  t_mat_old[1, ],
                                                  t_mat_old[2, ],
                                                  r1_young, 
                                                  r2_young,
                                                  r1_old, 
                                                  r2_old,
                                                  eta1_young,
                                                  eta2_young,
                                                  eta1_old,
                                                  eta2_old,
                                                  sig)
            } 
        }
    }
    
    #-----------------------------
    # Acceptance Probability
    #----------------------------
    keep.young <- lapply(keep.young, function(x) x / n.mcmc)
    keep.old <- lapply(keep.old, function(x) x / n.mcmc)
    
    keep.young.eta.1 <- keep.young.eta.1 / n.mcmc
    keep.young.eta.2 <- keep.young.eta.2 / n.mcmc
    keep.old.eta.1 <- keep.old.eta.1 / n.mcmc
    keep.old.eta.2 <- keep.old.eta.2 / n.mcmc
    
    #--------------
    # Write output
    #--------------
    # print(draws)
    return(list(draws = draws, 
                # accept_young = keep.young,
                # accept_old = keep.old,
                accept_young_eta_1 = keep.young.eta.1,
                accept_young_eta_2 = keep.young.eta.2,
                accept_old_eta_1 = keep.old.eta.1,
                accept_old_eta_2 = keep.old.eta.2,
                start = start.lst, 
                # tune_young = tune.lst.young, 
                # tune_old = tune.lst.old, 
                tune_young_eta_1 = tune.young.eta.1, 
                tune_young_eta_2 = tune.young.eta.2, 
                tune_old_eta_1 = tune.old.eta.1, 
                tune_old_eta_2 = tune.old.eta.2, 
                burn = burn, 
                thin = thin, 
                n.mcmc = n.mcmc, 
                sampleidx = sampleidx))
}