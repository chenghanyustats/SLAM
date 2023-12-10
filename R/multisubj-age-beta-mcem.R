## careful about M-step
# The algorithm fixes the number of groups and assumes known number of t's
# a_eta_young_1 = 1 / 2
# b_eta_young_1 = 1 / 2
# a_eta_young_2 = 1 / 2
# b_eta_young_2 = 1 / 2
# a_eta_old_1 = 1 / 2
# b_eta_old_1 = 1 / 2
# a_eta_old_2 = 1 / 2
# b_eta_old_2 = 1 / 2


multisubj_age_beta_r_ga_mcem <- function(multi_y_young, multi_y_old, x, H0, 
                                         theta_init = c(1, 1), 
                                         epsilon = 1e-4, D = 100, max_iter = 100,
                                         lower = c(0.001, 0.001),
                                         upper = c(1/0.001, 1/0.001), 
                                         ctrl = list(TOL = 1e-5, trace = 0),
                                         ga_shape = 6, ga_rate = 5,
                                         is.h.par = FALSE,
                                         n_mcmc_e = 1000, 
                                         burn_e = 10, 
                                         thin_e = 1, 
                                         n_mcmc_final = 1000, 
                                         burn_final = 10, 
                                         thin_final = 1, 
                                         start.lst, name.par, adapt = TRUE, 
                                         keep.young.eta.1, keep.young.eta.2,
                                         keep.old.eta.1, keep.old.eta.2,
                                         tune.young.eta.1, tune.young.eta.2,
                                         tune.old.eta.1, tune.old.eta.2,
                                         tune.len = 40, 
                                         target.accept.rate = 0.35,
                                         a_h = 1, b_h = 1,
                                         a_young_1, b_young_1, 
                                         a_young_2, b_young_2,
                                         a_young_3, b_young_3,
                                         a_old_1, b_old_1, 
                                         a_old_2, b_old_2,
                                         alpha_r_young_1 = 1, beta_r_young_1 = 1,
                                         alpha_r_young_2 = 1, beta_r_young_2 = 1,
                                         alpha_r_old_1 = 1, beta_r_old_1 = 1,
                                         alpha_r_old_2 = 1, beta_r_old_2 = 1,
                                         a_eta_young_1 = 1 / 2, 
                                         b_eta_young_1 = 1 / 2,
                                         a_eta_young_2 = 1 / 2, 
                                         b_eta_young_2 = 1 / 2,
                                         a_eta_old_1 = 1 / 2, 
                                         b_eta_old_1 = 1 / 2,
                                         a_eta_old_2 = 1 / 2, 
                                         b_eta_old_2 = 1 / 2,
                                         final_draw = TRUE) {
    
    ################################
    # Arguments:
    # Rewrite it #######
    #   multi_y: matrix of multi-subject observations
    #   x: univariate (time) inputs same for all subjects
    #   H0: difference matrix of x
    #   theta_init: initial value of theta = (sig, tau, h)
    #   espilon: tolerence level for convergence
    #   D: MCMC sample size in the E-step
    #   max_iter: maximum number of iterations of the stochastic EM algorithm
    #   lower: lower bound for theta in the optimization
    #   upper: upper bound for theta in the optimization
    #   ctrl: control parameter used in the optimization
    #   is.sig.par: put a IG prior on sigma^2 or not
    #   n.mcmc: number of MCMC iterations for the final samples of t and r
    #   start.lst: initial values of t and r
    #   burn: burn-in period
    #   thin: thining 
    #   name.par: names of parameters to be sampled in MCMC
    #   adapt: whether to adaptive sampling or not
    #   tune.lst: tuning variance in the random walk noraml proposal
    #   keep.lst: number of times the draws is accepted
    #   tune.len: how often to update the tunning variance
    #   target.accept.rate: targeted accepted rate of the metropolis 
    #   proposal_type: "normal", "uniform" or "transform_normal"
    #   intvl: closed interval of inputs
    ################################
    
    ## =========================================================================
    ### create objects for later use
    ## =========================================================================
    n_theta <- length(theta_init)
    n_subj_young <- ncol(multi_y_young)
    n_subj_old <- ncol(multi_y_old)
    n_subj <- n_subj_old + n_subj_young
    n_t <- n_subj * 2
    # n_r <- 2 * 2 ## two groups (young and old) and two erp components
    # a <- intvl[1]
    # b <- intvl[2]
    # a_1 <- max(a_young_1, a_old)
    a_1 <- max(a_young_1, a_old_1)
    b_1 <- min(b_young_1, b_old_1)
    a_2 <- max(a_young_2, a_old_2)
    b_2 <- min(b_young_2, b_old_2)
    
    eps <- epsilon + 1
    count <- 1
    
    # Initialization of parameters
    theta_mat <- matrix(0L, max_iter, n_theta)
    theta_mat[1, ] <- theta_k <- theta_init
    theta_new <- theta_k
    mar_post_new <- 0
    
    # alpha_eta_star <- (n_subj - 1) + 2 * alpha_eta
    
    
    ## =========================================================================
    ## Multi-subject (young and old) MCEM Algorithm for DGP
    ## =========================================================================
    print("multisubj_age_beta_r_ga_mcem (f_js, theta, t_js, r_j1, r_j2) begins")
    while(eps > epsilon) {
        
        if (count == 100) {
            print("Does not converge")
            break
        }
        
        # *****************************************
        # ******** MC E-step ************ #
        #  ****************************************
        ## MCMC algorithm for the MC E-step
        if (count == 1) {
            start.lst <- start.lst
        } else{
            ## rewrite for r1 and r2
            post_mean <- apply(mcmc_output$draws, 2, mean)
            start.lst <- list(t_young = matrix(runif(n_subj_young * 2, 
                                                     min = c(a_young_1, a_young_2), 
                                                     max = c(b_young_1, b_young_2)), 
                                               2, n_subj_young),
                              t_old = matrix(runif(n_subj_old * 2, 
                                                   min = c(a_old_1, a_old_2), 
                                                   max = c(b_old_1, b_old_2)), 
                                             2, n_subj_old),
                              r1_young = post_mean[n_t + 1],
                              r2_young = post_mean[n_t + 2],
                              r1_old = post_mean[n_t + 3],
                              r2_old = post_mean[n_t + 4], 
                              eta1_young = post_mean[n_t + 5],
                              eta2_young = post_mean[n_t + 6],
                              eta1_old = post_mean[n_t + 7],
                              eta2_old = post_mean[n_t + 8],
                              sigma = post_mean[n_t + 9])
            
            # start_lst <- list(t_young = matrix(runif(no_data * 2,
            #                                          min = c(a_sin_1, a_sin_2),
            #                                          max = c(b_sin_1, b_sin_2)), 2, no_data),
            #                   t_old = matrix(runif(no_data * 2,
            #                                        min = c(a_cos_1, a_cos_2),
            #                                        max = c(b_cos_1, b_cos_2)), 2, no_data),
            #                   r1_young = 0.3,
            #                   r2_young = 0.7,
            #                   r1_old = 0.3,
            #                   r2_old = 0.7,
            #                   eta1_young = 1,
            #                   eta2_young = 1,
            #                   eta1_old = 1,
            #                   eta2_old = 1,
            #                   sigma = 1)
        }
        
        print("E-step sampling for t, r, eta and sigma")
        # print(paste("startlist", dim(start.lst$t_young)))
        mcmc_output <- mcmc_Estep_age_beta_r_ga_sig(
            multi_y_young = multi_y_young, 
            multi_y_old = multi_y_old, 
            x_vec = x, 
            n_subj_young = n_subj_young, 
            n_subj_old = n_subj_old,
            theta_k = theta_k, H0 = H0, 
            start.lst = start.lst,
            n.mcmc = n_mcmc_e, 
            burn = burn_e, 
            thin = thin_e, 
            name.par = name.par, adapt = adapt, 
            keep.young.eta.1 = keep.young.eta.1, 
            keep.young.eta.2 = keep.young.eta.2,
            keep.old.eta.1 = keep.old.eta.1, 
            keep.old.eta.2 = keep.old.eta.2,
            tune.young.eta.1 = tune.young.eta.1, 
            tune.young.eta.2 = tune.young.eta.2,
            tune.old.eta.1 = tune.old.eta.1, 
            tune.old.eta.2 = tune.old.eta.2,
            tune.len = tune.len, 
            target.accept.rate = target.accept.rate,
            a_young_1 = a_young_1, 
            b_young_1 = b_young_1, 
            a_young_2 = a_young_2, 
            b_young_2 = b_young_2,
            a_old_1 = a_old_1, 
            b_old_1 = b_old_1, 
            a_old_2 = a_old_2, 
            b_old_2 = b_old_2,
            a_1 = a_1, 
            b_1 = b_1, 
            a_2 = a_2,
            b_2 = b_2,
            alpha_r_young_1 = alpha_r_young_1, 
            beta_r_young_1 = beta_r_young_1,
            alpha_r_young_2 = alpha_r_young_2, 
            beta_r_young_2 = beta_r_young_2,
            alpha_r_old_1 = alpha_r_old_1, 
            beta_r_old_1 = beta_r_old_1,
            alpha_r_old_2 = alpha_r_old_2,
            beta_r_old_2 = beta_r_old_2,
            a_eta_young_1 = a_eta_young_1, 
            b_eta_young_1 = b_eta_young_1,
            a_eta_young_2 = a_eta_young_2, 
            b_eta_young_2 = b_eta_young_2,
            a_eta_old_1 = a_eta_old_1, 
            b_eta_old_1 = b_eta_old_1,
            a_eta_old_2 = a_eta_old_2, 
            b_eta_old_2 = b_eta_old_2,
            ga_shape = ga_shape, ga_rate = ga_rate)
        
        sample_t_young <- mcmc_output$draws[, 1:(2*n_subj_young)]
        sample_t_old <- mcmc_output$draws[, (2*n_subj_young + 1):(2*n_subj_young + 2*n_subj_old)]
        sample_sig <- mcmc_output$draws[, ncol(mcmc_output$draws)]
        # format_for_M_step <- function(sample_t) {
        #     sample_t_lst <- list()
        #     D <- nrow(sample_t)
        #     S <- ncol(sample_t) / 2
        #     for (s in 1:S) {
        #         sample_t_lst[[s]] <- t(sample_t[, c(s, s + S)])
        #     }
        #     return(sample_t_lst)
        # }
        
        sample_t_young_mstep <- apply(sample_t_young, 2, sample, size = D)
        sample_t_old_mstep <- apply(sample_t_old, 2, sample, size = D)
        sample_sig_mstep <- sample(sample_sig, size = D)
        
        
        sample_t_young_lst <- format_for_M_step(sample_t_young_mstep)
        sample_t_old_lst <- format_for_M_step(sample_t_old_mstep)    
        ## the sample should be dimension nn(=2) X D
        
        
        
        # sample_t_mstep <- apply(sample_t, 2, sample, size = D)
        
        
        
        
        
        # *****************************************
        # ******** M-step for theta ************* #
        # *****************************************
        print("M-step updating kernel parameters")
        
        res <- Rsolnp::solnp(pars = theta_k,
                             fun = log_mar_lik_gp_der_mc_age_multi_sig,
                             LB = lower, control = ctrl,
                             multi_y_young = multi_y_young,
                             multi_y_old = multi_y_old,
                             x_vec = x, H0 = H0,
                             der_mc_mat_young_lst = sample_t_young_lst,
                             der_mc_mat_old_lst = sample_t_old_lst,
                             sample_sig = sample_sig_mstep,
                             is.h.par = is.h.par,
                             a_h = a_h, b_h = b_h)
        theta_k <- res$par
        mar_post_k <- res$values[length(res$values)]
        
        # res <- stats::optim(par = theta_k,
        #                     fn = log_mar_lik_gp_der_mc_age_multi,
        #                     lower = lower, control = ctrl,
        #                     multi_y_young = multi_y_young,
        #                     multi_y_old = multi_y_old,
        #                     x_vec = x, H0 = H0,
        #                     der_mc_mat_young_lst = sample_t_young_lst,
        #                     der_mc_mat_old_lst = sample_t_old_lst,
        #                     ga_shape = ga_shape, ga_rate = ga_rate,
        #                     is.sig.par = is.sig.par,
        #                     is.h.par = is.h.par,
        #                     a_h = a_h, b_h = b_h)
        # theta_k <- res$par
        # mar_post_k <- res$value
        
        print(paste("theta =", theta_k))
        print(paste("(increasing) marginal posterior =", -mar_post_k))
        
        # ******** Update epsilon and parameters ************ #
        eps <- sum((theta_new - theta_k) ^ 2)
        theta_new <- theta_k
        mar_post_new <- mar_post_k
        count <- count + 1
        print(paste("eps =", round(eps, 4), " count", count))
        theta_mat[count, ] <- theta_new
        
    }
    
    
    ## =========================================================================
    ## Final MCMC samples for t, r, eta, sigma
    ## =========================================================================
    if (final_draw) {
        print("Final sampling for t, r, eta, sigma")
        mcmc_output <- mcmc_Estep_age_beta_r_ga_sig(
            multi_y_young = multi_y_young, 
            multi_y_old = multi_y_old, 
            x_vec = x, 
            n_subj_young = n_subj_young, 
            n_subj_old = n_subj_old,
            theta_k = theta_k, H0 = H0, 
            start.lst = start.lst,
            n.mcmc = n_mcmc_final, 
            burn = burn_final, 
            thin = thin_final, 
            name.par = name.par, adapt = adapt, 
            
            keep.young.eta.1 = keep.young.eta.1, 
            keep.young.eta.2 = keep.young.eta.2,
            keep.old.eta.1 = keep.old.eta.1, 
            keep.old.eta.2 = keep.old.eta.2,
            tune.young.eta.1 = tune.young.eta.1, 
            tune.young.eta.2 = tune.young.eta.2,
            tune.old.eta.1 = tune.old.eta.1, 
            tune.old.eta.2 = tune.old.eta.2,
            tune.len = tune.len, 
            target.accept.rate = target.accept.rate,
            
            a_young_1 = a_young_1, 
            b_young_1 = b_young_1, 
            a_young_2 = a_young_2, 
            b_young_2 = b_young_2,
            a_old_1 = a_old_1, 
            b_old_1 = b_old_1, 
            a_old_2 = a_old_2, 
            b_old_2 = b_old_2,
            a_1 = a_1, 
            b_1 = b_1, 
            a_2 = a_2,
            b_2 = b_2,
            alpha_r_young_1 = alpha_r_young_1, 
            beta_r_young_1 = beta_r_young_1,
            alpha_r_young_2 = alpha_r_young_2, 
            beta_r_young_2 = beta_r_young_2,
            alpha_r_old_1 = alpha_r_old_1, 
            beta_r_old_1 = beta_r_old_1,
            alpha_r_old_2 = alpha_r_old_2,
            beta_r_old_2 = beta_r_old_2,
            a_eta_young_1 = a_eta_young_1, 
            b_eta_young_1 = b_eta_young_1,
            a_eta_young_2 = a_eta_young_2, 
            b_eta_young_2 = b_eta_young_2,
            a_eta_old_1 = a_eta_old_1, 
            b_eta_old_1 = b_eta_old_1,
            a_eta_old_2 = a_eta_old_2, 
            b_eta_old_2 = b_eta_old_2,
            ga_shape = ga_shape, ga_rate = ga_rate)
    }

    cat("\n")
    print("Done!")
    colnames(theta_mat) <- c("tau", "h")
    
    return(list(mcmc_output = mcmc_output, 
                theta_mat = theta_mat[1:count, ]))
    
}
