
# young: g1
# old: g2


Estep_one_way <- function(multi_y_g1, multi_y_g2, x_vec,
                          n_subj_g1, n_subj_g2, theta_k, H0,
                          start.lst,
                          n.mcmc = 1000, burn = 100, thin = 1,
                          name.par, adapt = TRUE,
                          tune.lst.g1, keep.lst.g1,
                          tune.lst.g2, keep.lst.g2,
                          keep.g1.eta.1, keep.g1.eta.2,
                          keep.g2.eta.1, keep.g2.eta.2,
                          tune.g1.eta.1, tune.g1.eta.2,
                          tune.g2.eta.1, tune.g2.eta.2,
                          keep.beta0.1,
                          keep.beta0.2,
                          keep.beta1.1,
                          keep.beta1.2,
                          tune.beta0.1,
                          tune.beta0.2,
                          tune.beta1.1,
                          tune.beta1.2,
                          tune.len = 50, target.accept.rate = 0.35,
                          a_1, b_1, a_2, b_2,
                          a_eta_g1_1,
                          b_eta_g1_1 ,
                          a_eta_g1_2,
                          b_eta_g1_2,
                          a_eta_g2_1,
                          b_eta_g2_1,
                          a_eta_g2_2,
                          b_eta_g2_2,
                          ga_shape, ga_rate,
                          mu0_1, sig0_1,
                          mu0_2, sig0_2,
                          mu1_1, sig1_1,
                          mu1_2, sig1_2,
                          link) {

  #-------------------------------
  # Adaptive tuning
  #-------------------------------
  Tb <- tune.len  # frequency of adaptive tuning

  keep.tmp.lst.g1 <- keep.lst.g1
  keep.tmp.lst.g2 <- keep.lst.g2

  # for eta
  keep.tmp.g1.eta.1 <- keep.g1.eta.1  # track MH acceptance rate for adaptive tuning
  keep.tmp.g1.eta.2 <- keep.g1.eta.2
  keep.tmp.g2.eta.1 <- keep.g2.eta.1
  keep.tmp.g2.eta.2 <- keep.g2.eta.2

  keep.tmp.beta0.1 <- keep.beta0.1
  keep.tmp.beta0.2 <- keep.beta0.2
  keep.tmp.beta1.1 <- keep.beta1.1
  keep.tmp.beta1.2 <- keep.beta1.2

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
  t_mat_g1 <- start.lst$t_g1
  t_mat_g2 <- start.lst$t_g2

  t_vec_1_g1 <- t_mat_g1[1, ]
  t_vec_2_g1 <- t_mat_g1[2, ]
  t_vec_1_g2 <- t_mat_g2[1, ]
  t_vec_2_g2 <- t_mat_g2[2, ]

  beta0_1 <- start.lst$beta0_1
  beta0_2 <- start.lst$beta0_2
  beta1_1 <- start.lst$beta1_1
  beta1_2 <- start.lst$beta1_2

  eta1_g1 <- start.lst$eta1_g1
  eta2_g1 <- start.lst$eta2_g1
  eta1_g2 <- start.lst$eta1_g2
  eta2_g2 <- start.lst$eta2_g2
  xi1_g1 <- log(eta1_g1)
  xi2_g1 <- log(eta2_g1)
  xi1_g2 <- log(eta1_g2)
  xi2_g2 <- log(eta2_g2)

  if (link == "logit") {
    r1_g1 <- 1 / (1 + exp(-beta0_1))
    r2_g1 <- 1 / (1 + exp(-beta0_2))
    r1_g2 <- 1 / (1 + exp(-(beta0_1 + beta1_1)))
    r2_g2 <- 1 / (1 + exp(-(beta0_2 + beta1_2)))
  }

  if (link == "probit") {
    r1_g1 <- pnorm(beta0_1)
    r2_g1 <- pnorm(beta0_2)
    r1_g2 <- pnorm(beta0_1 + beta1_1)
    r2_g2 <- pnorm(beta0_2 + beta1_2)
  }

  if (link == "cloglog") {
    r1_g1 <- 1 - exp(-exp(beta0_1))
    r2_g1 <- 1 - exp(-exp(beta0_2))
    r1_g2 <- 1 - exp(-exp(beta0_1 + beta1_1))
    r2_g2 <- 1 - exp(-exp(beta0_2 + beta1_2))
  }

  sig <- start.lst$sig

  #-------------------------------
  # M-H algorithm
  #-------------------------------
  for (i in 1:n.mcmc) {
    # if (i %% 5 == 0) cat("iter:", i, "\r")
    # flush.console()

    # ====================================================================
    #------------------------------------------------------------------
    # Sampling t given (r, theta^(i))
    #------------------------------------------------------------------
    # ====================================================================
    Sigma_lst_g1 <- lapply(1:n_subj_g1, function(s) {
      Kdf <- computeCovDer1(idx1 = t_mat_g1[, s],
                            idx2 = x_vec, tau = tau, h = h)

      Kdd <- computeCovDer2(idx1 = t_mat_g1[, s], tau = tau, h = h)
      Sigma <- Kff - quad.form.inv(Kdd, Kdf)
      return((Sigma + t(Sigma)) / 2)
    })

    Sigma_lst_g2 <- lapply(1:n_subj_g2, function(s) {
      Kdf <- computeCovDer1(idx1 = t_mat_g2[, s],
                            idx2 = x_vec, tau = tau, h = h)

      Kdd <- computeCovDer2(idx1 = t_mat_g2[, s], tau = tau, h = h)
      Sigma <- Kff - quad.form.inv(Kdd, Kdf)
      return((Sigma + t(Sigma)) / 2)
    })

    # Update tuning parameter
    #------------------------
    if (adapt == TRUE & i %% Tb == 0) {
      # Adaptive tuning
      for (s in 1:n_subj_g1) {
        keep.tmp.lst.g1[[s]] <- keep.tmp.lst.g1[[s]] / Tb
        tune.lst.g1[[s]] <- get.tune(tune.lst.g1[[s]],
                                        keep.tmp.lst.g1[[s]], k = i,
                                     target = target.accept.rate)
        keep.tmp.lst.g1[[s]] <- 0
      }
      for (s in 1:n_subj_g2) {
        keep.tmp.lst.g2[[s]] <- keep.tmp.lst.g2[[s]] / Tb
        tune.lst.g2[[s]] <- get.tune(tune.lst.g2[[s]],
                                      keep.tmp.lst.g2[[s]], k = i,
                                     target = target.accept.rate)
        keep.tmp.lst.g2[[s]] <- 0
      }

      # Adaptive tuning
      keep.tmp.beta0.1 <- keep.tmp.beta0.1 / Tb
      keep.tmp.beta0.2 <- keep.tmp.beta0.2 / Tb
      keep.tmp.beta1.1 <- keep.tmp.beta1.1 / Tb
      keep.tmp.beta1.2 <- keep.tmp.beta1.2 / Tb

      tune.beta0.1 <- get.tune(tune.beta0.1, keep.tmp.beta0.1,
                               k = i, target = target.accept.rate)
      tune.beta0.2 <- get.tune(tune.beta0.2, keep.tmp.beta0.2,
                               k = i, target = target.accept.rate)
      tune.beta1.1 <- get.tune(tune.beta1.1, keep.tmp.beta1.1,
                               k = i, target = target.accept.rate)
      tune.beta1.2 <- get.tune(tune.beta1.2, keep.tmp.beta1.2,
                               k = i, target = target.accept.rate)

      keep.tmp.beta0.1 <- 0
      keep.tmp.beta0.2 <- 0
      keep.tmp.beta1.1 <- 0
      keep.tmp.beta1.2 <- 0

      # Adaptive tuning
      keep.tmp.g1.eta.1 <- keep.tmp.g1.eta.1 / Tb
      keep.tmp.g1.eta.2 <- keep.tmp.g1.eta.2 / Tb
      keep.tmp.g2.eta.1 <- keep.tmp.g2.eta.1 / Tb
      keep.tmp.g2.eta.2 <- keep.tmp.g2.eta.2 / Tb

      tune.g1.eta.1 <- get.tune(tune.g1.eta.1, keep.tmp.g1.eta.1,
                                   k = i, target = target.accept.rate)
      tune.g1.eta.2 <- get.tune(tune.g1.eta.2, keep.tmp.g1.eta.2,
                                   k = i, target = target.accept.rate)
      tune.g2.eta.1 <- get.tune(tune.g2.eta.1, keep.tmp.g2.eta.1,
                                 k = i, target = target.accept.rate)
      tune.g2.eta.2 <- get.tune(tune.g2.eta.2, keep.tmp.g2.eta.2,
                                 k = i, target = target.accept.rate)

      keep.tmp.g1.eta.1 <- 0
      keep.tmp.g1.eta.2 <- 0
      keep.tmp.g2.eta.1 <- 0
      keep.tmp.g2.eta.2 <- 0
    }

    # ---------------
    ## update t_g1: M-H
    # ---------------
    for (s in 1:n_subj_g1) {

      ts_star_g1 <-
        truncnorm::rtruncnorm(n = 2, a = c(a_1, a_2), b = c(b_1, b_2),
                              mean = c(t_vec_1_g1[s], t_vec_2_g1[s]),
                              sd = tune.lst.g1[[s]])

      Kdf_new <- computeCovDer1(idx1 = ts_star_g1, idx2 = x_vec,
                                tau = tau, h = h)
      Kdd_new <- computeCovDer2(idx1 = ts_star_g1, tau = tau, h = h)

      Sigma_new <- Kff - quad.form.inv(Kdd_new, Kdf_new)

      log_rho_g1 <- metro_log_ratio_t(
        ys = multi_y_g1[, s],
        t_vec = t_mat_g1[, s],
        ts_new = ts_star_g1,
        Sigma_s = Sigma_lst_g1[[s]],
        Sigma_s_new = Sigma_new,
        sig = sig,
        r1 = r1_g1,
        r2 = r2_g1,
        eta1 = eta1_g1,
        eta2 = eta2_g1,
        a1 = a_1,
        a2 = a_2,
        b1 = b_1,
        b2 = b_2,
        sd_trun = tune.lst.g1[[s]]
      )

      if (log(runif(1)) < log_rho_g1) {
        t_mat_g1[, s] <- ts_star_g1
        keep.lst.g1[[s]] <- keep.lst.g1[[s]] + 1
        keep.tmp.lst.g1[[s]] <- keep.tmp.lst.g1[[s]] + 1
      }
    }

    for (s in 1:n_subj_g2) {
      ts_star_g2 <-
        truncnorm::rtruncnorm(n = 2, a = c(a_1, a_2), b = c(b_1, b_2),
                              mean = c(t_vec_1_g2[s], t_vec_2_g2[s]),
                              sd = tune.lst.g2[[s]])
      Kdf_new <- computeCovDer1(idx1 = ts_star_g2, idx2 = x_vec,
                                tau = tau, h = h)
      Kdd_new <- computeCovDer2(idx1 = ts_star_g2, tau = tau, h = h)
      Sigma_new <- Kff - emulator::quad.form.inv(Kdd_new, Kdf_new)

      log_rho_g2 <- metro_log_ratio_t(
        ys = multi_y_g2[, s],
        t_vec = t_mat_g2[, s],
        ts_new = ts_star_g2,
        Sigma_s = Sigma_lst_g2[[s]],
        Sigma_s_new = Sigma_new,
        sig = sig,
        r1 = r1_g2,
        r2 = r2_g2,
        eta1 = eta1_g2,
        eta2 = eta2_g2,
        a1 = a_1,
        a2 = a_2,
        b1 = b_1,
        b2 = b_2,
        sd_trun = tune.lst.g2[[s]]
      )

      if (log(runif(1)) < log_rho_g2) {
        t_mat_g2[, s] <- ts_star_g2
        keep.lst.g2[[s]] <- keep.lst.g2[[s]] + 1
        keep.tmp.lst.g2[[s]] <- keep.tmp.lst.g2[[s]] + 1
      }
    }

    # ====================================================================
    # Sampling r (regression coefficients beta) given (t, theta^(i))
    # ====================================================================

    t_vec_1_g1 <- t_mat_g1[1, ]
    t_vec_2_g1 <- t_mat_g1[2, ]
    t_vec_1_g2 <- t_mat_g2[1, ]
    t_vec_2_g2 <- t_mat_g2[2, ]



    # ---------
    # beta0_1: M-H (g1 is the baseline level (z = 0))
    # ---------

    beta0_1_star <- rnorm(1, beta0_1, tune.beta0.1)


    log_rho_beta0_1 <- metro_log_ratio_b0(
      t_vec_g1 = t_vec_1_g1,
      t_vec_g2 = t_vec_1_g2,
      beta0_new = beta0_1_star,
      beta0 = beta0_1,
      beta1 = beta1_1,
      eta_g1 = eta1_g1,
      eta_g2 = eta1_g2,
      a_t = a_1,
      b_t = b_1,
      mu0 = mu0_1,
      sig0 = sig0_1,
      link = link
    )

    if (log(runif(1)) < log_rho_beta0_1) {
      beta0_1 <- beta0_1_star
      keep.beta0.1 <- keep.beta0.1 + 1
      keep.tmp.beta0.1 <- keep.tmp.beta0.1 + 1
    }


    # ---------
    # beta0_2: M-H
    # ---------
    beta0_2_star <- rnorm(1, beta0_2, tune.beta0.2)

    log_rho_beta0_2 <- metro_log_ratio_b0(
      t_vec_g1 = t_vec_2_g1,
      t_vec_g2 = t_vec_2_g2,
      beta0_new = beta0_2_star,
      beta0 = beta0_2,
      beta1 = beta1_2,
      eta_g1 = eta2_g1,
      eta_g2 = eta2_g2,
      a_t = a_2,
      b_t = b_2,
      mu0 = mu0_2,
      sig0 = sig0_2,
      link
    )

    if (log(runif(1)) < log_rho_beta0_2) {
      beta0_2 <- beta0_2_star
      keep.beta0.2 <- keep.beta0.2 + 1
      keep.tmp.beta0.2 <- keep.tmp.beta0.2 + 1
    }

    # ---------
    # beta1_1: M-H
    # ---------
    beta1_1_star <- rnorm(1, beta1_1, tune.beta1.1)

    log_rho_beta1_1 <- metro_log_ratio_b1(
      t_vec_g2 = t_vec_1_g2,
      beta0 = beta0_1,
      beta1_new = beta1_1_star,
      beta1 = beta1_1,
      eta_g2 = eta1_g2,
      a_t = a_1, b_t = b_1,
      mu1 = mu1_1,
      sig1 = sig1_1,
      link = link
    )

    if (log(runif(1)) < log_rho_beta1_1) {
      beta1_1 <- beta1_1_star
      keep.beta1.1 <- keep.beta1.1 + 1
      keep.tmp.beta1.1 <- keep.tmp.beta1.1 + 1
    }

    # ---------
    # beta1_2: M-H
    # ---------
    beta1_2_star <- rnorm(1, beta1_2, tune.beta1.2)

    log_rho_beta1_2 <- metro_log_ratio_b1(
      t_vec_g2 = t_vec_2_g2,
      beta0 = beta0_2,
      beta1_new = beta1_2_star,
      beta1 = beta1_2,
      eta_g2 = eta1_g2,
      a_t = a_2, b_t = b_2,
      mu1 = mu1_2,
      sig1 = sig1_2,
      link = link
    )

    if (log(runif(1)) < log_rho_beta1_2) {
      beta1_2 <- beta1_2_star
      keep.beta1.2 <- keep.beta1.2 + 1
      keep.tmp.beta1.2 <- keep.tmp.beta1.2 + 1
    }

    if (link == "logit") {
      r1_g1 <- 1 / (1 + exp(-beta0_1))
      r2_g1 <- 1 / (1 + exp(-beta0_2))
      r1_g2 <- 1 / (1 + exp(-(beta0_1 + beta1_1)))
      r2_g2 <- 1 / (1 + exp(-(beta0_2 + beta1_2)))
    }

    if (link == "probit") {
      r1_g1 <- pnorm(beta0_1)
      r2_g1 <- pnorm(beta0_2)
      r1_g2 <- pnorm(beta0_1 + beta1_1)
      r2_g2 <- pnorm(beta0_2 + beta1_2)
    }

    if (link == "cloglog") {
      r1_g1 <- 1 - exp(-exp(beta0_1))
      r2_g1 <- 1 - exp(-exp(beta0_2))
      r1_g2 <- 1 - exp(-exp(beta0_1 + beta1_1))
      r2_g2 <- 1 - exp(-exp(beta0_2 + beta1_2))

    }

    # ====================================================================
    # Sampling eta given (t, r, theta^(i))
    # ====================================================================

    # eta1_g1
    xi1_g1_star <- rnorm(1, xi1_g1, tune.g1.eta.1)

    log_rho_1_g1_eta <- metro_log_ratio_eta(
      t_vec = t_vec_1_g1,
      # beta0 = beta0_1,
      # beta1 = 0,
      r = r1_g1,
      xi_new = xi1_g1_star,
      xi = xi1_g1,
      a = a_1,
      b = b_1,
      a_eta = a_eta_g1_1,
      b_eta = b_eta_g1_1)

    if (log(runif(1)) < log_rho_1_g1_eta) {
      xi1_g1 <- xi1_g1_star
      eta1_g1 <- exp(xi1_g1)
      keep.g1.eta.1 <- keep.g1.eta.1 + 1
      keep.tmp.g1.eta.1 <- keep.tmp.g1.eta.1 + 1
    }

    # eta2_g1
    xi2_g1_star <- rnorm(1, xi2_g1, tune.g1.eta.2)
    log_rho_2_g1_eta <- metro_log_ratio_eta(
      t_vec = t_vec_2_g1,
      r = r2_g1,
      # beta0 = beta0_2,
      # beta1 = 0,
      xi_new = xi2_g1_star,
      xi = xi2_g1,
      a = a_2,
      b = b_2,
      a_eta = a_eta_g1_2,
      b_eta = b_eta_g1_2)

    if (log(runif(1)) < log_rho_2_g1_eta) {
      xi2_g1 <- xi2_g1_star
      eta2_g1 <- exp(xi2_g1)
      keep.g1.eta.2 <- keep.g1.eta.2 + 1
      keep.tmp.g1.eta.2 <- keep.tmp.g1.eta.2 + 1
    }

    # eta1_g2
    xi1_g2_star <- rnorm(1, xi1_g2, tune.g2.eta.1)
    log_rho_1_g2_eta <- metro_log_ratio_eta(
      t_vec = t_vec_1_g2,
      r = r1_g2,
      xi_new = xi1_g2_star,
      xi = xi1_g2,
      a = a_1,
      b = b_1,
      a_eta = a_eta_g2_1,
      b_eta = b_eta_g2_1)

    if (log(runif(1)) < log_rho_1_g2_eta) {
      xi1_g2 <- xi1_g2_star
      eta1_g2 <- exp(xi1_g2)
      keep.g2.eta.1 <- keep.g2.eta.1 + 1
      keep.tmp.g2.eta.1 <- keep.tmp.g2.eta.1 + 1
    }

    # eta2_g2
    xi2_g2_star <- rnorm(1, xi2_g2, tune.g2.eta.2)
    log_rho_2_g2_eta <- metro_log_ratio_eta(
      t_vec = t_vec_2_g2,
      r = r2_g2,
      xi_new = xi2_g2_star,
      xi = xi2_g2,
      a = a_2,
      b = b_2,
      a_eta = a_eta_g2_2,
      b_eta = b_eta_g2_2)
    if (log(runif(1)) < log_rho_2_g2_eta) {
      xi2_g2 <- xi2_g2_star
      eta2_g2 <- exp(xi2_g2)
      keep.g2.eta.2 <- keep.g2.eta.2 + 1
      keep.tmp.g2.eta.2 <- keep.tmp.g2.eta.2 + 1
    }


    # ====================================================================
    # Sampling sigma^2 given (t, r, theta^(i))
    # ====================================================================

    sig <- sqrt(sample_sig2_slam(multi_y_g1 = multi_y_g1,
                                 multi_y_g2 = multi_y_g2,
                                 x_vec = x_vec, theta = theta_k,
                                 Sigma_g1 = Sigma_lst_g1,
                                 Sigma_g2 = Sigma_lst_g2,
                                 ga_shape = ga_shape,
                                 ga_rate = ga_rate))

    #  Save samples
    # -----------------------------------------------------
    if (i > burn) {
      if (i %% thin == 0) {
        draws[(i - burn) %/% thin, ] <- c(t_mat_g1[1, ],
                                          t_mat_g1[2, ],
                                          t_mat_g2[1, ],
                                          t_mat_g2[2, ],
                                          beta0_1,
                                          beta0_2,
                                          beta1_1,
                                          beta1_2,
                                          r1_g1,
                                          r2_g1,
                                          r1_g2,
                                          r2_g2,
                                          eta1_g1,
                                          eta2_g1,
                                          eta1_g2,
                                          eta2_g2,
                                          sig)
      }
    }
  }

  #-----------------------------
  # Acceptance Probability
  #----------------------------

  keep.lst.g1 <- lapply(keep.lst.g1, function(x) x / n.mcmc)
  keep.lst.g2 <- lapply(keep.lst.g2, function(x) x / n.mcmc)

  keep.g1.eta.1 <- keep.g1.eta.1 / n.mcmc
  keep.g1.eta.2 <- keep.g1.eta.2 / n.mcmc
  keep.g2.eta.1 <- keep.g2.eta.1 / n.mcmc
  keep.g2.eta.2 <- keep.g2.eta.2 / n.mcmc

  keep.beta0.1 <- keep.beta0.1 / n.mcmc
  keep.beta0.2 <- keep.beta0.2 / n.mcmc
  keep.beta1.1 <- keep.beta1.1 / n.mcmc
  keep.beta1.2 <- keep.beta1.2 / n.mcmc


  #--------------
  # Write output
  #--------------
  # print(draws)
  return(list(draws = draws,
              accept_g1 = keep.lst.g1,
              accept_g2 = keep.lst.g2,
              accept_g1_eta_1 = keep.g1.eta.1,
              accept_g1_eta_2 = keep.g1.eta.2,
              accept_g2_eta_1 = keep.g2.eta.1,
              accept_g2_eta_2 = keep.g2.eta.2,
              accept_beta0_1 = keep.beta0.1,
              accept_beta0_2 = keep.beta0.2,
              accept_beta1_1 = keep.beta1.1,
              accept_beta1_2 = keep.beta1.2,
              start = start.lst,
              tune_g1 = tune.lst.g1,
              tune_g2 = tune.lst.g2,
              tune_g1_eta_1 = tune.g1.eta.1,
              tune_g1_eta_2 = tune.g1.eta.2,
              tune_g2_eta_1 = tune.g2.eta.1,
              tune_g2_eta_2 = tune.g2.eta.2,
              tune_beta0_1 = tune.beta0.1,
              tune_beta0_2 = tune.beta0.2,
              tune_beta1_1 = tune.beta1.1,
              tune_beta1_2 = tune.beta1.2,
              burn = burn,
              thin = thin,
              n.mcmc = n.mcmc,
              sampleidx = sampleidx))
}
