# The algorithm fixes the number of groups at two and
# assumes known number of t's (M = 2)

# young: g1
# old: g2

mcem_slam <- function(multi_y_g1, multi_y_g2, x, H0, a_1, b_1, a_2, b_2,
                      start.lst, name.par,
                      theta_init = c(1, 1),
                      epsilon = 1e-3, D = 100, max_iter = 100,
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
                      adapt = TRUE,
                      tune.lst.g1 = rep(list(0.5), ncol(multi_y_g1)),
                      keep.lst.g1 = rep(list(0), ncol(multi_y_g1)),
                      tune.lst.g2 = rep(list(0.5), ncol(multi_y_g2)),
                      keep.lst.g2 = rep(list(0), ncol(multi_y_g2)),
                      keep.g1.eta.1 = 0,
                      keep.g1.eta.2 = 0,
                      keep.g2.eta.1 = 0,
                      keep.g2.eta.2 = 0,
                      tune.g1.eta.1 = 2,
                      tune.g1.eta.2 = 2,
                      tune.g2.eta.1 = 2,
                      tune.g2.eta.2 = 2,
                      keep.beta0.1 = 0,
                      keep.beta0.2 = 0,
                      keep.beta1.1 = 0,
                      keep.beta1.2 = 0,
                      tune.beta0.1 = 2,
                      tune.beta0.2 = 2,
                      tune.beta1.1 = 2,
                      tune.beta1.2 = 2,
                      tune.len = 40,
                      target.accept.rate = 0.35,
                      a_h = 1, b_h = 1,
                      mu0_1 = 0,
                      mu0_2 = 0,
                      mu1_1 = 0,
                      mu1_2 = 0,
                      sig0_1 = 1,
                      sig0_2 = 1,
                      sig1_1 = 1,
                      sig1_2 = 1,
                      a_eta_g1_1 = 1 / 2,
                      b_eta_g1_1 = 1 / 2,
                      a_eta_g1_2 = 1 / 2,
                      b_eta_g1_2 = 1 / 2,
                      a_eta_g2_1 = 1 / 2,
                      b_eta_g2_1 = 1 / 2,
                      a_eta_g2_2 = 1 / 2,
                      b_eta_g2_2 = 1 / 2,
                      final_draw = TRUE,
                      link = "logit") {

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
  n_subj_g1 <- ncol(multi_y_g1)
  n_subj_g2 <- ncol(multi_y_g2)
  n_subj <- n_subj_g2 + n_subj_g1
  n_t <- n_subj * 2

  eps <- epsilon + 1
  count <- 1

  # Initialization of parameters
  theta_mat <- matrix(0L, max_iter, n_theta)
  theta_mat[1, ] <- theta_k <- theta_init
  theta_new <- theta_k
  mar_post_new <- 0

  ## =========================================================================
  ## Multi-subject (g1 and g2) MCEM Algorithm for SLAM
  ## =========================================================================
  print("MCEM SLAM begins")
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
      start.lst <- list(t_g1 = matrix(runif(n_subj_g1 * 2,
                                               min = c(a_1, a_2),
                                               max = c(b_1, b_2)),
                                         2, n_subj_g1),
                        t_g2 = matrix(runif(n_subj_g2 * 2,
                                             min = c(a_1, a_2),
                                             max = c(b_1, b_2)),
                                       2, n_subj_g2),
                        beta0_1 = post_mean[n_t + 1],
                        beta0_2 = post_mean[n_t + 2],
                        beta1_1 = post_mean[n_t + 3],
                        beta1_2 = post_mean[n_t + 4],
                        eta1_g1 = post_mean[n_t + 5],
                        eta2_g1 = post_mean[n_t + 6],
                        eta1_g2 = post_mean[n_t + 7],
                        eta2_g2 = post_mean[n_t + 8],
                        sigma = post_mean[n_t + 9])

    }

    print("E-step sampling t, beta, eta and sigma")
    mcmc_output <- Estep_one_way(
      multi_y_g1 = multi_y_g1,
      multi_y_g2 = multi_y_g2,
      x_vec = x,
      n_subj_g1 = n_subj_g1,
      n_subj_g2 = n_subj_g2,
      theta_k = theta_k, H0 = H0,
      start.lst = start.lst,
      n.mcmc = n_mcmc_e,
      burn = burn_e,
      thin = thin_e,
      name.par = name.par,
      adapt = adapt,
      tune.lst.g1 = tune.lst.g1,
      keep.lst.g1 = keep.lst.g1,
      tune.lst.g2 = tune.lst.g2,
      keep.lst.g2 = keep.lst.g2,
      keep.g1.eta.1 = keep.g1.eta.1,
      keep.g1.eta.2 = keep.g1.eta.2,
      keep.g2.eta.1 = keep.g2.eta.1,
      keep.g2.eta.2 = keep.g2.eta.2,
      tune.g1.eta.1 = tune.g1.eta.1,
      tune.g1.eta.2 = tune.g1.eta.2,
      tune.g2.eta.1 = tune.g2.eta.1,
      tune.g2.eta.2 = tune.g2.eta.2,
      keep.beta0.1 = keep.beta0.1,
      keep.beta0.2 = keep.beta0.2,
      keep.beta1.1 = keep.beta1.1,
      keep.beta1.2 = keep.beta1.2,
      tune.beta0.1 = tune.beta0.1,
      tune.beta0.2 = tune.beta0.2,
      tune.beta1.1 = tune.beta1.1,
      tune.beta1.2 = tune.beta1.2,
      tune.len = tune.len,
      target.accept.rate = target.accept.rate,
      a_1 = a_1,
      b_1 = b_1,
      a_2 = a_2,
      b_2 = b_2,
      a_eta_g1_1 = a_eta_g1_1,
      b_eta_g1_1 = b_eta_g1_1,
      a_eta_g1_2 = a_eta_g1_2,
      b_eta_g1_2 = b_eta_g1_2,
      a_eta_g2_1 = a_eta_g2_1,
      b_eta_g2_1 = b_eta_g2_1,
      a_eta_g2_2 = a_eta_g2_2,
      b_eta_g2_2 = b_eta_g2_2,
      ga_shape = ga_shape, ga_rate = ga_rate,
      mu0_1 = mu0_1, sig0_1 = sig0_1,
      mu0_2 = mu0_2, sig0_2 = sig0_2,
      mu1_1 = mu1_1, sig1_1 = sig1_1,
      mu1_2 = mu1_2, sig1_2 = sig1_2,
      link = link)

    sample_t_g1 <- mcmc_output$draws[, 1:(2*n_subj_g1)]
    sample_t_g2 <- mcmc_output$draws[, (2*n_subj_g1 + 1):(2*n_subj_g1 + 2*n_subj_g2)]
    sample_sig <- mcmc_output$draws[, ncol(mcmc_output$draws)]

    sample_t_g1_mstep <- apply(sample_t_g1, 2, sample, size = D)
    sample_t_g2_mstep <- apply(sample_t_g2, 2, sample, size = D)
    sample_sig_mstep <- sample(sample_sig, size = D)


    sample_t_g1_lst <- format_for_M_step(sample_t_g1_mstep)
    sample_t_g2_lst <- format_for_M_step(sample_t_g2_mstep)


    # *****************************************
    # ******** M-step for theta ************* #
    # *****************************************
    print("M-step updating kernel parameters")

    res <- Rsolnp::solnp(pars = theta_k,
                         fun = log_mar_lik_slam,
                         LB = lower, control = ctrl,
                         multi_y_g1 = multi_y_g1,
                         multi_y_g2 = multi_y_g2,
                         x_vec = x, H0 = H0,
                         der_mc_mat_g1_lst = sample_t_g1_lst,
                         der_mc_mat_g2_lst = sample_t_g2_lst,
                         sample_sig = sample_sig_mstep,
                         is.h.par = is.h.par,
                         a_h = a_h, b_h = b_h)
    theta_k <- res$par
    mar_post_k <- res$values[length(res$values)]

    print(paste("theta =", theta_k))
    print(paste("marg post =", -mar_post_k))

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
    print("Final sampling for t, beta, eta, sigma")
    mcmc_output <- Estep_one_way(
      multi_y_g1 = multi_y_g1,
      multi_y_g2 = multi_y_g2,
      x_vec = x,
      n_subj_g1 = n_subj_g1,
      n_subj_g2 = n_subj_g2,
      theta_k = theta_k, H0 = H0,
      start.lst = start.lst,
      n.mcmc = n_mcmc_final,
      burn = burn_final,
      thin = thin_final,
      name.par = name.par, adapt = adapt,
      tune.lst.g1 = tune.lst.g1,
      keep.lst.g1 = keep.lst.g1,
      tune.lst.g2 = tune.lst.g2,
      keep.lst.g2 = keep.lst.g2,
      keep.g1.eta.1 = keep.g1.eta.1,
      keep.g1.eta.2 = keep.g1.eta.2,
      keep.g2.eta.1 = keep.g2.eta.1,
      keep.g2.eta.2 = keep.g2.eta.2,
      tune.g1.eta.1 = tune.g1.eta.1,
      tune.g1.eta.2 = tune.g1.eta.2,
      tune.g2.eta.1 = tune.g2.eta.1,
      tune.g2.eta.2 = tune.g2.eta.2,
      keep.beta0.1 = keep.beta0.1,
      keep.beta0.2 = keep.beta0.2,
      keep.beta1.1 = keep.beta1.1,
      keep.beta1.2 = keep.beta1.2,
      tune.beta0.1 = tune.beta0.1,
      tune.beta0.2 = tune.beta0.2,
      tune.beta1.1 = tune.beta1.1,
      tune.beta1.2 = tune.beta1.2,
      tune.len = tune.len,
      target.accept.rate = target.accept.rate,
      a_1 = a_1,
      b_1 = b_1,
      a_2 = a_2,
      b_2 = b_2,
      a_eta_g1_1 = a_eta_g1_1,
      b_eta_g1_1 = b_eta_g1_1,
      a_eta_g1_2 = a_eta_g1_2,
      b_eta_g1_2 = b_eta_g1_2,
      a_eta_g2_1 = a_eta_g2_1,
      b_eta_g2_1 = b_eta_g2_1,
      a_eta_g2_2 = a_eta_g2_2,
      b_eta_g2_2 = b_eta_g2_2,
      ga_shape = ga_shape, ga_rate = ga_rate,
      mu0_1 = mu0_1, sig0_1 = sig0_1,
      mu0_2 = mu0_2, sig0_2 = sig0_2,
      mu1_1 = mu1_1, sig1_1 = sig1_1,
      mu1_2 = mu1_2, sig1_2 = sig1_2,
      link = link)
  }

  cat("\n")
  print("Done!")
  colnames(theta_mat) <- c("tau", "h")

  return(list(mcmc_output = mcmc_output,
              theta_mat = theta_mat[1:count, ]))

}
