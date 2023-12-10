get_pi_t_sig <- function(yJ, x, x_test, idx_der, sample_sig, tau, h) {

  len_test <- length(x_test)
  mm <- dim(idx_der)[1]

  mean_var_lst <- lapply(1:mm, function(i) {
    tau_new <- tau*sample_sig[i]
    Kff <- compute_cov_1d(idx1 = x, idx2 = x, tau = tau_new, h = h)
    Kffnew <- compute_cov_1d(idx1 = x, idx2 = x_test, tau = tau_new, h = h)
    Kfnewfnew <- compute_cov_1d(idx1 = x_test, idx2 = x_test, tau = tau_new, h = h)
    K_joint <- compute_joint_cov(idx_obs = x, idx_der = idx_der[i, ],
                                 sig = sample_sig[i], tau = tau_new, h = h)
    Kdfnew <- computeCovDer1(idx1 = idx_der[i, ], idx2 = x_test,
                             tau = tau_new, h = h)
    return(cond_norm_mean_var(mean_vec_1 = 0, mean_vec_2 = 0, obs_2 = yJ,
                              cov_mat_1 = Kfnewfnew, cov_mat_2 = K_joint,
                              cov_mat_12 = t(rbind(Kffnew, Kdfnew))))
  })

  pred_f_der_mat <- sapply(mean_var_lst, function(x) {
    mvnfast::rmvn(1, mu = x$condMean,
                  sigma = x$condVar)})

  PI_f_der <- apply(pred_f_der_mat, 1, quantile, prob = c(0.025, 0.975))
  mu_test_der_mean <- apply(pred_f_der_mat, 1, mean)

  return(list(pred_f = pred_f_der_mat, ci_low = PI_f_der[1, ],
              ci_high = PI_f_der[2, ],
              mu_test = mu_test_der_mean))
}
