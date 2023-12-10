
# g1: g1
# g2: g2

log_mar_lik_slam <- function(theta,
                             multi_y_g1,
                             multi_y_g2,
                             x_vec,
                             der_mc_mat_g1_lst,
                             der_mc_mat_g2_lst,
                             H0,
                             sample_sig,
                             a_h = 1, b_h = 1,
                             is.h.par = FALSE) {

  ## der_mc_mat_g1_lst: a list of length n_subj_g1, with element
  ##                       a matrix of dim 2 x D

  # nn <- nrow(der_mc_mat)
  D_g1 <- ncol(der_mc_mat_g1_lst[[1]])
  D_g2 <- ncol(der_mc_mat_g2_lst[[1]])
  nn <- 2
  n <- nrow(multi_y_g1)
  n_subj_g1 <- ncol(multi_y_g1)
  n_subj_g2 <- ncol(multi_y_g2)

  mu <- rep(0L, n)
  Kff <- se_ker(H0 = H0, tau = theta[1], h = theta[2])
  # Psi <- Sigma + diag(theta[1] ^ 2, n)
  # Psi <- (t(Psi) + Psi) / 2


  log_all_g1 <- sapply(1:n_subj_g1, function(s) {
    # dd_g1 <- t(der_mc_mat_g1[, c(s, s + n_subj_g1)])
    Kdf_lst <- as.list(data.frame(apply(der_mc_mat_g1_lst[[s]], 2,
                                        computeCovDer1,
                                        idx2 = x_vec,
                                        tau = theta[1], h = theta[2])))
    Kdf_lst_new <- lapply(Kdf_lst, function(x) {
      matrix(x, nrow = nn, ncol = n)
    })

    Kdd_lst <- as.list(data.frame(apply(der_mc_mat_g1_lst[[s]], 2,
                                        computeCovDer2,
                                        tau = theta[1], h = theta[2])))
    Kdd_lst_new <- lapply(Kdd_lst, function(x){
      matrix(x, nrow = nn, ncol = nn)
    })


    den_value_lst_new <- lapply(1:D_g1, function(m) {
      Sigma <- Kff - quad.form.inv(Kdd_lst_new[[m]], Kdf_lst_new[[m]])
      # if(!is.positive.definite(Psi)) Psi <- nearPD(Psi)$mat
      Psi <- sample_sig[m] ^ 2 * (Sigma + diag(n + 1e-4))
      Psi <- (Psi + t(Psi)) / 2
      mvnfast::dmvn(multi_y_g1[, s], mu = mu, sigma = Psi)
    })


    return(log(sum(unlist(den_value_lst_new)) / D_g1))
  })


  log_all_g2 <- sapply(1:n_subj_g2, function(s) {
    # dd_g2 <- t(der_mc_mat_g2[, c(s, s + n_subj_g2)])
    Kdf_lst <- as.list(data.frame(apply(der_mc_mat_g2_lst[[s]], 2,
                                        computeCovDer1,
                                        idx2 = x_vec,
                                        tau = theta[1], h = theta[2])))
    Kdf_lst_new <- lapply(Kdf_lst, function(x) {
      matrix(x, nrow = nn, ncol = n)
    })

    Kdd_lst <- as.list(data.frame(apply(der_mc_mat_g2_lst[[s]], 2,
                                        computeCovDer2,
                                        tau = theta[1], h = theta[2])))
    Kdd_lst_new <- lapply(Kdd_lst, function(x){
      matrix(x, nrow = nn, ncol = nn)
    })


    den_value_lst_new <- lapply(1:D_g2, function(m) {
      Sigma <- Kff - quad.form.inv(Kdd_lst_new[[m]], Kdf_lst_new[[m]])
      Psi <- sample_sig[m] ^ 2 * (Sigma + diag(n + 1e-4))
      Psi <- (Psi + t(Psi)) / 2
      # if(!is.positive.definite(Psi)) Psi <- nearPD(Psi)$mat
      mvnfast::dmvn(multi_y_g2[, s], mu = mu, sigma = Psi)
    })

    return(log(sum(unlist(den_value_lst_new)) / D_g2))
  })

  log_summ <- sum(log_all_g1) + sum(log_all_g2)

  if (is.h.par) {
    C <- dgamma(theta[2], shape = a_h, rate = b_h, log = TRUE)
    return(-log_summ - C)
  } else {
    return(-log_summ)
  }
}
