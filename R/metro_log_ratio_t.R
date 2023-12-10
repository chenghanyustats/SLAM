metro_log_ratio_t <- function(ys, t_vec, ts_new,
                              Sigma_s, Sigma_s_new,
                              sig, r1, r2, eta1, eta2,
                              a1, b1, a2, b2,
                              sd_trun = NULL) {
  # conditional of ts (Metropolis ratio)
  # ============================================
  n <- length(ys)

  Psi <- (diag(n) + Sigma_s) * sig ^ 2
  Psi <- (Psi + t(Psi)) / 2
  if (!matrixcalc::is.positive.definite(Psi)) {
    Psi <- as.matrix(Matrix::nearPD(Psi)$mat)
  }

  Psi_new <- (diag(n) + Sigma_s_new) * sig ^ 2
  Psi_new <- (Psi_new + t(Psi_new)) / 2
  if (!matrixcalc::is.positive.definite(Psi_new)) {
    Psi <- as.matrix(Matrix::nearPD(Psi_new)$mat)
  }

  mu <- rep(0L, n)

  log_normal <- mvnfast::dmvn(ys, mu = mu, sigma = Psi, log = TRUE)
  log_normal_new <- mvnfast::dmvn(ys, mu = mu, sigma = Psi_new,
                                  log = TRUE)
  shape1_1 <- r1 * eta1
  shape2_1 <- (1 - r1) * eta1

  shape1_2 <- r2 * eta2
  shape2_2 <- (1 - r2) * eta2

  log_gbeta_t_new_1 <- log_dgbeta(x = ts_new[1], shape1 = shape1_1,
                                  shape2 = shape2_1, a = a1, b = b1)
  log_gbeta_t_new_2 <- log_dgbeta(x = ts_new[2], shape1 = shape1_2,
                                  shape2 = shape2_2, a = a2, b = b2)

  log_trun_t_new_1 <-
    log(truncnorm::dtruncnorm(x = ts_new[1], mean = t_vec[1], sd = sd_trun,
                              a = a1, b = b1))
  log_trun_t_new_2 <-
    log(truncnorm::dtruncnorm(x = ts_new[2], mean = t_vec[2], sd = sd_trun,
                              a = a2, b = b2))

  log_gbeta_t_1 <- log_dgbeta(x = t_vec[1], shape1 = shape1_1,
                              shape2 = shape2_1, a = a1, b = b1)
  log_gbeta_t_2 <- log_dgbeta(x = t_vec[2], shape1 = shape1_2,
                              shape2 = shape2_2, a = a2, b = b2)

  log_trun_t_1 <-
    log(truncnorm::dtruncnorm(x = t_vec[1], mean = ts_new[1], sd = sd_trun,
                              a = a1, b = b1))
  log_trun_t_2 <-
    log(truncnorm::dtruncnorm(x = t_vec[2], mean = ts_new[2], sd = sd_trun,
                              a = a2, b = b2))

  A_new <- log_normal_new + log_gbeta_t_new_1 + log_gbeta_t_new_2 -
    log_trun_t_new_1 - log_trun_t_new_2
  A <- log_normal + log_gbeta_t_1 + log_gbeta_t_2 - log_trun_t_1 - log_trun_t_2
  return(A_new - A)
}
