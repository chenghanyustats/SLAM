metro_log_ratio_eta <- function(t_vec, r, xi_new, xi, a, b,
                                                 a_eta, b_eta) {
  # conditional of r (Metropolis ratio)
  # ============================================
  if (r == 1) {
    r <- 1 - 1e-6
  }
  if (r == 0) {
    r <- 1e-6
  }

  eta_new <- exp(xi_new)
  eta <- exp(xi)
  shape1_new <- r * eta_new
  shape2_new <- (1 - r) * eta_new
  shape1 <- r * eta
  shape2 <- (1 - r) * eta

  sum_log_gbeta_t_new <- sum(log_dgbeta(x = t_vec, shape1 =  shape1_new,
                                        shape2 = shape2_new,
                                        a = a, b = b))
  sum_log_gbeta_t <- sum(log_dgbeta(x = t_vec, shape1 = shape1,
                                    shape2 = shape2,
                                    a = a, b = b))

  log_ga_eta_new <- stats::dgamma(x = eta_new, shape = a_eta,
                                  rate = b_eta, log = TRUE)
  log_ga_eta <- stats::dgamma(x = eta, shape = a_eta,
                              rate = b_eta, log = TRUE)

  A_new <- sum_log_gbeta_t_new + log_ga_eta_new

  A <- sum_log_gbeta_t + log_ga_eta

  return(A_new - A + xi_new - xi)
}
