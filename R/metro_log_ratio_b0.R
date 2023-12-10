
# young: g1
# old: g2

metro_log_ratio_b0 <- function(t_vec_g1, t_vec_g2,
                               beta0_new, beta0, beta1,
                               eta_g1, eta_g2,
                               a_t, b_t,
                               mu0, sig0,
                               link) {

  # conditional of beta0 (Metropolis ratio)
  # ============================================

  if (link == "logit") {
    r_g1_new <- 1 / (1 + exp(-beta0_new))
    r_g1 <- 1 / (1 + exp(-beta0))
    r_g2_new <- 1 / (1 + exp(-(beta0_new + beta1)))
    r_g2 <- 1 / (1 + exp(-(beta0 + beta1)))
  }

  if (link == "probit") {
    r_g1_new <- pnorm(beta0_new)
    r_g1 <- pnorm(beta0)
    r_g2_new <- pnorm(beta0_new + beta1)
    r_g2 <- pnorm(beta0 + beta1)
  }

  if (link == "cloglog") {
    r_g1_new <- 1 - exp(-exp(beta0_new))
    r_g1 <- 1 - exp(-exp(beta0))
    r_g2_new <- 1 - exp(-exp(beta0_new + beta1))
    r_g2 <- 1 - exp(-exp(beta0 + beta1))
  }

  shape1_g1_new <- r_g1_new * eta_g1
  shape2_g1_new <- (1 - r_g1_new) * eta_g1

  shape1_g1 <- r_g1 * eta_g1
  shape2_g1 <- (1 - r_g1) * eta_g1

  shape1_g2_new <- r_g2_new * eta_g2
  shape2_g2_new <- (1 - r_g2_new) * eta_g2

  shape1_g2 <- r_g2 * eta_g2
  shape2_g2 <- (1 - r_g2) * eta_g2



  sum_log_gbeta_t_g1_new <- sum(log_dgbeta(x = t_vec_g1,
                                              shape1 = shape1_g1_new,
                                              shape2 = shape2_g1_new,
                                              a = a_t,
                                              b = b_t))
  sum_log_gbeta_t_g1 <- sum(log_dgbeta(x = t_vec_g1,
                                          shape1 = shape1_g1,
                                          shape2 = shape2_g1,
                                          a = a_t,
                                          b = b_t))

  sum_log_gbeta_t_g2_new <- sum(log_dgbeta(x = t_vec_g2,
                                            shape1 = shape1_g2_new,
                                            shape2 = shape2_g2_new,
                                            a = a_t,
                                            b = b_t))
  sum_log_gbeta_t_g2 <- sum(log_dgbeta(x = t_vec_g2,
                                        shape1 = shape1_g2,
                                        shape2 = shape2_g2,
                                        a = a_t,
                                        b = b_t))
  sum_log_gbeta_t_new <- sum_log_gbeta_t_g1_new + sum_log_gbeta_t_g2_new
  sum_log_gbeta_t <- sum_log_gbeta_t_g1 + sum_log_gbeta_t_g2


  log_normal_beta0_new <- dnorm(x = beta0_new, mean = mu0, sd = sig0, log = TRUE)
  log_normal_beta0 <- dnorm(x = beta0, mean = mu0, sd = sig0, log = TRUE)

  A_new <- sum_log_gbeta_t_new + log_normal_beta0_new

  A <- sum_log_gbeta_t + log_normal_beta0

  return(A_new - A)
}
