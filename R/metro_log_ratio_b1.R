# young: g1
# old: g2

metro_log_ratio_b1 <- function(t_vec_g2,
                                       beta0, beta1_new, beta1,
                                       eta_g2,
                                       a_t, b_t,
                                       mu1, sig1, link) {

  # conditional of beta0 (Metropolis ratio)
  # ============================================

  r_g2_new <- 1 / (1 + exp(-(beta0 + beta1_new)))
  r_g2 <- 1 / (1 + exp(-(beta0 + beta1)))


  if (link == "logit") {
    r_g2_new <- 1 / (1 + exp(-(beta0 + beta1_new)))
    r_g2 <- 1 / (1 + exp(-(beta0 + beta1)))
  }

  if (link == "probit") {
    r_g2_new <- pnorm(beta0 + beta1_new)
    r_g2 <- pnorm(beta0 + beta1)
  }

  if (link == "cloglog") {
    r_g2_new <- 1 - exp(-exp(beta0 + beta1_new))
    r_g2 <- 1 - exp(-exp(beta0 + beta1))
  }

  shape1_g2_new <- r_g2_new * eta_g2
  shape2_g2_new <- (1 - r_g2_new) * eta_g2

  shape1_g2 <- r_g2 * eta_g2
  shape2_g2 <- (1 - r_g2) * eta_g2

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

  log_normal_beta1_new <- dnorm(x = beta1_new, mean = mu1, sd = sig1, log = TRUE)
  log_normal_beta1 <- dnorm(x = beta1, mean = mu1, sd = sig1, log = TRUE)

  A_new <- sum_log_gbeta_t_g2_new + log_normal_beta1_new

  A <- sum_log_gbeta_t_g2 + log_normal_beta1

  return(A_new - A)
}
