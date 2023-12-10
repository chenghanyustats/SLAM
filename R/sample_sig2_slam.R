# young: g1
# old: g2

sample_sig2_slam <- function(multi_y_g1, multi_y_g2, x_vec,
                             theta, Sigma_g1, Sigma_g2,
                             ga_shape = 5, ga_rate = 5) {
  n_subj_g1 <- ncol(multi_y_g1)
  n <- nrow(multi_y_g1)

  n_subj_g2 <- ncol(multi_y_g2)

  quadinv_lst_g1 <- lapply(1:n_subj_g1, function(s) {
    quad.form.inv(Sigma_g1[[s]] + diag(n), multi_y_g1[, s])
  })

  quadinv_lst_g2 <- lapply(1:n_subj_g2, function(s) {
    quad.form.inv(Sigma_g2[[s]] + diag(n), multi_y_g2[, s])
  })

  quadinvsum_g1 <- sum(unlist(quadinv_lst_g1))
  quadinvsum_g2 <- sum(unlist(quadinv_lst_g2))

  invgamma::rinvgamma(1, shape = ga_shape + (n * (n_subj_g1 + n_subj_g2)) / 2,
                      rate = ga_rate + (1 / 2) * (quadinvsum_g1 + quadinvsum_g2))
}
