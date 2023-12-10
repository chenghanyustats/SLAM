format_for_M_step <- function(sample_t) {
  sample_t_lst <- list()
  D <- nrow(sample_t)
  S <- ncol(sample_t) / 2
  for (s in 1:S) {
    sample_t_lst[[s]] <- t(sample_t[, c(s, s + S)])
  }
  return(sample_t_lst)
}
