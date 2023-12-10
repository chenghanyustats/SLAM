
get_amp_at_mean_lat <- function(pred_f, x_test, r1_sample, r2_sample) {
  ## fix r at posterior mean
  mm <- length(r1_sample)
  # amplitude_sample <- matrix(0, mm, 2)
  amplitude_sample_mean <- matrix(0, mm, 2)
  r1_post_mean <- mean(r1_sample)
  r2_post_mean <- mean(r2_sample)
  r1_mean_idx <- find_closest(x_test, r1_post_mean)
  r2_mean_idx <- find_closest(x_test, r2_post_mean)
  for (m in 1:mm) {
    # r1_idx <- find_closest(x_test, r1_sample[m])
    # r2_idx <- find_closest(x_test, r2_sample[m])
    # amplitude_sample[m, 1] <- abs(pred_f[r1_idx, m])
    # amplitude_sample[m, 2] <- abs(pred_f[r2_idx, m])
    amplitude_sample_mean[m, 1] <- abs(pred_f[r1_mean_idx, m])
    amplitude_sample_mean[m, 2] <- abs(pred_f[r2_mean_idx, m])
  }
  return(amplitude_sample_mean)
}
