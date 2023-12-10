
## mean peak over range of latency r (Mean Peak)
get_amp_mean <- function(pred_f, x_test, lat_sample_1, lat_sample_2,
                         is.abs = FALSE) {
  mm <- length(lat_sample_1)
  amplitude_sample <- matrix(0, mm, 2)

  a1 <- find_closest(x_test, range(lat_sample_1)[1])
  b1 <- find_closest(x_test, range(lat_sample_1)[2])
  for (m in 1:mm) {
    if(is.abs) {
      pred_r1 <- abs(pred_f[a1:b1, m])
    } else {
      pred_r1 <- pred_f[a1:b1, m]
    }
    amplitude_sample[m, 1] <- mean(pred_r1)
  }

  a2 <- find_closest(x_test, range(lat_sample_2)[1])
  b2 <- find_closest(x_test, range(lat_sample_2)[2])
  for (m in 1:mm) {
    if(is.abs) {
      pred_r2 <- abs(pred_f[a2:b2, m])
    } else {
      pred_r2 <- pred_f[a2:b2, m]
    }
    amplitude_sample[m, 2] <- mean(pred_r2)
  }
  return(amplitude_sample)
}
