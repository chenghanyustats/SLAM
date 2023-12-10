
## peak at mid latency r (Integral Peak)
get_amp_integrate <- function(pred_f, x_test, lat_sample_1, lat_sample_2,
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

    all_sum_pred_1 <- sum(pred_r1)
    sum_pred_1 <- pred_r1[1]
    for (i in 2:length(pred_r1)) {
      sum_pred_1 <- sum_pred_1 + pred_r1[i]
      if (abs(sum_pred_1) > abs(all_sum_pred_1 / 2)) break
    }
    amplitude_sample[m, 1] <- (pred_r1[i-1] + pred_r1[i]) / 2
  }

  a2 <- find_closest(x_test, range(lat_sample_2)[1])
  b2 <- find_closest(x_test, range(lat_sample_2)[2])
  for (m in 1:mm) {
    if(is.abs) {
      pred_r2 <- abs(pred_f[a2:b2, m])
    } else {
      pred_r2 <- pred_f[a2:b2, m]
    }
    all_sum_pred_2 <- sum(pred_r2)
    sum_pred_2 <- pred_r2[1]
    for (j in 2:length(pred_r2)) {
      sum_pred_2 <- sum_pred_2 + pred_r2[j]
      if (abs(sum_pred_2) > abs(all_sum_pred_2 / 2)) break
    }
    amplitude_sample[m, 2] <- (pred_r2[j-1] + pred_r2[j]) / 2
  }
  return(amplitude_sample)
}
