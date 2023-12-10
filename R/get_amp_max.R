
## max peak over range of latency r (Max Peak)
get_amp_max <- function(pred_f, x_test, lat_sample_1, lat_sample_2,
                        is.abs = FALSE, is.old = FALSE) {
  mm <- length(lat_sample_1)
  amplitude_sample <- matrix(0, mm, 2)

  a1 <- find_closest(x_test, range(lat_sample_1)[1])
  b1 <- find_closest(x_test, range(lat_sample_1)[2])
  for (m in 1:mm) {
    if(is.abs) {
      pred_r1 <- abs(pred_f[a1:b1, m])
      amplitude_sample[m, 1] <- max(pred_r1)
    } else {
      pred_r1 <- pred_f[a1:b1, m]
      amplitude_sample[m, 1] <- pred_r1[which.max(abs(pred_r1))]
    }
  }

  a2 <- find_closest(x_test, range(lat_sample_2)[1])
  b2 <- find_closest(x_test, range(lat_sample_2)[2])
  for (m in 1:mm) {
    if(is.abs) {
      pred_r2 <- abs(pred_f[a2:b2, m])
      amplitude_sample[m, 2] <- max(pred_r2)
    } else {
      pred_r2 <- pred_f[a2:b2, m]
      if(is.old) {
        amplitude_sample[m, 2] <- pred_r2[which.min(abs(pred_r2))]
      } else {
        amplitude_sample[m, 2] <- pred_r2[which.max(abs(pred_r2))]
      }
    }
  }
  return(amplitude_sample)
}

