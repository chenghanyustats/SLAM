regfcn_sin <- function(x, k) {
  reg <- -2 * sin(2 * pi * x + k/15 - 0.3)
  return(reg)
}
