find_closest <- function(vec, number) {
  which(abs(vec - number) == min(abs(vec - number)))
}
