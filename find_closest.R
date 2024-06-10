A <- c(1.1, 2.5, 3.8, 4.6)
B <- c(1.0, 2.0, 3.0, 4.0, 5.0)

find_closest <- function(a, B) {
  differences <- abs(B - a)
  min_index <- which.min(differences)
  return(B[min_index])
}

A_replaced <- sapply(A, find_closest, B)
