########################################
######## Randomly slice sample #########
########################################
sample_rows <- function(data, minimum = 5, maximum = 7) {
  n <- sample(minimum:maximum, 1)
  data |> slice_sample(n = n)
}
