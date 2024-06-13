################################## Missing Data Generator ####################################
########################### Suppose 50% missing - Missing at Random ##########################
########################### Suppose a reasonable weights #####################################
set.seed(123)
miss_pattern <- matrix(c(rep(1, 14), 0, rep(1, 4)), nrow = 1, byrow = TRUE)
ampute_test <- ampute(data = famx, prop = 0.5, patterns = miss_pattern)
md.pattern(ampute_test$amp, rotate.names = TRUE) # Missing pattern
ampute_test$mech
ampute_test$weights
reasonable_weights <- ampute_test$weights
reasonable_weights[1,] <- c(0, 0, 0, 0, 0, 0.5, 0, 0, 0, 0, 0, 1, 1, 2, 0, 0, 0, 0, 0)
ampute_test <- ampute(data = famx, prop = 0.5, patterns = miss_pattern, 
                      mech = "MAR", weights = reasonable_weights) # Missing proportions
miss50_famx <- ampute_test$amp

###############################################################################################
variable_to_ampute <- "newx"
n <- nrow(famx)
num_missing <- floor(0.5*n)
missing_indices <- sample(1:n, num_missing)
miss50_famx <- famx
miss50_famx[missing_indices, variable_to_ampute] <- NA

