set.seed(1234)

# define set of data with 3 covariates for testing
n <- 20
X3 <- matrix(rnorm(n * 3), nrow = n, ncol = 3)
X3 <- cbind(matrix(1, nrow = n, ncol = 1), X3) # add intercept
theta3_truth <- matrix(c(1, 0, 1, 0), ncol = 1)
theta3_truth_bool <- as.logical(theta3_truth)
theta3_truth_idx <- which(theta3_truth_bool)
y3 <- X3 %*% theta3_truth + rnorm(n)
