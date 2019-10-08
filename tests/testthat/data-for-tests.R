set.seed(1234)

# define set of data with 3 covariates for testing
n <- 20
X3 <- matrix(rnorm(n * 3), nrow = n, ncol = 3)
X3 <- cbind(matrix(1, nrow = n, ncol = 1), X3) # add intercept
theta3 <- matrix(c(1, 0, 1, 0), ncol = 1)
y3 <- X3 %*% theta3 + rnorm(n)
