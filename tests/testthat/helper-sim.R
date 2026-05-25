# Shared simulation setup for all CoSparTA tests.
# Sourced automatically by testthat before any test file runs.

gen_U1_with_covariates <- function(n, X, gamma_true, pi0_true, alpha_true, beta_true, seed) {
  set.seed(seed)
  lambda_i <- as.vector(exp(X %*% gamma_true))
  active   <- rbinom(n, 1, 1 - pi0_true)
  u        <- numeric(n)
  idx      <- which(active == 1)
  if (length(idx) > 0) {
    u[idx] <- rgamma(length(idx), shape = alpha_true, rate = beta_true / lambda_i[idx])
  }
  u / sqrt(sum(u^2))
}

# ---- DGP ----
set.seed(42)
K <- 2
n <- 100
p <- 20
w <- 10

X_cov       <- cbind(1, rbinom(n, 1, 0.5), runif(n, 20, 80))
gamma1_true <- c(0,  0.8,  0.15)
gamma2_true <- c(0, -0.5,  0.25)

U1_true <- cbind(
  gen_U1_with_covariates(n, X_cov, gamma1_true, pi0_true = 0.2,
                         alpha_true = 3,   beta_true = 2,   seed = 1),
  gen_U1_with_covariates(n, X_cov, gamma2_true, pi0_true = 0.2,
                         alpha_true = 1.5, beta_true = 2.5, seed = 2)
)

u2_k1   <- c(rep(1, p/2), rep(0, p/2))
u2_k2   <- c(rep(0, p/2), rep(1, p/2))
U2_true <- cbind(u2_k1 / sqrt(sum(u2_k1^2)),
                 u2_k2 / sqrt(sum(u2_k2^2)))

set.seed(3)
u3_k1   <- c(rgamma(w/2, shape = 2, rate = 1), rep(0, w/2))
u3_k2   <- c(rep(0, w/2), rgamma(w/2, shape = 2, rate = 1))
U3_true <- cbind(u3_k1 / sqrt(sum(u3_k1^2)),
                 u3_k2 / sqrt(sum(u3_k2^2)))

weight       <- c(2, 1.5)
scalenpw     <- sqrt(n) * sqrt(p) * sqrt(w)
sparse       <- 20
lambda_tensor <- (
  weight[1] * (U1_true[,1] %o% U2_true[,1] %o% U3_true[,1]) +
  weight[2] * (U1_true[,2] %o% U2_true[,2] %o% U3_true[,2])
) * scalenpw / sparse

set.seed(99)
X_obs <- array(rpois(n * p * w, lambda_tensor), dim = c(n, p, w))
