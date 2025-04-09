# Tentativa de implementar algorítimo RAM para caso unidimensional, considerando
# theta ~ U(0,1), y_i | theta ~ Ber(theta). Com isso, pode-se obter a posteriori
# dada por theta | y ~ Beta(sum(y) + 1, n - sum(y) + 1)

theta_true <- 0.6; n <- 100; l <- 300
y <- rbinom(n, 1, theta_true)
theta <- vector(mode='numeric', length=l); theta[1] <- 0.1
S <- vector(mode='numeric', length=l); S[1] <- 1
eta <- rbeta(l, 1, 1)
gamma <- 0.5

for (i in 2:l) {
  U_l <- rnorm(1)
  theta_star <- theta[i-1] + S[i-1] %*% U_l
  gamma_l <- min(1, dbeta(theta_star, sum(y) + 1, n - sum(y) + 1)/
                  dbeta(theta[i-1], sum(y) + 1, n - sum(y) + 1))
  A <- S[i-1] %*% (diag(1) + eta[i]*(gamma_l - gamma) %*%
                     U_l %*% t(U_l)/sum(U_l^2)) %*% t(S[i-1])
  S[i] <- t(chol(A))
  # S[i] <- sqrt(A)
  if (rbinom(1,1, gamma_l) == 1) theta[i] <- theta_star else theta[i] <- theta[i-1]
}

mean(theta)
