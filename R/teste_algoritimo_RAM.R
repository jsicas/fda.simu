# Testes com o algorítimo RAM.

# packages
require(mvtnorm)
require(wavethresh)
require(fda.simu)

# Caso Unidimensional -----------------------------------------------------
# Considerando theta ~ U(0,1), y_i | theta ~ Ber(theta). Com isso, pode-se obter
# a posteriori dada por
#             theta | y ~ Beta(sum(y) + 1, n - sum(y) + 1)

theta_true <- 0.6  # parâmetro verdadeiro
n <- 100           # número de obs. geradas
l <- 5000          # número de iterações
y <- rbinom(n, 1, theta_true)
theta <- vector(mode='numeric', length=l); theta[1] <- 0.2  # chuite inicial
S <- vector(mode='numeric', length=l); S[1] <- 1
gamma <- 0.7
eta <- seq(1, l)^(-gamma)

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

# mean(theta[1000:l])



# Algorítimo RAM para Normal Bivariada ------------------------------------
# Gerando de uma Normal variada com média (-1, 1)' e matriz de variâncias
# [5, 2]
# [2, 5]

l <- 4500  # quantidade de iterações
gamma <- 2/3

mu <- c(-1,1)
Sigma <- matrix(c(5, 2,
                  2, 5), 2)
theta <- matrix(0, l, 2)
theta[1,] <- c(0,0)  # chute inicial
S <- vector(mode='list', length=l)
S[[1]] <- matrix(c(2, 0,
                   3, 3), 2, byrow=T)

gamma_l <- vector(mode='numeric', l-1)


for (i in 2:l) {
  U_l <- t(rmvnorm(1, c(0,0), diag(2)))
  theta_star <- t(theta[i-1,] + S[[i-1]] %*% U_l)

  gamma_l[i-1] <- min(1, dmvnorm(theta_star, mu, Sigma)/
                        dmvnorm(theta[i-1,], mu, Sigma))

  if (rbinom(1, 1, gamma_l[i-1]) == 1) {
    theta[i,] <- theta_star
    } else theta[i,] <- theta[i-1,]

  A <- S[[i-1]] %*% (diag(2) + eta[i]*(gamma_l[i-1] - gamma) *
                       U_l %*% t(U_l)/as.vector(t(U_l) %*% U_l)) %*% t(S[[i-1]])
  S[[i]] <- t(chol(A))
}

# colMeans(theta[-(1:500),])



# Amostrando da Posteriori com Erro Gamma ---------------------------------
#' @title teste
#'
#' @export
#'
#' @importFrom mvtnorm rmvnorm
#' @importFrom wavethresh GenW
#'
#' @param theta_1 chute inicial.
#' @param L número de iterações.
#' @param S_1 matriz triangular inferior com diagonal positiva.
#' @param d coeficientes empíricos de ondaletas.
#' @param alpha hiperparâmetro.
#' @param tau hiperparâmetro.
#' @param lambda parâmetro da Gamma.
#' @param gamma

# ram <- function(theta_1, L, S_1, d, alpha, tau, lambda, gamma=2/3) {

# definindo parâmetros da função
n <- 32               # quantidade de pontos
theta_1 <- rep(20, n)  # chute inicial
S_1 <- diag(n)
alpha <- 0.8
tau <- 2
lambda <- 4
d <- wd(f_test(n)$bumps + dgamma(n, shape=4, rate=lambda),
        family='DaubExPhase', filter.number=5)  # coeficientes empíricos
L <- 5

dim <- length(theta_1)  # dimensão de theta
theta <- matrix(0, L, dim)
theta[1,] <- theta_1
S_l <- vector(mode='list', length=L)
S_l[[1]] <- S_1  # chute inicial
gamma_l <- vector(mode='numeric', L - 1)  # taxa de aprendizado
eta <- seq(1, L)^(-gamma)

for (i in 2:L) {
  # proposta
  U_l <- t(rmvnorm(1, rep(0, dim), diag(dim)))
  theta_star <- t(theta[i-1,] + S_l[[i-1]] %*% U_l)

  # taxa de aceitação
  gamma_l[i-1] <- min(1, post_gamma(theta_star, d, alpha, tau, lambda)/
                        post_gamma(theta[i-1,], d, alpha, tau, lambda))

  if (rbinom(1, 1, gamma_l[i-1]) == 1) {
    theta[i,] <- theta_star
  } else {theta[i,] <- theta[i-1,]}

  # atualizando S_l
  A <- S_l[[i-1]] %*% (diag(dim) + eta[i]*(gamma_l[i] - gamma) *
                         U_l %*% t(U_l) / as.vector(t(U_l) %*% U_l)) %*% t(S_l[[i-1]])
  S_l[[i]] <- t(chol(A))
}
# return(list('theta'=theta, 'gamma_l'=gamma_l))
# }




# ram(theta_1, L = 50, S_1, d, alpha, tau, lambda, gamma)


post_gamma(theta_star, d, alpha, tau, lambda)

