# Testes com o algorítimo RAM.

# packages
library(mvtnorm)
library(wavethresh)
library(fda.simu)

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
gamma <- 2/3
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

mean(theta[1000:l])



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

colMeans(theta[-(1:500),])



# Amostrando da Posteriori com Erro Gamma ---------------------------------

# definindo parâmetros da função
set.seed(282829)
n <- 8                # quantidade de pontos
L <- 500000                 # quantidade de iterações
theta_1 <- c(-2.2445997, -4.2212565, -12.4897425, 5.8349722,
             -1.5304893, 3.9793896, 0.7659587, -8.0117955)  # chute inicial
S_1 <- diag(n)         # chute inicial
alpha <- 0.8
tau <- 2
beta <- 10  # parâmetro gamma
lambda <- 40  # parâmetro gamma
# d <- wd(f_test(n)$bumps,
#         # + dgamma(n, shape=beta, rate=lambda),
#         family='DaubExPhase', filter.number=5)  # coeficientes empíricos

theta <- matrix(0, L, n); theta[1,] <- theta_1
S_l <- vector(mode='list', length=L); S_l[[1]] <- S_1
gamma_l <- vector(mode='numeric', L - 1)  # taxa de aprendizado
gamma <- 2/3
eta <- seq(1, L)^(-gamma)

for (i in 2:L) {if (i %% 5000 == 0) message(i)
  # proposta
  U_l <- t(rmvnorm(1, rep(0, n), diag(n)))
  theta_star <- t(theta[i-1,] + S_l[[i-1]] %*% U_l)

  # taxa de aceitação
  gamma_l[i-1] <- min(1, post_gamma(theta_star, d, alpha, beta, tau, lambda)/
                        post_gamma(theta[i-1,], d, alpha, beta, tau, lambda))

  if (rbinom(1, 1, gamma_l[i-1]) == 1) theta[i,] <- theta_star
  else theta[i,] <- theta[i-1,]

  # atualizando S_l
  A <- S_l[[i-1]] %*% (diag(n) + eta[i]*(gamma_l[i-1] - gamma) *
                         U_l %*% t(U_l) / as.vector(t(U_l) %*% U_l)) %*% t(S_l[[i-1]])
  S_l[[i]] <- t(chol(A))
}




post_gamma(theta_1, d, alpha, beta, tau, lambda)



# gerando w(d - theta) no suporte
set.seed(282829)
n <- 8
f <- f_test(n)$bumps
beta <- 10; lambda <- 40  # parâmetros da gamma
y <- f + rgamma(n, shape=beta, rate=lambda)
d <- wd(y, filter.number=5, family='DaubExPhase')
d_emp <- c(accessC(d, lev=0), d$D)

# gerando valores positivos
x <- c(4, 5, 3, 4, 2, 4, 4, 2)
W <- t(GenW(n, filter.number=5, family='DaubExPhase'))
theta_1 <- d_emp - W %*% x

t(W) %*% as.vector((d_emp - theta_1))  # wdt_i


# parametros da distribuição
alpha <- 0.8; tau <- 2; beta <- 10; lambda <- 40
post_gamma(theta_1, d, alpha, beta, tau, lambda)







