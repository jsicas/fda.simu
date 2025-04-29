# Neste arquivo são construidos exemplos para a utilização da regra de shrinkage
# com erro gamma.

# packages
require(fda.simu)
require(wavethresh)
require(mvtnorm)

# Exemplo 1 ---------------------------------------------------------------
# Considere o modelo
#         y = f + e
# onde
#     y: vetor de valores observados;
#     f: é a função verdadeira;
#     e: é o erro Gamma(beta, lambda).
#
# Daí, no domínio das ondaletas, temos:
#         d = theta + epsilon
# onde
#     d: vetor de coenficientes empíricos de ondaleta;
#     theta: vetor de coeficientes verdadeiros de ondaleta;
#     epsilon: é o erro com distribuição calculada através do método do
#       jacobiano.

set.seed(282829)

# definindo parâmetros
n <- 8
tau <- 5
gamma <- 2/3
n_ite <- 1000
beta <- 10; lambda <- 5  # parâmetros do erro Gamma
f <- f_test(n)$doppler
e <- rgamma(n, shape=beta, rate=lambda)
y <- f + e  # valores observados

# aplicando DWT
dwt <- wd(y, filter.number=5, family='DaubExPhase')  # coeficientes empíricos
theta_true <- wd(f, filter.number=5, family='DaubExPhase') |>
  (\(x) c(accessC(x, lev=0), x$D))()  # coeficientes verdadeiros

# criando objetos para simulação
theta <- matrix(0, n_ite, n)    # matriz contendo amostra de theta
S_l <- vector(mode='list', length=n_ite)      # lista para armazenar S_l
gamma_l <- vector(mode='numeric', n_ite - 1)  # vetor para armazenar gamma_l

# gerandop ponto inicial
theta_1 <- gera_ponto(n, y, filter.number=5, family='DaubExPhase')
post_gamma(theta_1, dwt, beta, lambda, tau, alpha)
theta[1,] <- theta_1             # chute inicial para theta
S_l[[1]] <- diag(n)              # chute inicial para S
eta <- seq(1, n_ite)^(-gamma)    # parâmetro do algorítimo RAM

for (i in 2:n_ite) {if (i %% 10000 == 0) message(i)
  # proposta
  U_l <- t(rmvnorm(1, rep(0, n), diag(n)))
  theta_star <- t(theta[i-1,] + S_l[[i-1]] %*% U_l)

  # taxa de aceitação
  gamma_l[i-1] <- min(1, post_gamma(theta_star, dwt, beta, tau, lambda, alpha)/
                        post_gamma(theta[i-1,], dwt, beta, tau, lambda, alpha))

  if (rbinom(1, 1, gamma_l[i-1]) == 1) theta[i,] <- theta_star
  else theta[i,] <- theta[i-1,]

  # atualizando S_l
  A <- S_l[[i-1]] %*% (diag(n) + eta[i]*(gamma_l[i-1] - gamma) *
                         U_l %*% t(U_l) / as.vector(t(U_l) %*% U_l)) %*% t(S_l[[i-1]])
  S_l[[i]] <- t(chol(A))
}

theta_true
colMeans(theta)


set.seed(282829)
aa <- ram_gamma(theta_1, S_1 = NULL, f, alpha=0.8, tau=2, beta, lambda,
                n_ite = 1000, gamma = 2/3,
                filter.number=5, family='DaubExPhase')
theta_true
colMeans(aa$theta)




set.seed(282829)
n <- 16
f <- f_test(n)$bumps
beta <- 10; lambda <- 80  # parâmetros da gamma
y <- f + rgamma(n, shape=beta, rate=lambda)
d <- wd(y, filter.number=5, family='DaubExPhase')
d_emp <- c(accessC(d, lev=0), d$D)

# gerando valores positivos
x <- c(4, 5, 3, 4, 2, 4, 4, 2, 4, 5, 3, 4, 2, 4, 4, 2)
W <- t(GenW(n, filter.number=5, family='DaubExPhase'))
theta_1 <- d_emp - W %*% x

t(W) %*% as.vector((d_emp - theta_1))  # wdt_i


# parametros da distribuição
alpha <- 0.8; tau <- 2; beta <- 10; lambda <- 40
post_gamma(theta_1, d, alpha, beta, tau, lambda)
