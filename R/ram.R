#' @title Gera Amostra da Posteriori com Erro Gamma
#'
#' @importFrom wavethresh wd
#' @importFrom mvtnorm rmvnorm
#'
#' @export
ram <- function(theta_1 = NULL, S_1 = NULL, y, alpha=0.8, tau=2,
                beta, lambda, n_ite = 50000, gamma = 2/3,
                filter.number=5, family='DaubExPhase') {
  # criando objetos
  n <- length(y)                  # quantidade de pontos por função
  theta <- matrix(0, n_ite, n)    # matriz contendo amostra de theta
  S_l <- vector(mode='list', length=n_ite)      # lista para armazenar S_l
  gamma_l <- vector(mode='numeric', n_ite - 1)  # vetor para armazenar gamma_l

  # gerando amostra
  dwt <- wd(y, filter.number=filter.number, family=family)  # coeficientes empíricos

  # definindo alguns parâmetors
  if (is.null(S_1)) S_1 <- diag(n)  # chute inicial para matriz S caso não definida
  if (is.null(theta_1)) theta_1 <- gera_ponto(n, y)
  theta[1,] <- theta_1             # chute inicial
  S_l[[1]] <- S_1                  # chute inicial
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

  return(list('theta'=theta, 'S'=S_l, 'gamma_l'=gamma_l,
              'parametros'=c('n_ite'=n_ite, 'alpha'=alpha, 'tau'=tau, 'beta'=beta,
                             'lambda'=lambda, 'gamma'=gamma,
                             'filter.number'=filter.number, 'family'=family)
  ))
}
