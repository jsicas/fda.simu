#' @title Gera Amostra da Posteriori com Erro Gamma por Meio do Algorítimo RAM
#'
#' @importFrom wavethresh wd
#' @importFrom mvtnorm rmvnorm
#'
#' @export
#'
#' @param theta_1 Chute inicial para \eqn{\theta}. Se não for definido, será
#'   gerado um ponto aleatório baseado na função `gera_ponto`.
#' @param S_1 Chute inicial para a matriz S. Se não for definida será utilziado
#'   a identidaide.
#' @param f Valores verdadeiros.
#' @param alpha Parâmetro da mistura da priori spike and slab utilizada.
#' @param tau Parâmetro da logística presente na priori.
#' @param beta,lambda parâmetros do erro \eqn{Gamma(\beta, \lambda)}.
#' @param n_ite número de iterações.
#' @param gamma Parâmetro do algorítimo RAM.
#' @inheritParams wavethresh::wd
#'
#' @references
#' Vihola, M. (2012). Robust adaptive Metropolis algorithm with coerced
#' acceptance rate. Stat Comput 22, 997–1008. DOI:
#' [10.1007/s11222-011-9269-5](https://doi.org/10.1007/s11222-011-9269-5).
#'
#' Sousa, A.R.S., Garcia, N.L. (2023). Wavelet shrinkage in
#' nonparametric regression models with positive noise. \emph{Journal of
#' Statistical Computation and Simulation}. DOI:
#' [10.1080/00949655.2023.2215372](https://doi.org/10.1080/00949655.2023.2215372).
#'
#' @examples
#' sample <- ram_gamma(f = f_test(n = 8)$bumps, beta = 10, lambda = 30, n_ite = 3000)
#' colMeans(sample$theta)

ram_gamma <- function(theta_1 = NULL, S_1 = NULL, f, alpha=0.8, tau=2,
                      beta, lambda, n_ite = 50000, gamma = 2/3,
                      filter.number=5, family='DaubExPhase') {

  # criando objetos
  n <- length(f)                  # quantidade de pontos por função
  theta <- matrix(0, n_ite, n)    # matriz contendo amostra de theta
  S_l <- vector(mode='list', length=n_ite)      # lista para armazenar S_l
  gamma_l <- vector(mode='numeric', n_ite - 1)  # vetor para armazenar gamma_l

  # gerando amostra
  y <- f + rgamma(n, shape=beta, rate=lambda)  # gerando valores observados
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
    gamma_l[i-1] <- min(1, post_gamma(theta_star, dwt, alpha, beta, tau, lambda)/
                          post_gamma(theta[i-1,], dwt, alpha, beta, tau, lambda))

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
