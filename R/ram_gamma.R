#' @title Gera Amostra da Posteriori com Erro Gamma por Meio do Algorítimo RAM
#'
#' @importFrom wavethresh wd
#' @importFrom wavethresh GenW
#' @importFrom mvtnorm rmvnorm
#'
#' @export
#'
#' @param theta_1 Chute inicial para \eqn{\theta}.
#' @param S_1 Chute inicial para a matriz S. Se não for definida será utilziado
#'   a matriz identidaide.
#' @param d coeficientes empíricos de ondaleta. Objeto da classe `wd`, ver
#'   função [wavethresh::wd].
#' @param alpha Parâmetro da mistura da priori spike and slab utilizada.
#' @param tau Parâmetro da logística presente na priori.
#' @param beta,lambda parâmetros do erro \eqn{Gamma(\beta, \lambda)}.
#' @param n_ite número de iterações (tamanho da cadeia gerada).
#' @param gamma Parâmetro do algorítimo RAM. Por padrão, utiliza-se
#'   \eqn{\frac{2}{3}}.
#' @inheritParams wavethresh::wd
#'
#' @references Vihola, M. (2012). Robust adaptive Metropolis algorithm with
#' coerced acceptance rate. Stat Comput 22, 997–1008. DOI:
#' [10.1007/s11222-011-9269-5](https://doi.org/10.1007/s11222-011-9269-5).
#'
#' Sousa, A.R.S., Garcia, N.L. (2023). Wavelet shrinkage in
#' nonparametric regression models with positive noise. \emph{Journal of
#' Statistical Computation and Simulation}. DOI:
#' [10.1080/00949655.2023.2215372](https://doi.org/10.1080/00949655.2023.2215372).
#'
#' @example examples/ram_gamma_exam.R

ram_gamma <- function(theta_1, S_1=NULL, d, n_ite=50000, alpha=0.8,
                      tau=2, beta, lambda, gamma=2/3,
                      filter.number=5, family='DaubExPhase') {
  # verificando ponto inicial detheta
  if (post_gamma(theta_1, d, beta, lambda, tau, alpha, filter.number, family) == 0)
    stop('Ponto inicial inválido, forneça um ponto com densidade maior que 0.')

  # criando objetos
  M <- 2^nlevelsWT(d)            # quantidade de pontos por função
  theta <- matrix(0, n_ite, M)   # matriz contendo amostra de theta
  gamma_l <- vector(mode='numeric', n_ite - 1)  # vetor para armazenar gamma_l
  eta <- seq(1, n_ite)^(-gamma)  # parâmetro do algorítimo RAM
  W <- t(GenW(M, filter.number=filter.number, family=family))
  theta_mudou <- TRUE            # indica quando theta muda

  # armazenando chute inicial
  if (is.null(S_1)) S_1 <- diag(M)  # chute inicial se S não definida
  theta[1,] <- theta_1              # chute inicial para theta
  S_l <- S_1                        # chute inicial para S

  for (i in 2:n_ite) {
    # proposta
    U_l <- matrix(rnorm(M))
    theta_star <- t(theta[i-1,] + S_l %*% U_l)

    if (theta_mudou == TRUE) {
      den <- post_gamma(theta[i-1,], d, beta, lambda, tau, alpha,
                        filter.number, family, W)
      theta_mudou <- FALSE
    }

    # taxa de aceitação
    gamma_l[i-1] <- min(1, post_gamma(theta_star, d, beta, lambda, tau, alpha,
                                        filter.number, family, W)/den)

    if (rbinom(1, 1, gamma_l[i-1]) == 1) {
      theta[i,] <- theta_star
      theta_mudou <- TRUE
    } else {
      theta[i,] <- theta[i-1,]
    }

    # atualizando S_l
    A <- S_l %*% (diag(M) + eta[i]*(gamma_l[i-1] - gamma) * U_l %*% t(U_l) /
                    as.vector(t(U_l) %*% U_l)) %*% t(S_l)
    S_l <- t(chol(A))
  }

  return(list('theta'=theta, 'gamma_l'=gamma_l,
              'parametros'=c('n_ite'=n_ite, 'alpha'=alpha, 'tau'=tau,
                             'beta'=beta, 'lambda'=lambda, 'gamma'=gamma,
                             'filter.number'=filter.number, 'family'=family)
  ))
}
