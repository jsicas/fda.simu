#' @title Função de Simulação para Dados Funcionais com Erro Gamma
#'
#' @importFrom wavethresh wd GenW
#' @importFrom foreach foreach
#' @importFrom doFuture %dofuture%
#' @import future
#'
#' @export
#'
#' @details
#' Essa função faz várias repetições.
#'
#' @param alpha_comp Matriz com as funções componentes, cada coluna é uma
#'   função.
#' @param rep número de replicações.
#' @param I Tamanho da amostra gerada, isto é, quantidade de observações.
#' @param bt Vetor contendo os índices de burn-in e thinning. Se não for
#'   especificado será utilizado uma queima de 15% das observações e um thinning
#'   de 3 observações.
#' @inheritParams ram_gamma
#' @inheritParams gera_ponto
#'
#' @example examples/exam_simu_ram_gamma.R

simu_ram_gamma <- function(alpha_comp, I, n_ite, rep, alpha, tau, beta, lambda,
                           gamma=2/3, bt,
                           filter.number=5, family='DaubExPhase') {
  # geral
  M <- nrow(alpha_comp)  # quantidade de pontos por função
  L <- ncol(alpha_comp)  # número de curvas componentes
  if (missing(bt)) bt <- seq(round(0.15 * n_ite), n_ite, 3)

  result <- foreach(j = 1:rep, .combine=rbind,
          .options.future=list(seed=TRUE)) %dofuture% {

    # gerando amostra
    y <- apply(matrix(runif(L*I), nrow=L), 2, \(col) col/sum(col))  # pesos
    f <- alpha_comp %*% y
    e <- matrix(rgamma(M*I, shape=beta, rate=lambda), ncol=I)       # erro
    A <- f + e  # amostra

    # DWT
    D <- apply(A, MARGIN=2, \(col) wd(col, filter.number, family))

    # gerando e verificando chutes iniciais para theta
    theta_1 <- matrix(0, nrow=M, ncol=I)  # armazena o chute em colunas

    for (i in 1:I) {  # verificando chutes
      theta_1[,i] <- gera_ponto(a=NULL, d=D[[i]], lim_sup=beta/lambda,
                                filter.number=filter.number, family=family)
      if (post_gamma(theta_1[,i], d=D[[i]], beta=beta, lambda=lambda, tau=tau,
                     alpha=alpha, filter.number=filter.number,
                     family=family) == 0) stop('Ponto inválido gerado.')
    }

    # RAM
    delta_D <- matrix(0, M, I)  # pre-alocando

    for (i in 1:I) {
      # d_i representa a i-esima colunca de delta_D
      d_i <- ram_gamma(theta_1[,i], S_1=NULL, d=D[[i]], n_ite=n_ite, alpha=alpha,
                       tau=tau, beta=beta, lambda=lambda, gamma=gamma,
                       filter.number=filter.number, family=family)
      delta_D[,i] <- colMeans(d_i$theta[bt,])
    }

    Gamma_hat <- delta_D %*% t(y) %*% solve(y %*% t(y))

    # IDWT
    alpha_hat <- GenW(M, filter.number, family) %*% Gamma_hat  # curvas estimadas

    # calculando erro
    MSE <- colMeans((alpha_hat - alpha_comp)^2)
    AMSE <- mean(MSE)

    # resultado
    return(c('MSE'=MSE, 'AMSE'=AMSE))
  }

  return(result)
}
