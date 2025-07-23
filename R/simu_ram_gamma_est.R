#' @title Função de Simulação para Dados Funcionais com Erro Gamma
#'
#' @importFrom wavethresh wd
#' @importFrom wavethresh GenW
#' @importFrom doFuture %dofuture%
#'
#' @export
#'
#' @details
#' Essa função faz apenas uma única repetição do método. E não será utilizada.
#'
#' @param alpha_comp Matriz com as funções componentes, cada coluna é uma
#'   função.
#' @param I Tamanho da amostra gerada, isto é, quantidade de observações.
#' @param bt Vetor contendo os índices de burn-in e thinning. Se não for
#'   especificado será utilizado uma queima de 15% das observações e um thinning
#'   de 3 observações.
#' @inheritParams ram_gamma
#' @inheritParams gera_ponto

simu_ram_gamma_est <- function(alpha_comp, I, n_ite, alpha, tau, beta, lambda,
                           gamma=2/3, filter.number=5, family='DaubExPhase',
                           lim_sup=15, bt) {
  # geral
  M <- nrow(alpha_comp)  # quantidade de pontos por função
  L <- ncol(alpha_comp)  # número de curvas componentes
  if (missing(bt)) bt <- seq(round(0.15 * n_ite), n_ite, 3)

  # gerando amostra
  y <- apply(matrix(runif(L*I), nrow=L), 2, \(col) col/sum(col))  # pesos
  f <- alpha_comp %*% y
  e <- matrix(rgamma(M*I, shape=beta, rate=lambda), ncol=I)       # erro
  A <- f + e  # amostra

  # DWT
  D <- apply(A, MARGIN=2, \(col) wd(col, filter.number=filter.number,
                                    family=family))

  # gerando e verificando chutes iniciais para theta
  theta_1 <- matrix(0, nrow=M, ncol=I)  # cada coluna é um chute

  for (i in 1:I) {
    theta_1[,i] <- gera_ponto(a=NULL, D[[i]], lim_sup, filter.number, family)
    if (post_gamma(theta_1[,i], D[[i]], beta, lambda, tau, alpha, filter.number,
                   family) == 0) stop('Ponto inválido gerado.')
  }

  # RAM
  delta_D <- foreach(i = 1:I, .combine=cbind,
                     .options.future=list(seed=TRUE)) %dofuture% {
    ram_gamma(theta_1[,i], S_1=NULL, D[[i]], n_ite, alpha, tau, beta,
              lambda, gamma, filter.number=5, family='DaubExPhase') |>
      (\(x) colMeans(x$theta[bt,]))()
  }

  Gamma_hat <- delta_D %*% t(y) %*% solve(y %*% t(y))

  # IDWT
  alpha_hat <- GenW(M, filter.number, family) %*% Gamma_hat  # curvas estimadas

  # calculando erro
  MSE <- colMeans((alpha_hat - alpha_comp)^2)
  AMSE <- mean(MSE)

  # resultado
  return(list('MSE'=MSE, 'AMSE'=AMSE))
}
