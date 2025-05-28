#' @title Função de Simulação para Dados Funcionais com Erro Gamma
#'
#' @importFrom wavethresh wd
#' @importFrom wavethresh GenW
#' @importFrom doFuture %dofuture%
#'
#' @export
#'
#' @param alpha_comp Matriz com as funções componentes, cada coluna é uma
#'   função.
#' @param I Tamanho da amostra gerada, isto é, quantidade de observações.
#' @inheritParams ram_gamma
#' @param bt Vetor contendo os índices de burn-in e thinning. Se não for
#'   especificado será utilizado uma queima de 5% das observações e um thinning
#'   de 3 observações.

simu_ram_gamma <- function(alpha_comp, I, n_ite, beta, lambda, alpha, tau,
                           gamma=2/3, filter.number=5, family='DaubExPhase',
                           lim_sup=15, bt=NULL) {
  # geral
  M <- nrow(alpha_comp)  # quantidade de pontos por função
  L <- ncol(alpha_comp)  # número de curvas componentes
  if (is.null(bt)) bt <- seq(round(0.05 * n_ite), n_ite, 3)

  # gerando amostra
  y <- apply(matrix(runif(L*I), nrow=L), 2, \(col) col/sum(col))  # pesos
  f <- alpha_comp %*% y
  e <- matrix(rgamma(M*I, shape=beta, rate=lambda), ncol=I)     # erro
  A <- f + e  # amostra

  # DWT
  D <- apply(A, MARGIN=2, \(col) wd(col, filter.number, family))

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
  alpha_hat <- GenW(M, filter.number, family) %*% delta_D  # curvas estimadas

  # calculando erro

}
