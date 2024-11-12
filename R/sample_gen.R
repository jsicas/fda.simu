#' @title Gera uma Amostra de Dados Funcionais
#'
#' @export
#'
#' @usage
#' sample_gen(fun_comp, snr, n=10, stand=T, signal=7)
#'
#' @param fun_comp objeto com as \eqn{L} funções componentes, onde cada linha é uma
#' função.
#' @param snr razão sinal-ruído.
#' @param n tamanho da amostra gerada.
#' @param stand Se TRUE (default) as funções componentes são padronizada para terem
#' razão sinal-ruído igual a \code{signal}.
#' @param signal desvio padrão dos dados após padronização.
#'
#' @description
#' Esta função gera uma amostra de dados funcionais agregados baseado na entrada
#' das funções componentes, razão sinal-ruído e tamanho da amostra.
#'
#' @returns Retorna uma \code{matriz} e cada linha representa uma observação.
#'
#' @examples
#' bumps <- f_test()$bumps
#' doppler <- f_test()$doppler
#' par(mfrow=c(2,1))
#' plot(bumps, type='l', main='Bumps'); plot(doppler, type='l', main='Doppler')
#'
#' fun_comp <- matrix(c(bumps, doppler), nrow=2, byrow=T)
#' sample <- sample_gen(fun_comp, snr=5, n=7)
#' plot(1:1024, ylim=c(min(sample), max(sample)), type='n')
#' lapply(1:nrow(sample), \(i) lines(sample[i,], col=i)) |> invisible()

sample_gen <- function(fun_comp, snr, n=10, stand=T, signal=7) {
  if (isTRUE(stand))
    fun_comp <- signal * fun_comp/apply(fun_comp, 1, sd)  # garantindo sd(sinal)=7
  L <- nrow(fun_comp)  # número de curvas componentes
  y <- apply(matrix(runif(L*n), nrow=L), 2,
             function(col) col/sum(col))  # soma dos pesos igual a 1
  fun_agr <- t(y) %*% fun_comp + rnorm(n*ncol(fun_comp), 0, 7/snr)
  return(fun_agr)
}
