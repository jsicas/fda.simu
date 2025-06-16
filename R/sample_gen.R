#' @title Gera uma Amostra de Dados Funcionais
#'
#' @export
#'
#' @param alpha_comp objeto com as \eqn{L} funções componentes, onde cada coluna é uma
#' função.
#' @param snr razão sinal-ruído.
#' @param I quantidade de observações que se deseja gerar.
#' @param beta,lambda Parâmetros para erro \eqn{Gamma(\beta, \lambda)}.
#' @param erro tipo de erro, podendo ser "normal" ou "gamma".
#'
#' @description
#' Esta função gera uma amostra de dados funcionais agregados baseado na entrada
#' das funções componentes, razão sinal-ruído e tamanho da amostra. É necessário
#' que \eqn{SD(fun\_comp) = 7}.
#'
#' @returns
#' Retorna uma \code{lista} contendo dois elementos: \code{fun}, onde cada
#' coluna é um elemento da amostra gerada; e \code{y}, os quais são os pesos
#' utilizados.
#'
#' @examples
#' # Exemplo 1: Erro normal, L = 2
#' x <- (1:1024)/1024
#' fun_comp <- matrix(c(bumps, doppler), ncol=2)
#' par(mfrow=c(2,1))
#' plot(x, fun_comp[,1], type='l', main='Bumps', ylab='y')
#' plot(x, fun_comp[,2], type='l', main='Doppler', ylab='y')
#' par(mfrow=c(1,1))
#' sample <- sample_gen(fun_comp, snr=5, I=7)
#' plot(x=x, y=x, main='Amostra Gerada (erro Normal)', ylab='', type='n',
#'      ylim=c(min(sample$fun), max(sample$fun)))
#' for (i in 1:ncol(sample$fun)) lines(x, sample$fun[,i], col=i)
#'
#' # Example 2: Erro gamma, L = 3
#' fun_comp <- cbind(fun_comp, f_test()$blocks)
#' plot(x, fun_comp[,3], type='l', main='blocks', ylab='y')
#' sample <- sample_gen(fun_comp, I=7, beta=2, lambda=2, erro='gamma')
#' plot(x=x, y=x, main='Amostra Gerada (erro Gamma)', ylab='', type='n',
#'      ylim=c(min(sample$fun), max(sample$fun)))
#' for (i in 1:ncol(sample$fun)) lines(x, sample$fun[,i], col=i)
#'
#' plot(x, fun_comp %*% sample$y[,1], type='l', main='Obs 1', ylab='y')
#' lines(x, sample$fun[,1], col='blue')
#' legend('topright', bty='n', lwd=2,
#'        legend=c('Real', 'Amostra'), col=c('black', 'blue'))

sample_gen <- function(alpha_comp, snr, I=10, beta, lambda,
                       erro=c('normal', 'gamma')) {
  erro <- match.arg(erro)
  # verificações
  if (all(apply(alpha_comp, 2, sd) == 7)) {
    stop('É necessário que sd(sinal) = 7 para todas as funções componentes.')
  }
  if ((missing(snr) & erro == 'normal') |
      ((missing(beta) | missing(lambda)) & erro == 'gamma'))
    stop('Parâmetros da função incorretos.')

  # função
  L <- ncol(alpha_comp)  # número de curvas componentes
  M <- nrow(alpha_comp)  # número de pontos por obs.
  y <- matrix(runif(L*I), nrow=L) |> # soma dos pesos igual a 1
    apply(2, \(col) col/sum(col))
  if (erro == 'normal') {
    signal <- 7
    fun_agr <- alpha_comp %*% y + rnorm(I*M, 0, signal/snr)
  } else if(erro == 'gamma') {
    fun_agr <- alpha_comp %*% y + rgamma(I*M, shape=beta, rate=lambda)
  }
  return(list('fun'=fun_agr, 'y'=y))
}
