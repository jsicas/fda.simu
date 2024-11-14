#' @title Função Densidade de Probabilidade da Exponencial Dupla (Laplace)
#'
#' @usage ddexp(x, mu, scale)
#'
#' @param x valor a ser avaliado.
#' @param mu parâmetro de locação.
#' @param scale parâmetro de escala.
#'
#' @details
#' Seja \eqn{X} uma variável aleatória com distribuição de Laplace, então, sua
#' função densidade de probabilidade é dada por:
#' \deqn{f(x) = \dfrac{1}{2\beta} \exp\left\{ -\dfrac{|x-\mu|}{\beta} \right\}}

ddexp <- function(x, mu, scale) {
  if (scale <= 0) stop('Parâmetro scale deve ser positivo.')
  return(exp(-abs(x-mu)/scale)/(2*scale))
}
