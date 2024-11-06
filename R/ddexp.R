#' @title Função Densidade de Probabilidade da Exponencial Dupla (Laplace)
#'
#' @usage ddexp(x, mu, lambda)
#'
#' @param x valor a ser avaliado.
#' @param mu parâmetro de locação.
#' @param lambda parâmetro de escala.
#'
#' @details
#' Seja \eqn{X} uma variável aleatória com distribuição de Laplace, então, sua função
#' densidade de probabilidade é dada por:
#' \deqn{f(x) = \dfrac{\exp\left\{ -\frac{|x-\mu|}{\beta} \right\}}{2\beta}}

ddexp <- function(x, mu, lambda) {
  return(exp(-abs(x-mu)/lambda)/(2*lambda))
}
