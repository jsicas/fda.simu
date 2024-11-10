#' @title Regra de Encolhimento Bayesiana com Priori Epanechnikov Generalizada
#'
#' @export
#'
#' @usage epanec_shrink(d, a, b, lambda)
#'
#' @param d coeficiente de ondaleta.
#' @param a parâmetro \eqn{\alpha} da mistura.
#' @param b parâmetro \eqn{\beta} da Epanechnikov Generalizada. Ver
#' \code{\link{ddexp}}.
#' @param lambda parâmetro da priori exponencial atribuida à \eqn{\sigma^2}.
#'
#' @returns Retorna um vetor de mesmo tamanho inicial com os coeficientes de
#' ondaleta estimados.
#'
#' @details
#' Uma variável aleatória \eqn{X} tem distribuição de Epanechnikov
#' Generalizada com parâmetro \eqn{\beta} (\eqn{DEG(\beta)}) se, e somente se, sua
#' densidade é dada por:
#' \deqn{g(x; \beta) = \frac{3}{4\beta^3} (\beta^2 - x^2)
#' \mathcal{I}^{(x)}_{(-\beta, \beta)}}
#' além disso, \eqn{\mathbb{E}(X) = 0} e \eqn{\mathbb{V}(X) = \frac{\beta^2}{5}}.
#'
#' Além disso, a priori utilizada é
#' \deqn{\pi(d; \alpha,\beta) = \alpha\delta_0(d) + (1-\alpha) g(d;\beta)}
#' onde \eqn{\alpha \in (0,1)}, \eqn{\delta_0} é o delta de Dirac centrado em
#' \eqn{0} e \eqn{g(d; \beta)} é a densidade de Epanechnikov Generalizada
#' com parâmetro \eqn{\beta}.
#'
#' Com isso, e assumindo a priori de \eqn{\sigma^2} sendo uma distribuição
#' exponencial com parâmetro \eqn{\lambda}, obtem-se a regra de shrinkage Bayesiana:
#' \deqn{\delta(d) = \dfrac{(1-\alpha) \frac{3 \sqrt{2\lambda}}{8 \beta^3}
#' \left[ \frac{2 \lambda\beta^2 + 3\beta\sqrt{2\lambda} + 3}{2\lambda^2}
#' \left( e^{-\sqrt{2\lambda}(\beta - d)} - e^{-\sqrt{2\lambda}(\beta+d)} \right)
#' +\frac{(\lambda\beta^2 - 3) d\sqrt{2\lambda} - d^3\lambda\sqrt{2\lambda}}{\lambda^2}
#' \right]}{\alpha \mathcal{ED}\left(0, \frac{1}{\sqrt{2\lambda}} \right) +
#' (1 - \alpha) \frac{3\sqrt{2\lambda}}{8\beta^3} \left[\frac{\beta}{\lambda}
#' \left( e^{-\sqrt{2\lambda}(\beta + d)} + e^{-\sqrt{2\lambda}(\beta-d)} \right)
#' + \frac{2}{\sqrt{2\lambda}} (\beta^2 - d^2 - \frac{1}{\lambda}) \right]}
#' }
#' onde \eqn{\mathcal{ED}\left(0, \frac{1}{\sqrt{2\lambda}} \right)} é a função
#' de densidade da distribuição exponencial dupla com média \eqn{0} e parâmetro
#' de escala \eqn{\frac{1}{\sqrt{2\lambda}}}.
#'
#' @references
#' NULL
#'
#' @examples
#' curve(epanec_shrink(x, 0.8, 4, 6), -6, 6)

epanec_shrink <- function(d, a, b, lambda) {
  exp1 <- (2*lambda*b^2 + 3*b*sqrt(2*lambda) + 3)/(2*lambda^2) *
    (exp(-(b-d)*sqrt(2*lambda)) - exp(-(b+d)*sqrt(2*lambda)))
  exp2 <- b/lambda*(exp(-(b+d)*sqrt(2*lambda)) + exp(-(b-d)*sqrt(2*lambda)))
  num <- (1-a)*3*sqrt(2*lambda)/(8*b^3) * (exp1 +
    ((lambda*b^2 - 3)*d*sqrt(2*lambda) - lambda*d^3*sqrt(2*lambda))/lambda^2)
  den <- a*ddexp(d, 0, 1/sqrt(2*lambda)) + ((1-a)*3*sqrt(2*lambda))/(8*b^3) *
    (exp2 + 2/sqrt(2*lambda)*(b^2 - d^2 - 1/lambda))
  return(num/den)
}
