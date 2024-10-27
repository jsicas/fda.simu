#' @title Regra de Encolhimento Bayesiana com Priori Logística
#'
#' @export
#'
#' @importFrom purrr map
#'
#' @usage logis_shrink_map(d, a, s, t)
#'
#' @param d vetor de coeficientes de ondaleta empírico.
#' @param a parâmetro \eqn{\alpha} da mistura.
#' @param s desvio padrão dos coeficientes de ondaleta.
#' @param t parâmetro \eqn{\tau} da logística.

logis_shrink_map <- function(d, a, s, t) {
  u <- rnorm(10000)
  delta <- map(d, ~(1-a) * mean((s*u + .) * dlogis(s*u + ., scale=t)) /
                  (a * dnorm(., sd=s)/s + (1-a) * mean(dlogis(s*u + ., scale=t))))
  return(unlist(delta))
}
