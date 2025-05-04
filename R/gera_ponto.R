#' @title Gera um Ponto no Suprote da Posteriori com Erro Gamma
#'
#' @export
#'
#' @importFrom wavethresh GenW
#' @importFrom wavethresh accessC
#'
#' @param y Valores observados.
#' @param lim_sup Valor máximo para uma coordenada gerado para o ponto gerado a.
#' @inheritParams wavethresh::wd
#'
#' @description
#' Gera um ponto no domínio do tempo da função original garantindo que esta
#' dentro do suporte da posteriori com erro Gamma [post_gamma()].
#'
#' @examples
#' y <- c(1,2,7,2,5,6,0,3)  # valores observados
#' d <- wd(y, filter.number = 5, family = 'DaubExPhase')
#' (theta <- gera_ponto(y))
#' post_gamma(theta, d, beta = 5, lambda = 1, tau = 4, alpha = 0.7)  # posteriori

gera_ponto <- function(y, lim_sup = 20, filter.number = 5,
                       family = 'DaubExPhase') {
  dwt <- wd(y, filter.number=filter.number, family=family)
  d <- c(accessC(dwt, lev=0), dwt$D)          # coeficientes empíricos
  a <- runif(length(y), 0.01, lim_sup)        # W(d - theta) = a >= 0
  W <- t(GenW(length(y), filter.number=filter.number, family=family))
  theta <- d - W %*% a
  return(theta)
}
