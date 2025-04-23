#' @title Gera um Ponto no Suprote da Posteriori com Erro Gamma
#'
#' @export
#'
#' @importFrom wavethresh GenW
#' @importFrom wavethresh accessC
#'
#' @param n Quantidade de pontos por função.
#' @param y Dados.
#' @param lim_sup Valor máximo gerado para o ponto gerado a.
#' @param ... configurações referênte à DWT, como `filter.number` e `family`.

gera_ponto <- function(n, y, lim_sup = 20, filter.number = 5,
                       family = 'DaubExPhase') {
  dwt <- wd(y, filter.number=filter.number, family=family)
  d <- c(accessC(dwt, lev=0), dwt$D)  # coeficientes empíricos
  a <- runif(n, 0.01, lim_sup)        # W(d - theta) = a >= 0
  W <- t(GenW(n, filter.number=filter.number, family=family))
  theta <- d - W %*% a
  return(theta)
}
