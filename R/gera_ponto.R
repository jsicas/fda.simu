#' @title Gera um Ponto no Suprote da Posteriori com Erro Gamma
#'
#' @export
#'
#' @importFrom wavethresh GenW
#' @importFrom wavethresh accessC
#' @importFrom wavethresh nlevelsWT
#'
#' @param a vetor com todas as entradas positivas. Se `NULL` é gerado
#'   aleatóriamente.
#' @param d Coeficientes empíricos de ondaletas.
#' @param lim_sup Valor máximo para uma coordenada gerado para o ponto gerado a.
#' @inheritParams wavethresh::wd
#'
#' @description Gera um ponto no domínio do tempo da função original garantindo
#'   que esta dentro do suporte da posteriori com erro Gamma, ver função
#'   [post_gamma].
#'
#' @examples
#' # library(wavethresh)
#' y <- c(1,2,7,2,5,6,0,3)  # valores observados
#' d <- wd(y, filter.number=5, family='DaubExPhase')
#' (theta <- gera_ponto(a=NULL, d))
#' post_gamma(theta, d, beta=5, lambda=1, tau=4, alpha=0.7)  # posteriori

gera_ponto <- function(a=NULL, d, lim_sup=10, filter.number=5,
                       family='DaubExPhase') {
  M <- 2^nlevelsWT(d)
  W <- t(GenW(M, filter.number, family))
  if (is.null(a)) a <- runif(M, 0.01, lim_sup)
  d_emp <- c(accessC(d, lev=0), d$D)  # coeficientes empíricos
  return(as.vector(d_emp - W %*% a))  # theta
}
