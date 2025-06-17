#' @title Desagregador Dados Funcionais
#'
#' @export
#'
#' @importFrom wavethresh wd
#' @importFrom wavethresh GenW
#' @importFrom wavethresh accessC
#' @importFrom wavethresh threshold
#'
#' @param data matriz contendo, em cada coluna, uma observações.
#' @param y matriz com os pesos conhecidos de cada funcional.
#' @param policy política para escolha de limiar. Os possíveis valores são
#'   "sure", "universal", "cv", "fdr", "logistica", "epanechnikov". Para mais
#'   detalhes, ver \code{\link[wavethresh]{wavethresh::threshold.wd}}.
#' @param filter.number controla a suavidade da ondaleta.
#' @param family Especifíca a família da ondaleta, e.g. "DaubExPhase",
#'   "DaubLeAsymm".
#' @inheritParams logis_shrink
#'
#' @returns Retorna uma matriz com a função recuperar no domínio do tempo de
#'   forma que cada coluna representa uma função recuperada.
#'
#' @references Sousa, A.R.S. (2024). A wavelet-based method in aggregated
#' functional data analysis. \emph{Monte Carlo Methods and Applications}, 30(1),
#' 19-30. DOI: [10.1515/mcma-2023-2016](https://doi.org/10.1515/mcma-2023-2016).
#'
#' @examples
#' set.seed(282829)
#' x <- (1:1024)/1024
#' alpha_comp <- matrix(c(f_test()$bumps, f_test()$doppler), ncol=2)
#' sample <- sample_gen(alpha_comp, snr=5, I=10)
#' plot(sample)
#' fun_recup <- desagrega(sample$fun, sample$y)  # recuperando funções
#'
#' # plot
#' par(mar=c(3, 4, 2, 2), mgp=c(2, 0.5, 0))
#' layout(matrix(c(1,1,2,3), 2, byrow=T))
#' plot(sample)
#' plot(x, fun_recup[,1], main='Bumps', ylab='y', col='blue', type='l')
#' lines(x, alpha_comp[,1])
#' legend('topright', lwd=2, bty='n', cex=0.85, col=c('black', 'blue'),
#'        legend=c('Função Verdadeira', 'Função Recuperada'))
#' plot(x, fun_recup[,2], main='Doppler', ylab='y', col='blue', type='l')
#' lines(x, alpha_comp[,2])
#' legend('bottomright', lwd=2, bty='n', cex=0.85, col=c('black', 'blue'),
#'        legend=c('Função Verdadeira', 'Função Recuperada'))

desagrega <- function(data, y, policy='sure', filter.number=10,
                      family='DaubExPhase', a,s,tau) {
  D <- apply(data, MARGIN=2, wd, filter.number, family)  # coeficientes empiricos
  if (policy %in% c('universal', 'sure', 'cv', 'fdr')) {
    D_shrink <- sapply(D, \(d) {
      thr <- threshold(d, policy=policy)
      c(accessC(thr, lev=0), thr$D)
    })
  } else if (policy == 'logistica') {
    if (missing(a) | missing(s) | missing(tau))
      stop('Especifique os hiperparâmetros da priori.')
    D_shrink <- sapply(D, \(x) c(accessC(x, lev=0), logis_shrink(x$D,a,s,tau)))
  } else if (policy == 'epanechnikov') {
    stop('Regra não implementada!')
  } else {
    stop('Politica de limiar mal especificada.')
  }
  gamma <- D_shrink %*% t(y) %*% solve(y %*% t(y))
  alpha <- GenW(nrow(data), filter.number=filter.number, family=family) %*%
    gamma  # funções recuperadas
  return(alpha)
}
