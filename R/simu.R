#' @title Simulação para Política de Limiar
#'
#' @export
#'
#' @importFrom furrr future_map
#' @importFrom furrr furrr_options
#' @importFrom future nbrOfWorkers
#'
#' @param rep quantidade de replicações da simulação.
#' @inheritParams sample_gen
#' @inheritParams desagrega
#'
#' @returns
#' Retorna um objeto da classe \code{data.frame} com o cálculo do MSE de cada
#' função e o AMSE.
#'
#' @details
#' A função assume, para o cálculo da razão sinal-ruído, que \eqn{sd(sinal)=7}.
#' Os pesos para cada função são gerados aleatóriamente seguindo uma
#' \eqn{\mathcal{U}(0,1)} e são modificados de forma a garantir que sua soma
#' seja \eqn{1}.
#'
#' @references
#' Sousa, A.R.S. (2024). A wavelet-based method in aggregated functional data
#' analysis. \emph{Monte Carlo Methods and Applications}, 30(1), 19-30. DOI:
#' [10.1515/mcma-2023-2016](https://doi.org/10.1515/mcma-2023-2016).
#'
#' @examples
#' x <- (1:1024)/1024
#' fun_comp <- matrix(c(f_test()$bumps, f_test()$doppler), ncol=2)
#' par(mfrow=c(2,1))
#' plot(x, fun_comp[,1], type='l', main='Bumps', ylab='y')
#' plot(x,fun_comp[,2], type='l', main='Doppler', ylab='y')
#'
#' # plan(multisession, workers=2)  # paralelizando execução
#' simu(fun_comp, rep=4, I=10, snr=5)

simu <- function(alpha_comp, snr, rep, I=10, policy='sure',
                 filter.number=5, family='DaubExPhase') {
  if (nbrOfWorkers() == 1)
    message('A simulação não está sendo paralelizada.')

  future_map(1:rep, ~{
    # gerando amostra
    sample <- sample_gen(alpha_comp, snr=snr, I=I, erro='normal')
    # wavelets
    alpha <- desagrega(sample$fun, sample$y, policy=policy,
                       filter.number=filter.number, family=family)
    # calculando erro
    MSE <- colMeans((alpha - alpha_comp)^2)
    c(MSE, mean(MSE))
  }, .options=furrr_options(seed=TRUE), .progress=T) |>
    do.call(rbind.data.frame, args=_) |>
    setNames(nm=c(paste0('MSE_', 1:ncol(alpha_comp)), 'AMSE'))
}
