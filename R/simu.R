#' @title Simulação para Política de Limiar
#'
#' @export
#'
#' @importFrom furrr future_map
#' @importFrom furrr furrr_options
#' @importFrom future nbrOfWorkers
#'
#' @usage
#' simu(fun_comp, snr, rep, n=10, policy='sure', filter.number=10,
#'    stand=F, signal=7)
#'
#' @param rep quantidade de replicações da simulação.
#' @inheritParams sample_gen
#' @inheritParams desagrega
#'
#' @returns Retorna um objeto da classe \code{data.frame} com o cálculo do MSE
#' de cada função e o AMSE.
#'
#' @details
#' A função assume que, para o cálculo da razão sinal-ruído, que
#' \eqn{sd(sinal)=7}.
#' Os pesos para cada função são gerados aleatóriamente seguindo uma
#' \eqn{\mathcal{U}(0,1)} e são modificados de forma a garantir que sua soma é
#' \eqn{1}.
#'
#' @references
#' Sousa, A.R.S. (2024). A wavelet-based method in aggregated functional data
#' analysis. \emph{Monte Carlo Methods and Applications}, 30(1), 19-30.
#' [doi:10.1515/mcma-2023-2016](https://doi.org/10.1515/mcma-2023-2016).
#'
#' @examples
#' bumps <- f_test()$bumps
#' doppler <- f_test()$doppler
#' par(mfrow=c(2,1))
#' plot(bumps, type='l', main='Bumps'); plot(doppler, type='l', main='Doppler')
#'
#' # plan(multisession, workers=2)  # paralelizando execução
#' fun_comp <- matrix(c(bumps, doppler), nrow=2, byrow=T)
#' simu(fun_comp, rep=4, n=10, snr=5)

simu <- function(fun_comp, snr, rep, n=10, policy='sure', filter.number=10,
                 stand=F, signal=7) {
  if(nbrOfWorkers() == 1)
    message('\nCuidado: a simulação não está sendo paralelizada.\n')
  future_map(1:rep, ~{
    # gerando amostra
    sample <- sample_gen(fun_comp, snr, n, stand=stand, signal=signal)

    # wavelets
    alpha <- desagrega(sample$fun, sample$y, policy=policy,
                       filter.number=filter.number)
    # calculando erro
    data.frame('MSE'=t(rowMeans((t(alpha) - fun_comp)^2)), 'AMSE'=mean(MSE))
  }, .options=furrr_options(seed=TRUE), .progress=T) |>
    do.call(rbind, args=_) |>
    setNames(nm=c(paste0('MSE_', 1:nrow(fun_comp)), 'AMSE'))
}
