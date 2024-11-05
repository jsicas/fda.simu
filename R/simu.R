#' @title Simulação para Política de Limiar
#'
#' @export
#'
#' @importFrom furrr future_map
#' @importFrom furrr furrr_options
#'
#' @param fun_comp objeto com as \eqn{l} funções componentes. Cada linha é uma
#' função.
#' @param rep quantidade de replicações da simulação.
#' @param n tamanho da amostra gerada.
#' @param snr razão sinal-ruído.
#' @inheritParams desagrega
#'
#' @returns Retorna um objeto da classe `data.frame` com o cálculo do MSE de cada
#' função e o AMSE.
#'
#' @details
#' A função assume que os dados de entrada tem o desvio padrão 7 para efetuar o
#' cálculo da razão sinal-ruído. Além disso, os pesos para cada função são gerados
#' aleatóriamente seguindo uma \eqn{\mathcal{U}(0,1)}.
#'
#' @references
#' Sousa, A.R.S. (2024). A wavelet-based method in aggregated functional data
#' analysis. \emph{Monte Carlo Methods and Applications}, 30(1), 19-30.
#' [https://doi.org/10.1515/mcma-2023-2016](https://doi.org/10.1515/mcma-2023-2016).
#'
#' @examples
#' bumps <- f_test()$bumps
#' doppler <- f_test()$doppler
#' par(mfrow=c(2,1))
#' plot(bumps, type='l', main='Bumps'); plot(doppler, type='l', main='Doppler')
#'
#' funcoes_comp <- matrix(c(bumps, doppler), nrow=2, byrow=T)
#' simu(funcoes_comp, rep=4, n=10, snr=5)

simu <- function(fun_comp, rep, n=10, snr, policy='sure', filter.number=10) {
  future_map(1:rep, ~{
    # gerando amostra
    y1 <- runif(n)  # pesos da curva 1
    y2 <- 1 - y1    # pesos da curva 2
    y <- matrix(c(y1, y2), nrow=2, byrow=T)
    fun_agr <- matrix(0, n, ncol(fun_comp))  # pre-alocando memoria
    for (i in 1:n) fun_agr[i,] <- y1[i]*fun_comp[1,] + y2[i]*fun_comp[2,] +
      rnorm(ncol(fun_comp), 0, 7/snr)

    # wavelets
    alpha <- desagrega(fun_agr, y, policy=policy, filter.number=filter.number)

    # calculando erro
    MSE_1 <- sum((alpha[,1] - fun_comp[1,])^2) / ncol(fun_comp)  # MSE da fç componente 1
    MSE_2 <- sum((alpha[,2] - fun_comp[2,])^2) / ncol(fun_comp)  # MSE da fç componente 2
    AMSE <- (MSE_1 + MSE_2) / 2
    data.frame('MSE_1'=MSE_1, 'MSE_2'=MSE_2, 'AMSE'=AMSE)
  }, .options=furrr_options(seed=TRUE), .progress=T) |>
    do.call(rbind, args=_)
}
