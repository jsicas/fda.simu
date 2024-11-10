#' @title DEPRECIADO: Simulação2 para Política de Limiar
#'
#' @importFrom furrr future_map
#' @importFrom furrr furrr_options
#' @importFrom future nbrOfWorkers

simu2 <- function(fun_comp, rep, n=10, snr, policy='sure', filter.number=10) {
  warning('Esta função está depreciada.')
  if(nbrOfWorkers() == 1) warning('Simulação não está sendo paralelizada.')
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
