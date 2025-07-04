#' @title Funções para Teste
#'
#' @export
#'
#' @usage f_test(n=1024, signal=7, snr=7, noise=FALSE)
#'
#' @param n números de pontos considerados da função.
#' @param signal desvio padrão da função verdadeira.
#' @param snr razão sinal-ruído.
#' @param noise indica se os pontos gerados devem apresentar um ruido gaussiano.
#'
#' @details
#' Todas as funções são normalizadas para apresentarem o mesmo desvio padrão
#' \code{signal}, cujo default é \eqn{7}.
#'
#' @returns Retorna uma lista com \code{n} pontos de cada função.
#'
#' @references
#' Donoho, D.L., Johnstone, I.M. (1994). Ideal spatial adaptation by wavelet
#' shrinkage. \emph{Biometrika}, 81, 425–455.

f_test <- function(n=1024, signal=7, snr=7, noise=FALSE)
{
  x <- seq(1, n)/n
  t <- c(0.1, 0.13, 0.15, 0.23, 0.25, 0.4, 0.44, 0.65, 0.76, 0.78, 0.81)

  # bumps
  h.bumps <- c(4, 5, 3, 4, 5, 4.2, 2.1, 4.3, 3.1, 5.1, 4.2)
  w <- c(0.005, 0.005, 0.006, 0.01, 0.01, 0.03, 0.01, 0.01, 0.005, 0.008, 0.005)
  bumps <- rep(0, n)
  for (i in seq(1, length(h.bumps))) {
    bumps <- bumps + h.bumps[i] * pmax(0, (1 - abs((x - t[i])/w[i])))^4
  }

  # blocks
  h.blocks <- c(4, -5, 3, -4, 5, -4.2, 2.1, 4.3, -3.1, 2.1, -4.2)
  blocks <- rep(0, n)
  for (i in seq(1, length(h.blocks))) {
    blocks <- blocks + (h.blocks[i] * (1 + sign(x - t[i])))/2
  }

  # doppler
  doppler <- sqrt(x * (1 - x)) * sin(2 * pi * (1 - 0.05) / (x + 0.05))  # wavethresh
  # doppler <- sqrt(x * (1 - x)) * sin(2 * pi * (1 + 0.05) / (x + 0.05))  # Alex

  # heavisine
  heavisine <- 4 * sin(4 * pi * x) - sign(x - 0.3) - sign(0.72 - x)

  # logit
  logit <- 1/(1 + exp(-20 * (x - 0.5)))

  # SpaHet
  spahet <- sqrt(x * (1 - x)) * sin(2 * pi * (1 + 2^(-0.6))/(x + 2^(-0.6)))

  # example.1 - piecewise polynomial
  example <- rep(0, length(x))
  i1 <- x <= 0.5
  i2 <- (x > 0.5) & (x <= 0.75)
  i3 <- x > 0.75
  example[i1] <- -16 * x[i1]^3 + 12 * x[i1]^2
  example[i2] <- (x[i2] * (16 * x[i2]^2 - 40 * x[i2] + 28))/3 - 1.5
  example[i3] <- (x[i3] * (16 * x[i3]^2 - 32 * x[i3] + 16))/3

  # padronizando
  bumps <- bumps/sd(bumps) * signal
  blocks <- blocks/sd(blocks) * signal
  doppler <- doppler/sd(doppler) * signal
  heavisine <- heavisine/sd(heavisine) * signal
  logit <- logit/sd(logit) * signal
  spahet <- spahet/sd(spahet) * signal
  example <- example/sd(example) * signal

  # noise
  if (noise == TRUE) {
    result <- list('bumps'=bumps + rnorm(n, 0, signal/snr),
                   'blocks'=blocks + rnorm(n, 0, signal/snr),
                   'doppler'=doppler + rnorm(n, 0, signal/snr),
                   'heavisine'=heavisine + rnorm(n, 0, signal/snr),
                   'logit'=logit + rnorm(n, 0, signal/snr),
                   'spahet'=spahet + rnorm(n, 0, signal/snr),
                   'example'=example + rnorm(length(example), 0, signal/snr))
  }
  else {
    result <- list('x'=x, 'bumps'=bumps, 'blocks'=blocks, 'doppler'=doppler,
                   'heavisine'=heavisine, 'logit'=logit, 'spahet'=spahet,
                   'example'=example)
  }
  return(result)
}
