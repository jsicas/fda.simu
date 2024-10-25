#' @title Funções para Teste
#' @export
#'
#' @param n números de pontos considerados da função.
#' @param signal desvio padrão da função verdadeira.
#' @param snr razão sinal-ruído.
#' @param noise indica se os pontos gerados devem apresentar um ruido gaussiano.
#'
#' @references
#' Donoho, D.L. and Johnstone, I.M. (1994), Ideal spatial adaptation by wavelet
#' shrinkage. \enph{Biometrika}, 81, 425–455.
#'
#' Sousa, A.R.S. (2020). Bayesian wavelet shrinkage with logistic prior.
#' \emph{Communications in Statistics - Simulation and Computation}.

f_test <- function(n=1024, signal=7, snr=7, noise=FALSE)
{
  x <- seq(1, n)/n
  t <- c(0.1, 0.13, 0.15, 0.23, 0.25, 0.4, 0.44, 0.65, 0.76, 0.78, 0.81)

  # Bumps
  h.bumps <- c(4, 5, 3, 4, 5, 4.2, 2.1, 4.3, 3.1, 5.1, 4.2)
  w <- c(0.005, 0.005, 0.006, 0.01, 0.01, 0.03, 0.01, 0.01, 0.005, 0.008, 0.005)
  bumps <- rep(0, n)
  for (i in seq(1, length(h.bumps))) {
    bumps <- bumps + h.bumps[i] * pmax(0, (1 - abs((x - t[i])/w[i])))^4
  }

  # Blocks
  h.blocks <- c(4, -5, 3, -4, 5, -4.2, 2.1, 4.3, -3.1, 2.1, -4.2)
  blocks <- rep(0, n)
  for (i in seq(1, length(h.blocks))) {
    blocks <- blocks + (h.blocks[i] * (1 + sign(x - t[i])))/2
  }

  # Doppler
  doppler <- sqrt(x * (1 - x)) * sin(2 * pi * (1 - 0.05) / (x + 0.05))  # wavethresh
  # doppler <- sqrt(x * (1 - x)) * sin(2 * pi * (1 + 0.05) / (x + 0.05))  # Alex

  # Heavisine
  heavisine <- 4 * sin(4 * pi * x) - sign(x - 0.3) - sign(0.72 - x)

  # Logit
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

  # Padronizando
  bumps <- bumps/sqrt(var(bumps)) * signal
  blocks <- blocks/sqrt(var(blocks)) * signal
  doppler <- doppler/sqrt(var(doppler)) * signal
  heavisine <- heavisine/sqrt(var(heavisine)) * signal
  logit <- logit/sqrt(var(logit)) * signal
  spahet <- spahet/sqrt(var(spahet)) * signal
  example <- example/sqrt(var(example)) * signal

  # Noise
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
