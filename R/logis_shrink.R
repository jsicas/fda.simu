#' @title Regra de Encolhimento Bayesiana com Priori Logística
#'
#' @export
#'
#' @usage logis_shrink(d, a, s, t)
#'
#' @param d vetor de coeficientes de ondaleta empírico.
#' @param a parâmetro \eqn{\alpha} da mistura.
#' @param s desvio padrão dos coeficientes de ondaleta.
#' @param t parâmetro \eqn{\tau} da logística.
#'
#' @returns Retorna um vetor de mesmo tamanho inicial com os coeficientes de
#' ondaleta estimados.
#'
#' @details
#' Seja \eqn{\theta} uma variável aleatória com distribuição logística, então,
#' sua densidade é dada por:
#' \deqn{g(\theta; \tau) = \dfrac{\exp\left\{-\frac{\theta}{\tau}\right\}}
#' {\tau\left(1 + \exp\left\{-\frac{\theta}{\tau}\right\}\right)^2} \;
#' \mathcal{I}_{\mathbb{R}}^{(\theta)}}
#' onde \eqn{\tau > 0} é um parâmetro da distribuição.
#'
#' Além disso, a priori utilizada é
#' \deqn{\pi(\theta; \tau, \alpha) = \alpha\delta_0(\theta) + (1-\alpha) g(\theta; \tau)}
#' onde \eqn{\alpha \in (0,1)}, \eqn{\delta_0} é o delta de Dirac centrado em
#' \eqn{0} e \eqn{g(\theta; \tau)} é a densidade da logística.
#'
#' Com isso, obtem-se a regra de shrinkage Bayesiana:
#' \deqn{\delta(d) = \dfrac{(1-\alpha)\int_\mathbb{R} (\sigma u + d)
#' g(\sigma u + d; \tau) \phi(u) du}{\frac{\alpha}{\sigma} \phi\left(
#' \frac{d}{\sigma}\right) + (1-\alpha) \int_\mathbb{R}g(\sigma u + d; \tau)
#' \phi(u) du}}
#' onde \eqn{\phi} é a densidade da distribuição normal padrão.
#'
#' @references
#' Sousa, A.R.S. (2020). Bayesian wavelet shrinkage with logistic prior.
#' \emph{Communications in Statistics - Simulation and Computation}, 51(8),
#' 4700–4714. [doi:10.1080/03610918.2020.1747076](https://doi.org/10.1080/03610918.2020.1747076).
#'
#' @examples
#' x <- seq(1,10)
#' logis_shrink(x, 0.8, 1, 10)
#'
#' x <- seq(-7, 7, 0.05)
#' t <- c(3, 5, 10, 20, 30, 40, 50)
#' plot(x=1, main='Variando t', type='n', xlab='', ylab='',
#'      xlim=c(-6,6), ylim=c(-6,6))
#' abline(v=c(-3, 0, 3), h=0, lty=c(1,9,1,9))
#' for (i in 1:length(t)) {
#'   lines(x, logis_shrink(x, 0.8, 1, t[i]), col=i, lwd=2)
#' }

logis_shrink <- function(d, a, s, t) {
  u <- rnorm(10000)
  delta <- vector(length(d), mode='double')
  for(i in 1:length(d)) {
    logis <- dlogis(s*u + d[i], scale=t)
    int1 <- mean((s*u + d[i]) * logis)
    int2 <- mean(logis)
    delta[i] <- (1-a) * int1/(a * dnorm(d[i], sd=s)/s + (1-a) * int2)
  }
  return(delta)
}
