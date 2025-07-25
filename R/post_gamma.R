#' @title Posteriori dos Coeficientes de Ondaletas Considerando Erro Gamma
#'
#' @export
#'
#' @importFrom wavethresh GenW
#'
#' @param theta vetor
#' @param d coeficientes empíricos.
#' @param beta,lambda parâmetro da \eqn{Gamma(\beta, \lambda)}.
#' @param tau parâmetro da logística.
#' @param alpha coeficiente da spike and slab (0.8 é o valor padrão).
#' @inheritParams wavethresh::wd
#' @param suppressMessage Se `TRUE` (default) suprime a mensagem "Fora do
#'   suporte." caso o ponto `theta` esteja fora do domínio.
#' @param W Matriz que carrega a DWT. Por padrão ela é gerada dentro da função,
#'   contudo, para uma simulação é interessante gerá-la previamente para
#'   melhorar o desemepenho.
#'
#' @details Partindo do modelo vetorial no dominio do tempo
#' \deqn{
#' \boldsymbol y = \boldsymbol f + \boldsymbol e
#' }
#'   com \eqn{e_i \sim Gamma(\beta, \lambda)}. Então, no domínio das ondaletas,
#'   temos que o modelo é dado por
#' \deqn{
#' \boldsymbol d = \boldsymbol\theta + \boldsymbol\varepsilon}
#'   onde \eqn{\boldsymbol d = W \boldsymbol y} são os coeficiente empíricos,
#'   \eqn{\boldsymbol\theta = W \boldsymbol f} são os coeficientes de ondaleta
#'   verdadeiros e \eqn{\boldsymbol \varepsilon = W \boldsymbol e} é o erro no
#'   domínio das ondaletas Além disso, atribuindo a priori dos coeficientes de
#'   ondaleta \eqn{\boldsymbol\theta} consiederado que seus elementos são
#'   independentes, tem-se:
#' \deqn{
#' \displaystyle\pi(\boldsymbol\theta; \alpha, \tau) = \prod_{i=1}^n[\alpha
#' \delta_0(\theta_i) + (1 - \alpha) g(\theta_i; \tau)]
#' }
#'   onde \eqn{\alpha \in (0,1)}, \eqn{\delta_0} é o delta de Dirac com massa em
#'   0 e \eqn{g(\theta; \tau)} é a função densidade de probabilidade logística,
#'   dada por
#' \deqn{
#' g(\theta_i; \tau) = \dfrac{\exp\left\{-\frac{\theta}{\tau}\right\}}{\tau
#' \left( 1 + \exp\left\{-\frac{\theta}{\tau}\right\} \right)^2} \;
#' \mathcal{I}_\mathbb{R}(\theta) \;, \hspace{0.45cm} \tau > 0
#' }
#'   Portanto, com isso é possível obter a verossimilhança e computar a
#'   posteriori, dada por
#' \deqn{
#' \displaystyle \pi(\boldsymbol\theta \mid \boldsymbol d) \propto \prod_{i=1}^n
#' \left[ \alpha \delta_0(\theta_i) + (1 - \alpha)
#' \frac{\exp\left\{-\frac{\theta_i}{\tau}\right\}}{\tau \left( 1 +
#' \exp\left\{-\frac{\theta_i}{\tau}\right\} \right)^2} \right] \\ \times
#' \exp\left\{-\lambda \sum_{i=1}^n\sum_{k=1}^n w_{ki} (d_k -\theta_k)\right\}
#' \left( \prod_{i=1}^n\sum_{k=1}^n w_{ki} (d_k -\theta_k) \right)^{\beta-1} \\
#' \times \prod_{i=1}^n \mathcal{I}_{(0,\infty)}\left(\sum_{k=1}^n w_{ki} (d_k
#' -\theta_k)\right)
#' }
#'
#' @references Sousa, A.R.S., Garcia, N.L (2023): Wavelet shrinkage in
#' nonparametric regression models with positive noise. \emph{Journal of
#' Statistical Computation and Simulation}, DOI:
#'   [10.1080/00949655.2023.2215372](https://doi.org/10.1080/00949655.2023.2215372).

post_gamma <- function(theta, d, beta, lambda, tau, alpha=0.8, filter.number=5,
                       family='DaubExPhase', W=NULL, suppressMessage=T) {
  if (class(d) == 'wd') d_emp <- c(accessC(d, lev=0), d$D)
  else d_emp <- d
  if (is.null(W)) W <- t(GenW(n=length(theta), filter.number, family))
  wdt_i <- t(W) %*% as.vector((d_emp - theta))
  if (all(wdt_i > 0)) {
    return(prod(ifelse(theta == 0, alpha, 0) + (1 - alpha) * dlogis(theta, scale=tau)) *
             exp(-lambda * sum(wdt_i)) * (prod(wdt_i))^(beta - 1))
  } else {
    if (suppressMessage == F) message('Fora do suporte.')
    return(0)
  }
}
