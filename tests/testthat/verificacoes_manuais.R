require(wavethresh)
require(fda.simu)


#---- sample_gen --------------------
set.seed(123)
fun_comp2 <- matrix(c(f_test()$bumps, f_test()$doppler), 2, byrow=T)
n <- 4; snr <- 5; L <- nrow(fun_comp2)
y <- apply(matrix(runif(n*L), nrow=L), 2, \(x) x/sum(x))
fun_agr2 <- t(y) %*% fun_comp2 + rnorm(n*ncol(fun_comp2), 0, 7/snr)

set.seed(123)
sample <- sample_gen(fun_comp2, snr, n, stand=F)

identical(sample$fun, fun_agr2)  # verificando função recuperada
identical(sample$y, y)           # verificando pesos


#---- desagrega --------------------
filter.number <- 3
D <- apply(fun_agr2, MARGIN=1, wd, filter.number, family='DaubExPhase')
D_shrink <- sapply(D, \(x = threshold(D, policy='sure')) c(accessC(x, lev=0), x$D))
gamma <- D_shrink %*% t(y) %*% solve(y %*% t(y))
alpha <- t(GenW(ncol(fun_agr2), filter.number, family='DaubExPhase') %*% gamma)

alpha_recup <- desagrega(fun_agr2, y, filter.number=filter.number)

identical(alpha, alpha_recup)


#---- simu --------------------
set.seed(123)
n <- 4; snr <- 5; L <- nrow(fun_comp2)
fun_comp2 <- matrix(c(f_test()$bumps, f_test()$doppler), 2, byrow=T)
sample <- sample_gen(fun_comp2, snr, n, stand=F)

set.seed(123)
val <- future_map(1:1, ~{
  # gerando amostra
  sample <- sample_gen(fun_comp2, snr, n, stand=F, 7)
  # wavelets
  alpha <- desagrega(sample$fun, sample$y, policy='sure',
                     filter.number=filter.number)
  # calculando erro
  MSE <- rowMeans((alpha - fun_comp2)^2)
  c(MSE, mean(MSE))
}, .options=furrr_options(seed=TRUE), .progress=T) |>
  do.call(rbind.data.frame, args=_) |>
  setNames(nm=c(paste0('MSE_', 1:nrow(fun_comp2)), 'AMSE'))

set.seed(123)
simu <- simu(fun_comp2, snr=snr, rep=1, n=4, filter.number=filter.number)

identical(simu, val)
