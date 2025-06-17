# setup -------------------------------------------------------------------
I <- 5     # número de observações
M <- 1024  # número de pontos por observação
snr <- 5   # erro normal
beta <- 2; lambda <- 2  # erro gamma
alpha_comp <- matrix(c(f_test(M)$bumps, f_test(M)$doppler, f_test(M)$blocks), M)
L <- ncol(alpha_comp)  # quantidade de curvas componentes
set.seed(123123)
y <- matrix(runif(I*L), nrow=L) |> apply(2, \(col) col/sum(col))
sample_norm <- alpha_comp %*% y + rnorm(I*M, 0, 7/snr)
set.seed(123123)
y <- matrix(runif(I*L), nrow=L) |> apply(2, \(col) col/sum(col))
sample_gamma <- alpha_comp %*% y + rgamma(I*M, shape=beta, rate=lambda)

D <- apply(sample_norm, MARGIN=2, wd, filter.number=5, family='DaubExPhase')
D_shrink <- matrix(nrow=M, ncol=I)
for (j in seq_along(D)) {
  thr <- threshold(D[[j]], policy='sure')
  D_shrink[,j] <- c(accessC(thr, lev=0), thr$D)
}
gamma <- D_shrink %*% t(y) %*% solve(y %*% t(y))
alpha <- GenW(M, filter.number=5, family='DaubExPhase') %*% gamma


# sample_gen --------------------------------------------------------------
test_that('sample_gen: erro normal', {
  set.seed(123123)
  sample <- sample_gen(alpha_comp, snr=5, I=5, erro='normal')
  expect_equal(sample$fun, sample_norm)
  expect_equal(sample$y, y)
})

test_that('sample_gen: erro gamma', {
  set.seed(123123)
  sample <- sample_gen(alpha_comp, I=5, beta=2, lambda=2, erro='gamma')
  expect_equal(sample$fun, sample_gamma)
  expect_equal(sample$y, y)
})


# desagrega ---------------------------------------------------------------
test_that('desagrega: sure policy', {
  set.seed(123123)
  sample <- sample_gen(alpha_comp, snr=5, I=5, erro='normal')
  des <- desagrega(sample$fun, sample$y, policy='sure', filter.number=5,
                   family='DaubExPhase')
  expect_equal(des, alpha)
})

