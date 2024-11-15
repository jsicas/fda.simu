#---- configurações iniciais --------------------

set.seed(282829)  # garantindo reprodutibilidade

## gerando objetos e amostra
snr <- 5; n <- 4
fun_comp_test <- matrix(c(f_test()$bumps, f_test()$doppler), nrow=2, byrow=T)
L <- nrow(fun_comp_test)  # número de curvas componentes
y <- apply(matrix(runif(L*n), nrow=L), 2, \(col) col/sum(col))  # soma dos pesos igual a 1
fun_agr_test <- t(y) %*% fun_comp_test + rnorm(n*ncol(fun_comp_test), 0, 7/snr)
D <- apply(fun_agr_test, MARGIN=1, wd, family='DaubExPhase')
D_shrink <- sapply(D, \(x = threshold(D, policy='sure')) c(accessC(x, lev=0), x$D))
gamma <- D_shrink %*% t(y) %*% solve(y %*% t(y))
alpha <- t(GenW(ncol(fun_agr_test), filter.number=10, family='DaubExPhase') %*% gamma)


#---- testes ------------------------------------
test_that('sample_gen: sem padronização', {
  set.seed(282829)
  sample <- sample_gen(fun_comp_test, snr, n, stand=F)

  expect_identical(sample$fun, fun_agr_test)
  expect_identical(sample$y, y)
})


test_that('sample_gen: com padronização', {
  set.seed(282829)
  sample <- sample_gen(fun_comp_test, snr, n, stand=T)

  expect_equal(sample$fun, fun_agr_test)
})


test_that('desagrega', {
  expect_equal(desagrega(fun_agr_test, y), alpha)
})


test_that('simu: testando duas curvas', {
  set.seed(282829)
  comp <- matrix(c(f_test(signal=3)$bumps, f_test()$doppler), nrow=2, byrow=T)
  msg1 <- '\nCuidado: a simulação não está sendo paralelizada.'
  msg2 <- 'As funcoes componentes foram normalizadas (sd(sinal) != signal).'

  expect_message(simu(comp, rep=3, snr=5), msg1)   # mensagem para falta de paralelização
  expect_snapshot_value(suppressMessages(simu(fun_agr_test, rep=3, snr=5)),
                        style='json2')  # funcionalidade geral
  expect_snapshot_value(suppressMessages(simu(comp, rep=3, snr=5)),
                        style='json2')  # funcionalidade para sd(sinal) != signal
})


test_that('simu: testando 4 curvas', {
  set.seed(282829)
  db_test <- matrix(c(f_test()$bumps, f_test()$doppler, f_test()$blocks,
                      f_test()$heavi), nrow=4, byrow=T)

  expect_snapshot_value(suppressMessages(simu(db_test, 4, snr=5)), style='json2')
})
