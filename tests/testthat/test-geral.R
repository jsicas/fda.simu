#---- configurações iniciais --------------------
set.seed(282829)
db_test <- matrix(c(f_test()$bumps, f_test()$doppler), nrow=2, byrow=T)

# gerando todos os objetos necessários
snr <- 5; n <- 4

## gerando amostra
L <- nrow(db_test)  # número de curvas componentes
y <- apply(matrix(runif(L*n), nrow=L), 2,
           function(col) col/sum(col))  # soma dos pesos igual a 1
fun_agr_test <- t(y) %*% db_test + rnorm(n*ncol(db_test), 0, 7/snr)

set.seed(282829)
teste <- sample_gen(db_test, snr, n, stand=F)


#---- testes ------------------------------------
test_that('simu: testando duas curvas', {
  expect_snapshot_value(simu(db_test[1:2,], 4, snr=5), style='json2')
})

# test_that('simu: testando 4 curvas', {
#   expect_snapshot_value(simu(db_test, 4, snr=5), style='json2')
# })

# Adicionar testes para:
#     sample_gen
#     simu




teste <- function(a, b, c=T) {
  if (c == T) {
    if (a == b) print('=')
    if (a != b) print('!=')
    print(a); print(b)
  } else print('AAAA')
}


t <- function(a, b, ...) {
  do.call(teste, list(b, a, ...))
}

t(2,3, F)







