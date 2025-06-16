# setup -------------------------------------------------------------------
set.seed(123123)

I <- 5     # número de observações
M <- 1024  # número de pontos por observação
alpha_comp <- matrix(c(f_test(M)$bumps, f_test(M)$doppler, f_test(M)$blocks), M)
L <- ncol(alpha_comp)  # quantidade de curvas componentes
y <- matrix(runif(I*L), nrow=L) |> apply(2, \(col) col/sum(col))



# sample_gen --------------------------------------------------------------
test_that('sample_gen', {

})


