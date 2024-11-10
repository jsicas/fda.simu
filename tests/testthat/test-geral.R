#---- configurações iniciais --------------------
set.seed(282829)
db_test <- matrix(c(f_test()$bumps, f_test()$doppler), nrow=2, byrow=T)


#---- testes ------------------------------------
test_that('simu: testando duas curvas', {
  expect_snapshot_value(simu(db_test[1:2,], 4, snr=5), style='json2')
})

# test_that('simu: testando 4 curvas', {
#   expect_snapshot_value(simu(db_test, 4, snr=5), style='json2')
# })
