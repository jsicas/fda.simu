# Este arquivo fará todo o processamento da simulação para o modelo de dados
# funcionais com erro Gamma(beta, lambda). Como parâmetros temos:
#
#   - M = 32 e 128 (pontos por curva);
#   - L = 2 e 4 (quantidade de funções componentes);
#   - SNR = 3 (beta = 36/49 e lambda = 18/49) e 9 (beta = 324/49 e lambda = 162/49);
#   - I = 50 (número de curvas);
#   - alpha = 0.75 e tau = 5 (parâmetros da mistura da priori);
#   - n_ite = 100000 (tamanho da cadeia);
#   - rep = 500 (quantidade de replicações).
#
# Além disso, para L = 2 serão utilizadas as funções Bumps e Doppler, para L = 4
# serão utilziadas todas as funções de testes de Donoho e Johnstone. Por fim,
# vale ressaltar que a esperança do erro será fixada em 2 para todos os
# cenários. Com relação às ondaletas, foi utilizado filter.number=5 e
# family='DaubExPhase'.
#
# Referências
#   Sousa, A.R.S. (2020). Bayesian wavelet shrinkage with logistic prior.
# Communications in Statistics - Simulation and Computation. DOI:
# 10.1080/03610918.2020.1747076.
#
#   Sousa, A.R.S. (2024). A wavelet-based method in aggregated functional data
# analysis. Monte Carlo Methods and Applications, 30(1), 19-30. DOI:
# 10.1515/mcma-2023-2016.


# setup -------------------------------------------------------------------
set.seed(1832323)
require(fda.simu)
require(future)
plan(multisession, workers=25)

# a estrutura dos objetos da simulação vão seguir: simu_SNR_M_L

# M = 32 ------------------------------------------------------------------
M <- 32
alpha_2 <- matrix(c(f_test(M)$bumps, f_test(M)$doppler), ncol=2)
alpha_4 <- matrix(c(f_test(M)$bumps, f_test(M)$doppler, f_test(M)$blocks,
                    f_test(M)$heavisine), ncol=4)

# SNR = 3 ===============================
simu_3_32_2 <- simu_ram_gamma(alpha_2, I=50, n_ite=100000, rep=500, alpha=0.75,
                              tau=5, beta=36/49, lambda=18/49)
simu_3_32_4 <- simu_ram_gamma(alpha_4, I=50, n_ite=100000, rep=500, alpha=0.75,
                              tau=5, beta=36/49, lambda=18/49)

print('SNR = 3 e M = 32: feito')

# SNR = 9 ===============================
simu_9_32_2 <- simu_ram_gamma(alpha_2, I=50, n_ite=100000, rep=500, alpha=0.75,
                              tau=5, beta=324/49, lambda=162/49)
simu_9_32_4 <- simu_ram_gamma(alpha_4, I=50, n_ite=100000, rep=500, alpha=0.75,
                              tau=5, beta=324/49, lambda=162/49)

print('SNR = 9 e M = 32: feito')

# M = 128 -----------------------------------------------------------------
M <- 128
alpha_2 <- matrix(c(f_test(M)$bumps, f_test(M)$doppler), ncol=2)
alpha_4 <- matrix(c(f_test(M)$bumps, f_test(M)$doppler, f_test(M)$blocks,
                    f_test(M)$heavisine), ncol=4)

# SNR = 3 ===============================
simu_3_32_2 <- simu_ram_gamma(alpha_2, I=50, n_ite=100000, rep=500, alpha=0.75,
                              tau=5, beta=36/49, lambda=18/49)
simu_3_32_4 <- simu_ram_gamma(alpha_4, I=50, n_ite=100000, rep=500, alpha=0.75,
                              tau=5, beta=36/49, lambda=18/49)

print('SNR = 3 e M = 128: feito')

# SNR = 9 ===============================
simu_9_32_2 <- simu_ram_gamma(alpha_2, I=50, n_ite=100000, rep=500, alpha=0.75,
                              tau=5, beta=324/49, lambda=162/49)
simu_9_32_4 <- simu_ram_gamma(alpha_4, I=50, n_ite=100000, rep=500, alpha=0.75,
                              tau=5, beta=324/49, lambda=162/49)

print('SNR = 9 e M = 128: feito')

save.image('simulacao2_FAPESP.RData')
