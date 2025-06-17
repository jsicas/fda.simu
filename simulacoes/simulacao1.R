# Código da simulação 1. Os cenários abordados foram:
#
#   - M: 32, 512, 2048 (quantidade de pontos por curva);
#   - snr: 1, 3 e 5;
#   - L:
#       2 (bumps e doppler);
#       4 (funções de DJ);
#       5 (funções de DJ e logit).
#
# Além disso, também é utilizado I = 50 (tamanho da amostra) e 500 replicações.
# Com relação às ondaletas, foi utilizado filter.number=5 e family='DaubExPhase'
# e a polistica de escolha de limiar 'sure'.


# setup -------------------------------------------------------------------
set.seed(678512)
require(fda.simu)
plan(multisession, workers=5)


# SNR = 1 -----------------------------------------------------------------
# M = 32 ==================================
M <- 32
alpha_comp_2 <- matrix(c(f_test(M)$bumps, f_test(M)$doppler), ncol=2)
alpha_comp_4 <- matrix(c(f_test(M)$bumps, f_test(M)$doppler, f_test(M)$blocks,
                         f_test(M)$heavisine), ncol=4)
alpha_comp_5 <- matrix(c(f_test(M)$bumps, f_test(M)$doppler, f_test(M)$blocks,
                         f_test(M)$heavisine, f_test(M)$logit), ncol=5)

# simulação
# a estrutura dos objetos da simulação vão seguir: simu_SNR_M_L
simu_1_32_2 <- simu(alpha_comp_2, snr=1, rep=500, I=50)  # L = 2
simu_1_32_4 <- simu(alpha_comp_4, snr=1, rep=500, I=50)  # L = 4
simu_1_32_5 <- simu(alpha_comp_5, snr=1, rep=500, I=50)  # L = 5

# M = 512 =================================
M <- 512
alpha_comp_2 <- matrix(c(f_test(M)$bumps, f_test(M)$doppler), ncol=2)
alpha_comp_4 <- matrix(c(f_test(M)$bumps, f_test(M)$doppler, f_test(M)$blocks,
                         f_test(M)$heavisine), ncol=4)
alpha_comp_5 <- matrix(c(f_test(M)$bumps, f_test(M)$doppler, f_test(M)$blocks,
                         f_test(M)$heavisine, f_test(M)$logit), ncol=5)

# simulacao
simu_1_512_2 <- simu(alpha_comp_2, snr=1, rep=500, I=50)  # L = 2
simu_1_512_4 <- simu(alpha_comp_4, snr=1, rep=500, I=50)  # L = 4
simu_1_512_5 <- simu(alpha_comp_5, snr=1, rep=500, I=50)  # L = 5

# M = 2048 ================================
M <- 2048
alpha_comp_2 <- matrix(c(f_test(M)$bumps, f_test(M)$doppler), ncol=2)
alpha_comp_4 <- matrix(c(f_test(M)$bumps, f_test(M)$doppler, f_test(M)$blocks,
                         f_test(M)$heavisine), ncol=4)
alpha_comp_5 <- matrix(c(f_test(M)$bumps, f_test(M)$doppler, f_test(M)$blocks,
                         f_test(M)$heavisine, f_test(M)$logit), ncol=5)

# simulacao
simu_1_2048_2 <- simu(alpha_comp_2, snr=1, rep=500, I=50)  # L = 2
simu_1_2048_4 <- simu(alpha_comp_4, snr=1, rep=500, I=50)  # L = 4
simu_1_2048_5 <- simu(alpha_comp_5, snr=1, rep=500, I=50)  # L = 5


# SNR = 3 -----------------------------------------------------------------
# M = 32 ==================================
M <- 32
alpha_comp_2 <- matrix(c(f_test(M)$bumps, f_test(M)$doppler), ncol=2)
alpha_comp_4 <- matrix(c(f_test(M)$bumps, f_test(M)$doppler, f_test(M)$blocks,
                         f_test(M)$heavisine), ncol=4)
alpha_comp_5 <- matrix(c(f_test(M)$bumps, f_test(M)$doppler, f_test(M)$blocks,
                         f_test(M)$heavisine, f_test(M)$logit), ncol=5)

# simulação
simu_3_32_2 <- simu(alpha_comp_2, snr=3, rep=500, I=50)  # L = 2
simu_3_32_4 <- simu(alpha_comp_4, snr=3, rep=500, I=50)  # L = 4
simu_3_32_5 <- simu(alpha_comp_5, snr=3, rep=500, I=50)  # L = 5

# M = 512 =================================
M <- 512
alpha_comp_2 <- matrix(c(f_test(M)$bumps, f_test(M)$doppler), ncol=2)
alpha_comp_4 <- matrix(c(f_test(M)$bumps, f_test(M)$doppler, f_test(M)$blocks,
                         f_test(M)$heavisine), ncol=4)
alpha_comp_5 <- matrix(c(f_test(M)$bumps, f_test(M)$doppler, f_test(M)$blocks,
                         f_test(M)$heavisine, f_test(M)$logit), ncol=5)

# simulacao
simu_3_512_2 <- simu(alpha_comp_2, snr=3, rep=500, I=50)  # L = 2
simu_3_512_4 <- simu(alpha_comp_4, snr=3, rep=500, I=50)  # L = 4
simu_3_512_5 <- simu(alpha_comp_5, snr=3, rep=500, I=50)  # L = 5

# M = 2048 ================================
M <- 2048
alpha_comp_2 <- matrix(c(f_test(M)$bumps, f_test(M)$doppler), ncol=2)
alpha_comp_4 <- matrix(c(f_test(M)$bumps, f_test(M)$doppler, f_test(M)$blocks,
                         f_test(M)$heavisine), ncol=4)
alpha_comp_5 <- matrix(c(f_test(M)$bumps, f_test(M)$doppler, f_test(M)$blocks,
                         f_test(M)$heavisine, f_test(M)$logit), ncol=5)

# simulacao
simu_3_2048_2 <- simu(alpha_comp_2, snr=3, rep=500, I=50)  # L = 2
simu_3_2048_4 <- simu(alpha_comp_4, snr=3, rep=500, I=50)  # L = 4
simu_3_2048_5 <- simu(alpha_comp_5, snr=3, rep=500, I=50)  # L = 5


# SNR = 5 -----------------------------------------------------------------
# M = 32 ==================================
M <- 32
alpha_comp_2 <- matrix(c(f_test(M)$bumps, f_test(M)$doppler), ncol=2)
alpha_comp_4 <- matrix(c(f_test(M)$bumps, f_test(M)$doppler, f_test(M)$blocks,
                         f_test(M)$heavisine), ncol=4)
alpha_comp_5 <- matrix(c(f_test(M)$bumps, f_test(M)$doppler, f_test(M)$blocks,
                         f_test(M)$heavisine, f_test(M)$logit), ncol=5)

# simulação
simu_5_32_2 <- simu(alpha_comp_2, snr=5, rep=500, I=50)  # L = 2
simu_5_32_4 <- simu(alpha_comp_4, snr=5, rep=500, I=50)  # L = 4
simu_5_32_5 <- simu(alpha_comp_5, snr=5, rep=500, I=50)  # L = 5

# M = 512 =================================
M <- 512
alpha_comp_2 <- matrix(c(f_test(M)$bumps, f_test(M)$doppler), ncol=2)
alpha_comp_4 <- matrix(c(f_test(M)$bumps, f_test(M)$doppler, f_test(M)$blocks,
                         f_test(M)$heavisine), ncol=4)
alpha_comp_5 <- matrix(c(f_test(M)$bumps, f_test(M)$doppler, f_test(M)$blocks,
                         f_test(M)$heavisine, f_test(M)$logit), ncol=5)

# simulacao
simu_5_512_2 <- simu(alpha_comp_2, snr=5, rep=500, I=50)  # L = 2
simu_5_512_4 <- simu(alpha_comp_4, snr=5, rep=500, I=50)  # L = 4
simu_5_512_5 <- simu(alpha_comp_5, snr=5, rep=500, I=50)  # L = 5

# M = 2048 ================================
M <- 2048
alpha_comp_2 <- matrix(c(f_test(M)$bumps, f_test(M)$doppler), ncol=2)
alpha_comp_4 <- matrix(c(f_test(M)$bumps, f_test(M)$doppler, f_test(M)$blocks,
                         f_test(M)$heavisine), ncol=4)
alpha_comp_5 <- matrix(c(f_test(M)$bumps, f_test(M)$doppler, f_test(M)$blocks,
                         f_test(M)$heavisine, f_test(M)$logit), ncol=5)

# simulacao
simu_5_2048_2 <- simu(alpha_comp_2, snr=5, rep=500, I=50)  # L = 2
simu_5_2048_4 <- simu(alpha_comp_4, snr=5, rep=500, I=50)  # L = 4
simu_5_2048_5 <- simu(alpha_comp_5, snr=5, rep=500, I=50)  # L = 5

