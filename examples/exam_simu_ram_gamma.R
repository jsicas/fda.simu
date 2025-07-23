I <- 5                  # número de observações
M <- 32                 # número de pontos por observação
beta <- 2; lambda <- 2  # erro gamma
alpha_comp <- matrix(c(f_test(M)$bumps, f_test(M)$doppler, f_test(M)$heavi), M)
L <- ncol(alpha_comp)   # quantidade de curvas componentes
simu_ram_gamma(alpha_comp, I, n_ite=300, rep=3, alpha=0.7, tau=3,
               beta=beta,lambda=lambda)
