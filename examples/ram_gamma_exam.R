# require(wavethresh)
set.seed(123)

# definindo parâmetros
M <- 16                      # quantidade de pontos por função
I <- 2                       # quantidade de observações
beta <- 196/25; lambda <- 2  # parâmetros do erro Gamma(beta, lambda)

# funções componentes: Bumps e Doppler
alpha <- matrix(c(f_test(M)$bumps, f_test(M)$heavisine), ncol = 2)
par(mfrow = c(1,2))
plot(alpha[,1], type = 'b', ylab='', main = 'Bumps')
plot(alpha[,2], type = 'b', ylab='', main = 'Heavisine')
par(mfrow = c(1,1))

# pesos
y1 <- runif(I)                               # pesos da curva 1 (bumps)
y <- matrix(c(y1, 1 - y1), nrow=2, byrow=T)  # matriz de pesos

# amostra
f <- alpha %*% y  # combinação das funções verdadeiras
e <- matrix(rgamma(M * I, shape = beta, rate = lambda), nrow = M,
            ncol = I, byrow = T)
A <- f + e

# DWT
D <- apply(A, MARGIN = 2, \(x)  # mesmo que W %*% A
           wd(x, filter.number = 5, family = 'DaubExPhase'))

# executando algorítmo RAM
tau <- 5; alpha_priori <- 0.8  # parâmetros da priori

## gerando chutes iniciais para theta
theta_1 <- matrix(0, I, M)
for (i in 1:I) {
  for (j in 1:100) {  # tentando 100 vezes para cada coluna de D
    theta_1[i,] <- gera_ponto(a = NULL, D[[i]], lim_sup = 15, filter.number = 5,
                              family = 'DaubExPhase')
    if (post_gamma(theta_1[i,], D[[i]], beta, tau, lambda, alpha_priori) != 0) break
    if (j == 100) stop('Nenhum ponto inicial válido para theta_1 encontrado.')
  }
}

delta_D <- matrix(0, nrow = M, ncol = I)
for (i in 1:I) {
  d_i <- ram_gamma(theta_1[i,], S_1 = NULL, D[[i]], n_ite = 3000,
                   alpha_priori, tau, beta, lambda, gamma = 2/3,
                   filter.number = 5, family = 'DaubExPhase')
  delta_D[,i] <- colMeans(d_i$theta[seq(50, 3000, 2),])  # burn-in e thining
}

# estimando Gamma
Gamma_hat <- delta_D %*% t(y) %*% solve(y %*% t(y))  # delta(D)y'(yy')^-1

# IDWT
alpha_hat <- GenW(M, filter.number = 5, family = 'DaubExPhase') %*% Gamma_hat

# resultados
knitr::kable(data.frame('theta_i' = 1:M, 'bumps' = alpha[,1],
                        'bumps_est' = alpha_hat[,1], 'doppler' = alpha[,2],
                        'doppler_est' = alpha_hat[,2]))

par(mfrow=c(1,2))
plot(alpha[,1], type = 'b', main = 'Bumps', ylab = '',
     ylim = c(min(alpha[,1]), max(alpha_hat[,1])))
lines(alpha_hat[,1], type = 'b', col = 'blue')
legend('topright', bty = 'n', lwd = 1, cex = 0.6,
       legend = c('Curva real', 'Cruva recuperada'), col = c('black', 'blue'))
plot(alpha[,2], type = 'b', main = 'Bumps',
     ylim = c(min(alpha[,2]), max(alpha_hat[,2])))
lines(alpha_hat[,2], type = 'b', col = 'blue')
legend('bottomright', bty = 'n', lwd = 1, cex = 0.6,
       legend=c('Curva real', 'Cruva recuperada'), col = c('black', 'blue'))
