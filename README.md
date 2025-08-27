
# fda.simu

O pacote `fda.simu` permite a utilização de técnincas da *Análise de
Dados Funcionais Agregados*. No geral, neste tipo de problema tem-se o
interesse em estimar curvas componentes $`\alpha_1(t)`$, , $\alpha_L(t)$
a partir de uma amostra agregada, isto é, cada elemento da amostra é uma
combinação linear das $L$ funções componentes. O modelo para esse tipo
de problema é dado por:

$$
A(t) = \sum_{l=1}^L y_l\alpha_l(t) + e(t)
$$

## Installation

A instalação do pacote é feita através de:

``` r
# install.packages("devtools")
devtools::install_github("jsicas/fda.simu")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(fda.simu)
## basic example code
```

## Referências

Sousa, A. R. S. (2024). A Wavelet-Based Method in Aggregated Functional
Data Analysis. *Monte Carlo Methods and Applications*, 30 (1): 19–30.
DOI: [10.1515/mcma-2023-2016](https://doi.org/10.1515/mcma-2023-2016).
