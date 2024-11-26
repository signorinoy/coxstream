
<!-- README.md is generated from README.Rmd. Please edit that file -->

# coxstream

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![CRAN
status](https://www.r-pkg.org/badges/version/coxstream)](https://CRAN.R-project.org/package=coxstream)
[![R-CMD-check](https://github.com/SignorinoY/coxstream/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/SignorinoY/coxstream/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/SignorinoY/coxstream/graph/badge.svg)](https://app.codecov.io/gh/SignorinoY/coxstream)
<!-- badges: end -->

The goal of coxstream is to …

## Installation

You can install the development version of coxstream like so:

``` r
# install.packages("pak")
pak::pak("SignorinoY/CoxStream")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(coxstream)
#> Loading required package: survival
## basic example code
formula <- Surv(time, status) ~ X1 + X2 + X3 + X4 + X5
fit <- coxstream(
  formula, sim[sim$batch_id == 1, ],
  degree = 4, boundary = c(0, 3), idx_col = "patient_id"
)
for (batch in 2:10) {
  fit <- update(fit, sim[sim$batch_id == batch, ])
}
summary(fit)
#> Call:
#> coxstream(formula = formula, data = sim[sim$batch_id == 1, ], 
#>     degree = 4, boundary = c(0, 3), idx_col = "patient_id")
#> 
#>       coef exp(coef)      se     z      p    
#> X1 0.96860   2.63426 0.05107 18.97 <2e-16 ***
#> X2 0.99686   2.70976 0.05355 18.61 <2e-16 ***
#> X3 0.99161   2.69557 0.05368 18.47 <2e-16 ***
#> X4 0.88254   2.41704 0.05205 16.96 <2e-16 ***
#> X5 1.01817   2.76811 0.05092 19.99 <2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#>    exp(coef) exp(-coef) lower .95 upper .95
#> X1 2.6343    0.3796     2.3833    2.9116   
#> X2 2.7098    0.3690     2.4398    3.0096   
#> X3 2.6956    0.3710     2.4264    2.9946   
#> X4 2.4170    0.4137     2.1826    2.6766   
#> X5 2.7681    0.3613     2.5052    3.0587
```
