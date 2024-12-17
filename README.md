
<!-- README.md is generated from README.Rmd. Please edit that file -->

# CoxStream

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
pak::pak("SignorinoY/coxstream")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(coxstream)
## basic example code
formula <- survival::Surv(time, status) ~ X1 + X2 + X3 + X4 + X5
fit <- coxstream(
  formula, sim[sim$batch_id == 1, ],
  degree = 1, boundary = c(0, 3), idx_col = "patient_id"
)
for (batch in 2:10) {
  fit <- update(fit, sim[sim$batch_id == batch, ])
}
summary(fit)
#> Call:
#> coxstream(formula = formula, data = sim[sim$batch_id == 1, ], 
#>     degree = 1, boundary = c(0, 3), idx_col = "patient_id")
#> 
#>       coef exp(coef)      se     z      p    
#> X1 0.96244   2.61808 0.04949 19.45 <2e-16 ***
#> X2 0.98791   2.68562 0.05140 19.22 <2e-16 ***
#> X3 0.99161   2.69558 0.05260 18.85 <2e-16 ***
#> X4 0.88978   2.43459 0.05029 17.69 <2e-16 ***
#> X5 0.99063   2.69292 0.04812 20.59 <2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#>    exp(coef) exp(-coef) lower .95 upper .95
#> X1 2.6181    0.3820     2.3760    2.8848   
#> X2 2.6856    0.3724     2.4283    2.9703   
#> X3 2.6956    0.3710     2.4316    2.9883   
#> X4 2.4346    0.4107     2.2061    2.6868   
#> X5 2.6929    0.3713     2.4506    2.9592
```

Besides the estimated coefficients, we can also obtain the estimated
baseline hazard function. The following code plots the estimated
cumulative baseline hazard function and the true cumulative baseline
hazard function.

``` r
time <- seq(0, 3, length.out = 100)
basehaz_pred <- basehaz(fit, time)
basehaz_true <- cbind(time, 0.5 * time^2)
plot(
  time, basehaz_pred[, 2],
  type = "l", lty = 2,
  ylim = range(basehaz_pred[, 4], basehaz_pred[, 5]),
  xlab = expression(t), ylab = expression(Lambda[0](t)),
  main = "Cumulative Baseline Hazard (Estimate vs True)"
)
polygon(
  c(time, rev(time)), c(basehaz_pred[, 4], rev(basehaz_pred[, 5])),
  col = rgb(0.5, 0.5, 0.5, 0.4), border = NA
)
lines(basehaz_true[, 1], basehaz_true[, 2])
legend("topleft", legend = c("Estimate", "True"), lty = c(2, 1))
```

<img src="man/figures/README-evaluation-1.png" alt="Cumulative Baseline Hazard (Estimate vs True)" width="100%" />
