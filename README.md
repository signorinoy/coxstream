
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
formula <- Surv(time, status) ~ X1 + X2 + X3 + X4 + X5
fit <- coxstream(
  formula, sim[sim$batch_id == 1, ],
  n_basis = 3, boundary = c(0, 3), idx_col = "patient_id"
)
for (batch in 2:10) {
  fit <- update(fit, sim[sim$batch_id == batch, ])
}
summary(fit)
#> Call:
#> coxstream(formula = formula, data = sim[sim$batch_id == 1, ], 
#>     n_basis = 3, boundary = c(0, 3), idx_col = "patient_id")
#> 
#> Number of basis:  7 
#>       coef exp(coef)      se     z      p    
#> X1 0.98898   2.68848 0.05123 19.30 <2e-16 ***
#> X2 1.01042   2.74675 0.05364 18.84 <2e-16 ***
#> X3 1.01112   2.74869 0.05431 18.62 <2e-16 ***
#> X4 0.89296   2.44234 0.05207 17.15 <2e-16 ***
#> X5 1.03635   2.81892 0.05106 20.30 <2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#>    exp(coef) exp(-coef) lower .95 upper .95
#> X1 2.6885    0.3720     2.4317    2.9724   
#> X2 2.7467    0.3641     2.4727    3.0512   
#> X3 2.7487    0.3638     2.4711    3.0574   
#> X4 2.4423    0.4094     2.2054    2.7048   
#> X5 2.8189    0.3547     2.5505    3.1156
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
