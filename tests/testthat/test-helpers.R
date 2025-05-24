# Test prox_forward
test_that("prox_forward returns identity for p == q", {
  expect_equal(prox_forward(3, 3, 4), diag(7))
})

test_that("prox_forward errors for invalid input", {
  expect_error(prox_forward(3, 2, 1))
  expect_error(prox_forward(-1, 2, 1))
  expect_error(prox_forward(2.5, 3, 1))
  expect_error(prox_forward(2, 3, -1))
})

test_that("prox_forward returns correct shape", {
  mat <- prox_forward(2, 4, 1)
  expect_equal(dim(mat), c(5, 3))
})

# Test bernstein
test_that("bernstein returns correct shape and values", {
  x <- seq(0, 1, length.out = 5)
  b <- bernstein(x, 3, c(0, 1))
  expect_equal(dim(b), c(5, 3))
  expect_true(all(b >= 0))
  expect_true(all(abs(rowSums(b) - 1) < 1e-8))
})

test_that("bernstein errors for out-of-bound x", {
  expect_error(bernstein(-0.1, 3, c(0, 1)))
  expect_error(bernstein(1.1, 3, c(0, 1)))
})

# Test basehaz_ode, basehaz_grad, basehaz_hess
test_that("basehaz_ode returns correct value", {
  parms <- list(alpha = c(0.1, 0.2, 0.3), n_basis = 3, boundary = c(0, 1))
  val <- basehaz_ode(0.5, NULL, parms)
  expect_type(val, "list")
  expect_true(val[[1]] > 0)
})

test_that("basehaz_grad returns correct gradient", {
  parms <- list(alpha = c(0.1, 0.2, 0.3), n_basis = 3, boundary = c(0, 1))
  grad <- basehaz_grad(0.5, NULL, parms)
  expect_type(grad, "list")
  expect_equal(length(grad[[1]]), 3)
})

test_that("basehaz_hess returns correct Hessian", {
  parms <- list(alpha = c(0.1, 0.2, 0.3), n_basis = 3, boundary = c(0, 1))
  hess <- basehaz_hess(0.5, NULL, parms)
  expect_type(hess, "list")
  expect_equal(dim(hess[[1]]), c(3, 3))
})

# Test objective: gradient and hessian correctness
test_that("objective gradient and hessian are correct (numDeriv)", {
  set.seed(1)
  n <- 10
  n_features <- 2
  n_basis <- 3
  x <- matrix(rnorm(n * n_features), n, n_features)
  time <- sort(runif(n, 0, 1))
  status <- sample(0:1, n, replace = TRUE)
  cluster <- rep(1:2, length.out = n)
  boundary <- c(0, 1)
  par <- c(runif(n_basis), runif(n_features))
  theta_prev <- par
  hess_prev <- diag(length(par))
  time_int <- NULL

  # function for numDeriv
  fn <- function(par) {
    res <- objective(
      par, time, status, x, cluster, boundary, theta_prev, hess_prev, time_int
    )
    as.numeric(res)
  }
  gr <- function(par) {
    res <- objective(
      par, time, status, x, cluster, boundary, theta_prev, hess_prev, time_int
    )
    attr(res, "gradient")
  }
  hs <- function(par) {
    res <- objective(
      par, time, status, x, cluster, boundary, theta_prev, hess_prev, time_int
    )
    hess <- attr(res, "hessian")
    as.matrix(hess)
  }

  grad_num <- numDeriv::grad(fn, par)
  grad_ana <- gr(par)
  expect_equal(grad_ana, grad_num, tolerance = 1e-4)

  hess_num <- numDeriv::hessian(fn, par)
  hess_ana <- hs(par)
  dimnames(hess_num) <- NULL
  dimnames(hess_ana) <- NULL
  expect_equal(hess_ana, hess_num, tolerance = 1e-3)
})
