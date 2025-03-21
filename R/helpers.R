#' Survival Function Wrapper
#' 
#' This alias provides direct access to the survival::Surv function.
#' @inheritParams survival::Surv
#' @export
Surv <- survival::Surv # nolint: object_name_linter.

#' Return odd number
#'
#' If the number is even, adds 1 to convert it into an odd number.
#'
#' @param x A numeric vector.
#' @return A numeric vector with all elements as odd numbers.
make_odd <- function(x) x + (x %% 2 == 0)

#' Generic function for basehaz
#'
#' @param object Any object.
#' @param ... Additional arguments.
#'
#' @return A numeric vector of baseline hazard.
#' @export
basehaz <- function(object, ...) {
  UseMethod("basehaz")
}

#' Baseline Hazard Function for ODE Evaluation
#'
#' @param t A non-negative real number specifying the time point at which to
#' evaluate the baseline hazard function.
#' @param y The initial (state) values of the ODE system, which are not used in
#' Cox proportional hazards models.
#' @param parms A list of parameters containing the following components:
#' \describe{
#'   \item{alpha}{The coefficients of the Bernstein polynomial.}
#'   \item{n_basis}{A positive odd integer specifying the number of fourier
#' basis functions.}
#'  \item{period}{The period of the Fourier basis functions.}
#' }
#'
#' @return The baseline hazard function evaluated at the specified time point.
basehaz_ode <- function(t, y, parms) {
  b <- fda::fourier(t, nbasis = parms$n_basis, period = parms$period)
  basehaz <- exp(sum(parms$alpha * b))
  list(basehaz)
}

#' Gradient of the Baseline Hazard Function
#'
#' @param t A non-negative real number specifying the time point at which to
#' evaluate the gradient of the baseline hazard function.
#' @param y The initial (state) values of the ODE system, which are not used in
#' Cox proportional hazards models.
#' @param parms A list of parameters containing the following components:
#' \describe{
#'   \item{alpha}{The coefficients of the Bernstein polynomial.}
#'   \item{n_basis}{A positive odd integer specifying the number of fourier
#' basis functions.}
#'  \item{period}{The period of the Fourier basis functions.}
#' }
#'
#' @return The gradient with respect to the coefficients of the baseline hazard
#' function evaluated at the specified time point.
basehaz_grad <- function(t, y, parms) {
  b <- fda::fourier(t, nbasis = parms$n_basis, period = parms$period)
  list(exp(sum(parms$alpha * b)) * b)
}

#' Hessian Matrix of the Baseline Hazard Function
#'
#' @param t A non-negative real number specifying the time point at which to
#' evaluate the Hessian matrix of the baseline hazard function.
#' @param y The initial (state) values of the ODE system, which are not used in
#' Cox proportional hazards models.
#' @param parms A list of parameters containing the following components:
#' \describe{
#'   \item{alpha}{The coefficients of the Bernstein polynomial.}
#'   \item{n_basis}{A positive odd integer specifying the number of fourier
#' basis functions.}
#'  \item{period}{The period of the Fourier basis functions.}
#' }
#'
#' @return The Hessian matrix with respect to the coefficients of the baseline
#' hazard function evaluated at the specified time point.
basehaz_hess <- function(t, y, parms) {
  b <- fda::fourier(t, nbasis = parms$n_basis, period = parms$period)
  basehaz_hess <- t(b) %*% b * exp(sum(parms$alpha * b))
  list(basehaz_hess)
}

objective <- function(
    theta, x, time, delta, n_basis_cur, n_basis_pre, boundary,
    theta_prev, hess_prev, time_int) {

  n_params <- length(theta)
  period <- abs(boundary[2] - boundary[1])

  alpha <- c(theta[seq_len(n_basis_cur)], rep(0, n_basis_pre - n_basis_cur))
  beta <- theta[(n_basis_pre + 1):n_params]
  parms <- list(alpha = alpha, n_basis = n_basis_pre, period = period)

  loss1 <- (theta - theta_prev) %*% hess_prev %*% (theta - theta_prev) / 2
  grad1 <- as.vector(hess_prev %*% (theta - theta_prev))
  hess1 <- hess_prev

  u <- unique(time)
  b_uniq <- fda::fourier(u, n_basis_pre, period)
  b <- b_uniq[match(time, u), ]
  cbh_uniq <- as.matrix(deSolve::ode(
    y = 0, times = c(0, u), func = basehaz_ode, parms = parms, method = "ode45"
  ))[-1, -1, drop = FALSE]
  cbh <- cbh_uniq[match(time, u), , drop = FALSE]
  dcbh_uniq <- as.matrix(deSolve::ode(
    y = rep(0, n_basis_pre), times = c(0, u), func = basehaz_grad,
    parms = parms, method = "ode45"
  ))[-1, -1, drop = FALSE]
  dcbh <- dcbh_uniq[match(time, u), , drop = FALSE]
  d2cbh_uniq <- as.matrix(deSolve::ode(
    y = rep(0, n_basis_pre**2), times = c(0, u), func = basehaz_hess,
    parms = parms, method = "ode45"
  ))[-1, -1, drop = FALSE]
  d2cbh <- d2cbh_uniq[match(time, u), , drop = FALSE]
  eta <- x %*% beta
  dalpha <- t(dcbh) %*% exp(eta) - t(b) %*% delta
  dbeta <- t(x) %*% (cbh * exp(eta) - delta)
  d2alpha <- matrix(colSums(sweep(d2cbh, 1, exp(eta), "*")), n_basis_pre)
  d2beta <- t(x) %*% (as.vector(cbh * exp(eta)) * x)
  dalph_dabeta <- t(dcbh) %*% (as.vector(exp(eta)) * x)

  loss2 <- t(cbh) %*% exp(eta) - t(b %*% alpha + eta) %*% delta
  grad2 <- c(dalpha, dbeta)
  hess2 <- cbind(
    rbind(d2alpha, t(dalph_dabeta)),
    rbind(dalph_dabeta, d2beta)
  )

  if (length(time_int) > 0) {
    u_int <- unique(time_int)
    cbh_int_uniq <- as.matrix(deSolve::ode(
      y = 0, times = c(0, u_int), func = basehaz_ode, parms = parms,
      method = "ode45"
    ))[-1, -1, drop = FALSE]
    cbh_int <- cbh_int_uniq[match(time_int, u_int), , drop = FALSE]
    dcbh_int_uniq <- as.matrix(deSolve::ode(
      y = rep(0, n_basis_pre), times = c(0, u_int), func = basehaz_grad,
      parms = parms, method = "ode45"
    ))[-1, -1, drop = FALSE]
    dcbh_int <- dcbh_int_uniq[match(time_int, u_int), , drop = FALSE]
    d2cbh_int_uniq <- as.matrix(deSolve::ode(
      y = rep(0, n_basis_pre**2), times = c(0, u_int), func = basehaz_hess,
      parms = parms, method = "ode45"
    ))[-1, -1, drop = FALSE]
    d2cbh_int <- d2cbh_int_uniq[match(time_int, u_int), , drop = FALSE]
    x_int <- x[names(time_int), , drop = FALSE]
    eta_int <- x_int %*% beta

    dalpha <- t(dcbh_int) %*% exp(eta_int)
    dbeta <- t(x_int) %*% (cbh_int * exp(eta_int))
    d2alpha_int <- matrix(
      colSums(sweep(d2cbh_int, 1, exp(eta_int), "*")), n_basis_pre
    )
    d2beta_int <- t(x_int) %*% (as.vector(cbh_int * exp(eta_int)) * x_int)
    dalph_dabeta_int <- t(dcbh_int) %*% (as.vector(exp(eta_int)) * x_int)

    loss3 <- -t(cbh_int) %*% exp(eta_int)
    grad3 <- -c(dalpha, dbeta)
    hess3 <- -cbind(
      rbind(d2alpha_int, t(dalph_dabeta_int)),
      rbind(dalph_dabeta_int, d2beta_int)
    )
  } else {
    loss3 <- 0
    grad3 <- 0
    hess3 <- 0
  }

  loss <- loss1 + loss2 + loss3
  grad <- grad1 + grad2 + grad3
  hess <- hess1 + hess2 + hess3
  list(value = loss, gradient = grad, hessian = hess)
}
