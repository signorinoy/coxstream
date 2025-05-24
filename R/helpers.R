#' Proximal Operator (Forward)
#'
#' Constructs a transformation matrix to map the coefficients of Bernstein
#' polynomial bases from one dimension to another, ensuring
#' consistency with the properties of Bernstein polynomial bases. Note that
#' the original dimension must be less than or equal to the target dimension.
#'
#' @param p Integer. The number of basis functions in the original dimension.
#' @param q Integer. The number of basis functions in the target dimension.
#' @param d Integer. The number of features.
#'
#' @return A matrix representing the transformation from the original dimension
#'         to the target dimension.
#'
#' @details This function iteratively constructs the transformation matrix by
#'          expanding the basis functions from \code{p} to \code{q}, while
#'          preserving the mathematical properties of Bernstein polynomials.
prox_forward <- function(p, q, d) {
  if (
    !is.numeric(p) || !is.numeric(q) || !is.numeric(d) ||
      p != as.integer(p) || q != as.integer(q) || d != as.integer(d)
  ) {
    stop("p, q, and d must be integers")
  }
  if (p <= 0 || q <= 0 || d <= 0) stop("p, q, and d must be positive integers")
  if (p > q) stop("p must be less than or equal to q")
  if (p == q) {
    return(diag(p + d))
  }
  prox_forward <- diag(p)
  for (i in p:(q - 1)) {
    if (i == 1) {
      prox <- matrix(1, 2, 1)
    } else {
      prox <- matrix(0, i + 1, i)
      diag(prox) <- (i - 0:(i - 1)) / i
      diag(prox[-1, ]) <- (1:i) / i
    }
    prox_forward <- prox %*% prox_forward
  }
  prox_forward <- rbind(
    cbind(prox_forward, matrix(0, q, d)),
    cbind(matrix(0, d, p), diag(d))
  )
  prox_forward
}

bernstein <- function(x, n_basis, boundary) {
  if (n_basis < 0 && n_basis != as.integer(n_basis)) {
    stop("n_basis must be a non-negative integer")
  }
  if (length(boundary) != 2) stop("boundary must be a vector of length 2")
  if (any(x < boundary[1]) || any(x > boundary[2])) {
    stop("x must be within the range defined by boundary")
  }
  b <- if (n_basis == 1) {
    matrix(1, nrow = length(x), ncol = 1)
  } else {
    splines2::bpoly(
      x, n_basis - 1,
      intercept = TRUE, Boundary.knots = boundary
    )
  }
  b
}

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
#'   \item{alpha}{The coefficients of the Bernstein polynomial bases.}
#'   \item{n_basis}{A positive odd integer specifying the number of Bernstein
#' polynomial bases.}
#'  \item{boundary}{A numeric vector of length 2 that defines the lower and
#' upper bounds for normalization.}
#' }
#'
#' @return The baseline hazard function evaluated at the specified time point.
basehaz_ode <- function(t, y, parms) {
  b <- bernstein(t, parms$n_basis, parms$boundary)
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
#'   \item{alpha}{The coefficients of the bernstein polynomial.}
#'   \item{n_basis}{A positive odd integer specifying the number of bernstein
#' basis functions.}
#'  \item{boundary}{A numeric vector of length 2 that defines the lower and
#' upper bounds for normalization.}
#' }
#'
#' @return The gradient with respect to the coefficients of the baseline hazard
#' function evaluated at the specified time point.
basehaz_grad <- function(t, y, parms) {
  b <- bernstein(t, parms$n_basis, parms$boundary)
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
#'   \item{alpha}{The coefficients of the bernstein polynomial.}
#'   \item{n_basis}{A positive odd integer specifying the number of bernstein
#' basis functions.}
#'  \item{boundary}{A numeric vector of length 2 that defines the lower and
#' upper bounds for normalization.}
#' }
#'
#' @return The Hessian matrix with respect to the coefficients of the baseline
#' hazard function evaluated at the specified time point.
basehaz_hess <- function(t, y, parms) {
  b <- bernstein(t, parms$n_basis, parms$boundary)
  basehaz_hess <- t(b) %*% b * exp(sum(parms$alpha * b))
  list(basehaz_hess)
}

#' Objective Function for Optimization
#'
#' Computes the objective value, gradient, and Hessian for the penalized
#' likelihood in the Cox model with a Bernstein polynomial bases baseline
#' hazard. This function is designed for use in optimization routines.
#'
#' @param par Numeric vector. Parameters to be optimized, including the
#'   Bernstein polynomial coefficients of dimension \code{n_basis} (alpha)
#'   and regression coefficients (beta).
#' @param time Numeric vector. Observed event or censoring times.
#' @param status Numeric vector. Event indicators (1 = event, 0 = censored).
#' @param x Numeric matrix. Covariate data for each observation.
#' @param subject Vector. Subject identifiers for repeated measurements.
#' @param boundary Numeric vector of length 2. The lower and upper bounds for
#'  the Bernstein polynomial basis functions.
#' @param theta_prev Numeric vector. Previous parameter estimates used in the
#'   proximal term, including the Bernstein polynomial coefficients of dimension
#'   \code{n_basis_prev} (alpha) and regression coefficients (beta). The number
#'   of basis functions in the previous estimate must be greater than or equal
#'   to the number of basis functions in the current estimate.
#' @param hess_prev Numeric matrix. Previous Hessian matrix used in the proximal
#'   term, which is a square matrix of size \code{length(theta_prev)}.
#' @param time_int Numeric vector. Additional time points stored in a named
#'   vector, where the names correspond to the subject identifiers. These time
#'   points are used for evaluating the baseline hazard function and its
#'   derivatives.
#'
#' @return The objective value with attributes "gradient" and "hessian".
#'
#' @importFrom utils head
#' @importFrom stats ave
objective <- function(
    par, time, status, x, subject, boundary, theta_prev, hess_prev, time_int) {
  n_parameters <- length(par)
  n_features <- ncol(x)
  n_basis <- n_parameters - n_features

  n_parameters_pre <- length(theta_prev)
  n_basis_pre <- n_parameters_pre - n_features

  # Extract the parameters from the input vector
  alpha <- par[seq_len(n_basis)]
  beta <- par[(n_basis + 1):n_parameters]
  parms <- list(
    alpha = alpha, n_basis = n_basis, boundary = boundary
  )

  prox <- prox_forward(n_basis, n_basis_pre, n_features)
  theta_prev <- qr.solve(prox, theta_prev)
  hess1 <- t(prox) %*% hess_prev %*% prox
  loss1 <- as.numeric((par - theta_prev) %*% hess1 %*% (par - theta_prev) / 2)
  grad1 <- as.vector(hess1 %*% (par - theta_prev))

  time <- c(time, time_int)
  subject <- c(subject, names(time_int))
  order_idx <- order(subject, time)
  time <- time[order_idx]
  subject <- subject[order_idx]

  idx_lag <- ave(
    seq_along(subject), subject,
    FUN = function(x) c(NA, head(x, -1))
  )

  u <- sort(unique(time))
  b_uniq <- bernstein(u, n_basis, boundary)
  b <- b_uniq[match(time, u), ]
  cbh_uniq <- as.matrix(deSolve::ode(
    y = 0, times = c(0, u), func = basehaz_ode, parms = parms, method = "ode45"
  ))[-1, -1, drop = FALSE]
  dcbh_uniq <- as.matrix(deSolve::ode(
    y = rep(0, n_basis), times = c(0, u), func = basehaz_grad,
    parms = parms, method = "ode45"
  ))[-1, -1, drop = FALSE]
  d2cbh_uniq <- as.matrix(deSolve::ode(
    y = rep(0, n_basis**2), times = c(0, u), func = basehaz_hess,
    parms = parms, method = "ode45"
  ))[-1, -1, drop = FALSE]

  cbh <- cbh_uniq[match(time, u), , drop = FALSE]
  cbh <- cbh - ifelse(is.na(idx_lag), 0, cbh[idx_lag, , drop = FALSE])
  dcbh <- dcbh_uniq[match(time, u), , drop = FALSE]
  dcbh[!is.na(idx_lag), ] <- dcbh[!is.na(idx_lag), ] -
    dcbh[idx_lag[!is.na(idx_lag)], , drop = FALSE]
  d2cbh <- d2cbh_uniq[match(time, u), , drop = FALSE]
  d2cbh[!is.na(idx_lag), ] <- d2cbh[!is.na(idx_lag), ] -
    d2cbh[idx_lag[!is.na(idx_lag)], , drop = FALSE]

  restore_idx <- order(order_idx)[seq_len(nrow(x))]
  b <- b[restore_idx, , drop = FALSE]
  cbh <- cbh[restore_idx, , drop = FALSE]
  dcbh <- dcbh[restore_idx, , drop = FALSE]
  d2cbh <- d2cbh[restore_idx, , drop = FALSE]

  eta <- x %*% beta
  dalpha <- t(dcbh) %*% exp(eta) - t(b) %*% status
  dbeta <- t(x) %*% (cbh * exp(eta) - status)
  d2alpha <- matrix(colSums(sweep(d2cbh, 1, exp(eta), "*")), n_basis)
  d2beta <- t(x) %*% (as.vector(cbh * exp(eta)) * x)
  dalph_dabeta <- t(dcbh) %*% (as.vector(exp(eta)) * x)

  loss2 <- t(cbh) %*% exp(eta) - t(b %*% alpha + eta) %*% status
  grad2 <- c(dalpha, dbeta)
  hess2 <- cbind(
    rbind(d2alpha, t(dalph_dabeta)),
    rbind(dalph_dabeta, d2beta)
  )

  loss <- loss1 + loss2
  attr(loss, "gradient") <- grad1 + grad2
  attr(loss, "hessian") <- hess1 + hess2
  loss
}
