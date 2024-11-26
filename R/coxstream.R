#' Streaming Cox model with repeated observations
#'
#' `r lifecycle::badge('experimental')`
#' @param formula A formula expression as for regression models, of the form
#' \code{response ~ predictors}. The response must be a survival object as
#' returned by the \code{\link{Surv}} function.
#' @param data A data frame containing the variables in the model.
#' @param degree The degree of the Bernstein polynomial.
#' @param boundary A vector of length 2 containing the boundary knots.
#' @param idx_col The column name of the index column, which is used to
#' distinguish different patients.
#'
#' @return An object of class \code{coxstream}.
#' @export
#'
#' @examples
#' library(coxstream)
#' formula <- Surv(time, status) ~ X1 + X2 + X3 + X4 + X5
#' fit <- coxstream(
#'   formula, sim[sim$batch_id == 1, ],
#'   degree = 6, boundary = c(0, 3), idx_col = "patient_id"
#' )
#' for (batch in 2:10) {
#'   fit <- update(fit, sim[sim$batch_id == batch, ])
#' }
#' summary(fit)
coxstream <- function(formula, data, degree, boundary, idx_col) {
  mf <- stats::model.frame(formula, data)
  y <- stats::model.response(mf)
  time <- y[, 1]
  delta <- y[, 2]
  names(time) <- names(delta) <- data[[idx_col]]
  x <- stats::model.matrix(formula, data)[, -1]
  rownames(x) <- data[[idx_col]]

  # Store the survival times for each patient
  time_stored <- sort(time[delta == 0])
  n_passes <- length(time_stored)

  sorted <- order(time)
  x <- x[sorted, ]
  time <- time[sorted]
  delta <- delta[sorted]

  p <- ncol(x)
  q <- degree + 1
  theta_prev <- numeric(p + q)
  hess_prev <- matrix(0, p + q, p + q)
  time_int <- c()

  res <- stats::optim(
    par = theta_prev, fn = fn, gr = gr, method = "BFGS",
    x = x, time = time, delta = delta, degree = degree, boundary = boundary,
    theta_prev = theta_prev, hess_prev = hess_prev, time_int = time_int
  )
  coef <- res$par
  hess <- hess(
    coef,
    x = x, time = time, delta = delta, degree = degree, boundary = boundary,
    theta_prev = theta_prev, hess_prev = hess_prev, time_int = time_int
  )

  coef_names <- c(paste0("Basis ", 1:(degree + 1)), colnames(x))
  names(coef) <- coef_names
  colnames(hess) <- rownames(hess) <- coef_names

  fit <- list(
    degree = degree,
    boundary = boundary,
    theta_prev = coef,
    hess_prev = hess,
    time_stored = time_stored,
    n_passes = n_passes,
    formula = formula,
    idx_col = idx_col,
    call = match.call()
  )
  class(fit) <- "coxstream"
  return(fit)
}
