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
#' @param verbose A logical indicating whether to print the optimization
#' messages. Default is \code{FALSE}.
#'
#' @return An object of class \code{coxstream}.
#' @export
#'
#' @examples
#' library(coxstream)
#' formula <- survival::Surv(time, status) ~ X1 + X2 + X3 + X4 + X5
#' fit <- coxstream(
#'   formula, sim[sim$batch_id == 1, ],
#'   degree = 3, boundary = c(0, 3), idx_col = "patient_id"
#' )
#' for (batch in 2:10) {
#'   fit <- update(fit, sim[sim$batch_id == batch, ])
#' }
#' summary(fit)
coxstream <- function(
    formula, data, degree, boundary, idx_col, verbose = FALSE) {
  mf <- stats::model.frame(formula, data)
  y <- stats::model.response(mf)
  time <- y[, 1]
  delta <- y[, 2]
  names(time) <- names(delta) <- data[[idx_col]]
  x <- stats::model.matrix(formula, data)[, -1]
  rownames(x) <- data[[idx_col]]

  # Store the survival times for each patient
  time_stored <- sort(time[delta == 0])
  time_unique <- unique(time)

  sorted <- order(time)
  x <- x[sorted, ]
  time <- time[sorted]
  delta <- delta[sorted]

  p <- ncol(x)
  q <- degree + 1
  theta_prev <- numeric(p + q)
  hess_prev <- matrix(0, p + q, p + q)
  time_int <- c()

  res <- stats::nlm(
    f = objective, p = theta_prev, hessian = TRUE,
    x = x, time = time, delta = delta, degree = degree, boundary = boundary,
    theta_prev = theta_prev, hess_prev = hess_prev, time_int = time_int,
    print.level = ifelse(verbose, 2, 0)
  )
  if (res$code > 3) {
    stop("The optimization did not converge.")
  }
  coef <- res$estimate
  hess <- res$hessian

  coef_names <- c(paste0("Basis ", 1:(degree + 1)), colnames(x))
  names(coef) <- coef_names
  colnames(hess) <- rownames(hess) <- coef_names

  fit <- list(
    degree = degree,
    boundary = boundary,
    theta_prev = coef,
    hess_prev = hess,
    time_stored = time_stored,
    time_unique = time_unique,
    formula = formula,
    idx_col = idx_col,
    verbose = verbose,
    call = match.call()
  )
  class(fit) <- "coxstream"
  return(fit)
}
