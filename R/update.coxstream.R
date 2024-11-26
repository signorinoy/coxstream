#' Update the \code{coxstream} with new data.
#'
#' `r lifecycle::badge('experimental')`
#' @param object A \code{coxstream} object.
#' @param data A data frame containing the variables in the model.
#' @param ... Additional arguments (not used).
#'
#' @return An object of class \code{coxstream}.
#' @export
update.coxstream <- function(object, data, ...) {
  idx_col <- object$idx_col
  formula <- object$formula
  time_stored <- object$time_stored

  mf <- stats::model.frame(formula, data)
  y <- stats::model.response(mf)
  time <- y[, 1]
  delta <- y[, 2]
  names(time) <- names(delta) <- data[[idx_col]]
  x <- stats::model.matrix(formula, data)[, -1]
  rownames(x) <- data[[idx_col]]

  # The patients that are already stored and the new patients
  patients_stored <- names(time_stored)
  patients_int <- intersect(patients_stored, data[[idx_col]])
  patients_new <- data[[idx_col]][!data[[idx_col]] %in% patients_int]

  # Select the stored survival times that int with the new patients
  time_int <- time_stored[patients_int]

  # Update the stored survival times
  time_stored[patients_int] <- time[patients_int]
  # Select the survival patients in the inted stored survival times
  delta_int <- delta[patients_int]
  patients_remove <- names(delta_int[delta_int == 1])
  time_stored <- time_stored[!names(time_stored) %in% patients_remove]
  time_new <- time[data[[idx_col]] %in% patients_new & delta == 0]
  time_stored <- c(time_stored, time_new)
  time_stored <- sort(time_stored)

  sorted <- order(time)
  x <- x[sorted, ]
  time <- time[sorted]
  delta <- delta[sorted]

  theta_prev <- object$theta_prev
  hess_prev <- object$hess_prev

  sr <- stats::optim(
    par = theta_prev, fn = fn, gr = gr, method = "BFGS",
    x = x, time = time, delta = delta, degree = object$degree,
    boundary = object$boundary, theta_prev = theta_prev, hess_prev = hess_prev,
    time_int = time_int
  )

  object$theta_prev <- sr$par
  object$hess_prev <- hess(
    sr$par,
    x = x, time = time, delta = delta, degree = object$degree,
    boundary = object$boundary, theta_prev = theta_prev, hess_prev = hess_prev,
    time_int = time_int
  )
  object$time_stored <- time_stored
  object$n_passes <- object$n_passes + length(patients_new)
  class(object) <- "coxstream"
  return(object)
}
