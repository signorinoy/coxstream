#' Update the \code{coxstream} with new data.
#'
#' `r lifecycle::badge('experimental')`
#' @param object A \code{coxstream} object.
#' @param data A data frame containing the variables in the model.
#' @param degree An integer representing the degree of the Bernstein polynomial,
#' or a string specifying the degree of the Bernstein polynomial as a function
#' of the number of patients, "auto", which is calculated as
#' \eqn{[\alpha N^{\nu}]}, where \eqn{N} is the number of the unique survival
#' times.
#' @param alpha A numeric value used to calculate the degree of the Bernstein
#' polynomial. Default is 2.
#' @param nu A numeric value used to calculate the degree of the Bernstein
#' polynomial. Default is 0.2.
#' @param ... Additional arguments (not used).
#'
#' @return An object of class \code{coxstream}.
#' @export
update.coxstream <- function(
    object, data, degree = "auto", alpha = 2, nu = 0.2, ...) {
  formula <- object$formula
  boundary <- object$boundary
  idx_col <- object$idx_col
  time_stored <- object$time_stored

  mf <- stats::model.frame(formula, data)
  y <- stats::model.response(mf)
  time <- y[, 1]
  delta <- y[, 2]
  names(time) <- names(delta) <- data[[idx_col]]
  x <- stats::model.matrix(formula, data)[, -1]
  rownames(x) <- data[[idx_col]]

  time_unique <- unique(c(object$time_unique, time))
  n_features <- length(object$theta_prev) - object$degree - 1

  if (degree == "auto") {
    degree <- max(
      object$degree, round(alpha * length(time_unique)^nu)
    )
    degree <- min(degree, object$degree + 1)
  } else if (is.numeric(degree)) {
    degree <- as.integer(degree)
  } else {
    stop("The degree must be an integer or 'auto'.")
  }

  while (degree > object$degree) {
    object$degree <- object$degree + 1
    prox <- matrix(0, object$degree + 1, object$degree)
    diag(prox) <- (object$degree - 0:(object$degree - 1)) / object$degree
    diag(prox[-1, ]) <- (1:object$degree) / object$degree
    prox <- rbind(
      cbind(prox, matrix(0, object$degree + 1, n_features)),
      cbind(matrix(0, n_features, object$degree), diag(n_features))
    )
    object$theta_prev <- as.vector(prox %*% object$theta_prev)
    object$hess_prev <- prox %*% object$hess_prev %*% t(prox)
  }

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

  res <- trust::trust(
    objfun = objective, parinit = theta_prev, rinit = 1, rmax = 10,
    x = x, time = time, delta = delta, degree = degree, boundary = boundary,
    theta_prev = theta_prev, hess_prev = hess_prev, time_int = time_int
  )
  if (!res$converged) {
    warning("The optimization did not converge.")
  }
  object$theta_prev <- res$argument
  object$hess_prev <- res$hessian

  coef_names <- c(paste0("Basis ", 1:(object$degree + 1)), colnames(x))
  names(object$theta_prev) <- coef_names
  colnames(object$hess_prev) <- rownames(object$hess_prev) <- coef_names

  object$time_stored <- time_stored
  object$time_unique <- time_unique
  class(object) <- "coxstream"
  return(object)
}
