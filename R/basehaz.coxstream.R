#' Predict the baseline hazard function for \code{coxstream} objects
#'
#' @param object An object of class \code{coxstream}.
#' @param newdata A numeric vector of time points at which to predict the
#' baseline hazard function. If \code{NULL}, the function will predict the
#' baseline hazard function at 100 equally spaced time points between the
#' boundaries of the \code{coxstream} object.
#' @param conf.int A numeric value between 0 and 1 indicating the confidence
#' level for the confidence interval. Default is 0.95.
#' @param ... Additional arguments (not unused).
#'
#' @return A matrix with one row for each time point, and columns containing
#' the baseline hazard function, standard error, and confidence interval.
#' @export
basehaz.coxstream <- function(
    object, newdata = NULL, conf.int = 0.95, ...) {
  if (!inherits(object, "coxstream")) {
    stop("object must be of class 'coxstream'")
  }

  degree <- object$degree
  boundary <- object$boundary

  if (!is.null(newdata)) {
    time <- newdata
  } else {
    time <- seq.int(boundary[1], boundary[2], length.out = 100)
  }
  q <- degree + 1
  alpha <- object$theta_prev[1:q]
  u <- unique(time)
  b_pre <- splines2::bpoly(u, degree, TRUE, boundary)
  b <- b_pre[match(time, u), ]
  basehaz <- exp(b %*% alpha)

  vcov_alpha <- vcov(object)[1:q, 1:q]
  vcov_lp <- b %*% vcov_alpha %*% t(b)
  vcov_basehaz <- diag(as.vector(basehaz^2)) %*% vcov_lp
  se_basehaz <- sqrt(diag(vcov_basehaz))
  z <- stats::qnorm((1 + conf.int) / 2)

  basehaz_matrix <- cbind(
    time, basehaz, se_basehaz,
    basehaz - z * se_basehaz, basehaz + z * se_basehaz
  )
  dimnames(basehaz_matrix) <- list(
    NULL, c(
      "time", "hazard", "se(hazard)",
      paste("lower .", round(100 * conf.int, 2), sep = ""),
      paste("upper .", round(100 * conf.int, 2), sep = "")
    )
  )
  return(basehaz_matrix)
}
