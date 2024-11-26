#' Extract variance-covariance matrix from a \code{coxstream} bbject
#'
#' @param object An object of class \code{coxstream}.
#' @param ... Additional arguments (not unused).
#'
#' @return A matrix representing the variance-covariance matrix of the
#' coefficients.
#' @export
vcov.coxstream <- function(object, ...) {
  if (!inherits(object, "coxstream")) {
    stop("The provided object must be of class 'coxstream'.")
  }
  return(solve(object$hess_prev))
}
