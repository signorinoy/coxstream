#' Extract the coefficients from a \code{coxstream} object
#'
#' @param object An object of class \code{coxstream}.
#' @param ... Additional arguments (not used).
#'
#' @return A numeric vector of coefficients.
#' @export
coef.coxstream <- function(object, ...) {
  if (!inherits(object, "coxstream")) {
    stop("The object must be of class 'coxstream'.")
  }
  return(object$theta_prev)
}
