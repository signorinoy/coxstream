#' Summary method for \code{coxstream} object
#'
#' @param object An object of class \code{coxstream}.
#' @param conf.int A numeric value between 0 and 1 indicating the confidence
#' level for the confidence interval. Default is 0.95.
#' @param ... Additional arguments (not unused).
#'
#' @return An object of class \code{"summary.coxstream"}, containing the following
#' components:
#' \describe{
#'   \item{call}{The matched call.}
#'   \item{coefficients}{A matrix with one row for each coefficient, and
#'   columns containing the coefficient, the hazard ratio (exp(coef)), standard
#'   error, Wald statistic, and P value.}
#'   \item{conf.int}{A matrix with one row for each coefficient, containing the
#'   confidence limits for exp(coef).}
#' }
#' @export
summary.coxstream <- function(object, conf.int = 0.95, ...) {
  if (!inherits(object, "coxstream")) {
    stop("object must be of class 'coxstream'")
  }
  coef <- coef(object)
  vcov <- vcov(object)

  se <- sqrt(diag(vcov))
  z_scores <- coef / se
  p_values <- stats::pchisq(z_scores^2, df = 1, lower.tail = FALSE)
  coef_matrix <- cbind(coef, exp(coef), se, z_scores, p_values)
  dimnames(coef_matrix) <- list(
    names(coef), c("coef", "exp(coef)", "se", "z", "p")
  )

  z <- stats::qnorm((1 + conf.int) / 2)
  conf_int_matrix <- cbind(
    exp(coef), exp(-coef), exp(coef - z * se), exp(coef + z * se)
  )
  dimnames(conf_int_matrix) <- list(
    names(coef), c(
      "exp(coef)", "exp(-coef)",
      paste("lower .", round(100 * conf.int, 2), sep = ""),
      paste("upper .", round(100 * conf.int, 2), sep = "")
    )
  )

  summary_list <- list(
    call = object$call, coefficients = coef_matrix, conf.int = conf_int_matrix,
    degree = object$degree, boundary = object$boundary
  )
  class(summary_list) <- "summary.coxstream"
  return(summary_list)
}
