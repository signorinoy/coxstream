#' Prediction method for \code{coxtrans} objects.
#' @param object An object of class \code{coxtrans}.
#' @param newdata Optional new data for making predictions. If omitted,
#'   predictions are made using the data used for fitting the model.
#' @param type The type of prediction to perform. Options include:
#'   \describe{
#'     \item{\code{"lp"}}{The linear predictor.}
#'     \item{\code{"terms"}}{The components of the linear predictor.}
#'     \item{\code{"risk"}}{The risk score \eqn{\exp(\text{lp})}.}
#'     \item{\code{"expected"}}{The expected number of events, given the
#'              covariates and follow-up time.}
#'     \item{\code{"survival"}}{The survival probability, given the covariates
#'             and follow-up time.}
#'   }
#' @param ... Additional arguments (not unused).
#' @return A numeric vector of predictions.
#' @export
predict.coxstream <- function(
    object, newdata = NULL,
    type = c("lp", "terms", "risk", "expected", "survival"), ...) {
  type <- match.arg(type)
  x <- stats::model.matrix(object$formula, newdata)[, -1]
  p <- ncol(x)
  q <- object$degree + 1
  beta <- object$theta_prev[(q + 1):(q + p)]
  lp <- x %*% beta
  if (type == "lp") {
    return(lp)
  } else if (type == "risk") {
    return(exp(lp))
  } else {
    stop("type must be one of 'lp' or 'risk'")
  }
}
