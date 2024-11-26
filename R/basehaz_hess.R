#' Hessian Matrix of the Baseline Hazard Function
#'
#' @param t A non-negative real number specifying the time point at which to
#' evaluate the Hessian matrix of the baseline hazard function.
#' @param y The initial (state) values of the ODE system, which are not used in
#' Cox proportional hazards models.
#' @param parms A list of parameters containing the following components:
#' \describe{
#'   \item{alpha}{The coefficients of the Bernstein polynomial.}
#'   \item{degree}{A non-negative integer specifying the degree of the
#' Bernstein polynomial.}
#'   \item{boundary}{Boundary knots at which to anchor the Bernstein
#' polynomial.}
#' }
#'
#' @return The Hessian matrix with respect to the coefficients of the baseline
#' hazard function evaluated at the specified time point.
basehaz_hess <- function(t, y, parms) {
  b <- splines2::bpoly(t, parms$degree, TRUE, parms$boundary)
  basehaz_hess <- t(b) %*% b * exp(sum(parms$alpha * b))
  return(list(basehaz_hess))
}
