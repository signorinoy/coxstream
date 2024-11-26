hess <- function(
    theta, x, time, delta, degree, boundary,
    theta_prev, hess_prev, time_int) {
  p <- ncol(x)
  q <- degree + 1
  alpha <- theta[1:q]
  beta <- theta[(q + 1):(q + p)]
  parms <- list(alpha = alpha, degree = degree, boundary = boundary)

  hess1 <- hess_prev

  u <- unique(time)
  cbh_pre <- as.matrix(deSolve::ode(
    y = 0, times = c(0, u), func = basehaz, parms = parms, method = "ode45"
  ))[-1, -1, drop = FALSE]
  cbh <- cbh_pre[match(time, u), , drop = FALSE]
  dcbh_pre <- as.matrix(deSolve::ode(
    y = rep(0, q), times = c(0, u), func = basehaz_grad, parms = parms,
    method = "ode45"
  ))[-1, -1, drop = FALSE]
  dcbh <- dcbh_pre[match(time, u), , drop = FALSE]
  d2cbh_pre <- as.matrix(deSolve::ode(
    y = rep(0, q * q), times = c(0, u), func = basehaz_hess, parms = parms,
    method = "ode45"
  ))[-1, -1, drop = FALSE]
  d2cbh <- d2cbh_pre[match(time, u), , drop = FALSE]
  eta <- x %*% beta
  d2alpha <- matrix(colSums(sweep(d2cbh, 1, exp(eta), "*")), q)
  d2beta <- t(x) %*% (as.vector(cbh * exp(eta)) * x)
  dalph_dabeta <- t(dcbh) %*% (as.vector(exp(eta)) * x)
  hess2 <- cbind(
    rbind(d2alpha, t(dalph_dabeta)),
    rbind(dalph_dabeta, d2beta)
  )

  if (length(time_int) > 0) {
    u_int <- unique(time_int)
    cbh_int_pre <- as.matrix(deSolve::ode(
      y = 0, times = c(0, u_int), func = basehaz, parms = parms,
      method = "ode45"
    ))[-1, -1, drop = FALSE]
    cbh_int <- cbh_int_pre[match(time_int, u_int), , drop = FALSE]
    dcbh_int_pre <- as.matrix(deSolve::ode(
      y = rep(0, q), times = c(0, u_int), func = basehaz_grad,
      parms = parms, method = "ode45"
    ))[-1, -1, drop = FALSE]
    dcbh_int <- dcbh_int_pre[match(time_int, u_int), , drop = FALSE]
    d2cbh_int_pre <- as.matrix(deSolve::ode(
      y = rep(0, q * q), times = c(0, u_int), func = basehaz_hess,
      parms = parms, method = "ode45"
    ))[-1, -1, drop = FALSE]
    d2cbh_int <- d2cbh_int_pre[match(time_int, u_int), , drop = FALSE]
    eta_int <- eta[names(time_int), , drop = FALSE]
    d2alpha_int <- matrix(colSums(sweep(d2cbh_int, 1, exp(eta_int), "*")), q)
    x_int <- x[names(time_int), , drop = FALSE]
    d2beta_int <- t(x_int) %*% (as.vector(exp(eta_int) * cbh_int) * x_int)
    dalph_dabeta_int <- t(dcbh_int) %*% (as.vector(exp(eta_int)) * x_int)
    hess3 <- -cbind(
      rbind(d2alpha_int, t(dalph_dabeta_int)),
      rbind(dalph_dabeta_int, d2beta_int)
    )
  } else {
    hess3 <- 0
  }

  hess <- hess1 + hess2 + hess3
  return(hess)
}
