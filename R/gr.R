gr <- function(
    theta, x, time, delta, degree, boundary,
    theta_prev, hess_prev, time_int) {
  p <- ncol(x)
  q <- degree + 1
  alpha <- theta[1:q]
  beta <- theta[(q + 1):(q + p)]
  parms <- list(alpha = alpha, degree = degree, boundary = boundary)

  dtheta1 <- as.vector(hess_prev %*% (theta - theta_prev))

  u <- unique(time)
  b_pre <- splines2::bpoly(u, degree, TRUE, boundary)
  b <- b_pre[match(time, u), ]
  cbh_pre <- as.matrix(deSolve::ode(
    y = 0, times = c(0, u), func = basehaz, parms = parms, method = "ode45"
  ))[-1, -1, drop = FALSE]
  cbh <- cbh_pre[match(time, u), , drop = FALSE]
  dcbh_pre <- as.matrix(deSolve::ode(
    y = rep(0, q), times = c(0, u), func = basehaz_grad, parms = parms,
    method = "ode45"
  ))[-1, -1, drop = FALSE]
  dcbh <- dcbh_pre[match(time, u), , drop = FALSE]
  eta <- (x %*% beta)
  dalpha <- t(dcbh) %*% exp(eta) - t(b) %*% delta
  dbeta <- t(x) %*% (cbh * exp(eta) - delta)
  dtheta2 <- c(dalpha, dbeta)

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
    eta_int <- eta[names(time_int), , drop = FALSE]
    x_int <- x[names(time_int), , drop = FALSE]
    dalpha <- t(dcbh_int) %*% exp(eta_int)
    dbeta <- t(x_int) %*% (cbh_int * exp(eta_int))
    dtheta3 <- -c(dalpha, dbeta)
  } else {
    dtheta3 <- 0
  }

  dtheta1 + dtheta2 + dtheta3
}
