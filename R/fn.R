fn <- function(
    theta, x, time, delta, degree, boundary,
    theta_prev, hess_prev, time_int) {
  p <- ncol(x)
  q <- degree + 1
  alpha <- theta[1:q]
  beta <- theta[(q + 1):(q + p)]
  parms <- list(alpha = alpha, degree = degree, boundary = boundary)

  l1 <- (theta - theta_prev) %*% hess_prev %*% (theta - theta_prev) / 2

  u <- unique(time)
  b_pre <- splines2::bpoly(u, degree, TRUE, boundary)
  b <- b_pre[match(time, u), ]
  cbh_pre <- as.matrix(deSolve::ode(
    y = 0, times = c(0, u), func = basehaz, parms = parms, method = "ode45"
  ))[-1, -1, drop = FALSE]
  cbh <- cbh_pre[match(time, u), , drop = FALSE]
  eta <- x %*% beta
  l2 <- t(cbh) %*% exp(eta) - t(b %*% alpha + eta) %*% delta

  if (length(time_int) > 0) {
    u_int <- unique(time_int)
    cbh_int_pre <- as.matrix(deSolve::ode(
      y = 0, times = c(0, u_int), func = basehaz, parms = parms,
      method = "ode45"
    ))[-1, -1, drop = FALSE]
    cbh_int <- cbh_int_pre[match(time_int, u_int), , drop = FALSE]
    eta_int <- eta[names(time_int), ]
    l3 <- -sum(cbh_int * exp(eta_int))
  } else {
    l3 <- 0
  }

  l1 + l2 + l3
}
