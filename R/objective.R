objective <- function(
    theta, x, time, delta, degree, boundary,
    theta_prev, hess_prev, time_int) {
  p <- ncol(x)
  q <- degree + 1
  alpha <- theta[1:q]
  beta <- theta[(q + 1):(q + p)]
  parms <- list(alpha = alpha, degree = degree, boundary = boundary)

  loss1 <- (theta - theta_prev) %*% hess_prev %*% (theta - theta_prev) / 2
  grad1 <- as.vector(hess_prev %*% (theta - theta_prev))
  hess1 <- hess_prev

  u <- unique(time)
  b_pre <- splines2::bpoly(u, degree, TRUE, boundary)
  b <- b_pre[match(time, u), ]
  cbh_pre <- as.matrix(deSolve::ode(
    y = 0, times = c(0, u), func = basehaz_ode, parms = parms, method = "ode45"
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
  dalpha <- t(dcbh) %*% exp(eta) - t(b) %*% delta
  dbeta <- t(x) %*% (cbh * exp(eta) - delta)
  d2alpha <- matrix(colSums(sweep(d2cbh, 1, exp(eta), "*")), q)
  d2beta <- t(x) %*% (as.vector(cbh * exp(eta)) * x)
  dalph_dabeta <- t(dcbh) %*% (as.vector(exp(eta)) * x)

  loss2 <- t(cbh) %*% exp(eta) - t(b %*% alpha + eta) %*% delta
  grad2 <- c(dalpha, dbeta)
  hess2 <- cbind(
    rbind(d2alpha, t(dalph_dabeta)),
    rbind(dalph_dabeta, d2beta)
  )

  if (length(time_int) > 0) {
    u_int <- unique(time_int)
    cbh_int_pre <- as.matrix(deSolve::ode(
      y = 0, times = c(0, u_int), func = basehaz_ode, parms = parms,
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
    x_int <- x[names(time_int), , drop = FALSE]
    eta_int <- eta[names(time_int), ]

    dalpha <- t(dcbh_int) %*% exp(eta_int)
    dbeta <- t(x_int) %*% (cbh_int * exp(eta_int))
    d2alpha_int <- matrix(colSums(sweep(d2cbh_int, 1, exp(eta_int), "*")), q)
    d2beta_int <- t(x_int) %*% (as.vector(exp(eta_int) * cbh_int) * x_int)
    dalph_dabeta_int <- t(dcbh_int) %*% (as.vector(exp(eta_int)) * x_int)

    loss3 <- -sum(cbh_int * exp(eta_int))
    grad3 <- -c(dalpha, dbeta)
    hess3 <- -cbind(
      rbind(d2alpha_int, t(dalph_dabeta_int)),
      rbind(dalph_dabeta_int, d2beta_int)
    )
  } else {
    loss3 <- 0
    grad3 <- 0
    hess3 <- 0
  }

  loss <- loss1 + loss2 + loss3
  grad <- grad1 + grad2 + grad3
  hess <- hess1 + hess2 + hess3
  attr(loss, "gradient") <- grad
  attr(loss, "hessian") <- hess
  return(loss)
}
