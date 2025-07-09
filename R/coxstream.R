#' @title Alias for Survival Function
#' @description This function creates an alias for the `Surv` function from the
#' `survival` package, allowing it to be used without explicitly referencing
#' the package namespace.
#'
#' @param time The follow-up time for the subject.
#' @param time2 The ending time for the interval (if applicable).
#' @param event The status indicator, normally 0=alive, 1=dead.
#'              Other values are allowed.
#' @param type Character string specifying the type of censoring.
#'             Possible values are "right", "left", "interval", or "counting".
#' @param origin The origin for counting time, used when `type = "counting"`.
#'
#' @return A `Surv` object representing the survival data.
#'
#' @importFrom survival Surv
#' @seealso \code{\link[survival]{Surv}}
#' @export
Surv <- survival::Surv # nolint: object_name_linter.

#' Streaming Cox model with repeated observations
#'
#' `r lifecycle::badge('experimental')`
#'
#' Fits a streaming Cox proportional hazards model that accommodates repeated
#' observations for subjects.
#'
#' @param formula A formula expression as for regression models, of the form
#' \code{response ~ predictors}. The response must be a survival object as
#' returned by the \code{\link{Surv}} function.
#' @param data A data frame containing the variables in the model.
#' @param n_basis A positive integer specifying the number of basis functions.
#' @param boundary A vector of length 2 containing the boundary knots.
#' @param subject_col The column name of the index column, which is used to
#' distinguish different patients.
#' @param scale The scaling factor for the pre-estimated coefficients.
#' @param ... Additional arguments (not used).
#' @return Additional arguments passed to the optimization function.
#'
#' @return An object of class \code{coxstream} containing the fitted model.
#' The object contains the following components:
#' \describe{
#'   \item{logLik}{The log-likelihood of the fitted model.}
#'   \item{theta}{The estimated coefficients.}
#'   \item{hessian}{The Hessian matrix of the fitted model.}
#'   \item{n_basis}{The number of basis functions used.}
#'   \item{n_features}{The number of features in the model.}
#'   \item{n_samples}{The number of samples in the model.}
#'   \item{formula}{The formula used to fit the model.}
#'   \item{boundary}{The boundary knots used in the model.}
#'   \item{subject_col}{The index column used to distinguish different
#'                       patients.}
#'   \item{scale}{The scaling factor used in the model.}
#'   \item{time_stored}{The survival times stored for each patient.}
#'   \item{call}{The matched call.}
#' }
#'
#' @details This function employs the `nlm` optimization method to estimate
#'   the model parameters. During the pre-estimation stage, parameters are
#'   projected to mitigate bias introduced in the early stage.
#'
#' @examples
#' formula <- Surv(time, status) ~ X1 + X2 + X3 + X4 + X5
#' fit <- coxstream(
#'   formula, sim[sim$batch_id == 1, ],
#'   n_basis = 5, boundary = c(0, 3), subject_col = "patient_id"
#' )
#' for (batch in 2:10) {
#'   fit <- update(fit, sim[sim$batch_id == batch, ])
#' }
#' summary(fit)
#'
#' @seealso \code{\link{update.coxstream}}, \code{\link{summary.coxstream}}.
#'
#' @importFrom stats model.frame model.matrix model.response nlm
#'
#' @export
coxstream <- function(
    formula, data, n_basis, boundary, subject_col, scale = 2.0, ...) {
  # Input validation
  if (!inherits(formula, "formula")) stop("formula must be a formula object")
  if (!inherits(data, "data.frame")) stop("data must be a data frame")
  if (!is.numeric(n_basis) || length(n_basis) != 1 || n_basis <= 0) {
    stop("n_basis must be a positive integer")
  }
  if (
    !is.numeric(boundary) || length(boundary) != 2 || boundary[1] >= boundary[2]
  ) {
    stop("boundary must be a numeric vector of length 2 with increasing values")
  }
  if (!is.data.frame(data) || !all(subject_col %in% colnames(data))) {
    stop("subject_col must be a column name in the data frame")
  }
  if (!is.numeric(scale) || length(scale) != 1 || scale <= 0) {
    stop("scale must be a positive numeric value")
  }

  # Prepare model frame and extract response and predictors
  mf <- stats::model.frame(formula, data)
  y <- stats::model.response(mf)
  time <- y[, 1, drop = FALSE]
  status <- y[, 2, drop = FALSE]
  x <- stats::model.matrix(formula, data)[, -1]
  subject <- data[[subject_col]]

  # Initialize parameters
  n_observations <- nrow(data)
  n_features <- ncol(x)

  n_basis_pre <- round(n_basis * scale)
  n_parameters <- n_features + n_basis_pre

  # Store the survival times for each patient
  time_stored <- tapply(time, subject, max, na.rm = TRUE)
  is_event <- tapply(status, subject, max, na.rm = TRUE)
  time_stored <- time_stored[!subject[is_event == 1] %in% names(time_stored)]

  # Initialize the Parameters
  theta_init <- numeric(n_basis + n_features)
  theta <- numeric(n_parameters)
  hess <- matrix(0, n_parameters, n_parameters)
  time_int <- c()

  ## 1. Dynamic Activation Stage
  res <- nlm(
    f = objective, p = theta_init,
    time = time, status = status, x = x, subject = subject,
    boundary = boundary, theta_prev = theta, hess_prev = hess,
    time_int = time_int, ...
  )

  # 2. Pre-estimation stage
  prox <- prox_forward(n_basis, n_basis_pre, n_features)
  theta <- as.vector(prox %*% res$estimate)
  names(theta) <- c(paste0("Basis ", seq_len(n_basis_pre)), colnames(x))
  obj <- objective(
    theta, time, status, x, subject, boundary, theta, hess, time_int
  )

  object <- list(
    logLik = as.numeric(-obj),
    theta = theta,
    hessian = attr(obj, "hessian"),
    n_basis = c(n_basis),
    n_features = n_features,
    n_observations = n_observations,
    formula = formula,
    boundary = boundary,
    subject_col = subject_col,
    scale = scale,
    time_stored = time_stored,
    call = match.call()
  )
  class(object) <- "coxstream"
  object
}

#' Update the \code{coxstream} with new data.
#'
#' `r lifecycle::badge('experimental')`
#'
#' This function updates the \code{coxstream} object with new data, allowing
#' the dynamic activation of the newly added spline basis functions.
#'
#' @param object A \code{coxstream} object.
#' @param newdata A \code{data.frame} containing the new data to update the
#'   model. It must include the variables specified in the original model
#'   formula.
#' @param n_basis Either a numeric value specifying the number of basis
#'   functions or \code{"auto"} to automatically determine the number of basis
#'   functions based on the sample size and scaling parameters.
#'   Defaults to \code{"auto"}, which means the number of basis functions will
#'   be determined by n_basis = round(alpha * n_samples^nu).
#' @param alpha A positive numeric value controlling the growth rate of the
#'   number of basis functions when \code{n_basis = "auto"}. Defaults to 1.0.
#' @param nu A positive numeric value controlling the scaling exponent for
#'   determining the number of basis functions when \code{n_basis = "auto"}.
#'   Defaults to 0.2.
#' @param ... Additional arguments passed to the optimization function.
#'
#' @details
#'
#' @return An object of class \code{coxstream}.
#'
#' @importFrom stats update
#' @importFrom stats model.frame model.response model.matrix nlm
#' @importFrom utils tail
#'
#' @export
update.coxstream <- function(
    object, newdata, n_basis = "auto", alpha = 1.0, nu = 0.2, ...) {
  if (!inherits(object, "coxstream")) {
    stop("object must be of class 'coxstream'")
  }
  if (!inherits(newdata, "data.frame")) stop("newdata must be a data frame")
  if (!is.numeric(n_basis) && n_basis != "auto") {
    stop("n_basis must be numeric or 'auto'")
  }
  if (!is.numeric(alpha) || length(alpha) != 1 || alpha <= 0) {
    stop("alpha must be a positive numeric value")
  }
  if (!is.numeric(nu) || length(nu) != 1 || nu <= 0) {
    stop("nu must be a positive numeric value")
  }

  # Retrieve Parameters
  theta <- object$theta
  hess <- object$hessian

  n_observations <- object$n_observations
  n_features <- object$n_features
  n_basis_last <- tail(object$n_basis, 1)
  n_basis_pre_last <- length(theta) - n_features

  formula <- object$formula
  boundary <- object$boundary
  scale <- object$scale
  subject_col <- object$subject_col
  time_stored <- object$time_stored

  mf <- stats::model.frame(formula, newdata)
  y <- stats::model.response(mf)
  time <- y[, 1, drop = FALSE]
  status <- y[, 2, drop = FALSE]
  x <- stats::model.matrix(formula, newdata)[, -1]
  subject <- newdata[[subject_col]]

  # Update Parameters
  n_observations <- n_observations + nrow(x)
  if (n_basis == "auto") {
    n_basis <- max(n_basis_last, round(alpha * n_observations^nu))
  } else if (is.numeric(n_basis)) {
    n_basis <- max(n_basis_last, n_basis)
  } else {
    stop("n_basis must be numeric or 'auto'")
  }
  object$n_observations <- n_observations
  object$n_basis <- c(object$n_basis, n_basis)

  time_int <- time_stored[intersect(names(time_stored), subject)]

  time_max <- tapply(time, subject, max, na.rm = TRUE)
  is_event <- tapply(status, subject, max, na.rm = TRUE)
  time_stored <- c(time_stored, time_max)
  time_stored <- time_stored[!duplicated(names(time_stored), fromLast = TRUE)]
  object$time_stored <- time_stored[
    !names(time_stored) %in% subject[is_event == 1]
  ]

  ## 1. Dynamic Activation Stage
  n_basis_pre <- round(n_basis * scale)
  if (n_basis_pre > n_basis_pre_last) {
    prox <- prox_forward(n_basis_pre_last, n_basis_pre, n_features)
    theta <- as.vector(prox %*% theta)
    prox_inv <- MASS::ginv(prox)
    hess <- t(prox_inv) %*% hess %*% prox_inv
  }
  prox <- prox_forward(n_basis, n_basis_pre, n_features)
  theta_init <- qr.solve(prox, theta)
  res <- nlm(
    f = objective, p = theta_init,
    time = time, status = status, x = x, subject = subject,
    boundary = boundary, theta_prev = theta, hess_prev = hess,
    time_int = time_int, ...
  )

  # 2. Pre-estimation stage
  prox <- prox_forward(n_basis, n_basis_pre, n_features)
  theta <- as.vector(prox %*% res$estimate)
  names(theta) <- c(paste0("Basis ", seq_len(n_basis_pre)), colnames(x))
  obj <- objective(
    theta, time, status, x, subject, boundary, theta, hess, time_int
  )

  object$logLik <- as.numeric(-obj)
  object$theta <- theta
  object$hessian <- attr(obj, "hessian")
  object$call <- match.call()
  object
}

#' Compute Log-Likelihood for 'coxstream' Objects
#'
#' This function retrieves the log-likelihood value from an object of class
#' 'coxstream'.
#'
#' @param object An object of class 'coxstream'.
#' @param ... Additional arguments (currently unused).
#'
#' @return The log-likelihood value stored in the 'coxstream' object.
#'
#' @importFrom stats logLik
#'
#' @export
logLik.coxstream <- function(object, ...) {
  if (!inherits(object, "coxstream")) {
    stop("object must be of class 'coxstream'")
  }
  object$logLik
}

#' Calculate the Akaike Information Criterion (AIC) for a 'coxstream' object
#'
#' This function computes the AIC for a given 'coxstream' object. The AIC is a
#' measure of the relative quality of a statistical model for a given dataset.
#' It is calculated as: \eqn{-2 * logLik + 2 * k}, where \code{k} is the number
#' of parameters in the model.
#'
#' @param object An object of class \code{"coxstream"} for which the AIC is to
#'   be calculated.
#' @param ... Additional arguments (currently not used).
#'
#' @return A numeric value representing the AIC of the model.
#'
#' @importFrom stats AIC
#'
#' @seealso \code{\link{logLik.coxstream}}
#'
#' @export
AIC.coxstream <- function(object, ...) {
  if (!inherits(object, "coxstream")) {
    stop("object must be of class 'coxstream'")
  }
  n_parameters <- tail(object$n_basis, 1) + object$n_features
  -2 * logLik(object) + 2 * n_parameters
}

#' Extract the coefficients from a \code{coxstream} object
#'
#' @param object An object of class \code{coxstream}.
#' @param ... Additional arguments (not used).
#'
#' @return A numeric vector of coefficients.
#' @export
coef.coxstream <- function(object, ...) {
  if (!inherits(object, "coxstream")) {
    stop("object must be of class 'coxstream'")
  }
  n_basis_pre <- length(object$theta) - object$n_features
  n_basis <- tail(object$n_basis, 1)
  prox <- prox_forward(n_basis, n_basis_pre, object$n_features)
  theta <- qr.solve(prox, object$theta)

  names(theta) <- c(
    paste0("Basis ", seq_len(n_basis)),
    names(object$theta[(n_basis_pre + 1):length(object$theta)])
  )
  theta
}

#' Extract variance-covariance matrix from a \code{coxstream} object
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
  n_basis_pre <- length(object$theta) - object$n_features
  n_basis <- tail(object$n_basis, 1)
  prox <- prox_forward(n_basis, n_basis_pre, object$n_features)

  hessian <- object$hessian
  hess <- t(prox) %*% hessian %*% prox
  solve(hess)
}

#' Summary method for \code{coxstream} object
#'
#' @param object An object of class \code{coxstream}.
#' @param conf.int A numeric value between 0 and 1 indicating the confidence
#' level for the confidence interval. Default is 0.95.
#' @param ... Additional arguments (not unused).
#'
#' @return An object of class \code{"summary.coxstream"}, containing the
#' following components:
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
    n_basis = tail(object$n_basis, 1), boundary = object$boundary
  )
  class(summary_list) <- "summary.coxstream"
  return(summary_list)
}

#' Print method for \code{summary.coxstream} object
#'
#' @param x A summary object produced from a fitted \code{coxstream} model. This
#' object contains information such as model coefficients, the log-likelihood,
#' BIC, and confidence intervals.
#' @param digits An integer controlling the number of significant digits to
#' print for numeric values.
#' @param signif.stars Logical; if \code{TRUE}, significance stars are printed
#' along with the p-values.
#' @param ... Additional arguments (not used).
#' @return The function prints the summary of the \code{coxstream} model and
#' returns the object \code{x} invisibly.
#'
#' @details The function provides a formatted output that includes:
#' \itemize{
#'   \item \strong{Call:} The original function call that produced the model.
#'   \item \strong{Coefficients:} The regression coefficients, their standard
#'         errors, z-values, and p-values, formatted in a table. Significance
#'         stars are shown next to p-values if \code{signif.stars} is
#'         \code{TRUE}.
#'   \item \strong{Confidence intervals:} The exponentiated coefficients along
#'         with their confidence intervals.
#' }
#'
#' @importFrom stats printCoefmat
#'
#' @export
print.summary.coxstream <- function(
    x, digits = max(3, getOption("digits") - 3),
    signif.stars = getOption("show.signif.stars"), ...) {
  cat("Call:\n")
  print(x$call)
  cat("\nNumber of basis functions: ", x$n_basis, "\n\n")

  idx <- !grepl("^Basis", rownames(x$coefficients))
  printCoefmat(
    x$coefficients[idx, ],
    digits = digits, signif.stars = signif.stars,
    cs.ind = 1:3, tst.ind = 4, P.values = TRUE, has.Pvalue = TRUE
  )
  print(
    format(x$conf.int[idx, ], digits = digits),
    quote = FALSE
  )

  invisible(x)
}


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
  beta <- object$theta_prev[(object$n_basis_pre + 1):length(object$theta_prev)]
  lp <- x %*% beta
  if (type == "lp") {
    return(lp)
  } else if (type == "risk") {
    return(exp(lp))
  } else {
    stop("type must be one of 'lp' or 'risk'")
  }
}

#' Predict the cumulative baseline hazard function for \code{coxstream} objects
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
#' the cumulative baseline hazard function, standard error, and confidence
#' interval.
#' @export
basehaz.coxstream <- function(
    object, newdata = NULL, conf.int = 0.95, ...) {
  if (!inherits(object, "coxstream")) {
    stop("object must be of class 'coxstream'")
  }

  n_basis <- tail(object$n_basis, 1)
  boundary <- object$boundary

  if (!is.null(newdata)) {
    time <- newdata
  } else {
    time <- seq.int(boundary[1], boundary[2], length.out = 100)
  }
  alpha <- coef(object)[seq_len(n_basis)]
  parms <- list(alpha = alpha, n_basis = n_basis, boundary = boundary)

  u <- unique(time)
  b_uniq <- bernstein(u, n_basis, boundary)
  b <- b_uniq[match(time, u), ]
  cbh_uniq <- as.matrix(deSolve::ode(
    y = 0, times = c(0, u), func = basehaz_ode, parms = parms, method = "ode45"
  ))[-1, -1, drop = FALSE]
  cbh <- cbh_uniq[match(time, u), , drop = FALSE]

  vcov_alpha <- vcov(object)[seq_len(n_basis), seq_len(n_basis)]
  vcov_lp <- b %*% vcov_alpha %*% t(b)
  vcov_basehaz <- diag(as.vector(cbh^2)) %*% vcov_lp
  se_basehaz <- sqrt(diag(vcov_basehaz))
  z <- stats::qnorm((1 + conf.int) / 2)

  basehaz_matrix <- cbind(
    time, cbh, se_basehaz,
    cbh - z * se_basehaz, cbh + z * se_basehaz
  )
  dimnames(basehaz_matrix) <- list(
    NULL, c(
      "time", "hazard", "se(hazard)",
      paste("lower .", round(100 * conf.int, 2), sep = ""),
      paste("upper .", round(100 * conf.int, 2), sep = "")
    )
  )
  basehaz_matrix
}
