#' Streaming Cox model with repeated observations
#'
#' `r lifecycle::badge('experimental')`
#' @param formula A formula expression as for regression models, of the form
#' \code{response ~ predictors}. The response must be a survival object as
#' returned by the \code{\link{Surv}} function.
#' @param data A data frame containing the variables in the model.
#' @param n_basis A positive integer specifying the number of basis functions.
#' @param boundary A vector of length 2 containing the boundary knots.
#' @param idx_col The column name of the index column, which is used to
#' distinguish different patients.
#' @param scale The scaling factor for the pre-estimated coefficients.
#' @param alpha A numeric value used to calculate the number of basis functions.
#' Default is 2.0.
#' @param nu A numeric value used to calculate the number of basis functions.
#' Default is 0.2.
#' @param ... Additional arguments (not used).
#' @return An object of class \code{coxstream} is returned.
#' @export
#'
#' @examples
#' library(coxstream)
#' formula <- Surv(time, status) ~ X1 + X2 + X3 + X4 + X5
#' fit <- coxstream(
#'   formula, sim[sim$batch_id == 1, ],
#'   n_basis = 3, boundary = c(0, 3), idx_col = "patient_id"
#' )
#' for (batch in 2:10) {
#'   fit <- update(fit, sim[sim$batch_id == batch, ])
#' }
#' summary(fit)
coxstream <- function(
    formula, data, n_basis, boundary, idx_col,
    scale = 2, alpha = 2.0, nu = 0.2, ...) {
  if (n_basis %% 1 != 0 || n_basis < 1 ) {
    stop("The number of basis functions must be an positive integer.")
  }
  mf <- stats::model.frame(formula, data)
  y <- stats::model.response(mf)
  time <- y[, 1]
  delta <- y[, 2]
  names(time) <- names(delta) <- data[[idx_col]]
  x <- stats::model.matrix(formula, data)[, -1]
  rownames(x) <- data[[idx_col]]

  # Store the survival times for each patient
  time_stored <- sort(time[delta == 0])
  time_unique <- unique(time)

  sorted <- order(time)
  x <- x[sorted, ]
  time <- time[sorted]
  delta <- delta[sorted]

  n_basis_cur <- n_basis
  n_basis_pre <- round(n_basis * scale)

  n_features <- ncol(x)
  n_params <- n_features + n_basis_pre
  theta_prev <- numeric(n_params)
  hess_prev <- matrix(0, n_params, n_params)
  time_int <- c()

  cur_idx <- c(seq_len(n_basis_cur), (n_basis_pre + 1):n_params)
  res <- trust::trust(
    objfun = objective, parinit = theta_prev[cur_idx], rinit = 1, rmax = 10,
    x = x, time = time, delta = delta, 
    n_basis = n_basis_cur, boundary = boundary,
    theta_prev = theta_prev[cur_idx], hess_prev = hess_prev[cur_idx, cur_idx],
    time_int = time_int
  )
  if (!res$converged) {
    warning("The optimization did not converge.")
  }
  theta_prev[cur_idx] <- res$argument
  res <- objective(
    theta_prev, x, time, delta, n_basis_pre, boundary, theta_prev, 
    hess_prev, time_int
  )
  loss <- res$value
  hess <- res$hessian

  coef_names <- c(paste0("Basis ", seq_len(n_basis_pre)), colnames(x))
  names(theta_prev) <- coef_names
  colnames(hess) <- rownames(hess) <- coef_names

  fit <- list(
    n_basis_cur = n_basis_cur,
    n_basis_pre = n_basis_pre,
    boundary = boundary,
    logLik = loss,
    theta_prev = theta_prev,
    hess_prev = hess,
    time_stored = time_stored,
    time_unique = time_unique,
    formula = formula,
    idx_col = idx_col,
    alpha = alpha,
    nu = nu,
    scale = scale,
    call = match.call()
  )
  class(fit) <- "coxstream"
  fit
}

#' Update the \code{coxstream} with new data.
#'
#' `r lifecycle::badge('experimental')`
#' @param object A \code{coxstream} object.
#' @param data A data frame containing the variables in the model.
#' @param n_basis Either a positive integer specifying the number of basis
#' functions, or the string "auto" to automatically determine this number based
#' \eqn{[\alpha N^{\nu}]}, where \eqn{N} is the number of unique survival times.
#' @param ... Additional arguments (not used).
#'
#' @return An object of class \code{coxstream}.
#' @export
update.coxstream <- function(object, data, n_basis = "auto", ...) {
  formula <- object$formula
  boundary <- object$boundary
  idx_col <- object$idx_col
  scale <- object$scale
  time_stored <- object$time_stored

  mf <- stats::model.frame(formula, data)
  y <- stats::model.response(mf)
  time <- y[, 1]
  delta <- y[, 2]
  names(time) <- names(delta) <- data[[idx_col]]
  x <- stats::model.matrix(formula, data)[, -1]
  rownames(x) <- data[[idx_col]]

  time_unique <- unique(c(object$time_unique, time))

  if (n_basis == "auto") {
    n_basis_cur <- max(
      object$n_basis_cur, round(object$alpha * length(time_unique)^object$nu)
    )
  } else if (is.numeric(n_basis)) {
    n_basis_cur <- as.integer(n_basis)
  } else {
    stop("The n_basis must be an integer or 'auto'.")
  }
  n_basis_pre <- round(n_basis_cur * scale)
  n_features <- ncol(x)
  n_params <- n_basis_pre + n_features

  if (n_basis_pre > object$n_basis_pre) {
    idx <- c(
      seq_len(object$n_basis_pre), n_basis_pre + seq_len(n_features)
    )
    theta_prev <- rep(0, n_params)
    theta_prev[idx] <- object$theta_prev
    hess_prev <- matrix(0, n_params, n_params)
    hess_prev[idx, idx] <- object$hess_prev
  } else {
    theta_prev <- object$theta_prev
    hess_prev <- object$hess_prev
  }

  # The patients that are already stored and the new patients
  patients_stored <- names(time_stored)
  patients_int <- intersect(patients_stored, data[[idx_col]])
  patients_new <- data[[idx_col]][!data[[idx_col]] %in% patients_int]
  # Select the stored survival times that int with the new patients
  time_int <- time_stored[patients_int]
  # Update the stored survival times
  time_stored[patients_int] <- time[patients_int]
  # Select the survival patients in the inted stored survival times
  delta_int <- delta[patients_int]
  patients_remove <- names(delta_int[delta_int == 1])
  time_stored <- time_stored[!names(time_stored) %in% patients_remove]
  time_new <- time[data[[idx_col]] %in% patients_new & delta == 0]
  time_stored <- c(time_stored, time_new)
  time_stored <- sort(time_stored)

  sorted <- order(time)
  x <- x[sorted, ]
  time <- time[sorted]
  delta <- delta[sorted]

  cur_idx <- c(seq_len(n_basis_cur), (n_basis_pre + 1):n_params)
  res <- trust::trust(
    objfun = objective, parinit = theta_prev[cur_idx], rinit = 1, rmax = 10,
    x = x, time = time, delta = delta, 
    n_basis = n_basis_cur, boundary = boundary,
    theta_prev = theta_prev[cur_idx], hess_prev = hess_prev[cur_idx, cur_idx],
    time_int = time_int
  )
  if (!res$converged) {
    warning("The optimization did not converge.")
  }
  theta_prev[cur_idx] <- res$argument
  object$theta_prev <- theta_prev

  res <- objective(
    theta_prev, x, time, delta, n_basis_pre, boundary, theta_prev, 
    hess_prev, time_int
  )
  object$hess_prev <- res$hessian

  if (n_basis_cur < n_basis_pre) {
    object$theta_prev[(n_basis_cur + 1):n_basis_pre] <- 0
  }
  coef_names <- c(paste0("Basis ", seq_len(n_basis_pre)), colnames(x))
  names(object$theta_prev) <- coef_names
  colnames(object$hess_prev) <- rownames(object$hess_prev) <- coef_names

  object$time_stored <- time_stored
  object$time_unique <- time_unique
  object$n_basis_cur <- n_basis_cur
  object$n_basis_pre <- n_basis_pre
  class(object) <- "coxstream"
  object
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
    stop("The object must be of class 'coxstream'.")
  }
  idx <- c(
    seq_len(object$n_basis_cur),
    (object$n_basis_pre + 1):length(object$theta_prev)
  )
  object$theta_prev[idx]
}

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

  idx <- c(
    seq_len(object$n_basis_cur),
    (object$n_basis_pre + 1):length(object$theta_prev)
  )
  solve(object$hess_prev[idx, idx])
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
    n_basis_cur = object$n_basis_cur, boundary = object$boundary
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
#' @export
print.summary.coxstream <- function(
    x, digits = max(3, getOption("digits") - 3),
    signif.stars = getOption("show.signif.stars"), ...) {
  cat("Call:\n")
  print(x$call)
  cat("\n")
  cat("Number of basis: ", x$n_basis_cur, "\n")

  n_basis <- x$n_basis_cur
  stats::printCoefmat(x$coefficients[(n_basis + 1):nrow(x$coefficients), ],
    digits = digits, signif.stars = signif.stars,
    cs.ind = 1:3, tst.ind = 4, P.values = TRUE, has.Pvalue = TRUE
  )

  print(
    format(x$conf.int[(n_basis + 1):nrow(x$conf.int), ], digits = digits),
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

  n_basis <- object$n_basis_cur
  boundary <- object$boundary

  if (!is.null(newdata)) {
    time <- newdata
  } else {
    time <- seq.int(boundary[1], boundary[2], length.out = 100)
  }
  alpha <- object$theta_prev[seq_len(n_basis)]
  parms <- list(alpha = alpha, n_basis = n_basis, boundary = boundary)

  u <- unique(time)
  b_pre <- chebyshev(u, n_basis, boundary)
  b <- b_pre[match(time, u), ]
  cbh_pre <- as.matrix(deSolve::ode(
    y = 0, times = c(0, u), func = basehaz_ode, parms = parms, method = "ode45"
  ))[-1, -1, drop = FALSE]
  cbh <- cbh_pre[match(time, u), , drop = FALSE]

  vcov_alpha <- vcov(object)[1:n_basis, 1:n_basis]
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

#' Cross-validation for \code{coxstream} objects
#' @param formula A formula expression as for regression models, of the form
#' \code{response ~ predictors}. The response must be a survival object as
#' returned by the \code{\link{Surv}} function.
#' @param data A data frame containing the variables in the model.
#' @param boundary A vector of length 2 containing the boundary knots.
#' @param idx_col The column name of the index column, which is used to
#' distinguish different patients.
#' @param nu A numeric value used to calculate the number of basis functions.
#' Default is 0.2.
#' @param n_folds A positive integer specifying the number of folds for
#' cross-validation. Default is 10.
#' @param n_alphas A positive integer specifying the number of alpha values to
#' test. Default is 5
#' @param n_basis_min A positive integer specifying the minimum number of
#' basis functions. Default is 3.
#' @param seed An integer specifying the random seed for reproducibility.
#' Default is 0.
#' @return An object of class \code{cv.coxstream} is returned.
#' @export
cv.coxstream <- function(
    formula, data, boundary, idx_col,
    nu = 0.2, n_folds = 10, n_alphas = 5, n_basis_min = 3, seed = 0) {
  set.seed(seed)
  folds <- sample(seq_len(n_folds), nrow(data), replace = TRUE)
  n_basises <- seq_len(n_alphas) + n_basis_min - 1

  cv_score <- matrix(0, nrow = n_alphas, ncol = n_folds)
  for (fold in seq_len(n_folds)) {
    data_train <- data[folds != fold, ]

    data_test <- data[folds == fold, ]
    y_test <- stats::model.response(stats::model.frame(formula, data_test))

    for (i in seq_along(n_basises)) {
      fit <- coxstream(
        formula, data_train,
        n_basis = n_basises[i],
        boundary = boundary, idx_col = idx_col, scale = 1.0
      )
      df_test <- data.frame(
        time = y_test[, 1], status = y_test[, 2],
        pred = -predict(fit, newdata = data_test, type = "lp")
      )
      cv_score[i, fold] <-
        concordance(Surv(time, status) ~ pred, df_test)$concordance
    }
  }
  n_uniques <- length(unique(data$time)) * (n_folds - 1) / n_folds
  alpha <- n_basises / n_uniques^nu
  cvm <- apply(cv_score, 1, mean)
  cvsd <- apply(cv_score, 1, stats::sd)
  index <- c(which.max(cvm), which.max(cvm - cvsd))

  object <- list(
    alpha = alpha,
    cvm = cvm,
    cvsd = cvsd,
    cvup = cvm + cvsd,
    cvlo = cvm - cvsd,
    nbasis = n_basises,
    alpha.min = alpha[index[1]],
    alpha.1se = alpha[index[2]]
  )
  class(object) <- "cv.coxstream"
  object
}
