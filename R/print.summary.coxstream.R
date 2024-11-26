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

  q <- x$degree + 1
  stats::printCoefmat(x$coefficients[(q + 1):nrow(x$coefficients), ],
    digits = digits, signif.stars = signif.stars,
    cs.ind = 1:3, tst.ind = 4, P.values = TRUE, has.Pvalue = TRUE
  )

  print(
    format(x$conf.int[(q + 1):nrow(x$conf.int), ], digits = digits),
    quote = FALSE
  )

  invisible(x)
}
