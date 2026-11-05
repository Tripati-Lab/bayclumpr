#' Posterior draws from a Bayesian calibration
#'
#' Returns the posterior draws of a Bayesian calibration as a \code{data.frame}, in the same shape as the results of the bootstrap calibrations (\code{cal.ols()} and others). The result can be passed to \code{rec.bayesian()}, \code{rec.clumped()}, or \code{cal.ci()}.
#'
#' @param fit A result from \code{cal.bayesian()} or \code{cal.ion.bayesian()}, or a single \code{stanfit} from it.
#' @param model The model to use when \code{fit} holds more than one. Use \code{"errors"}, \code{"noerrors"}, or \code{"mixed"}, or the name of the element (such as \code{"BLM1_fit"}).
#'
#' @return A \code{data.frame} with one row per draw and the columns \code{alpha}, \code{beta}, and \code{sigma}. Models with a second predictor also have \code{gamma}. Mixed models have one set of draws per material, stacked, with a \code{Material} column.
#'
#' @export

cal.draws <- function(fit, model = NULL) {
  if (inherits(fit, "bayclump_cal")) {
    if (is.null(model)) {
      if (length(fit) > 1) {
        stop("Choose a model with the model argument: ", paste(names(fit), collapse = ", "), ".")
      }
      model <- names(fit)[1]
    }
    key <- if (model %in% names(.model_names)) .model_names[[model]] else model
    if (!key %in% names(fit)) {
      stop("The model ", model, " was not fit.")
    }
    fit <- fit[[key]]
  }
  if (!inherits(fit, "stanfit")) {
    stop("fit must be a result from cal.bayesian() or cal.ion.bayesian(), or a stanfit from it.")
  }

  pars <- intersect(c("alpha", "beta", "gamma", "sigma"), fit@sim$pars_oi)
  p <- rstan::extract(fit, pars = pars)

  # Mixed models have one column of draws per material
  if (is.matrix(p$alpha)) {
    mats <- attr(fit, "Material")
    if (is.null(mats)) mats <- seq_len(ncol(p$alpha))
    draws <- lapply(seq_along(mats), function(j) {
      data.frame(alpha = p$alpha[, j], beta = p$beta[, j], sigma = as.vector(p$sigma), Material = mats[j])
    })
    return(do.call(rbind, draws))
  }

  as.data.frame(lapply(p, as.vector))
}
