#' Posterior draws from a Bayesian calibration
#'
#' Returns the posterior draws of a Bayesian calibration as a \code{data.frame}, in the same shape as the results of the bootstrap calibrations (\code{cal.ols()} and others). The result can be passed to \code{rec.bayesian()}, \code{rec.clumped()}, or \code{cal.ci()}.
#'
#' @param fit A result from \code{cal.bayesian()} or \code{cal.ion.bayesian()}, or a single \code{stanfit} from it.
#' @param model The model to use when \code{fit} holds more than one. Use \code{"errors"}, \code{"noerrors"}, \code{"mixed_errors"}, or \code{"mixed_noerrors"}, or the name of the element (such as \code{"BLM1_fit"}).
#' @param newMaterial For mixed models, whether to return draws for a material that is not in the calibration. These draws combine the uncertainty in the shared means with the spread among materials.
#'
#' @return A \code{data.frame} with one row per draw and the columns \code{alpha}, \code{beta}, and \code{sigma}. Models with a second predictor also have \code{gamma}. Mixed models have one set of draws per material, stacked, with a \code{Material} column. With \code{newMaterial = TRUE}, there is a single set of draws and no \code{Material} column.
#'
#' @examples
#' \donttest{
#' ds <- cal.dataset(error = "S1", nobs = 50, seed = 1)
#' fit <- cal.bayesian(ds, models = "errors", numSavedSteps = 1000)
#' head(cal.draws(fit))
#' }
#'
#' @export

cal.draws <- function(fit, model = NULL, newMaterial = FALSE) {
  if (inherits(fit, "bayclump_cal")) {
    if (is.null(model)) {
      if (length(fit) > 1) {
        stop("Choose a model with the model argument: ", paste(names(fit), collapse = ", "), ".")
      }
      model <- names(fit)[1]
    }
    if (model == "mixed") model <- "mixed_noerrors"
    key <- if (model %in% names(.model_names)) .model_names[[model]] else model
    if (!key %in% names(fit)) {
      stop("The model ", model, " was not fit.")
    }
    fit <- fit[[key]]
  }
  if (!inherits(fit, "stanfit")) {
    stop("fit must be a result from cal.bayesian() or cal.ion.bayesian(), or a stanfit from it.")
  }

  # Draws for a material that is not in the calibration
  if (newMaterial) {
    if (!"alpha_new" %in% fit@sim$pars_oi) {
      stop("newMaterial = TRUE needs a mixed model fitted with this version of bayclumpr.")
    }
    p <- rstan::extract(fit, pars = c("alpha_new", "beta_new", "sigma"))
    return(data.frame(alpha = as.vector(p$alpha_new), beta = as.vector(p$beta_new), sigma = as.vector(p$sigma)))
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
