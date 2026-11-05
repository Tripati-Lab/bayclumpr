#' Average replicate measurements of D47
#'
#' Averages the D47 replicates of each sample and returns the standard error of each mean. The result can be used as a calibration or reconstruction dataset.
#'
#' @param data A \code{data.frame} with one row per replicate and the columns \code{Sample} and \code{D47}.
#' @param pooledSD Optional SD used for every sample instead of the SD of its own replicates, such as the long-term SD of standards measured in the same lab.
#'
#' @details \code{D47error} is the standard error of the mean, SD / sqrt(N). By default, the SD comes from the replicates of each sample. Samples with a single replicate then have no SE, so a warning is given. Other columns are kept when they have a single value within every sample, such as \code{Material} or \code{Temperature}.
#'
#' @return A \code{data.frame} with one row per sample and the columns \code{Sample}, \code{D47} (mean), \code{D47sd} (SD of the replicates), \code{N} (number of replicates), and \code{D47error} (standard error of the mean), followed by the columns that were kept.
#'
#' @examples
#' reps <- data.frame(
#'   Sample = rep(c("A", "B"), each = 3),
#'   D47 = c(0.601, 0.598, 0.605, 0.702, 0.695, 0.699),
#'   Material = "calcite"
#' )
#' avg.replicates(reps)
#' avg.replicates(reps, pooledSD = 0.008)
#'
#' @importFrom stats sd
#'
#' @export

avg.replicates <- function(data, pooledSD = NULL) {
  if (!all(c("Sample", "D47") %in% names(data))) {
    stop("data needs the columns Sample and D47.")
  }
  if (anyNA(data$Sample) || anyNA(data$D47)) {
    stop("Sample and D47 can't have missing values.")
  }
  if (!is.null(pooledSD) && !(length(pooledSD) == 1 && pooledSD > 0)) {
    stop("pooledSD must be a single positive number.")
  }

  # Keep samples in the order they first appear
  samples <- factor(data$Sample, levels = unique(data$Sample))

  # Columns with a single value within every sample are kept
  other <- setdiff(names(data), c("Sample", "D47", "D47error"))
  keep <- Filter(function(v) all(tapply(data[[v]], samples, function(z) length(unique(z)) == 1)), other)

  out <- do.call(rbind, lapply(split(data, samples), function(d) {
    n <- nrow(d)
    s <- if (is.null(pooledSD)) sd(d$D47) else pooledSD
    cbind(
      data.frame(Sample = d$Sample[1], D47 = mean(d$D47), D47sd = sd(d$D47), N = n, D47error = s / sqrt(n)),
      d[1, keep, drop = FALSE]
    )
  }))
  row.names(out) <- NULL

  if (anyNA(out$D47error)) {
    warning("Samples with a single replicate have no SE. Set pooledSD to use a common SD for all samples.")
  }
  out
}
