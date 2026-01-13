#' Fit a Trend-MoMo model with Poisson GLM
#' @param model A \code{"trendMoMo"} object from \code{\link{trendMoMo}}.
#' @param data Data frame with columns \code{x} (age), \code{t} (year),
#'   \code{D} (deaths), \code{E} (exposure). Additional columns are created internally.
#' @return A list with components: \code{fit} (glm), \code{model}, \code{ages.fit}, \code{years.fit}.
#' @export
trendModelFit <- function(model, data){
  stopifnot(inherits(model, "trendMoMo"))
  stopifnot(all(c("x", "t", "D", "E") %in% names(data)))

  #Add offset if not present
  if (!"o" %in% names(data)) {
    data$o <- 0
  }
  #add data$c if cohort effect is present
  if (!is.null(model$cohortAgeFun)) {
    data$c <- data$t - data$x
  }
  #add data$tj if trend changes are present
  if (!is.null(model$trendChanges)) {
    for (j in seq_along(model$trendChanges)) {
      change_year <- model$trendChanges[j]
      data[[paste0("t", j)]] <- pmax(data$t - change_year, 0)
    }
  }
  fit <- glm(formula = model$gnmFormula,
             family = poisson(link = "log"), data = data)
  ages.fit <- sort(unique(data$x))
  years.fit <- sort(unique(data$t))
  out <- list(fit = fit, model = model, ages.fit = ages.fit,
              years.fit = years.fit)
  out
}
