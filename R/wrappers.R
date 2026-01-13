#' Fit–extract–constrain (sequential constraints)
#' @inheritParams trendModelFit
#' @return A list with \code{modelFit}, raw \code{coef}, and adjusted \code{coef_adj}.
#' @export
fitTrendMoMo <- function(model, data){
  mf <- trendModelFit(model, data)
  cf <- getCoefficients(mf)
  ca <- applyIdentifiabilityConstraints(cf, mf)
  out <- list(modelFit = mf, coef = cf, coef_adj = ca)
  class(out) <- "trendMoMoFit"
  out
}

#' Fit–extract–constrain (projection constraints)
#' @inheritParams trendModelFit
#' @return A list with \code{modelFit}, raw \code{coef}, and adjusted \code{coef_adj}.
#' @export
fitTrendMoMo2 <- function(model, data){
  mf <- trendModelFit(model, data)
  cf <- getCoefficients(mf)
  ca <- applyIdentifiabilityConstraints2(cf, mf)
  out <- list(modelFit = mf, coef = cf, coef_adj = ca)
  class(out) <- "trendMoMoFit"
  out
}
