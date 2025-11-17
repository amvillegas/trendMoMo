#' Trend-MoMo model constructor (prototype)
#'
#' Builds a StMoMo-style mortality model that allows for:
#' (i) an age-specific linear improvement term \eqn{d_x (t - \bar{t})}, and
#' (ii) one or more **trend breaks** at specified calendar years
#' through terms \eqn{d^{(i)}_x (t - t_i)_+}. Period and cohort factors are
#' inherited from \pkg{StMoMo} via \code{periodAgeFun} and \code{cohortAgeFun}.
#'
#' @details
#' The general formulation implemented in this model is:
#' \deqn{
#' \log m_{x,t} = \alpha_x
#' + \kappa^{(1)}_t
#' + \sum_{j=1}^M \beta^{(j)}_x \kappa^{(j)}_t
#' + d_x (t - \bar{t})
#' + \sum_{i=1}^B d^{(i)}_x (t - t_i)_+
#' + \gamma_{t-x},
#' }
#' where:
#' \itemize{
#'   \item \eqn{\alpha_x} is the age-specific intercept;
#'   \item \eqn{\kappa^{(1)}_t} and \eqn{\beta^{(j)}_x \kappa^{(j)}_t} are
#'   standard period terms from \pkg{StMoMo};
#'   \item \eqn{d_x (t - \bar{t})} is a baseline linear improvement term;
#'   \item \eqn{d^{(i)}_x (t - t_i)_+} are post-break adjustments for each trend change
#'   at calendar year \eqn{t_i};
#'   \item \eqn{\gamma_{t-x}} represents the cohort effect.
#' }
#'
#' Identifiability is subsequently enforced via
#' \code{\link{applyIdentifiabilityConstraints}} or
#' \code{\link{applyIdentifiabilityConstraints2}},
#' following the transformations in the referenced note
#' (e.g. zero mean, zero linear trend, and orthogonality to broken-line components).
#'
#' **Note:** At this stage, the implementation has only been tested for
#' the case of \eqn{\kappa^{(1)}_t} — that is, for *trend-break extensions
#' of the APCI model*.  Other combinations of period and cohort factors
#'
#' remain experimental and unverified.#' The \code{trendAgeFun} argument controls whether each break’s slope change is
#' age-specific (\code{"NP"}) or age-constant (\code{"1"}).
#'
#' @param link Link function passed to \code{\link[StMoMo]{StMoMo}} (default \code{"log"}).
#' @param staticAgeFun Logical; include an \eqn{\alpha_x} term (default \code{TRUE}).
#' @param periodAgeFun Period-age interaction as in \pkg{StMoMo}; e.g. \code{"NP"} or \code{"1"}.
#' @param cohortAgeFun Cohort-age interaction (or \code{NULL} for no cohort term).
#' @param trend Logical; include a baseline linear improvement term \eqn{d_x (t - \bar{t})}? (default \code{FALSE}).
#' @param trendChanges Numeric vector of calendar years \eqn{t_i} at which trend breaks occur (default \code{NULL}).
#' @param trendAgeFun Character vector of length 1 or \code{length(trendChanges)} with entries
#'   \code{"NP"} (age-varying \eqn{d^{(i)}_x}) or \code{"1"} (age-constant \eqn{d^{(i)}}).
#' @param constFun Optional function to post-process model parameters (for compatibility with \pkg{StMoMo} constructors).
#'
#' @return An object of classes \code{"trendMoMo"} and \code{"StMoMo"} containing:
#'   \item{gnmFormula}{a \code{glm}-style formula with terms for age, period, cohort, baseline trend, and broken trends;}
#'   \item{textFormula}{a textual representation of the model;}
#'   \item{meta}{list elements \code{trend}, \code{trendChanges}, \code{trendAgeFun}, \code{N}, \code{M}, etc.}
#'
#' @references
#' Villegas, A. M., & Arik, Ayse. (2025).
#' *Further Properties of APC Models with Trend Breaks.*
#' UNSW Sydney, School of Risk & Actuarial Studies. :contentReference[oaicite:1]{index=1}
#'
#' @seealso \code{\link[StMoMo]{StMoMo}},
#' \code{\link{trendModelFit}},
#' \code{\link{applyIdentifiabilityConstraints}},
#' \code{\link{applyIdentifiabilityConstraints2}}
#'
#' @examples
#' mod <- trendMoMo(
#'   trend = TRUE,
#'   trendChanges = c(1976, 1999, 2012),
#'   trendAgeFun = c("NP", "1", "NP"),
#'   cohortAgeFun = "NP"
#' )
#' @export

trendMoMo <- function(link = "log", staticAgeFun = TRUE, periodAgeFun = "NP",
                      cohortAgeFun = NULL,
                      trend = FALSE,
                      trendChanges = NULL,
                      trendAgeFun = c("NP", "1"),
                      constFun = function(ax, bx, kt, b0x, gc, wxt, ages)
                        list(ax = ax, bx = bx, kt = kt, b0x = b0x, gc = gc)) {

  #Create estimation model
  model <- StMoMo::StMoMo(link = link,
                          staticAgeFun = staticAgeFun,
                          periodAgeFun = periodAgeFun, cohortAgeFun = cohortAgeFun)
  textFormula <- model$textFormula
  #Add constant improvement rates if necessary
  if (trend){
    model$gnmFormula <- paste(model$gnmFormula, "factor(x):t", sep = " + ")
    textFormula <- paste(textFormula, " + d[x]t", sep = "")
  }
  #Add trend changes if necessary
  M <- 0
  if (!is.null(trendChanges)) {
    M <- length(trendChanges)
    trendAgeFun <- match.arg(trendAgeFun, c("NP", "1"), several.ok = TRUE)
    if (length(trendAgeFun) != M & length(trendAgeFun) != 1) {
      stop("trendAgeFun must be of length 1 or equal to the number of trend changes")
    }
    if (length(trendAgeFun) == 1) {
      trendAgeFun <- rep(trendAgeFun, M)
    }
    for (j in 1:M) {
      if (trendAgeFun[j] == "NP") {
        model$gnmFormula <- paste(model$gnmFormula, " + factor(x):t", j, sep = "")
        textFormula <- paste(textFormula, " + d", j, "[x]max(t-t",j,", 0)", sep = "")
      } else if (trendAgeFun[j] == "1") {
        model$gnmFormula <- paste(model$gnmFormula, " + t", j, sep = "")
        textFormula <- paste(textFormula, " + d", j, "max(t-t",j,", 0)", sep = "")
      }

    }

  }
  out <- list(link = link,
              staticAgeFun = staticAgeFun,
              periodAgeFun = model$periodAgeFun,
              cohortAgeFun = model$cohortAgeFun,
              trend = trend,
              trendChanges = trendChanges,
              trendAgeFun = trendAgeFun,
              N = model$N,
              M = M,
              textFormula = textFormula,
              constFun = constFun,
              gnmFormula = model$gnmFormula)


  #Estimation type
  class(out) <-  c("trendMoMo", "StMoMo")
  out
}
