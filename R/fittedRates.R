#' Function to extract fitted rates from a fitted trend MoMo model
#'
#' This function extracts the fitted rates from a fitted trend MoMo model object.
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
#' @param fittedModel A fitted trend MoMo model object.
#' @return A matrix of fitted rates with ages in rows and years in columns.
#' @export
fittedRatesTrendMoMo <- function (fittedModel)
{

  stopifnot(inherits(fittedModel, "trendMoMoFit"))

  #predict rates using predictLinkTrendMoMo

  logRates <- predictLinkTrendMoMo(
    ax = fittedModel$coef_adj$ax,
    kt1 = fittedModel$coef_adj$kt1,
    dx = fittedModel$coef_adj$dx,
    dx_list = fittedModel$coef_adj$dx_list,
    gc = fittedModel$coef_adj$gc,
    ages = fittedModel$modelFit$ages.fit,
    years = fittedModel$modelFit$years.fit,
    trendChanges = fittedModel$modelFit$model$trendChanges
  )
  fitted_rates <- exp(logRates)
  rownames(fitted_rates) <- fittedModel$modelFit$ages.fit
  colnames(fitted_rates) <- fittedModel$modelFit$years.fit
  return(fitted_rates)
  #
  # coef_adj <- fittedModel$coef_adj
  # modelFit <- fittedModel$modelFit
  # years.fit <- modelFit$years.fit
  # ages.fit <- modelFit$ages.fit
  # tbar <- mean(years.fit)
  #
  # matrix_ax <- matrix(coef_adj$ax, nrow = length(ages.fit), ncol = length(years.fit),
  #                     byrow = FALSE)
  #
  # matrix_dx <- matrix(coef_adj$dx, nrow = length(ages.fit), ncol = length(years.fit),
  #                     byrow = FALSE)
  # matrix_year_effect <- matrix(years.fit, nrow = length(ages.fit), ncol = length(years.fit),
  #                              byrow = TRUE)
  # matrix_dx_term <- matrix_dx * (matrix_year_effect - tbar)
  #
  # #d_x^(i) terms
  # dx_i_term <- matrix(0, nrow = length(ages.fit), ncol = length(years.fit))
  # if (length(coef_adj$dx_list) > 0) {
  #   for (i in 1:length(coef_adj$dx_list)) {
  #     dx_i <- coef_adj$dx_list[[i]]
  #     ti <- fittedModel$modelFit$model$trendChanges[i]
  #     matrix_dx_i <- matrix(dx_i, nrow = length(ages.fit), ncol = length(years.fit),
  #                           byrow = FALSE)
  #     dx_i_term <- dx_i_term + matrix_dx_i * pmax(0, matrix_year_effect - ti)
  #   }
  # }
  # #kt1 term
  # kt1_term <- matrix(0, nrow = length(ages.fit), ncol = length(years.fit))
  # if (!is.null(coef_adj$kt1)) {
  #   matrix_kt1 <- matrix(coef_adj$kt1, nrow = length(ages.fit), ncol = length(years.fit),
  #                           byrow = TRUE)
  #   kt1_term <- matrix_kt1
  # }
  # #gc term (cohort effect)
  # gc_term <- matrix(0, nrow = length(ages.fit), ncol = length(years.fit))
  # nAges <- length(ages.fit)
  # nYears <- length(years.fit)
  # if (!is.null(coef_adj$gc)) {
  #   for (i in 1:nAges) {
  #     for (j in 1:nYears) {
  #       cohort_index <- years.fit[j] - ages.fit[i]
  #                 gc_term[i, j] <- coef_adj$gc[as.character(cohort_index)]
  #     }
  #   }
  #
  # }
  #
  # fitted_log_rates <- matrix_ax + matrix_dx_term + dx_i_term + kt1_term + gc_term
  # fitted_rates <- exp(fitted_log_rates)
  # rownames(fitted_rates) <- ages.fit
  # colnames(fitted_rates) <- years.fit
  # return(fitted_rates)
}

#'Function to predict log rates given a_x, k_t, d_x, d_xi, gc
#'
#'This function predicts rates given the components of the trend MoMo model.
#'@details
#'The general formulation implemented in this model is:
#'\deqn{
#'\log m_{x,t} = \alpha_x + \kappa^{(1)}_t +
#'\sum_{j=1}^M \beta^{(j)}_x \kappa^{(j)}_t +
#'d_x (t - \bar{t}) + \sum_{i=1}^B d^{(i)}_x (t - t_i)_+ + \gamma_{t-x}}
#'@param ax Numeric vector of age-specific intercepts.
#'@param kt1 Numeric vector of period factor.
#'@param dx Numeric vector of baseline trend components.
#'@param dx_list List of numeric vectors of trend break components.
#'@param gc Numeric vector of cohort effects.
#'@param trendChanges Numeric vector of calendar years at which trend breaks occur.
#'@param ages Numeric vector of ages.
#'@param years Numeric vector of years.
#'@param years.fit Optional numeric vector of years to calcuate tbar (if NULL, uses years).
#'@return A matrix of predicted rates with ages in rows and years in columns.
#'@export
predictLinkTrendMoMo <- function(ax, kt1, dx, dx_list, gc, ages, years, trendChanges, years.fit = NULL) {
  #calculate tbar
  if (is.null(years.fit)) {
    years.fit <- years
  }
  tbar <- mean(years.fit)
  nAges <- length(ages)
  nYears <- length(years)

  #check lengths for ax only it is not NULL
  if (!is.null(ax) & length(ax) != nAges) {
    stop("Length of ax does not match number of ages.")
  }
  #check lengths for kt1 only it is not NULL
  if (!is.null(kt1) & length(kt1) != nYears) {
    stop("Length of kt1 does not match number of years.")
  }
  #check lengths for dx only it is not NULL
  if (!is.null(dx) & length(dx) != nAges) {
    stop("Length of dx does not match number of ages.")
  }
  #check lengths for gc only it is not NULL
  if (!is.null(gc)) {
    cohort_indices <- years - min(ages)
    if (length(gc) < length(unique(cohort_indices))) {
      stop("Length of gc does not match number of cohorts.")
    }
  }

  #ax term if ax is NULL set to 0
  if (is.null(ax)) {
    ax <- rep(0, nAges)
  }
  #dx term if dx is NULL set to 0
  if (is.null(dx)) {
    dx <- rep(0, nAges)
  }
  matrix_ax <- matrix(ax, nrow = length(ages), ncol = length(years),
                      byrow = FALSE)

  matrix_dx <- matrix(dx, nrow = length(ages), ncol = length(years),
                      byrow = FALSE)
  matrix_year_effect <- matrix(years, nrow = length(ages), ncol = length(years),
                               byrow = TRUE)
  matrix_dx_term <- matrix_dx * (matrix_year_effect - tbar)

  #d_x^(i) terms
  dx_i_term <- matrix(0, nrow = length(ages), ncol = length(years))
  if (length(dx_list) > 0) {
    for (i in 1:length(dx_list)) {
      dx_i <- dx_list[[i]]
      ti <- trendChanges[i]
      matrix_dx_i <- matrix(dx_i, nrow = length(ages), ncol = length(years),
                            byrow = FALSE)
      dx_i_term <- dx_i_term + matrix_dx_i * pmax(0, matrix_year_effect - ti)
    }
  }
  #kt1 term
  kt1_term <- matrix(0, nrow = length(ages), ncol = length(years))
  if (!is.null(kt1)) {
    matrix_kt1 <- matrix(kt1, nrow = length(ages), ncol = length(years),
                         byrow = TRUE)
    kt1_term <- matrix_kt1
  }
  #gc term (cohort effect)
  gc_term <- matrix(0, nrow = length(ages), ncol = length(years))

  if (!is.null(gc)) {
    for (i in 1:nAges) {
      for (j in 1:nYears) {
        cohort_index <- years[j] - ages[i]
        gc_term[i, j] <- gc[as.character(cohort_index)]
      }
    }

  }

  predicted_log_rates <- matrix_ax + matrix_dx_term + dx_i_term + kt1_term + gc_term
  return(predicted_log_rates)
}
