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
  coef_adj <- fittedModel$coef_adj
  modelFit <- fittedModel$modelFit
  years.fit <- modelFit$years.fit
  ages.fit <- modelFit$ages.fit
  tbar <- mean(years.fit)

  matrix_ax <- matrix(coef_adj$ax, nrow = length(ages.fit), ncol = length(years.fit),
                      byrow = FALSE)

  matrix_dx <- matrix(coef_adj$dx, nrow = length(ages.fit), ncol = length(years.fit),
                      byrow = FALSE)
  matrix_year_effect <- matrix(years.fit, nrow = length(ages.fit), ncol = length(years.fit),
                               byrow = TRUE)
  matrix_dx_term <- matrix_dx * (matrix_year_effect - tbar)

  #d_x^(i) terms
  dx_i_term <- matrix(0, nrow = length(ages.fit), ncol = length(years.fit))
  if (length(coef_adj$dx_list) > 0) {
    for (i in 1:length(coef_adj$dx_list)) {
      dx_i <- coef_adj$dx_list[[i]]
      ti <- fittedModel$modelFit$model$trendChanges[i]
      matrix_dx_i <- matrix(dx_i, nrow = length(ages.fit), ncol = length(years.fit),
                            byrow = FALSE)
      dx_i_term <- dx_i_term + matrix_dx_i * pmax(0, matrix_year_effect - ti)
    }
  }
  #kt1 term
  kt1_term <- matrix(0, nrow = length(ages.fit), ncol = length(years.fit))
  if (!is.null(coef_adj$kt1)) {
    matrix_kt1 <- matrix(coef_adj$kt1, nrow = length(ages.fit), ncol = length(years.fit),
                            byrow = TRUE)
    kt1_term <- matrix_kt1
  }
  #gc term (cohort effect)
  gc_term <- matrix(0, nrow = length(ages.fit), ncol = length(years.fit))
  nAges <- length(ages.fit)
  nYears <- length(years.fit)
  if (!is.null(coef_adj$gc)) {
    for (i in 1:nAges) {
      for (j in 1:nYears) {
        cohort_index <- years.fit[j] - ages.fit[i]
                  gc_term[i, j] <- coef_adj$gc[as.character(cohort_index)]
      }
    }

  }

  fitted_log_rates <- matrix_ax + matrix_dx_term + dx_i_term + kt1_term + gc_term
  fitted_rates <- exp(fitted_log_rates)
  rownames(fitted_rates) <- ages.fit
  colnames(fitted_rates) <- years.fit
  return(fitted_rates)
}
