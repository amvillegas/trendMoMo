#' Plot adjusted Trend-MoMo coefficients
#' @param fittedModel The object returned by \code{\link{fitTrendMoMo}} or \code{\link{fitTrendMoMo2}}.
#' @export
plotTrendMoMo <- function(fittedModel){
  coef <- fittedModel$coef
  coef_adj <- fittedModel$coef_adj
  modelFit <- fittedModel$modelFit
  years.fit <- modelFit$years.fit
  ages.fit <- modelFit$ages.fit
  tbar <- mean(years.fit)

  #Check if there are cohort effects to plot to know if we need to adjust the par(mfrow)
  n_plots <- 0
  if (!is.null(coef$ax)) n_plots <- n_plots + 1
  if (!is.null(coef$dx)) n_plots <- n_plots + 1
  if (!is.null(coef$kt1)) n_plots <- n_plots + 1
  if (length(coef$dx_list) > 0) n_plots <- n_plots + 1
  if (!is.null(coef$gc)) n_plots <- n_plots + 1
  if (n_plots == 0) {
    stop("No coefficients to plot")
  } else if (n_plots == 1) {
    par(mfrow = c(1,1))
  } else if (n_plots == 2) {
    par(mfrow = c(1,2))
  } else if (n_plots <= 4) {
    par(mfrow = c(2,2))
  } else {
    par(mfrow = c(3,2))
  }


  if (!is.null(coef$ax)) {
    plot(ages.fit, coef_adj$ax, type = "l", main = "ax", xlab = "Age", ylab = "ax")
  }
  if (!is.null(coef$dx)) {
    plot(ages.fit, coef_adj$dx, type = "l", main = "dx", xlab = "Age", ylab = "dx",
         ylim = range(c(0, coef_adj$dx)))
    #add dashed horizontal line at y=0
    abline(h = 0, lty = 2)
  }

  #plot all dx_list in one subplot
  if (length(coef$dx_list) > 0) {
    matplot(ages.fit, do.call(cbind, coef_adj$dx_list), type = "l", main = "dx(i)", xlab = "Age", ylab = "dx_list")
    #add dashed horizontal line at y=0
    abline(h = 0, lty = 2)
    #Remove box around legend
    legend("topright", legend = paste0("dx", seq_along(coef$dx_list)), col = 1:length(coef$dx_list), lty = 1, bty = "n")
  }
  if (!is.null(coef$kt1)) {
    plot(years.fit, coef_adj$kt1, type = "l", main = "kt1", xlab = "Year", ylab = "kt1")
    #add dashed horizontal line at y=0
    abline(h = 0, lty = 2)
  }
  if (!is.null(coef$gc)) {
    plot(as.numeric(names(coef_adj$gc)), coef_adj$gc, type = "l", main = "gc", xlab = "Cohort", ylab = "gc")
    #add dashed horizontal line at y=0
    abline(h = 0, lty = 2)
  }

  par(mfrow = c(1,1))

}
