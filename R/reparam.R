#' Reparameterise to deviations from a fixed long-term improvement rate (Eq. 11)
#'
#' Transforms the fitted Trend-MoMo parameters so that the baseline age-specific
#' improvement \eqn{d_x} equals a chosen long-term assumption \eqn{\tilde d_x}, and the
#' broken-trend components become **segment deviations** \eqn{\tilde d^{(j)}_x} between
#' consecutive break dates \eqn{t_j}.
#'
#' @details
#' This function implements the reparameterisation described by Equation (11) in the note
#' *Further Properties of APC Models with Trend Breaks*. In terms of mortality improvement
#' rates \eqn{\Delta \log m_{x,t}}, the new parametrisation is:
#'
#' \deqn{
#' \Delta \log m_{x,t}
#' = \tilde d_x
#' + \sum_{j=0}^{M} \tilde d^{(j)}_x \, \mathbf{1}(t_j < t \le t_{j+1})
#' + \Delta \kappa^{(1)}_t
#' + \Delta \gamma_{t-x},
#' }
#'
#' where:
#' \itemize{
#'   \item \eqn{\tilde d_x} is the **long-term** age-specific improvement rate (user-chosen);
#'   \item \eqn{\tilde d^{(j)}_x} is the **deviation** from \eqn{\tilde d_x} on the segment \eqn{(t_j, t_{j+1}]};
#'   \item \eqn{\Delta \kappa^{(1)}_t} captures **short-term shocks**; and
#'   \item \eqn{\Delta \gamma_{t-x}} captures **cohort effects** in improvements.
#' }
#'
#' The deviations are obtained from the fitted (post-constraint) components as:
#' \deqn{
#' \tilde d^{(0)}_x = d_x - \tilde d_x,\qquad
#' \tilde d^{(j)}_x = d_x + \sum_{i=1}^{j} d^{(i)}_x - \tilde d_x,\;\; j=1,\dots,M,
#' }
#' so that each segment’s mean improvement equals \eqn{\tilde d_x + \tilde d^{(j)}_x}.
#' This decomposition cleanly separates (i) **long-term level** \eqn{\tilde d_x},
#' (ii) **temporal deviations** \eqn{\tilde d^{(j)}_x}, (iii) **shocks** \eqn{\Delta\kappa^{(1)}_t},
#' and (iv) **cohort** effects \eqn{\Delta\gamma_{t-x}} for projection purposes.
#'
#' @note
#' This reparameterisation has been developed and tested in the setting with a single period factor
#' \eqn{\kappa^{(1)}_t} (i.e., trend-break extensions of the APCI model).
#'
#' @param fittedModel Output of \code{\link{fitTrendMoMo}} or \code{\link{fitTrendMoMo2}}
#'   containing \code{$coef_adj} = \{a_x, d_x, d^{(j)}_x, \kappa^{(1)}_t, \gamma_{t-x}\} after
#'   identifiability constraints are applied.
#' @param dx_new Numeric vector \eqn{\tilde d_x} (length equal to the number of ages) giving the
#'   **long-term** age-specific improvement assumption to anchor to.
#'
#' @return
#' The input \code{fittedModel} with \code{$coef_adj} replaced by the reparameterised lists:
#' \itemize{
#'   \item \code{dx} set to \eqn{\tilde d_x};
#'   \item \code{dx_list} augmented with \code{dx0} (for \eqn{\tilde d^{(0)}_x}) and updated
#'         \code{dx1, ..., dxM} holding \eqn{\tilde d^{(j)}_x}, \eqn{j=1,\dots,M};
#'   \item other components (e.g., \code{ax}, \code{kt1}, \code{gc}) unchanged.
#' }
#'
#' @examples
#' # After fitting and constraining:
#' #   fit <- fitTrendMoMo(mod, df)
#' # Choose a long-term improvement curve (e.g., flat 1% p.a. across ages)
#' dx_anchor <- rep(0.01, length(fit$modelFit$ages.fit))
#' fit_dev   <- getTrendDeviationParam(fit, dx_anchor)
#' # Now fit_dev$coef_adj$dx is the long-term rate (tilde d_x)
#' # and fit_dev$coef_adj$dx_list$dx0, dx1, ... are segment deviations (tilde d_x^(j))
#'
#' @seealso
#' \code{\link{fitTrendMoMo}}, \code{\link{fitTrendMoMo2}},
#' \code{\link{applyIdentifiabilityConstraints}},
#' \code{\link{applyIdentifiabilityConstraints2}}
#'
#' @references
#' Villegas, A. M., & Arik, Ayse (2025).
#' *Further Properties of APC Models with Trend Breaks* (Equation 11).
#'
#' @aliases getTrendDeviationParam getTrenDeviationParam
#' @export

getTrendDeviationParam <- function(fittedModel, dx_new){
  coef_adj <- fittedModel$coef_adj
  stopifnot(length(dx_new) == length(coef_adj$dx))

  Mbreaks <- length(coef_adj$dx_list)

  coef_new <- coef_adj
  coef_new$dx <- dx_new
  coef_new$dx_list$dx0 <- coef_adj$dx - coef_new$dx

  if (Mbreaks > 0) {
    for (j in seq_len(Mbreaks)) {
      nm <- paste0("dx", j)
      coef_new$dx_list[[nm]] <- coef_adj$dx - coef_new$dx
      for (i in seq_len(j)) {
        coef_new$dx_list[[nm]] <- coef_new$dx_list[[nm]] + coef_adj$dx_list[[paste0("dx", i)]]
      }
    }
  }
  fittedModel$coef_adj <- coef_new
  fittedModel
}
