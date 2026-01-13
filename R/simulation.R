#' Function to simulat d^{(i)}_x terms
#'
#' This funcion simulates future values of the age-specific trend deviation terms
#' d^{(i)}_x for a Trend-MoMo model with breakpoints
#'
#' @details
#' The d^{(i)}_x terms represent age-specific deviations from the linear trend
#' after each breakpoint t_i. This function simulates future values of these terms
#' based on a given model specification. In the simulation it assumes that the
#' future d^{(i)}_x terms do not depend on age, i.e., they are constant across ages.
#'
#' @param fittedModel The object returned by \code{\link{fitTrendMoMo}} or \code{\link{fitTrendMoMo2}}.
#' @param h Numeric; the number of years ahead to simulate.
#' @param nSim Numeric; the number of simulation paths to generate.
#' @param method_size  Character; the simulation method to use. Options are:
#' \code{"normal"} for normal distribution, \code{"AR1"} for autoregressive model, possibly
#' with box cox transformation, more to be implemented.
#' @param method_size_params A list of additional parameters for the size simulation method.
#' For \code{method_size = "normal"}, this can include:
#' \itemize{
#'   \item \code{mean}: Numeric; the mean of the normal distribution (default 0).
#'   \item \code{sd}: Numeric; the standard deviation of the normal distribution (default 0.01).
#' }
#' For \code{method_size = "AR1"}, this can include:
#' \itemize{
#'   \item \code{phi}: Numeric; the AR(1) coefficient (default 0.5).
#'   \item \code{sd}: Numeric; the standard deviation of the white noise (default 0.01).
#'   \item \code{init}: Numeric; the initial value of the AR(1) process (default 0).
#'   \item \code{lambda} : Numeric; the box cox transformation parameter (default NULL, no transformation).
#' }
#' @param method_duration Character; the simulation method for duration. Options are:
#' \code{"exponential"} for exponential distribution, \code{"gamma"} for gamma distribution,
#'  more to be implemented.
#' @param method_duration_params A list of additional parameters for the duration simulation method.
#' For \code{method_duration = "exponential"}, this can include:
#' \itemize{
#'   \item \code{mean}: Numeric; the mean number of year between trend changes (default 10).
#' }
#' For \code{method_duration = "gamma"}, this can include:
#' \itemize{
#'   \item \code{shape}: Numeric; the shape parameter of the gamma distribution (default 2).
#'   \item \code{rate}: Numeric; the rate parameter of the gamma distribution (default 0.2).
#' }
#' @return A list where each element corresponds to a simulation of the d^{(i)}_x terms. It has two
#' #' components:
#' \item{dx_list}{A list of length equal to the number of breakpoints in the model. Each element is a numeric
#' vector representing the simulated d^{(i)}_x values for that breakpoint.}
#' \item{breakpoints}{A numeric vector of the breakpoint years in the model.}
#'
#' @export
simulateTrendDeviations <- function(fittedModel, h = 20, nSim = 1,
                                    method_size = c("normal", "AR1"),
                                    method_size_params = list(),
                                    method_duration = c("exponential", "gamma"), method_duration_params = list(),
                                    ...) {
  stopifnot(inherits(fittedModel, "trendMoMoFit"))
  method_size <- match.arg(method_size)
  method_duration <- match.arg(method_duration)

  # Extract model information
  trendChanges <- fittedModel$modelFit$model$trendChanges #Current times of breakpoints
  dx_list <- fittedModel$coef_adj$dx_list #Current d^{(i)}_x terms
  nBreaks <- length(trendChanges)
  ages <- fittedModel$modelFit$ages.fit
  nAges <- length(ages)

  #Simulate new times of breakpoints
  simulated_breakpoints <- list()
  for (sim in 1:nSim) {
    sim_breakpoints <- c()
    if (nBreaks > 0) {
      current_time <- max(fittedModel$modelFit$years.fit)
      while (current_time < max(fittedModel$modelFit$years.fit) + h) {
        #Simulate duration until next breakpoint
        if (method_duration == "exponential") {
          mean_duration <- method_duration_params$mean
          if (is.null(mean_duration)) {
            mean_duration <- 10
          }
          duration <- ceiling(rexp(1, rate = 1/mean_duration))
        } else if (method_duration == "gamma") {
          shape_param <- method_duration_params$shape
          rate_param <- method_duration_params$rate
          if (is.null(shape_param)) {
            shape_param <- 2
          }
          if (is.null(rate_param)) {
            rate_param <- 0.2
          }
          duration <- ceiling(rgamma(1, shape = shape_param, rate = rate_param))
        }

        else {
          stop("Unsupported method_duration")
        }
        current_time <- current_time + duration
        if (current_time < max(fittedModel$modelFit$years.fit) + h) {
          sim_breakpoints <- c(sim_breakpoints, current_time)
        }
      }
    }
    simulated_breakpoints[[sim]] <- sim_breakpoints
  }
  #Simulate d^{(i)}_x terms for each simulation

  #get lambda for box cox transformation if provided
  lambda <- NULL
  if(!is.null(method_size_params$lambda)) {
    lambda <- method_size_params$lambda
    if (lambda == 0) {
      lambda <- NULL
    }
  }

  simulated_dx_list <- list()
  for (sim in 1:nSim) {
    #check if there are any breakpoints
    if (length(simulated_breakpoints) < sim) {
      sim_breakpoints <- c()
    } else {
    sim_breakpoints <- simulated_breakpoints[[sim]]
    }

    sim_dx_list <- list()
    #Initial current size for AR(1) based on init parameter
    size_1 <- NULL
    if (method_size == "AR1") {
      size_1 <- method_size_params$init
      if (is.null(size_1)) {
        size_1 <- 0
      }
      size_1_tilde <- size_1
      if (!is.null(lambda)) {
        size_1_tilde <-(exp(size_1 * lambda) - 1)/lambda
      }
    }
    for (i in seq_along(sim_breakpoints)) {
      #Simulate size of d^{(i)}_x
      if (method_size == "normal") {
        mean_size <- method_size_params$mean
        sd_size <- method_size_params$sd
        if (is.null(mean_size)) {
          mean_size <- 0
        }
        if (is.null(sd_size)) {
          sd_size <- 0.01
        }
        size <- rnorm(1, mean = mean_size, sd = sd_size)
      } else if (method_size == "AR1") {
        phi <- method_size_params$phi
        sd_size <- method_size_params$sd
        init <- method_size_params$init
        if (is.null(phi)) {
          phi <- 0.5
        }
        if (is.null(sd_size)) {
          sd_size <- 0.01
        }
        #Simulate one step of AR(1) depending on box-cox transformation
        size_tilde <- phi * size_1_tilde + rnorm(1, mean = 0, sd = sd_size)
        size_1_tilde <- size_tilde
        if (is.null(lambda)){
          size <- size_tilde

        } else {
          size <- log(lambda * size_tilde + 1)/lambda
        }

      } else {
        stop("Unsupported method_size")
      }
      #Assume d^{(i)}_x is constant across ages
      sim_dx_list[[i]] <- rep(size, nAges)
    }
    simulated_dx_list[[sim]] <- sim_dx_list
  }
  return(list(dx_list = simulated_dx_list, breakpoints = simulated_breakpoints))
}


#' Simulate kt1 future paths
#'
#' This function simulates future paths of the period effect kt1 for a Trend-MoMo model.
#'
#' @param fittedModel The object returned by \code{\link{fitTrendMoMo}} or \code{\link{fitTrendMoMo2}}.
#' @param h Numeric; the number of years ahead to simulate.
#' @param nSim Numeric; the number of simulation paths to generate.
#' @param method_kt Character; the simulation method to use. Options are:
#' \code{"rw"} for random walk, \code{"normal"} for normal distribution, more to be implemented.
#' @param method_params_kt A list of additional parameters for the simulation method.
#' For \code{method_kt = "rw"}, this can include:
#' \itemize{
#'   \item \code{sd}: Numeric; the standard deviation of the random walk increments (default estimated from the data).
#'   \item \code{last_kt1}: Numeric; the last observed value of kt1 to start the simulation from (default taken from the fitted model).
#' }
#' For \code{method_kt = "normal"}, this can include:
#' \itemize{
#'   \item \code{mean}: Numeric; the mean of the normal distribution (default 0).
#'   \item \code{sd}: Numeric; the standard deviation of the normal distribution (default estimated from the data).
#'   }
#' @return A matrix of dimension \code{h x nSim} where each column corresponds to a simulation of the kt1 path.
#' @export
simulateKt1 <- function(fittedModel, h = 20, nSim = 1, method_kt = c("normal", "rw"),
                        method_params_kt = list(), ...) {
  stopifnot(inherits(fittedModel, "trendMoMoFit"))
  # Extract model information
  kt1_current <- fittedModel$coef_adj$kt1
  last_kt1 <- as.numeric(tail(kt1_current, n = 1))
  simulated_kt1 <- matrix(0, nrow = h, ncol = nSim)

  #Default parameters if not provided depending on method
  if (method_kt == "rw") {
    #Estimate sd of increments from historical kt1
    kt1_diff <- diff(kt1_current)
    sd_est <- sd(kt1_diff, na.rm = TRUE)
    if (is.na(sd_est) || sd_est == 0) {
      sd_est <- 0.02
    }
    if (is.null(method_params_kt$sd)) {
      method_params_kt$sd <- sd_est
    }
    if (is.null(method_params_kt$last_kt1)) {
      method_params_kt$last_kt1 <- last_kt1
    }
  }
  if (method_kt == "normal") {
    #Estimate sd from historical kt1
    sd_est <- sd(kt1_current, na.rm = TRUE)
    if (is.na(sd_est) || sd_est == 0) {
      sd_est <- 0.02
    }
    if (is.null(method_params_kt$mean)) {
      method_params_kt$mean <- 0
    }
    if (is.null(method_params_kt$sd)) {
      method_params_kt$sd <- sd_est
    }
  }

  #Simulate kt1 paths
  for (sim in 1:nSim) {
    if (method_kt == "rw") {
      increments <- rnorm(h, mean = 0, sd = method_params_kt$sd)
      kt1_path <- method_params_kt$last_kt1 + cumsum(increments)
      simulated_kt1[, sim] <- kt1_path
    } else if (method_kt == "normal") {
      kt1_path <- rnorm(h, mean = method_params_kt$mean, sd = method_params_kt$sd)
      simulated_kt1[, sim] <- kt1_path
    } else {
      stop("Unsupported method")
    }

  }
  return(simulated_kt1)
}

#' Simulate future mortality rates
#'
#' This function simulates future mortality rates based on a fitted Trend-MoMo model.
#'
#' @param fittedModel The object returned by \code{\link{fitTrendMoMo}} or \code{\link{fitTrendMoMo2}}.
#' @param h Numeric; the number of years ahead to simulate.
#' @param nSim Numeric; the number of simulation paths to generate.
#' @param ... Additional parameters to pass to \code{\link{simulateKt1}} and \code{\link{simulateTrendDeviations}}.
#' @return A list where each element corresponds to a simulation of future mortality rates.
#' Each element is a matrix of dimension \code{length(ages) x h} with ages in rows and years in columns.
#' @export
simulateFutureRates <- function(fittedModel, h = 20, nSim = 1, ...) {
  stopifnot(inherits(fittedModel, "trendMoMoFit"))
  ages <- fittedModel$modelFit$ages.fit
  nAges <- length(ages)
  years_fit <- fittedModel$modelFit$years.fit
  last_year <- max(years_fit)
  future_years <- (last_year + 1):(last_year + h)
  simulated_rates <- list()
  #Simulate kt1 paths

  #check if method_kt is provided in ...
  args_list <- list(...)
  if (!"method_kt" %in% names(args_list)) {
    args_list$method_kt <- "normal"
  }


  simulated_kt1 <- simulateKt1(fittedModel, h = h, nSim = nSim, method_kt = args_list$method_kt,
                              method_params_kt = args_list$method_params_kt)
  #Simulate d^{(i)}_x terms
  simulated_dx_list <- simulateTrendDeviations(fittedModel, h = h, nSim = nSim, ...)
  for (sim in 1:nSim) {
    kt1_path <- simulated_kt1[, sim]
    dx_list_sim <- simulated_dx_list$dx_list[[sim]]
    #Build full dx_list including original breakpoints
    full_dx_list <- list()
    original_breaks <- fittedModel$modelFit$model$trendChanges
    all_breaks <- c(original_breaks, simulated_dx_list$breakpoints[[sim]])
    all_breaks_sorted <- sort(all_breaks)
    for (i in seq_along(all_breaks_sorted)) {
      brk <- all_breaks_sorted[i]
      if (brk %in% original_breaks) {
        idx <- which(original_breaks == brk)
        full_dx_list[[i]] <- fittedModel$coef_adj$dx_list[[idx]]
      } else {
        idx_sim <- which(simulated_dx_list$breakpoints[[sim]] == brk)
        #print all breaks
        # print("All breaks:")
        # print(simulated_dx_list$breakpoints[[sim]])
        # #print sim amd break sim
        # print(paste("Simulation:", sim, "Breakpoint:", brk))
        # print(paste("Index in simulated dx_list:", idx_sim))

        full_dx_list[[i]] <- dx_list_sim[[idx_sim]]
        # print(paste("dx value:", dx_list_sim[[idx_sim]]))
        # print("/n")
      }
    }
    # Get cohort effects for future years
    gc <- fittedModel$coef_adj$gc
    if (!is.null(gc)) {
      cohort_years <- (min(future_years)-max(ages)):(max(future_years)-min(ages))
      gc_future <- rep(0, length(cohort_years))
      names(gc_future) <- as.character(cohort_years)
      #replace gc_future values with gc values where available
      index_gc <- names(gc_future)[names(gc_future) %in% names(gc)]
      gc_future[index_gc] <- gc[index_gc]
    } else {
      gc_future <- NULL
    }
    #Predict future rates
    logRates_future <- predictLinkTrendMoMo(
      ax = fittedModel$coef_adj$ax,
      kt1 = kt1_path,
      dx = fittedModel$coef_adj$dx,
      dx_list = full_dx_list,
      gc = gc_future,
      ages = ages,
      years = future_years,
      trendChanges = all_breaks_sorted,
      years.fit = years_fit
    )
    rates_future <- exp(logRates_future)
    rownames(rates_future) <- ages
    colnames(rates_future) <- future_years
    simulated_rates[[sim]] <- rates_future
  }
  # return rates and simulated kt1 and dx_list
list(rates = simulated_rates,
       kt1 = simulated_kt1,
       dx_list = simulated_dx_list)
}

