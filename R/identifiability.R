#' Apply identifiability constraints (sequential reallocation)
#' @param coef Output of \code{\link{getCoefficients}}.
#' @param modelFit The object returned by \code{\link{trendModelFit}}.
#' @return A list mirroring \code{coef} but adjusted.
#' @export
applyIdentifiabilityConstraints <- function(coef, modelFit){
  years.fit <- modelFit$years.fit
  ages.fit <- modelFit$ages.fit
  tbar <- mean(years.fit)

  ax <- coef$ax
  kt1 <- coef$kt1
  dx <- coef$dx
  dx_list <- coef$dx_list
  gc <- coef$gc

  #Set ax to be centered in the average year if both ax and dx are not NULL
  if (!is.null(ax) & !is.null(dx)) {
    ax <- ax + dx * tbar
  }

  # plot(coef$ax + coef$dx*(2020) + coef$kt1["2020"] - (ax + dx*(2020-2001) + kt1["2020"]))

  #Remove quadratic trends from gc if gc, kt, ax, dx are not NULL
  if (!is.null(gc) & !is.null(kt1) & !is.null(ax) & !is.null(dx)) {
    n <- length(gc)
    c <- (min(years.fit)-max(ages.fit)):(max(years.fit)-min(ages.fit))
    phi <- coef(lm(gc ~ c + I(c^2)))
    gc <- gc - phi[1] - phi[2] * c - phi[3] * c^2
    ax <- ax + phi[1] - phi[2]*ages.fit + phi[3]*ages.fit^2 - 2*phi[3]*ages.fit*tbar
    dx <- dx - 2*phi[3]*ages.fit
    kt1 <- kt1 + phi[2]*years.fit + phi[3]*years.fit^2
  }


  #Make sum of kt1 = 0 if kt1 is not NULL and ax is not NULL
  if (!is.null(kt1) & !is.null(ax)) {
    kt1Mean <- mean(kt1)
    kt1 <- kt1 - kt1Mean
    ax <- ax + kt1Mean
  }



  #Remove linear trends from kt1 if kt1 is not NULL and ax and dx are not NULL
  if (!is.null(kt1) & !is.null(ax) & !is.null(dx)) {
    n <- length(kt1)
    tt <- years.fit - tbar
    phi <- coef(lm(kt1 ~ tt))
    kt1 <- kt1 - phi[1] - phi[2] * tt
    ax <- ax + phi[1]
    dx <- dx + phi[2]
  }

  #Remove trends associated wtih trend changes if dx_list is not NULL and ax and dx are not NULL
  # and kt1 is not NULL
  #Do for each of the trend changes
  if (!is.null(dx_list) & !is.null(ax) & !is.null(dx) & !is.null(kt1)) {
    for (i in seq_along(dx_list)) {
      dx_i <- dx_list[[i]]
      if (!is.null(dx_i)) {
        change_year <- modelFit$model$trendChanges[i]
        tt <- pmax(years.fit - change_year, 0)
        phi <- coef(lm(kt1 ~ tt))
        kt1 <- kt1 - phi[1] - phi[2] * tt
        ax <- ax + phi[1]
        dx_list[[i]] <- dx_i + phi[2]
      }
    }
  }
  return(list(ax = ax, kt1 = kt1, dx = dx, dx_list = dx_list, gc = gc))
}

#' Apply identifiability constraints (one-shot projection)
#' @inheritParams applyIdentifiabilityConstraints
#' @return A list of adjusted components \code{ax}, \code{kt1}, \code{dx}, \code{dx_list}, \code{gc}.
#' @export
applyIdentifiabilityConstraints2 <- function(coef, modelFit){

  years.fit <- modelFit$years.fit
  ages.fit <- modelFit$ages.fit
  tbar <- mean(years.fit)

  ax <- coef$ax
  kt1 <- coef$kt1
  dx <- coef$dx
  dx_list <- coef$dx_list
  gc <- coef$gc

  # plot(coef$ax + coef$dx*(2020) + coef$kt1["2020"] - (ax + dx*(2020) + kt1["2020"]))

  #Set ax to be centered in the average year if both ax and dx are not NULL
  if (!is.null(ax) & !is.null(dx)) {
    ax <- ax + dx * tbar
  }

  # plot(coef$ax + coef$dx*(2020) + coef$kt1["2020"] - (ax + dx*(2020-2001) + kt1["2020"]))



  #Remove quadratic trends from gc if gc, kt, ax, dx are not NULL
  if (!is.null(gc) & !is.null(kt1) & !is.null(ax) & !is.null(dx)) {
    n <- length(gc)
    c <- (min(years.fit)-max(ages.fit)):(max(years.fit)-min(ages.fit))
    phi <- coef(lm(gc ~ c + I(c^2)))
    gc <- gc - phi[1] - phi[2] * c - phi[3] * c^2
    ax <- ax + phi[1] - phi[2]*ages.fit + phi[3]*ages.fit^2 - 2*phi[3]*ages.fit*tbar
    dx <- dx - 2*phi[3]*ages.fit
    kt1 <- kt1 + phi[2]*years.fit + phi[3]*years.fit^2
  }



  #Make sum of kt1 = 0 if kt1 is not NULL and ax is not NULL
  if (!is.null(kt1) & !is.null(ax)) {
    kt1Mean <- mean(kt1)
    kt1 <- kt1 - kt1Mean
    ax <- ax + kt1Mean
  }

  # plot(coef$ax + coef$dx*(2020) + coef$kt1["2020"] - (ax + dx*(2020-2001) + kt1["2020"]))


  #Remove linear trends from kt if dx_list, kt, ax, dx are not NULL

  if (length(dx_list) > 0 & !is.null(kt1) & !is.null(ax) & !is.null(dx)){
    #Create design matrix for linear trends and trend breaks
    design_matrix <- matrix(1, nrow = length(years.fit), ncol = 1 + 1 + length(dx_list))
    colnames(design_matrix) <- c("Intercept", "LinearTrend", paste0("Break", seq_along(dx_list)))
    design_matrix[, "LinearTrend"] <- years.fit - tbar
    for (i in seq_along(dx_list)) {
      change_year <- modelFit$model$trendChanges[i]
      design_matrix[, paste0("Break", i)] <- pmax(years.fit - change_year, 0)
    }
    #Fit linear model to kt1

    phi <- coef(lm(kt1 ~ design_matrix - 1)) #-1 to remove intercept
    #Remove linear trend and trend breaks from kt1
    kt1 <- kt1 - design_matrix %*% phi
    #Adjust ax and dx_list accordingly
    ax <- ax + phi["design_matrixIntercept"]
    dx <- dx + phi["design_matrixLinearTrend"]
    for (i in seq_along(dx_list)) {
      dx_list[[i]] <- dx_list[[i]] + phi[paste0("design_matrixBreak", i)]
    }

  } else if (!is.null(kt1) & !is.null(ax) & !is.null(dx)) {
    n <- length(kt1)
    tt <- years.fit - tbar
    phi <- coef(lm(kt1 ~ tt))
    kt1 <- kt1 - phi[1] - phi[2] * tt
    ax <- ax + phi[1]
    dx <- dx + phi[2]
  }

  return(list(ax = ax, kt1 = kt1, dx = dx, dx_list = dx_list, gc = gc))
}


#' Apply identifiability constraints (remove trends up - to a given year)
#' @inheritParams applyIdentifiabilityConstraints
#' @param lastYear The last year up to which to remove trends. If NULL, use all years in \code{modelFit$years.fit}.
#' @return A list of adjusted components \code{ax}, \code{kt1}, \code{dx}, \code{dx_list}, \code{gc}.
#' @export
applyIdentifiabilityConstraintsTrends <- function(coef, modelFit, lastYear = NULL){
  tbar <- mean(years.fit)
  if (is.null(lastYear)) {
    years.fit <- modelFit$years.fit
  } else {
    years.fit <- modelFit$years.fit[modelFit$years.fit <= lastYear]
  }
  ages.fit <- modelFit$ages.fit


  ax <- coef$ax
  kt1 <- coef$kt1
  dx <- coef$dx
  dx_list <- coef$dx_list
  gc <- coef$gc
  #index of years.fit in the original kt1
  year_indices <- match(years.fit, modelFit$years.fit)


  #Make sum of kt1 = 0 if kt1 is not NULL and ax is not NULL
  if (!is.null(kt1) & !is.null(ax)) {
    kt1Mean <- mean(kt1)
    kt1 <- kt1 - kt1Mean
    ax <- ax + kt1Mean
  }

  # plot(coef$ax + coef$dx*(2020) + coef$kt1["2020"] - (ax + dx*(2020-2001) + kt1["2020"]))


  #Remove linear trends from kt if dx_list, kt, ax, dx are not NULL

  if (length(dx_list) > 0 & !is.null(kt1) & !is.null(ax) & !is.null(dx)){
    #Create design matrix for linear trends and trend breaks

    #find number of breaks that are within years.fit
    dx_list_sub <- dx_list[seq_along(modelFit$model$trendChanges)[modelFit$model$trendChanges < max(years.fit)]]

    design_matrix <- matrix(1, nrow = length(years.fit), ncol = 1 + 1 + length(dx_list_sub))
    colnames(design_matrix) <- c("Intercept", "LinearTrend", paste0("Break", seq_along(dx_list_sub)))
    design_matrix[, "LinearTrend"] <- years.fit - tbar
    for (i in seq_along(dx_list_sub)) {
      change_year <- modelFit$model$trendChanges[i]
      design_matrix[, paste0("Break", i)] <- pmax(years.fit - change_year, 0)
    }
    #Fit linear model to kt1



    phi <- coef(lm(kt1[year_indices] ~ design_matrix - 1)) #-1 to remove intercept
    #Remove linear trend and trend breaks from kt1
    kt1[year_indices] <- kt1[year_indices] - design_matrix %*% phi

    #Adjust for years beyond lastYear
    design_matrixLastYear <- matrix(1, nrow = length(modelFit$years.fit) - length(years.fit), ncol = 1 + 1 + length(dx_list_sub))
    colnames(design_matrixLastYear) <- c("Intercept", "LinearTrend", paste0("Break", seq_along(dx_list_sub)))
    design_matrixLastYear[, "LinearTrend"] <- modelFit$years.fit[modelFit$years.fit > max(years.fit)] - tbar
    for (i in seq_along(dx_list_sub)) {
      change_year <- modelFit$model$trendChanges[i]
      design_matrixLastYear[, paste0("Break", i)] <- pmax(modelFit$years.fit[modelFit$years.fit > max(years.fit)] - change_year, 0)
    }
    kt1[-year_indices] <- kt1[-year_indices] - design_matrixLastYear %*% phi

    #Adjust ax and dx_list accordingly
    ax <- ax + phi["design_matrixIntercept"]
    dx <- dx + phi["design_matrixLinearTrend"]
    for (i in seq_along(dx_list_sub)) {
      dx_list[[i]] <- dx_list[[i]] + phi[paste0("design_matrixBreak", i)]
    }

  }

  return(list(ax = ax, kt1 = kt1, dx = dx, dx_list = dx_list, gc = gc))
}
