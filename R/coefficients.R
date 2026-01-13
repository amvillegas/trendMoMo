#' Extract structured coefficients from a fitted Trend-MoMo
#' @param modelFit The object returned by \code{\link{trendModelFit}}.
#' @return A list with elements \code{ax}, \code{kt1}, \code{dx}, \code{dx_list}, \code{gc}.
#' @export
getCoefficients <- function(modelFit){
  modelCoef <- modelFit$fit$coefficients

  #ax
  ax <- NULL
  if (modelFit$model$staticAgeFun){
    ax <- modelCoef[grep(pattern = "^factor[(]x[)][[:digit:]]+$",
                         names(modelCoef))]
    names(ax) <- modelFit$ages.fit
  }

  #kt1
  kt1 <- NULL
  if (!is.null(modelFit$model$periodAgeFun)) {
    kt1 <- modelCoef[grep(pattern = "^factor[(]t[)][-]?[[:digit:]]+$",
                          names(modelCoef))]
    #Repleace "factor(t)1981" "factor(t)1982" "factor(t)1983" ... with "1981" "1982" "1983"
    names(kt1) <- sub("factor\\(t\\)", "", names(kt1))
    #check if all years in ages.fit are present in names(kt1)
    missing_years <- setdiff(as.character(modelFit$years.fit), names(kt1))
    if (length(missing_years) > 0) {
      kt1 <- c(kt1, setNames(rep(0, length(missing_years)), missing_years))
    }
    #reorder kt1 according to years.fit
    kt1 <- kt1[as.character(modelFit$years.fit)]
    #Check if there are NAs in kt1 and replace with 0
    kt1[is.na(kt1)] <- 0
  }
  #Trend term
  dx <- NULL
  if (modelFit$model$trend) {
    dx <- modelCoef[grep(pattern = ":t$", names(modelCoef))]
    names(dx) <- modelFit$ages.fit
    #check if all ages in ages.fit are present in names(dx)
    missing_ages <- setdiff(as.character(modelFit$ages.fit), names(dx))
    if (length(missing_ages) > 0) {
      dx <- c(dx, setNames(rep(0, length(missing_ages)), missing_ages))
    }
    #reorder dx according to ages.fit
    dx <- dx[as.character(modelFit$ages.fit)]
    #Check if there are NAs in dx and replace with 0
    dx[is.na(dx)] <- 0
  }

  #Trend change terms
  dx_list <- list()
  if (!is.null(modelFit$model$trendChanges)) {
    for (i in seq_along(modelFit$model$trendChanges)) {
      change_year <- modelFit$model$trendChanges[i]
      #case 1: trendAgeFun = "1"
      if (modelFit$model$trendAgeFun[i] == "1") {

        term_pattern <- paste0("t", i, "$")
        dx_i <- modelCoef[grep(pattern = term_pattern, names(modelCoef))]
        #In this case it should be only one coefficient. Check if length(dx_i) != 1
        if (length(dx_i) != 1) {
          stop(paste("Error: Expected one coefficient for trend change", i,
                     "but found", length(dx_i)))
        }
        dx_i <- rep(dx_i, length(modelFit$ages.fit))
        names(dx_i) <- modelFit$ages.fit
        #Replace NAs with 0
        dx_i[is.na(dx_i)] <- 0
        dx_list[[paste0("dx", i)]] <- dx_i
      } else if (modelFit$model$trendAgeFun[i] == "NP") {
        #case 2: trendAgeFun = "NP"
        term_pattern <- paste0("factor[(]x[)][[:digit:]]+:t", i, "$")
        dx_i <- modelCoef[grep(pattern = term_pattern, names(modelCoef))]
        names(dx_i) <- sub("factor\\(x\\)", "", names(dx_i))
        # Replace 20:t1"  "21:t1"  "22:t1"  "23:t1"  "24:t1" , ... with "20" "21" "22" "23" "24"
        names(dx_i) <- sub(paste0(":t", i), "", names(dx_i))

        #check if all ages in ages.fit are present in names(dx_i)
        missing_ages <- setdiff(as.character(modelFit$ages.fit), names(dx_i))
        if (length(missing_ages) > 0) {
          dx_i <- c(dx_i, setNames(rep(0, length(missing_ages)), missing_ages))
        }
        #reorder dx_i according to ages.fit
        dx_i <- dx_i[as.character(modelFit$ages.fit)]
        #Check if there are NAs in dx_i and replace with 0
        dx_i[is.na(dx_i)] <- 0
        dx_list[[paste0("dx", i)]] <- dx_i
      }
    }
  }
  #cohort effect
  gc <- NULL
  if (!is.null(modelFit$model$cohortAgeFun)) {
    #Get the cohort effect coefficients: They are of the form factor(c)XXXX
    gc <- modelCoef[grep(pattern = "^factor[(]c[)][-]?[[:digit:]]+$",
                         names(modelCoef))]
    names(gc) <- sub("factor\\(c\\)", "", names(gc))
    #check if all cohorts in data are present in names(gc): cohorts are all
    # possible values of t - x with t in modelFit$ages.fit and x in modelFit$years.fit

    cohorts <- sort(unique(modelFit$fit$data$t - modelFit$fit$data$x))

    missing_cohorts <- setdiff(as.character(cohorts), names(gc))
    if (length(missing_cohorts) > 0) {
      gc <- c(gc, setNames(rep(0, length(missing_cohorts)), missing_cohorts))
    }
    #reorder gc according to cohorts
    gc <- gc[as.character(cohorts)]
    #Check if there are NAs in gc and replace with 0
    gc[is.na(gc)] <- 0
  }
  return(list(ax = ax, kt1 = kt1, dx = dx, dx_list = dx_list, gc = gc))
}
