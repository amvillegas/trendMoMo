# trendMoMo (prototype)

This is a minimal prototype that extends `StMoMo` with:
- a baseline age-specific linear improvement rate,
- multiple broken-trend components at user-specified calendar years,
- identifiability constraint handlers (sequential and one-shot),
- coefficient extraction and plotting helpers.

## Install locally

```r
# in R
install.packages("devtools")
devtools::install_local("trendMoMo")
```

## Minimal example

```r
library(StMoMo)
library(trendMoMo)

mod <- trendMoMo(trend = TRUE,
                 trendChanges = c(1976, 1999, 2012),
                 trendAgeFun = c("NP", "1", "NP"),
                 cohortAgeFun = "NP")

# df with columns: x (age), t (year), D (deaths), E (exposure)
fit1 <- fitTrendMoMo(mod, df)
plotTrendMoMo(fit1)

dx_anchor <- rep(0.01, length(fit1$modelFit$ages.fit))
fit1_dev <- getTrenDeviationParam(fit1, dx_anchor)
```

> Note: This is a research prototype; formulas and constraints are kept simple for clarity.