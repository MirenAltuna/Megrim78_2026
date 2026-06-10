#' Stock input data
#'
#' Input data related to the stock as an FLStock object
#'
#' @name stock
#' @format RData file
#' @tafOriginator WGBIE
#' @tafYear 2026
#' @tafAccess Public
#' @tafSource script

# AI (03/11/2020): Scrit updated by Ane Iriondo for WKTADSA 2020.

library(icesTAF)
library(xlsx)
library(dplyr)
library(FLCore)

setwd("C:/Users/maltuna/OneDrive - AZTI/Oilarra/WGBIE/WGBIE_2026/2.Checks/Input")

end.yr <- 2025

# utility functions
get_sheet <- function(sheet) {
  read.xlsx2(
    taf.boot.path("initial", "data", "Inputdata_MEGRIM_a4a_WGBIE2026.xlsx"),
    sheetName = sheet,
    check.names = FALSE,
    colClasses = numeric()
  )
}

get_age_quant <- function(sheet, units = NA) {
  get_sheet(sheet) %>%
    taf2long(names = c("year", "age", "data")) %>%
    as.FLQuant(units = units)
}

# collate data into an FLStock

# catch.n
catch.n <- get_age_quant("catches.n", "10^3")

# landings.n
landings.n <-
  get_age_quant("landings.n", "10^3") %>%
  window(start = 1984, end = end.yr)


# discards.n
discards.n <-
  get_sheet("discards.n") %>%
  mutate('9' = 0, '10' = 0) %>%
  taf2long(names = c("year", "age", "data")) %>%
  as.FLQuant(units = "10^3") %>%
  window(start = 1984, end= end.yr)


# catches
catch <-
  get_sheet("meg78_caton") %>%
  rename(data = catches) %>%
  select(year, data) %>%
  as.FLQuant(units = "t", quant = "age")

#landings
landings <-
  get_sheet("meg78_caton") %>%
  rename(data = landings) %>%
  select(year, data) %>%
  as.FLQuant(units = "t", quant = "age")

#discards
discards <-
  get_sheet("meg78_caton") %>%
  rename(data = discards) %>%
  select(year, data) %>%
  as.FLQuant(units = "t", quant = "age")

#catches.wt
catch.wt <- get_age_quant("catches.wt", "kg")

#landings.wt
landings.wt <- get_age_quant("landings.wt", "kg")

#discards.wt
discards.wt <- get_age_quant("discards.wt", "kg")


#stock.wt
stock.wt <- get_age_quant("catches.wt", "kg")


#mortality
m <- replace(stock.wt, TRUE, 0.2)
units(m) <- "m"

#maturity
mat <- get_age_quant("mat", "")

#harvest.spwn, harvest before spawning
harvest.spwn <- replace(mat, TRUE, 0)

#m.spwn, mortality before spawning
m.spwn <- replace(mat, TRUE, 0)

# FLSTOCK
stock <-
  FLStock(
    name = "meg78abd", desc = "run1",
    catch.n = catch.n, landings.n = landings.n, discards.n = discards.n,
    catch.wt = catch.wt, landings.wt = landings.wt, discards.wt = discards.wt, stock.wt = stock.wt,
    catch = catch, landings = landings, discards = discards,
    m = m, mat = mat,
    harvest.spwn = harvest.spwn, m.spwn = m.spwn
  )


range(stock, c("minfbar", "maxfbar")) <- c(3, 6)

save(stock, file="bootstrap/data/stock/meg78_stock.RData")
