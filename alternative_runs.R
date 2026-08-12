
setwd("C:/USE/GitHub/Megrim78_2026")

## libraries
# library(remotes)
# install_github("flr/FLa4a")
library(FLa4a)
library(a4adiags)
library(FLCore)
library(tidyr)
library(ggplot2)
library(dplyr)
library(icesAdvice)
library(gridExtra)
library(openxlsx)

set.seed(1234)

## Functions ----
reduce_pg_index <- function(ix, pg = 7) {
  
  ages <- as.numeric(dimnames(ix@index)$age)
  maxage <- max(ages)
  
  # 1. Sumar edades pg:maxage en la edad pg
  ix@index[ac(pg), ] <- apply(ix@index[ac(pg:maxage), , drop = FALSE], 2, sum, na.rm = TRUE)
  
  # 2. Recortar edades
  ix <- trim(ix, age = 1:pg)
  
  # 3. Actualizar plusgroup
  ix@range["plusgroup"] <- pg
  
  return(ix)
}

source('./docs/retro_analysis_f.R') # retrospective function

## load the data
load("Input/bootstrap/data/stock/meg78_stock_sop_BE_corrected.RData")
load("Input/IGFS_EVHOE_index/index_sep.RData")


## prepare the data  ----

###### Put index weight as stock.wt

df <- read.xlsx("C:/Users/maltuna/OneDrive - AZTI/Oilarra/WGBIE/WGBIE_2026/0.Original_country_data/Accession&email/Surveys/Weight/Hans/MegIndexMeanLen.xlsx")

df2 <- df %>%
  mutate(MeanWeight = MeanWeight / 1000) # Pasar a kg

mean_0305 <- df2 %>%
  filter(Year %in% c(2003, 2004, 2005)) %>% # Media por edad de los tres primeros años para los datos históricos (<2003)
  group_by(Age) %>%
  summarise(MeanWeight = mean(MeanWeight, na.rm = TRUE),
            .groups = "drop")

hist_wt <- expand.grid(
  Age = unique(mean_0305$Age),
  Year = 1984:2002 # Expandir 1984-2002 usando dicha media index data starts in 2003
) %>%
  left_join(mean_0305, by = "Age")

wt_all <- bind_rows(
  hist_wt,
  df2 %>% select(Age, Year, MeanWeight) # Juntar históricos + observados
) %>%
  arrange(Age, Year)

wt.mat <- wt_all %>%
  pivot_wider(names_from = Year, values_from = MeanWeight) %>%
  arrange(Age) # Matriz edad x año

ages <- wt.mat$Age
years <- names(wt.mat)[-1]

wt.mat <- as.matrix(wt.mat[, -1])

arr <- array(
  wt.mat,
  dim = c(
    length(ages),
    length(years),
    1, 1, 1, 1
  ),
  dimnames = list(
    age    = as.character(ages),
    year   = as.character(years),
    unit   = "unique",
    season = "all",
    area   = "unique",
    iter   = "1"
  )
)

stock.wt <- FLQuant(arr)
units(stock.wt) <- "kg"
stock@stock.wt <- stock.wt

###### Put plus group at age 7
stk7 <- setPlusGroup(stock, 7)
stk7@catch.n['1',as.character(1984:2000)] <- NA # We do not really believe that the increase in 1-year-olds in the catch is real so we shouldn't formulate a model that treats it as real.

idx7 <- index_sep
idx7[[1]] <- reduce_pg_index(idx7[[1]], 7)
idx7[[2]] <- reduce_pg_index(idx7[[2]], 7)
idx7[[3]] <- reduce_pg_index(idx7[[3]], 7)
idx7[[3]]@index[,"2017"] <- NA

idx7[[1]]@index[ac(1:3),ac(2015:2021)] <- NA # Porcupine
idx7[[2]]@index[ac(1:3),ac(2015:2021)] <- NA # IGFS
idx7[[3]]@index[ac(1:3),ac(2015:2021)] <- NA # EVHOE

## RUNS (Ernesto) ----
fit00 <- sca(stk7, idx7)
res00 <- residuals(fit00, stk7, idx7)
plot(res00)
plot(res00, by = "age")

fmod <- ~te(age, year, k = c(5, 10), bs = "tp", by=as.numeric(year>2000)) + s(age, k = 5) + s(year, k=10) + s(year, k=5, by=as.numeric(age==7))
srmod <- ~I(as.numeric(year<=1996)) + s(year, k=10, by=as.numeric(year>=1997))
qmod <- list(~factor(age),~factor(age), ~factor(age))
fit01 <- sca(stk7, idx7, fmodel=fmod, srmodel=srmod, qmodel=qmod)
res01 <- residuals(fit01, stk7, idx7)
plot(res01)

cthDg01 <- computeCatchDiagnostics(fit01, stk7)
plot(cthDg01)
plot(cthDg01, type="prediction", probs=c(0.025, 0.975))

n <- 4
# list to hold data for retrospective fits
nret <- as.list(1:n)
stks <- FLStocks(lapply(nret, function(x){window(stk7, end=(range(stk7)["maxyear"]-x))}))
idxs <- lapply(nret, function(x){window(idx7, end=(range(idx7)["maxyear"]-x))})
# fit to each list element, note scas can be paralelized
fits01 <- scas(stks, idxs, fmodel=list(fmod), srmodel=list(srmod), qmodel=list(qmod), workers=n)
# update stock object with fit
stks <- stks + fits01
# add candidate fit
stks[[5]] <- stk7 + simulate(fit01, 250)
plot(window(stks, start=2000)) + theme(legend.position = "none") + scale_colour_manual(values = rep("black", n+1))

plot(stk7 + simulate(fit01, 250))


## RUNS (Miren) ----

# option 1
srmod <- ~ bevholt(CV = 0.3)
fmod <- ~s(age, k = 3, by = breakpts(year, 2013)) + factor(year)
qmod <- list(~factor(age),~factor(age), ~factor(age))
n1mod <- ~s(age, k = 3) 
vmod <- list(~s(age, k = 3), ~1, ~1, ~1) 

# Alternative options for srmodel
# # option 2 (all the index data available from 2003 onwards)
# srmod <- ~ factor(ifelse(year < 2003, "pre2003", as.character(year)))
# 
# # option 3 (one index data available from 1997 onwards)
# srmod <- ~ factor(ifelse(year < 1997, "pre1998", as.character(year)))
#
# # option 4 (no use factor use spline)
# srmod <- ~ s(year, k=10, by=as.numeric(year>2003))
#
# # option 5
# srmod <- ~I(as.numeric(year<=1996)) + s(year, k=10, by=as.numeric(year>=1997))
# 
# Alternative options for fmodel
# fmod <- ~ti(age, year, k = c(3, 25))
# fmod <- ~te(age, year, k = c(5, 26))

# RUN
fit01 <- sca(stk7, idx7, fmodel = fmod, qmodel = qmod, srmodel = srmod, vmodel = vmod, n1model = n1mod)
stk01 <- stk7 + fit01

# Residuals
res01 <- residuals(fit01, stk7, idx7)
plot(res01)

# Retro analysis plot with monrho values
results <- run_retro_analysis(stk7, idx7, fit01, fmod, qmod, srmod, vmod, n1mod)
results$rho_table <- results$rho_table %>% mutate(x = 2025, y = 0)
results$rho_table$qname <- c("F" = "F", "SSB" = "SB", "Recruitment" = "Rec", "Catch" = "C")

new_names <- c("Rec" = "Recruitment", "SB" = "SSB", "C" = "Catch", "F" = "F")
plot(FLStocks(results$retro), col = 1, lwd = 1) +
  facet_wrap(~qname, scales = 'free_y', labeller = labeller(qname = new_names)) +
  geom_text(
    data = results$rho_table,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust = 1, vjust = 0) +
  theme_bw() +
  labs(color = "N years removed")

### Selectivity and catchability plot
fitted <- predict(pars(fit01))

a  <- xyplot(data~age,groups=year,stk01@harvest,type='b',ylim=c(0,1),ylab='F',main='Fishing mortality')
a1 <- xyplot(data~age,groups=year,data=fitted$qmodel[1],type='b',ylab='Catchability',main="Porcupine")
a2 <- xyplot(data~age,groups=year,data=fitted$qmodel[2],type='b',ylab='Catchability',main="IGFS")
a3 <- xyplot(data~age,groups=year,data=fitted$qmodel[3],type='b',ylab='Catchability',main="EVHOE")
grid.arrange(a,a1,a2, a3, ncol=2)
