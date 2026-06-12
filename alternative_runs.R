
setwd("C:/USE/GitHub/Megrim78_2026")

## load the data
load("Input/bootstrap/data/stock/meg78_stock_sop_BE_corrected.RData")
load("Input/IGFS_EVHOE_index/index_sep.RData")

## libraries
# library(remotes)
# install_github("flr/FLa4a")
library(FLa4a)
library(a4adiags)
library(FLCore)
library(ggplot2)
library(dplyr)
library(icesAdvice)
library(gridExtra)

set.seed(1234)

## functions
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

## prepare the data

stk7 <- setPlusGroup(stock, 7)
stk7@catch.n['1',as.character(1984:2000)] <- NA # We do not really believe that the increase in 1-year-olds in the catch is real so we shouldn't formulate a model that treats it as real.

idx7 <- index_sep
idx7[[1]] <- reduce_pg_index(idx7[[1]], 7)
idx7[[2]] <- reduce_pg_index(idx7[[2]], 7)
idx7[[2]]@index[,"2017"] <- NA
idx7[[3]] <- reduce_pg_index(idx7[[3]], 7)
idx7[[3]]@index[,"2017"] <- NA

idx7[[1]]@index[ac(1:3),ac(2015:2021)] <- NA # Porcupine
idx7[[2]]@index[ac(1:3),ac(2015:2021)] <- NA # IGFS
idx7[[3]]@index[ac(1:3),ac(2015:2021)] <- NA # EVHOE
# idx7[["EVHOE"]]@index[,ac(2019:2020)] <- NA # We do not have arguments to remove these years from EVHOE

## RUNS (Ernesto) ----
fit00 <- sca(stk7, idx7)
res00 <- residuals(fit00, stk7, idx7)
plot(res00)
plot(res00, by = "age")

fmod <- ~te(age, year, k = c(5, 10), bs = "tp", by=as.numeric(year>2000)) + s(age, k = 5) + s(year, k=10) + s(year, k=5, by=as.numeric(age==7))
srmod <- ~factor(replace(year, year<1999, 1999))
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
fits01 <- scas(stks, idxs, fmodel=list(fmod), srmodel=list(srmod), workers=n)
# update stock object with fit
stks <- stks + fits01
# add candidate fit
stks[[5]] <- stk7 + simulate(fit01, 250)
plot(window(stks, start=2000)) + theme(legend.position = "none") + scale_colour_manual(values = rep("black", n+1))

plot(stk7 + simulate(fit01, 250))


## RUNS (Miren tries) ----

fmod <- ~te(age, year, k = c(3, 25))
# srmod <- ~ bevholt(CV = 0.3)
srmod <- ~factor(replace(year, year<2001, 2001)) # I put 2001 insted of 1999 because one of the index starts in 2001 and the catch data for age 1 in 2000
qmod <- list(~factor(age),~factor(age), ~factor(age))
n1mod <- ~s(age, k = 3)
vmod <- list(~s(age, k = 3), ~1, ~1,~1)

fit01 <- sca(stk7, idx7, fmodel = fmod, qmodel = qmod, srmodel = srmod, vmodel = vmod, n1model = n1mod)
stk01 <- stk7 + fit01
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
fits01 <- scas(stks, idxs, fmodel=list(fmod), qmodel=list(qmod), srmodel=list(srmod), vmodel=list(vmod), n1model=list(n1mod), workers=n)
# update stock object with fit
stks <- stks + fits01
# add candidate fit
stks[[5]] <- stk7 + simulate(fit01, 250)
plot(window(stks, start=2000)) + theme(legend.position = "none") + scale_colour_manual(values = rep("black", n+1))

plot(stk7 + simulate(fit01, 250))

# Retro analysis plot with monrho values
source('./docs/retro_analysis_f.R') # retrospective function

results <- run_retro_analysis(stock = stk7, index = idx7, fit1 = fit01, fmod, qmod, srmod, vmod, n1mod)
results$rho_table <- results$rho_table %>% mutate(x = 2025, y = 0)
results$rho_table$qname <- c("F" = "F", "SSB" = "SSB", "Recruitment" = "Rec", "Catch" = "Catch")

new_names <- c("Rec" = "Recruitment", "SSB" = "SSB", "Catch" = "Catch", "F" = "F")
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
