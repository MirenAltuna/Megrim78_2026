
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
library(FLBRP)

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
iterMedians <- function(x) {
  # x es un FLQuant con iter
  apply(x, c(1,2,3,4,5), median, na.rm = TRUE) |>
    FLQuant(dimnames = c(dimnames(x)[1:5], list(iter = "1")),
            units    = units(x))
}
fitMedian <- function(fit) {
  
  # 1. Mediana del stock (FLStock)
  iterMedians <- function(x) {
    apply(x, c(1,2,3,4,5), median, na.rm = TRUE) |>
      FLQuant(dimnames = c(dimnames(x)[1:5], list(iter = "1")),
              units    = units(x))
  }
  
  stock_med <- qapply(fit, iterMedians)
  
  # 2. Predecir todas las iteraciones
  fitted_all <- predict(pars(fit))
  
  # 3. Mediana de cada FLQuant dentro de predict()
  fitted_med <- lapply(fitted_all, function(sublist) {
    lapply(sublist, function(q) {
      arr <- array(q[], dim = dim(q), dimnames = dimnames(q))
      med <- apply(arr, c(1,2,3,4,5), median, na.rm = TRUE)
      dim(med) <- c(dim(q)[1:5], 1)
      dn <- dimnames(q)
      dn$iter <- "1"
      FLQuant(med, dimnames = dn, units = units(q))
    })
  })
  
  # 4. Reconstruir un objeto SCA con stock + fitted_med
  out <- list(
    stock   = stock_med,
    stkmodel = fitted_med$stkmodel,
    qmodel   = fitted_med$qmodel,
    vmodel   = fitted_med$vmodel
  )
  
  return(out)
}
residualsMedian <- function(res) {
  
  # Función interna para colapsar un FLQuant a la mediana
  iterMedians <- function(x) {
    arr <- array(x[], dim = dim(x), dimnames = dimnames(x))
    med <- apply(arr, c(1,2,3,4,5), median, na.rm = TRUE)
    
    # reconstruir FLQuant con iter = 1
    dim(med) <- c(dim(x)[1:5], 1)
    dn <- dimnames(x)
    dn$iter <- "1"
    
    FLQuant(med, dimnames = dn, units = units(x))
  }
  
  # Aplicar a cada FLQuant dentro del objeto a4aFitResiduals
  res_med_list <- lapply(res@.Data, iterMedians)
  
  # Reconstruir un objeto a4aFitResiduals
  out <- new("a4aFitResiduals")
  out@.Data <- res_med_list
  out@names <- res@names
  out@desc  <- paste(res@desc, "(median collapsed)")
  out@lock  <- FALSE
  
  return(out)
}

source('./docs/retro_analysis_f.R') # retrospective function simple
source('./docs/retro_analysis_f_mcmc.R') # retrospective function mcmc

## Load the data ----
load("Input/bootstrap/data/stock/meg78_stock_sop_BE_corrected.RData")
load("Input/IGFS_EVHOE_index/index_sep.RData")

## Prepare the data ----
# Put index weight as stock.wt

df <- read.xlsx("Input/IGFS_EVHOE_index/MegIndexMeanLen.xlsx")

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

# Put plus group at age 7
stk7 <- setPlusGroup(stock, 7)
stk7@catch.n['1',as.character(1984:2000)] <- NA # We do not really believe that the increase in 1-year-olds in the catch is real so we shouldn't formulate a model that treats it as real.

idx7 <- index_sep
idx7[[1]] <- reduce_pg_index(idx7[[1]], 7)
idx7[[2]] <- reduce_pg_index(idx7[[2]], 7)
idx7[[3]] <- reduce_pg_index(idx7[[3]], 7)
idx7[[3]]@index[,"2017"] <- NA

#index(idx7[[3]]) <- replaceZeros(index(idx7)[[3]], frac=0.1)

# idx7[[1]]@index[ac(1:3),ac(2015:2021)] <- NA # Porcupine
# idx7[[2]]@index[ac(1:3),ac(2015:2021)] <- NA # IGFS
# idx7[[3]]@index[ac(1:3),ac(2015:2021)] <- NA # EVHOE

## The best setting until now ----

srmod <- ~ bevholt(CV = 0.1)
fmod <- ~s(replace(age, age>6,6), k = 6, by = breakpts(year, 2013)) + s(year, k=20)# +ti(age, year, k=c(6,10))
qmod <- list(~s(age, k = 5), ~s(age, k = 5), ~s(age, k = 5))
n1mod <- ~s(age, k = 5)
vmod <- list(~s(age, k = 3), ~1, ~1, ~1) 

### Simple RUN ----
#idx7[[3]] <- window(idx7[[3]], 2003)
#stk7 <- window(stk7, 2001)
# ctn <- catch.n(stk7)
# ctn[] <- 0.5
# catch.n(stk7) <- FLQuantDistr(catch.n(stk7), ctn)

fit01 <- sca(stk7, idx7, fmodel = fmod, qmodel = qmod, srmodel = srmod, vmodel = vmod, n1model = n1mod)
stk01 <- stk7 + fit01

#### Residuals ----
res01 <- residuals(fit01, stk7, idx7)
plot(res01)

plot(computeCatchDiagnostics(fit01, stk7))

#### Selectivity and catchability plot  ----
fitted <- predict(pars(fit01))

a  <- xyplot(data~age,groups=year,stk01@harvest,type='b',ylim=c(0,1),ylab='F',main='Fishing mortality')
a1 <- xyplot(data~age,groups=year,data=fitted$qmodel[1],type='b',ylab='Catchability',main="Porcupine")
a2 <- xyplot(data~age,groups=year,data=fitted$qmodel[2],type='b',ylab='Catchability',main="IGFS")
a3 <- xyplot(data~age,groups=year,data=fitted$qmodel[3],type='b',ylab='Catchability',main="EVHOE")
grid.arrange(a,a1,a2, a3, ncol=2)

#### Observed and predicted catches plot  ----

fits <- simulate(fit01, 1000)
stks <- stk7 + fits

pred <- catch(stks)   # FLQuant con iter = 1000

pred_median <- apply(pred, 2, median, na.rm = TRUE)
pred_p5     <- apply(pred, 2, quantile, 0.05, na.rm = TRUE)
pred_p95    <- apply(pred, 2, quantile, 0.95, na.rm = TRUE)

df <- tibble(
  Year      = as.numeric(dimnames(pred)$year),
  Observed  = as.numeric(catch(stock)),
  Median    = as.numeric(pred_median),
  P5        = as.numeric(pred_p5),
  P95       = as.numeric(pred_p95)
)

ggplot(df, aes(x = Year)) +
  geom_ribbon(aes(ymin = P5, ymax = P95),
              fill = "steelblue", alpha = 0.25) +
  geom_line(aes(y = Median, color = "Predicted"), size = 1.2) +
  geom_line(aes(y = Observed, color = "Observed"), size = 1.2) +
  scale_color_manual(values = c("Observed" = "black",
                                "Predicted" = "steelblue4")) +
  labs(x = "Year",
       y = "Catch (tonnes)",
       color = "Type") +
  theme_bw()

#### Retro analysis plot with monrho values  ----
results <- run_retro_analysis(stk7, idx7, fit01, fmod, qmod, srmod, vmod, n1mod)

fit01a <- sca(window(stk7, end=2024), window(idx7, end=2024), fmodel = fmod, qmodel = qmod, srmodel = srmod, vmodel = vmod, n1model = n1mod)


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

### MCMC RUN ----

## original
# mcmc_ctrl <- SCAMCMC(
#   mcmc   = 100000, 
#   mcsave = 200, mcseed = 10)

## Proposed by Ernesto
mcmc_ctrl <- SCAMCMC(
  mcmc   = 1000000, 
  mcsave = 2000, mcseed = 10)      

fits <- sca(stk7, idx7, fmodel = fmod, qmodel = qmod, srmodel = srmod, vmodel = vmod, n1model = n1mod,
            fit = "MCMC", mcmc = mcmc_ctrl)
stks <- stk7 + fits

fit1 <- fitMedian(fits)
stk1 <- qapply(stks, iterMedians)

#### Residuals ----
res <- residuals(fits, stk7, idx7)

res_median <- residualsMedian(res)
plot(res_median)

#### Selectivity and catchability plot ----
a  <- xyplot(data~age,groups=year,stk1@harvest,type='b',ylim=c(0,1),ylab='F',main='Fishing mortality')
a1 <- xyplot(data~age,groups=year,data=fit1$qmodel[1],type='b',ylab='Catchability',main="Porcupine")
a2 <- xyplot(data~age,groups=year,data=fit1$qmodel[2],type='b',ylab='Catchability',main="IGFS")
a3 <- xyplot(data~age,groups=year,data=fit1$qmodel[3],type='b',ylab='Catchability',main="EVHOE")
grid.arrange(a,a1,a2, a3, ncol=2)

#### Observed and predicted catches plot  ----
pred <- catch(stks)

pred_median <- apply(pred, 2, median, na.rm = TRUE)
pred_p5     <- apply(pred, 2, quantile, 0.05, na.rm = TRUE)
pred_p95    <- apply(pred, 2, quantile, 0.95, na.rm = TRUE)

df <- tibble(
  Year      = as.numeric(dimnames(pred)$year),
  Observed  = as.numeric(catch(stock)),
  Median    = as.numeric(pred_median),
  P5        = as.numeric(pred_p5),
  P95       = as.numeric(pred_p95)
)

ggplot(df, aes(x = Year)) +
  geom_ribbon(aes(ymin = P5, ymax = P95),
              fill = "steelblue", alpha = 0.25) +
  geom_line(aes(y = Median, color = "Predicted"), size = 1.2) +
  geom_line(aes(y = Observed, color = "Observed"), size = 1.2) +
  scale_color_manual(values = c("Observed" = "black",
                                "Predicted" = "steelblue4")) +
  labs(x = "Year",
       y = "Catch (tonnes)",
       color = "Type") +
  theme_bw()

#### Retro  ----
results <- run_retro_analysis_mcmc(stk7, idx7,
                                   fits, fmod, qmod,
                                   srmod, vmod, n1mod,
                                   mcmc_ctrl,
                                   back = 5)

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

#### MCMC diagnostic  ----
library(coda)

fits_b <- burnin(fits, 400)

##### Traceplot  ----
fitmc01.mc <- FLa4a::as.mcmc(fits_b)
traceplot(mcmc.list(mc01=fitmc01.mc[,1]), lwd=1.5, col=c(2,4), lty=1)

##### Autocorrelation and crosscorrelation analysis  ----
acfplot(fitmc01.mc[,1], lwd=3, ylim=c(-1, 1))
crosscorr.plot(fitmc01.mc)

##### Geweke diagnostic  ----
geweke.plot(fitmc01.mc[,1])

##### Cumulative means  ----
cm01 <- fitmc01.mc[,1]
cm01 <- cumsum(cm01) / seq_along(cm01)
plot(cm01, type="l", xlab="samples", ylab="mean")

##### Distribution density  ----
densplot(fitmc01.mc[,1])

##### Acceptance rate  ----
data.frame(hessian_scale=c(fitSumm(fits)), row.names=rownames(fitSumm(fits)))

data.frame(hessian_scale=c(autocorr.diag(fitmc01.mc[,1])),
           row.names=c(0, 1, 5, 10, 50))
