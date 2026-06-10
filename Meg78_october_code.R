
library(a4adiags)
library(gridExtra)
library(grid)
library(lattice)
library(RColorBrewer)
library(dplyr)
library(FLCore)
library(icesAdvice)
library(reshape2)
library(ggplot2)
library(readxl)
library(tidyverse)
library(tibble)
library(icesDatras)
library(proxy)
library(stringr)
library(FLa4a)

set.seed(1234)

load("model/MegFit_Final.Rdata")
load("model/SubModels_Final.Rdata")

## functions ----
source('./docs/retro_analysis_f_mcmc.R') # retrospective function

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

mcmc_ctrl <- SCAMCMC(
  mcmc   = 5000,  # total iterations
  mcsave = 1, mcdiag = TRUE)     # save every 1th

## Base case model configuration

qmod <- list(~s(age, k = 4), ~s(age, k = 4))
fmod <- ~te(age, year, k = c(5, 25))
srmod <- ~ bevholt(CV=0.3)
n1mod <- ~s(age, k = 3)
vmod <- list(~s(age, k = 3), ~1, ~1)

# alternatives to improve the indices residuals for age 7+
# 
# qmod <- list(
#   ~s(pmin(age, 5), k = 4),
#   ~s(pmin(age, 4), k = 3)
# )
# 
# fmod <- ~te(pmin(age, 6), year, k = c(3, 28))

## The fit ----

fits <- sca(stock, index, fmodel = fmod, qmodel = qmod,
            srmodel = srmod, vmod = vmod, n1mod = n1mod,
            fit = "MCMC", mcmc = mcmc_ctrl)

# fitSumm(fits)

stks <- stock + fits

fit1 <- fitMedian(fits) # Objective: acurate between 0.25 and 0.40
stk1 <- qapply(stks, iterMedians)

### Retrospective analysis plot

results <- run_retro_analysis_mcmc(stock, index,
                                   fits, fmod, qmod,
                                   srmod, vmod, n1mod,
                                   mcmc_ctrl,
                                   back = 5)

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
a  <- xyplot(data~age,groups=year,stk1@harvest,type='b',ylim=c(0,1),ylab='F',main='Fishing mortality')
a1 <- xyplot(data~age,groups=year,data=fit1$qmodel[1],type='b',ylab='Catchability',main="Porcupine")
a2 <- xyplot(data~age,groups=year,data=fit1$qmodel[2],type='b',ylab='Catchability',main="IRLFR")
grid.arrange(a,a1,a2,ncol=2)

### Residuals plot

res <- residuals(fits, stock, index)

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

res_median <- residualsMedian(res)
plot(res_median)

### Observed and predicted catches plot

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
