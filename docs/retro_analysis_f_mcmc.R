
run_retro_analysis_mcmc <- function(stock, index, fit1, fmod, qmod, srmod, vmod, n1mod, mcmc_ctrl, back = 5) {

  # iterMedians function
  iterMedians <- function(x) {
    # x es un FLQuant con iter
    apply(x, c(1,2,3,4,5), median, na.rm = TRUE) |>
      FLQuant(dimnames = c(dimnames(x)[1:5], list(iter = "1")),
              units    = units(x))
  }
  
  # Crear objeto retro
  retro <- split(1:back, 1:back)
  retro <- lapply(retro, function(x) {
    yr <- range(stock)["maxyear"] - x
    stk <- window(stock, end = yr)
    idx <- window(index, end = yr)
    stk + sca(stk, idx, fmodel = fmod, qmodel = qmod, srmodel = srmod, vmodel = vmod, n1model = n1mod,
              fit   = "MCMC",
              mcmc  = mcmc_ctrl)
  })
  
  retro_median <- lapply(retro, function(stock_obj) {
    qapply(stock_obj, iterMedians)
  })
  
  retro_median$"0" <- stock + fit1
  
  # Calcular Mohn's rho para F
  Retro_F <- data.frame(
    Y0 = c(fbar(retro_median$`0`)),
    Y1 = c(fbar(retro_median$`1`), NA),
    Y2 = c(fbar(retro_median$`2`), NA, NA),
    Y3 = c(fbar(retro_median$`3`), NA, NA, NA),
    Y4 = c(fbar(retro_median$`4`), NA, NA, NA, NA),
    Y5 = c(fbar(retro_median$`5`), NA, NA, NA, NA, NA)
  )
  rho_f <- mohn(Retro_F)
  
  # Calcular Mohn's rho para SSB
  Retro_SSB <- data.frame(
    Y0 = c(ssb(retro_median$`0`)),
    Y1 = c(ssb(retro_median$`1`), NA),
    Y2 = c(ssb(retro_median$`2`), NA, NA),
    Y3 = c(ssb(retro_median$`3`), NA, NA, NA),
    Y4 = c(ssb(retro_median$`4`), NA, NA, NA, NA),
    Y5 = c(ssb(retro_median$`5`), NA, NA, NA, NA, NA)
  )
  rho_ssb <- mohn(Retro_SSB)
  
  # Calcular Mohn's rho para Recruitment
  recr <- function(x) x@stock.n[1, ]
  Retro_R <- data.frame(
    Y0 = c(recr(retro_median$`0`)),
    Y1 = c(recr(retro_median$`1`), NA),
    Y2 = c(recr(retro_median$`2`), NA, NA),
    Y3 = c(recr(retro_median$`3`), NA, NA, NA),
    Y4 = c(recr(retro_median$`4`), NA, NA, NA, NA),
    Y5 = c(recr(retro_median$`5`), NA, NA, NA, NA, NA)
  )
  rho_rec <- mohn(Retro_R)
  
  # Crear tabla resumen
  rho_table <- data.frame(
    Metric = c("F", "SSB", "Recruitment", "Catch"),
    Mohns_Rho = c(rho_f, rho_ssb, rho_rec, NA)
  )
  
  rho_table <- rho_table %>%
    mutate(label = ifelse(is.na(Mohns_Rho), "", paste0("rho = ", round(as.numeric(Mohns_Rho), 3))))
  
  list(
    retro = retro_median,
    rho_table = rho_table
  )
  
}
