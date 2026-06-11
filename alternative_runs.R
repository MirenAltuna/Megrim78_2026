
load("Input/bootstrap/data/stock/meg78_stock_sop_BE_corrected.RData")

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

# plusgroup at 7
stk7 <- setPlusGroup(stock, 7)

idx7 <- index_sep
idx7[[1]] <- reduce_pg_index(idx7[[1]], 7)
idx7[[2]] <- reduce_pg_index(idx7[[2]], 7)
idx7[[2]]@index[,"2017"] <- NA
idx7[[3]] <- reduce_pg_index(idx7[[3]], 7)
idx7[[3]]@index[,"2017"] <- NA
idx7[[2]]@index[ac(1:3),ac(2015:2021)] <- NA # IRLFR
idx7[[1]]@index[ac(1:3),ac(2015:2021)] <- NA # Porcupine
idx7[["EVHOE"]]@index[,ac(2019:2020)] <- NA

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
