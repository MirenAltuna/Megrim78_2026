
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
library(matrixStats)

set.seed(1234)

# Working directory
setwd("C:/USE/GitHub/Megrim78_2026")

## load the data
load("Input/bootstrap/data/stock/meg78_stock_sop_BE_corrected.RData")

## functions ----
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


## Put index weight as stock.wt  ----

df <- read.xlsx("C:/Users/maltuna/OneDrive - AZTI/Oilarra/WGBIE/WGBIE_2026/0.Original_country_data/Accession&email/Surveys/Weight/Hans/MegIndexMeanLen.xlsx")


ggplot(df, aes(x = Year, y = MeanWeight, color = factor(Age), group = Age)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  labs(
    x = "Year",
    y = "Mean weight",
    color = "Age",
    title = "Mean weight at age by year"
  ) +
  theme_bw()


# Pasar a kg
df2 <- df %>%
  mutate(MeanWeight = MeanWeight / 1000)

# Media por edad de los tres primeros años
mean_0305 <- df2 %>%
  filter(Year %in% c(2003, 2004, 2005)) %>%
  group_by(Age) %>%
  summarise(MeanWeight = mean(MeanWeight, na.rm = TRUE),
            .groups = "drop")

# Expandir 1984-2002 usando dicha media # 1984tik 2025era. Indizea 2003tik aurrera.
hist_wt <- expand.grid(
  Age = unique(mean_0305$Age),
  Year = 1984:2002
) %>%
  left_join(mean_0305, by = "Age")

# Juntar históricos + observados
wt_all <- bind_rows(
  hist_wt,
  df2 %>% select(Age, Year, MeanWeight)
) %>%
  arrange(Age, Year)

# Matriz edad x año
wt.mat <- wt_all %>%
  pivot_wider(names_from = Year, values_from = MeanWeight) %>%
  arrange(Age)

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

## plus group at 7  ----
stock <- setPlusGroup(stock, 7)

###### ---------------------------------------------------------------------------- AGE ----

## comparative plots  ----

sw <- as.data.frame(stock@stock.wt)
cw <- as.data.frame(stock@catch.wt)

comp <- left_join(
  sw[, c("age","year","data")],
  cw[, c("age","year","data")],
  by = c("age","year"),
  suffix = c(".stock", ".catch")
)

# Líneas por edad
ggplot(comp, aes(year)) +
  geom_line(aes(y = data.stock, colour = "stock.wt")) +
  geom_line(aes(y = data.catch, colour = "catch.wt")) +
  facet_wrap(~ age, scales = "free_y") +
  labs(y = "Weight (kg)", colour = "") +
  theme_bw()


## Diferencia absoluta
comp$diff <- comp$data.catch - comp$data.stock

ggplot(comp,
       aes(x = year, y = age, fill = diff)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0
  ) +
  labs(fill = "Catch - Stock") +
  theme_bw()

## Ratio catch/stock
comp$ratio <- comp$data.catch / comp$data.stock

ggplot(comp,
       aes(x = year, y = age, fill = ratio)) +
  geom_tile() +
  scale_fill_viridis_c() +
  labs(fill = "Catch/Stock") +
  theme_bw()

## Comparación 1:1
ggplot(comp,
       aes(x = data.stock,
           y = data.catch,
           colour = factor(age))) +
  geom_point(alpha = 0.7) +
  geom_abline(slope = 1, intercept = 0,
              linetype = 2) +
  labs(
    x = "stock.wt",
    y = "catch.wt",
    colour = "Age"
  ) +
  theme_bw()

## Diferencia media por edad
comp %>%
  group_by(age) %>%
  summarise(
    stock = mean(data.stock, na.rm = TRUE),
    catch = mean(data.catch, na.rm = TRUE)
  ) %>%
  tidyr::pivot_longer(c(stock, catch)) %>%
  ggplot(aes(age, value, colour = name)) +
  geom_line(linewidth = 1.2) +
  geom_point() +
  theme_bw() +
  labs(y = "Mean weight")


###### ---------------------------------------------------------------------------- LENGHT ----

## INDEX ----

df <- read.csv("C:/Users/maltuna/OneDrive - AZTI/Oilarra/WGBIE/WGBIE_2026/0.Original_country_data/Accession&email/Surveys/Weight/Hans/MegBio.csv")

index_nlen <- df %>%
  dplyr::count(Year, LngtClassCm)

ggplot(index_nlen, aes(x = LngtClassCm, y = n, colour = factor(Year))) +
  geom_line(linewidth = 1) +
  geom_point(size = 1.5) +
  labs(
    x = "Length class (cm)",
    y = "Number of individuals",
    colour = "Year"
  ) +
  theme_bw()

index_nlen$Source <- "Index"

index_prop <- index_nlen %>%
  group_by(Year) %>%
  mutate(prop = n / sum(n),
         Year = as.character(Year)) %>%
  ungroup() %>%
  select(Year, LngtClassCm, prop) %>%
  mutate(Source = "Index") %>% filter(Year > 2021)

## CATCH ----

ld_country <- read.xlsx("Input/Ctry_info/TLD_cntry.xlsx")

# 2025 ----
fr <- read.xlsx("C:/Users/maltuna/OneDrive - AZTI/Oilarra/WGBIE/WGBIE_2026/0.Original_country_data/Accession&email/Catch_length/FRA_2025_WGBIE_meg.27.7b-k8abd.xlsx", sheet = "SD")
sp <- read.xlsx("C:/Users/maltuna/OneDrive - AZTI/Oilarra/WGBIE/WGBIE_2026/0.Original_country_data/Accession&email/Catch_length/SP_2025_WGBIE_meg.27.7b-k8abd.xlsx", sheet = "SD")
uk <- read.xlsx("C:/Users/maltuna/OneDrive - AZTI/Oilarra/WGBIE/WGBIE_2026/0.Original_country_data/Accession&email/Catch_length/UK_2025_WGBIE_meg.27.7b-k8abd.xlsx", sheet = "SD")

colnames(uk) <- colnames(sp) <- colnames(fr)
unique(sp$UnitAgeOrLength)
unique(fr$UnitAgeOrLength)
unique(uk$UnitAgeOrLength)

fr$AgeLength <- as.numeric(fr$AgeLength) /10
fr$UnitAgeOrLength <- "cm"

fr_sp25 <- rbind(fr, sp)
fr_sp_uk25 <- rbind(fr_sp25, uk)

fr_sp_uk25$Year <- "2025" 

# 2024  ----
fr <- read.xlsx("C:/Users/maltuna/OneDrive - AZTI/Oilarra/WGBIE/WGBIE_2025/0.Original_country_data/Accession&email/Catch_length/FRA_2024_WGBIE_meg.27.7b-k8abd.xlsx", sheet = "SD")
sp <- read.xlsx("C:/Users/maltuna/OneDrive - AZTI/Oilarra/WGBIE/WGBIE_2025/0.Original_country_data/Accession&email/Catch_length/SP_2024_WGBIE_meg.27.7b-k8abd.xlsx", sheet = "SD")
uk <- read.xlsx("C:/Users/maltuna/OneDrive - AZTI/Oilarra/WGBIE/WGBIE_2025/0.Original_country_data/Accession&email/Catch_length/UK_2024_WGBIE_meg.27.7b-k8abd.xlsx", sheet = "SD")

colnames(uk) <- colnames(sp) <- colnames(fr)
unique(sp$UnitAgeOrLength)
unique(fr$UnitAgeOrLength)
unique(uk$UnitAgeOrLength)

fr$AgeLength <- as.numeric(fr$AgeLength) /10
fr$UnitAgeOrLength <- "cm"

fr_sp24 <- rbind(fr, sp)
fr_sp_uk24 <- rbind(fr_sp24, uk)

fr_sp_uk24$Year <- "2024" 

# 2023  ----
fr <- read.xlsx("C:/Users/maltuna/OneDrive - AZTI/Oilarra/WGBIE/WGBIE_2024_lcitores/0.Original_country_data/Accession&email/Catch_length/FRA_2023_WGBIE_meg.27.7b-k8abd.xlsx", sheet = "SD")
sp <- read.xlsx("C:/Users/maltuna/OneDrive - AZTI/Oilarra/WGBIE/WGBIE_2024_lcitores/0.Original_country_data/Accession&email/Catch_length/SP_2023_WGBIE_meg.27.7b-k8abd.xlsx", sheet = "SD")
uk <- read.xlsx("C:/Users/maltuna/OneDrive - AZTI/Oilarra/WGBIE/WGBIE_2024_lcitores/0.Original_country_data/Accession&email/Catch_length/UK_2023_WGBIE_meg.27.7b-k8abd.xlsx", sheet = "SD")

colnames(uk) <- colnames(sp) <- colnames(fr)
unique(sp$UnitAgeOrLength)
unique(fr$UnitAgeOrLength)
unique(uk$UnitAgeOrLength)

fr$AgeLength <- as.numeric(fr$AgeLength) /10
fr$UnitAgeOrLength <- "cm"

fr_sp23 <- rbind(fr, sp)
fr_sp_uk23 <- rbind(fr_sp23, uk)

fr_sp_uk23$Year <- "2023" 

# 2022  ----
fr <- read.xlsx("C:/Users/maltuna/OneDrive - AZTI/Oilarra/WGBIE/Previous_years_Airiondo/WGBIE_2023/0.Original_country_data/Accession&email/Catch_length/FRA_2022_WGBIE_meg.27.7b-k8abd.xlsx", sheet = "SD")
sp <- read.xlsx("C:/Users/maltuna/OneDrive - AZTI/Oilarra/WGBIE/Previous_years_Airiondo/WGBIE_2023/0.Original_country_data/Accession&email/Catch_length/SP_2022_WGBIE_meg.27.7b-k8abd.xlsx", sheet = "SD")
uk <- read.xlsx("C:/Users/maltuna/OneDrive - AZTI/Oilarra/WGBIE/Previous_years_Airiondo/WGBIE_2023/0.Original_country_data/Accession&email/Catch_length/UK_2022_WGBIE_meg.27.7b-k8abd.xlsx", sheet = "SD")

colnames(uk) <- colnames(sp) <- colnames(fr)
unique(sp$UnitAgeOrLength)
unique(fr$UnitAgeOrLength)
unique(uk$UnitAgeOrLength)

fr$AgeLength <- as.numeric(fr$AgeLength) /10
fr$UnitAgeOrLength <- "cm"

fr_sp22 <- rbind(fr, sp)
fr_sp_uk22 <- rbind(fr_sp22, uk)

fr_sp_uk22$Year <- "2022" 

# 2021  ----
fr <- read.xlsx("C:/Users/maltuna/OneDrive - AZTI/Oilarra/WGBIE/Previous_years_Airiondo/WGBIE_2022/0.Original_country_data/Accessions&email/Catch_length/FRA_2021_WGBIE_meg.27.7b-k8abd.xlsx", sheet = "SD")
sp <- read.xlsx("C:/Users/maltuna/OneDrive - AZTI/Oilarra/WGBIE/Previous_years_Airiondo/WGBIE_2022/0.Original_country_data/Accessions&email/Catch_length/SP_2021_WGBIE_meg.27.7b-k8abd.xlsx", sheet = "SD")
uk <- read.xlsx("C:/Users/maltuna/OneDrive - AZTI/Oilarra/WGBIE/Previous_years_Airiondo/WGBIE_2022/0.Original_country_data/Accessions&email/Catch_length/UK_2021_WGBIE_meg.27.7b-k8abd.xlsx", sheet = "SD")

colnames(uk) <- colnames(sp) <- colnames(fr)
unique(sp$UnitAgeOrLength)
unique(fr$UnitAgeOrLength)
unique(uk$UnitAgeOrLength)

fr$AgeLength <- as.numeric(fr$AgeLength) /10
fr$UnitAgeOrLength <- "cm"

fr_sp21 <- rbind(fr, sp)
fr_sp_uk21 <- rbind(fr_sp21, uk)

fr_sp_uk21$Year <- "2021" 

# Join all the years data  ----

fr_sp_uk2524 <- rbind(fr_sp_uk25, fr_sp_uk24)
fr_sp_uk2523 <- rbind(fr_sp_uk2524, fr_sp_uk23)
fr_sp_uk2522 <- rbind(fr_sp_uk2523, fr_sp_uk22)
fr_sp_uk2521 <- rbind(fr_sp_uk2522, fr_sp_uk21)

# Put in a correct format  ----
# Columnas donde -9 actúa como valor ausente
cols_con_minus9 <- c("AgeLength", "NumberCaught", "MeanWeight", "MeanLength",
                     "SampledCatch", "NumSamplesLngt", "NumLngtMeas",
                     "NumSamplesAge", "NumAgeMeas", "PlusGroup")

fr_sp_uk2521 <- fr_sp_uk2521 %>%
  mutate(across(all_of(cols_con_minus9), ~ ifelse(. == -9, NA, .)))

# Suma NumberCaught 
catch_nlen <- fr_sp_uk2521 %>%
  filter(!is.na(NumberCaught),         # elimina filas sin dato
         !is.na(AgeLength)) %>%         # elimina filas sin longitud
  mutate(LngtClassCm = AgeLength) %>% 
  group_by(Year, LngtClassCm, CatchCategory, Country) %>%
  summarise(
    n = sum(NumberCaught, na.rm = TRUE),  # suma total por grupo
    .groups = "drop"
  ) %>%
  arrange(LngtClassCm)

# proporciones
len_prop <- catch_nlen %>% 
  group_by(Year, Country, CatchCategory) %>%
  mutate(prop = n / sum(n))

# preparar la tabla de capturas
ld_country2 <- ld_country %>% 
  mutate(
    Year = as.character(year),
    Country = recode(country,
                     "SP" = "ES"),
    CatchCategory = recode(catch.cat,
                           "landings" = "L",
                           "discards" = "D")
  )

# unir capturas y proporciones
len_w <- len_prop %>% 
  left_join(
    ld_country2 %>%
      select(Year, Country, CatchCategory, total),
    by = c("Year","Country","CatchCategory")
  )

# ponderar
len_w <- len_w %>% 
  mutate(weighted = prop * total)

# composición final por año
catch_len_year <- len_w %>% 
  group_by(Year, LngtClassCm) %>%
  summarise(value = sum(weighted, na.rm = TRUE),
            .groups = "drop") %>%
  group_by(Year) %>%
  mutate(prop = value / sum(value))

catch_prop <- catch_len_year %>%
  ungroup() %>%
  select(Year, LngtClassCm, prop) %>%
  mutate(Source = "Catch") %>% filter(Year > 2021)

## CATCH and INDEX ----
len_compare <- bind_rows(catch_prop, index_prop) 

## Frecuencia relativa (%)

ggplot(len_compare,
       aes(x = LngtClassCm,
           y = prop,
           colour = Source)) +
  geom_line(linewidth = 1) +
  facet_wrap(~Year) +
  scale_y_continuous(labels = scales::percent) +
  labs(
    y = "Relative frequency",
    x = "Length class (cm)"
  ) +
  theme_bw()

## Comparar distribuciones acumuladas
df_cum <- len_compare %>%
  group_by(Year, Source) %>%
  arrange(LngtClassCm) %>%
  mutate(
    cumprop = cumsum(prop)
  )

ggplot(df_cum,
       aes(x = LngtClassCm,
           y = cumprop,
           colour = Source)) +
  geom_line(linewidth = 1) +
  facet_wrap(~Year) +
  scale_y_continuous(labels = scales::percent) +
  labs(
    y = "Cumulative frequency",
    x = "Length class (cm)"
  ) +
  theme_bw()
