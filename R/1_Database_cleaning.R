# ENVIRONMENTAL DATABASE (FROM AST SCRIPTS)
## NOTE: The relevant output here is: ./output/localities_total_unique_ch10m_Bio(v14).csv

# Libraries

library(readxl)
library(dplyr)
library(ggplot2)
library(raster)
library(tidyr)
library(sp)

# Cleaning NAs and repeated occurrences in pixels

localities <- read.csv("./data/occs/Localidades_Total_v14.csv")

dim(localities)
species <- unique(localities$NAME1)
occs <- localities %>%
  dplyr::select(XCOOR, YCOOR) %>%
  SpatialPoints()

# rasters

files <- list.files("./data/env/chelsa10m/bio", full.names = T)

predictors <- stack(files)
names(predictors)

# extract

env.tot <- raster::extract(predictors, occs, cellnumber = T, df = T)
env.tot.pt <- occs %>% data.frame() %>%
  cbind(localities$NAME1) %>%
  cbind(env.tot) %>%
  rename(NAME1 = `localities$NAME1`)

names(env.tot.pt)
dim(env.tot.pt)
ncol(env.tot.pt)

sum(is.na(env.tot.pt))
localities.clean <- env.tot.pt %>%
  #drop_na(6:11) %>%
  drop_na(6:24) %>%
  distinct(NAME1, cells, .keep_all = T)

head(localities.clean)
dim(localities.clean)
# Chelsa10m: after cleaning NAs and suppressing pixel duplicates, from 28763 records we conserved 18443 (Localidades_v14)

# save again:
write.csv(localities.clean, file = "./output/localities_total_unique_ch10m_Bio(v14).csv")

