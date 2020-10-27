# BROENNIMANN ECOSPAT CLIM-GRID Object for Overlap and Equivalence-Similarity tests
## This script produces the pca.grid object needed to calculate the overlap and do the equivalence and similarity tests
## Two pca.grid objects are produced: one only for the endemic species and the other for all the species plus AoE's environment.
## Ecospat requires that species have at least 5 points (or 5 valid rellocations) to calculate the kernel.
## Only 328 spp out of 386 spp of Bignonieae fulfil this condition. All the 166 endemics fulfill this condition.
## We have grids for all endemics: 166 spp.

# Libraries

library(tidyverse)
library(vegan)
library(ecospat)

### 0. Loading data ####
locs <-  read.csv("./output/localities_total_unique_ch10m_Bio(v14).csv", row.names = 1, stringsAsFactors = FALSE) # Localidades_v14
#locs <- read.csv("./output/2_Environment2.csv", stringsAsFactors = FALSE)[-1] # For Environment 2: Spp + AoE env space

# Loading PCA object
#load("./output/2_PCA_join_spp_AoE.rdata") 
load("./output/1_pca_amb_sel(v14).rdata")
pca_env.sel

# AoE and endemic species
ca_db <- read.csv("./data/biogeography/spp_listby_AoE.csv", stringsAsFactors = FALSE)[-1]
ca_db <- ca_db[str_detect(ca_db$Area_name, pattern = "D2.*"), ]

### 1. Extracting pca scores ####

# For endemic species 
sco.pca <- scores(pca_env.sel, scaling = -1)
scores.sp <- vector(mode = "list", length = length(unique(ca_db$Species)))
scores.clim <- sco.pca$sites

# ## For AoE (Using Environment2 - Load first: "./output/2_PCA_join_spp_AoE.rdata")
# scores.aoe <- vector(mode = "list", length(unique(ca_db$Area_name)))
# scores.clim <- sco.pca$sites

### 2. Extracting names ####

# For endemic species
locs_end <- locs %>%
  filter(NAME1 %in% unique(ca_db$Species))
spp_names <- unique(locs_end$NAME1)

# # For areas of endemism (Read first: locs <- "./output/2_Environment2.csv")
# locs_ppx_aoe <- locs %>%
#   filter(xname %in% unique(ca_db$Area_name))
# aoe_names <- unique(locs_ppx_aoe$xname)

### 3. Creating clim-grid object ####

# For endemics 
for (i in 1:166) scores.sp[[i]] <- sco.pca$sites[locs_end$NAME1 == spp_names[i],]

# For AoEs
#for (i in 1:28) scores.aoe[[i]] <- sco.pca$sites[locs_ppx_aoe$xname == aoe_names[i],]

# Density grids for all endemic species and AoE

# For endemics
pca.grid.spp <- vector(mode = "list", length = 166)

# # For AoE
# pca.grid.aoe <- vector(mode = "list", length = 28)

# For Endemics
a <- Sys.time()
k <- 0
nogrid_index <- 0
for (i in 1:166) { ## For endemics 1 in 1:166
  
  if (length(scores.sp[[i]])/2 >= 5) {
    
    pca.grid.spp[[i]] <- ecospat.grid.clim.dyn(scores.clim, scores.clim, scores.sp[[i]], R = 100)
    
  } else {
    
    pca.grid[[i]] <- "<5 relocations"
    k <- k + 1
    nogrid_index <- c(nogrid_index, i)
    
  }
  
}
b <- Sys.time()
a - b # This loop took -2.137236 mins with 166 spp using PCA for overlap (Spp only)
# This loop took -8.188191 mins with 166 spp using PCA for Environment2(Spp + AoEs env)
k
# k = 0/166. THE GRIDS FOR ALL THE ENDEMIC SPECIES ARE AVAILABLE.

# # For AoE
# a <- Sys.time()
# k <- 0
# nogrid_index <- 0
# for (i in 1:28) { ## For endemics 1 in 1:166
#   
#   if (length(scores.aoe[[i]])/2 >= 5) {
#     
#     pca.grid.aoe[[i]] <- ecospat.grid.clim.dyn(scores.clim, scores.clim, scores.aoe[[i]], R = 100)
#     
#   } else {
#     
#     pca.grid[[i]] <- "<5 relocations"
#     k <- k + 1
#     nogrid_index <- c(nogrid_index, i)
#     
#   }
#   
# }
# b <- Sys.time()
# a - b # This loop took  -1.2908 mins with 28 aoe in Linux-pc 
# k 
# # K = 0/28. THE GRIDS FOR ALL THE AOES ARE AVAILABLE

# # In case of species with no grid:
# nogrid_index <- nogrid_index[-1]
# # <Write here the indexes of spp with no grids>
# # Exclude these species from pca.grid object:
# pca.grid <- pca.grid[-nogrid_index]

# Saving clean pca.grid object
saveRDS(pca.grid.spp, file = "./output/Grid_clim/1_ecospat_climgrid_endemicspp.rds")
pca.grid.spp <- readRDS(file = "./output/Grid_clim/1_ecospat_climgrid_endemicspp.rds")


# For clim.grid objects for spp only and spp+Aoe using Environment2:
# saveRDS(pca.grid.spp, file = "./output/Grid_clim/2_ecospat_climgrid_endemicspp.rds")
# pca.grid.spp <- readRDS(file = "./output/Grid_clim/2_ecospat_climgrid_endemicspp.rds")

# saveRDS(pca.grid.aoe, file = "./output/Grid_clim/2_ecospat_climgrid_AoE.rds")
# pca.grid.aoe <- readRDS(file = "./output/Grid_clim/2_ecospat_climgrid_AoE.rds")

length(pca.grid.spp) # 166 endemic spp
length(pca.grid.aoe) # 28 AoEs



