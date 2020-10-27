# Processing the ecospat results from SDrumont CCC: ENDEMIC SPECIES
## In this script all the results from ecospat overlap and test of equivalency and similarity are integrated 
## in one data frame. Besides, new information is added such as Geographical overlap between AoEs and 
## classes of niche overlap among others. The final output contains a table that integrates all the 
## ecospat analysis with information about AoEs.

library(dplyr)
library(ggplot2)
library(ggpubr)
library(grid)
library(gridExtra)
library(stringr)

### 0. Results of equivalency test for endemic species ####

# Results:

# NOTE: Load first one result only (ex. Ecospat EQUIVALENCE LOWER: ./end_niche_eq_lower.csv) 
# and run step by step the code below. Save the sugested objects along the script.
# Once the script is finished, erase all objects from the Global Environment, and start again
# with the second result.

eq_res <- read.csv("./data/ecospat_results/end_niche_eq_lower.csv", stringsAsFactors = FALSE) # 1. Ecospat EQUIVALENCE LOWER
eq_res <- read.csv("./data/ecospat_results/end_niche_sim_lower.csv", stringsAsFactors = FALSE) # 2. Ecospat SIMILARITY LOWER
eq_res <- read.csv("./data/ecospat_results/end_niche_eq_greater.csv", stringsAsFactors = FALSE) # 3. Ecospat EQUIVALENCE GREATER
eq_res <- read.csv("./data/ecospat_results/end_niche_sim_greater.csv", stringsAsFactors = FALSE) # 4. Ecospat SIMILARITY GREATER

# AoE species 
ca_db <- read.csv("./data/biogeography/spp_listby_AoE.csv", stringsAsFactors = FALSE)[-1]
ca_db <- ca_db[str_detect(ca_db$Area_name, pattern = "D2.*"), ]

# Species names from env data
locs <-  read.csv("./output/localities_total_unique_ch10m_Bio(v14).csv", row.names = 1, stringsAsFactors = FALSE) # Localidades_v14
locs_end <- locs %>% filter(NAME1 %in% unique(ca_db$Species))

# Identity code for SDumont results: Species 1 to 166
spp_code <- tibble(spp_names = unique(locs_end$NAME1), code = 1:166)
spp_code2 <- tibble(spp_names = unique(locs_end$NAME1), ID = 1:166)
#write.csv(spp_code2, "./output/AoE/Table_endemicspp_code_SDumont.csv")

# Inclusion of code in ca_db
code_vec <- vector(mode = "numeric", length = length(ca_db$Species))
for (i in seq_along(ca_db$Species)) {
  for (j in seq_along(spp_code$spp_names)){
    
    if (ca_db$Species[i] == spp_code$spp_names[j]) {
      
      code_vec[i] <- spp_code$code[j]
      
    }
    
  }
  
}
ca_db$Sp.1 <- code_vec
ca_db$Sp.2 <- code_vec

### 1. Classifying niche ovelarp degrees in classes following Rödder and Engler 2011 ####

eq_res_jp <- eq_res %>%
  mutate(D_classes = case_when(
    D <= 0.2 ~ "Limited",
    D > 0.2 & D <= 0.4 ~ "Low",
    D > 0.4 & D <= 0.6 ~ "Moderate",
    D > 0.6 & D <= 0.8 ~ "High",
    D > 0.8 ~ "Very_high"
  ),
  I_classes = case_when(
    I <= 0.2 ~ "Limited",
    I > 0.2 & I <= 0.4 ~ "Low",
    I > 0.4 & I <= 0.6 ~ "Moderate",
    I > 0.6 & I <= 0.8 ~ "High",
    I > 0.8 ~ "Very_high"
  ),
  Dp_class = ifelse(D_p <= 0.049999, "P<0.05", "P>=0.05"),
  Ip_class = ifelse(I_p <= 0.049999, "P<0.05", "P>=0.05")
  )

head(eq_res_jp)
head(ca_db)

#################### Geographical overlap endemics ##################################
## NOTE: The object can be load at the end of this section. This code is slow to run. 
### 2. Calculating geographical intersection area between endemic species ####

# library(sf)
# library(units)
# library(purrr)
# 
# eq_res_jp
# head(eq_res_jp)
# 
# SPP_1_vec <- eq_res_jp$Sp.1
# SPP_2_vec <- eq_res_jp$Sp.2
# 
# 
# # Checking for species with 2 or less unique occurrence points for which convex hulls cannot be calculated
# 
# locs_end %>%
#   group_by(NAME1) %>%
#   summarise()
# 
# less2vec <- locs_end %>%
#   group_by(NAME1) %>%
#   summarise(total = n()) %>%
#   filter(total <= 2) %>%
#   pull(NAME1)
# 
# # 21 species have 2 or less occurrence points.
# # To resolve the problem, I decided to put a buffer over them.
# # When applying buffer of 1° radius to estimate the area, 5/21 species have occurrence points very far a away between each other: 5, 8, 14, 16, 21
# 
# ### Function spp_geointersection
# ## This function calculates the intersection area in km^2 between two species from the endemic species loality-env dataframe.
# ## It is necessary to calculate previously all the possible combinations between two species taken from the total number of species.
# ## To ease the process species must have a numeric code that represents its position inside a vector of unique values obtained from the dataframe.
# ## Each element of the species pair must be separated in a different vector.
# ## A buffer of 1 degree is applied to species with 2 or less occurrences points. 
# ## Arguments:
# ## 1) vec1 and vec2: Two numerical vectors representing each element of the pair of numbers that identified the two species.
# ## 2) df_locs: endemic species env dataframe with column NAME1 containing species names.
# 
# 
# spp_geointersection <- function(vec1, vec2, df_locs) {
#   
#   if (length(vec1) != length(vec2)) {
#     
#     stop("vec1 and vec2 must have the same lenght")
#     
#   }
#   
#   library(sf)
#   library(dplyr)
#   library(units)
#   
#   geo_inter_vec <- vector(mode = "numeric")
#   
#   spp_nvec <- unique(df_locs$NAME1)
#   
#   less2vec <- df_locs %>%
#     group_by(NAME1) %>%
#     summarise(total = n()) %>%
#     filter(total <= 2) %>%
#     pull(NAME1)
#   
#   if (spp_nvec[vec1] %in% less2vec) {
#     
#     ch1 <- df_locs %>%
#       filter(NAME1 == spp_nvec[vec1]) %>%
#       st_as_sf(coords = c("XCOOR", "YCOOR"), crs = 4326) %>%
#       summarise(geometry = st_combine(geometry)) %>%
#       st_buffer(dist = 1)
#     
#   } else {
#     
#     ch1 <- df_locs %>%
#       filter(NAME1 == spp_nvec[vec1]) %>%
#       st_as_sf(coords = c("XCOOR", "YCOOR"), crs = 4326) %>%
#       summarise(geometry = st_combine(geometry)) %>%
#       st_convex_hull()
#     
#   }
#   
#   if (spp_nvec[vec2] %in% less2vec) {
#     
#     ch2 <- df_locs %>%
#       filter(NAME1 == spp_nvec[vec2]) %>%
#       st_as_sf(coords = c("XCOOR", "YCOOR"), crs = 4326) %>%
#       summarise(geometry = st_combine(geometry)) %>%
#       st_buffer(dist = 1)
#     
#   } else {
#     
#     ch2 <- df_locs %>%
#       filter(NAME1 == spp_nvec[vec2]) %>%
#       st_as_sf(coords = c("XCOOR", "YCOOR"), crs = 4326) %>%
#       summarise(geometry = st_combine(geometry)) %>%
#       st_convex_hull()
#     
#   }
#   
#   intch12 <- st_intersection(ch1, ch2)
#   
#   geo_inter_vec <- set_units(st_area(intch12), "km^2") %>% as.vector()
#   
#   if (identical(geo_inter_vec, numeric(0)) == TRUE) {
#     
#     geo_inter_vec <- 0
#     
#   }
#   
#   return(geo_inter_vec)
#   
# } 
# spp_geointersection2 <- function(vec1, vec2, df_locs) {
#   
#   if (length(vec1) != length(vec2)) {
#     
#     stop("vec1 and vec2 must have the same lenght")
#     
#   }
#   
#   library(sf)
#   library(dplyr)
#   library(units)
#   
#   geo_inter_vec <- vector(mode = "numeric")
#   
#   spp_nvec <- unique(df_locs$NAME1)
#   
#   ch1 <- df_locs %>%
#     filter(NAME1 == spp_nvec[vec1]) %>%
#     st_as_sf(coords = c("XCOOR", "YCOOR"), crs = 4326) %>%
#     summarise(geometry = st_combine(geometry)) %>%
#     st_convex_hull()
#   
#   ch2 <- df_locs %>%
#     filter(NAME1 == spp_nvec[vec2]) %>%
#     st_as_sf(coords = c("XCOOR", "YCOOR"), crs = 4326) %>%
#     summarise(geometry = st_combine(geometry)) %>%
#     st_convex_hull()
#   
#   intch12 <- st_intersection(ch1, ch2)
#   
#   geo_inter_vec <- set_units(st_area(intch12), "km^2") %>% as.vector()
#   
#   if (identical(geo_inter_vec, numeric(0)) == TRUE) {
#     
#     geo_inter_vec <- 0
#     
#   }
#   
#   return(geo_inter_vec)
#   
# } ## WITHOUT BUFFER
# 
# # Test with one pair of species
# a <- spp_geointersection(vec1 = 1, vec2 = 1, df_locs = locs_end)
# 
# # Calculating interceptions for all species
# # a <- Sys.time()
# # interceptions_geo <- purrr::pmap(list(vec1 = SPP_1_vec, vec2 = SPP_2_vec),
# #                                  .f = spp_geointersection,
# #                                  df_locs = locs_end)
# # b <- Sys.time()
# # a - b # ime difference of -8.394915 mins
# 
# index <- which(unlist(interceptions_geo) != unlist(interceptions_geo2))
# 
# View(cbind(unlist(interceptions_geo2)[index], unlist(interceptions_geo)[index]))
# # 1083/27394 comparisons have a geographical overlap equal to 0 before applying the buffer. 
# 
# # Test with parallel processing (Not significant gain. The process seems to be not adequate for parallel processing)
# # library(furrr)
# # a <- Sys.time()
# # numCores <- detectCores()
# # cl <- makeCluster(4)
# # interceptions_geo <- future_pmap(list(vec1 = SPP_1_vec, vec2 = SPP_2_vec),
# #             .f = spp_geointersection,
# #             df_locs = locs_end)
# # 
# # stopCluster(cl)
# # b <- Sys.time()
# # a - b

#interceptions_geo_unlist <-  unlist(interceptions_geo)
#save(interceptions_geo_unlist, file = "./output/AoE/AoE_end_spp_geooverlap.rdata")

### 2.1. Classifying geographical area overlap size ####

# Assigning intersections to eq_res_jp dataframe

# Geo overlap object
load("./output/Geo_overlap/AoE_end_spp_geooverlap.rdata")

glimpse(interceptions_geo_unlist)
length(interceptions_geo_unlist)  # 27390

eq_res_jp$Geo_overlap <- interceptions_geo_unlist 
glimpse(eq_res_jp)

### K-means for overlap area classification

set.seed(1)
kme_df <- eq_res_jp %>%
  filter(Geo_overlap != 0) %>%
  dplyr::select(Geo_overlap) %>%
  kmeans(centers = 4)

# Looking at centers
kme_df$centers

eq_res_jp %>%
  filter(Geo_overlap != 0) %>%
  ggplot(aes(x = Geo_overlap)) +
  geom_histogram() +
  geom_vline(xintercept = c(36085.89, 269101.92, 905311.11, 4939384.04), color = c("blue", "red", "green", "black"))

# Intermediate dataframe with clusters and geo_overlap classes
interm_df_kme <- eq_res_jp %>%
  filter(Geo_overlap != 0) %>%
  dplyr::select(X, Geo_overlap)

interm_df_kme <- interm_df_kme %>%
  mutate(Geo_cluster = factor(kme_df$cluster, levels = c("1", "2", "3", "4"), ordered = TRUE),
         Geo_classes = case_when(
           Geo_cluster == 1 ~ "Low",
           Geo_cluster == 2 ~ "Medium",
           Geo_cluster == 3 ~ "High",
           Geo_cluster == 4 ~ "Very_high"
         ))

glimpse(interm_df_kme)

# Feeding main dataframe with Geo_classes
eq_res_jp <- eq_res_jp %>%
  left_join(interm_df_kme, by = c("X", "Geo_overlap")) %>%
  mutate(Geo_classes = ifelse(is.na(Geo_cluster), "Disjunct", Geo_classes)) %>%
  dplyr::select(-Geo_cluster)

eq_res_jp$Geo_classes <- factor(eq_res_jp$Geo_classes, levels = c("Disjunct", "Low", "Medium", "High", "Very_high"), 
                                ordered = TRUE)

glimpse(eq_res_jp)

### 3. Creating columns to identify the analsysis ####

# NOTE: This code identifies the analysis (Similatity = "Sim", Equivalence = "Eq") 
# and the alternative hypotheses in ecospat (Lower = "L", Greater = "G").
# Execute the line corresponding to the result that is currently running.

# 1. Ecospat EQUIVALENCE LOWER
# eq_res_jp$Analysis <- "Eq"
# eq_res_jp$Alternative <- "L"
# 
# # 2. Ecospat SIMILARITY LOWER
# eq_res_jp$Analysis <- "Sim"
# eq_res_jp$Alternative <- "L"
# 
# # 3. Ecospat EQUIVALENCE GREATER
# eq_res_jp$Analysis <- "Eq"
# eq_res_jp$Alternative <- "G"
# 
# # 4. Ecospat SIMILARITY GREATER
eq_res_jp$Analysis <- "Sim"
eq_res_jp$Alternative <- "G"

glimpse(eq_res_jp)

### Saving eq_res_jp
#write.csv(eq_res_jp, "./output/AoE/Table_eqivalence_lower_JPNG.csv") # 1. Ecospat EQUIVALENCE LOWER
#write.csv(eq_res_jp, "./output/AoE/Table_similarity_lower_JPNG.csv") # 2. Ecospat SIMILARITY LOWER
#write.csv(eq_res_jp, "./output/AoE/Table_eqivalence_greater_JPNG.csv") # 3. Ecospat EQUIVALENCE GREATER
#write.csv(eq_res_jp, "./output/AoE/Table_similarity_greater_JPNG.csv") # 4. Ecospat SIMILARITY GREATER


### 4. Dataframe with pair-wise comparisons among all endemic species ####

# NOTE: Run this sections once all the dataframes in the section three are saved.

alleqlow <- read.csv("./output/AoE/Table_eqivalence_lower_JPNG.csv", stringsAsFactors = FALSE)[, -c(1:2)] # EQUIVALENCE LOWER
allsimlow <- read.csv("./output/AoE/Table_similarity_lower_JPNG.csv", stringsAsFactors = FALSE)[, -c(1:2)] # SIMILARITY LOWER
alleqgre <- read.csv("./output/AoE/Table_eqivalence_greater_JPNG.csv", stringsAsFactors = FALSE)[, -c(1:2)] # EQUIVALENCE GREATER
allsimgre <- read.csv("./output/AoE/Table_similarity_greater_JPNG.csv", stringsAsFactors = FALSE)[, -c(1:2)] # SIMILARITY GREATER

### 4.1. Unifying dataframes from all analyses (All pair-wise comp): ####
ecospat_res_all <- dplyr::bind_rows(alleqlow, allsimlow, alleqgre, allsimgre)
glimpse(ecospat_res_all)

#write.csv(ecospat_res_all, "./output/AoE/AoE_ECOSPAT_RES_allcomp.csv")

### 5. Dataframes containing only pair-wise comparisons inside each AoE ####

unique(ca_db$Area_name)

aoe_builtdf <- function(x) {
  
  aoe_spp <- ca_db %>% 
    filter(Area_name == x) %>%
    dplyr::select(Sp.1) %>%
    pull()
  
  aoe_df <- eq_res_jp %>%
    filter(Sp.1 %in% aoe_spp, Sp.2 %in% aoe_spp) %>%
    mutate(AoE_name = x)
  
  return(aoe_df)
  
}

AoE_comp_df <- purrr::map(.x = unique(ca_db$Area_name), aoe_builtdf) %>%
  bind_rows()

AoE_comp_df$D_classes <- factor(AoE_comp_df$D_classes, levels = c("Limited", "Low", "Moderate", "High", "Very_high"), ordered = TRUE)
AoE_comp_df$I_classes <- factor(AoE_comp_df$I_classes, levels = c("Limited", "Low", "Moderate", "High", "Very_high"), ordered = TRUE)
AoE_comp_df$AoE_name <- factor(AoE_comp_df$AoE_name, levels = unique(AoE_comp_df$AoE_name), ordered = FALSE)

# Including regional classification of AoEs.
AoE_regions <- read.csv("./data/biogeography/AoE_listby_Reg.csv", stringsAsFactors = FALSE)[-1] %>%
  rename(Area_name = AoE_name, Region = AoE_regs)

name_vec <- as.character(AoE_comp_df$AoE_name)
region_vec <- vector(mode = "character", length = length(name_vec))
for (i in seq_along(name_vec)) {
  
  for (j in seq_along(AoE_regions$Area_name)) {
    
    if (name_vec[i] == AoE_regions$Area_name[j]) {
      
      region_vec[i] <- AoE_regions$Region[j]
      
    }
  }
}
rm(AoE_regions)
AoE_comp_df$AoE_regs <- region_vec

glimpse(AoE_comp_df)

### Saving AoE_comp_df
#write.csv(AoE_comp_df, "./output/AoE/Table_endemicsspp_comparison_eq_lower.csv") # Ecospat EQUIVALENCE LOWER
#write.csv(AoE_comp_df, "./output/AoE/Table_endemicsspp_comparison_sim_lower.csv") # Ecospat SIMILARITY LOWER
#write.csv(AoE_comp_df, "./output/AoE/Table_endemicsspp_comparison_eq_greater.csv") # Ecospat EQUIVALENCE GREATER
#write.csv(AoE_comp_df, "./output/AoE/Table_endemicsspp_comparison_sim_greater.csv") # Ecospat EQUIVALENCE GREATER

### 5.1. Unifying dataframes from all analyses (Pair-wise comp inside AoEs): ####

# NOTE: Run this sections once all the dataframes in the section three are saved.

eqlower <- read.csv("./output/AoE/Table_endemicsspp_comparison_eq_lower.csv", stringsAsFactors = FALSE)[, -c(1:2)] # Ecospat EQUIVALENCE LOWER
simlower <- read.csv("./output/AoE/Table_endemicsspp_comparison_sim_lower.csv", stringsAsFactors = FALSE)[, -c(1:2)] # Ecospat SIMILARITY LOWER
eqgreater <- read.csv("./output/AoE/Table_endemicsspp_comparison_eq_greater.csv", stringsAsFactors = FALSE)[, -c(1:2)] # Ecospat EQUIVALENCE GREATER
simgreater <- read.csv("./output/AoE/Table_endemicsspp_comparison_sim_greater.csv", stringsAsFactors = FALSE)[, -c(1:2)] # Ecospat SIMILARITY GREATER

ecospat_res <- dplyr::bind_rows(eqlower, simlower, eqgreater, simgreater)
glimpse(ecospat_res)

#write.csv(ecospat_res, "./output/AoE/AoE_ECOSPAT_RES_areasonly.csv")