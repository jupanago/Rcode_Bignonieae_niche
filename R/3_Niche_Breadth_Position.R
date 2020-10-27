## Calculating niche position and breadth

## Breadth will be calculated as the mean distance to the centroid of each species pc scores in the space
## defined by the first two PCs.

## Postion will be taken as the explicit point defined by the medians in each Principal Component.
## This calculation will be done only for PC1 and PC2 because they are important for visualization.
## The important thing here is to be able to compare positions given the grouping by AoEs (or some 
## other trait). Therefore, the aggregation of positions (centroids) give a clue of how similar
## or different the endemic species are. 

## I will use the usedist package to do the distance measures. 
## Note that the Adonis function of the vegan package is also able to measure distance between centroids. This
## was explained in the usedist documentation. 

### Libraries

library(tidyverse)
library(usedist)
library(vegan)

### 0. Loading data ####
big_db <- read.csv("./data/occs/Localidades_Total_v14.csv")
ca_db <- read.csv("./data/biogeography/spp_listby_AoE.csv", stringsAsFactors = FALSE)[45:243, -1]
AoE_data <- read.csv("./data/biogeography/AoE_raster/AoE_Chelsa10m_bio.csv", stringsAsFactors = FALSE) # AoE Environment
locs <- read.csv("./output/localities_total_unique_ch10m_Bio(v14).csv", row.names = 1, stringsAsFactors = FALSE) # Localidades_v14
ambiente <- locs %>% dplyr::select(6:24) # Spp bioclim (Dangerous variables: Bio8,9,18,19)

### 1. Preparing PCA summary ####
load("./output/2_PCA_join_spp_AoE.rdata")
pca_summ <- summary(pca_env.sel)

# Sites pc scores
sites_sco <- as_tibble(pca_summ$sites)

# Names for spp and AoEs from the environmental databases use to produce the PCA.
sites_sco$xname <- c(locs$NAME1, AoE_data$AoE_name)
# # NOTE: This should be identical to:
# x <- read.csv("./output/2_Environment2.csv")
# sum(x$xname == sites_sco$xname ) 
# # If identical, this sum is [1] 133930

### 2. Calculating medians and centroid distances (Position) ####

### 2.1. Medians dataframe ####
# Five number summary for each principal component (Max, IQ1, Median, IQ3, max).
# For now, only the two first principal components are used. 

medians_pcs <- sites_sco %>%
  group_by(xname) %>%
  group_modify(~{
    .x %>%
      purrr::map_dfc(fivenum) %>%
      mutate(pcvar = c("min", "Q1", "median", "Q3", "max"))
  }) %>%
  pivot_wider(names_from = pcvar, values_from = c(PC1, PC2, PC3, PC4, PC5, PC6)) %>%
  #dplyr::select(xname, PC1_median, PC2_median) %>% # using only PC1 and PC2 medians
  dplyr::select(xname, contains("median")) %>% # Using all the PCs medians
  ungroup() %>%
  mutate(spp_pts = paste0("pt_", row_number())) # For usedist; this identifies each species point. 

### 2.2. AoE dataframe with medians' information (Centroid) ####

pc_median_AoE_df <- ca_db %>%
  left_join(medians_pcs, by = c("Species" = "xname"))

AoE_regs <- read.csv("./data/biogeography/AoE_listby_Reg.csv")[-1]

pc_median_AoE_df <- left_join(pc_median_AoE_df , AoE_regs, by = c("Area_name" = "AoE_name"))

### 2.2.1. Visualization of position in environmental space ####

sites_sco %>% 
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(size = 0.1, color = "gray") +
  geom_point(data = pc_median_AoE_df, aes(x = PC1_median, y = PC2_median, colour = Area_name, shape = AoE_regs), size = 1) +
  geom_vline(aes(xintercept = 0), linetype = "dashed") +
  geom_hline(aes(yintercept = 0), linetype = "dashed") +
  # geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2), arrow = arrow(length = unit(0.2,"cm")),
  #              color="blue") +
  # geom_path(data = cec_df, aes(x = x, y = y)) +
  theme_classic() +
 # theme(legend.position = "none") +
  labs(x = paste("PC1", paste0(format(as.data.frame(pca_summ$cont)[2, 1]*100, digits = 4), "%"), sep = " "),
       y = paste("PC2", paste0(format(as.data.frame(pca_summ$cont)[2, 2]*100, digits = 4), "%"), sep = " "))

### 2.3. Distance matrix between each Centroid [Defined as: (PC1_median, PC2_median) point] ####
pts_matrix <- as.matrix(pc_median_AoE_df[,c("PC1_median", "PC2_median", "PC3_median", "PC4_median", "PC5_median", "PC6_median")])
rownames(pts_matrix) <- pc_median_AoE_df$spp_pts
pts_distances <- dist(pts_matrix)


### 2.4. Distance within and between AoE's species elements ####

## NOTE: Species are identified by the pt_x name under the Item columns. 
## This table shows the distance between points within the same AoE and among different AoEs.
## The Columns Group1 and 2, indicates among which pair of areas the comparison is made.
## The column Group_comp indicates the kind of comparison made: within or within. 
## The column distance shows the distance between the species points compared as indicated by
## the columns Item1 and 2. 

groups_dist <- dist_groups(d = pts_distances, g = pc_median_AoE_df$Area_name) %>%
  mutate(Group_comp = str_extract_all(Label, pattern = "Within|Between")) %>%
  dplyr::select(1, 2, 3, 4, 7, 6)

# groups_dist  <- apply(groups_dist, 2, as.character)
# write.csv(groups_dist, file = "./output/Position/2_groups_distance_wb.csv")

### 2.5. Distance to the centroid within and between groups ####

## NOTE: This table shows the distance between each species point to the centroid of an AoE.
## Both the distance to the centroid of the AoE to which it belongs and the distance to the centroid
## of other different AoEs are calculated.
## The column Item_ID shows the name of the species. The Item_group indicates the AoE to which it belongs.
## The Centroid_Group column shows the group centroid to which the species distance is calculated.

centroid_dist_allpts <- dist_to_centroids(d = pts_distances, g = pc_median_AoE_df$Area_name) %>%
  group_by(CentroidGroup) %>%
  mutate(Item_ID = pc_median_AoE_df$Species,
         Item_group = pc_median_AoE_df$Area_name,
         Centroid_Group = str_replace(CentroidGroup, pattern = "([a-z]{4}).*\\s([a-z]{5}).*", replacement = "\\1_\\2"),
         Group_comp = ifelse(Item_group == CentroidGroup, "Within", "Between")) %>%
  ungroup() %>%
  dplyr::select(1, 4, 5, 6, 7, 3)

# write.csv(centroid_dist_allpts, "./output/Position/2_sppcentroid_to_aoes_centroid_wb.csv")

### 2.5.1. Meand distance to centroid: species centroids to AoE centroids. ####
position_df <- centroid_dist_allpts %>%
  group_by(Group_comp, Item_group, Centroid_Group) %>%
  summarise(Mean_dist_cent = mean(CentroidDistance), 
            Sd_dist_cent = sd(CentroidDistance),
            Med_dist_cent = median(CentroidDistance),
            IQR_dist_cent = IQR(CentroidDistance)) 

position_fivenum_df <- centroid_dist_allpts %>%
  group_by(Group_comp, Item_group, Centroid_Group) %>%
  dplyr::select(CentroidDistance) %>%
  group_modify(~ {
    .x %>%
      purrr::map_dfc(fivenum) %>%
      mutate(pcvar = c("min", "Q1", "median", "Q3", "max"))
  }) %>%
  pivot_wider(names_from = pcvar, values_from = CentroidDistance)

position_df <- left_join(position_df, position_fivenum_df, by = c("Group_comp", "Item_group", "Centroid_Group"))

position_df %>%
  filter(Group_comp == "Within")

# write.csv(position_df, "./output/Position/2_sppcentroid_mean_dist_to_AoEcentroid_wb.csv")

### 2.6. Distance between centroids of specific AoEs ####

## NOTE: This table show the distance between the centroid of the specified AoEs.

comb_mat <- combn(unique(pc_median_AoE_df$Area_name), 2)
aoe1_vec <- comb_mat[1, ]
aoe2_vec <- comb_mat[2, ]

aoe_btn_cent_dist <- function(pcmedians_df, item_dist_df, aoe_1, aoe_2) {
  
  library(usedist)
  
  aoe_g1 <- pcmedians_df %>%
    filter(Area_name == aoe_1) %>%
    pull(spp_pts)
  
  aoe_g2 <- pcmedians_df %>%
    filter(Area_name == aoe_2) %>%
    pull(spp_pts)
 
  bet_dist_vec <- dist_between_centroids(d = item_dist_df, idx1 = aoe_g1, idx2 = aoe_g2)
  
  return(bet_dist_vec)
  
}
# aoe_btn_cent_dist(pcmedians_df = pc_median_AoE_df, item_dist_df = pts_distances, aoe_1 = "D2CA0", aoe_2 = "D2CA24")

btn_cent_dist_vec <- purrr::pmap(.l = list(aoe_1 = aoe1_vec,
                                           aoe_2 = aoe2_vec), 
                                 .f = aoe_btn_cent_dist, 
                                 pcmedians_df = pc_median_AoE_df, 
                                 item_dist_df = pts_distances)

cent_between_dist_df <- tibble(AoE_1 = aoe1_vec, 
                               AoE_2 = aoe2_vec, 
                               Dist_between = unlist(btn_cent_dist_vec))

# write.csv(cent_between_dist_df, "./output/Position/2_aoespecific_between_centroid_dist.csv")

### 3. Breadth: preparing dataframes ####

## Format for "usedist"
# pts <- data.frame(
#   x = c(-1, 0, 0, 1, 2, 3, 3, 4), => PC1
#   y = c(0, 1, -1, 0, 0, 1, -1, 0), => PC2
#   Item = LETTERS[1:8], => spp_pts
#   Group = rep(c("Control", "Treatment"), each=4)) => xname

aoe_vec <- unique(ca_db$Area_name)

pcs_points <- sites_sco %>%
  filter(xname %in% aoe_vec) %>%
  group_by(xname) %>%
  sample_n(size = 1000, replace = TRUE) %>%
  mutate(spp_pts = paste0("pt_", row_number())) 
  

# Distance between points
pts_matrix <- as.matrix(pcs_points[,c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6")])
rownames(pts_matrix) <- pcs_points$spp_pts
pts_distances <- dist(pts_matrix)

# Distance within and between groups' elements
groups_dist <- dist_groups(d = pts_distances, g = pcs_points$xname) %>%
  mutate(Group_comp = str_extract_all(Label, pattern = "Within|Between")) %>%
  dplyr::select(1, 2, 3, 4, 7, 6)

# Distance to Centroid within and between groups
centroid_dist_allpts <- dist_to_centroids(pts_distances, pcs_points$xname) %>%
  group_by(CentroidGroup) %>%
  mutate(Item_ID = pcs_points$xname,
         Centroid_Group = str_replace(CentroidGroup, pattern = "([a-z]{4}).*\\s([a-z]{5}).*", replacement = "\\1_\\2"),
         Group_comp = ifelse(Item_ID == CentroidGroup, "Within", "Between")) %>%
  ungroup() %>%
  dplyr::select(1, 4, 5, 6, 3)

breadth_df <- centroid_dist_allpts %>%
  group_by(Group_comp, Item_ID, Centroid_Group) %>%
  summarise(meandist_centroid = mean(CentroidDistance), 
         sddist_centroid = sd(CentroidDistance),
         median_centroid = median(CentroidDistance),
         IQR_centroid = IQR(CentroidDistance)) 

breadth_fivenum_df <- centroid_dist_allpts %>%
  group_by(Group_comp, Item_ID, Centroid_Group) %>%
  dplyr::select(CentroidDistance) %>%
  group_modify(~ {
    .x %>%
      purrr::map_dfc(fivenum) %>%
      mutate(pcvar = c("min", "Q1", "median", "Q3", "max"))
  }) %>%
  pivot_wider(names_from = pcvar, values_from = CentroidDistance)
  
breadth_df <- left_join(breadth_df, breadth_fivenum_df)

# Distance between centroids of specific groups.
# Groups must be explicit by indicating their elements. 
dist_between_centroids(d = pts_distances, idx1 = pcs_points$spp_pts[1:4], idx2 = pcs_points$spp_pts[5:12])

### 3.1. Species breadth function ####

pcs_points <- sites_sco %>%
  mutate(spp_pts = paste0("pt_", row_number()))

# Species with only one point. No breadth was calculated for these species. 
one_point_only <- pcs_points %>%
  group_by(xname) %>%
  summarise(n = length(xname)) %>%
  filter(n == 1) %>%
  pull(xname)

spp_distto_cent <- function(pcs_df, g_name) {
  
  library(usedist)
  library(tidyr)
  library(dplyr)
  library(purrr)
  
  subset_pts <- pcs_df %>%
    filter(xname == g_name)
  
  pts_matrix <- as.matrix(subset_pts[,1:6])
  rownames(pts_matrix) <- subset_pts$spp_pts
  pts_distances <- dist(pts_matrix)
  
  cent_dist_df <- dist_to_centroids(pts_distances, subset_pts$xname) %>%
    group_by(CentroidGroup) %>%
    mutate(Item_ID = subset_pts$xname,
           Centroid_Group = str_replace(CentroidGroup, pattern = "([a-z]{4}).*\\s([a-z]{5}).*", replacement = "\\1_\\2"),
           Group_comp = ifelse(Item_ID == CentroidGroup, "Within", "Between")) %>%
    ungroup() %>%
    dplyr::select(1, 4, 5, 6, 3)
  
  breadth_df <- cent_dist_df %>%
    group_by(Group_comp, Item_ID, Centroid_Group) %>%
    summarise(Mean_dist_cent = mean(CentroidDistance), 
              Sd_dist_cent = sd(CentroidDistance),
              Med_dist_cent = median(CentroidDistance),
              IQR_dist_cent = IQR(CentroidDistance)) 
  
  breadth_fivenum_df <- cent_dist_df %>%
    group_by(Group_comp, Item_ID, Centroid_Group) %>%
    dplyr::select(CentroidDistance) %>%
    group_modify(~ {
      .x %>%
        purrr::map_dfc(fivenum) %>%
        mutate(pcvar = c("min", "Q1", "median", "Q3", "max"))
    }) %>%
    pivot_wider(names_from = pcvar, values_from = CentroidDistance)
  
  breadth_df <- left_join(breadth_df, breadth_fivenum_df, by = c("Group_comp", "Item_ID", "Centroid_Group"))
  
  return(breadth_df)
  
}

g_names_vec <- unique(pcs_points$xname) 
g_names_vec <- g_names_vec[!g_names_vec %in% one_point_only]

spp_distto_cent(pcs_df = pcs_points, g_name = g_names_vec[365])

# This was calculated in Azure. Firts for speceis only, and then for AoEs. AoEs didn't work
# a <- Sys.time()
# mean_dist_list <- purrr::map(g_names_vec, spp_distto_cent, pcs_df = pcs_points) 
# spp_centr_meandist_df <- do.call("rbind", mean_dist_list)
# b <- Sys.time()
# a - b

# From AZURE speceis niche bradth
breadth_spp_df <- read.csv("./output/Breadth/2_Breadth_AoEspp_joint_AZURE.csv")

### Reverse Engeneering: ####

d <- pts_distances
g <- pcs_points$xname
squared <- FALSE 

d <- stats::as.dist(d)
g <- as.factor(g)
dsize <- attr(d, "Size")
if (length(g) != dsize) {
  stop("Length of grouping vector (g) must equal number of observations in ", 
       "dist object (d)")
}
dlabels <- attr(d, "Labels")
idxs <- utils::combn(dsize, 2)
idx1 <- idxs[1, ]
idx2 <- idxs[2, ]
level1 <- levels(g)[pmin(as.numeric(g[idx1]), as.numeric(g[idx2]))]
level2 <- levels(g)[pmax(as.numeric(g[idx1]), as.numeric(g[idx2]))]
data.frame(Item1 = if (is.null(dlabels)) 
  idx1
  else dlabels[idx1], Item2 = if (is.null(dlabels)) 
    idx2
  else dlabels[idx2], Group1 = g[idx1], Group2 = g[idx2], 
  Label = factor(ifelse(level1 == level2, "Within", "Between")), 
  Distance = dist_get(d, idx1, idx2), stringsAsFactors = FALSE)
