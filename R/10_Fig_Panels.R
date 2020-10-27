### Script for drawing panels of niche properties among AoEs Regional Environment 

library(patchwork)

## 1. Map of the AoE ####

library(tidyverse)
library(sp)
library(rgdal)
library(gridExtra)
library(grid)

## 1.1. Loading  data ####

AoE_regs <- read.csv("./output/AoE/AoE_ECOSPAT_RES_areasonly.csv", stringsAsFactors = FALSE) %>%  
  filter(Analysis == "Eq", Alternative == "G") %>%
  dplyr::select(AoE_name, AoE_regs) %>%
  distinct() %>%
  mutate(region_name = case_when(
    AoE_regs == "ESA_THDD" ~ "ESA-THDD", 
    AoE_regs == "ESA_Atlantic_Forest" ~ "ESA-Atlantic Forest", 
    AoE_regs == "ESA_Dry_Diagonal" ~ "ESA-Dry Diagonal", 
    AoE_regs == "Guiana_centered" ~ "Guiana Centered", 
    AoE_regs == "Amazonia_centered" ~ "Amazonia Centered", 
    AoE_regs == "Mesoamerica" ~ "Mesoamerica", 
    AoE_regs == "Northern_Andes" ~ "Northern Andes"
  )) 

ca_db <- read.csv("./data/biogeography/spp_listby_AoE.csv", stringsAsFactors = FALSE)[45:243,] %>%
  dplyr::select(-1) %>%
  rename(AoE_name = Area_name) %>%
  group_by(AoE_name) %>%
  summarise(Number_spp = n()) %>%
  mutate(id = as.numeric(str_replace(AoE_name, pattern = "D2CA([0-9]+)", replacement = "\\1"))) %>%
  arrange(id) %>%
  dplyr::select(-3)

aoe_info_df <- left_join(ca_db, AoE_regs, by = "AoE_name")
rm(AoE_regs, ca_db)

layers_regions_ca <- list(
  Mesoamerica = c("2dgD_CA12_NWYuc-grid"),
  Northern_Andes = c("2dgD_CA14_NWColVen-grid", "2dgD_CA17_NWCol-grid", "2dgD_CA22_NWCol2-grid", "2dgD_CA26_NWCol3-grid"),
  Guiana_centred = c("2dgD_CA3_AM-grid", "2dgD_CA10_AM-grid", "2dgD_CA25_AM-grid"),
  Amazonia_centred = c("2dgD_CA4_AM-grid", "2dgD_CA7_AM-grid", "2dgD_CA16_AM-grid", "2dgD_CA18_AM-grid", "2dgD_CA20_AM-grid", "2dgD_CA23_AM-grid"),
  ESA_DD = c("2dgD_CA2_DD-grid","2dgD_CA5_DD-grid", "2dgD_CA19_DD-grid"),
  ESA_THDD = c("2dgD_CA0_AF-grid", "2dgD_CA6_AFDD-grid", "2dgD_CA8_AF-grid", "2dgD_CA9_AF-grid", "2dgD_CA11_THDD-grid", "2dgD_CA27_DDb-grid"),
  ESA_Atlantic = c("2dgD_CA1_AF-grid", "2dgD_CA13_AF-grid", "2dgD_CA15_AF-grid", "2dgD_CA21_AF-grid", "2dgD_CA24_AF-grid")
)


## 1.2. Creating AoE polygons in ggplot format ####

ndm_geompol <- function(shp_dns, layer_names, alfa = 0.3, fill_ca = "viridis", colour_ca = "black", size_border = 0.8) {
  
  # Libraries
  library(ggplot2)
  library(viridis)
  library(wesanderson)
  library(rgdal)
  library(sf)
  
  # fortyfing function
  fortifying_ca <- function(x) {
    x <- st_as_sf(x)
    ca_pol <- st_union(x)
    ca_pol <- as_Spatial(ca_pol)
    ca_pol_df <- fortify(ca_pol)
    return(ca_pol_df)
  }
  
  # naming function
  naming_ca <- function(x) {
    
    if (grepl(x, pattern = ".*1dgD.*")) {
      x <- gsub(x, pattern = ".*1dgD_(.*)_.*-grid", replacement = "\\1")
      y <- paste("D1", x, sep = "")
      return(y)
    }
    
    if (grepl(x, pattern = ".*2dgD.*")) {
      x <- gsub(x, pattern = ".*2dgD_(.*)_.*-grid", replacement = "\\1")
      y <- paste("D2", x, sep = "")
      return(y)
    }
    
    if (grepl(x, pattern = ".*3dgD.*")) {
      x <- gsub(x, pattern = ".*3dgD_(.*)_.*-grid", replacement = "\\1")
      y <- paste("D3", x, sep = "")
      return(y)
    }
    
    if (grepl(x, pattern = ".*1dgT.*")) {
      x <- gsub(x, pattern = ".*1dgT_(.*)_.*-grid", replacement = "\\1")
      y <- paste("T1", x, sep = "")
      return(y)
    }
    
    if (grepl(x, pattern = ".*2dgT.*")) {
      x <- gsub(x, pattern = ".*2dgT_(.*)_.*-grid", replacement = "\\1")
      y <- paste("T2", x, sep = "")
      return(y)
    }
    
    if (grepl(x, pattern = ".*3dgT.*")) {
      x <- gsub(x, pattern = ".*3dgT_(.*)_.*-grid", replacement = "\\1")
      y <- paste("T3", x, sep = "")
      return(y)
    }
    
    if (grepl(x, pattern = ".*1dgS.*")) {
      x <- gsub(x, pattern = ".*1dgS_(.*)_.*-grid", replacement = "\\1")
      y <- paste("S1", x, sep = "")
      return(y)
    }
    
    if (grepl(x, pattern = ".*2dgS.*")) {
      x <- gsub(x, pattern = ".*2dgD_(.*)_.*-grid", replacement = "\\1")
      y <- paste("D2", x, sep = "")
      return(y)
    }
    
    if (grepl(x, pattern = ".*3dgS.*")) {
      x <- gsub(x, pattern = ".*3dgS_(.*)_.*-grid", replacement = "\\1")
      y <- paste("S3", x, sep = "")
      return(y)
    }
    
    if (grepl(x, pattern = ".*1dgH.*")) {
      x <- gsub(x, pattern = ".*1dgH_(.*)_.*-grid", replacement = "\\1")
      y <- paste("H1", x, sep = "")
      return(y)
    }
    
    if (grepl(x, pattern = ".*2dgH.*")) {
      x <- gsub(x, pattern = ".*2dgH_(.*)_.*-grid", replacement = "\\1")
      y <- paste("H2", x, sep = "")
      return(y)
    }
    
    if (grepl(x, pattern = ".*3dgH.*")) {
      x <- gsub(x, pattern = ".*3dgH_(.*)_.*-grid", replacement = "\\1")
      y <- paste("H3", x, sep = "")
      return(y)
    }
    
  }
  
  # Defining projection
  wgs84.proj4 <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
  
  # Defining color pallette
  
  if (fill_ca == "viridis") {
    
    fillingv_ca <- viridis(n = length(layer_names), alpha = alfa)
    
  }
  
  if (fill_ca == "cavalcanti") {
    
    #jpn_pal <- wes_palette("Cavalcanti1", type = "continuous", n = 7)
    jpn_pal <- c("#D8B70A", "#496715", "#376138", "#A2A475", "#8CA685", "#887F65", "#972D15")
    rep_n <- lapply(layers_regions_ca, length)
    
    fillingv_ca <- c(
      rep(jpn_pal[1], rep_n[[1]]),
      rep(jpn_pal[2], rep_n[[2]]),
      rep(jpn_pal[3], rep_n[[3]]),
      rep(jpn_pal[4], rep_n[[4]]),
      rep(jpn_pal[5], rep_n[[5]]),
      rep(jpn_pal[6], rep_n[[6]]),
      rep(jpn_pal[7], rep_n[[7]])
    )
    
    rm(rep_n, jpn_pal)
  }
  
  if (fill_ca == "darjeeling") {
    
    #jpn_pal <- wes_palette("Darjeeling1", type = "continuous", n = 7)
    jpn_pal <- c("#FF0000", "#556A5B", "#50A45C", "#F2AD00", "#6D5B4E", "#C49647", "#5BBCD6")
    rep_n <- lapply(layers_regions_ca, length)
    
    fillingv_ca <- c(
      rep(jpn_pal[1], rep_n[[1]]),
      rep(jpn_pal[2], rep_n[[2]]),
      rep(jpn_pal[3], rep_n[[3]]),
      rep(jpn_pal[4], rep_n[[4]]),
      rep(jpn_pal[5], rep_n[[5]]),
      rep(jpn_pal[6], rep_n[[6]]),
      rep(jpn_pal[7], rep_n[[7]])
    )
    
    rm(rep_n, jpn_pal)
  }
  
  
  # Initiating list
  gg_ca <- vector(mode = "list", length = length(layer_names))
  
  # Creating geoms list
  for (i in seq_along(layer_names)) {
    
    ca <- readOGR(dsn = shp_dns, layer = layer_names[i], p4s = wgs84.proj4)
    ca_poldf <- fortifying_ca(ca)
    
    gg_ca[[i]][[1]] <- geom_polygon(data = ca_poldf, 
                                    aes(x = long, y = lat, group = group),
                                    fill = ifelse(fill_ca == "viridis", fillingv_ca[i],
                                                  ifelse(fill_ca == "cavalcanti", fillingv_ca[i], 
                                                         ifelse(fill_ca == "darjeeling", fillingv_ca[i], fill_ca))),
                                    colour = colour_ca,
                                    size = size_border,
                                    alpha = alfa)
    gg_ca[[i]][[2]] <- naming_ca(layer_names[i])
    
  }
  
  return(gg_ca)
  
}

x <- ndm_geompol(shp_dns = "./data/biogeography/AoE_shapes/", layer_names = unlist(layers_regions_ca), fill_ca = "darjeeling", alfa = 0.3)
# y <- ndm_geompol(shp_dns = "./data/biogeography/AoE_shapes/", layer_names = unlist(layers_regions_ca)[1], fill_ca = "cavalcanti", alfa = 0.8)
# 
# pattinfo <- aoe_info_df %>%
#   filter(AoE_name == y[[1]][[2]]) %>%
#   unlist()
# 
# base_map + 
#   y[[1]][[1]] +
#   theme(plot.title = element_text(size = 15, face = "bold")) +
#   labs(title = y[[1]][[2]],
#        subtitle = paste(paste("Geographic region:", pattinfo[5], sep = " "), "\n",  
#                         paste("Pattern type:", pattinfo[2], "\nSpp:", pattinfo[3]), sep = ""))

## 1.3. Creating panel maps ####
wgs84.proj4 <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
Amer <- readOGR(dsn = "./data/GIS/", layer = "Amer_land", p4s = wgs84.proj4)
Amer_df <- fortify(Amer)

base_map <- ggplot(Amer_df, aes(x = long, y = lat, group = group)) + 
  geom_polygon(fill = "white", 
               col = "black") +
  coord_cartesian(xlim = c(-111, -34),
                  ylim = c(-37, 39),
                  expand = TRUE) +
  labs(x = NULL, y = NULL) +
  theme_linedraw() +
  theme(panel.grid = element_blank())

maps_list <- vector(mode = "list", length = length(x))

for (i in seq_along(maps_list)) {
  
  pattinfo <- aoe_info_df %>%
    filter(AoE_name == x[[i]][[2]]) %>%
    unlist()
  
  maps_list[[i]] <- base_map + 
    x[[i]][[1]] +
    theme(plot.title = element_text(size = 15, face = "bold")) +
    labs(title = x[[i]][[2]],
         subtitle = paste(paste("GR:", pattinfo[4], sep = " "), "\n",  
                          paste("Spp:", pattinfo[2]), sep = ""))
  
}

#dir.create("./figs/Panel_figs")
#saveRDS(maps_list, file = "./figs/Panel_figs/A_AoEmaps.rds")

maps_list <- readRDS(file = "./figs/Panel_figs/A_AoEmaps.rds")
maps_list[[1]] 

### 1.3.1. Changing area names code ####
#NOTE: From 7_Overlap_stats.r

AoE_comp_df <- read.csv("./output/AoE/AoE_ECOSPAT_RES_areasonly.csv") %>%  filter(Analysis == "Eq", Alternative == "G")

AoE_comp_df$AoE_regs <- factor(AoE_comp_df$AoE_regs, levels = c("Mesoamerica", "Northern_Andes", "Guiana_centered", "Amazonia_centered", 
                                                                "ESA_Dry_Diagonal", "ESA_THDD", "ESA_Atlantic_Forest"), ordered = TRUE)
AoE_comp_df <- AoE_comp_df %>% arrange(AoE_regs)

new_area_code <- data.frame(unique(cbind(as.character(AoE_comp_df$AoE_name), as.character(AoE_comp_df$AoE_regs))))
colnames(new_area_code) <- c("aoe", "reg")
new_area_code$reg <- as.character(new_area_code$reg)
new_area_code$aoe <- as.character(new_area_code$aoe)

area_code <- new_area_code %>%
  mutate(reg2 = case_when(
    reg == "Mesoamerica" ~ "MesoAme",
    reg == "Northern_Andes" ~ "N_Andes",
    reg == "Guiana_centered" ~ "Guiana",
    reg == "Amazonia_centered" ~ "Amazonia",
    reg == "ESA_Dry_Diagonal" ~ "D_Diagonal",
    reg == "ESA_THDD" ~ "DD_AF",
    reg == "ESA_Atlantic_Forest" ~ "Atl_Forest"
  )) %>%
  group_by(reg2) %>%
  add_tally() %>%
  mutate(aoen = 1:n,
         areacode = paste(reg2, aoen, sep = "_")) %>%
  pull(areacode)

maps_list[[13]]

for (i in seq_along(area_code)) {
  
  maps_list[[i]] <- maps_list[[i]] +
    labs(title = area_code[i])
  
}

#saveRDS(maps_list, file = "./figs/Panel_figs/A_AoEmaps.rds")

## 1.4. Map of all areas ####

# Labels
test_labels <- c("Mesoamerica", "Northern Andes", "Guiana centred", "Amazonia centred",
                 "Eastern South America: Dry Diagonal", "Eastern South America: Throughout Dry Diagonal", "Eastern South America: Atlantic Forest" )

legend_position <- sort(seq(0.5, 1, 0.07142857)[-1], decreasing = TRUE)

text_legend <- textGrob(label = test_labels, x = 0.6, 
                        y = legend_position, just = "left",
                        gp = gpar(fontsize = 10), check = TRUE) 

rect_position_y <- legend_position
rect_position_x <- rep(0.55, length(rect_position_y))

# jpn_palette <- brewer.pal(n = 7, name = "Accent")
jpn_palette <- c("#FF0000", "#556A5B", "#50A45C", "#F2AD00", "#6D5B4E", "#C49647", "#5BBCD6")

square_legend <- rectGrob(x = rect_position_x, y = rect_position_y, 
                          width = 0.05, height = 0.05,
                          gp = gpar(fill = jpn_palette, alpha = 0.2, col = jpn_palette, lwd = 2))

grob_legend <- grobTree(text_legend, square_legend, name = "jp_legendgrob")

ndm_geompol2 <- function(shp_dns, layer_names, alfa = 0.3, fill_ca = "viridis", colour_ca = "black", size_border = 0.3) {
  
  # Libraries
  library(ggplot2)
  library(viridis)
  library(wesanderson)
  library(rgdal)
  library(sf)
  
  # fortyfing function
  fortifying_ca <- function(x) {
    x <- st_as_sf(x)
    ca_pol <- st_union(x)
    ca_pol <- as_Spatial(ca_pol)
    ca_pol_df <- fortify(ca_pol)
    return(ca_pol_df)
  }
  
  # naming function
  naming_ca <- function(x) {
    
    if (grepl(x, pattern = ".*1dgD.*")) {
      x <- gsub(x, pattern = ".*1dgD_(.*)_.*-grid", replacement = "\\1")
      y <- paste("D1", x, sep = "")
      return(y)
    }
    
    if (grepl(x, pattern = ".*2dgD.*")) {
      x <- gsub(x, pattern = ".*2dgD_(.*)_.*-grid", replacement = "\\1")
      y <- paste("D2", x, sep = "")
      return(y)
    }
    
    if (grepl(x, pattern = ".*3dgD.*")) {
      x <- gsub(x, pattern = ".*3dgD_(.*)_.*-grid", replacement = "\\1")
      y <- paste("D3", x, sep = "")
      return(y)
    }
    
    if (grepl(x, pattern = ".*1dgT.*")) {
      x <- gsub(x, pattern = ".*1dgT_(.*)_.*-grid", replacement = "\\1")
      y <- paste("T1", x, sep = "")
      return(y)
    }
    
    if (grepl(x, pattern = ".*2dgT.*")) {
      x <- gsub(x, pattern = ".*2dgT_(.*)_.*-grid", replacement = "\\1")
      y <- paste("T2", x, sep = "")
      return(y)
    }
    
    if (grepl(x, pattern = ".*3dgT.*")) {
      x <- gsub(x, pattern = ".*3dgT_(.*)_.*-grid", replacement = "\\1")
      y <- paste("T3", x, sep = "")
      return(y)
    }
    
    if (grepl(x, pattern = ".*1dgS.*")) {
      x <- gsub(x, pattern = ".*1dgS_(.*)_.*-grid", replacement = "\\1")
      y <- paste("S1", x, sep = "")
      return(y)
    }
    
    if (grepl(x, pattern = ".*2dgS.*")) {
      x <- gsub(x, pattern = ".*2dgD_(.*)_.*-grid", replacement = "\\1")
      y <- paste("D2", x, sep = "")
      return(y)
    }
    
    if (grepl(x, pattern = ".*3dgS.*")) {
      x <- gsub(x, pattern = ".*3dgS_(.*)_.*-grid", replacement = "\\1")
      y <- paste("S3", x, sep = "")
      return(y)
    }
    
    if (grepl(x, pattern = ".*1dgH.*")) {
      x <- gsub(x, pattern = ".*1dgH_(.*)_.*-grid", replacement = "\\1")
      y <- paste("H1", x, sep = "")
      return(y)
    }
    
    if (grepl(x, pattern = ".*2dgH.*")) {
      x <- gsub(x, pattern = ".*2dgH_(.*)_.*-grid", replacement = "\\1")
      y <- paste("H2", x, sep = "")
      return(y)
    }
    
    if (grepl(x, pattern = ".*3dgH.*")) {
      x <- gsub(x, pattern = ".*3dgH_(.*)_.*-grid", replacement = "\\1")
      y <- paste("H3", x, sep = "")
      return(y)
    }
    
  }
  
  # Defining projection
  wgs84.proj4 <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
  
  # Defining color pallette
  
  if (fill_ca == "viridis") {
    
    fillingv_ca <- viridis(n = length(layer_names), alpha = alfa)
    
  }
  
  if (fill_ca == "cavalcanti") {
    
    #jpn_pal <- wes_palette("Cavalcanti1", type = "continuous", n = 7)
    jpn_pal <- c("#D8B70A", "#496715", "#376138", "#A2A475", "#8CA685", "#887F65", "#972D15")
    rep_n <- lapply(layers_regions_ca, length)
    
    fillingv_ca <- c(
      rep(jpn_pal[1], rep_n[[1]]),
      rep(jpn_pal[2], rep_n[[2]]),
      rep(jpn_pal[3], rep_n[[3]]),
      rep(jpn_pal[4], rep_n[[4]]),
      rep(jpn_pal[5], rep_n[[5]]),
      rep(jpn_pal[6], rep_n[[6]]),
      rep(jpn_pal[7], rep_n[[7]])
    )
    
    rm(rep_n, jpn_pal)
  }
  
  if (fill_ca == "darjeeling") {
    
    #jpn_pal <- wes_palette("Darjeeling1", type = "continuous", n = 7)
    jpn_pal <- c("#FF0000", "#556A5B", "#50A45C", "#F2AD00", "#6D5B4E", "#C49647", "#5BBCD6")
    rep_n <- lapply(layers_regions_ca, length)
    
    fillingv_ca <- c(
      rep(jpn_pal[1], rep_n[[1]]),
      rep(jpn_pal[2], rep_n[[2]]),
      rep(jpn_pal[3], rep_n[[3]]),
      rep(jpn_pal[4], rep_n[[4]]),
      rep(jpn_pal[5], rep_n[[5]]),
      rep(jpn_pal[6], rep_n[[6]]),
      rep(jpn_pal[7], rep_n[[7]])
    )
    
    rm(rep_n, jpn_pal)
  }
  
  
  # Initiating list
  gg_ca <- vector(mode = "list", length = length(layer_names))
  
  # Creating geoms list
  for (i in seq_along(layer_names)) {
    
    ca <- readOGR(dsn = shp_dns, layer = layer_names[i], p4s = wgs84.proj4)
    ca_poldf <- fortifying_ca(ca)
    
    gg_ca[[i]][[1]] <- geom_polygon(data = ca_poldf, 
                                    aes(x = long, y = lat, group = group),
                                    fill = ifelse(fill_ca == "viridis", fillingv_ca[i],
                                                  ifelse(fill_ca == "cavalcanti", fillingv_ca[i], 
                                                         ifelse(fill_ca == "darjeeling", fillingv_ca[i], fill_ca))),
                                    colour = ifelse(colour_ca == "black", black, fillingv_ca[i]),
                                    size = size_border,
                                    alpha = alfa)
    gg_ca[[i]][[2]] <- naming_ca(layer_names[i])
    
  }
  
  return(gg_ca)
  
}
#x <- ndm_geompol2(shp_dns = "./data/biogeography/AoE_shapes/", layer_names = unlist(layers_regions_ca), fill_ca = "darjeeling", alfa = 0.2, colour_ca = "darjeeling")

all_areas_plot <- base_map +
  x[[1]][[1]] +
  x[[10]][[1]] +
  x[[11]][[1]] +
  x[[12]][[1]] +
  x[[13]][[1]] +
  x[[14]][[1]] +
  x[[15]][[1]] +
  x[[2]][[1]] +
  x[[3]][[1]] +
  x[[4]][[1]] +
  x[[5]][[1]] +
  x[[6]][[1]] +
  x[[7]][[1]] +
  x[[8]][[1]] +
  x[[9]][[1]] +
  x[[16]][[1]] +
  x[[17]][[1]] +
  x[[18]][[1]] +
  x[[19]][[1]] +
  x[[20]][[1]] +
  x[[21]][[1]] +
  x[[22]][[1]] +
  x[[23]][[1]] +
  x[[24]][[1]] +
  x[[25]][[1]] +
  x[[26]][[1]] +
  x[[27]][[1]] +
  x[[28]][[1]]  +
  annotation_custom(grob = grob_legend, 
                    xmin = -150, xmax = -80, 
                    ymin = -75, ymax = -10) +
  labs(title = "Areas of endemism of Bignonieae",
       subtitle = "Approximated geographic regions") +
  theme(plot.title = element_text(face = "bold"))

png("./figs/ALL_AoEs_map.png", width = 1100, height = 1000, res = 130)
print(all_areas_plot)
dev.off()

## 2. AoE hull in environmental space Regional ####

# Libraries
library(tidyverse)
library(vegan)
library(ggordiplots)

## 2.1. Load data ####
ca_db <- read.csv("./data/biogeography/spp_listby_AoE.csv", stringsAsFactors = FALSE)[45:243,]
load("./output/2_PCA_join_spp_AoE.rdata")
comp_env <- read.csv("./output/2_Environment2.csv")

## 2.2. Defining Color palette ####
layers_regions_ca <- list(
  Mesoamerica = c("2dgD_CA12_NWYuc-grid"),
  Northern_Andes = c("2dgD_CA14_NWColVen-grid", "2dgD_CA17_NWCol-grid", "2dgD_CA22_NWCol2-grid", "2dgD_CA26_NWCol3-grid"),
  Guiana_centred = c("2dgD_CA3_AM-grid", "2dgD_CA10_AM-grid", "2dgD_CA25_AM-grid"),
  Amazonia_centred = c("2dgD_CA4_AM-grid", "2dgD_CA7_AM-grid", "2dgD_CA16_AM-grid", "2dgD_CA18_AM-grid", "2dgD_CA20_AM-grid", "2dgD_CA23_AM-grid"),
  ESA_DD = c("2dgD_CA2_DD-grid","2dgD_CA5_DD-grid", "2dgD_CA19_DD-grid"),
  ESA_THDD = c("2dgD_CA0_AF-grid", "2dgD_CA6_AFDD-grid", "2dgD_CA8_AF-grid", "2dgD_CA9_AF-grid", "2dgD_CA11_THDD-grid", "2dgD_CA27_DDb-grid"),
  ESA_Atlantic = c("2dgD_CA1_AF-grid", "2dgD_CA13_AF-grid", "2dgD_CA15_AF-grid", "2dgD_CA21_AF-grid", "2dgD_CA24_AF-grid")
) %>% 
  lapply(str_replace_all, pattern = ".*_(CA[0-9]+)_.*", replacement = "D2\\1")

# From wesanderson package "Darjeeling1" for 7 colors. Fifth color from the same palette but with n = 9
jpn_pal <- c("#FF0000", "#556A5B", "#50A45C", "#F2AD00", "#6D5B4E", "#C49647", "#5BBCD6")
rep_n <- lapply(layers_regions_ca, length)

fillingv_ca <- c(
  rep(jpn_pal[1], rep_n[[1]]),
  rep(jpn_pal[2], rep_n[[2]]),
  rep(jpn_pal[3], rep_n[[3]]),
  rep(jpn_pal[4], rep_n[[4]]),
  rep(jpn_pal[5], rep_n[[5]]),
  rep(jpn_pal[6], rep_n[[6]]),
  rep(jpn_pal[7], rep_n[[7]])
)

rm(rep_n, jpn_pal)

## 2.3. Intermediate objects necessary to produce the plots ####

# Summary to obtain percentage of variation explained by PCs
pca_summ <- summary(pca_env.sel)
# AoE names' vector
aoe_name_vec <- unique(ca_db$Area_name)
# Points in environmental space
sites_sco <- scores(pca_env.sel, display = "sites", scaling = 2)
sites_sco <- as_tibble(sites_sco)
sites_sco$xname <- comp_env$xname

# Calculating hull in environmetal space (PC1 x PC2) 
hull <- sites_sco %>%
  filter(xname %in% aoe_name_vec[1]) %>%
  group_by(xname) %>%
  slice(chull(PC1, PC2))

# Ellipse (Just in case it is required)
# ellip <- sites_sco %>%
#   filter(xname %in% aoe_name_vec)

## 2.4. Plot of environmental space ####

sites_sco %>%
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(size = 0.1, alpha = 0.5, color = "gray") +
  geom_vline(aes(xintercept = 0), linetype = "dashed") +
  geom_hline(aes(yintercept = 0), linetype = "dashed") +
  geom_polygon(data = hull, fill = fillingv_ca[1], colour = "black", alpha = 0.3, size = 0.2) +
  coord_cartesian(xlim = c(-0.3748672, 0.4475069), ylim = c(-0.4032099, 0.1944293)) +
  geom_vline(aes(xintercept = 0), linetype = "dashed") +
  geom_hline(aes(yintercept = 0), linetype = "dashed") +         
  #stat_ellipse(data = ellip, aes(color = factor(xname)), level = 0.95) +
  theme_classic() +
  labs(title = paste("Environmental space of", "xname", sep = " "),
       x = paste("PC1", paste0(format(as.data.frame(pca_summ$cont)[2, 1]*100, digits = 4), "%"), sep = " "),
       y = paste("PC2", paste0(format(as.data.frame(pca_summ$cont)[2, 2]*100, digits = 4), "%"), sep = " ")) + 
  theme(legend.position = "none",
        plot.title = element_text(size = 12, face = "bold"))

## 2.5. Function to plot the environment ####

env_plot <- function(aoe_pca, aoe_env, aoe_name, env_fill, alfa = 0.3) {
  
  pca_summ <- summary(aoe_pca)
  # AoE names' vector
  aoe_name_vec <- unique(ca_db$Area_name)
  # Points in environmental space
  sites_sco <- scores(pca_env.sel, display = "sites", scaling = 2)
  sites_sco <- as_tibble(sites_sco)
  sites_sco$xname <- aoe_env$xname
  
  hull <- sites_sco %>%
    filter(xname %in% aoe_name) %>%
    group_by(xname) %>%
    slice(chull(PC1, PC2))
  
  aoeenvplot <- sites_sco %>%
    ggplot(aes(x = PC1, y = PC2)) +
    geom_point(size = 0.1, alpha = 0.5, color = "gray") +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_hline(aes(yintercept = 0), linetype = "dashed") +
    geom_polygon(data = hull, fill = env_fill, colour = "black", alpha = alfa, size = 0.2) +
    coord_cartesian(xlim = c(-0.3748672, 0.4475069), ylim = c(-0.4032099, 0.1944293)) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_hline(aes(yintercept = 0), linetype = "dashed") +         
    #stat_ellipse(data = ellip, aes(color = factor(xname)), level = 0.95) +
    theme_classic() +
    labs(subtitle = "Environmental space",
         x = paste("PC1", paste0(format(as.data.frame(pca_summ$cont)[2, 1]*100, digits = 4), "%"), sep = " "),
         y = paste("PC2", paste0(format(as.data.frame(pca_summ$cont)[2, 2]*100, digits = 4), "%"), sep = " ")) + 
    theme(legend.position = "none")
  
  return(aoeenvplot)
  
}

aoe_name_vec <- unlist(layers_regions_ca)

env_plot(aoe_pca = pca_env.sel, aoe_env = comp_env, aoe_name = aoe_name_vec[1], env_fill = fillingv_ca[1], alfa = 0.3)

env_list <- purrr::map2(.x = aoe_name_vec, .y = fillingv_ca, env_plot, aoe_pca = pca_env.sel, aoe_env = comp_env, alfa = 0.3)

#dir.create("./figs/paper/Panel_regionalenv")
saveRDS(env_list, file = "./figs/Panel_figs/B_AoEenv.rds")

env_list  <- readRDS("./figs/paper/Panel_regionalenv/B_AoEenv.rds")
env_list[2]

# library(patchwork)
i <- 1
(maps_list[[i]] | env_list[[i]])  + plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face = "bold"))

## 2.6 Plot of all environments at once ####

pca_summ <- summary(pca_env.sel)
sites_sco <- scores(pca_env.sel, display = "sites", scaling = 2)
sites_sco <- as_tibble(sites_sco)

env_basemap <-  sites_sco %>%
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(size = 0.1, alpha = 0.5, color = "gray") +
  geom_vline(aes(xintercept = 0), linetype = "dashed") +
  geom_hline(aes(yintercept = 0), linetype = "dashed") +
  #coord_cartesian(xlim = c(-0.3748672, 0.4475069), ylim = c(-0.4032099, 0.1944293)) +
  geom_vline(aes(xintercept = 0), linetype = "dashed") +
  geom_hline(aes(yintercept = 0), linetype = "dashed") +         
  #stat_ellipse(data = ellip, aes(color = factor(xname)), level = 0.95) +
  theme_classic() +
  labs(title = "Regional environmental space",
       x = paste("PC1", paste0(format(as.data.frame(pca_summ$cont)[2, 1]*100, digits = 4), "%"), sep = " "),
       y = paste("PC2", paste0(format(as.data.frame(pca_summ$cont)[2, 2]*100, digits = 4), "%"), sep = " ")) + 
  theme(legend.position = "none", 
        plot.title = element_text(face = "bold"))

# Labels
test_labels <- c("Mesoamerica", "Northern Andes", "Guiana centred", "Amazonia centred",
                 "ESA: Dry Diagonal", "ESA: THDD", "ESA: Atlantic Forest" )

legend_position <- sort(seq(0.5, 1, 0.07142857)[-1], decreasing = TRUE)

text_legend <- textGrob(label = test_labels, x = 0.6, 
                        y = legend_position, just = "left",
                        gp = gpar(fontsize = 10), check = TRUE) 

rect_position_y <- legend_position
rect_position_x <- rep(0.55, length(rect_position_y))

# jpn_palette <- brewer.pal(n = 7, name = "Accent")
jpn_palette <- c("#FF0000", "#556A5B", "#50A45C", "#F2AD00", "#6D5B4E", "#C49647", "#5BBCD6")

square_legend <- rectGrob(x = rect_position_x, y = rect_position_y, 
                          width = 0.05, height = 0.05,
                          gp = gpar(fill = jpn_palette, alpha = 0.2, col = jpn_palette, lwd = 2))

grob_legend <- grobTree(text_legend, square_legend, name = "jp_legendgrob")


env_geomhull <- function (aoe_pca, aoe_env, aoe_name, fill_ca, colour_ca = "black", alfa = 0.3, size_border = 0.2) {
  
  # Libraries
  library(ggplot2)
  library(viridis)
  library(wesanderson)
  library(vegan)
  
  layers_regions_ca <- list(
    Mesoamerica = c("2dgD_CA12_NWYuc-grid"),
    Northern_Andes = c("2dgD_CA14_NWColVen-grid", "2dgD_CA17_NWCol-grid", "2dgD_CA22_NWCol2-grid", "2dgD_CA26_NWCol3-grid"),
    Guiana_centred = c("2dgD_CA3_AM-grid", "2dgD_CA10_AM-grid", "2dgD_CA25_AM-grid"),
    Amazonia_centred = c("2dgD_CA4_AM-grid", "2dgD_CA7_AM-grid", "2dgD_CA16_AM-grid", "2dgD_CA18_AM-grid", "2dgD_CA20_AM-grid", "2dgD_CA23_AM-grid"),
    ESA_DD = c("2dgD_CA2_DD-grid","2dgD_CA5_DD-grid", "2dgD_CA19_DD-grid"),
    ESA_THDD = c("2dgD_CA0_AF-grid", "2dgD_CA6_AFDD-grid", "2dgD_CA8_AF-grid", "2dgD_CA9_AF-grid", "2dgD_CA11_THDD-grid", "2dgD_CA27_DDb-grid"),
    ESA_Atlantic = c("2dgD_CA1_AF-grid", "2dgD_CA13_AF-grid", "2dgD_CA15_AF-grid", "2dgD_CA21_AF-grid", "2dgD_CA24_AF-grid")
  )
  
  if (fill_ca == "viridis") {
    
    fillingv_ca <- viridis(n = length(layers_regions_ca), alpha = alfa)
    
  }
  
  if (fill_ca == "cavalcanti") {
    
    #jpn_pal <- wes_palette("Cavalcanti1", type = "continuous", n = 7)
    jpn_pal <- c("#D8B70A", "#496715", "#376138", "#A2A475", "#8CA685", "#887F65", "#972D15")
    rep_n <- lapply(layers_regions_ca, length)
    
    fillingv_ca <- c(
      rep(jpn_pal[1], rep_n[[1]]),
      rep(jpn_pal[2], rep_n[[2]]),
      rep(jpn_pal[3], rep_n[[3]]),
      rep(jpn_pal[4], rep_n[[4]]),
      rep(jpn_pal[5], rep_n[[5]]),
      rep(jpn_pal[6], rep_n[[6]]),
      rep(jpn_pal[7], rep_n[[7]])
    )
    
    rm(rep_n, jpn_pal)
  }
  
  if (fill_ca == "darjeeling") {
    
    #jpn_pal <- wes_palette("Darjeeling1", type = "continuous", n = 7)
    jpn_pal <- c("#FF0000", "#556A5B", "#50A45C", "#F2AD00", "#6D5B4E", "#C49647", "#5BBCD6")
    rep_n <- lapply(layers_regions_ca, length)
    
    fillingv_ca <- c(
      rep(jpn_pal[1], rep_n[[1]]),
      rep(jpn_pal[2], rep_n[[2]]),
      rep(jpn_pal[3], rep_n[[3]]),
      rep(jpn_pal[4], rep_n[[4]]),
      rep(jpn_pal[5], rep_n[[5]]),
      rep(jpn_pal[6], rep_n[[6]]),
      rep(jpn_pal[7], rep_n[[7]])
    )
    
    rm(rep_n, jpn_pal)
  }
  
  pca_summ <- summary(aoe_pca)
  
  # Points in environmental space
  sites_sco <- scores(pca_env.sel, display = "sites", scaling = 2)
  sites_sco <- as_tibble(sites_sco)
  sites_sco$xname <- aoe_env$xname
  
  # Initiating list
  gg_ca <- vector(mode = "list", length = length(aoe_name))
  
  # Creating geoms list
  for (i in seq_along(aoe_name)) {
    
    hull <- sites_sco %>%
      filter(xname %in% aoe_name[i]) %>%
      group_by(xname) %>%
      slice(chull(PC1, PC2))
    
    gg_ca[[i]][[1]] <- geom_polygon(data = hull, 
                                    fill = ifelse(fill_ca == "viridis", fillingv_ca[i],
                                                  ifelse(fill_ca == "cavalcanti", fillingv_ca[i], 
                                                         ifelse(fill_ca == "darjeeling", fillingv_ca[i], fill_ca))),
                                    colour = ifelse(colour_ca == "black", "black", fillingv_ca[i]),
                                    size = size_border,
                                    alpha = alfa)
    gg_ca[[i]][[2]] <- aeo_names_vec[i]
    
  }
  
  return(gg_ca)
  
}

y <- env_geomhull(aoe_pca = pca_env.sel, aoe_env = comp_env, aoe_name = aeo_names_vec[1], fill_ca = "darjeeling", colour_ca = "darjeeling", alfa = 0.2)

aeo_names_vec <- unlist(layers_regions_ca)

y <-  env_geomhull(aoe_pca = pca_env.sel, aoe_env = comp_env, aoe_name = aeo_names_vec, fill_ca = "darjeeling", colour_ca = "darjeeling", alfa = 0.1)
str(y)

allenv_plot <- env_basemap +
  y[[10]][[1]] +
  y[[11]][[1]] +
  y[[12]][[1]] +
  y[[13]][[1]] +
  y[[14]][[1]] +
  y[[15]][[1]] +
  y[[1]][[1]] +
  y[[2]][[1]] +
  y[[3]][[1]] +
  y[[4]][[1]] +
  y[[5]][[1]] +
  y[[6]][[1]] +
  y[[7]][[1]] +
  y[[8]][[1]] +
  y[[9]][[1]] +
  y[[16]][[1]] +
  y[[17]][[1]] +
  y[[18]][[1]] +
  y[[19]][[1]] +
  y[[20]][[1]] +
  y[[21]][[1]] +
  y[[22]][[1]] +
  y[[23]][[1]] +
  y[[24]][[1]] +
  y[[25]][[1]] +
  y[[26]][[1]] +
  y[[27]][[1]] +
  y[[28]][[1]]  +
  annotation_custom(grob = grob_legend, 
                    xmin = 0, xmax = 0.55, 
                    ymin = -0.23, ymax = 0.2)

png("./figs/ALL_AoEs_Env.png", width = 1100, height = 1000, res = 130)
print(allenv_plot)
dev.off()

library(patchwork)

png("./figs/ALL_AoEs_and_Env.png", width = 1100, height = 2000, res = 130)
all_areas_plot / allenv_plot
dev.off()

## 3. Spp position in environmental space ####

# NOTE: Position was calculated in all the PCs, which is not possible to depict. 
# Therefore, we show the centroid in the PC1 and PC2. It is calculated as the mean pc score plus its standard deviation.
# This is only a way to intuitively show the possition of species in the env space.

## 3.1. Load data ####
ca_db <- read.csv("./data/biogeography/spp_listby_AoE.csv", stringsAsFactors = FALSE)[45:243,]
load("./output/2_PCA_join_spp_AoE.rdata")
comp_env <- read.csv("./output/2_Environment2.csv")

## 3.2. Defining plotting order vector ####
layers_regions_ca <- list(
  Mesoamerica = c("2dgD_CA12_NWYuc-grid"),
  Northern_Andes = c("2dgD_CA14_NWColVen-grid", "2dgD_CA17_NWCol-grid", "2dgD_CA22_NWCol2-grid", "2dgD_CA26_NWCol3-grid"),
  Guiana_centred = c("2dgD_CA3_AM-grid", "2dgD_CA10_AM-grid", "2dgD_CA25_AM-grid"),
  Amazonia_centred = c("2dgD_CA4_AM-grid", "2dgD_CA7_AM-grid", "2dgD_CA16_AM-grid", "2dgD_CA18_AM-grid", "2dgD_CA20_AM-grid", "2dgD_CA23_AM-grid"),
  ESA_DD = c("2dgD_CA2_DD-grid","2dgD_CA5_DD-grid", "2dgD_CA19_DD-grid"),
  ESA_THDD = c("2dgD_CA0_AF-grid", "2dgD_CA6_AFDD-grid", "2dgD_CA8_AF-grid", "2dgD_CA9_AF-grid", "2dgD_CA11_THDD-grid", "2dgD_CA27_DDb-grid"),
  ESA_Atlantic = c("2dgD_CA1_AF-grid", "2dgD_CA13_AF-grid", "2dgD_CA15_AF-grid", "2dgD_CA21_AF-grid", "2dgD_CA24_AF-grid")
) %>% 
  lapply(str_replace_all, pattern = ".*_(CA[0-9]+)_.*", replacement = "D2\\1")

aoe_name_vec <- unlist(layers_regions_ca)

## 3.3. Function to plot species niche position ####

# Modified from 05_Vis_niche_properties.R
sppaoe_position_pca <- function (.pca_env.sel, .comp_env, .ca_db, area_name, save_plot = "NO") {
  
  sites_sco <- scores(.pca_env.sel, display = "sites", scaling = 2)
  sites_sco <- as_tibble(sites_sco)
  sites_sco$xname <- .comp_env$xname
  
  # Breadth and position in pca mean/sd dataframe
  
  table_mean <- data.frame(pc1_mean = vector(mode = "numeric", length = length(unique(sites_sco$xname))),
                           pc1_sd = vector(mode = "numeric", length = length(unique(sites_sco$xname))),
                           pc1_lower = vector(mode = "numeric", length = length(unique(sites_sco$xname))),
                           pc1_upper = vector(mode = "numeric", length = length(unique(sites_sco$xname))),
                           pc2_mean = vector(mode = "numeric", length = length(unique(sites_sco$xname))),
                           pc2_sd = vector(mode = "numeric", length = length(unique(sites_sco$xname))),
                           pc2_lower = vector(mode = "numeric", length = length(unique(sites_sco$xname))),
                           pc2_upper = vector(mode = "numeric", length = length(unique(sites_sco$xname))))
  
  for (i in seq_along(unique(sites_sco$xname))) {
    
    j <- sites_sco %>% filter(xname == unique(sites_sco$xname)[i]) %>% pull(PC1) %>% length()
    
    if (j >= 2) {
      table_mean[i, ] <- sites_sco %>%
        filter(xname == unique(sites_sco$xname)[i]) %>%
        summarise(pc1_mean = mean(PC1),
                  pc1_sd = sd(PC1),
                  pc1_lower = pc1_mean - pc1_sd,
                  pc1_upper = pc1_mean + pc1_sd,
                  pc2_mean = mean(PC2),
                  pc2_sd = sd(PC2),
                  pc2_lower = pc2_mean - pc2_sd,
                  pc2_upper = pc2_mean + pc2_sd)
      
    } else {
      
      table_mean[i, ] <- sites_sco %>%
        filter(xname == unique(sites_sco$xname)[i]) %>%
        summarise(pc1_mean = mean(PC1),
                  pc1_sd = pc1_mean,
                  pc1_lower = pc1_mean,
                  pc1_upper = pc1_mean,
                  pc2_mean = mean(PC2),
                  pc2_sd = pc2_mean,
                  pc2_lower = pc2_mean,
                  pc2_upper = pc2_mean)
      
    }
    
  }
  
  table_med50 <- data.frame(Q0 = vector(mode = "numeric", length = length(unique(sites_sco$xname))),
                            Q1 = vector(mode = "numeric", length = length(unique(sites_sco$xname))),
                            Q2 = vector(mode = "numeric", length = length(unique(sites_sco$xname))),
                            Q3 = vector(mode = "numeric", length = length(unique(sites_sco$xname))),
                            Q4 = vector(mode = "numeric", length = length(unique(sites_sco$xname))))
  
  for (i in seq_along(unique(sites_sco$xname))) {
    
    table_med50[i, ] <- sites_sco %>%
      filter(xname == unique(sites_sco$xname)[i]) %>%
      pull(PC1) %>%
      quantile()
    
  }
  
  table_mean$xname <- unique(sites_sco$xname)
  table_med50$xname <- unique(sites_sco$xname)
  stat_df <- left_join(table_mean, table_med50, by = "xname")
  stat_df <- stat_df %>% arrange(pc1_mean, pc2_mean)
  
  
  imp_pc1 <- summary(.pca_env.sel)$cont$importance[2, 1] * 100
  imp_pc2 <- summary(.pca_env.sel)$cont$importance[2, 2] * 100
  
  
  ### Plot 
  
  stat_df$xname <- factor(stat_df$xname, levels = unique(stat_df$xname), ordered = TRUE)
  
  aoe_spp_vec <- .ca_db %>%
    filter(Area_name == area_name) %>%
    dplyr::select(Species) %>% 
    pull()
  
  crux_df <- stat_df %>%
    filter(xname %in% aoe_spp_vec) %>%
    mutate(PC1 = pc1_mean,
           PC2 = pc2_mean)
  
  hull <- sites_sco %>%
    filter(xname %in% area_name) %>%
    group_by(xname) %>%
    slice(chull(PC1, PC2))
  
  plot_sppaoe <- sites_sco %>%
    ggplot(aes(x = PC1, y = PC2)) +
    #geom_point(size = 0.2, alpha = 0.5, color = "gray") +
    geom_segment(data = crux_df, aes(x = pc1_lower, xend = pc1_upper, yend = PC2, color = factor(xname))) +
    geom_segment(data = crux_df, aes(y = pc2_lower, yend = pc2_upper, xend = PC1, color = factor(xname))) +
    geom_point(data = crux_df, aes(x = PC1, y = PC2, color = factor(xname)), shape = 21, fill = "black", size = 1) +
    geom_polygon(data = hull, fill = NA, colour = "black", alpha = 0.3, size = 0.2) +
    # geom_errorbarh(data = crux_df, aes(xmin = pc1_lower, xmax = pc1_upper, color = factor(xname))) +
    # geom_errorbar(data = crux_df, aes(ymin = pc2_lower, ymax = pc2_upper, color = factor(xname))) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_hline(aes(yintercept = 0), linetype = "dashed") +
    coord_cartesian(xlim = c(-0.3748672, 0.4475069), ylim = c(-0.4032099, 0.1944293)) +
    scale_color_grey() +
    theme_classic() +
    theme(legend.position = "none") +
    labs(x = paste("PC1", format(imp_pc1, digits = 4), "%", sep = " "), 
         y = paste("PC2", format(imp_pc2, digits = 4), "%", sep = " "),
         subtitle = "Niche position")
  
  if (save_plot == "YES") {
    
    png(paste0("./figs/AoE/2_pca_sppaoe_position_", area_name, ".png"), width = 800, height = 800, res = 120)
    print(plot_sppaoe)
    dev.off()
    
  } else {
    
    return(plot_sppaoe)
    
  }
  
}

jaja <- sppaoe_position_pca(.pca_env.sel = pca_env.sel, .comp_env = comp_env, .ca_db = ca_db, area_name = "D2CA0", save_plot = "NO")

pos_list <- purrr::map(.x = aoe_name_vec, sppaoe_position_pca, .pca_env.sel = pca_env.sel, .comp_env = comp_env, .ca_db = ca_db, save_plot = "NO")

#saveRDS(pos_list, file = "./figs/Panel_figs/C_SppPos.rds")

pos_list <- readRDS("./figs/Panel_figs/C_SppPos.rds")
pos_list[1]

i <- 2
(maps_list[[i]] | env_list[[i]]) / (pos_list[[i]] | plot_spacer()) + plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face = "bold"))

## 4. Spp breadth in environmental space ####

library(tidyverse)
library(scales)

## 4.1. Loading data ####

ca_db <- read.csv("./data/biogeography/spp_listby_AoE.csv", stringsAsFactors = FALSE)[45:243,]

niche_properties_allspp <- read.csv("./output/Breadth/2_Niche_properties_Allspp_PCA.csv")[-1] %>%
  mutate(pos = Mean_dist_cent + Sd_dist_cent,
         pos2 = Mean_dist_cent - Sd_dist_cent,
         sd_range = pos - pos2)

# Organizing plot by the range between max and min sd to define niche class
# a <- niche_properties_allspp %>%
#   ggplot(aes(x = 1, y = sd_range)) +
#   geom_boxplot()
# 
# a <- ggplot_build(a)
# a
# # lower = 0.07392203
# # middle = 0.1029429
# # upper = 0.1324147
# # min = 0
# # max = 1.375

# Organizing plots by Mean_dist to define niche class
# a <- niche_properties_allspp %>%
#   ggplot(aes(x = 1, y = Mean_dist_cent)) +
#   geom_boxplot()
# 
# a <- ggplot_build(a)
# a
# # lower = 0.1006569
# # middle = 0.1278217
# # upper = 0.1577236
# # min = 0.01772329
# # max = 0.2414453

niche_properties_allspp <- niche_properties_allspp %>%
  mutate(niche_category = case_when(
    sd_range <= 0.07392203 ~ "Narrow", # lower = 0.07392203
    sd_range > 0.07392203 &  sd_range <= 0.1324147  ~ "Medium", # upper = 0.1324147
    sd_range > 0.1324147 ~ "Wide" # upper = 0.1324147
  ), 
  niche_category2 = case_when(
    Mean_dist_cent <= 0.1006569 ~ "Narrow", # lower = 0.1006569
    Mean_dist_cent > 0.1006569 &  Mean_dist_cent <= 0.1577236  ~ "Medium", # upper = 0.1324147
    Mean_dist_cent > 0.1577236 ~ "Wide" # upper = 0.1324147
  )) %>%
  dplyr::select(-c("pos", "pos2", "sd_range"))

glimpse(niche_properties_allspp)

## 4.1.1. AoE breadth ####

# NOTE: The idea of including the breadth of the AoE was abandoned. 
# Given the higher number of points (one per pixel) the dispersion is low.
# It is even lower than the dispersion for some endemic species of the same AoE.

# aoe_breadth_df <- read.csv("./output/AoE/Breadth_AoE_SPPJOINT/2_Breadth_AoEspp_joint_aoeOnly.csv")[-1] 
# glimpse(aoe_breadth_df)
# 
# # a <- aoe_breadth_df %>%
# #   ggplot(aes(y = Mean_dist_cent)) +
# #   geom_boxplot()
# # a <- ggplot_build(a)
# # # lower = 0.08935668
# # # middle = 0.1472739
# # # upper = 0.1657443
# # # min = 0.08935668
# # # max = 0.2204663
# 
# aoe_breadth_df %>%
#   mutate(niche_category = case_when(
#     Mean_dist_cent <= 0.08935668 ~ "Narrow", # lower = 0.1006569
#     Mean_dist_cent > 0.08935668 &  Mean_dist_cent <= 0.1657443  ~ "Medium", # upper = 0.1324147
#     Mean_dist_cent > 0.1657443 ~ "Wide" # upper = 0.1324147
#   )) %>%
#   ggplot(aes(x = Item_ID, y = Mean_dist_cent)) +
#   geom_point(size = 1.5, aes(colour = niche_category)) +
#   geom_errorbar(aes(ymin = Mean_dist_cent - Sd_dist_cent, ymax = Mean_dist_cent + Sd_dist_cent, colour = niche_category), width = 0.2) +
#   ylim(0, 0.4) +
#   coord_flip() +
#   scale_color_manual(values = c("Narrow" = "#333333", "Medium" = "#989898", "Wide" = "#CCCCCC")) +
#   #scale_color_grey() +
#   theme_classic() +
#   theme(axis.title.y = element_blank(),
#         axis.text.y = element_text(face = "italic"), 
#         legend.position = "right") +
#   labs(subtitle = paste("Niche breadth of endemic species in", "aoe_name", sep = " "),
#        y = "Mean distance to centroid in PCA",
#        colour = "Niche class")


## 4.2. Defining plotting order vector ####
layers_regions_ca <- list(
  Mesoamerica = c("2dgD_CA12_NWYuc-grid"),
  Northern_Andes = c("2dgD_CA14_NWColVen-grid", "2dgD_CA17_NWCol-grid", "2dgD_CA22_NWCol2-grid", "2dgD_CA26_NWCol3-grid"),
  Guiana_centred = c("2dgD_CA3_AM-grid", "2dgD_CA10_AM-grid", "2dgD_CA25_AM-grid"),
  Amazonia_centred = c("2dgD_CA4_AM-grid", "2dgD_CA7_AM-grid", "2dgD_CA16_AM-grid", "2dgD_CA18_AM-grid", "2dgD_CA20_AM-grid", "2dgD_CA23_AM-grid"),
  ESA_DD = c("2dgD_CA2_DD-grid","2dgD_CA5_DD-grid", "2dgD_CA19_DD-grid"),
  ESA_THDD = c("2dgD_CA0_AF-grid", "2dgD_CA6_AFDD-grid", "2dgD_CA8_AF-grid", "2dgD_CA9_AF-grid", "2dgD_CA11_THDD-grid", "2dgD_CA27_DDb-grid"),
  ESA_Atlantic = c("2dgD_CA1_AF-grid", "2dgD_CA13_AF-grid", "2dgD_CA15_AF-grid", "2dgD_CA21_AF-grid", "2dgD_CA24_AF-grid")
) %>% 
  lapply(str_replace_all, pattern = ".*_(CA[0-9]+)_.*", replacement = "D2\\1")

aoe_name_vec <- unlist(layers_regions_ca)

## 4.3. Function to plot breadth ####

# Modified from 05_Vis_Niche_properties.R
niche_breadth_plot <- function (breadth_df, aoe_df, aoe_name, save_plot = "NO") {
  
  library(ggplot2)
  library(wesanderson)
  
  endspp <- aoe_df %>%
    filter(Area_name == aoe_name) %>%
    dplyr::select(Species) %>%
    pull()
  
  interm_df <- breadth_df %>%
    filter(NAME1 %in% endspp) %>%
    mutate(pos = Mean_dist_cent + Sd_dist_cent,
           pos2 = Mean_dist_cent - Sd_dist_cent,
           sd_range = pos - pos2) %>%
    arrange(Mean_dist_cent) %>%
    filter(Mean_dist_cent != 0)
  
  interm_df$NAME1 <- factor(interm_df$NAME1, levels = interm_df$NAME1, ordered = TRUE)
  interm_df$niche_category2 <- factor(interm_df$niche_category2, levels = c("Wide", "Medium", "Narrow"), ordered = TRUE)
  
  breadth_plot <- interm_df %>% 
    ggplot(aes(x = NAME1, y = Mean_dist_cent)) +
    geom_point(size = 1, aes(colour = niche_category2)) +
    #geom_point(size = 1) +
    geom_errorbar(aes(ymin = Mean_dist_cent - Sd_dist_cent, ymax = Mean_dist_cent + Sd_dist_cent, colour = niche_category2), width = 0.2) +
    #geom_errorbar(aes(ymin = Mean_dist_cent - Sd_dist_cent, ymax = Mean_dist_cent + Sd_dist_cent), width = 0.2) +
    ylim(0, 0.3) +
    coord_flip() +
    geom_hline(yintercept = 0.1006569, lty = 2) +
    geom_hline(yintercept = 0.1577236, lty = 2) +
    scale_color_manual(values = c("Narrow" = "#333333", "Medium" = "#989898", "Wide" = "#CCCCCC")) +
    #scale_color_grey() +
    annotate(geom = "text", x = 0.5, y = 0.097, label = "Narrow", size = 2, hjust = 1) +
    annotate(geom = "text", x = 0.5, y = 0.13, label = "Medium", size = 2, hjust = 0.5) +
    annotate(geom = "text", x = 0.5, y = 0.16, label = "Wide", size = 2, hjust = 0) +
    theme_classic() +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_text(face = "italic"), 
          legend.position = "none") +
    labs(subtitle = "Niche breadth",
         y = "Mean distance to centroid in PCA",
         colour = "Niche class")
  
  if (aoe_name == "D2CA0") {
    
    breadth_plot <- breadth_plot +
      theme(axis.text.y = element_text(size = 1))
    
  }
  
  if (save_plot == "YES") {
    
    png(paste0("./figs/AoE/2_breadth_", aoe_name, ".png"), width = 800, height = 800, res = 120)
    print(breadth_plot)
    dev.off()
    
  } else {
    
    return(breadth_plot)
    
  }
}

niche_breadth_plot(breadth_df = niche_properties_allspp, aoe_df = ca_db, aoe_name = "D2CA0")

bre_list <- purrr::map(.x = aoe_name_vec, niche_breadth_plot,breadth_df = niche_properties_allspp, aoe_df = ca_db, save_plot = "NO")
bre_list[[28]]

saveRDS(bre_list, file = "./figs/Panel_figs/D_SppBre.rds")

bre_list <- readRDS("./figs/paper/Panel_regionalenv/D_SppBre.rds")

i <- 28
(maps_list[[i]] | env_list[[i]]) / (pos_list[[i]] | bre_list[[i]]) + plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face = "bold"))

## 5. Niche overlap ####

library(tidyverse)

# 5.1. Loading data ####

locs <-  read.csv("./output/localities_total_unique_ch10m_Bio(v14).csv", row.names = 1, stringsAsFactors = FALSE) # Localidades_v14
ca_db <- read.csv("./data/biogeography/spp_listby_AoE.csv", stringsAsFactors = FALSE)[45:243,]

locs_end <- locs %>%
  filter(NAME1 %in% unique(ca_db$Species))
spp_names <- unique(locs_end$NAME1)

AoE_comp_df <- read.csv("./output/AoE/AoE_ECOSPAT_RES_areasonly.csv") %>%  filter(Analysis == "Eq", Alternative == "G")

AoE_comp_df$AoE_regs <- factor(AoE_comp_df$AoE_regs, levels = c("Mesoamerica", "Northern_Andes", "Guiana_centered", "Amazonia_centered", 
                                                                "ESA_Dry_Diagonal", "ESA_THDD", "ESA_Atlantic_Forest"), ordered = TRUE)
AoE_comp_df$D_classes <- factor(AoE_comp_df$D_classes, levels = c("Limited", "Low", "Moderate", "High", "Very_high"), ordered = TRUE)
AoE_comp_df$I_classes <- factor(AoE_comp_df$I_classes, levels = c("Limited", "Low", "Moderate", "High", "Very_high"), ordered = TRUE)
AoE_comp_df$Geo_classes <- factor(AoE_comp_df$Geo_classes, levels = c("Disjunct", "Low", "Medium", "High", "Very_high"), ordered = TRUE)
AoE_comp_df <- AoE_comp_df %>% arrange(AoE_regs)

Spname1 <- vector(mode = "character", length = length(AoE_comp_df$Sp.1))
Spname2 <- vector(mode = "character", length = length(AoE_comp_df$Sp.2))
for (i in seq_along(AoE_comp_df$Sp.1)) {
  
  Spname1[i] <- spp_names[AoE_comp_df$Sp.1[i]]
  Spname2[i] <- spp_names[AoE_comp_df$Sp.2[i]]
  
}
AoE_comp_df$Sp1name <- Spname1
AoE_comp_df$Sp2name <- Spname2

AoE_comp_df <- AoE_comp_df %>%
  mutate(Sp1namec = str_replace_all(string = Sp1name, pattern = "^([A-z]{5}).*\\s([A-z]{3}).*$", replacement = paste("\\1", "\\2", sep = "_")),
         Sp2namec = str_replace_all(string = Sp2name, pattern = "^([A-z]{5}).*\\s([A-z]{3}).*$", replacement = paste("\\1", "\\2", sep = "_")))

jpn_pal <- wes_palette("Zissou1", n = 6, type = "continuous")

## 5.2. Defining plotting order vector ####
layers_regions_ca <- list(
  Mesoamerica = c("2dgD_CA12_NWYuc-grid"),
  Northern_Andes = c("2dgD_CA14_NWColVen-grid", "2dgD_CA17_NWCol-grid", "2dgD_CA22_NWCol2-grid", "2dgD_CA26_NWCol3-grid"),
  Guiana_centred = c("2dgD_CA3_AM-grid", "2dgD_CA10_AM-grid", "2dgD_CA25_AM-grid"),
  Amazonia_centred = c("2dgD_CA4_AM-grid", "2dgD_CA7_AM-grid", "2dgD_CA16_AM-grid", "2dgD_CA18_AM-grid", "2dgD_CA20_AM-grid", "2dgD_CA23_AM-grid"),
  ESA_DD = c("2dgD_CA2_DD-grid","2dgD_CA5_DD-grid", "2dgD_CA19_DD-grid"),
  ESA_THDD = c("2dgD_CA0_AF-grid", "2dgD_CA6_AFDD-grid", "2dgD_CA8_AF-grid", "2dgD_CA9_AF-grid", "2dgD_CA11_THDD-grid", "2dgD_CA27_DDb-grid"),
  ESA_Atlantic = c("2dgD_CA1_AF-grid", "2dgD_CA13_AF-grid", "2dgD_CA15_AF-grid", "2dgD_CA21_AF-grid", "2dgD_CA24_AF-grid")
) %>% 
  lapply(str_replace_all, pattern = ".*_(CA[0-9]+)_.*", replacement = "D2\\1")

aoe_name_vec <- unlist(layers_regions_ca)

## 5.3. Function to plot overlap ####

over_heatmap_plot <- function(df_end_spp, aoe_name, analysis = "Eq", alternative = "G", plot_save = "NO") {
  
  library(ggplot2)
  library(dplyr)
  library(wesanderson)
  
  #jpn_pal <- wes_palette("Zissou1", n = 6, type = "continuous")
  # breaks_jpng <- seq(0, 1, by = 0.1)
  # jpn_pal <- colorRampPalette(c("#3B9AB2", "#9EBE91", "#E4B80E", "#F21A00"))(length(breaks) - 1)
  
  heatmap_plot <- df_end_spp %>%
    filter(AoE_name == aoe_name, Analysis == analysis, Alternative == alternative) %>%
    distinct(D, .keep_all = TRUE) %>%
    ggplot(aes(x = Sp1namec, y = Sp2namec, fill = D)) +
    geom_tile(alpha = 0.8) +
    #scale_fill_gradientn(colours = jpn_pal, breaks = breaks_jpng, label = breaks_jpng) +
    #scale_fill_gradientn(colors = jpn_pal, values = c(0, 0.2, 0.4, 0.6, 0.8, 1), breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)) +
    # scale_fill_gradient2(midpoint = 0.5) +
    #scale_fill_viridis_c() +
    #scale_fill_gradient(low = "purple", high = "yellow") +
    scale_fill_gradient(low = "#E6E6E6", high = "#4D4D4D") +
    theme_light()  +
    theme(axis.text = element_text(face = "italic"),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          axis.title = element_blank(),
          legend.position = "right", # c(0.85,0.25),
          legend.title = element_text(size = 10, face = "bold"),
          axis.line = element_line(colour = "black"),
          axis.ticks = element_line(colour = "black")) +
    labs(fill = "D",
         subtitle = "Niche overlap")
  
  
  if (plot_save == "YES") {
    
    png(paste0("./figs/AoE/2_Overlap_heatmap_", aoe_name, ".png"), width = 1000, height = 1000, res = 80)
    print(heatmap_plot)
    dev.off()
    
  }
  
  if (aoe_name == "D2CA0") {
    
    heatmap_plot <- heatmap_plot +
      theme(axis.text = element_text(size = 2))
    
  }
  
  return(heatmap_plot)
  
}

over_heatmap_plot(df_end_spp = AoE_comp_df, aoe_name = "D2CA23", plot_save = "NO")

over_list <-purrr::map(aoe_name_vec, over_heatmap_plot, df_end_spp = AoE_comp_df, plot_save = "NO")
over_list[[14]]

# saveRDS(over_list, file = "./figs/Panel_figs/E_SppOver.rds")

over_list <- readRDS("./figs/Panel_figs/E_SppOver.rds")

i <- 7
(maps_list[[i]] | env_list[[i]]) / (pos_list[[i]] | bre_list[[i]]) / (over_list[[i]] | plot_spacer()) + 
  plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face = "bold"))

## 6. Equivalence and similarity tests ####

## 6.1. Load data ####

locs <-  read.csv("./output/localities_total_unique_ch10m_Bio(v14).csv", row.names = 1, stringsAsFactors = FALSE) # Localidades_v14
ca_db <- read.csv("./data/biogeography/spp_listby_AoE.csv", stringsAsFactors = FALSE)[45:243,]
locs_end <- locs %>%
  filter(NAME1 %in% unique(ca_db$Species))
spp_names <- unique(locs_end$NAME1)

AoE_comp_df <- read.csv("./output/AoE/AoE_ECOSPAT_RES_areasonly.csv") #%>%  filter(Analysis == "Eq", Alternative == "G")

AoE_comp_df$AoE_regs <- factor(AoE_comp_df$AoE_regs, levels = c("Mesoamerica", "Northern_Andes", "Guiana_centered", "Amazonia_centered", 
                                                                "ESA_Dry_Diagonal", "ESA_THDD", "ESA_Atlantic_Forest"), ordered = TRUE)
AoE_comp_df$D_classes <- factor(AoE_comp_df$D_classes, levels = c("Limited", "Low", "Moderate", "High", "Very_high"), ordered = TRUE)
AoE_comp_df$I_classes <- factor(AoE_comp_df$I_classes, levels = c("Limited", "Low", "Moderate", "High", "Very_high"), ordered = TRUE)
AoE_comp_df$Geo_classes <- factor(AoE_comp_df$Geo_classes, levels = c("Disjunct", "Low", "Medium", "High", "Very_high"), ordered = TRUE)
AoE_comp_df <- AoE_comp_df %>% arrange(AoE_regs)

Spname1 <- vector(mode = "character", length = length(AoE_comp_df$Sp.1))
Spname2 <- vector(mode = "character", length = length(AoE_comp_df$Sp.2))
for (i in seq_along(AoE_comp_df$Sp.1)) {
  
  Spname1[i] <- spp_names[AoE_comp_df$Sp.1[i]]
  Spname2[i] <- spp_names[AoE_comp_df$Sp.2[i]]
  
}
AoE_comp_df$Sp1name <- Spname1
AoE_comp_df$Sp2name <- Spname2

AoE_comp_df <- AoE_comp_df %>%
  mutate(Sp1namec = str_replace_all(string = Sp1name, pattern = "^([A-z]{5}).*\\s([A-z]{3}).*$", replacement = paste("\\1", "\\2", sep = "_")),
         Sp2namec = str_replace_all(string = Sp2name, pattern = "^([A-z]{5}).*\\s([A-z]{3}).*$", replacement = paste("\\1", "\\2", sep = "_")))


## 6.2. Defining plotting order vector ####
layers_regions_ca <- list(
  Mesoamerica = c("2dgD_CA12_NWYuc-grid"),
  Northern_Andes = c("2dgD_CA14_NWColVen-grid", "2dgD_CA17_NWCol-grid", "2dgD_CA22_NWCol2-grid", "2dgD_CA26_NWCol3-grid"),
  Guiana_centred = c("2dgD_CA3_AM-grid", "2dgD_CA10_AM-grid", "2dgD_CA25_AM-grid"),
  Amazonia_centred = c("2dgD_CA4_AM-grid", "2dgD_CA7_AM-grid", "2dgD_CA16_AM-grid", "2dgD_CA18_AM-grid", "2dgD_CA20_AM-grid", "2dgD_CA23_AM-grid"),
  ESA_DD = c("2dgD_CA2_DD-grid","2dgD_CA5_DD-grid", "2dgD_CA19_DD-grid"),
  ESA_THDD = c("2dgD_CA0_AF-grid", "2dgD_CA6_AFDD-grid", "2dgD_CA8_AF-grid", "2dgD_CA9_AF-grid", "2dgD_CA11_THDD-grid", "2dgD_CA27_DDb-grid"),
  ESA_Atlantic = c("2dgD_CA1_AF-grid", "2dgD_CA13_AF-grid", "2dgD_CA15_AF-grid", "2dgD_CA21_AF-grid", "2dgD_CA24_AF-grid")
) %>% 
  lapply(str_replace_all, pattern = ".*_(CA[0-9]+)_.*", replacement = "D2\\1")

aoe_name_vec <- unlist(layers_regions_ca)


## 6.3. Function to plot equivalence and similiraty tests ####

# Modified from 05_Vis_niche_properties.R
tests_plot <- function(df_end_spp, aoe_name, save_plot = "NO") {
  
  overlap_D <- df_end_spp %>%
    filter(AoE_name == aoe_name) %>%
    mutate(Alternative = case_when(
      Alternative == "G" ~ "Greater", 
      Alternative == "L" ~ "Lower"
    ),
    Analysis = case_when(
      Analysis == "Eq" ~ "Equivalence",
      Analysis == "Sim" ~ "Similarity"
    ),
    Dp_class = case_when(
      Dp_class == "P<0.05" ~ "Rejected",
      Dp_class == "P>=0.05" ~ "Not Rejected"
    ))
  
  testplot <- overlap_D %>%
    ggplot(aes(x = Dp_class)) +
    geom_bar(position = "dodge", width = 0.5) +
    facet_grid(Alternative ~ Analysis) +
    theme_linedraw() +
    theme(strip.text = element_text(face = "bold"),
          axis.text.x = element_text(hjust = 1, vjust = 1, angle = 45)) +
    labs(x = "Result of the significance test (alpha = 0.05)",
         y = "Species pair-wise comparisons",
         subtitle = "Niche equivalence and similarity tests")
  
  if (save_plot == "YES") {
    
    png(paste0("./figs/AoE/2_aoe_ecospat_tests_", AoE_name, ".png"), width = 700, height = 700, res = 150)
    print(testplot)
    dev.off()
    
  } else {
    
    return(testplot)
    
  }
  
}

tests_plot(df_end_spp = AoE_comp_df, aoe_name = "D2CA23", save_plot = "NO")


tests_list <-purrr::map(aoe_name_vec, tests_plot, df_end_spp = AoE_comp_df, save_plot = "NO")
tests_list[14]

saveRDS(tests_list, file = "./figs/Panel_figs/F_EqSimtests.rds")

tests_list <- readRDS("./figs/Panel_figs/F_EqSimtests.rds")

##  7. Composing the panel ####
library(patchwork)

maps_list <- readRDS(file = "./figs/Panel_figs/A_AoEmaps.rds")
env_list  <- readRDS("./figs/Panel_figs/B_AoEenv.rds")
pos_list <- readRDS("./figs/Panel_figs/C_SppPos.rds")
bre_list <- readRDS("./figs/Panel_figs/D_SppBre.rds")
over_list <- readRDS("./figs/Panel_figs/E_SppOver.rds")
tests_list <- readRDS("./figs/Panel_figs/F_EqSimtests.rds")

x <- seq(1, 28, 1)#[-14]

for (i in seq_along(x)) {
  
  png(paste("./figs/Panel_figs/", x[i],"_28.png", sep ="_"), width = 1400, height = 1800, res = 190)
  print((maps_list[[x[i]]] | env_list[[x[i]]]) / (pos_list[[x[i]]] | bre_list[[x[i]]]) / (over_list[[x[i]]] | tests_list[[x[i]]]) + 
          plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face = "bold")))
  dev.off()
  
}

# # Panel for D2CA23
# i <- 14
# png("./figs/paper/Panel_regionalenv/Panel1__14.png", width = 1400, height = 1800, res = 190)
# (maps_list[[i]] | env_list[[i]]) / (pos_list[[i]] | bre_list[[i]]) / (plot_spacer() | plot_spacer()) + 
#   plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face = "bold"))
# dev.off()

# end_spp <- ca_db %>%
#   filter(Area_name == "D2CA23") %>%
#   pull()
# 
# overlap_D <- all_spp_ecospat %>%
#   filter(Sp1name %in% end_spp, Sp2name %in% end_spp) %>%
#   mutate(Alternative = case_when(
#     Alternative == "G" ~ "Greater", 
#     Alternative == "L" ~ "Lower"
#   ),
#   Analysis = case_when(
#     Analysis == "Eq" ~ "Equivalence",
#     Analysis == "Sim" ~ "Similarity"
#   ),
#   Dp_class = case_when(
#     Dp_class == "P<0.05" ~ "Rejected",
#     Dp_class == "P>=0.05" ~ "Not Rejected"
#   ))


