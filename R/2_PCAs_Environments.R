## Building environmental spaces.
## Two pca objects are produced: one only for the endemic species (1_) [Used for overlap]
## and the other for all the species plus AoE's environment (2_). [Used for breadt and position]
## The environmental space is defined by the uncorrelated bioclim variables obtained from species occurrence data.
## Variables correlational structure is examined over all the species including non-endemics.

# Libraries

library(caret)
library(vegan)
library(dplyr)
library(stringr)
library(ggplot2)
library(ggordiplots)
library(viridis)
library(tidyverse)

### 0. Loading data ####
big_db <- read.csv("./data/occs/Localidades_Total_v14.csv")
ca_db <- read.csv("./data/biogeography/spp_listby_AoE.csv", stringsAsFactors = FALSE)[45:243,]
AoE_data <- read.csv("./data/biogeography/AoE_raster/AoE_Chelsa10m_bio.csv", stringsAsFactors = FALSE) # AoE Environment
locs <- read.csv("./output/localities_total_unique_ch10m_Bio(v14).csv", row.names = 1, stringsAsFactors = FALSE) # Localidades_v14
ambiente <- locs %>% dplyr::select(6:24) # Spp bioclim (Dangerous variables: Bio8,9,18,19)

### 1. Exclusion of correlated variables in species environment (Avoiding redundancy) ####
corr <- cor(scale(ambiente))
spp_env_sel <- data.frame(ambiente[, -findCorrelation(corr, cutoff = 0.7)])
names(spp_env_sel)
#[1] "bio_10" "bio_12" "bio_14" "bio_2" "bio_5" "bio_7" (Species)
# [1] "bio_10" "bio_11" "bio_12" "bio_13" "bio_2"  "bio_4" (Areas of Endemism)

# Pearson Correlation Table

cor_env_df <- cor(spp_env_sel)
corrplot::corrplot(cor_env_df, type = "upper", diag = FALSE, addCoef.col = TRUE, number.cex = 1)
#write.csv(cor_env_df, "./output/SPP/2_Bioclim_corr07.csv")

# Using species uncorrelated variables to define AoE environmenatal space: #[1] "bio_10" "bio_12" "bio_14" "bio_2" "bio_5" "bio_7"
aoe_env_sel <- dplyr::select(AoE_data, c("bio_10", "bio_12", "bio_14", "bio_2", "bio_5", "bio_7"))
glimpse(aoe_env_sel)
glimpse(spp_env_sel)

### 2. Joining environmental dataframes ####

comp_env <- rbind(spp_env_sel, aoe_env_sel)
names_ssp_aoe <- c(locs$NAME1, AoE_data$AoE_name)
comp_env$xname <- names_ssp_aoe
comp_env$class_id <- c(rep("spp", length(locs$NAME1)), rep("aoe", length(AoE_data$AoE_name)))
glimpse(comp_env)

# Saving environment
#write.csv(comp_env, "./output/2_Environment2.csv")

### 3. Join species-AoE PCA for beadth and position ####

pca_env.sel <- rda(scale(comp_env[-c(7, 8)]))
plot(pca_env.sel, scaling = -1)
#save(pca_env.sel,file = "./output/2_PCA_join_spp_AoE.rdata") 
load("./output/2_PCA_join_spp_AoE.rdata")

pca_summ <- summary(pca_env.sel)

# Table of scores
#write.csv(pca_summ$species, "./output/2_PCA_spp_scores.csv")
### 4. Figure of env space ####

PCA_plot <- function (area_name, areas_df, species_df, areas_env_df, save.pcaplot = "NO") {
  
  # Libraries
  library(caret)
  library(vegan)
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(ggordiplots)
  library(viridis)
  library(ecospat)
  library(ggrepel)
  library(gtable)
  library(grid)
  library(gridExtra)
  library(patchwork)
  
  if (area_name != "Regional") {
    
    # Prepearing data
    end_vec <- areas_df %>%
      filter(Area_name == area_name) %>%
      pull(Species)
    
    spp_df <- species_df %>%
      filter(NAME1 %in% end_vec) %>%
      dplyr::select(-c(XCOOR, YCOOR))
    
    aoe_df <- areas_env_df %>% 
      filter(AoE_name == area_name)
    
    names(spp_df) <- c("xname", "ID", "cells", "bio_1", "bio_10", "bio_11", "bio_12", "bio_13", "bio_14", "bio_15", "bio_16", "bio_17",
                       "bio_18", "bio_19", "bio_2", "bio_3", "bio_4", "bio_5", "bio_6", "bio_7", "bio_8", "bio_9")
    
    names(aoe_df) <- c("xname", "ID", "cells", "bio_1", "bio_10", "bio_11", "bio_12", "bio_13", "bio_14", "bio_15", "bio_16", "bio_17",
                       "bio_18", "bio_19", "bio_2", "bio_3", "bio_4", "bio_5", "bio_6", "bio_7", "bio_8", "bio_9")
    
    comp_env <- rbind(aoe_df, spp_df) %>%
      dplyr::select(-c(bio_8, bio_9, bio_18, bio_19))
    
    ambiente <- comp_env %>% dplyr::select(-c(xname, ID, cells)) # Spp bioclim (Dangerous variables: Bio8,9,18,19)
    
    ### 1. Exclusion of correlated variables in species environment (Avoiding redundancy) ####
    corr <- cor(scale(ambiente))
    env_sel <- data.frame(ambiente[, -findCorrelation(corr, cutoff = 0.7)])
    
  } else {
    
    ambiente <- species_df %>% dplyr::select(6:24) # Spp bioclim (Dangerous variables: Bio8,9,18,19)
    corr <- cor(scale(ambiente))
    spp_env_sel <- data.frame(ambiente[, -findCorrelation(corr, cutoff = 0.7)])
    aoe_env_sel <- dplyr::select(areas_env_df, c("bio_10", "bio_12", "bio_14", "bio_2", "bio_5", "bio_7"))
    env_sel <- rbind(spp_env_sel, aoe_env_sel)
    # names_ssp_aoe <- c(species_df$NAME1, areas_env_df$AoE_name)
    # env_sel$xname <- names_ssp_aoe
    
  }
  
  ### 2. Principal component analysis ####
  
  pca_env.sel <- rda(scale(env_sel))
  pca_summ <- summary(pca_env.sel)
  
  ## Saving pca data, pca scores, and importance
  
  pcscores_table <- format(pca_env.sel$CA$v, digits = 2)
  pcimportance_table <- format(as.data.frame(pca_summ$cont),digits = 3)
  colnames(pcimportance_table) <- str_replace_all(string = colnames(pcimportance_table), pattern = "importance\\.(PC[0-9]+)", replacement = "\\1") 
  
  ### 3. Env space and var plots ####
  
  spp_sco1 <- scores(pca_env.sel, display = "species")
  spp_sco1 <- as_tibble(spp_sco1)
  sites_sco1 <- scores(pca_env.sel, display = "sites", scaling = 2)
  sites_sco1 <- as_tibble(sites_sco1)
  
  # Note that tibbles do not contain rownames so we need to add them
  spp_tbl <- mutate(spp_sco1, vgntxt = rownames(scores(pca_env.sel, display = "species")))
  
  plt_vars <- ggplot(data = spp_tbl, aes(x = PC1, y = PC2, label = vgntxt)) +
    geom_text_repel(seed = 123, color = "black") +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_hline(aes(yintercept = 0), linetype = "dashed") +
    geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2), arrow = arrow(length = unit(0.2,"cm")),
                 color="blue") +
    theme_classic()
  
  plot_env_space <- ggplot(data = sites_sco1, aes(x = PC1, y = PC2)) +
    geom_point(size = 0.1, color = "gray") +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_hline(aes(yintercept = 0), linetype = "dashed") +
    theme_classic() +
    labs(x = paste("PC1", paste0(format(as.data.frame(pca_summ$cont)[2, 1]*100, digits = 4), "%"), sep = " "),
         y = paste("PC2", paste0(format(as.data.frame(pca_summ$cont)[2, 2]*100, digits = 4), "%"), sep = " "))
  
  
  importance_g <- tableGrob(pcimportance_table, theme = ttheme_minimal())
  
  importance_g <- gtable_add_grob(importance_g,
                                  grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                                  t = 1, l = 1, r = ncol(importance_g))
  
  importance_g <- gtable_add_grob(importance_g,
                                  grobs = segmentsGrob( # line across the bottom
                                    x0 = unit(0,"npc"),
                                    y0 = unit(0,"npc"),
                                    x1 = unit(1,"npc"),
                                    y1 = unit(0,"npc"),
                                    gp = gpar(lwd = 2.0)),
                                  t = nrow(importance_g),  nrow(importance_g), l = 1, r =  ncol(importance_g))
  
  scores_g <- tableGrob(pcscores_table, theme = ttheme_minimal())
  
  scores_g <- gtable_add_grob(scores_g,
                              grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                              t = 1, l = 1, r = ncol(scores_g))
  
  scores_g <- gtable_add_grob(scores_g,
                              grobs = segmentsGrob( # line across the bottom
                                x0 = unit(0,"npc"),
                                y0 = unit(0,"npc"),
                                x1 = unit(1,"npc"),
                                y1 = unit(0,"npc"),
                                gp = gpar(lwd = 2.0)),
                              t = nrow(scores_g),  nrow(scores_g), l = 1, r =  ncol(scores_g))
  
  degs <- "
  1122
  3333
  4444
  "
  
  # degs <- "
  # 11
  # 11
  # 22
  # 33
  # "
  
  p <- plt_vars + plot_env_space +
    importance_g + scores_g +
    plot_layout(design = degs) +
    plot_annotation(title = area_name,
                    tag_levels = "A") &
    theme(plot.tag = element_text(face = "bold"),
          plot.title = element_text(size = 15, face = "bold"))
  
  # p <- plot_env_space + 
  #   importance_g + scores_g + 
  #   plot_layout(design = degs) + 
  #   plot_annotation(#title = area_name,
  #                   tag_levels = "A") &
  #   theme(plot.tag = element_text(face = "bold"),
  #         plot.title = element_text(size = 15, face = "bold"))
  
  if (save.pcaplot == "YES") {
    
    png(paste0("./figs/", area_name, ".png"), width = 850, height = 1000, res = 120)
    print(p)
    dev.off()
    
  } else {
    
    return(p)
    
  }
  
}

x <- PCA_plot(area_name = "D2CA0", areas_df = ca_db, species_df = locs, areas_env_df = AoE_data, save.pcaplot = "NO")
x

PCA_plot(area_name = "Regional", areas_df = ca_db, species_df = locs, areas_env_df = AoE_data, save.pcaplot = "YES")


### 5. PCA for overlap ####

# NOTE: This PCA was built using only the species occurrence points to ease computational load.

big_db <- read.csv("./data/occs/Localidades_Total_v14.csv")
locs <- read.csv("./output/localities_total_unique_ch10m_Bio(v14).csv", row.names = 1) # Localidades_v14
ambiente <- locs %>% dplyr::select(6:24) # (Dangerous variables: Bio8,9,18,19)

corr <- cor(scale(ambiente))
ambiente.sel <- data.frame(ambiente[, -findCorrelation(corr, cutoff = 0.7)])
names(ambiente.sel)
#[1] "bio_10" "bio_12" "bio_14" "bio_2" "bio_5" "bio_7"
# [1] "bio_10" "bio_11" "bio_12" "bio_13" "bio_2"  "bio_4" (Areas of Endemism)
#write.csv(ambiente.sel, "./output/1_species_environment_uncorr.csv")

pca_env.sel <- rda(scale(ambiente.sel))
#save(pca_env.sel, file = "./output/1_pca_amb_sel(v14).rdata") #para overlap
load("./output/v14/1_pca_amb_sel(v14).rdata")

