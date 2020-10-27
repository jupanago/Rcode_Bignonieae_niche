## Script for statistical analysis of: Niche breadth with Environment 2 (Spp + AoEs Env)
## Main objective: look for patterns of niche properties similarities and differences among
## endemic species inside each Areas of Endemism.
## This file is done just to unify the code for the statistical analysis on niche breadth.

# Libraries
library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(purrr)
library(onewaytests)
library(gridExtra)
library(grid)

############## BREADTH FROM PCA ##########################
### 0. Loading data ####

## NOTE: Data frames origin:
# Data come from 03_Niche_properties_table 1 and 2 for environment 2 (Old Analyses): 
# Niche breadth from betadisper are conserved under the columns mean_distance & sd_distance.
# Creating new niche properties table
# stat_df <- read.csv("./data/OLD_ANALYSES/Allspp_Table_nicheproperties.csv")[c(10:11, 17:19,23:26)]
# Adding Breadth data from AZURE results using usedist package:
# breadth_spp_df <- read.csv("./output/Breadth/2_Breadth_AoEspp_joint_AZURE.csv")[-1]
# breadth_spp_df <- left_join(stat_df, breadth_spp_df, by = c("NAME1" = "Item_ID"))
# #write.csv(breadth_spp_df, "./output/Breadth/2_Niche_properties_Allspp_PCA.csv")

# Loading niche properties table
niche_properties_allspp <- read.csv("./output/Breadth/2_Niche_properties_Allspp_PCA.csv")[-1]

ca_db <- read.csv("./data/biogeography/spp_listby_AoE.csv", stringsAsFactors = FALSE)[45:243,] %>%
  dplyr::select(-1)

# Geographical regions
AoE_regs <- read.csv("./data/biogeography/AoE_listby_Reg.csv", stringsAsFactors = FALSE)[-1] 
ca_db <- left_join(ca_db, AoE_regs, by = c("Area_name" = "AoE_name"))

# Preparing AoE data frame
aoe_extrac_df <- function(x) {
  
  spp_aoe_vec <- ca_db %>%
    filter(Area_name == x) %>%
    dplyr::select(Species) %>%
    pull()
  
  regs_vec <- ca_db %>%
    filter(Area_name == x) %>%
    dplyr::select(AoE_regs) %>%
    pull()
  
  aoe_df <- niche_properties_allspp %>%
    filter(NAME1 %in% spp_aoe_vec) %>%
    mutate(AoE_name = x,
           AoE_regs = regs_vec)
  
  return(aoe_df)
  
}

aoe_extrac_df(x = "D2CA1")

aoe_names <- ca_db %>% 
  dplyr::select(Area_name) %>%
  pull() %>%
  unique()

aoe_dfs <- purrr::map(.x = aoe_names, .f = aoe_extrac_df) 
aoe_dfs <- do.call("rbind", aoe_dfs)

aoe_dfs$AoE_regs <- factor(aoe_dfs$AoE_regs, levels = c("Mesoamerica", "Northern_Andes", "Guiana_centered", "Amazonia_centered", 
                                                        "ESA_Dry_Diagonal", "ESA_THDD", "ESA_Atlantic_Forest"), ordered = TRUE)

### 1. Statistical summary ####
glimpse(aoe_dfs)

# By Aoe
aoe_dfs %>%
  group_by(AoE_name) %>%
  summarise(
    count = n(),
    betmean_beadth = mean(mean_distance, na.rm = TRUE),
    betsd_breadth = sd(mean_distance, na.rm = TRUE),
    pcamean_breadth = mean(Mean_dist_cent, na.rm = TRUE),
    pcasd_breadth = sd(Mean_dist_cent, na.rm = TRUE),
    median_breadth = median(Med_dist_cent , na.rm = TRUE),
    IQR_breadth = IQR(Med_dist_cent, na.rm = TRUE)
  )

# By region
aoe_dfs %>%
  group_by(AoE_regs) %>%
  summarise(
    count = n(),
    betmean_beadth = mean(mean_distance, na.rm = TRUE),
    betsd_breadth = sd(mean_distance, na.rm = TRUE),
    pcamean_breadth = mean(Mean_dist_cent, na.rm = TRUE),
    pcasd_breadth = sd(Mean_dist_cent, na.rm = TRUE),
    median_breadth = median(Med_dist_cent , na.rm = TRUE),
    IQR_breadth = IQR(Med_dist_cent, na.rm = TRUE)
  )

### 2. Statistical tests summary of niche breadth between areas among regions ####

stat_summ_aoe <- function(aoe_df, reg_name = "All") {
  
  library(car)
  library(nortest)
  
  if (reg_name == "All") {
    
    aoe_data <- aoe_df 
    
  } else {
    
    aoe_data <- aoe_df %>% 
      filter(AoE_regs == reg_name)
    
  }
  
  x <- leveneTest(Mean_dist_cent ~ AoE_name, data = aoe_data)
  y <- ad.test(aoe_data$Mean_dist_cent)
  z <- kruskal.test(aoe_data$Mean_dist_cent ~ aoe_data$AoE_name)
  
  stat_summ <- data.frame(
    Region = reg_name,
    Data = rep("Mean-distance-to-centroid", 3),
    Test = c("Anderson-Darling normality test",
             "Levene's Test for Homogeneity of Variance (center = median)",
             "Kruskal-Wallis rank sum test"),
    Statistic = c("A", "F", "Chi-squared"),
    Value = c(format(z$statistic, digits = 5), format(x$`F value`[1], digits = 5), format(z$statistic, digits = 5)),
    p_0.05 = c(format(z$p.value, digits = 5), format(x$`Pr(>F)`[1], digits = 5), format(z$p.value, digits = 5)),
    Df = c(NA, paste(x$Df[1], x$Df[2], sep = "-"), z$parameter),
    Result = c(
      ifelse(z$p.value < 0.05, "Normality rejected", "Normality accepted"),
      ifelse(x$`Pr(>F)`[1] < 0.05, "Homoscedasticity rejected", "Homoscedasticity accepted"),
      ifelse(z$p.value < 0.05, "Significant differences between groups", "No significant differences between groups")
    ))
  
  return(stat_summ)
}

stat_summ_aoe(aoe_df = aoe_dfs, reg_name = "All")

# Not possible for Mesoamerica as it has only one AoE.
stat_summ_aoe(aoe_df = aoe_dfs, reg_name = as.character(unique(ca_db$AoE_regs))[6])

regs_vec <- c("All", unique(as.character(aoe_dfs$AoE_regs))[-6])
stat_summ_df <- purrr::map(regs_vec, stat_summ_aoe, aoe_df = aoe_dfs)
stat_summ_df <- do.call("rbind", stat_summ_df)

### 3. Kruskal wallis bar plots ####

aoe_dfs <- aoe_dfs %>% arrange(AoE_regs)
aoe_dfs$AoE_name <- factor(aoe_dfs$AoE_name, levels = unique(aoe_dfs$AoE_name), ordered = FALSE)

pv <- aoe_dfs %>% 
  group_by(AoE_name) %>%
  summarize(p.value = kruskal.test(aoe_dfs$Mean_dist_cent ~ aoe_dfs$AoE_name)$p.value)

aoe_dfs %>%
  ggplot(aes(x = AoE_name, y = Mean_dist_cent)) +
  #geom_point(aes(y = D), size = 0.5, alpha = 0.2, position = position_jitter(width = 0.1, height = 0.1)) +
  geom_boxplot(alpha = 0.5, aes(fill = AoE_regs)) +
  geom_text(data = pv, aes(x = 5, y = 0.3, label = paste0("Kruskal-Wallis p=", format(p.value, digits = 3)))) +
  scale_y_continuous(breaks = seq(0, 0.3, 0.05)) +
  #coord_cartesian(ylim = c(0,1)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  labs(y = "Mean-distance-to-centroid", x = "Areas of Endemism", fill = "Regions")

kruskal_aoe_plot <- function(x) {
  
  pv <- aoe_dfs %>% 
    filter(AoE_regs == x) %>%
    group_by(AoE_name) %>%
    summarize(p.value = kruskal.test(aoe_dfs$Mean_dist_cent ~ aoe_dfs$AoE_name)$p.value)
  
  aoe_dfs %>%
    filter(AoE_regs == x)  %>%
    ggplot(aes(x = AoE_name, y = Mean_dist_cent)) +
    #geom_point(aes(y = D), size = 0.5, alpha = 0.2, position = position_jitter(width = 0.1, height = 0.1)) +
    geom_boxplot(alpha = 0.5, aes(fill = AoE_regs)) +
    geom_text(data = pv, aes(x = 3, y = 0.3, label = paste0("Kruskal-Wallis p=", format(p.value, digits = 3)))) +
    scale_y_continuous(breaks = seq(0, 3, 0.05)) +
    #coord_cartesian(ylim = c(0,1)) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    labs(y = "Mean-distance-to-centroid", x = "Areas of Endemism", fill = "Regions")
  
}

regs_vec <- unique(aoe_dfs$AoE_regs)

kruskal_aoe_plot(x = regs_vec[6])


### 4. Correlation between betadisp_meandist and PCA_meandist ####
cor_meandist <- cor.test(aoe_dfs$mean_distance, aoe_dfs$Mean_dist_cent, method = "kendall")
ggpubr::ggscatter(aoe_dfs, x = "mean_distance", y = "Mean_dist_cent", 
                  add = "reg.line", add.params = list(color = "black"), conf.int = TRUE, 
                  cor.coef = TRUE, cor.method = "kendall",
                  xlab = "Betadisp_meandist", ylab = "PCA_meandist",
                  color = "grey", alpha = 0.5)

### 5. Statistical analysis ####

# Kruskall wallis is not applicable to our data. See Overlap_stats.R 

### 5.0 Loading data ####

# Loading niche properties table
niche_properties_allspp <- read.csv("./output/Breadth/2_Niche_properties_Allspp_PCA.csv")[-1]

ca_db <- read.csv("./data/biogeography/spp_listby_AoE.csv", stringsAsFactors = FALSE)[45:243,] %>%
  dplyr::select(-1)

# Geographical regions
AoE_regs <- read.csv("./data/biogeography/AoE_listby_Reg.csv", stringsAsFactors = FALSE)[-1] 
ca_db <- left_join(ca_db, AoE_regs, by = c("Area_name" = "AoE_name"))


# Preparing AoE data frame
aoe_extrac_df <- function(x) {
  
  spp_aoe_vec <- ca_db %>%
    filter(Area_name == x) %>%
    dplyr::select(Species) %>%
    pull()
  
  regs_vec <- ca_db %>%
    filter(Area_name == x) %>%
    dplyr::select(AoE_regs) %>%
    pull()
  
  aoe_df <- niche_properties_allspp %>%
    filter(NAME1 %in% spp_aoe_vec) %>%
    mutate(AoE_name = x,
           AoE_regs = regs_vec)
  
  return(aoe_df)
  
}

aoe_extrac_df(x = "D2CA1")

aoe_names <- ca_db %>% 
  dplyr::select(Area_name) %>%
  pull() %>%
  unique()

aoe_dfs <- purrr::map(.x = aoe_names, .f = aoe_extrac_df) 
aoe_dfs <- do.call("rbind", aoe_dfs)

aoe_dfs$AoE_regs <- factor(aoe_dfs$AoE_regs, levels = c("Mesoamerica", "Northern_Andes", "Guiana_centered", "Amazonia_centered", 
                                                        "ESA_Dry_Diagonal", "ESA_THDD", "ESA_Atlantic_Forest"), ordered = TRUE)

aoe_dfs$AoE_name <- factor(aoe_dfs$AoE_name, levels = unique(aoe_dfs$AoE_name), ordered = FALSE)

aoe_dfs <- aoe_dfs %>% arrange(AoE_regs)
### 5.1. Description ####

total <- describe(Mean_dist_cent ~ AoE_name, data = aoe_dfs)
regions <- describe(Mean_dist_cent ~ AoE_name, data = aoe_dfs)

### 5.2. Assumptions: Normality and Homoscedasticity ####

### 5.2.1. Regions ####

regs_vec <- c("All", unique(as.character(aoe_dfs$AoE_regs)))[-2] # Mesoamerica not showed; only one AoE

stat_summ_aoe2 <- function(aoe_df, reg_name = "All") {
  
  library(car)
  library(nortest)
  
  if (reg_name == "All") {
    
    aoe_data <- aoe_df 
    
  } else {
    
    aoe_data <- aoe_df %>% 
      filter(AoE_regs == reg_name)
    
  }
  
  x <- leveneTest(Mean_dist_cent ~ AoE_name, data = aoe_data)
  y <- ad.test(aoe_data$Mean_dist_cent)
  
  stat_summ <- data.frame(
    Region = reg_name,
    Data = rep("Breadth: Mean distance to centroid", 2),
    Test = c("Anderson-Darling normality test",
             "Levene's Test for Homogeneity of Variance (center = median)"),
    Statistic = c("A", "F"),
    Value = c(format(y$statistic[[1]], digits = 5), format(x$`F value`[1], digits = 5)),
    p_0.05 = c(format(y$p.value, digits = 5), format(x$`Pr(>F)`[1], digits = 5)),
    Df = c(NA, paste(x$Df[1], x$Df[2], sep = "-")),
    Result = c(
      ifelse(y$p.value < 0.05, "Normality rejected", "Normality accepted"),
      ifelse(x$`Pr(>F)`[1] < 0.05, "Homoscedasticity rejected", "Homoscedasticity accepted"))
  )
  
  return(stat_summ)
  
}
#stat_summ_aoe2(aoe_df = AoE_comp_df, reg_name = "Mesoamerica")
stat_summ_df <- purrr::map(regs_vec, stat_summ_aoe2, aoe_df = aoe_dfs)
stat_summ_df <- do.call("rbind", stat_summ_df)

### 5.2.2. AoEs ####

stat_summ_aoe3 <- function(aoe_df, aoe_name = "All") {
  
  library(car)
  library(nortest)
  
  if (aoe_name == "All") {
    
    aoe_data <- aoe_df 
    
  } else {
    
    aoe_data <- aoe_df %>% 
      filter(AoE_name == aoe_name)
    
  }
  
  y <- ad.test(aoe_data$Mean_dist_cent)
  
  stat_summ <- data.frame(
    AoE = aoe_name,
    Data = "Breadth: Mean distance to centroid",
    Test = c("Anderson-Darling normality test"),
    Statistic = c("A"),
    Value = c(format(y$statistic[[1]], digits = 5)),
    p_0.05 = c(format(y$p.value, digits = 5)),
    Result = c(ifelse(y$p.value < 0.05, "Normality rejected", "Normality accepted"))
  )
  
  return(stat_summ)
  
}

aoe_vec <- as.character(unique(ca_db$Area_name))

# Normality
aoe_dfs$AoE_name <- factor(aoe_dfs$AoE_name, levels = unique(aoe_dfs$AoE_name), ordered = FALSE)
nor.test(Mean_dist_cent ~ AoE_name, data = aoe_dfs, method = "SW")

# Levene homogeneity of variance
aoe_levene <- leveneTest(Mean_dist_cent ~ AoE_name, data = aoe_dfs)

### 5.3. Welch's heteroscedastic F with trimmed means and Winsorized variances ####

### AoEs
welchaoe <- welch.test(Mean_dist_cent ~ AoE_name, data = aoe_dfs, rate = 0.1)
# Pair comparison
paircomparison_aoe <- paircomp(welchaoe, adjust.method = "bonferroni", na.rm = TRUE)
colnames(paircomparison_aoe) <- str_replace_all(colnames(paircomparison_aoe), pattern = " ", replacement = "_")
colnames(paircomparison_aoe)[4] <- "No_difference"

diff_rejected <- paircomparison_aoe %>%
  filter(No_difference == "Reject")

### Regions
welchregs <- welch.test(Mean_dist_cent ~ AoE_regs, data = aoe_dfs, rate = 0.1, na.rm = TRUE)
# Pair comparison
paircomparison_regs <- paircomp(welchregs, adjust.method = "bonferroni")
# Pair comparison
colnames(paircomparison_regs) <- str_replace_all(colnames(paircomparison_regs), pattern = " ", replacement = "_")
colnames(paircomparison_regs)[4] <- "No_difference"

diff_rejected_regs <- paircomparison_regs %>%
  filter(No_difference == "Reject")

grid.table(diff_rejected_regs, theme)

grid.arrange(tableGrob(diff_rejected, theme = ttheme_minimal()))
grid.arrange(tableGrob(diff_rejected_regs, theme = ttheme_minimal()))

welchANOVA <- function(aoe_df, reg_name = "All") {
  
  if (reg_name == "All") {
    
    aoe_data <- aoe_df 
    
  } else {
    
    aoe_data <- aoe_df %>% 
      filter(AoE_regs == reg_name)
    
  }
  
  x <- welch.test(Mean_dist_cent ~ AoE_name, data = aoe_data, rate = 0.1, na.rm = TRUE)
  
  stat_summ <- data.frame(
    Region = reg_name,
    Data = "Breadth: Mean distance to centroid",
    Test = "Welch's ANOVA",
    Statistic = "F",
    Value = x$statistic,
    num_df = x$parameter[1],
    denom_df = x$parameter[2],
    p_0.05 = x$p.value,
    Result = ifelse(x$p.value < 0.05, "Difference is statistically significant.", "Difference is not statistically significant.")
  )
  
  return(stat_summ)
  
  
}

aoe_dfs %>% 
  filter(AoE_regs == regs_vec[3])

regs_vec <- c("All", unique(as.character(aoe_dfs$AoE_regs)))[-2]

wANOVA_summ_df <- purrr::map(regs_vec, welchANOVA, aoe_df = aoe_dfs)
wANOVA_summ_df <- do.call("rbind", wANOVA_summ_df)

# bread_df <- wANOVA_summ_df %>%
#   dplyr::select(2, 1, 3:9) %>%
#   rename(Niche_property = Data)
# write.csv(bread_df, "./output/Breadth/2_ANOVAS_table_breadth.csv")

# a <- format(bind_rows(list(pos_df, bread_df, over_df)), digits = 4)
# write.csv(a, "./output/AoE/2_ANOVAS_table_ALL.csv")


### 5.4. Boxplot ####

## Drawing stats box

# Labels
labels_vec <- vector(mode = "character", length = 7)

for (i in seq_along(wANOVA_summ_df$Region)) {
  
  labels_vec[i] <- paste(paste0(pull(wANOVA_summ_df[i,][1]), ":"), 
                         paste0(paste(paste(pull(wANOVA_summ_df[i,][4]), wANOVA_summ_df[i,][6], sep = ","), pull(format(wANOVA_summ_df[i,][5], digits = 5)), sep = " = "), ";"), 
                         paste("p.value", format(wANOVA_summ_df[i,][8], digits = 5), sep = " = "), sep = " ")
  
}

test_labels <- c(
  #"Welch's heteroscedastic F with trimmed means and Winsorized variances (alfa = 0.05)",
  labels_vec
)

legend_position <- sort(seq(0.5, 1, 0.07142857)[-1], decreasing = TRUE)

text_legend <- textGrob(label = test_labels, x = 0.6, 
                        y = legend_position, just = "left",
                        gp = gpar(fontsize = 10), check = TRUE) 

rect_position_y <- legend_position
rect_position_x <- rep(0.55, length(rect_position_y))

# jpn_palette <- brewer.pal(n = 7, name =  "Accent")
jpn_palette <- c("#FF0000", "#556A5B", "#50A45C", "#F2AD00", "#6D5B4E", "#C49647", "#5BBCD6")

square_legend <- rectGrob(x = rect_position_x, y = rect_position_y, 
                          width = 0.05, height = 0.05,
                          gp = gpar(fill = c("white", jpn_palette[-1]), alpha = 0.8))

grob_legend <- grobTree(text_legend, square_legend, name = "jp_legendgrob")

aoe_dfs$AoE_regs <- factor(aoe_dfs$AoE_regs, levels = c("Mesoamerica", "Northern_Andes", "Guiana_centered", "Amazonia_centered", 
                                                        "ESA_Dry_Diagonal", "ESA_THDD", "ESA_Atlantic_Forest"), ordered = TRUE)

aoe_dfs <- aoe_dfs %>% arrange(AoE_regs)

# x <- aoe_dfs %>%
#   group_by(AoE_regs, AoE_name) %>%
#   summarise(aoe_medD = median(Mean_dist_cent, na.rm = TRUE)) %>%
#   group_by(AoE_regs) %>%
#   mutate(ordering = as.numeric(str_replace(AoE_name, pattern = "D2CA([0-9]+)", replacement = "\\1"))) %>%
#   arrange(AoE_regs, ordering, aoe_medD) %>%
#   pull(AoE_name)
# 
# aoe_dfs$AoE_name <- factor(aoe_dfs$AoE_name, levels = x, ordered = FALSE)

x <- AoE_comp_df %>%
  group_by(AoE_regs, AoE_name) %>%
  summarise(aoe_medD = median(D)) %>%
  mutate(ordering = as.numeric(str_replace(AoE_name, pattern = "D2CA([0-9]+)", replacement = "\\1"))) %>%
  group_by(AoE_regs) %>%
  arrange(AoE_regs, ordering, aoe_medD) %>%
  pull(AoE_name)

aoe_dfs$AoE_name <- factor(aoe_dfs$AoE_name, levels = x, ordered = FALSE)

new_area_code <- aoe_dfs[, c(21, 22)] %>% distinct(AoE_name, AoE_regs)
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

breadth_anova_plot <- aoe_dfs %>%
  ggplot(aes(x = AoE_name, y = Mean_dist_cent)) +
  geom_boxplot(alpha = 0.8, aes(fill = AoE_regs)) +
  scale_fill_manual(values =  jpn_palette) +
  #geom_jitter(alpha = 0.2, size = 0.2, width = 0.2) +
  #geom_text(data = pv, aes(x = 3, y = 0.9, label = paste0("Kruskal-Wallis p=", format(p.value, digits = 3)))) +
  #geom_text(aes(x = 3, y = 0.9, label = paste0("Kruskal-Wallis p=", format(pv$p.value[1], digits = 3)))) +
  # coord_cartesian(ylim = c(0,0.35)) +
  scale_y_continuous(breaks = seq(0, 3, 0.05)) +
  scale_x_discrete(labels = area_code) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        legend.position = "none") +
  labs(y = "Niche breadth:\nMean distance to centroid", x = "Areas of Endemism", fill = "Regions") +
  annotation_custom(grob = grob_legend, 
                    xmin = -9, xmax = 10, 
                    ymin = 0.15, ymax = 0.35)

png("./figs/Boxplot_Breadth.png", height = 800, width = 1000, res = 120)
breadth_anova_plot 
dev.off()

#saveRDS(breadth_anova_plot, "./figs/Boxplot_Breadth.rds")
breadth_anova_plot <- readRDS("./figs/Boxplot_Breadth.rds")

### 5.5. Auxiliary figure: categories of niche breadth ####

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

breadth_df <- left_join(ca_db, niche_properties_allspp, by = c("Species" = "NAME1"))[-1] %>%
  dplyr::select(1,2,13,14,23)

area_code <- area_code[, c(1, 6)]
area_code$areacode <- factor(area_code$areacode, levels = unique(area_code$areacode), ordered = TRUE)

breadth_df <- left_join(breadth_df, area_code, by = c("Area_name" = "aoe"))
breadth_df$niche_category2 <- factor(breadth_df$areacode, levels = unique(breadth_df$niche_category2))

categ_aux <- breadth_df %>% 
  filter(!is.na(niche_category2)) %>% 
  ggplot(aes(x = areacode, fill = factor(niche_category2, levels = c("Narrow", "Medium", "Wide"), ordered = TRUE))) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels=scales::percent) +
  scale_fill_grey(start = 0.3) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  labs(fill = "Breadth",
       x = "Areas of Endemism",
       y = "Spp")

png("./figs/Categ_Breadth.png", height = 400, width = 1000, res = 120)
categ_aux
dev.off()

#saveRDS(categ_aux, "./figs/Auxiliary_Breadth.rds")
categ_aux <- readRDS("./figs/Auxiliary_Breadth.rds")

############### BEADTH FROM BETADISPER ####################
### 0. Loading data ####
# Data come from 03_Niche_properties_table 1 and 2 for environment 2
stat_df <- read.csv("./output/SPP/Allspp_Table_nicheproperties.csv")
ca_db <- read.csv("./data/biogeography/spp_listby_AoE.csv", stringsAsFactors = FALSE)[45:243,]

# Data come from 09_Ecospat_Results_ENDEMICSPP.R
# This is only to extract the column of geographical regions.
AoE_regs <- read.csv("./output/AoE/AoE_ECOSPAT_RES_areasonly.csv") %>%  filter(Analysis == "Eq", Alternative == "G") %>%
  dplyr::select(AoE_name, AoE_regs) %>%
  distinct()

ca_db <- left_join(ca_db, AoE_regs, by = c("Area_name" = "AoE_name"))

# Prepearing AoE data frame
aoe_extrac_df <- function(x) {
  
  spp_aoe_vec <- ca_db %>%
    filter(Area_name == x) %>%
    dplyr::select(Species) %>%
    pull()
  
  regs_vec <- ca_db %>%
    filter(Area_name == x) %>%
    dplyr::select(AoE_regs) %>%
    pull()
  
  aoe_df <- stat_df %>%
    filter(NAME1 %in% spp_aoe_vec) %>% 
    dplyr::select(NAME1, mean_distance, sd_distance) %>%
    mutate(AoE_name = x,
           AoE_regs = regs_vec)
  
  return(aoe_df)
  
}

aoe_extrac_df(x = "D2CA1")

aoe_names <- ca_db %>% 
  dplyr::select(Area_name) %>%
  pull() %>%
  unique()

aoe_dfs <- purrr::map(.x = aoe_names, .f = aoe_extrac_df) 
aoe_dfs <- do.call("rbind", aoe_dfs)

aoe_dfs$AoE_regs <- factor(aoe_dfs$AoE_regs, levels = c("Mesoamerica", "Northern_Andes", "Guiana_centered", "Amazonia_centered", 
                                                        "ESA_Dry_Diagonal", "ESA_THDD", "ESA_Atlantic_Forest"), ordered = TRUE)


### 1. Statistical summary ####

# By Aoe
aoe_dfs %>%
  group_by(AoE_name) %>%
  summarise(
    count = n(),
    mean_beadth = mean(mean_distance, na.rm = TRUE),
    sd_breadth = sd(mean_distance, na.rm = TRUE),
    median_breadth = median(mean_distance, na.rm = TRUE),
    IQR_breadth = IQR(mean_distance, na.rm = TRUE)
  )

# By region
aoe_dfs %>%
  group_by(AoE_regs) %>%
  summarise(
    count = n(),
    mean_beadth = mean(mean_distance, na.rm = TRUE),
    sd_breadth = sd(mean_distance, na.rm = TRUE),
    median_breadth = median(mean_distance, na.rm = TRUE),
    IQR_breadth = IQR(mean_distance, na.rm = TRUE)
  )

### 2. Statistical tests summary of niche breadth between areas among regions ####

stat_summ_aoe <- function(aoe_df, reg_name = "All") {
  
  library(car)
  library(nortest)
  
  if (reg_name == "All") {
    
    aoe_data <- aoe_df 
    
  } else {
    
    aoe_data <- aoe_df %>% 
      filter(AoE_regs == reg_name)
    
  }
  
  x <- leveneTest(mean_distance ~ AoE_name, data = aoe_data)
  y <- ad.test(aoe_data$mean_distance)
  z <- kruskal.test(aoe_data$mean_distance ~ aoe_data$AoE_name)
  
  stat_summ <- data.frame(
    Region = reg_name,
    Data = rep("Mean-distance-to-centroid", 3),
    Test = c("Anderson-Darling normality test",
             "Levene's Test for Homogeneity of Variance (center = median)",
             "Kruskal-Wallis rank sum test"),
    Statistic = c("A", "F", "Chi-squared"),
    Value = c(format(z$statistic, digits = 5), format(x$`F value`[1], digits = 5), format(z$statistic, digits = 5)),
    p_0.05 = c(format(z$p.value, digits = 5), format(x$`Pr(>F)`[1], digits = 5), format(z$p.value, digits = 5)),
    Df = c(NA, paste(x$Df[1], x$Df[2], sep = "-"), z$parameter),
    Result = c(
      ifelse(z$p.value < 0.05, "Normality rejected", "Normality accepted"),
      ifelse(x$`Pr(>F)`[1] < 0.05, "Homoscedasticity rejected", "Homoscedasticity accepted"),
      ifelse(z$p.value < 0.05, "Significant differences between groups", "No significant differences between groups")
    ))
  
  return(stat_summ)
  
  
}

stat_summ_aoe(aoe_df = aoe_dfs, reg_name = "All")

# Not possible for Mesoamerica as it has only one AoE.
stat_summ_aoe(aoe_df = aoe_dfs, reg_name = as.character(unique(ca_db$AoE_regs))[6])

regs_vec <- c("All", unique(as.character(aoe_dfs$AoE_regs))[-6])
stat_summ_df <- purrr::map(regs_vec, stat_summ_aoe, aoe_df = aoe_dfs)
stat_summ_df <- do.call("rbind", stat_summ_df)

### 3. Kruskal wallis bar plots ####

aoe_dfs$AoE_name <- factor(aoe_dfs$AoE_name, levels = unique(aoe_dfs$AoE_name), ordered = FALSE)
aoe_dfs <- aoe_dfs %>% arrange(AoE_regs)

pv <- aoe_dfs %>% 
  group_by(AoE_name) %>%
  summarize(p.value = kruskal.test(aoe_dfs$mean_distance ~ aoe_dfs$AoE_name)$p.value)

aoe_dfs %>%
  ggplot(aes(x = AoE_name, y = mean_distance)) +
  #geom_point(aes(y = D), size = 0.5, alpha = 0.2, position = position_jitter(width = 0.1, height = 0.1)) +
  geom_boxplot(alpha = 0.5, aes(fill = AoE_regs)) +
  geom_text(data = pv, aes(x = 5, y = 2.9, label = paste0("Kruskal-Wallis p=", format(p.value, digits = 3)))) +
  scale_y_continuous(breaks = seq(0, 3, 0.5)) +
  #coord_cartesian(ylim = c(0,1)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  labs(y = "Mean-distance-to-centroid", x = "Areas of Endemism", fill = "Regions")

kruskal_aoe_plot <- function(x) {
  
  pv <- aoe_dfs %>% 
    filter(AoE_regs == x) %>%
    group_by(AoE_name) %>%
    summarize(p.value = kruskal.test(aoe_dfs$mean_distance ~ aoe_dfs$AoE_name)$p.value)
  
  aoe_dfs %>%
    filter(AoE_regs == x)  %>%
    ggplot(aes(x = AoE_name, y = mean_distance)) +
    #geom_point(aes(y = D), size = 0.5, alpha = 0.2, position = position_jitter(width = 0.1, height = 0.1)) +
    geom_boxplot(alpha = 0.5, aes(fill = AoE_regs)) +
    geom_text(data = pv, aes(x = 3, y = 3, label = paste0("Kruskal-Wallis p=", format(p.value, digits = 3)))) +
    scale_y_continuous(breaks = seq(0, 3, 0.5)) +
    #coord_cartesian(ylim = c(0,1)) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    labs(y = "Mean-distance-to-centroid", x = "Areas of Endemism", fill = "Regions")
  
}

regs_vec <- unique(aoe_dfs$AoE_regs)

kruskal_aoe_plot(x = regs_vec[7])

