## Script for statistical analysis of: Niche position with Environment 2 (Spp + AoEs Env)
## Main objective: look for patterns of niche properties similarities and differences among
## endemic species inside each Areas of Endemism.
## This file is done just to unify the code for the statistical analysis on niche position.

# Libraries
library(dplyr)
library(ggplot2)
library(ggpubr)
library(purrr)
library(RColorBrewer)
library(onewaytests)

### 0. Loading data ####

# position_df <- read.csv("./output/Position/2_sppcentroid_mean_dist_to_AoEcentroid_wb.csv")
# glimpse(position_df)
# 
# position_df %>%
#   filter(Group_comp == "Within")

position_df2 <- read.csv("./output/Position/2_sppcentroid_to_aoes_centroid_wb.csv")[-1] %>%
  filter(Group_comp == "Within")
glimpse(position_df2)

ca_db <- read.csv("./data/biogeography/spp_listby_AoE.csv", stringsAsFactors = FALSE)[45:243,] %>%
  dplyr::select(-1)

# Geographical regions
AoE_regs <- read.csv("./data/biogeography/AoE_listby_Reg.csv", stringsAsFactors = FALSE)[-1] 
ca_db <- left_join(ca_db, AoE_regs, by = c("Area_name" = "AoE_name"))

aoe_dfs <- left_join(position_df2, AoE_regs, by = c("Item_group" = "AoE_name"))

aoe_dfs$AoE_regs <- factor(aoe_dfs$AoE_regs, levels = c("Mesoamerica", "Northern_Andes", "Guiana_centered", "Amazonia_centered", 
                                                        "ESA_Dry_Diagonal", "ESA_THDD", "ESA_Atlantic_Forest"), ordered = TRUE)


### 1. Statistical summary ####
glimpse(aoe_dfs)

# By Aoe
aoe_dfs %>%
  group_by(Centroid_Group) %>%
  summarise(
    count = n(),
    pcamean_position = mean(CentroidDistance, na.rm = TRUE),
    pcasd_position = sd(CentroidDistance, na.rm = TRUE),
    median_position = median(CentroidDistance, na.rm = TRUE),
    IQR_position = IQR(CentroidDistance, na.rm = TRUE)
  )

# By region
aoe_dfs %>%
  group_by(AoE_regs) %>%
  summarise(
    count = n(),
    pcamean_position = mean(CentroidDistance, na.rm = TRUE),
    pcasd_position = sd(CentroidDistance, na.rm = TRUE),
    median_position = median(CentroidDistance, na.rm = TRUE),
    IQR_position = IQR(CentroidDistance, na.rm = TRUE)
  )

### 2. Statistical tests summary of niche position between areas among regions ####

stat_summ_aoe <- function(aoe_df, reg_name = "All") {
  
  library(car)
  library(nortest)
  
  if (reg_name == "All") {
    
    aoe_data <- aoe_df 
    
  } else {
    
    aoe_data <- aoe_df %>% 
      filter(AoE_regs == reg_name)
    
  }
  
  x <- leveneTest(CentroidDistance ~ Centroid_Group, data = aoe_data)
  y <- ad.test(aoe_data$CentroidDistance)
  z <- kruskal.test(aoe_data$CentroidDistance ~ aoe_data$Centroid_Group)
  
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

stat_summ_aoe(aoe_df = aoe_dfs, reg_name = as.character(unique(ca_db$AoE_regs))[-6])

regs_vec <- c("All", unique(as.character(aoe_dfs$AoE_regs))[-6])
stat_summ_df <- purrr::map(regs_vec, stat_summ_aoe, aoe_df = aoe_dfs)
stat_summ_df <- do.call("rbind", stat_summ_df)

### 3. Kruskal wallis bar plots ####

aoe_dfs$AoE_regs <- factor(aoe_dfs$AoE_regs, levels = c("Mesoamerica", "Northern_Andes", "Guiana_centered", "Amazonia_centered", 
                                                        "ESA_Dry_Diagonal", "ESA_THDD", "ESA_Atlantic_Forest"), ordered = TRUE)

aoe_dfs <- aoe_dfs %>% arrange(AoE_regs)

aoe_dfs$Centroid_Group <- factor(aoe_dfs$Centroid_Group, levels = unique(aoe_dfs$Centroid_Group), ordered = FALSE)


pv <- aoe_dfs %>% 
  group_by(Centroid_Group) %>%
  summarize(p.value = kruskal.test(aoe_dfs$CentroidDistance ~ aoe_dfs$Centroid_Group)$p.value)

aoe_dfs %>%
  ggplot(aes(x = Centroid_Group, y = CentroidDistance)) +
  #geom_point(aes(y = D), size = 0.5, alpha = 0.2, position = position_jitter(width = 0.1, height = 0.1)) +
  geom_boxplot(alpha = 0.5, aes(fill = AoE_regs)) +
  geom_text(data = pv, aes(x = 5, y = 0.4, label = paste0("Kruskal-Wallis p=", format(p.value, digits = 3)))) +
  scale_y_continuous(breaks = seq(0, 0.5, 0.05)) +
  #coord_cartesian(ylim = c(0,1)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  labs(y = "Position: Mean-distance-to-centroid", x = "Areas of Endemism", fill = "Regions")

kruskal_aoe_plot <- function(x) {
  
  pv <- aoe_dfs %>% 
    filter(AoE_regs == x) %>%
    group_by(Centroid_Group) %>%
    summarize(p.value = kruskal.test(aoe_dfs$CentroidDistance ~ aoe_dfs$Centroid_Group)$p.value)
  
  aoe_dfs %>%
    filter(AoE_regs == x)  %>%
    ggplot(aes(x = Centroid_Group, y = CentroidDistance)) +
    #geom_point(aes(y = D), size = 0.5, alpha = 0.2, position = position_jitter(width = 0.1, height = 0.1)) +
    geom_boxplot(alpha = 0.5, aes(fill = AoE_regs)) +
    geom_text(data = pv, aes(x = 3, y = 0.4, label = paste0("Kruskal-Wallis p=", format(p.value, digits = 3)))) +
    scale_y_continuous(breaks = seq(0, 5, 0.05)) +
    #coord_cartesian(ylim = c(0,1)) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    labs(y = "Mean-distance-to-centroid", x = "Areas of Endemism", fill = "Regions")
  
}

regs_vec <- unique(aoe_dfs$AoE_regs)

kruskal_aoe_plot(x = regs_vec[6])

### 4. Statistical analysis ####

# Kruskall wallis is not applicable to our data. See Overlap_stats.R 

### 4.0 Loading data ####

position_df <- read.csv("./output/Position/2_sppcentroid_to_aoes_centroid_wb.csv")[-1] %>%
  filter(Group_comp == "Within")
glimpse(position_df)

ca_db <- read.csv("./data/biogeography/spp_listby_AoE.csv", stringsAsFactors = FALSE)[45:243,] %>%
  dplyr::select(-1)

# Geographical regions
AoE_regs <- read.csv("./data/biogeography/AoE_listby_Reg.csv", stringsAsFactors = FALSE)[-1] 
ca_db <- left_join(ca_db, AoE_regs, by = c("Area_name" = "AoE_name"))

# position_df2 <- position_df %>%
#   group_by(Group_comp, Item_group, Centroid_Group) %>%
#   summarise(Mean_dist_cent = mean(CentroidDistance), 
#             Sd_dist_cent = sd(CentroidDistance),
#             Med_dist_cent = median(CentroidDistance),
#             IQR_dist_cent = IQR(CentroidDistance)) 
# 
# position_fivenum_df <- position_df %>%
#   group_by(Group_comp, Item_group, Centroid_Group) %>%
#   dplyr::select(CentroidDistance) %>%
#   group_modify(~ {
#     .x %>%
#       purrr::map_dfc(fivenum) %>%
#       mutate(pcvar = c("min", "Q1", "median", "Q3", "max"))
#   }) %>%
#   pivot_wider(names_from = pcvar, values_from = CentroidDistance)
# 
# position_df2 <- left_join(position_df2, position_fivenum_df, by = c("Group_comp", "Item_group", "Centroid_Group"))
# 
# position_df2 <- position_df2 %>%
#   filter(Group_comp == "Within")
# 
# aoe_dfs <- left_join(position_df2, AoE_regs, by = c("Item_group" = "AoE_name"))
# colnames(aoe_dfs)[3] <- "AoE_name"

aoe_dfs <- left_join(position_df, AoE_regs, by = c("Item_group" = "AoE_name"))
colnames(aoe_dfs)[3] <- "AoE_name"

aoe_dfs$AoE_regs <- factor(aoe_dfs$AoE_regs, levels = c("Mesoamerica", "Northern_Andes", "Guiana_centered", "Amazonia_centered", 
                                                        "ESA_Dry_Diagonal", "ESA_THDD", "ESA_Atlantic_Forest"), ordered = TRUE)

aoe_dfs <- aoe_dfs %>% arrange(AoE_regs)

aoe_dfs$AoE_name <- factor(aoe_dfs$AoE_name, levels = unique(aoe_dfs$AoE_name), ordered = FALSE)

### 4.1. Description ####

total <- describe(CentroidDistance ~ AoE_name, data = aoe_dfs)
regions <- describe(CentroidDistance ~ AoE_regs, data = aoe_dfs)

### 4.2. Assumptions: Normality and Homoscedasticity ####

### 4.2.1. Regions ####

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
  
  x <- leveneTest(CentroidDistance ~ AoE_name, data = aoe_data)
  y <- ad.test(aoe_data$CentroidDistance)
  
  stat_summ <- data.frame(
    Region = reg_name,
    Data = rep("Position: Spp-AoE centroids distance", 2),
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

### 4.2.2. AoEs ####

stat_summ_aoe3 <- function(aoe_df, aoe_name = "All") {
  
  library(car)
  library(nortest)
  
  if (aoe_name == "All") {
    
    aoe_data <- aoe_df 
    
  } else {
    
    aoe_data <- aoe_df %>% 
      filter(AoE_name == aoe_name)
    
  }
  
  y <- ad.test(aoe_data$CentroidDistance)
  
  stat_summ <- data.frame(
    AoE = aoe_name,
    Data = "Position: Spp-AoE centroids distance",
    Test = c("Anderson-Darling normality test"),
    Statistic = c("A"),
    Value = c(format(y$statistic[[1]], digits = 5)),
    p_0.05 = c(format(y$p.value, digits = 5)),
    Result = c(ifelse(y$p.value < 0.05, "Normality rejected", "Normality accepted"))
  )
  
  return(stat_summ)
  
}
aoe_vec <- as.character(unique(aoe_dfs$AoE_name))

# Normality
aoe_dfs$AoE_name <- factor(aoe_dfs$AoE_name, levels = unique(aoe_dfs$AoE_name), ordered = FALSE)
nor.test(CentroidDistance ~ AoE_name, data = aoe_dfs, method = "SW")

# Levene homogeneity of variance
aoe_levene <- leveneTest(CentroidDistance ~ AoE_name, data = aoe_dfs)

### 4.3. Welch's heteroscedastic F with trimmed means and Winsorized variances ####

welchaoe <- welch.test(CentroidDistance ~ AoE_name, data = aoe_dfs, rate = 0.1)
paircomparison_aoe <- paircomp(welchaoe, adjust.method = "bonferroni", na.rm = TRUE)
colnames(paircomparison_aoe) <- str_replace_all(colnames(paircomparison_aoe), pattern = " ", replacement = "_")
colnames(paircomparison_aoe)[4] <- "No_difference"
diff_rejected <- paircomparison_aoe %>%
  filter(No_difference == "Reject")
grid.arrange(tableGrob(diff_rejected, theme = ttheme_minimal()))


welchregs <- welch.test(CentroidDistance ~ AoE_regs, data = aoe_dfs, rate = 0.1, na.rm = TRUE)
paircomparison_regs <- paircomp(welchregs, adjust.method = "bonferroni")
colnames(paircomparison_regs) <- str_replace_all(colnames(paircomparison_regs), pattern = " ", replacement = "_")
colnames(paircomparison_regs)[4] <- "No_difference"
diff_rejected_regs <- paircomparison_regs %>%
  filter(No_difference == "Reject")
grid.arrange(tableGrob(diff_rejected_regs, theme = ttheme_minimal()))

welchANOVA <- function(aoe_df, reg_name = "All") {
  
  if (reg_name == "All") {
    
    aoe_data <- aoe_df 
    
  } else {
    
    aoe_data <- aoe_df %>% 
      filter(AoE_regs == reg_name)
    
  }
  
  x <- welch.test(CentroidDistance ~ AoE_name, data = aoe_data, rate = 0.1)
  
  stat_summ <- data.frame(
    Region = reg_name,
    Data = "Position: Spp-AoE centroids distance",
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

regs_vec <- c("All", unique(as.character(aoe_dfs$AoE_regs)))[-2]

wANOVA_summ_df <- purrr::map(regs_vec, welchANOVA, aoe_df = aoe_dfs)
wANOVA_summ_df <- do.call("rbind", wANOVA_summ_df)

# pos_df <- wANOVA_summ_df %>%
#   dplyr::select(2, 1, 3:9) %>%
#   rename(Niche_property = Data)
# #write.csv(pos_df, "./output/Position/2_ANOVAS_table_position.csv")

## Complete ANOVAS table ####

# pos_df <- read.csv("./output/Position/2_ANOVAS_table_position.csv")[-1]
# bread_df <- read.csv("./output/Breadth/2_ANOVAS_table_breadth.csv")[-1]
# over_df <- read.csv("./output/Niche_overlap/2_ANOVAS_table_overlap.csv")[-1]
# 
# a <- format(bind_rows(list(pos_df, bread_df, over_df)), digits = 4)
# write.csv(a, "./output/AoE/2_ANOVAS_table_ALLREGIONS.csv")

### 4.4. Boxplot ####

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

# jpn_palette <- brewer.pal(n = 7, name = "Accent")
jpn_palette <- c("#FF0000", "#556A5B", "#50A45C", "#F2AD00", "#6D5B4E", "#C49647", "#5BBCD6")

square_legend <- rectGrob(x = rect_position_x, y = rect_position_y, 
                          width = 0.05, height = 0.05,
                          gp = gpar(fill = c("white", jpn_palette[-1]), alpha = 0.8))

grob_legend <- grobTree(text_legend, square_legend, name = "jp_legendgrob")

aoe_dfs$AoE_regs <- factor(aoe_dfs$AoE_regs, levels = c("Mesoamerica", "Northern_Andes", "Guiana_centered", "Amazonia_centered", 
                                                        "ESA_Dry_Diagonal", "ESA_THDD", "ESA_Atlantic_Forest"), ordered = TRUE)

aoe_dfs <- aoe_dfs %>% arrange(AoE_regs)

new_area_code <- aoe_dfs[, c(3, 7)] %>% distinct(AoE_name, AoE_regs)
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

position_anova_plot <- aoe_dfs %>%
  ggplot(aes(x = AoE_name, y = CentroidDistance)) +
  geom_boxplot(aes(fill = AoE_regs), alpha = 0.8) +
  scale_fill_manual(values =  jpn_palette) +
  #geom_jitter(alpha = 0.2, size = 0.2, width = 0.2) +
  #geom_text(data = pv, aes(x = 3, y = 0.9, label = paste0("Kruskal-Wallis p=", format(p.value, digits = 3)))) +
  #geom_text(aes(x = 3, y = 0.9, label = paste0("Kruskal-Wallis p=", format(pv$p.value[1], digits = 3)))) +
  #coord_cartesian(ylim = c(0,0.5)) +
  scale_y_continuous(breaks = seq(0, 3, 0.05)) +
  scale_x_discrete(labels = area_code) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        legend.position = "none") +
  labs(y = "Niche Position:\nSpp-AoE centroids' distance", x = "Areas of Endemism", fill = "Regions") +
  annotation_custom(grob = grob_legend, 
                    xmin = -9, xmax = 10, 
                    ymin = 0.25, ymax = 0.5)


# g <- ggplot_build(position_anova_plot)
# pull(unique(g$data[[1]]["fill"]))
# [1] "#7FC97F" "#BEAED4" "#FDC086" "#FFFF99" "#386CB0" "#F0027F" "#BF5B17"

png("./figs/Boxplot_Position.png", height = 800, width = 1000, res = 120)
position_anova_plot 
dev.off()

#saveRDS(position_anova_plot, "./figs/Boxplot_Position.rds")
position_anova_plot <- readRDS("./figs/Boxplot_Position.rds")


### 4.5. Auxiliary figure: position categories ####

a <- aoe_dfs  %>%
  ggplot(aes(x = 1, y = CentroidDistance)) +
  geom_boxplot()

a <- ggplot_build(a)
a[1]
# lower = 0.04992205
# middle = 0.09053683
# upper = 0.1339504 
# min = 0.01078003
# max = 0.216511

aoe_dfs <- aoe_dfs %>%
  mutate(position_category = case_when(
    CentroidDistance <= 0.04992205 ~ "Near", # lower = 0.04992205
    CentroidDistance > 0.04992205 &  CentroidDistance <= 0.1339504  ~ "Middling", # upper = 0.1339504
    CentroidDistance > 0.1339504 ~ "Far" # upper = 0.1339504
  ))

aoe_dfs <- left_join(aoe_dfs, area_code, by = c("AoE_name" = "aoe"))
aoe_dfs$position_category <- factor(aoe_dfs$position_category, levels = c("Near", "Middling", "Far"), ordered = TRUE)
aoe_dfs$areacode <- factor(aoe_dfs$areacode, levels = unique(aoe_dfs$areacode), ordered = TRUE)

categ_pos <- aoe_dfs %>% 
  #filter(!is.na(niche_category2)) %>% 
  ggplot(aes(x = areacode, fill = position_category)) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels=scales::percent) +
  scale_fill_grey(start = 0.3) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  labs(fill = "Position",
       x = "Areas of Endemism",
       y = "Distance from Spp centroid\nto AoE centroid")

png("./figs/Categ_Breadth.png", height = 400, width = 1000, res = 120)
categ_pos
dev.off()

#saveRDS(categ_pos, "./figs/Auxiliary_Position.rds")
categ_pos <- readRDS("./figs/Auxiliary_Position.rds")

### ALL Boxplots FIGURE ####

library(patchwork)

breadth_anova_plot <- readRDS("./figs/Boxplot_Breadth.rds")
position_anova_plot <- readRDS("./figs/Boxplot_Position.rds")
overlap_anova_plot <- readRDS("./figs/Boxplot_OverlapD.rds")

png("./figs/Boxplot_ALL_ANOVAS.png", width = 1500, height = 1500, res = 150)
print((breadth_anova_plot / position_anova_plot / overlap_anova_plot) +
        plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face = "bold")))
dev.off()

### All Auxiliary figure ####

library(patchwork)

categ_aux <- readRDS("./figs/Auxiliary_Breadth.rds")
categ_pos <- readRDS("./figs/Auxiliary_Position.rds")
categ_over <- readRDS("./figs/Auxiliary_Overlap.rds")

png("./figs/Auxiliary_ALL.png", width = 1500, height = 1500, res = 150)
print((categ_aux / categ_pos / categ_over) +
        plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face = "bold")))
dev.off()

