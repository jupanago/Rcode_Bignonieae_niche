## Script for statistical analysis of: OVERLAP with Environment 1 (Spp only)
## Main objective: look for patterns of niche properties similarities and differences among
## endemic species inside each Areas of Endemism.

# Libraries

library(tidyverse)
library(onewaytests)
library(ggplot2)
library(grid)
library(gridExtra)
library(viridis)

### 0. Data ####
# Only the data from the Equivalence-Greater test table are used to analyse the overlap D.

AoE_comp_df <- read.csv("./output/AoE/AoE_ECOSPAT_RES_areasonly.csv") %>%  filter(Analysis == "Eq", Alternative == "G")

AoE_comp_df$AoE_regs <- factor(AoE_comp_df$AoE_regs, levels = c("Mesoamerica", "Northern_Andes", "Guiana_centered", "Amazonia_centered", 
                                                                "ESA_Dry_Diagonal", "ESA_THDD", "ESA_Atlantic_Forest"), ordered = TRUE)
AoE_comp_df$D_classes <- factor(AoE_comp_df$D_classes, levels = c("Limited", "Low", "Moderate", "High", "Very_high"), ordered = TRUE)
AoE_comp_df$I_classes <- factor(AoE_comp_df$I_classes, levels = c("Limited", "Low", "Moderate", "High", "Very_high"), ordered = TRUE)
AoE_comp_df$Geo_classes <- factor(AoE_comp_df$Geo_classes, levels = c("Disjunct", "Low", "Medium", "High", "Very_high"), ordered = TRUE)
AoE_comp_df <- AoE_comp_df %>% arrange(AoE_regs)

### 1. Statistical summary of AoE D overlap ####

AoE_comp_df %>%
  group_by(AoE_name) %>%
  summarise(
    count = n(),
    mean_D = mean(D, na.rm = TRUE),
    sd_D = sd(D, na.rm = TRUE),
    median_D = median(D, na.rm = TRUE),
    IQR_D = IQR(D, na.rm = TRUE),
    mean_G = mean(Geo_overlap, na.rm = TRUE),
    sd_G = sd(Geo_overlap, na.rm = TRUE),
    median_G = median(Geo_overlap, na.rm = TRUE),
    IQR_G = IQR(Geo_overlap, na.rm = TRUE)
  )


### 2. Statistical summary of Overlap between areas among regions ####

stat_summ_aoe <- function(aoe_df, reg_name = "All") {
  
  library(car)
  library(nortest)
  
  if (reg_name == "All") {
    
    aoe_data <- aoe_df 
    
  } else {
    
    aoe_data <- aoe_df %>% 
      filter(AoE_regs == reg_name)
    
  }
  
  x <- leveneTest(D ~ AoE_name, data = aoe_data)
  y <- ad.test(aoe_data$D)
  z <- kruskal.test(aoe_data$D ~ aoe_data$AoE_name)
  
  stat_summ <- data.frame(
    Region = reg_name,
    Data = rep("Schoener D", 3),
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
stat_summ_aoe(aoe_df = AoE_comp_df, reg_name = "All")

regs_vec <- c("All", unique(as.character(AoE_comp_df$AoE_regs)))[-2]

stat_summ_df <- purrr::map(regs_vec, stat_summ_aoe, aoe_df = AoE_comp_df)
stat_summ_df <- do.call("rbind", stat_summ_df)

### 3. Kruskal wallis bar plots ####

# KRUSKAL-WALLIS REJECTED: Assumption of homocedasticity not fulfilled. The shapes of the D values distributions are different.
# Ratio between the greatest and smallest variance = 72.8. The rule forbids to apply non-parametric equivalents of ANOVA
# when this ratio is greater than 4:1. Checking variance:
x <- AoE_comp_df %>%
  group_by(AoE_name) %>%
  summarise(variance = var(D)) %>%
  arrange(desc(variance)) %>% 
  dplyr::select(AoE_name, variance)
x$variance[1] / x$variance[28]

# Refer to table above: stat_sum_df and check that the levene test rejects the null hypothesis of homogeneity of variances for ALL
# and D2CA0.
# Check here how the shape of the distribution of D values differs for each area. 
AoE_comp_df %>% 
  filter(AoE_name != "D2CA0") %>%
  ggplot(mapping = aes(x = D)) + 
  geom_freqpoly(mapping = aes(color = AoE_name), binwidth = 0.05) +
  facet_wrap(~AoE_regs)

# Solution: Use Welch’s heteroscedastic F test with trimmed means and Winsorized variances (Dag et al 2018; Keldman et al 2008)


### 4. Statiscal Analysis ####

# IMPORTANT NOTE: Ambiguous species will be mantained. If they are eliminated from the analysis, seven areas are lost:
# [1] "D2CA17" "D2CA19" "D2CA6"  "D2CA8"  "D2CA9"  "D2CA13" "D2CA21"
{ca_db <- read.csv("./data/biogeography/spp_listby_AoE.csv", stringsAsFactors = FALSE)[45:243,] 
dup <- ca_db$Species[duplicated(ca_db$Species)]

spp_code <- read.csv("./output/AoE/Table_endemicspp_code_SDumont.csv")

dup_spp_code <- spp_code %>%
  filter(spp_names %in% dup) %>%
  dplyr::select(ID) %>%
  pull()

ca_db %>%
  filter(Species %in% dup) %>%
  group_by(Species) %>%
  summarise(count = n()) %>%
  arrange(desc(count))

a <- AoE_comp_df %>%
  group_by(AoE_name) %>%
  summarise(count_conf = n()) %>%
  filter(count_conf <= 7)

b <- AoE_comp_df %>%
  filter(!Sp.1 %in% dup_spp_code, !Sp.2 %in% dup_spp_code) %>%
  group_by(AoE_name) %>%
  summarise(count_noconf = n()) %>%
  filter(count_noconf <= 7)

full_join(b, a) %>%
  mutate(difference = count_conf - count_noconf)

c <- AoE_comp_df %>%
  filter(!Sp.1 %in% dup_spp_code, !Sp.2 %in% dup_spp_code) %>%
  dplyr::select(AoE_name) %>%
  unique() %>%
  pull() %>%
  as.character()

d <- AoE_comp_df %>%
  dplyr::select(AoE_name) %>%
  unique() %>%
  pull() %>%
  as.character()

length(c)
length(d)

d[!d %in% c] # Areas lost if conflicting species are eliminated
} 

library(onewaytests)

### 4.1. Description ####

total <- describe(D ~ AoE_name, data = AoE_comp_df)
regions <- describe(D ~ AoE_regs, data = AoE_comp_df)

### 4.2. Assumptions: Normality and Homoscedasticity ####

### 4.2.1. Regions ####

regs_vec <- c("All", unique(as.character(AoE_comp_df$AoE_regs)))[-2] # Mesoamerica not showed; only one AoE

stat_summ_aoe2 <- function(aoe_df, reg_name = "All") {
  
  library(car)
  library(nortest)
  
  if (reg_name == "All") {
    
    aoe_data <- aoe_df 
    
  } else {
    
    aoe_data <- aoe_df %>% 
      filter(AoE_regs == reg_name)
    
  }
  
  x <- leveneTest(D ~ AoE_name, data = aoe_data)
  y <- ad.test(aoe_data$D)
  
  stat_summ <- data.frame(
    Region = reg_name,
    Data = rep("Schoener D", 2),
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
stat_summ_df <- purrr::map(regs_vec, stat_summ_aoe2, aoe_df = AoE_comp_df)
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
  
  y <- ad.test(aoe_data$D)
  
  stat_summ <- data.frame(
    AoE = aoe_name,
    Data = "Schoener D",
    Test = c("Anderson-Darling normality test"),
    Statistic = c("A"),
    Value = c(format(y$statistic[[1]], digits = 5)),
    p_0.05 = c(format(y$p.value, digits = 5)),
    Result = c(ifelse(y$p.value < 0.05, "Normality rejected", "Normality accepted"))
  )
  
  return(stat_summ)
  
}

aoe_vec <- as.character(unique(AoE_comp_df$AoE_name))

# Adtest of normality 
aoe_noadtest <- AoE_comp_df %>% # Areas with less than n<7 cannot be submitet to ad.test()
  group_by(AoE_name) %>%
  summarise(count_conf = n()) %>%
  filter(count_conf <= 7) %>%
  dplyr::select(AoE_name) %>%
  pull() %>%
  as.character()

aoe_vec <- aoe_vec[!aoe_vec %in% aoe_nolevene]
aoe_vec <- c("All", aoe_vec)

stat_summ_aoe3(aoe_df = AoE_comp_df, aoe_name = "All")

stat_summ_df2 <- purrr::map(aoe_vec, stat_summ_aoe3, aoe_df = AoE_comp_df)
stat_summ_df2 <- do.call("rbind", stat_summ_df2) %>%
  mutate(n = as.numeric(ifelse(AoE == "All", -0.5, str_replace(AoE, pattern = "D2CA([0-9]+)", replacement = "\\1")))) %>%
  arrange(n) %>%
  dplyr::select(-n)

nor.test(D ~ AoE_name, data = AoE_comp_df, method = "LT")

# Levene homogeneity of variance
aoe_levene <- leveneTest(D ~ AoE_name, data = AoE_comp_df)

### 4.3. Welch's heteroscedastic F with trimmed means and Winsorized variances ####

welchaoe <- welch.test(D ~ AoE_name, data = AoE_comp_df, rate = 0.2)
paircomparison_aoe <- paircomp(welchaoe, adjust.method = "bonferroni")
colnames(paircomparison_aoe) <- str_replace_all(colnames(paircomparison_aoe), pattern = " ", replacement = "_")
colnames(paircomparison_aoe)[4] <- "No_difference"
diff_rejected <- paircomparison_aoe %>%
  filter(No_difference == "Reject")
grid.arrange(tableGrob(head(diff_rejected, 15), theme = ttheme_minimal()))


welchregs <- welch.test(D ~ AoE_regs, data = AoE_comp_df, rate = 0.2)
paircomparison_regs <- paircomp(welchregs, adjust.method = "bonferroni")

welchANOVA <- function(aoe_df, reg_name = "All") {
  
  if (reg_name == "All") {
    
    aoe_data <- aoe_df 
    
  } else {
    
    aoe_data <- aoe_df %>% 
      filter(AoE_regs == reg_name)
    
  }
  
  x <- welch.test(D ~ AoE_name, data = aoe_data, rate = 0.2)
  
  stat_summ <- data.frame(
    Region = reg_name,
    Data = "Schoener D",
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

regs_vec <- c("All", unique(as.character(AoE_comp_df$AoE_regs)))[-2]

wANOVA_summ_df <- purrr::map(regs_vec, welchANOVA, aoe_df = AoE_comp_df)
wANOVA_summ_df <- do.call("rbind", wANOVA_summ_df)

# over_df <- wANOVA_summ_df %>%
#   dplyr::select(2, 1, 3:9) %>%
#   rename(Niche_property = Data)
# #write.csv(over_df, "./output/Niche_overlap/2_ANOVAS_table_overlap.csv")

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
                          gp = gpar(fill = c("white",jpn_palette[-1]), alpha = 0.8))

grob_legend <- grobTree(text_legend, square_legend, name = "jp_legendgrob")

x <- AoE_comp_df %>%
  group_by(AoE_regs, AoE_name) %>%
  summarise(aoe_medD = median(D)) %>%
  mutate(ordering = as.numeric(str_replace(AoE_name, pattern = "D2CA([0-9]+)", replacement = "\\1"))) %>%
  group_by(AoE_regs) %>%
  arrange(AoE_regs, ordering, aoe_medD) %>%
  pull(AoE_name)

AoE_comp_df$AoE_name <- factor(AoE_comp_df$AoE_name, levels = x, ordered = FALSE)

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


overlap_anova_plot <- AoE_comp_df %>%
  ggplot(aes(x = AoE_name, y = D)) +
  geom_boxplot(alpha = 0.8, aes(fill = AoE_regs)) +
  scale_fill_manual(values =  jpn_palette) +
  #geom_jitter(alpha = 0.2, size = 0.2, width = 0.2) +
  #geom_text(data = pv, aes(x = 3, y = 0.9, label = paste0("Kruskal-Wallis p=", format(p.value, digits = 3)))) +
  #geom_text(aes(x = 3, y = 0.9, label = paste0("Kruskal-Wallis p=", format(pv$p.value[1], digits = 3)))) +
  coord_cartesian(ylim = c(0,1)) +
  scale_y_continuous(breaks = seq(0, 3, 0.2)) +
  scale_x_discrete(labels = area_code) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        legend.position = "none") +
  labs(y = "Schoener D", x = "Areas of Endemism", fill = "Regions") +
  annotation_custom(grob = grob_legend, 
                    xmin = -9, xmax = 10, 
                    ymin = 0.4, ymax = 1)

png("./figs/Boxplot_OverlapD.png", height = 800, width = 1000, res = 120)
overlap_anova_plot 
dev.off()

#saveRDS(overlap_anova_plot, "./figs/Boxplot_OverlapD.rds")
overlap_anova_plot <- readRDS("./figs/Boxplot_OverlapD.rds")

### 4.5. Auxiliary figure: categories of niche overlap ####

AoE_comp_df <- read.csv("./output/AoE/AoE_ECOSPAT_RES_areasonly.csv") %>%  filter(Analysis == "Eq", Alternative == "G")

AoE_comp_df$AoE_regs <- factor(AoE_comp_df$AoE_regs, levels = c("Mesoamerica", "Northern_Andes", "Guiana_centered", "Amazonia_centered", 
                                                                "ESA_Dry_Diagonal", "ESA_THDD", "ESA_Atlantic_Forest"), ordered = TRUE)
AoE_comp_df$D_classes <- factor(AoE_comp_df$D_classes, levels = c("Limited", "Low", "Moderate", "High", "Very_high"), ordered = TRUE)
AoE_comp_df$I_classes <- factor(AoE_comp_df$I_classes, levels = c("Limited", "Low", "Moderate", "High", "Very_high"), ordered = TRUE)
AoE_comp_df$Geo_classes <- factor(AoE_comp_df$Geo_classes, levels = c("Disjunct", "Low", "Medium", "High", "Very_high"), ordered = TRUE)
AoE_comp_df <- AoE_comp_df %>% arrange(AoE_regs)

x <- AoE_comp_df %>%
  group_by(AoE_regs, AoE_name) %>%
  summarise(aoe_medD = median(D)) %>%
  mutate(ordering = as.numeric(str_replace(AoE_name, pattern = "D2CA([0-9]+)", replacement = "\\1"))) %>%
  group_by(AoE_regs) %>%
  arrange(AoE_regs, ordering, aoe_medD) %>%
  pull(AoE_name)

AoE_comp_df$AoE_name <- factor(AoE_comp_df$AoE_name, levels = x, ordered = FALSE)

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
  ungroup() %>% 
  dplyr::select(1, 6)

overlap_df <- AoE_comp_df %>%
  dplyr::select(AoE_name, AoE_regs, D_classes) %>%
  left_join(area_code, by = c("AoE_name" = "aoe"))

overlap_df$areacode <- factor(overlap_df$areacode, levels = unique(overlap_df$areacode), ordered = TRUE)

categ_over <- overlap_df %>% 
  #filter(!is.na(niche_category2)) %>% 
  ggplot(aes(x = areacode, fill = D_classes)) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels=scales::percent) +
  scale_fill_grey() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  labs(fill = "Overlap",
       x = "Areas of Endemism",
       y = "Pair-wise Spp comparisons")

png("./figs/Categ_Overlap.png", height = 400, width = 1000, res = 120)
categ_over
dev.off()

#saveRDS(categ_over, "./figs/Auxiliary_Overlap.rds")
categ_over <- readRDS("./figs/Auxiliary_Overlap.rds")


### 4.6. Esospat all (166*166 spp) tests figure #####

AoE_comp_df <- read.csv("./output/AoE/AoE_ECOSPAT_RES_areasonly.csv")[, -1]

AoE_comp_df$D_classes <- factor(AoE_comp_df$D_classes, levels = c("Limited", "Low", "Moderate", "High", "Very_high"), ordered = TRUE)
AoE_comp_df$I_classes <- factor(AoE_comp_df$I_classes, levels = c("Limited", "Low", "Moderate", "High", "Very_high"), ordered = TRUE)
AoE_comp_df$AoE_name <- factor(AoE_comp_df$AoE_name, levels = unique(AoE_comp_df$AoE_name), ordered = FALSE)
AoE_comp_df$Geo_classes <- factor(AoE_comp_df$Geo_classes, levels = c("Disjunct", "Low", "Medium", "High", "Very_high"), ordered = TRUE)


glimpse(AoE_comp_df)
166*166

dp_eqlow_hist <- AoE_comp_df %>%
  filter(Analysis == "Eq", Alternative == "L") %>%
  ggplot(aes(x = D_p, fill = D_classes)) +
  geom_histogram() +
  scale_fill_grey() +
  theme_classic() +
  theme(legend.position = "none") +
  labs(x = paste0("Schoener D\nEquivalence Lower p-value\n", "\u03B1", "=0.05"),
       y = "Endemic species\npair comparisons")

dp_eqlow_col <- AoE_comp_df %>%
  filter(Analysis == "Eq", Alternative == "L") %>%
  ggplot(aes(x = Dp_class, fill = D_classes)) +
  geom_bar(width = 0.5) +
  scale_x_discrete(labels = c("Rejected", "Not Rejected")) +
  scale_fill_grey() +
  theme_classic() +
  theme(legend.position = "none") +
  labs(fill = "Niche overlap",
       x = NULL,
       y = NULL)

dp_eqgre_hist <- AoE_comp_df %>%
  filter(Analysis == "Eq", Alternative == "G") %>%
  ggplot(aes(x = D_p, fill = D_classes)) +
  geom_histogram() +
  theme_classic() +
  scale_fill_grey() +
  theme(legend.position = "none") +
  labs(x = paste0("Schoener D\nEquivalence Greater p-value\n", "\u03B1", "=0.05"),
       y = "Endemic species\npair comparisons")

dp_eqgre_col <- AoE_comp_df %>%
  filter(Analysis == "Eq", Alternative == "G") %>%
  ggplot(aes(x = Dp_class, fill = D_classes)) +
  geom_bar(width = 0.5) +
  scale_x_discrete(labels = c("Rejected", "Not Rejected")) +
  theme_classic() +
  scale_fill_grey() +
  theme(legend.position = "none") + 
  labs(fill = "Niche overlap",
       x = NULL,
       y = NULL)

dp_simlow_hist <- AoE_comp_df %>%
  filter(Analysis == "Sim", Alternative == "L") %>%
  ggplot(aes(x = D_p, fill = D_classes)) +
  geom_histogram() +
  theme_classic() +
  scale_fill_grey() +
  theme(legend.position = "none") +
  labs(x = paste0("Schoener D\nSimilarity Lower p-value\n", "\u03B1", "=0.05"),
       y = "Endemic species\npair comparisons")

dp_simlow_col <- AoE_comp_df %>%
  filter(Analysis == "Sim", Alternative == "L") %>%
  ggplot(aes(x = Dp_class, fill = D_classes)) +
  geom_bar(width = 0.5) +
  scale_x_discrete(labels = c("Not Rejected")) +
  theme_classic() +
  scale_fill_grey() +
  theme(legend.position = "none") +
  labs(fill = "Niche overlap",
       x = NULL,
       y = NULL)

dp_simgre_hist <- AoE_comp_df %>%
  filter(Analysis == "Sim", Alternative == "G") %>%
  ggplot(aes(x = D_p, fill = D_classes)) +
  geom_histogram() +
  theme_classic() +
  scale_fill_grey() +
  theme(legend.position = "none") +
  labs(x = paste0("Schoener D\nSimilarity Greater p-value\n", "\u03B1", "=0.05"),
       y = "Endemic species\npair comparisons")

dp_simgre_col <- AoE_comp_df %>%
  filter(Analysis == "Sim", Alternative == "G") %>%
  ggplot(aes(x = Dp_class, fill = D_classes)) +
  geom_bar(width = 0.5) +
  scale_x_discrete(labels = c("Rejected", "Not Rejected")) +
  theme_classic() +
  scale_fill_grey() +
  theme(legend.position = "none") + 
  theme_classic() +
  labs(fill = "Niche overlap", 
       x = NULL,
       y = NULL)

leg <- get_legend(dp_simgre_col)

dp_simgre_col <- dp_simgre_col + theme(legend.position = "none")

lay <- rbind(c(1, 1, 2, 3, 3, 4, NA),
             c(5, 5, 6, 7, 7, 8, 9))

ecospat_tests <- marrangeGrob(grobs = list(dp_eqlow_hist, dp_eqlow_col, dp_eqgre_hist, dp_eqgre_col,
                                           dp_simlow_hist, dp_simlow_col, dp_simgre_hist, dp_simgre_col, leg),
                              layout_matrix = lay, top = NULL)

# png("./figs/paper/ALLEND_ECOSPATTests.png", width = 1000, height = 750, res = 150)
# print(ecospat_tests)
# dev.off()

library(patchwork)

layout <- "AAB
CCD
EEF
GGH"
ecospat_tests <-  dp_eqgre_hist + dp_eqgre_col + dp_eqlow_hist + dp_eqlow_col +
  dp_simlow_hist + dp_simlow_col + dp_simgre_hist + dp_simgre_col +
  plot_layout(design = layout, guides = "collect")  +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold"))

png("./figs/paper/ALLEND_ECOSPATTests.png", width = 1000, height = 1200, res = 130)
print(ecospat_tests)
dev.off()


library(ecospat)
ecospat.niche.similarity.test()
### 4.7. Ecospat tests for end spp inside each AoE ####

{
  locs <-  read.csv("./output/localities_total_unique_ch10m_Bio(v14).csv", row.names = 1, stringsAsFactors = FALSE) # Localidades_v14
  AoE_comp_df <- read.csv("./output/AoE/AoE_ECOSPAT_RES_areasonly.csv") #%>%  filter(Analysis == "Eq", Alternative == "G")
  
  AoE_comp_df$AoE_regs <- factor(AoE_comp_df$AoE_regs, levels = c("Mesoamerica", "Northern_Andes", "Guiana_centered", "Amazonia_centered", 
                                                                  "ESA_Dry_Diagonal", "ESA_THDD", "ESA_Atlantic_Forest"), ordered = TRUE)
  AoE_comp_df$D_classes <- factor(AoE_comp_df$D_classes, levels = c("Limited", "Low", "Moderate", "High", "Very_high"), ordered = TRUE)
  AoE_comp_df$I_classes <- factor(AoE_comp_df$I_classes, levels = c("Limited", "Low", "Moderate", "High", "Very_high"), ordered = TRUE)
  AoE_comp_df$Geo_classes <- factor(AoE_comp_df$Geo_classes, levels = c("Disjunct", "Low", "Medium", "High", "Very_high"), ordered = TRUE)
  AoE_comp_df <- AoE_comp_df %>% arrange(AoE_regs)
  
  locs_end <- locs %>%
    filter(NAME1 %in% unique(ca_db$Species))
  spp_names <- unique(locs_end$NAME1)
  
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
  
  x <- AoE_comp_df %>%
    group_by(AoE_regs, AoE_name) %>%
    summarise(aoe_medD = median(D)) %>%
    mutate(ordering = as.numeric(str_replace(AoE_name, pattern = "D2CA([0-9]+)", replacement = "\\1"))) %>%
    group_by(AoE_regs) %>%
    arrange(AoE_regs, ordering, aoe_medD) %>%
    pull(AoE_name)
  
  AoE_comp_df$AoE_name <- factor(AoE_comp_df$AoE_name, levels = x, ordered = FALSE)
  
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
    ungroup() %>% 
    dplyr::select(1, 6)
  
  overlap_df <- AoE_comp_df %>%
    left_join(area_code, by = c("AoE_name" = "aoe")) %>% 
    mutate(HO_rejection = case_when(
      Dp_class == "P<0.05" ~ "Rejected",
      Dp_class == "P>=0.05" ~ "Not Rejected"
    ),
    test_state = paste(Analysis, Alternative, HO_rejection, sep = "_"))
  
  overlap_df$areacode <- factor(overlap_df$areacode, levels = unique(overlap_df$areacode), ordered = TRUE)
}

table(overlap_df$areacode, overlap_df$test_state)


plot_testaux <- function(df, analysis = "Eq", alternative = "L", fill_t = "H0: Equivalence Lower") {
  
  df %>% 
    filter(Analysis == analysis, Alternative == alternative) %>%
    ggplot(aes(x = areacode, fill = HO_rejection)) +
    geom_bar(position = "fill") +
    scale_y_continuous(labels=scales::percent) +
    scale_fill_grey(start = 0.4) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    labs(fill = fill_t,
         x = "Areas of Endemism",
         y = "Pair-wise Spp comparisons")
  
  
}

eql_plot <- plot_testaux(df = overlap_df, analysis = "Eq", alternative = "L", fill_t = "H0: Equivalence Lower")
eqg_plot <- plot_testaux(df = overlap_df, analysis = "Eq", alternative = "G", fill_t = "H0: Equivalence Greater")
siml_plot <- plot_testaux(df = overlap_df, analysis = "Sim", alternative = "L", fill_t = "H0: Similarity Lower")
simg_plot <- plot_testaux(df = overlap_df, analysis = "Sim", alternative = "G", fill_t = "H0: Similarity Greater")

png("./figs/Auxiliary_EqSimtests.png", width = 1500, height = 1800, res = 130)
print((eqg_plot / eql_plot / simg_plot / siml_plot) +
        plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face = "bold")))
dev.off()

overlap_df %>% 
  filter(Analysis == "Eq", Alternative == "L") %>%
  # summarise(largo = n()) %>% 
  distinct(Analysis, Alternative, areacode, Sp1name, .keep_all = TRUE) %>%
  summarise(largo = n()) %>% 
  ggplot(aes(x = areacode, fill = HO_rejection)) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels=scales::percent) +
  scale_fill_grey(start = 0.4) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  labs(fill = "H0: Equivalence Lower",
       x = "Areas of Endemism",
       y = "Pair-wise Spp comparisons")


### 4.8. Ecospat tests: review of species ####

{
  locs <-  read.csv("./output/localities_total_unique_ch10m_Bio(v14).csv", row.names = 1, stringsAsFactors = FALSE) # Localidades_v14
  AoE_comp_df <- read.csv("./output/AoE/AoE_ECOSPAT_RES_areasonly.csv") #%>%  filter(Analysis == "Eq", Alternative == "G")
  
  AoE_comp_df$AoE_regs <- factor(AoE_comp_df$AoE_regs, levels = c("Mesoamerica", "Northern_Andes", "Guiana_centered", "Amazonia_centered", 
                                                                  "ESA_Dry_Diagonal", "ESA_THDD", "ESA_Atlantic_Forest"), ordered = TRUE)
  AoE_comp_df$D_classes <- factor(AoE_comp_df$D_classes, levels = c("Limited", "Low", "Moderate", "High", "Very_high"), ordered = TRUE)
  AoE_comp_df$I_classes <- factor(AoE_comp_df$I_classes, levels = c("Limited", "Low", "Moderate", "High", "Very_high"), ordered = TRUE)
  AoE_comp_df$Geo_classes <- factor(AoE_comp_df$Geo_classes, levels = c("Disjunct", "Low", "Medium", "High", "Very_high"), ordered = TRUE)
  AoE_comp_df <- AoE_comp_df %>% arrange(AoE_regs)
  
  locs_end <- locs %>%
    filter(NAME1 %in% unique(ca_db$Species))
  spp_names <- unique(locs_end$NAME1)
  
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
  
  x <- AoE_comp_df %>%
    group_by(AoE_regs, AoE_name) %>%
    summarise(aoe_medD = median(D)) %>%
    mutate(ordering = as.numeric(str_replace(AoE_name, pattern = "D2CA([0-9]+)", replacement = "\\1"))) %>%
    group_by(AoE_regs) %>%
    arrange(AoE_regs, ordering, aoe_medD) %>%
    pull(AoE_name)
  
  AoE_comp_df$AoE_name <- factor(AoE_comp_df$AoE_name, levels = x, ordered = FALSE)
  
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
    ungroup() %>% 
    dplyr::select(1, 6)
  
  overlap_df <- AoE_comp_df %>%
    left_join(area_code, by = c("AoE_name" = "aoe")) %>% 
    mutate(HO_rejection = case_when(
      Dp_class == "P<0.05" ~ "Rejected",
      Dp_class == "P>=0.05" ~ "Not Rejected"
    ),
    test_state = paste(Analysis, Alternative, HO_rejection, sep = "_"))
  
  overlap_df$areacode <- factor(overlap_df$areacode, levels = unique(overlap_df$areacode), ordered = TRUE)
}

overlap_df <- overlap_df %>% 
  mutate(pairspp = paste(Sp.1, Sp.2, sep = "_"))

eq_reject_spp <- overlap_df %>%
  filter(HO_rejection == "Rejected", Analysis == "Eq", Alternative == "L") %>%
  pull(pairspp)

sim_reject_spp <- overlap_df %>%
  filter(HO_rejection == "Rejected", Analysis == "Sim") %>%
  pull(pairspp)

sum(sim_reject_spp %in% eq_reject_spp)
index_rr <- which(eq_reject_spp %in% sim_reject_spp)

rr_df <- overlap_df %>%
  filter(HO_rejection == "Rejected", Analysis == "Eq", Alternative == "L")

rr_df <- rr_df[index_rr, ]

unique(rr_df$areacode)
# [1] Amazonia_3   Amazonia_6   DD_AF_1      DD_AF_5      Atl_Forest_1

# overlap_df %>%
#   filter(HO_rejection == "Rejected", areacode == "Amazonia_3")
# 
# overlap_df %>%
#   filter(Analysis == "Sim", areacode == "Amazonia_3")

