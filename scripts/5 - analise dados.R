################################################################################
################## Comparison between metrics #################################
################################################################################

# Loading tables of calculated components
dados_comp_esp <- read.delim("dissi_comparar_esp.csv", sep = ";", header = TRUE)
dados_comp_tem <- read.delim("dissi_comparar_tem.csv", sep = ";", header = TRUE)

# Calculate differences for temporal data
dif_ST_tem <- dados_comp_tem$ST_Pois - dados_comp_tem$ST_Com
dif_OS_tem <- dados_comp_tem$OS_Pois - dados_comp_tem$OS_Com
dif_WN_tem <- dados_comp_tem$WN_Pois - dados_comp_tem$WN_Com

# Normality test (Shapiro-Wilk) for differences
cat("=== NORMALITY OF DIFFERENCES (TEMPORAL) ===\n")
cat("\nST (Species Turnover):\n")
shapiro.test(dif_ST_tem)

cat("\nOS (Rewiring):\n")
shapiro.test(dif_OS_tem)

cat("\nWN (Total dissimilarity):\n")
shapiro.test(dif_WN_tem)

# Visualization (optional, but helpful)
par(mfrow = c(1, 3))
qqnorm(dif_ST_tem, main = "ST - Differences")
qqline(dif_ST_tem, col = "red")
qqnorm(dif_OS_tem, main = "OS - Differences")
qqline(dif_OS_tem, col = "red")
qqnorm(dif_WN_tem, main = "WN - Differences")
qqline(dif_WN_tem, col = "red")

# Comparing temporal metrics using paired Wilcoxon test
wilcox.test(dados_comp_tem$ST_Pois, dados_comp_tem$ST_Com, paired = TRUE)
t.test(dados_comp_tem$ST_Pois, dados_comp_tem$ST_Com, paired = TRUE)

wilcox.test(dados_comp_tem$OS_Pois, dados_comp_tem$OS_Com, paired = TRUE)
t.test(dados_comp_tem$OS_Pois, dados_comp_tem$OS_Com, paired = TRUE)

wilcox.test(dados_comp_tem$WN_Pois, dados_comp_tem$WN_Com, paired = TRUE)
t.test(dados_comp_tem$WN_Pois, dados_comp_tem$WN_Com, paired = TRUE)

# Organizing data into long format and resolving duplicates
dados_long <- dados_comp_tem %>%
  select(ST_Pois, ST_Com, OS_Pois, OS_Com, WN_Pois, WN_Com) %>%
  pivot_longer(
    cols = everything(),
    names_to = c("Component", "Metric"),
    names_sep = "_",
    values_to = "Value"
  ) %>%
  group_by(Component, Metric) %>%
  mutate(ID = row_number()) %>%  # Creates a unique ID per row
  pivot_wider(
    names_from = Metric,
    values_from = Value
  ) %>%
  ungroup()

# Checking if it worked
glimpse(dados_long)

# Plotting
scatter_plot <- ggplot(dados_long, aes(x = Pois, y = Com)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", color = "blue", se = FALSE) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  facet_wrap(~Component, scales = "free") +
  theme_classic() +
  labs(
    x = "Poisot Metric",
    y = "Commondenom Metric",
    title = "A"
  )

# Data in long format for plotting
box_data <- dados_comp_tem %>%
  select(ST_Pois, ST_Com, OS_Pois, OS_Com, WN_Pois, WN_Com) %>%
  pivot_longer(
    cols = everything(),
    names_to = c("Component", "Metric"),
    names_sep = "_",
    values_to = "Value"
  )

# Plotting boxplot
box_plot <- ggplot(box_data, aes(x = Metric, y = Value, fill = Metric)) +
  geom_boxplot(alpha = 0.7) +
  stat_summary(fun = mean, geom = "point", shape = 21, size = 3, 
               fill = "black", color = "white") +
  facet_wrap(~Component, scales = "free") +
  theme_classic() +
  scale_fill_manual(values = c("gray70", "gray30")) +
  labs(
    x = "Metric",
    y = "Dissimilarity Value",
    title = "B"
  ) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold", size = 12),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  )

# Combine both with patchwork
scatter_plot / box_plot

ggsave("figure_dissimilarity_comparison.png", scatter_plot / box_plot,
       width = 10, height = 12, dpi = 300)

#### Now doing spatial comparisons
# Comparing spatial metrics using paired Wilcoxon test
wilcox.test(dados_comp_esp$ST_Poi, dados_comp_esp$ST_Com, paired = TRUE)
t.test(dados_comp_esp$ST_Poi, dados_comp_esp$ST_Com, paired = TRUE)

wilcox.test(dados_comp_esp$OS_Poi, dados_comp_esp$OS_Com, paired = TRUE)
t.test(dados_comp_esp$OS_Poi, dados_comp_esp$OS_Com, paired = TRUE)

wilcox.test(dados_comp_esp$WN_Poi, dados_comp_esp$WN_Com, paired = TRUE)
t.test(dados_comp_esp$WN_Poi, dados_comp_esp$WN_Com, paired = TRUE)

# Organizing data into long format and resolving duplicates
dados_long_esp <- dados_comp_esp %>%
  select(ST_Poi, ST_Com, OS_Poi, OS_Com, WN_Poi, WN_Com) %>%
  pivot_longer(
    cols = everything(),
    names_to = c("Component", "Metric"),
    names_sep = "_",
    values_to = "Value"
  ) %>%
  group_by(Component, Metric) %>%
  mutate(ID = row_number()) %>%
  pivot_wider(
    names_from = Metric,
    values_from = Value
  ) %>%
  ungroup()

# Checking if it worked
glimpse(dados_long_esp)

# Plotting
scatter_plot_esp <- ggplot(dados_long_esp, aes(x = Poi, y = Com)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", color = "blue", se = FALSE) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  facet_wrap(~Component, scales = "free") +
  theme_classic() +
  labs(
    x = "Poisot Metric",
    y = "Commondenom Metric",
    title = "A"
  )

# Data in long format for plotting
box_data_esp <- dados_comp_esp %>%
  select(ST_Poi, ST_Com, OS_Poi, OS_Com, WN_Poi, WN_Com) %>%
  pivot_longer(
    cols = everything(),
    names_to = c("Component", "Metric"),
    names_sep = "_",
    values_to = "Value"
  )

# Plotting boxplot
box_plot_esp <- ggplot(box_data_esp, aes(x = Metric, y = Value, fill = Metric)) +
  geom_boxplot(alpha = 0.7) +
  stat_summary(fun = mean, geom = "point", shape = 21, size = 3, 
               fill = "black", color = "white") +
  facet_wrap(~Component, scales = "free") +
  theme_classic() +
  scale_fill_manual(values = c("gray70", "gray30")) +
  labs(
    x = "Metric",
    y = "Dissimilarity Value",
    title = "B"
  ) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold", size = 12),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  )

scatter_plot_esp / box_plot_esp

ggsave("figure_dissimilarity_spatial.png", scatter_plot_esp / box_plot_esp,
       width = 10, height = 12, dpi = 300)

### MRM
# Loading data to calculate heterogeneity 
dados_heter <- read.delim("hetero.csv", sep = ";", header = TRUE)
dados_para_distancia <- dados_heter[, c("simpson", "patch.richness")]

# Calculate distance matrix
matriz_distancia <- dist(dados_para_distancia, method = "euclidean")
print(matriz_distancia)

# Convert 'dist' object to a full matrix
matriz_heterogeneidade_completa <- as.matrix(matriz_distancia)
df_heterogeneidade_para_csv <- as.data.frame(matriz_heterogeneidade_completa)

# Save the dataframe to a CSV file
write.csv(df_heterogeneidade_para_csv, "distance_matrix_heterogeneity.csv", 
          row.names = TRUE)

# Now load distance data
xy <- read.table("xy.csv", header = TRUE, sep = ";", row.names = 1)
distgeo.s <- dist(xy, method = "euclidean")
matriz_geografica_completa <- as.matrix(distgeo.s)
df_geografica_para_csv <- as.data.frame(matriz_geografica_completa)
write.csv(df_geografica_para_csv, "distance_matrix_geographic.csv", 
          row.names = TRUE)

# Transforming variables
dist_env <- list(matriz_distancia, distgeo.s)
sapply(dist_env, attr, "Size")
attr(matrizos, "Size")

class(matrizos)  # Should be "dist"
lapply(dist_env, class)  # All should be "dist"

##### Now everything together
# Beta Poisot
regres_beta_os <- ecodist::MRM(matrizos ~ distgeo.s + matriz_distancia, 
                               nperm = 1000, method = "linear")
regres_beta_os

regres_beta_st <- ecodist::MRM(matrizst ~ distgeo.s + matriz_distancia, 
                               nperm = 2000, method = "linear")
regres_beta_st

regres_beta_wt <- ecodist::MRM(matrizwt ~ distgeo.s + matriz_distancia, 
                               nperm = 1000, method = "linear")
regres_beta_wt

# Beta Common
regres_beta_os_com <- ecodist::MRM(matrizos_com ~ distgeo.s + matriz_distancia, 
                                   nperm = 1000, method = "linear")
regres_beta_os_com

regres_beta_st_com <- ecodist::MRM(matrizst_com ~ distgeo.s + matriz_distancia, 
                                   nperm = 1000, method = "linear")
regres_beta_st_com

regres_beta_wt_com <- ecodist::MRM(matrizwt_com ~ distgeo.s + matriz_distancia, 
                                   nperm = 1000, method = "linear")
regres_beta_wt_com

## Subcomponents
regres_beta_st_L <- ecodist::MRM(matrizST.l_com ~ distgeo.s + matriz_distancia, 
                                 nperm = 1000, method = "linear")
regres_beta_st_L

regres_beta_st_h <- ecodist::MRM(matrizST.h_com ~ distgeo.s + matriz_distancia, 
                                 nperm = 1000, method = "linear")
regres_beta_st_h

regres_beta_st_lh <- ecodist::MRM(matrizST.lh_com ~ distgeo.s + matriz_distancia, 
                                  nperm = 1000, method = "linear")
regres_beta_st_lh

regres_beta_wn_rep <- ecodist::MRM(matrizwn.rep_com ~ distgeo.s + matriz_distancia, 
                                   nperm = 1000, method = "linear")
regres_beta_wn_rep

regres_beta_wn_ric <- ecodist::MRM(matrizwn.ric_com ~ distgeo.s + matriz_distancia, 
                                   nperm = 1000, method = "linear")
regres_beta_wn_ric

### Plotting

# First Poisot
# Loading data
dados_plot_poisot <- read.csv("plot_relacao_line_poisot.csv", 
                              header = TRUE, sep = ";")

# Preparing data for plotting
dados_plot_poisot$dissimilarity <- dados_plot_poisot$bt
dados_plot_poisot$Heterogeneity_2 <- dados_plot_poisot$hetero
dados_plot_poisot$Beta_pois <- as.factor(dados_plot_poisot$Beta_pois) 

# Generate Plots - Distance vs. Dissimilarity
plot1_dist_poi <- ggplot(dados_plot_poisot, aes(x = Dist_km, y = dissimilarity,
                                                color = Beta_pois)) +
  geom_smooth(method = "lm", se = TRUE) +
  theme_minimal() +
  labs(title = NULL, x = "Distance (km)", y = "Dissimilarity") +
  theme(legend.position = "none")

# Heterogeneity vs. Dissimilarity
plot1_het_poi <- ggplot(dados_plot_poisot, aes(x = Heterogeneity_2, y = dissimilarity,
                                               color = Beta_pois)) +
  geom_smooth(method = "lm", se = TRUE) +
  theme_minimal() +
  labs(title = NULL, x = "Heterogeneity (Euclidean)", y = "") +
  theme(legend.position = "right")

# Combine and save
print(plot1_dist_poi + plot1_het_poi)
ggsave("figure_distance_heterogeneity.png", plot1_dist_poi + plot1_het_poi,
       width = 10, height = 12, dpi = 300)

# Now Commondenom
dados_plot_original <- read.csv("plot_relacao_line_common.csv", 
                                header = TRUE, sep = ";")

# Preparing data for plotting
dados_plot_original$dissimilarity <- dados_plot_original$bt
dados_plot_original$Heterogeneity_2 <- dados_plot_original$hetero
dados_plot_original$Beta_com <- as.factor(dados_plot_original$Beta_com) 

# Define categories
main_categories <- c("WN", "OS", "ST")
subcomponent_categories <- c("ST_L", "ST_H", "ST_LH", "WN_repl", "WN_rich")

# Filter dataframes
dados_plot_principais <- dados_plot_original %>%
  filter(Beta %in% main_categories)

dados_plot_subcomponentes <- dados_plot_original %>%
  filter(Beta %in% subcomponent_categories)

# GROUP 1: WN, OS, ST
plot1_dist <- ggplot(dados_plot_principais, aes(x = Dist_km, y = dissimilarity,
                                                color = Beta)) +
  geom_smooth(method = "lm", se = TRUE) +
  theme_minimal() +
  labs(title = NULL, x = "Distance (km)", y = "Dissimilarity") +
  theme(legend.position = "none")

plot1_het <- ggplot(dados_plot_principais, aes(x = Heterogeneity_2, y = dissimilarity,
                                               color = Beta)) +
  geom_smooth(method = "lm", se = TRUE) +
  theme_minimal() +
  labs(title = NULL, x = "Heterogeneity (Euclidean)", y = "") +
  theme(legend.position = "right")

# GROUP 2: Subcomponents
plot2_dist <- ggplot(dados_plot_subcomponentes, aes(x = Dist_km, y = dissimilarity,
                                                    color = Beta)) +
  geom_smooth(method = "lm", se = TRUE) +
  theme_minimal() +
  labs(title = NULL, x = "Distance (km)", y = "Dissimilarity") +
  theme(legend.position = "none")

plot2_het <- ggplot(dados_plot_subcomponentes, aes(x = Heterogeneity_2, y = dissimilarity,
                                                   color = Beta)) +
  geom_smooth(method = "lm", se = TRUE) +
  theme_minimal() +
  labs(title = NULL, x = "Heterogeneity (Euclidean)", y = "") +
  theme(legend.position = "right")

# Display the Plots
cat("### Group of Plots: WN, OS, ST ###\n")
print(plot1_dist + plot1_het)

cat("\n### Group of Plots: ST_L, ST_H, ST_LH, WN_repl, WN_rich ###\n")
print(plot2_dist + plot2_het)

# Save images
ggsave("figure_distance_heterogeneity_group1.png", plot1_dist + plot1_het,
       width = 10, height = 12, dpi = 300)

ggsave("figure_distance_heterogeneity_group2.png", plot2_dist + plot2_het,
       width = 10, height = 12, dpi = 300)