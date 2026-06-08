################################################################################
################# Dissimilarity Poisot et al. 2012 #############################
################################################################################

# Calculating dissimilarity between networks by dry and rainy seasons
# First step: calculate network dissimilarity (Poisot et al. 2012)
rede1 <- betalinkr(webs2array(dados1a, dados1b), 
                   index = "jaccard", function.dist = "vegdist", partitioning = "poisot",
                   partition.st = T, partition.rr = F)
rede1

rede2 <- betalinkr(webs2array(dados2a, dados2b), 
                   index = "jaccard", function.dist = "vegdist", partitioning = "poisot")
rede2

rede3 <- betalinkr(webs2array(dados3a, dados3b), 
                   index = "jaccard", function.dist = "vegdist", partitioning = "poisot")
rede3

rede4 <- betalinkr(webs2array(dados4a, dados4b), 
                   index = "jaccard", function.dist = "vegdist", partitioning = "poisot")
rede4

rede5 <- betalinkr(webs2array(dados5a, dados5b), 
                   index = "jaccard", function.dist = "vegdist", partitioning = "poisot")
rede5

rede6 <- betalinkr(webs2array(dados6a, dados6b), 
                   index = "jaccard", function.dist = "vegdist", partitioning = "poisot")
rede6

rede7 <- betalinkr(webs2array(dados7a, dados7b), 
                   index = "jaccard", function.dist = "vegdist", partitioning = "poisot")
rede7

rede8 <- betalinkr(webs2array(dados8a, dados8b), 
                   index = "jaccard", function.dist = "vegdist", partitioning = "poisot")
rede8

rede9 <- betalinkr(webs2array(dados9a, dados9b), 
                   index = "jaccard", function.dist = "vegdist", partitioning = "poisot")
rede9

rede10 <- betalinkr(webs2array(dados10a, dados10b), 
                    index = "jaccard", function.dist = "vegdist", partitioning = "poisot")
rede10

rede11 <- betalinkr(webs2array(dados11a, dados11b), 
                    index = "jaccard", function.dist = "vegdist", partitioning = "poisot")
rede11

rede12 <- betalinkr(webs2array(dados12a, dados12b), 
                    index = "jaccard", function.dist = "vegdist", binary = T, partitioning = "poisot")
rede12

# Calculate means and standard deviations
season <- read.delim("dissi_estac.csv", sep = ";", header = TRUE, row.names = 1) 
str(season)
mean(season$OS); var(season$OS); sd(season$OS)
mean(season$WN); var(season$WN); sd(season$WN)
mean(season$ST); var(season$ST); sd(season$ST)

# Calculating dissimilarity between networks using betalinkr_multi function
# which calculates dissimilarity for multiple networks
tab <- betalinkr_multi(webs2array(dados1, dados2, dados3, dados4, dados5, dados6,
                                  dados7, dados8, dados9, dados10, dados11, dados12),
                       index = "jaccard", function.dist = "vegdist", binary = T,
                       partitioning = "poisot")
tab

# Export table to CSV
write.csv2(tab, "dissi2.csv", row.names = FALSE) # Jaccard

# Where:
# ST = Changes in species composition (species turnover)
# OS = Spatial rewiring of interactions
# WN = Interaction beta-diversity (total network dissimilarity)
# S = Species composition dissimilarity

# Now load dissimilarity matrices
matrizos <- read.delim("matriz_os.csv", sep = ";", header = TRUE, row.names = 1) 
matrizos <- as.dist(matrizos, diag = FALSE, upper = FALSE)
matrizos

matrizst <- read.delim("matriz_st.csv", sep = ";", header = TRUE, row.names = 1) 
matrizst <- as.dist(matrizst, diag = FALSE, upper = FALSE)
matrizst

matrizwt <- read.delim("matriz_wt.csv", sep = ";", header = TRUE, row.names = 1) 
matrizwt <- as.dist(matrizwt, diag = FALSE, upper = FALSE)
matrizwt

# Running dbRDA
# 500 m buffer
cobertura_500a <- decostand(x = cobertura_500a, method = "standardize") 
dbrda.os_500 <- dbrda(matrizos ~ AREA + NP + ED + ENM, 
                      data = cobertura_500a)
dbrda.os_500
plot(dbrda.os_500)
summary(dbrda.os_500)
anova.cca(dbrda.os_500, by = "terms", permutations = how(nperm = 9999))
RsquareAdj(dbrda.os_500)
vif.cca(dbrda.os_500)

# Plotting
p1 <- ggord.dbrda(dbrda.os_500, ptslab = TRUE, size = 1, addsize = 3, repel = TRUE) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  theme(legend.position = 'right') +
  labs(title = expression(beta[OS])) 
p1

##
dbrda.st_500 <- dbrda(matrizst ~ AREA + NP + ED + ENM, data = cobertura_500a)
anova.cca(dbrda.st_500, by = "terms", permutations = how(nperm = 9999))
RsquareAdj(dbrda.st_500)
summary(dbrda.st_500)

# Plotting
p3 <- ggord.dbrda(dbrda.st_500, ptslab = TRUE, size = 1, addsize = 3, repel = TRUE, title = "St") +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  theme(legend.position = 'right') +
  labs(title = expression(beta[ST]))  
p3

##
dbrda.wt_500 <- dbrda(matrizwt ~ AREA + NP + ED + ENM, data = cobertura_500a)
summary(dbrda.wt_500)
anova.cca(dbrda.wt_500, by = "terms", permutations = how(nperm = 9999))
RsquareAdj(dbrda.wt_500)

# Plotting
p4 <- ggord.dbrda(dbrda.wt_500, ptslab = TRUE, size = 1, addsize = 3, repel = TRUE) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  theme(legend.position = 'right') +
  labs(title = expression(beta[WN])) 
p4

# Combine plots
combined_plot <- (p1 + p3 + p4) +
  plot_annotation(tag_levels = "a") +
  plot_layout(nrow = 1)

# Save
ggsave(
  filename = "Figure_3.jpg",
  plot = combined_plot,
  width = 18,
  height = 10,
  units = "cm",
  dpi = 300,
  device = "jpeg"
)

# 1500 m buffer
##
dbrda.os_1500 <- dbrda(matrizos ~ AREA_1500 + NP_1500 + ED_1500 + ENM_1500, 
                       data = cobertura_1500)

anova.cca(dbrda.os_1500, by = "terms", permutations = how(nperm = 9999))
RsquareAdj(dbrda.os_1500)
summary(dbrda.os_1500)
vif.cca(dbrda.os_1500)

# Plotting
p5 <- ggord.dbrda(dbrda.os_1500, ptslab = TRUE, size = 1, addsize = 3, repel = TRUE) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  theme(legend.position = 'right') +
  labs(title = expression(beta[OS])) 
p5

##
dbrda.st_1500 <- dbrda(matrizst ~ AREA_1500 + NP_1500 + ED_1500 + ENM_1500, 
                       data = cobertura_1500)
anova.cca(dbrda.st_1500, by = "terms", permutations = how(nperm = 9999))
RsquareAdj(dbrda.st_1500)
summary(dbrda.st_1500)

# Plotting
p7 <- ggord.dbrda(dbrda.st_1500, ptslab = TRUE, size = 1, addsize = 3, repel = TRUE) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  theme(legend.position = 'right') +
  labs(title = expression(beta[ST])) 
p7

##
dbrda.wt_1500 <- dbrda(matrizwt ~ AREA_1500 + NP_1500 + ED_1500 + ENM_1500, 
                       data = cobertura_1500)
anova.cca(dbrda.wt_1500, by = "terms", permutations = how(nperm = 9999))
RsquareAdj(dbrda.wt_1500)
summary(dbrda.wt_1500)

# Plotting
p8 <- ggord.dbrda(dbrda.wt_1500, ptslab = TRUE, size = 1, addsize = 3, repel = TRUE) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  theme(legend.position = 'right') +
  labs(title = expression(beta[WN])) 
p8

# Combine plots
combined_plot2 <- (p5 + p7 + p8) +
  plot_annotation(tag_levels = "a") +
  plot_layout(nrow = 1)

# Save
ggsave(filename = "Figure_S3.jpg", plot = combined_plot2, 
       width = 18, height = 10, units = "cm", dpi = 300, device = "jpeg")

# 3000 m buffer
cobertura_3000 <- decostand(x = cobertura_3000, method = "standardize")

dbrda.os_3000 <- dbrda(matrizos ~ AREA_3000 + NP_3000 + ED_3000 + ENM_3000, 
                       data = cobertura_3000)
anova.cca(dbrda.os_3000, by = "terms", permutations = how(nperm = 9999))
RsquareAdj(dbrda.os_3000)
vif.cca(dbrda.os_3000)
summary(dbrda.os_3000)

# Plotting
p9 <- ggord.dbrda(dbrda.os_3000, ptslab = TRUE, size = 1, addsize = 3, repel = TRUE) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  theme(legend.position = 'right') +
  labs(title = expression(beta[OS])) 
p9

##
dbrda.st_3000 <- dbrda(matrizst ~ AREA_3000 + NP_3000 + ED_3000 + ENM_3000, 
                       data = cobertura_3000)
anova.cca(dbrda.st_3000, by = "terms", permutations = how(nperm = 9999))
RsquareAdj(dbrda.st_3000)
summary(dbrda.st_3000)

# Plotting
p11 <- ggord.dbrda(dbrda.st_3000, ptslab = TRUE, size = 1, addsize = 3, repel = TRUE) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  theme(legend.position = 'right') +
  labs(title = expression(beta[ST])) 
p11

##
dbrda.wt_3000 <- dbrda(matrizwt ~ AREA_3000 + NP_3000 + ED_3000 + ENM_3000, 
                       data = cobertura_3000)
anova.cca(dbrda.wt_3000, by = "terms", permutations = how(nperm = 9999))
RsquareAdj(dbrda.wt_3000)
summary(dbrda.wt_3000)

# Plotting
p12 <- ggord.dbrda(dbrda.wt_3000, ptslab = TRUE, size = 1, addsize = 3, repel = TRUE) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  theme(legend.position = 'right') +
  labs(title = expression(beta[WN])) 
p12

# Combine plots
combined_plot3 <- (p9 + p11 + p12) +
  plot_annotation(tag_levels = "a") +
  plot_layout(nrow = 1)

# Save
ggsave(filename = "Figure_S4.jpg", plot = combined_plot3, 
       width = 18, height = 10, units = "cm", dpi = 300, device = "jpeg")
