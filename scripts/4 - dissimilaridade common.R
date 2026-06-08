################################################################################
################# Dissimilarity Commondenom ####################################
################################################################################

# Calculate Commondenom
redea <- betalinkr(webs2array(dados1a, dados1b), 
                   index = "jaccard", function.dist = "vegdist", partitioning = "commondenom",
                   partition.st = TRUE, partition.rr = TRUE)
redea

redeb <- betalinkr(webs2array(dados2a, dados2b), 
                   index = "jaccard", function.dist = "vegdist", partitioning = "commondenom",
                   partition.st = TRUE, partition.rr = TRUE)
redeb

redec <- betalinkr(webs2array(dados3a, dados3b), 
                   index = "jaccard", function.dist = "vegdist", partitioning = "commondenom",
                   partition.st = TRUE, partition.rr = TRUE)
redec

reded <- betalinkr(webs2array(dados4a, dados4b), 
                   index = "jaccard", function.dist = "vegdist", partitioning = "commondenom",
                   partition.st = TRUE, partition.rr = TRUE)
reded

redee <- betalinkr(webs2array(dados5a, dados5b), 
                   index = "jaccard", function.dist = "vegdist", partitioning = "commondenom",
                   partition.st = TRUE, partition.rr = TRUE)
redee

redef <- betalinkr(webs2array(dados6a, dados6b), 
                   index = "jaccard", function.dist = "vegdist", partitioning = "commondenom",
                   partition.st = TRUE, partition.rr = TRUE)
redef

redeg <- betalinkr(webs2array(dados7a, dados7b), 
                   index = "jaccard", function.dist = "vegdist", partitioning = "commondenom",
                   partition.st = TRUE, partition.rr = TRUE)
redeg

redeh <- betalinkr(webs2array(dados8a, dados8b), 
                   index = "jaccard", function.dist = "vegdist", partitioning = "commondenom",
                   partition.st = TRUE, partition.rr = TRUE)
redeh

redei <- betalinkr(webs2array(dados9a, dados9b), 
                   index = "jaccard", function.dist = "vegdist", partitioning = "commondenom",
                   partition.st = TRUE, partition.rr = TRUE)
redei

redej <- betalinkr(webs2array(dados10a, dados10b), 
                   index = "jaccard", function.dist = "vegdist", partitioning = "commondenom",
                   partition.st = TRUE, partition.rr = TRUE)
redej

redel <- betalinkr(webs2array(dados11a, dados11b), 
                   index = "jaccard", function.dist = "vegdist", partitioning = "commondenom",
                   partition.st = TRUE, partition.rr = TRUE)
redel

redem <- betalinkr(webs2array(dados12a, dados12b), 
                   index = "jaccard", function.dist = "vegdist", partitioning = "commondenom",
                   partition.st = TRUE, partition.rr = TRUE)
redem

season_common <- read.delim("dissi_estac_common.csv", sep = ";", header = TRUE) 
str(season_common)

mean(season_common$OS); var(season_common$OS); sd(season_common$OS)
mean(season_common$WN); var(season_common$WN); sd(season_common$WN)
mean(season_common$ST); var(season_common$ST); sd(season_common$ST)
mean(season_common$ST_L); var(season_common$ST_L); sd(season_common$ST_L)
mean(season_common$ST_h); var(season_common$ST_h); sd(season_common$ST_h)
mean(season_common$ST_LH); var(season_common$ST_LH); sd(season_common$ST_LH)
mean(season_common$WN.repl); var(season_common$WN.repl); sd(season_common$WN.repl)
mean(season_common$WN.rich); var(season_common$WN.rich); sd(season_common$WN.rich)

# Calculate dissimilarity between networks using betalinkr_multi function
# which calculates dissimilarity for multiple networks
tab_common <- betalinkr_multi(webs2array(dados1, dados2, dados3, dados4, dados5, dados6,
                                         dados7, dados8, dados9, dados10, dados11, dados12),
                              index = "jaccard", function.dist = "vegdist", binary = T,
                              partitioning = "commondenom", partition.st = TRUE, partition.rr = TRUE)
tab_common

# Export table to CSV
write.csv2(tab_common, "dissi_common2.csv", row.names = FALSE)

# Now load dissimilarity matrices
#
matrizos_com <- read.delim("matriz_os_com.csv", sep = ";", header = TRUE, row.names = 1) 
matrizos_com <- as.dist(matrizos_com, diag = FALSE, upper = FALSE)
matrizos_com

#
matrizst_com <- read.delim("matriz_st_com.csv", sep = ";", header = TRUE, row.names = 1) 
matrizst_com <- as.dist(matrizst_com, diag = FALSE, upper = FALSE)
matrizst_com

#
matrizwt_com <- read.delim("matriz_wt_com.csv", sep = ";", header = TRUE, row.names = 1) 
matrizwt_com <- as.dist(matrizwt_com, diag = FALSE, upper = FALSE)
matrizwt_com

#
matrizST.l_com <- read.delim("matriz_ST.l_com.csv", sep = ";", header = TRUE, row.names = 1) 
matrizST.l_com <- as.dist(matrizST.l_com, diag = FALSE, upper = FALSE)
matrizST.l_com

#
matrizST.h_com <- read.delim("matriz_ST.H_com.csv", sep = ";", header = TRUE, row.names = 1) 
matrizST.h_com <- as.dist(matrizST.h_com, diag = FALSE, upper = FALSE)
matrizST.h_com

#
matrizST.lh_com <- read.delim("matriz_ST.LH_com.csv", sep = ";", header = TRUE, row.names = 1) 
matrizST.lh_com <- as.dist(matrizST.lh_com, diag = FALSE, upper = FALSE)
matrizST.lh_com

# WN.repl
matrizwn.rep_com <- read.delim("matriz_WN_repl_com.csv", sep = ";", header = TRUE, row.names = 1) 
matrizwn.rep_com <- as.dist(matrizwn.rep_com, diag = FALSE, upper = FALSE)
matrizwn.rep_com

# WN.rich
matrizwn.ric_com <- read.delim("matriz_WN_rich_com.csv", sep = ";", header = TRUE, row.names = 1) 
matrizwn.ric_com <- as.dist(matrizwn.ric_com, diag = FALSE, upper = FALSE)
matrizwn.ric_com

# Run dbRDA
# 500 m buffer
cobertura_500a <- decostand(x = cobertura_500a, method = "standardize") 
dbrda.os_com_500 <- dbrda(matrizos_com ~ AREA + NP + ED + ENM, 
                          data = cobertura_500a)
dbrda.os_com_500
plot(dbrda.os_com_500)
summary(dbrda.os_com_500)
anova.cca(dbrda.os_com_500, by = "terms", permutations = how(nperm = 9999))
RsquareAdj(dbrda.os_com_500)
vif.cca(dbrda.os_com_500)

# Plotting
p1_c <- ggord.dbrda(dbrda.os_com_500, ptslab = TRUE, size = 1, addsize = 3, repel = TRUE) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  theme(legend.position = 'right') +
  labs(title = expression(beta[OS])) 
p1_c

##
dbrda.st_com_500 <- dbrda(matrizst_com ~ AREA + NP + ED + ENM, data = cobertura_500a)
anova.cca(dbrda.st_com_500, by = "terms", permutations = how(nperm = 9999))
RsquareAdj(dbrda.st_com_500)
summary(dbrda.st_com_500)

# Plotting
p3_c <- ggord.dbrda(dbrda.st_com_500, ptslab = TRUE, size = 1, addsize = 3, repel = TRUE, title = "St") +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  theme(legend.position = 'right') +
  labs(title = expression(beta[ST]))  
p3_c

##
dbrda.wt_com_500 <- dbrda(matrizwt_com ~ AREA + NP + ED + ENM, 
                          data = cobertura_500a)
summary(dbrda.wt_com_500)
anova.cca(dbrda.wt_com_500, by = "terms", permutations = how(nperm = 9999))
RsquareAdj(dbrda.wt_com_500)

# Plotting
p4_c <- ggord.dbrda(dbrda.wt_com_500, ptslab = TRUE, size = 1, addsize = 3, repel = TRUE) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  theme(legend.position = 'right') +
  labs(title = expression(beta[WN])) 
p4_c

##
dbrda.ST.l_com_500 <- dbrda(matrizST.l_com ~ AREA + NP + ED + ENM, 
                            data = cobertura_500a)
summary(dbrda.ST.l_com_500)
anova.cca(dbrda.ST.l_com_500, by = "terms", permutations = how(nperm = 9999))
RsquareAdj(dbrda.ST.l_com_500)

# Plotting
p5_c <- ggord.dbrda(dbrda.ST.l_com_500, ptslab = TRUE, size = 1, 
                    addsize = 3, repel = TRUE) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  theme(legend.position = 'right') +
  labs(title = expression(beta[ST_l])) 
p5_c

##
dbrda.ST.h_com_500 <- dbrda(matrizST.h_com ~ AREA + NP + ED + ENM, 
                            data = cobertura_500a)
summary(dbrda.ST.h_com_500)
anova.cca(dbrda.ST.h_com_500, by = "terms", permutations = how(nperm = 9999))
RsquareAdj(dbrda.ST.h_com_500)

# Plotting
p6_c <- ggord.dbrda(dbrda.ST.h_com_500, ptslab = TRUE, size = 1, 
                    addsize = 3, repel = TRUE) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  theme(legend.position = 'right') +
  labs(title = expression(beta[ST_h])) 
p6_c

##
dbrda.ST.lh_com_500 <- dbrda(matrizST.lh_com ~ AREA + NP + ED + ENM, 
                             data = cobertura_500a)
summary(dbrda.ST.lh_com_500)
anova.cca(dbrda.ST.lh_com_500, by = "terms", permutations = how(nperm = 9999))
RsquareAdj(dbrda.ST.lh_com_500)

# Plotting
p7_c <- ggord.dbrda(dbrda.ST.lh_com_500, ptslab = TRUE, size = 1, 
                    addsize = 3, repel = TRUE) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  theme(legend.position = 'right') +
  labs(title = expression(beta[ST_lh])) 
p7_c

## WN.repl
dbrda.WN.rep_com_500 <- dbrda(matrizwn.rep_com ~ AREA + NP + ED + ENM, 
                              data = cobertura_500a)
summary(dbrda.WN.rep_com_500)
anova.cca(dbrda.WN.rep_com_500, by = "terms", permutations = how(nperm = 9999))
RsquareAdj(dbrda.WN.rep_com_500)

# Plotting
p8_c <- ggord.dbrda(dbrda.WN.rep_com_500, ptslab = TRUE, size = 1, 
                    addsize = 3, repel = TRUE) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  theme(legend.position = 'right') +
  labs(title = expression(beta[WN_repl])) 
p8_c

## WN.rich
dbrda.WN.ric_com_500 <- dbrda(matrizwn.ric_com ~ AREA + NP + ED + ENM, 
                              data = cobertura_500a)
summary(dbrda.WN.ric_com_500)
anova.cca(dbrda.WN.ric_com_500, by = "terms", permutations = how(nperm = 9999))
RsquareAdj(dbrda.WN.ric_com_500)

# Plotting
p10_c <- ggord.dbrda(dbrda.WN.ric_com_500, ptslab = TRUE, size = 1, 
                     addsize = 3, repel = TRUE) +
  geom_hline(yintercept = 0, linetype = 4) +
  geom_vline(xintercept = 0, linetype = 4) +
  theme(legend.position = 'right') +
  labs(title = expression(beta[WN_ric])) 
p10_c

# Combine plots
combined_plot <- (p1_c + p3_c + p4_c) +
  plot_annotation(tag_levels = "a") +
  plot_layout(nrow = 1)

# Display
print(combined_plot)

# Save
ggsave(filename = "Figure_4.jpg", plot = combined_plot, width = 18,
       height = 10, units = "cm", dpi = 300, device = "jpeg")

# Supplementary figure
combined_plot6 <- (p5_c + p6_c + p7_c + p8_c + p10_c) +
  plot_annotation(tag_levels = "a") +
  plot_layout(nrow = 2)

print(combined_plot6)

# Save
ggsave(filename = "Figure_S5.jpg", plot = combined_plot6, width = 22,
       height = 15, units = "cm", dpi = 300, device = "jpeg")

# 1500 m buffer
##
dbrda.os_com_1500 <- dbrda(matrizos_com ~ AREA_1500 + NP_1500 + ED_1500 + ENM_1500, 
                           data = cobertura_1500)
anova.cca(dbrda.os_com_1500, by = "terms", permutations = how(nperm = 9999))
RsquareAdj(dbrda.os_com_1500)
summary(dbrda.os_com_1500)
vif.cca(dbrda.os_com_1500)

# Plotting
p12_c <- ggord.dbrda(dbrda.os_com_1500, ptslab = TRUE, size = 1, 
                     addsize = 3, repel = TRUE) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  theme(legend.position = 'right') +
  labs(title = expression(beta[OS])) 
p12_c

##
dbrda.st_com_1500 <- dbrda(matrizst_com ~ AREA_1500 + NP_1500 + ED_1500 + ENM_1500, 
                           data = cobertura_1500)
anova.cca(dbrda.st_com_1500, by = "terms", permutations = how(nperm = 9999))
anova.cca(dbrda.st_com_1500, by = "margin", permutations = how(nperm = 9999))
RsquareAdj(dbrda.st_com_1500)
summary(dbrda.st_com_1500)

# Plotting
p14_c <- ggord.dbrda(dbrda.st_com_1500, ptslab = TRUE, size = 1, addsize = 3, repel = TRUE) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  theme(legend.position = 'right') +
  labs(title = expression(beta[ST])) 
p14_c

##
dbrda.wt_com_1500 <- dbrda(matrizwt_com ~ AREA_1500 + NP_1500 + ED_1500 + ENM_1500, 
                           data = cobertura_1500)
anova.cca(dbrda.wt_com_1500, by = "terms", permutations = how(nperm = 9999))
anova.cca(dbrda.wt_com_1500, by = "margin", permutations = how(nperm = 9999))
RsquareAdj(dbrda.wt_com_1500)
summary(dbrda.wt_com_1500)

# Plotting
p15_c <- ggord.dbrda(dbrda.wt_com_1500, ptslab = TRUE, size = 1, addsize = 3, repel = TRUE) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  theme(legend.position = 'right') +
  labs(title = expression(beta[WN])) 
p15_c

##
dbrda.ST.l_com_1500 <- dbrda(matrizST.l_com ~ AREA_1500 + NP_1500 + ED_1500 + ENM_1500, 
                             data = cobertura_1500)
summary(dbrda.ST.l_com_1500)
anova.cca(dbrda.ST.l_com_1500, by = "terms", permutations = how(nperm = 9999))
RsquareAdj(dbrda.ST.l_com_1500)

# Plotting
p16_c <- ggord.dbrda(dbrda.ST.l_com_1500, ptslab = TRUE, size = 1, 
                     addsize = 3, repel = TRUE) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  theme(legend.position = 'right') +
  labs(title = expression(beta[ST_l])) 
p16_c

##
dbrda.ST.h_com_1500 <- dbrda(matrizST.h_com ~ AREA_1500 + NP_1500 + ED_1500 + ENM_1500, 
                             data = cobertura_1500)
summary(dbrda.ST.h_com_1500)
anova.cca(dbrda.ST.h_com_1500, by = "terms", permutations = how(nperm = 9999))
RsquareAdj(dbrda.ST.h_com_1500)

# Plotting
p17_c <- ggord.dbrda(dbrda.ST.h_com_1500, ptslab = TRUE, size = 1, 
                     addsize = 3, repel = TRUE) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  theme(legend.position = 'right') +
  labs(title = expression(beta[ST_h])) 
p17_c

##
dbrda.ST.lh_com_1500 <- dbrda(matrizST.lh_com ~ AREA_1500 + NP_1500 + ED_1500 + ENM_1500, 
                              data = cobertura_1500)
summary(dbrda.ST.lh_com_1500)
anova.cca(dbrda.ST.lh_com_1500, by = "terms", permutations = how(nperm = 9999))
RsquareAdj(dbrda.ST.lh_com_1500)

# Plotting
p18_c <- ggord.dbrda(dbrda.ST.lh_com_1500, ptslab = TRUE, size = 1, 
                     addsize = 3, repel = TRUE) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  theme(legend.position = 'right') +
  labs(title = expression(beta[ST_lh])) 
p18_c

## WN.repl
dbrda.WN.rep_com_1500 <- dbrda(matrizwn.rep_com ~ AREA_1500 + NP_1500 + ED_1500 + ENM_1500, 
                               data = cobertura_1500)
summary(dbrda.WN.rep_com_1500)
anova.cca(dbrda.WN.rep_com_1500, by = "terms", permutations = how(nperm = 9999))
RsquareAdj(dbrda.WN.rep_com_1500)

# Plotting
p19_c <- ggord.dbrda(dbrda.WN.rep_com_1500, ptslab = TRUE, size = 1, 
                     addsize = 3, repel = TRUE) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  theme(legend.position = 'right') +
  labs(title = expression(beta[WN_repl])) 
p19_c

## WN.rich
dbrda.WN.ric_com_1500 <- dbrda(matrizwn.ric_com ~ AREA_1500 + NP_1500 + ED_1500 + ENM_1500, 
                               data = cobertura_1500)
summary(dbrda.WN.ric_com_1500)
anova.cca(dbrda.WN.ric_com_1500, by = "terms", permutations = how(nperm = 9999))
RsquareAdj(dbrda.WN.ric_com_1500)

# Plotting
p21_c <- ggord.dbrda(dbrda.WN.ric_com_1500, ptslab = TRUE, size = 1, 
                     addsize = 3, repel = TRUE) +
  geom_hline(yintercept = 0, linetype = 4) +
  geom_vline(xintercept = 0, linetype = 4) +
  theme(legend.position = 'right') +
  labs(title = expression(beta[WN_ric])) 
p21_c

# Supplementary figure
combined_plot7 <- (p16_c + p17_c + p18_c + p19_c + p21_c) +
  plot_annotation(tag_levels = "a") +
  plot_layout(nrow = 2)

print(combined_plot7)

# Save
ggsave(filename = "Figure_S6.jpg", plot = combined_plot7, width = 22,
       height = 15, units = "cm", dpi = 300, device = "jpeg")

# 3000 m buffer
cobertura_3000 <- decostand(x = cobertura_3000, method = "standardize")

dbrda.os_com_3000 <- dbrda(matrizos_com ~ AREA_3000 + NP_3000 + ED_3000 + ENM_3000, 
                           data = cobertura_3000)
anova.cca(dbrda.os_com_3000, by = "terms", permutations = how(nperm = 9999))
RsquareAdj(dbrda.os_com_3000)
vif.cca(dbrda.os_com_3000)
summary(dbrda.os_com_3000)

# Plotting
p23_c <- ggord.dbrda(dbrda.os_com_3000, ptslab = TRUE, size = 1, 
                     addsize = 3, repel = TRUE) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  theme(legend.position = 'right') +
  labs(title = expression(beta[OS])) 
p23_c

##
dbrda.st_com_3000 <- dbrda(matrizst_com ~ AREA_3000 + NP_3000 + ED_3000 + ENM_3000, 
                           data = cobertura_3000)
anova.cca(dbrda.st_com_3000, by = "terms", permutations = how(nperm = 9999))
RsquareAdj(dbrda.st_com_3000)
summary(dbrda.st_com_3000)

# Plotting
p25_c <- ggord.dbrda(dbrda.st_com_3000, ptslab = TRUE, size = 1, 
                     addsize = 3, repel = TRUE) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  theme(legend.position = 'right') +
  labs(title = expression(beta[ST])) 
p25_c

##
dbrda.wt_com_3000 <- dbrda(matrizwt_com ~ AREA_3000 + NP_3000 + ED_3000 + ENM_3000, 
                           data = cobertura_3000)
anova.cca(dbrda.wt_com_3000, by = "terms", permutations = how(nperm = 9999))
RsquareAdj(dbrda.wt_com_3000)
summary(dbrda.wt_com_3000)

# Plotting
p26_c <- ggord.dbrda(dbrda.wt_com_3000, ptslab = TRUE, size = 1, 
                     addsize = 3, repel = TRUE) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  theme(legend.position = 'right') +
  labs(title = expression(beta[WN])) 
p26_c

##
dbrda.ST.l_com_3000 <- dbrda(matrizST.l_com ~ AREA_3000 + NP_3000 + ED_3000 + ENM_3000, 
                             data = cobertura_3000)
summary(dbrda.ST.l_com_3000)
anova.cca(dbrda.ST.l_com_3000, by = "terms", permutations = how(nperm = 9999))
RsquareAdj(dbrda.ST.l_com_3000)

# Plotting
p27_c <- ggord.dbrda(dbrda.ST.l_com_3000, ptslab = TRUE, size = 1, 
                     addsize = 3, repel = TRUE) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  theme(legend.position = 'right') +
  labs(title = expression(beta[ST_l])) 
p27_c

##
dbrda.ST.h_com_3000 <- dbrda(matrizST.h_com ~ AREA_3000 + NP_3000 + ED_3000 + ENM_3000, 
                             data = cobertura_3000)
summary(dbrda.ST.h_com_3000)
anova.cca(dbrda.ST.h_com_3000, by = "terms", permutations = how(nperm = 9999))
RsquareAdj(dbrda.ST.h_com_3000)

# Plotting
p28_c <- ggord.dbrda(dbrda.ST.h_com_3000, ptslab = TRUE, size = 1, 
                     addsize = 3, repel = TRUE) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  theme(legend.position = 'right') +
  labs(title = expression(beta[ST_h])) 
p28_c

##
dbrda.ST.lh_com_3000 <- dbrda(matrizST.lh_com ~ AREA_3000 + NP_3000 + ED_3000 + ENM_3000, 
                              data = cobertura_3000)
summary(dbrda.ST.lh_com_3000)
anova.cca(dbrda.ST.lh_com_3000, by = "terms", permutations = how(nperm = 9999))
RsquareAdj(dbrda.ST.lh_com_3000)

# Plotting
p29_c <- ggord.dbrda(dbrda.ST.lh_com_3000, ptslab = TRUE, size = 1, 
                     addsize = 3, repel = TRUE) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  theme(legend.position = 'right') +
  labs(title = expression(beta[ST_lh])) 
p29_c

## WN.repl
dbrda.WN.rep_com_3000 <- dbrda(matrizwn.rep_com ~ AREA_3000 + NP_3000 + ED_3000 + ENM_3000, 
                               data = cobertura_3000)
summary(dbrda.WN.rep_com_3000)
anova.cca(dbrda.WN.rep_com_3000, by = "terms", permutations = how(nperm = 9999))
RsquareAdj(dbrda.WN.rep_com_3000)

# Plotting
p30_c <- ggord.dbrda(dbrda.WN.rep_com_3000, ptslab = TRUE, size = 1, 
                     addsize = 3, repel = TRUE) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  theme(legend.position = 'right') +
  labs(title = expression(beta[WN_repl])) 
p30_c

## WN.rich
dbrda.WN.ric_com_3000 <- dbrda(matrizwn.ric_com ~ AREA_3000 + NP_3000 + ED_3000 + ENM_3000, 
                               data = cobertura_3000)
summary(dbrda.WN.ric_com_3000)
anova.cca(dbrda.WN.ric_com_3000, by = "terms", permutations = how(nperm = 9999))
RsquareAdj(dbrda.WN.ric_com_3000)

# Plotting
p32_c <- ggord.dbrda(dbrda.WN.ric_com_3000, ptslab = TRUE, size = 1, 
                     addsize = 3, repel = TRUE) +
  geom_hline(yintercept = 0, linetype = 4) +
  geom_vline(xintercept = 0, linetype = 4) +
  theme(legend.position = 'right') +
  labs(title = expression(beta[WN_ric])) 
p32_c

# Supplementary figure
combined_plot8 <- (p27_c + p28_c + p29_c + p30_c + p32_c) +
  plot_annotation(tag_levels = "a") +
  plot_layout(nrow = 2)

print(combined_plot8)

# Save
ggsave(filename = "Figure_S7.jpg", plot = combined_plot8, width = 22,
       height = 15, units = "cm", dpi = 300, device = "jpeg")
