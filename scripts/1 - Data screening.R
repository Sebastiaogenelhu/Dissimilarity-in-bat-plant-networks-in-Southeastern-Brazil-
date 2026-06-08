################################################################################
################# Script Network Analysis - Pains Project ######################
################# Dissimilarity in Interaction Networks ########################
################################################################################

## Loading initial packages
library(bipartite); library(amap); require(ggnetwork);
library(vegan); library(ggplot2); library(ggExtra); library(ggrepel);
library(car); library(visreg); library(multcomp); library(MuMIn);
library(betapart); library(adespatial); library(fields);
library(reshape); library(ade4); library(spdep); library(stats);
library(performance); library(ape); library(corrplot); library(ncf); 
library(ggstatsplot); 
library(tidyr); library(dplyr); library(patchwork);
library(tidyverse)  

## Setting working directory
dir()

# Remove previously created objects
rm(list = ls())

## First step: calculate network dissimilarity (Poisot et al. 2012)
# Creating work matrices
dados1 <- read.delim("Brega.txt", row.names = 1, head = TRUE)
dados1a <- read.delim("Brega_chuva.txt", row.names = 1, head = TRUE)
dados1b <- read.delim("Brega_seca.txt", row.names = 1, head = TRUE)

dados2 <- read.delim("coruja.txt", row.names = 1, head = TRUE)
dados2a <- read.delim("coruja_chuva.txt", row.names = 1, head = TRUE)
dados2b <- read.delim("coruja_seca.txt", row.names = 1, head = TRUE)

dados3 <- read.delim("Davi.txt", row.names = 1, head = TRUE)
dados3a <- read.delim("Davi_chuva.txt", row.names = 1, head = TRUE)
dados3b <- read.delim("Davi_seca.txt", row.names = 1, head = TRUE)

dados4 <- read.table("faroeste.txt", row.names = 1, head = TRUE)
dados4a <- read.table("faroeste_chuva.txt", row.names = 1, head = TRUE)
dados4b <- read.table("faroeste_seca.txt", row.names = 1, head = TRUE)

dados5 <- read.table("gruta_indio.txt", row.names = 1, head = TRUE)
dados5a <- read.table("gruta_indio_chuva.txt", row.names = 1, head = TRUE)
dados5b <- read.table("gruta_indio_seca.txt", row.names = 1, head = TRUE)
dados5b <- read.delim("gruta_indio_seca.csv", sep = ";", header = TRUE, row.names = 1)

dados6 <- read.table("Mastodonte.txt", row.names = 1, head = TRUE)
dados6a <- read.table("Mastodonte_chuva.txt", row.names = 1, head = TRUE)
dados6b <- read.table("Mastodonte_seca.txt", row.names = 1, head = TRUE)

dados7 <- read.table("nadinho.txt", row.names = 1, head = TRUE)
dados7a <- read.table("nadinho_chuva.txt", row.names = 1, head = TRUE)
dados7b <- read.table("nadinho_seca.txt", row.names = 1, head = TRUE)

dados8 <- read.table("olhodagua.txt", row.names = 1, head = TRUE)
dados8a <- read.table("olhodagua_chuva.txt", row.names = 1, head = TRUE)
dados8b <- read.table("olhodagua_seca.txt", row.names = 1, head = TRUE)

dados9 <- read.table("paca.txt", row.names = 1, head = TRUE)
dados9a <- read.table("paca_chuva.txt", row.names = 1, head = TRUE)
dados9b <- read.table("paca_seca.txt", row.names = 1, head = TRUE)

dados10 <- read.table("santoantonio.txt", row.names = 1, head = TRUE)
dados10a <- read.table("santoantonio_chuva.txt", row.names = 1, head = TRUE)
dados10b <- read.table("santoantonio_seca.txt", row.names = 1, head = TRUE)

dados11 <- read.table("sumidouro.txt", row.names = 1, head = TRUE)
dados11a <- read.table("sumidouro_chuva.txt", row.names = 1, head = TRUE)
dados11b <- read.table("sumidouro_seca.txt", row.names = 1, head = TRUE)

dados12 <- read.table("zizinho.txt", row.names = 1, head = TRUE)
dados12a <- read.table("zizinho_chuva.txt", row.names = 1, head = TRUE)
dados12b <- read.table("zizinho_seca.txt", row.names = 1, head = TRUE)

# Creating work matrices (transposed)
dados1 <- as.matrix(t(read.delim("Brega.txt", row.names = 1, head = TRUE)))
dados1a <- as.matrix(t(read.delim("Brega_chuva.txt", row.names = 1, head = TRUE)))
dados1b <- as.matrix(t(read.delim("Brega_seca.txt", row.names = 1, head = TRUE)))

dados2 <- as.matrix(t(read.delim("coruja.txt", row.names = 1, head = TRUE)))
dados2a <- as.matrix(t(read.delim("coruja_chuva.txt", row.names = 1, head = TRUE)))
dados2b <- as.matrix(t(read.delim("coruja_seca.txt", row.names = 1, head = TRUE)))

dados3 <- as.matrix(t(read.delim("Davi.txt", row.names = 1, head = TRUE)))
dados3a <- as.matrix(t(read.delim("Davi_chuva.txt", row.names = 1, head = TRUE)))
dados3b <- as.matrix(t(read.delim("Davi_seca.txt", row.names = 1, head = TRUE)))

dados4 <- as.matrix(t(read.table("faroeste.txt", row.names = 1, head = TRUE)))
dados4a <- as.matrix(t(read.table("faroeste_chuva.txt", row.names = 1, head = TRUE)))
dados4b <- as.matrix(t(read.table("faroeste_seca.txt", row.names = 1, head = TRUE)))

dados5 <- as.matrix(t(read.table("gruta_indio.txt", row.names = 1, head = TRUE)))
dados5a <- as.matrix(t(read.table("gruta_indio_chuva.txt", row.names = 1, head = TRUE)))
dados5b <- as.matrix(t(read.delim("gruta_indio_seca.csv", sep = ";", header = TRUE, row.names = 1)))

dados6 <- as.matrix(t(read.table("Mastodonte.txt", row.names = 1, head = TRUE)))
dados6a <- as.matrix(t(read.table("Mastodonte_chuva.txt", row.names = 1, head = TRUE)))
dados6b <- as.matrix(t(read.table("Mastodonte_seca.txt", row.names = 1, head = TRUE)))

dados7 <- as.matrix(t(read.table("nadinho.txt", row.names = 1, head = TRUE)))
dados7a <- as.matrix(t(read.delim("nadinho_chuva.txt", row.names = 1, head = TRUE)))
dados7b <- as.matrix(t(read.delim("nadinho_seca.txt", row.names = 1, head = TRUE)))

dados8 <- as.matrix(t(read.table("olhodagua.txt", row.names = 1, head = TRUE)))
dados8a <- as.matrix(t(read.table("olhodagua_chuva.txt", row.names = 1, head = TRUE)))
dados8b <- as.matrix(t(read.table("olhodagua_seca.txt", row.names = 1, head = TRUE)))

dados9 <- as.matrix(t(read.table("paca.txt", row.names = 1, head = TRUE)))
dados9a <- as.matrix(t(read.table("paca_chuva.txt", row.names = 1, head = TRUE)))
dados9b <- as.matrix(t(read.table("paca_seca.txt", row.names = 1, head = TRUE)))

dados10 <- as.matrix(t(read.table("santoantonio.txt", row.names = 1, head = TRUE)))
dados10a <- as.matrix(t(read.table("santoantonio_chuva.txt", row.names = 1, head = TRUE)))
dados10b <- as.matrix(t(read.table("santoantonio_seca.txt", row.names = 1, head = TRUE)))

dados11 <- as.matrix(t(read.table("sumidouro.txt", row.names = 1, head = TRUE)))
dados11a <- as.matrix(t(read.table("sumidouro_chuva.txt", row.names = 1, head = TRUE)))
dados11b <- as.matrix(t(read.table("sumidouro_seca.txt", row.names = 1, head = TRUE)))

dados12 <- as.matrix(t(read.table("zizinho.txt", row.names = 1, head = TRUE)))
dados12a <- as.matrix(t(read.table("zizinho_chuva.txt", row.names = 1, head = TRUE)))
dados12b <- as.matrix(t(read.table("zizinho_seca.txt", row.names = 1, head = TRUE)))
