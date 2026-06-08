################################################################################
################# Completeness of Interaction Networks (Chao2) #################
################################################################################

library(iNEXT)

# List of the 12 complete networks
redes_lista <- list(dados1, dados2, dados3, dados4, dados5, dados6,
                    dados7, dados8, dados9, dados10, dados11, dados12)
names(redes_lista) <- c("Brega", "Coruja", "Davi", "Faroeste", "Gruta_Indio", 
                        "Mastodonte", "Nadinho", "Olho_dagua", "Paca", 
                        "Santo_Antonio", "Sumidouro", "Zizinho")

# Function to calculate completeness using Chao2
calc_completude <- function(matriz) {
  mat_bin <- ifelse(matriz > 0, 1, 0)
  S_obs <- sum(mat_bin)
  
  # Count rare interactions (singletons and doubletons)
  # Based on interaction frequency per plant species (or bat)
  freq_plant <- colSums(mat_bin)
  singletons <- sum(freq_plant == 1)
  doubletons <- sum(freq_plant == 2)
  
  if (doubletons > 0) {
    S_est <- S_obs + (singletons^2) / (2 * doubletons)
  } else if (singletons > 0) {
    S_est <- S_obs + (singletons * (singletons - 1)) / 2
  } else {
    S_est <- S_obs
  }
  
  return(c(S_obs = S_obs, S_est = S_est, completude = S_obs / S_est))
}

# Apply to all networks
resultados <- t(sapply(redes_lista, calc_completude))
resultados <- as.data.frame(resultados)
resultados$Rede <- rownames(resultados)

print(resultados)
cat("\nMean completeness:", mean(resultados$completude))
cat("\nStandard deviation:", sd(resultados$completude))

################################################################################
################# Completeness of Interaction Networks #########################
################################################################################

# Load required packages
library(iNEXT)
library(bipartite)

################################################################################
# Create a list with all interaction matrices (complete networks)
################################################################################

# List of 12 networks (each network contains interactions from both seasons)
redes_completas <- list(
  Brega = dados1,
  Coruja = dados2,
  Davi = dados3,
  Faroeste = dados4,
  Gruta_Indio = dados5,
  Mastodonte = dados6,
  Nadinho = dados7,
  Olho_dagua = dados8,
  Paca = dados9,
  Santo_Antonio = dados10,
  Sumidouro = dados11,
  Zizinho = dados12
)

# List of 24 seasonal networks (wet and dry seasons separated)
redes_sazonais <- list(
  Brega_wet = dados1a, Brega_dry = dados1b,
  Coruja_wet = dados2a, Coruja_dry = dados2b,
  Davi_wet = dados3a, Davi_dry = dados3b,
  Faroeste_wet = dados4a, Faroeste_dry = dados4b,
  Gruta_Indio_wet = dados5a, Gruta_Indio_dry = dados5b,
  Mastodonte_wet = dados6a, Mastodonte_dry = dados6b,
  Nadinho_wet = dados7a, Nadinho_dry = dados7b,
  Olho_dagua_wet = dados8a, Olho_dagua_dry = dados8b,
  Paca_wet = dados9a, Paca_dry = dados9b,
  Santo_Antonio_wet = dados10a, Santo_Antonio_dry = dados10b,
  Sumidouro_wet = dados11a, Sumidouro_dry = dados11b,
  Zizinho_wet = dados12a, Zizinho_dry = dados12b
)

################################################################################
# Function to calculate completeness of a binary network (presence-absence)
################################################################################

calcular_completude <- function(matriz) {
  # Convert to presence-absence (if not already)
  matriz_pa <- ifelse(matriz > 0, 1, 0)
  
  # Observed number of interactions
  obs_interacoes <- sum(matriz_pa)
  
  # Using Chao2 estimator manually
  # For interaction networks, we can treat as an incidence matrix
  # where rows = bats, columns = plants
  
  n_bats <- nrow(matriz_pa)
  n_plants <- ncol(matriz_pa)
  
  # Number of observed interactions
  S_obs <- obs_interacoes
  
  # Simple method: use Chao2 estimator based on species frequency
  # Count how many species (bats or plants) have only 1 interaction (singletons) and 2 interactions (doubletons)
  
  # Interaction frequency per bat
  freq_bat <- rowSums(matriz_pa)
  singletons_bat <- sum(freq_bat == 1)
  doubletons_bat <- sum(freq_bat == 2)
  
  # Interaction frequency per plant
  freq_plant <- colSums(matriz_pa)
  singletons_plant <- sum(freq_plant == 1)
  doubletons_plant <- sum(freq_plant == 2)
  
  # Chao2 estimator for the entire matrix (as an approximation)
  Q1 <- singletons_bat + singletons_plant  # Rare species (occur in 1 interaction)
  Q2 <- doubletons_bat + doubletons_plant   # Species with 2 interactions
  
  if (Q2 > 0) {
    S_est <- S_obs + (Q1^2) / (2 * Q2)
  } else if (Q1 > 0) {
    S_est <- S_obs + (Q1 * (Q1 - 1)) / 2
  } else {
    S_est <- S_obs
  }
  
  completude <- S_obs / S_est
  
  return(data.frame(
    Network = deparse(substitute(matriz)),
    Observed_interactions = S_obs,
    Estimated_interactions = round(S_est, 2),
    Completeness = round(completude, 3),
    Q1 = Q1,
    Q2 = Q2
  ))
}

################################################################################
# Apply function to all networks
################################################################################

# Results for complete networks
resultados_completas <- data.frame()
for (nome in names(redes_completas)) {
  res <- data.frame(
    Network = nome,
    Observed_interactions = sum(redes_completas[[nome]] > 0),
    Estimated_interactions = NA,
    Completeness = NA
  )
  resultados_completas <- rbind(resultados_completas, res)
}

# Results for seasonal networks
resultados_sazonais <- data.frame()
for (nome in names(redes_sazonais)) {
  matriz <- redes_sazonais[[nome]]
  obs <- sum(matriz > 0)
  res <- data.frame(
    Network = nome,
    Observed_interactions = obs,
    Estimated_interactions = NA,
    Completeness = NA
  )
  resultados_sazonais <- rbind(resultados_sazonais, res)
}

################################################################################
# Alternative: Use bootstrap to calculate completeness
################################################################################

calcular_completude_bootstrap <- function(matriz, n_sim = 100) {
  matriz_pa <- ifelse(matriz > 0, 1, 0)
  obs <- sum(matriz_pa)
  
  # Simulate interaction accumulation curves
  n_interacoes_possiveis <- nrow(matriz_pa) * ncol(matriz_pa)
  
  # Bootstrap estimation
  estimativas <- numeric(n_sim)
  for (i in 1:n_sim) {
    # Sampling with replacement from matrix cells
    amostra <- sample(matriz_pa, n_interacoes_possiveis, replace = TRUE)
    estimativas[i] <- sum(amostra)
  }
  
  S_est <- mean(estimativas)
  completude <- obs / S_est
  
  return(data.frame(
    Observed_interactions = obs,
    Estimated_interactions = round(S_est, 2),
    Completeness = round(completude, 3),
    Completeness_sd = round(sd(estimativas) / S_est, 3)
  ))
}

################################################################################
# Simplified version using adapted `vegan::specpool` approach
################################################################################

# Transform each network into a list of species (interactions as "species")

calcular_completude_interacoes <- function(matriz) {
  # Convert to binary
  mat_bin <- ifelse(matriz > 0, 1, 0)
  
  # Observed number of interactions
  obs_interacoes <- sum(mat_bin)
  
  # Transpose matrix: samples = bats, species = plants
  incid_matrix <- t(mat_bin)
  
  # Check if iNEXT package is available
  if (require(iNEXT, quietly = TRUE)) {
    # Estimate interaction richness using iNEXT
    # Note: iNEXT works best with abundance data, but can be adapted
    # For presence-absence, we use the incidence method
    tryCatch({
      # Use Chao2 estimator
      S_obs <- obs_interacoes
      
      # Species with 1 occurrence (interactions occurring with only 1 bat)
      singletons <- sum(rowSums(incid_matrix) == 1)
      
      # Species with 2 occurrences
      doubletons <- sum(rowSums(incid_matrix) == 2)
      
      if (doubletons > 0) {
        S_est <- S_obs + (singletons^2) / (2 * doubletons)
      } else if (singletons > 0) {
        S_est <- S_obs + (singletons * (singletons - 1)) / 2
      } else {
        S_est <- S_obs
      }
      
      completude <- S_obs / S_est
      
      return(data.frame(
        Observed_interactions = S_obs,
        Estimated_interactions = round(S_est, 2),
        Completeness = round(completude, 3),
        Singletons = singletons,
        Doubletons = doubletons
      ))
    }, error = function(e) {
      return(data.frame(
        Observed_interactions = obs_interacoes,
        Estimated_interactions = NA,
        Completeness = NA,
        Singletons = NA,
        Doubletons = NA
      ))
    })
  } else {
    return(data.frame(
      Observed_interactions = obs_interacoes,
      Estimated_interactions = NA,
      Completeness = NA,
      Singletons = NA,
      Doubletons = NA
    ))
  }
}

################################################################################
# Run for all networks
################################################################################

# Complete networks
cat("\n========== COMPLETENESS OF COMPLETE NETWORKS ==========\n")
for (i in 1:length(redes_completas)) {
  nome <- names(redes_completas)[i]
  cat("\n", nome, ":\n")
  resultado <- calcular_completude_interacoes(redes_completas[[i]])
  print(resultado)
}

# Seasonal networks
cat("\n\n========== COMPLETENESS OF SEASONAL NETWORKS ==========\n")
resultados_sazonais_df <- data.frame()
for (i in 1:length(redes_sazonais)) {
  nome <- names(redes_sazonais)[i]
  cat("\n", nome, ":\n")
  resultado <- calcular_completude_interacoes(redes_sazonais[[i]])
  resultado$Network <- nome
  resultados_sazonais_df <- rbind(resultados_sazonais_df, resultado)
  print(resultado)
}

# Descriptive statistics
cat("\n\n========== SUMMARY STATISTICS ==========\n")
cat("\nMean completeness (complete networks):\n")
# Calculate mean

cat("\nMean completeness (seasonal networks):\n")
mean(resultados_sazonais_df$Completeness, na.rm = TRUE)
sd(resultados_sazonais_df$Completeness, na.rm = TRUE)

# Identify networks with low completeness (< 70%)
cat("\nNetworks with completeness < 70%:\n")
resultados_sazonais_df[resultados_sazonais_df$Completeness < 0.7, ]