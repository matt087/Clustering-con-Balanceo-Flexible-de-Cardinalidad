###############################################
# LIBRERÍAS
###############################################

packages <- c("flexclust", "lpSolve", "mclust", "aricode", 
              "cluster", "proxy", "ggplot2", "readr", "dplyr")
installed <- rownames(installed.packages())
to_install <- setdiff(packages, installed)
if (length(to_install) > 0) install.packages(to_install)

library(lpSolve)
library(proxy)
library(cluster)
library(ggplot2)
library(dplyr)

###############################################
# UTILIDADES
###############################################

prepare_data <- function(dataset) {
  dataset <- as.data.frame(dataset)
  X <- dataset
  list(X = X)
}

set.seed(22)
###############################################
# MILP
###############################################

solve_milp_assignment <- function(data, centroids, size_constraints) {
  n <- nrow(data)
  k <- nrow(centroids)
  
  cost_matrix <- proxy::dist(data, centroids, method = "cosine")
  cost_vec <- as.vector(t(as.matrix(cost_matrix)))
  
  constr1 <- matrix(0, n, n*k)
  for (i in 1:n) {
    constr1[i, ((i - 1) * k + 1):(i * k)] <- 1
  }
  
  constr2 <- matrix(0, k, n*k)
  for (j in 1:k) {
    constr2[j, seq(j, n*k, by = k)] <- 1
  }
  
  f.con <- rbind(constr1, constr2)
  f.dir <- c(rep("=", n), rep("=", k))
  f.rhs <- c(rep(1, n), size_constraints)
  
  result <- lp("min", cost_vec, f.con, f.dir, f.rhs, all.bin = TRUE)
  
  if (result$status != 0) {
    stop("MILP no encontró solución factible para estas cardinalidades.")
  }
  
  x_opt <- matrix(result$solution, nrow = n, byrow = TRUE)
  p <- apply(x_opt, 1, which.max)
  
  list(p = p)
}

###############################################
# KM-MILP
###############################################

clustering_with_size_constraints <- function(data, size_constraints, max_iter = 100, tol = 1e-6) {
  data <- as.matrix(data)
  n <- nrow(data)
  d <- ncol(data)
  k <- length(size_constraints)
  
  idx <- sample(1:n, k)
  centroids <- data[idx, , drop = FALSE]
  
  for (iter in 1:max_iter) {
    milp_res <- solve_milp_assignment(data, centroids, size_constraints)
    p <- milp_res$p
    
    new_centroids <- matrix(0, nrow = k, ncol = d)
    for (j in 1:k) {
      pts <- data[p == j, , drop = FALSE]
      if (nrow(pts) > 0) new_centroids[j, ] <- colMeans(pts)
      else new_centroids[j, ] <- centroids[j, ]
    }
    
    if (max(abs(centroids - new_centroids)) < tol) break
    centroids <- new_centroids
  }
  
  list(p = p, centroids = centroids)
}

###############################################
# FUNCIONES DE CARDINALIDAD
###############################################

compute_cardinality_ranges <- function(target, delta) {
  k <- length(target)
  
  if (length(delta) == 1) delta <- rep(delta, k)
  if (length(delta) != k) stop("delta debe tener longitud 1 o k.")
  
  lower <- floor(target * (1 - delta))
  upper <- round(target * (1 + delta))
  lower[lower < 1] <- 1
  
  list(lower = lower, upper = upper)
}

generate_cardinalities <- function(lower, upper, total_n) {
  k <- length(lower)
  res <- list()
  current <- integer(k)
  
  rec <- function(pos, remaining) {
    if (pos == k) {
      if (remaining >= lower[k] && remaining <= upper[k]) {
        current[k] <<- remaining
        res[[length(res) + 1]] <<- current
      }
      return()
    }
    
    min_rest <- sum(lower[(pos + 1):k])
    max_rest <- sum(upper[(pos + 1):k])
    
    for (v in lower[pos]:upper[pos]) {
      rest <- remaining - v
      if (rest < min_rest || rest > max_rest) next
      current[pos] <<- v
      rec(pos + 1, rest)
    }
  }
  
  rec(1, total_n)
  res
}

###############################################
# EVALUACIÓN (SILUETA + CSVI)
###############################################

evaluate_cardinality_solution <- function(X, y, size_constraints, target_cardinality) {
  resultat <- clustering_with_size_constraints(X, size_constraints)
  labels <- resultat$p
  
  dist_mat <- proxy::dist(X, method = "cosine")
  sil <- silhouette(labels, dist_mat)
  sil_mean <- mean(sil[, 3])
  
  card_real <- as.numeric(table(factor(labels, levels = 1:length(size_constraints))))
  
  ilvc <- sum(abs(card_real - target_cardinality))
  clvc <- sum(card_real != target_cardinality)
  
  n_total <- sum(card_real)
  k_total <- length(card_real)
  
  alpha <- 0.5
  norm_ilvc <- ilvc / n_total
  norm_clvc <- clvc / k_total
  csvi <- alpha * norm_ilvc + (1 - alpha) * norm_clvc
  
  list(
    silhouette = sil_mean,
    ilvc = ilvc,
    clvc = clvc,
    csvi = csvi
  )
}

###############################################
# MULTIOBJETIVO: PARETO + MO-SCORE (A MAXIMIZAR)
###############################################

pareto_dominates <- function(a, b) {
  (a$silhouette >= b$silhouette && a$CSVI <= b$CSVI) &&
    (a$silhouette >  b$silhouette || a$CSVI   < b$CSVI)
}

find_pareto_front <- function(df) {
  n <- nrow(df)
  dominated <- rep(FALSE, n)
  
  for (i in 1:n) {
    for (j in 1:n) {
      if (i != j) {
        if (pareto_dominates(df[j,], df[i,])) {
          dominated[i] <- TRUE
          break
        }
      }
    }
  }
  
  df[!dominated, ]
}

multiobjective_score <- function(sil, csvi, 
                                 w_sil, 
                                 w_csvi,
                                 sil_range = c(-1, 1)) {
  sil_min <- sil_range[1]
  sil_max <- sil_range[2]
  sil_norm <- (sil - sil_min) / (sil_max - sil_min)
  sil_norm <- max(min(sil_norm, 1), 0)  
  
  # Normalizar CSVI
  csvi_norm <- max(min(csvi, 1), 0)
  
  score <- w_sil * sil_norm + w_csvi * (1 - csvi_norm)
  return(score)
}

###############################################
# EXPLORACIÓN COMPLETA DEL ESPACIO
###############################################

run_flexible_cardinality_search <- function(dataset, 
                                            target_cardinality, 
                                            delta, 
                                            dataset_name = "dataset",
                                            w_sil = 0.5,
                                            w_csvi = 0.5) {
  
  data_prep <- prepare_data(dataset)
  X <- data_prep$X
  
  n <- nrow(X)
  k <- length(target_cardinality)
  
  ranges <- compute_cardinality_ranges(target_cardinality, delta)
  lower <- ranges$lower
  upper <- ranges$upper
  
  card_list <- generate_cardinalities(lower, upper, n)
  cat("Número de configuraciones posibles:", length(card_list), "\n")
  
  results <- list()
  
  for (i in seq_along(card_list)) {
    if (i %% 10 == 0 || i == 1) cat("Evaluando config", i, "de", length(card_list), "\n")
    
    eval <- evaluate_cardinality_solution(X, NULL, card_list[[i]], target_cardinality)
    
    results[[i]] <- data.frame(
      dataset = dataset_name,
      config_id = i,
      silhouette = eval$silhouette,
      ILVC = eval$ilvc,
      CLVC = eval$clvc,
      CSVI = eval$csvi,
      cardinality = paste(card_list[[i]], collapse = "-")
    )
  }
  
  df <- bind_rows(results)
  
  # MO_Score a maximizar 
  df$MO_Score <- mapply(
    multiobjective_score,
    sil  = df$silhouette,
    csvi = df$CSVI,
    MoreArgs = list(
      w_sil = w_sil,
      w_csvi = w_csvi,
      sil_range = c(-1, 1) 
    )
  )
  
  df
}

###############################################
# GRÁFICA PARETO CON PARETO FRONT
###############################################

plot_pareto <- function(df) { 
  df$Solution_ID <- factor(1:nrow(df))
  
  max_sil <- max(df$silhouette)
  max_csvi <- max(df$CSVI)
  scale_factor <- max_sil / max_csvi
  
  ggplot(df, aes(x = Solution_ID)) +
    geom_col(
      aes(y = silhouette),
      fill = "#5B9BD5", color = "#1F4E79", alpha = 0.75
    ) +
    geom_line(
      aes(y = CSVI * scale_factor, group = 1),
      color = "#ED7D31", size = 1
    ) +
    geom_point(
      aes(y = CSVI * scale_factor),
      color = "#ED7D31", size = 2
    ) +
    scale_y_continuous(
      name = "Silueta",
      sec.axis = sec_axis(~ . / scale_factor, name = "Violación CSVI")
    ) +
    labs(
      title = "Gráfica de Pareto",
      x = "Soluciones",
      subtitle = "Barras = Silueta, Línea = CSVI"
    ) +
    theme_minimal()
}

###############################################
# EJEMPLO DE USO
###############################################


#wine <- read_csv("tesis/wine.data",
              #   col_names = FALSE)
data(iris)
iris <- iris[,-5]
target_cardinality <- c(50, 50, 50)
delta <- c(0.05, 0.05, 0.05) #tolerancia de la restrigción
w_sil <- 0.9 #peso del coeficiente de silueta
w_csvi <- 0.1 #peso de la violación de tamaño (CSVI)

res <- run_flexible_cardinality_search(iris, 
                                       target_cardinality, 
                                       delta, 
                                       dataset_name = "iris",
                                       w_sil, 
                                       w_csvi)

pareto_front <- find_pareto_front(res)
best_solution <- pareto_front[which.max(pareto_front$MO_Score), ]

cat("\n==============================\n")
cat(" FRENTE PARETO\n")
cat("==============================\n")
print(pareto_front)

cat("\n==============================\n")
cat(" MEJOR SOLUCIÓN MULTIOBJETIVO (MO_Score MÁXIMO)\n")
cat("==============================\n")
print(best_solution)

print(plot_pareto(res))