###############################################
# LIBRERÍAS
###############################################

packages <- c("flexclust", "lpSolve", "mclust", "aricode", 
              "cluster", "proxy", "ggplot2")
installed <- rownames(installed.packages())
to_install <- setdiff(packages, installed)
if (length(to_install) > 0) install.packages(to_install)

library(lpSolve)
library(proxy)
library(cluster)
library(mclust)
library(aricode)
library(ggplot2)

###############################################
# UTILIDADES
###############################################
prepare_data <- function(dataset) {
  dataset <- as.data.frame(dataset)
  X <- dataset
  
  list(X = X)
}

###############################################
# MILP
###############################################

solve_milp_assignment <- function(data, centroids, size_constraints) {
  n <- nrow(data)
  k <- nrow(centroids)
  
  cost_matrix <- proxy::dist(data, centroids, method = "cosine")
  cost_vec <- as.vector(t(as.matrix(cost_matrix)))
  
  # Restricción 1: cada punto debe asignarse a un solo clúster
  constr1 <- matrix(0, n, n*k)
  for (i in 1:n) {
    constr1[i, ((i - 1) * k + 1):(i * k)] <- 1
  }
  
  # Restricción 2: cada clúster debe tener tamaño EXACTO
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
# KM-MILP COMPLETO
###############################################

clustering_with_size_constraints <- function(data, size_constraints, max_iter = 100, tol = 1e-6) {
  data <- as.matrix(data)
  n <- nrow(data)
  d <- ncol(data)
  k <- length(size_constraints)
  
  # Inicialización de centroides
  idx <- sample(1:n, k)
  centroids <- data[idx, , drop = FALSE]
  
  for (iter in 1:max_iter) {
    milp_res <- solve_milp_assignment(data, centroids, size_constraints)
    p <- milp_res$p
    
    new_centroids <- matrix(0, nrow = k, ncol = d)
    for (j in 1:k) {
      pts <- data[p == j, , drop = FALSE]
      new_centroids[j, ] <- colMeans(pts)
    }
    
    if (max(abs(centroids - new_centroids)) < tol) break
    centroids <- new_centroids
  }
  
  list(p = p, centroids = centroids)
}

###############################################
# GENERAR RANGOS DE CARDINALIDAD CON DELTA
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

###############################################
# GENERAR TODAS LAS CARDINALIDADES POSIBLES
# QUE SUMAN N Y RESPETAN LOS RANGOS
###############################################

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
# EVALUAR UNA CONFIGURACIÓN DE CARDINALIDAD
###############################################

evaluate_cardinality_solution <- function(X, y, size_constraints, target_cardinality) {
  resultat <- clustering_with_size_constraints(X, size_constraints)
  labels <- resultat$p
  
  # Silueta
  dist_mat <- proxy::dist(X, method = "cosine")
  sil <- silhouette(labels, dist_mat)
  sil_mean <- mean(sil[, 3])
  
  # Cardinalidad real
  card_real <- as.numeric(table(factor(labels, levels = 1:length(size_constraints))))
  
  # Violación
  violation <- sum(abs(card_real - target_cardinality))
  
  list(
    silhouette = sil_mean,
    size_violation = violation,
    card_real = card_real
  )
}

###############################################
# EXPLORACIÓN DEL ESPACIO CON DELTA
###############################################

run_flexible_cardinality_search <- function(dataset, target_cardinality, delta, dataset_name = "dataset") {
  data_prep <- prepare_data(dataset)
  X <- data_prep$X
  y <- data_prep$y
  
  n <- nrow(X)
  k <- length(target_cardinality)
  
  ranges <- compute_cardinality_ranges(target_cardinality, delta)
  lower <- ranges$lower
  upper <- ranges$upper
  
  card_list <- generate_cardinalities(lower, upper, n)
  cat("Número de configuraciones:", length(card_list), "\n")
  
  results <- vector("list", length(card_list))
  
  for (i in seq_along(card_list)) {
    cat("Evaluando config", i, "de", length(card_list), "\n")
    
    eval <- evaluate_cardinality_solution(X, y, card_list[[i]], target_cardinality)
    
    results[[i]] <- data.frame(
      dataset = dataset_name,
      config_id = i,
      silhouette = eval$silhouette,
      size_violation = eval$size_violation,
      cardinality = paste(card_list[[i]], collapse = "-"),
      stringsAsFactors = FALSE
    )
  }
  
  results_df <- do.call(rbind, results)
  results_df
}

###############################################
# FRENTE DE PARETO
###############################################

compute_pareto_front <- function(df) {
  n <- nrow(df)
  is_pareto <- rep(TRUE, n)
  
  for (i in 1:n) {
    for (j in 1:n) {
      if (i == j) next
      
      if (df$silhouette[j] >= df$silhouette[i] &&
          df$size_violation[j] <= df$size_violation[i] &&
          (df$silhouette[j] > df$silhouette[i] ||
           df$size_violation[j] < df$size_violation[i])) {
        is_pareto[i] <- FALSE
        break
      }
    }
  }
  
  df$is_pareto <- is_pareto
  df
}

###############################################
# GRÁFICA DE PARETO
###############################################

plot_pareto <- function(df) {
  df <- compute_pareto_front(df)
  pareto <- df[df$is_pareto, ]
  pareto <- pareto[order(pareto$size_violation), ]
  
  ggplot(df, aes(x = size_violation, y = silhouette)) +
    geom_point(alpha = 0.6) +
    geom_point(data = pareto, color = "red", size = 3) +
    geom_line(data = pareto, aes(x = size_violation, y = silhouette), color = "red") +
    xlab("Violación de tamaño") +
    ylab("Coeficiente de silueta") +
    ggtitle("Frente de Pareto: violación vs. silueta") +
    theme_minimal()
}

###############################################
# EJEMPLO COMPLETO
###############################################

#data(iris)
#iris$class <- iris$Species
#iris$Species <- NULL
library(readr)
wine <- read_csv("C:/Users/emont/OneDrive/Escritorio/Tesis/wine/wine.data", 
                 col_names = FALSE)
wine=wine[,-1]
target_cardinality <- c(59, 71, 48)
delta <- c(0.05, 0.05, 0.05)

res <- run_flexible_cardinality_search(wine, target_cardinality, delta, dataset_name = "wine")

###############################################
# IMPRIMIR CARDINALIDAD ORIGINAL Y LÍMITES
###############################################

# Calcular límites con delta
ranges <- compute_cardinality_ranges(target_cardinality, delta)
lower <- ranges$lower
upper <- ranges$upper

cat("\n========================================\n")
cat(" CARDINALIDAD ORIGINAL Y LÍMITES DELTA\n")
cat("========================================\n")

cat("Cardinalidad original:      ", paste(target_cardinality, collapse = " - "), "\n")
cat("Límites inferiores (delta): ", paste(lower, collapse = " - "), "\n")
cat("Límites superiores (delta): ", paste(upper, collapse = " - "), "\n")

###############################################
# IMPRIMIR LAS 3 MEJORES SOLUCIONES POR SILUETA
###############################################

cat("\n========================================\n")
cat(" TOP 5 MEJORES SOLUCIONES (por silueta)\n")
cat("========================================\n")

best5 <- res[order(-res$silhouette), ]
print(best5)

# Graficar Pareto
print(plot_pareto(res))