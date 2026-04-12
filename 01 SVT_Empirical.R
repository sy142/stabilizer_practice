rm(list=ls())
suppressPackageStartupMessages({
  library(readxl)
  library(lavaan)
  library(dplyr)
  library(metafor)
  library(boot)
})

set.seed(9186)
data <- read_excel("C:/Users/Salim/Desktop/Yildiz Teknik YL Tez/Makale/EC Revs/Rev2/Structural Equation Modeling Companion/Datasets/analizliksonAVE.xlsx")
data$BRIAN <- data$BRO_Toplam

cfa_model <- '
  MEQ_Factor =~ MEQ1 + MEQ2 + MEQ8 + MEQ9 + MEQ10 + MEQ11 + MEQ15 + MEQ17 + MEQ18 + MEQ19
  MEQ1 ~~ MEQ2
  MEQ2 ~~ MEQ10
  MEQ11 ~~ MEQ15
  MEQ17 ~~ MEQ18
  MEQ2 ~~ MEQ8
'

moderators <- c(
  "YasGruplari", "Cinsiyet", "EgitimDuzeyi", "MedeniDurum",
  "GelirDuzeyi", "BKISiniflamasi", "SigaraKullanma",
  "AlkolKullanma", "KronikHastalikVarligi", "TibbiBeslenmeTedavisiAlma"
)

compute_adaptive_weights <- function(data, moderators, cfa_model) {
  cat("\n[ADAPTIVE CALIBRATION] Computing data-driven weight system...\n")
  
  all_deltas_metric <- list()
  all_deltas_scalar <- list()
  all_invariance <- c()
  moderator_names <- c()
  
  for(mod in moderators) {
    if(!(mod %in% names(data))) next
    
    d <- data[!is.na(data[[mod]]), ]
    d[[mod]] <- factor(d[[mod]])
    levs <- levels(d[[mod]])
    
    if(length(levs) < 2) next
    
    fit_config <- tryCatch(
      cfa(cfa_model, data=d, group=mod, std.lv=TRUE, estimator="MLR", missing="fiml"),
      error=function(e) NULL
    )
    if(is.null(fit_config)) next
    
    fit_metric <- tryCatch(
      cfa(cfa_model, data=d, group=mod, group.equal="loadings",
          std.lv=TRUE, estimator="MLR", missing="fiml"),
      error=function(e) NULL
    )
    if(is.null(fit_metric)) next
    
    fit_scalar <- tryCatch(
      cfa(cfa_model, data=d, group=mod, 
          group.equal=c("loadings", "intercepts"),
          std.lv=TRUE, estimator="MLR", missing="fiml"),
      error=function(e) NULL
    )
    if(is.null(fit_scalar)) next
    
    fm_config <- fitmeasures(fit_config)
    fm_metric <- fitmeasures(fit_metric)
    fm_scalar <- fitmeasures(fit_scalar)
    
    get_delta <- function(fm1, fm2, index) {
      val1 <- as.numeric(fm1[index])
      val2 <- as.numeric(fm2[index])
      if(is.na(val1) || is.na(val2)) return(NA)
      return(abs(val1 - val2))
    }
    
    delta_cfi_metric <- get_delta(fm_config, fm_metric, "cfi.scaled")
    delta_tli_metric <- get_delta(fm_config, fm_metric, "tli.scaled")
    delta_rmsea_metric <- get_delta(fm_config, fm_metric, "rmsea.scaled")
    delta_srmr_metric <- get_delta(fm_config, fm_metric, "srmr")
    
    delta_cfi_scalar <- get_delta(fm_metric, fm_scalar, "cfi.scaled")
    delta_tli_scalar <- get_delta(fm_metric, fm_scalar, "tli.scaled")
    delta_rmsea_scalar <- get_delta(fm_metric, fm_scalar, "rmsea.scaled")
    delta_srmr_scalar <- get_delta(fm_metric, fm_scalar, "srmr")
    
    if(any(is.na(c(delta_cfi_metric, delta_tli_metric, delta_rmsea_metric, delta_srmr_metric,
                   delta_cfi_scalar, delta_tli_scalar, delta_rmsea_scalar, delta_srmr_scalar)))) next
    
    all_deltas_metric[[mod]] <- c(CFI = delta_cfi_metric, TLI = delta_tli_metric, 
                                  RMSEA = delta_rmsea_metric, SRMR = delta_srmr_metric)
    all_deltas_scalar[[mod]] <- c(CFI = delta_cfi_scalar, TLI = delta_tli_scalar,
                                  RMSEA = delta_rmsea_scalar, SRMR = delta_srmr_scalar)
    
    threshold_cfi <- 0.010
    threshold_tli <- 0.010
    threshold_rmsea <- 0.015
    threshold_srmr <- 0.030
    
    is_invariant <- (delta_cfi_metric <= threshold_cfi && delta_rmsea_metric <= threshold_rmsea &&
                       delta_srmr_metric <= threshold_srmr && delta_tli_metric <= threshold_tli &&
                       delta_cfi_scalar <= threshold_cfi && delta_rmsea_scalar <= threshold_rmsea &&
                       delta_srmr_scalar <= threshold_srmr && delta_tli_scalar <= threshold_tli)
    
    all_invariance <- c(all_invariance, is_invariant)
    moderator_names <- c(moderator_names, mod)
  }
  
  if(length(all_deltas_metric) < 3) {
    cat("[ADAPTIVE] Not enough data. Using default weights.\n")
    return(c(0.40, 0.10, 0.30, 0.20))
  }
  
  delta_metric_matrix <- do.call(rbind, all_deltas_metric)
  delta_scalar_matrix <- do.call(rbind, all_deltas_scalar)
  delta_max_matrix <- pmax(delta_metric_matrix, delta_scalar_matrix)
  
  cat("[ADAPTIVE] Analyzing", nrow(delta_max_matrix), "moderators\n")
  cat("[ADAPTIVE] Invariant:", sum(all_invariance), "Non-invariant:", sum(!all_invariance), "\n\n")
  
  correlation_matrix <- cor(delta_max_matrix)
  cat("Correlation Matrix of Delta Values:\n")
  print(round(correlation_matrix, 3))
  
  normalize_deltas <- function(deltas, threshold) {
    return(pmin(deltas / threshold, 1))
  }
  
  norm_matrix <- cbind(
    CFI = normalize_deltas(delta_max_matrix[, "CFI"], 0.010),
    TLI = normalize_deltas(delta_max_matrix[, "TLI"], 0.010),
    RMSEA = normalize_deltas(delta_max_matrix[, "RMSEA"], 0.015),
    SRMR = normalize_deltas(delta_max_matrix[, "SRMR"], 0.030)
  )
  
  cv_values <- apply(norm_matrix, 2, function(x) sd(x) / (mean(x) + 0.001))
  
  cat("\nCoefficient of Variation for each index:\n")
  cat(" CFI:", round(cv_values["CFI"], 3), "\n")
  cat(" TLI:", round(cv_values["TLI"], 3), "\n")
  cat(" RMSEA:", round(cv_values["RMSEA"], 3), "\n")
  cat(" SRMR:", round(cv_values["SRMR"], 3), "\n")
  
  mean_cor <- apply(abs(correlation_matrix), 1, function(x) mean(x[x < 1]))
  redundancy_penalty <- 1 / (1 + mean_cor)
  
  discriminant_power <- numeric(4)
  names(discriminant_power) <- c("CFI", "TLI", "RMSEA", "SRMR")
  
  for(idx in c("CFI", "TLI", "RMSEA", "SRMR")) {
    inv_mean <- mean(norm_matrix[all_invariance, idx])
    if(sum(!all_invariance) > 0) {
      non_inv_mean <- mean(norm_matrix[!all_invariance, idx])
      raw_diff <- non_inv_mean - inv_mean
      discriminant_power[idx] <- max(0, raw_diff)
    } else {
      discriminant_power[idx] <- sd(norm_matrix[, idx])
    }
  }
  
  cat("\nDiscriminant Power (separation between groups):\n")
  cat(" CFI:", round(discriminant_power["CFI"], 3), "\n")
  cat(" TLI:", round(discriminant_power["TLI"], 3), "\n")
  cat(" RMSEA:", round(discriminant_power["RMSEA"], 3), "\n")
  cat(" SRMR:", round(discriminant_power["SRMR"], 3), "\n")
  
  variability_score <- 1 / (1 + cv_values)
  raw_weights <- discriminant_power * redundancy_penalty * variability_score
  final_weights <- raw_weights / sum(raw_weights)
  
  cat("\n[ADAPTIVE] Final Calibrated Weights:\n")
  cat(" CFI:", round(final_weights["CFI"], 4), "\n")
  cat(" TLI:", round(final_weights["TLI"], 4), "\n")
  cat(" RMSEA:", round(final_weights["RMSEA"], 4), "\n")
  cat(" SRMR:", round(final_weights["SRMR"], 4), "\n")
  
  return(as.numeric(final_weights))
}

calibrated_weights <- compute_adaptive_weights(data, moderators, cfa_model)

test_measurement_invariance_adaptive <- function(data, moderator, cfa_model, weights = calibrated_weights) {
  d <- data[!is.na(data[[moderator]]), ]
  d[[moderator]] <- factor(d[[moderator]])
  
  if(length(unique(d[[moderator]])) < 2) return(NULL)
  
  fit_config <- tryCatch(
    cfa(cfa_model, data=d, group=moderator, std.lv=TRUE, estimator="MLR", missing="fiml"),
    error=function(e) NULL
  )
  
  if(is.null(fit_config)) return(NULL)
  
  fit_metric <- tryCatch(
    cfa(cfa_model, data=d, group=moderator, 
        group.equal="loadings",
        std.lv=TRUE, estimator="MLR", missing="fiml"),
    error=function(e) NULL
  )
  
  fit_scalar <- tryCatch(
    cfa(cfa_model, data=d, group=moderator, 
        group.equal=c("loadings", "intercepts"),
        std.lv=TRUE, estimator="MLR", missing="fiml"),
    error=function(e) NULL
  )
  
  fit_strict <- tryCatch(
    cfa(cfa_model, data=d, group=moderator,
        group.equal=c("loadings", "intercepts", "residuals"),
        std.lv=TRUE, estimator="MLR", missing="fiml"),
    error=function(e) NULL
  )
  
  if(is.null(fit_metric) || is.null(fit_scalar) || is.null(fit_strict)) return(NULL)
  
  fm_config <- fitMeasures(fit_config, c("cfi.scaled", "tli.scaled", "rmsea.scaled", "srmr"))
  fm_metric <- fitMeasures(fit_metric, c("cfi.scaled", "tli.scaled", "rmsea.scaled", "srmr"))
  fm_scalar <- fitMeasures(fit_scalar, c("cfi.scaled", "tli.scaled", "rmsea.scaled", "srmr"))
  fm_strict <- fitMeasures(fit_strict, c("cfi.scaled", "tli.scaled", "rmsea.scaled", "srmr"))
  
  get_delta <- function(fm1, fm2, index) {
    val1 <- as.numeric(fm1[index])
    val2 <- as.numeric(fm2[index])
    if(is.na(val1) || is.na(val2)) return(NA)
    return(val1 - val2)
  }
  
  delta_cfi_metric <- get_delta(fm_config, fm_metric, "cfi.scaled")
  delta_tli_metric <- get_delta(fm_config, fm_metric, "tli.scaled")
  delta_rmsea_metric <- -get_delta(fm_config, fm_metric, "rmsea.scaled")
  delta_srmr_metric <- -get_delta(fm_config, fm_metric, "srmr")
  
  delta_cfi_scalar <- get_delta(fm_metric, fm_scalar, "cfi.scaled")
  delta_tli_scalar <- get_delta(fm_metric, fm_scalar, "tli.scaled")
  delta_rmsea_scalar <- -get_delta(fm_metric, fm_scalar, "rmsea.scaled")
  delta_srmr_scalar <- -get_delta(fm_metric, fm_scalar, "srmr")
  
  delta_cfi_strict <- get_delta(fm_scalar, fm_strict, "cfi.scaled")
  delta_tli_strict <- get_delta(fm_scalar, fm_strict, "tli.scaled")
  delta_rmsea_strict <- -get_delta(fm_scalar, fm_strict, "rmsea.scaled")
  delta_srmr_strict <- -get_delta(fm_scalar, fm_strict, "srmr")
  
  if(any(is.na(c(delta_cfi_metric, delta_tli_metric, delta_rmsea_metric, delta_srmr_metric,
                 delta_cfi_scalar, delta_tli_scalar, delta_rmsea_scalar, delta_srmr_scalar,
                 delta_cfi_strict, delta_tli_strict, delta_rmsea_strict, delta_srmr_strict)))) return(NULL)
  
  threshold_cfi <- 0.010
  threshold_tli <- 0.010
  threshold_rmsea <- 0.015
  threshold_srmr <- 0.030
  
  norm_cfi <- min(max(pmax(0, delta_cfi_metric), pmax(0, delta_cfi_scalar), pmax(0, delta_cfi_strict)) / threshold_cfi, 1)
  norm_tli <- min(max(pmax(0, delta_tli_metric), pmax(0, delta_tli_scalar), pmax(0, delta_tli_strict)) / threshold_tli, 1)
  norm_rmsea <- min(max(pmax(0, delta_rmsea_metric), pmax(0, delta_rmsea_scalar), pmax(0, delta_rmsea_strict)) / threshold_rmsea, 1)
  norm_srmr <- min(max(pmax(0, delta_srmr_metric), pmax(0, delta_srmr_scalar), pmax(0, delta_srmr_strict)) / threshold_srmr, 1)
  
  mi_score <- weights[1] * norm_cfi + weights[2] * norm_tli + weights[3] * norm_rmsea + weights[4] * norm_srmr
  
  is_invariant <- (delta_cfi_metric <= threshold_cfi && delta_rmsea_metric <= threshold_rmsea &&
                     delta_srmr_metric <= threshold_srmr && delta_tli_metric <= threshold_tli &&
                     delta_cfi_scalar <= threshold_cfi && delta_rmsea_scalar <= threshold_rmsea &&
                     delta_srmr_scalar <= threshold_srmr && delta_tli_scalar <= threshold_tli &&
                     delta_cfi_strict <= threshold_cfi && delta_rmsea_strict <= threshold_rmsea &&
                     delta_srmr_strict <= threshold_srmr && delta_tli_strict <= threshold_tli)
  
  groups <- levels(d[[moderator]])
  group_contributions <- list()
  
  for(stage in c("metric", "scalar", "strict")) {
    if(stage == "metric") {
      fit_base <- fit_config
      fit_restricted <- fit_metric
      stage_deltas <- c(CFI = delta_cfi_metric, TLI = delta_tli_metric,
                        RMSEA = delta_rmsea_metric, SRMR = delta_srmr_metric)
    } else if(stage == "scalar") {
      fit_base <- fit_metric
      fit_restricted <- fit_scalar
      stage_deltas <- c(CFI = delta_cfi_scalar, TLI = delta_tli_scalar,
                        RMSEA = delta_rmsea_scalar, SRMR = delta_srmr_scalar)
    } else {
      fit_base <- fit_scalar
      fit_restricted <- fit_strict
      stage_deltas <- c(CFI = delta_cfi_strict, TLI = delta_tli_strict,
                        RMSEA = delta_rmsea_strict, SRMR = delta_srmr_strict)
    }
    
    if(abs(stage_deltas["CFI"]) > threshold_cfi || abs(stage_deltas["TLI"]) > threshold_tli ||
       abs(stage_deltas["RMSEA"]) > threshold_rmsea || abs(stage_deltas["SRMR"]) > threshold_srmr) {
      
      for(g in groups) {
        group_data <- d[d[[moderator]] == g, ]
        
        fit_g <- tryCatch(
          cfa(cfa_model, data=group_data, std.lv=TRUE, estimator="MLR", missing="fiml"),
          error=function(e) NULL
        )
        
        if(!is.null(fit_g)) {
          fm_g <- fitMeasures(fit_g, c("cfi.scaled", "tli.scaled", "rmsea.scaled", "srmr"))
          
          if(!is.null(group_contributions[[g]])) {
            group_contributions[[g]][[stage]] <- list(
              cfi = fm_g["cfi.scaled"],
              tli = fm_g["tli.scaled"],
              rmsea = fm_g["rmsea.scaled"],
              srmr = fm_g["srmr"],
              n = nrow(group_data)
            )
          } else {
            group_contributions[[g]] <- list()
            group_contributions[[g]][[stage]] <- list(
              cfi = fm_g["cfi.scaled"],
              tli = fm_g["tli.scaled"],
              rmsea = fm_g["rmsea.scaled"],
              srmr = fm_g["srmr"],
              n = nrow(group_data)
            )
          }
        }
      }
    }
  }
  
  for(g in names(group_contributions)) {
    contrib <- group_contributions[[g]]
    group_mi_score <- 0
    
    for(stage in names(contrib)) {
      stage_contrib <- contrib[[stage]]
      stage_score <- 0
      
      if(!is.na(stage_contrib$cfi)) {
        cfi_dev <- abs(0.95 - stage_contrib$cfi)
        stage_score <- stage_score + weights[1] * min(cfi_dev / threshold_cfi, 1)
      }
      if(!is.na(stage_contrib$tli)) {
        tli_dev <- abs(0.95 - stage_contrib$tli)
        stage_score <- stage_score + weights[2] * min(tli_dev / threshold_tli, 1)
      }
      if(!is.na(stage_contrib$rmsea)) {
        rmsea_dev <- max(0, stage_contrib$rmsea - 0.06)
        stage_score <- stage_score + weights[3] * min(rmsea_dev / threshold_rmsea, 1)
      }
      if(!is.na(stage_contrib$srmr)) {
        srmr_dev <- max(0, stage_contrib$srmr - 0.08)
        stage_score <- stage_score + weights[4] * min(srmr_dev / threshold_srmr, 1)
      }
      
      group_contributions[[g]][[stage]]$stage_mi_score <- stage_score
      group_mi_score <- group_mi_score + stage_score
    }
    
    group_contributions[[g]]$total_mi_score <- group_mi_score
  }
  
  return(list(
    moderator = moderator,
    mi_score = mi_score,
    is_invariant = is_invariant,
    delta_cfi_metric = delta_cfi_metric,
    delta_tli_metric = delta_tli_metric,
    delta_rmsea_metric = delta_rmsea_metric,
    delta_srmr_metric = delta_srmr_metric,
    delta_cfi_scalar = delta_cfi_scalar,
    delta_tli_scalar = delta_tli_scalar,
    delta_rmsea_scalar = delta_rmsea_scalar,
    delta_srmr_scalar = delta_srmr_scalar,
    delta_cfi_strict = delta_cfi_strict,
    delta_tli_strict = delta_tli_strict,
    delta_rmsea_strict = delta_rmsea_strict,
    delta_srmr_strict = delta_srmr_strict,
    group_contributions = group_contributions
  ))
}

calculate_stabilization_metrics <- function(data, moderator) {
  d <- data[!is.na(data[[moderator]]), ]
  d[[moderator]] <- factor(d[[moderator]])
  groups <- levels(d[[moderator]])
  if (length(groups) < 2) return(NULL)
  
  model_without <- '
    MEQ_Factor =~ MEQ1 + MEQ2 + MEQ8 + MEQ9 + MEQ10 + MEQ11 + MEQ15 + MEQ17 + MEQ18 + MEQ19
    YJIPAQ ~ c*MEQ_Factor
    MEQ1 ~~ MEQ2
    MEQ2 ~~ MEQ10
    MEQ11 ~~ MEQ15
    MEQ17 ~~ MEQ18
    MEQ2 ~~ MEQ8
  '
  
  model_with <- '
    MEQ_Factor =~ MEQ1 + MEQ2 + MEQ8 + MEQ9 + MEQ10 + MEQ11 + MEQ15 + MEQ17 + MEQ18 + MEQ19
    BRIAN ~ a*MEQ_Factor
    YJIPAQ ~ cprime*MEQ_Factor + b*BRIAN
    MEQ1 ~~ MEQ2
    MEQ2 ~~ MEQ10
    MEQ11 ~~ MEQ15
    MEQ17 ~~ MEQ18
    MEQ2 ~~ MEQ8
  '
  
  model_meq_brian <- '
    MEQ_Factor =~ MEQ1 + MEQ2 + MEQ8 + MEQ9 + MEQ10 + MEQ11 + MEQ15 + MEQ17 + MEQ18 + MEQ19
    BRIAN ~ meqb*MEQ_Factor
    MEQ1 ~~ MEQ2
    MEQ2 ~~ MEQ10
    MEQ11 ~~ MEQ15
    MEQ17 ~~ MEQ18
    MEQ2 ~~ MEQ8
  '
  
  c_paths <- c(); cprime_paths <- c(); meqb_paths <- c()
  group_ns <- c(); group_names <- c()
  
  for (g in groups) {
    group_data <- d[d[[moderator]] == g, ]
    if (nrow(group_data) < 30) {
      cat("  [Warning] Moderator:", moderator, "| Group:", g, "| n =", nrow(group_data), "< 30 (excluded from analysis)\n")
      next
    }
    
    fit_without <- tryCatch(
      sem(model_without, data = group_data, std.lv = TRUE, estimator = "MLR", missing = "fiml"),
      error = function(e) NULL
    )
    fit_with <- tryCatch(
      sem(model_with, data = group_data, std.lv = TRUE, estimator = "MLR", missing = "fiml"),
      error = function(e) NULL
    )
    fit_meqb <- tryCatch(
      sem(model_meq_brian, data = group_data, std.lv = TRUE, estimator = "MLR", missing = "fiml"),
      error = function(e) NULL
    )
    
    if (!is.null(fit_without) && !is.null(fit_with) && !is.null(fit_meqb)) {
      std_without <- standardizedSolution(fit_without)
      std_with    <- standardizedSolution(fit_with)
      std_meqb    <- standardizedSolution(fit_meqb)
      
      c_path      <- std_without[std_without$label == "c", "est.std"][1]
      cprime_path <- std_with[std_with$label == "cprime", "est.std"][1]
      meqb_path   <- std_meqb[std_meqb$label == "meqb", "est.std"][1]
      
      if (!is.na(c_path) && !is.na(cprime_path) && !is.na(meqb_path)) {
        c_paths      <- c(c_paths, c_path)
        cprime_paths <- c(cprime_paths, cprime_path)
        meqb_paths   <- c(meqb_paths, meqb_path)
        group_ns     <- c(group_ns, nrow(group_data))
        group_names  <- c(group_names, g)
      }
    }
  }
  
  if (length(c_paths) < 2) return(NULL)
  
  wmean <- function(x, w) weighted.mean(x, w)
  wsd   <- function(x, w) sqrt(sum(w * (x - wmean(x, w))^2) / sum(w))
  
  mean_c      <- wmean(c_paths,      group_ns)
  mean_cprime <- wmean(cprime_paths, group_ns)
  mean_meqb   <- wmean(meqb_paths,   group_ns)
  
  sd_c      <- wsd(c_paths,      group_ns)
  sd_cprime <- wsd(cprime_paths, group_ns)
  sd_meqb   <- wsd(meqb_paths,   group_ns)
  
  epsilon_mu <- 0.0001
  mean_c_protected <- sign(mean_c) * max(abs(mean_c), epsilon_mu)
  mean_cprime_protected <- sign(mean_cprime) * max(abs(mean_cprime), epsilon_mu)
  mean_meqb_protected <- sign(mean_meqb) * max(abs(mean_meqb), epsilon_mu)
  
  cv_without <- (sd_c / abs(mean_c_protected)) * 100
  cv_with    <- (sd_cprime / abs(mean_cprime_protected)) * 100
  cv_meqb    <- (sd_meqb / abs(mean_meqb_protected)) * 100
  
  eps_log <- 1e-8
  delta_log <- log(cv_without + eps_log) - log(cv_with + eps_log)
  pct_reduction <- 100 * (1 - exp(-delta_log))
  
  I_g <- numeric(length(group_names))
  for(i in seq_along(group_names)) {
    dist_before <- abs(c_paths[i] - mean_c)
    dist_after <- abs(cprime_paths[i] - mean_cprime)
    I_g[i] <- ifelse(dist_after < dist_before, 1, 0)
  }
  
  orientation_consistency_ratio <- mean(I_g)
  
  # Dlog_mu: log(|mu1|) - log(|mu0|)  -> positive values indicate mean alignment
  # Dlog_sigma: log(sigma0) - log(sigma1)  -> positive values indicate variance reduction
  # Dl = Dlog_sigma + Dlog_mu  (exact decomposition)
  # OS = max(0, Dlog_mu) / (max(0, Dlog_mu) + max(0, Dlog_sigma) + eps)
  eps_os <- 1e-8
  delta_log_mu    <- log(abs(mean_cprime) + eps_os) - log(abs(mean_c) + eps_os)
  delta_log_sigma <- log(sd_c + eps_os) - log(sd_cprime + eps_os)
  orientation_share <- max(0, delta_log_mu) / (max(0, delta_log_mu) + max(0, delta_log_sigma) + eps_os)

  group_contributions <- list()
  for (i in seq_along(group_names)) {
    group_contributions[[group_names[i]]] <- list(
      c_path = c_paths[i],
      cprime_path = cprime_paths[i],
      meqb_path = meqb_paths[i],
      c_deviation = abs(c_paths[i] - mean_c),
      cprime_deviation = abs(cprime_paths[i] - mean_cprime),
      meqb_deviation = abs(meqb_paths[i] - mean_meqb),
      stabilization_contribution = (abs(c_paths[i] - mean_c) - abs(cprime_paths[i] - mean_cprime)),
      moderation_strength = abs(meqb_paths[i]),
      I_g = I_g[i],
      n = group_ns[i]
    )
  }
  
  return(list(
    moderator = moderator,
    n_groups = length(c_paths),
    total_n = sum(group_ns),
    cv_without = cv_without,
    cv_with = cv_with,
    cv_meqb = cv_meqb,
    delta_log = delta_log,
    pct_reduction = pct_reduction,
    orientation_consistency_ratio = orientation_consistency_ratio,
    orientation_share = orientation_share,
    delta_log_mu = delta_log_mu,
    delta_log_sigma = delta_log_sigma,
    group_contributions = group_contributions
  ))
}

stabilization_variable_test <- function(data, moderators) {
  cat("\n==================== STABILIZATION VARIABLE TEST (SVT) ====================\n")
  cat("Reference: Yilmaz & Cene (2025), Mathematics, Algorithm C1\n")
  cat("\nTesting whether BRIAN acts as a stabilizer (structurally independent, MI-linked)\n")
  cat("H0: E[Delta_l] <= 0 (no stabilization)\n")
  cat("H1: E[Delta_l]  > 0 (stabilization present)\n")
  cat("Decision: Reject H0 iff p_boot < alpha AND p_binom < alpha (dual-criterion)\n")
  cat("-----------------------------------------------------------------------------\n\n")
  
  all_results <- list()
  
  for (mod in moderators) {
    if (!(mod %in% names(data))) next
    
    mi_result  <- test_measurement_invariance_adaptive(data, mod, cfa_model, calibrated_weights)
    stab_result <- calculate_stabilization_metrics(data, mod)
    
    if (!is.null(mi_result) && !is.null(stab_result)) {
      all_results[[mod]] <- list(
        moderator = mod,
        mi_score = mi_result$mi_score,
        is_invariant = mi_result$is_invariant,
        delta_log = stab_result$delta_log,
        pct_reduction = stab_result$pct_reduction,
        cv_without = stab_result$cv_without,
        cv_with = stab_result$cv_with,
        cv_meqb = stab_result$cv_meqb,
        n_groups = stab_result$n_groups,
        total_n = stab_result$total_n,
        orientation_consistency_ratio = stab_result$orientation_consistency_ratio,
        orientation_share = stab_result$orientation_share,           # EKLENEN
        delta_log_mu = stab_result$delta_log_mu,                     # EKLENEN
        delta_log_sigma = stab_result$delta_log_sigma,               # EKLENEN
        group_contributions = stab_result$group_contributions,
        mi_group_contributions = mi_result$group_contributions
      )
      
      cat(sprintf("%-25s: MI=%s (%.3f), Dl=%6.3f, CV: %.1f%%->%.1f%% (%+.1f%%), OCR=%.2f, OS=%.3f\n",
                  mod,
                  ifelse(mi_result$is_invariant, "OK ", "BAD"),
                  mi_result$mi_score,
                  stab_result$delta_log,
                  stab_result$cv_without,
                  stab_result$cv_with,
                  stab_result$pct_reduction,
                  stab_result$orientation_consistency_ratio,
                  stab_result$orientation_share))
    }
  }
  
  if (length(all_results) == 0) {
    cat("\nNo valid results obtained.\n")
    return(NULL)
  }
  
  cat("\n==================== STEP 2: STABILIZATION QUANTIFICATION ====================\n")
  
  delta_vec <- sapply(all_results, function(x) x$delta_log)
  weights   <- sapply(all_results, function(x) x$total_n)
  mi_scores <- sapply(all_results, function(x) x$mi_score)
  is_invariant <- sapply(all_results, function(x) x$is_invariant)
  os_vec    <- sapply(all_results, function(x) x$orientation_share)
  ocr_vec   <- sapply(all_results, function(x) x$orientation_consistency_ratio)
  
  # Criterion 1: Bootstrap test (Mathematics Section 4.1.4, Eq. 42-45)
  boot_mean <- function(d, idx) {
    w <- weights[idx] / sum(weights[idx])
    sum(d[idx] * w)
  }
  boot_res <- boot(delta_vec, boot_mean, R = 1000)
  mean_delta <- boot_res$t0
  se_delta   <- sd(boot_res$t)
  
  wvar <- sum(weights * (delta_vec - mean_delta)^2) / sum(weights)
  d_cohen <- if (is.finite(wvar) && wvar > 0) mean_delta / sqrt(wvar) else NA_real_
  
  z_stat <- mean_delta / se_delta
  p_boot <- pnorm(z_stat, lower.tail = FALSE)
  
  ci_perc <- boot.ci(boot_res, type = "perc")$percent[4:5]
  
  cat("\n--- Criterion 1: Bootstrap Test for Mean Stabilization (Eq. 42-45) ---\n")
  cat("Weighted mean Delta_l:", round(mean_delta, 3), "\n")
  cat("SE (bootstrap):       ", round(se_delta, 3), "\n")
  cat("z-statistic:          ", round(z_stat, 3), "\n")
  cat("p_boot (one-tailed):  ", round(p_boot, 4), "\n")
  cat("95% CI (percentile):  [", round(ci_perc[1], 3), ", ", round(ci_perc[2], 3), "]\n", sep = "")
  cat("Cohen's d:            ", round(d_cohen, 3), "\n")
  
  # Criterion 2: Binomial test (Mathematics Section 4.1.4, Eq. 46-47, Remark 11)
  # Per Remark 11: count moderators where (Dl_m > 0) OR (OCR_m >= 0.5)
  criterion2_pass_per_mod <- (delta_vec > 0) | (ocr_vec >= 0.5)
  S_obs <- sum(criterion2_pass_per_mod)
  M_total <- length(criterion2_pass_per_mod)
  
  btest <- binom.test(S_obs, M_total, p = 0.5, alternative = "greater")
  p_binom <- btest$p.value
  
  cat("\n--- Criterion 2: Binomial Test for Directional Consistency (Eq. 46-47) ---\n")
  cat("Per moderator: stabilized if Delta_l > 0 OR OCR >= 0.5 (Remark 11)\n")
  cat("Stabilized moderators: ", S_obs, "/", M_total, "\n")
  cat("  Delta_l > 0:  ", sum(delta_vec > 0), "/", M_total, "\n")
  cat("  OCR >= 0.5:   ", sum(ocr_vec >= 0.5), "/", M_total, "\n")
  cat("p_binom (one-tailed): ", round(p_binom, 4), "\n")
  
  # MI Complementarity (C1 diagnostic, Remark 1)
  cat("\n==================== C1 DIAGNOSTICS (Remark 1) ====================\n")
  n_mi_bad  <- sum(!is_invariant)
  n_mi_good <- sum(is_invariant)
  cat("\nMI Status Distribution:\n")
  cat("  MI Satisfied:", n_mi_good, "moderators\n")
  cat("  MI Violated :", n_mi_bad,  "moderators\n")
  
  if (n_mi_bad > 0) {
    stab_mi_bad  <- delta_vec[!is_invariant]
    stab_mi_good <- delta_vec[ is_invariant]
    cat("\nMean Delta_l by MI status:\n")
    cat("  MI violated : ", round(mean(stab_mi_bad), 3), "\n", sep = "")
    cat("  MI satisfied: ", round(mean(stab_mi_good), 3), "\n", sep = "")
  }
  
  r_mi <- suppressWarnings(cor(mi_scores, delta_vec))
  cat("\nMI-Complementarity: cor(S_MI, Delta_l) = ", round(r_mi, 3), "\n")
  cat("Cross-moderator consistency: ", sum(delta_vec > 0), "/", M_total, " moderators with Delta_l > 0\n")
  
  # Orientation Share (Mathematics Eq. 25, Section 3.4)
  mean_os <- weighted.mean(os_vec, weights)
  cat("\n==================== ORIENTATION SHARE (Eq. 25) ====================\n")
  cat("\nPer-moderator OS values:\n")
  for(i in seq_along(os_vec)) {
    cat(sprintf("  %-25s: OS = %.3f\n", names(os_vec)[i], os_vec[i]))
  }
  cat("\nWeighted mean OS:", round(mean_os, 3), "\n")
  cat("Interpretation: OS < 0.3 = Type A (Variance Purification)\n")
  cat("                OS > 0.7 = Type B (Directional Alignment)\n")
  cat("                0.3-0.7  = Type AB (Combined)\n")
  
  # ==================== DECISION (Algorithm C1, Section 4.1.4) ====================
  cat("\n==================== SVT DECISION (Algorithm C1) ====================\n\n")
  
  bootstrap_pass <- (p_boot < 0.05 && mean_delta > 0)
  binom_pass     <- (p_binom < 0.05)
  
  summary_df <- data.frame(
    moderator    = names(all_results),
    mi_satisfied = sapply(all_results, function(x) x$is_invariant),
    mi_score     = sapply(all_results, function(x) x$mi_score),
    delta_log    = sapply(all_results, function(x) x$delta_log),
    pct_reduction= sapply(all_results, function(x) x$pct_reduction),
    cv_without   = sapply(all_results, function(x) x$cv_without),
    cv_with      = sapply(all_results, function(x) x$cv_with),
    cv_meqb      = sapply(all_results, function(x) x$cv_meqb),
    orientation_consistency = sapply(all_results, function(x) x$orientation_consistency_ratio),
    orientation_share = sapply(all_results, function(x) x$orientation_share),
    n_groups     = sapply(all_results, function(x) x$n_groups),
    total_n      = sapply(all_results, function(x) x$total_n)
  )
  
  # Decision rule: p_boot < alpha AND p_binom < alpha (Eq. 48)
  if (bootstrap_pass && binom_pass) {
    
    # Mechanism classification via OS (Section 3.4)
    if (mean_os < 0.3) {
      mechanism_classification <- "Type A (Variance Purification)"
    } else if (mean_os > 0.7) {
      mechanism_classification <- "Type B (Directional Alignment)"
    } else {
      mechanism_classification <- "Type AB (Combined)"
    }
    
    cat("*** STABILIZER DETECTED ***\n\n")
    cat("Decision: BRIAN is classified as a STABILIZER VARIABLE (Definition 1)\n\n")
    cat("Criterion 1 (Bootstrap):  p_boot  = ", round(p_boot, 4), " < 0.05  [PASS]\n", sep = "")
    cat("Criterion 2 (Binomial):   p_binom = ", round(p_binom, 4), " < 0.05  [PASS]\n", sep = "")
    cat("Weighted mean Delta_l:    ", round(mean_delta, 3), "\n")
    cat("Cohen's d:                ", round(d_cohen, 3), "\n")
    cat("95% CI:                   [", round(ci_perc[1], 3), ", ", round(ci_perc[2], 3), "]\n", sep = "")
    cat("\nMechanism (OS = ", round(mean_os, 3), "): ", mechanism_classification, "\n", sep = "")
    
  } else {
    mechanism_classification <- "None (Not a stabilizer)"
    
    cat("*** NO STABILIZATION DETECTED ***\n\n")
    cat("Decision: BRIAN does NOT satisfy stabilizer variable conditions.\n\n")
    cat("Criterion 1 (Bootstrap):  p_boot  = ", round(p_boot, 4),
        ifelse(bootstrap_pass, " < 0.05  [PASS]", " >= 0.05  [FAIL]"), "\n", sep = "")
    cat("Criterion 2 (Binomial):   p_binom = ", round(p_binom, 4),
        ifelse(binom_pass, " < 0.05  [PASS]", " >= 0.05  [FAIL]"), "\n", sep = "")
  }
  
  svt_results <- list(
    # Primary test statistics
    test_statistic = z_stat,
    p_boot = p_boot,
    p_binom = p_binom,
    mean_delta_log = mean_delta,
    effect_size_d = d_cohen,
    ci_lower = ci_perc[1],
    ci_upper = ci_perc[2],
    # Mechanism
    orientation_share = mean_os,
    mean_orientation_consistency = mean(ocr_vec),
    mechanism = mechanism_classification,
    # Decision
    bootstrap_pass = bootstrap_pass,
    binom_pass = binom_pass,
    decision = ifelse(bootstrap_pass && binom_pass, "Stabilizer", "Not a stabilizer"),
    # Details
    summary = summary_df,
    detailed_results = all_results,
    calibrated_weights = calibrated_weights
  )
  
  saveRDS(svt_results, "C:/Users/Salim/Desktop/Yildiz Teknik YL Tez/Makale/EC Revs/Rev2/Structural Equation Modeling Companion/RDS/SVT_results_adaptive.rds")
  write.csv(summary_df, "C:/Users/Salim/Desktop/Yildiz Teknik YL Tez/Makale/EC Revs/Rev2/Structural Equation Modeling Companion/Outputs/SVT_summary_adaptive.csv", row.names = FALSE)
  
  cat("\nResults saved:\n")
  cat("  RDS:  .../Structural Equation Modeling Companion/RDS/SVT_results_adaptive.rds\n")
  cat("  CSV:  .../Structural Equation Modeling Companion/Outputs/SVT_summary_adaptive.csv\n")
  return(svt_results)
}

svt_results <- stabilization_variable_test(data, moderators)