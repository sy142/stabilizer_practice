rm(list = ls())
suppressPackageStartupMessages({
  library(readxl)
  library(lavaan)
  library(dplyr)
  library(stringr)
  library(boot)
  library(ggplot2)
  library(gridExtra)
})

set.seed(9186)

base_path <- "C:/Users/Salim/Desktop/Yildiz Teknik YL Tez/Makale/EC Revs/Rev2/Structural Equation Modeling Companion"
data_path   <- file.path(base_path, "Datasets")
output_path <- file.path(base_path, "Outputs")
rds_path    <- file.path(base_path, "RDS")
figure_path <- file.path(base_path, "Figures")

for (p in c(output_path, rds_path, figure_path)) {
  if (!dir.exists(p)) dir.create(p, recursive = TRUE)
}

turkce_ascii <- function(metin) {
  if (is.character(metin)) {
    r <- c("\u015f"="s","\u015e"="S","\u011f"="g","\u011e"="G",
           "\u00fc"="u","\u00dc"="U","\u00f6"="o","\u00d6"="O",
           "\u00e7"="c","\u00c7"="C","\u0131"="i","\u0130"="I")
    metin <- str_replace_all(metin, r)
  }
  metin
}

data <- read_excel(file.path(data_path, "analizliksonAVE.xlsx"))
data <- data %>% mutate(across(where(is.character), turkce_ascii))
data$BRIAN <- data$BRO_Toplam

cat("Multi-Path Negative Control SVT Analysis\n")
cat("Dataset: analizliksonAVE.xlsx, N =", nrow(data), "\n\n")

cat("Design:\n")
cat("  Positive control: BRIAN (satisfies C1 and C2)\n")
cat("  Negative controls: PSQI, KIDMED (structurally connected but no C1)\n\n")
cat("  Focal paths tested:\n")
cat("    1. MEQ -> YJIPAQ  (primary, same as main SVT analysis)\n")
cat("    2. MEQ -> PSF12   (PSQI->PSF12 exists at beta=-0.300)\n")
cat("    3. MEQ -> KIDMED  (BRIAN->KIDMED=-0.174, PSQI->KIDMED=-0.041)\n\n")

psqi_col   <- grep("^PSQI$|PSQI_Toplam|PUKIToplam|PUKI_Toplam", names(data), value = TRUE, ignore.case = TRUE)[1]
kidmed_col <- grep("^KIDMED$|KIDMED_Toplam|KIDMEDToplam", names(data), value = TRUE, ignore.case = TRUE)[1]
psf12_col  <- grep("^PSF12$|Physical_SF|PCS|FizikselBilesen", names(data), value = TRUE, ignore.case = TRUE)[1]

cat("Detected columns:\n")
cat("  PSQI   ->", ifelse(is.na(psqi_col), "NOT FOUND", psqi_col), "\n")
cat("  KIDMED ->", ifelse(is.na(kidmed_col), "NOT FOUND", kidmed_col), "\n")
cat("  PSF12  ->", ifelse(is.na(psf12_col), "NOT FOUND", psf12_col), "\n\n")

if (is.na(psqi_col) || is.na(kidmed_col) || is.na(psf12_col)) {
  cat("Available columns:\n")
  cat(paste(names(data), collapse = ", "), "\n")
  stop("One or more columns not found. Set manually above.")
}

data$PSQI_ctrl   <- data[[psqi_col]]
data$KIDMED_ctrl <- data[[kidmed_col]]
data$PSF12_ctrl  <- data[[psf12_col]]

cfa_model <- '
  MEQ_Factor =~ MEQ1 + MEQ2 + MEQ8 + MEQ9 + MEQ10 + MEQ11 + MEQ15 + MEQ17 + MEQ18 + MEQ19
  MEQ1 ~~ MEQ2
  MEQ2 ~~ MEQ10
  MEQ11 ~~ MEQ15
  MEQ17 ~~ MEQ18
  MEQ2 ~~ MEQ8
'

moderators <- c("YasGruplari","Cinsiyet","EgitimDuzeyi","MedeniDurum",
                "GelirDuzeyi","BKISiniflamasi","SigaraKullanma",
                "AlkolKullanma","KronikHastalikVarligi","TibbiBeslenmeTedavisiAlma")

compute_adaptive_weights <- function(data, moderators, cfa_model) {
  cat("[ADAPTIVE CALIBRATION] Computing data-driven weight system...\n")
  
  all_deltas_metric <- list()
  all_deltas_scalar <- list()
  all_invariance <- c()
  
  for (mod in moderators) {
    if (!(mod %in% names(data))) next
    d <- data[!is.na(data[[mod]]), ]
    d[[mod]] <- factor(d[[mod]])
    if (length(levels(d[[mod]])) < 2) next
    
    fit_config <- tryCatch(cfa(cfa_model, data=d, group=mod, std.lv=TRUE, estimator="MLR", missing="fiml"), error=function(e) NULL)
    if (is.null(fit_config)) next
    fit_metric <- tryCatch(cfa(cfa_model, data=d, group=mod, group.equal="loadings", std.lv=TRUE, estimator="MLR", missing="fiml"), error=function(e) NULL)
    if (is.null(fit_metric)) next
    fit_scalar <- tryCatch(cfa(cfa_model, data=d, group=mod, group.equal=c("loadings","intercepts"), std.lv=TRUE, estimator="MLR", missing="fiml"), error=function(e) NULL)
    if (is.null(fit_scalar)) next
    
    fm_config <- fitmeasures(fit_config)
    fm_metric <- fitmeasures(fit_metric)
    fm_scalar <- fitmeasures(fit_scalar)
    
    gd <- function(f1, f2, idx) {
      v1 <- as.numeric(f1[idx]); v2 <- as.numeric(f2[idx])
      if (is.na(v1) || is.na(v2)) NA else abs(v1 - v2)
    }
    
    dm <- c(CFI=gd(fm_config,fm_metric,"cfi.scaled"), TLI=gd(fm_config,fm_metric,"tli.scaled"),
            RMSEA=gd(fm_config,fm_metric,"rmsea.scaled"), SRMR=gd(fm_config,fm_metric,"srmr"))
    ds <- c(CFI=gd(fm_metric,fm_scalar,"cfi.scaled"), TLI=gd(fm_metric,fm_scalar,"tli.scaled"),
            RMSEA=gd(fm_metric,fm_scalar,"rmsea.scaled"), SRMR=gd(fm_metric,fm_scalar,"srmr"))
    
    if (any(is.na(c(dm, ds)))) next
    
    all_deltas_metric[[mod]] <- dm
    all_deltas_scalar[[mod]] <- ds
    
    th <- c(CFI=0.010, TLI=0.010, RMSEA=0.015, SRMR=0.030)
    is_inv <- all(dm <= th) && all(ds <= th)
    all_invariance <- c(all_invariance, is_inv)
    names(all_invariance)[length(all_invariance)] <- mod
  }
  
  if (length(all_deltas_metric) < 3) {
    cat("  Not enough data. Using default weights.\n")
    return(c(0.40, 0.10, 0.30, 0.20))
  }
  
  dm_mat <- do.call(rbind, all_deltas_metric)
  ds_mat <- do.call(rbind, all_deltas_scalar)
  dmax_mat <- pmax(dm_mat, ds_mat)
  
  cormat <- cor(dmax_mat)
  norm_mat <- cbind(CFI=pmin(dmax_mat[,"CFI"]/0.010,1), TLI=pmin(dmax_mat[,"TLI"]/0.010,1),
                    RMSEA=pmin(dmax_mat[,"RMSEA"]/0.015,1), SRMR=pmin(dmax_mat[,"SRMR"]/0.030,1))
  
  cv_v <- apply(norm_mat, 2, function(x) sd(x)/(mean(x)+0.001))
  mcor <- apply(abs(cormat), 1, function(x) mean(x[x<1]))
  rp <- 1/(1+mcor)
  
  dp <- numeric(4); names(dp) <- c("CFI","TLI","RMSEA","SRMR")
  for (idx in names(dp)) {
    if (sum(all_invariance)>0 && sum(!all_invariance)>0) {
      dp[idx] <- max(0, mean(norm_mat[!all_invariance,idx]) - mean(norm_mat[all_invariance,idx]))
    } else dp[idx] <- sd(norm_mat[,idx])
  }
  
  vs <- 1/(1+cv_v)
  raw_w <- dp*rp*vs
  final_w <- raw_w/sum(raw_w)
  
  cat("  Weights: CFI=", round(final_w["CFI"],4), " TLI=", round(final_w["TLI"],4),
      " RMSEA=", round(final_w["RMSEA"],4), " SRMR=", round(final_w["SRMR"],4), "\n\n")
  
  return(as.numeric(final_w))
}

test_mi_adaptive <- function(data, moderator, cfa_model, weights) {
  d <- data[!is.na(data[[moderator]]), ]
  d[[moderator]] <- factor(d[[moderator]])
  if (length(unique(d[[moderator]])) < 2) return(NULL)
  
  fc <- tryCatch(cfa(cfa_model, data=d, group=moderator, std.lv=TRUE, estimator="MLR", missing="fiml"), error=function(e) NULL)
  if (is.null(fc)) return(NULL)
  fm <- tryCatch(cfa(cfa_model, data=d, group=moderator, group.equal="loadings", std.lv=TRUE, estimator="MLR", missing="fiml"), error=function(e) NULL)
  fs <- tryCatch(cfa(cfa_model, data=d, group=moderator, group.equal=c("loadings","intercepts"), std.lv=TRUE, estimator="MLR", missing="fiml"), error=function(e) NULL)
  ft <- tryCatch(cfa(cfa_model, data=d, group=moderator, group.equal=c("loadings","intercepts","residuals"), std.lv=TRUE, estimator="MLR", missing="fiml"), error=function(e) NULL)
  if (is.null(fm) || is.null(fs) || is.null(ft)) return(NULL)
  
  fmc <- fitMeasures(fc, c("cfi.scaled","tli.scaled","rmsea.scaled","srmr"))
  fmm <- fitMeasures(fm, c("cfi.scaled","tli.scaled","rmsea.scaled","srmr"))
  fms <- fitMeasures(fs, c("cfi.scaled","tli.scaled","rmsea.scaled","srmr"))
  fmt <- fitMeasures(ft, c("cfi.scaled","tli.scaled","rmsea.scaled","srmr"))
  
  gds <- function(f1, f2, idx) { v1<-as.numeric(f1[idx]); v2<-as.numeric(f2[idx]); if(is.na(v1)||is.na(v2)) NA else v1-v2 }
  
  dc_m <- gds(fmc,fmm,"cfi.scaled"); dt_m <- gds(fmc,fmm,"tli.scaled")
  dr_m <- -gds(fmc,fmm,"rmsea.scaled"); ds_m <- -gds(fmc,fmm,"srmr")
  dc_s <- gds(fmm,fms,"cfi.scaled"); dt_s <- gds(fmm,fms,"tli.scaled")
  dr_s <- -gds(fmm,fms,"rmsea.scaled"); ds_s <- -gds(fmm,fms,"srmr")
  dc_t <- gds(fms,fmt,"cfi.scaled"); dt_t <- gds(fms,fmt,"tli.scaled")
  dr_t <- -gds(fms,fmt,"rmsea.scaled"); ds_t <- -gds(fms,fmt,"srmr")
  
  if (any(is.na(c(dc_m,dt_m,dr_m,ds_m,dc_s,dt_s,dr_s,ds_s,dc_t,dt_t,dr_t,ds_t)))) return(NULL)
  
  th_c <- 0.010; th_t <- 0.010; th_r <- 0.015; th_s <- 0.030
  
  nc <- min(max(pmax(0,dc_m),pmax(0,dc_s),pmax(0,dc_t))/th_c, 1)
  nt <- min(max(pmax(0,dt_m),pmax(0,dt_s),pmax(0,dt_t))/th_t, 1)
  nr <- min(max(pmax(0,dr_m),pmax(0,dr_s),pmax(0,dr_t))/th_r, 1)
  ns <- min(max(pmax(0,ds_m),pmax(0,ds_s),pmax(0,ds_t))/th_s, 1)
  
  mi_score <- weights[1]*nc + weights[2]*nt + weights[3]*nr + weights[4]*ns
  
  is_inv <- (dc_m<=th_c && dt_m<=th_t && dr_m<=th_r && ds_m<=th_s &&
               dc_s<=th_c && dt_s<=th_t && dr_s<=th_r && ds_s<=th_s &&
               dc_t<=th_c && dt_t<=th_t && dr_t<=th_r && ds_t<=th_s)
  
  list(moderator=moderator, mi_score=mi_score, is_invariant=is_inv)
}

cat("Computing adaptive weights and MI assessment...\n")
calibrated_weights <- compute_adaptive_weights(data, moderators, cfa_model)

mi_per_moderator <- list()
for (mod in moderators) {
  mi_res <- test_mi_adaptive(data, mod, cfa_model, calibrated_weights)
  if (!is.null(mi_res)) mi_per_moderator[[mod]] <- mi_res
}

cat("MI Status:\n")
for (mod in names(mi_per_moderator)) {
  mr <- mi_per_moderator[[mod]]
  cat(sprintf("  %-25s: MI=%s (S_MI=%.3f)\n", mod,
              ifelse(mr$is_invariant, "OK ", "BAD"), mr$mi_score))
}
cat("\n")

build_sem_models <- function(outcome_var, stab_var) {
  m0 <- paste0('
    MEQ_Factor =~ MEQ1 + MEQ2 + MEQ8 + MEQ9 + MEQ10 + MEQ11 + MEQ15 + MEQ17 + MEQ18 + MEQ19
    ', outcome_var, ' ~ c*MEQ_Factor
    MEQ1 ~~ MEQ2
    MEQ2 ~~ MEQ10
    MEQ11 ~~ MEQ15
    MEQ17 ~~ MEQ18
    MEQ2 ~~ MEQ8
  ')
  m1 <- paste0('
    MEQ_Factor =~ MEQ1 + MEQ2 + MEQ8 + MEQ9 + MEQ10 + MEQ11 + MEQ15 + MEQ17 + MEQ18 + MEQ19
    ', stab_var, ' ~ a*MEQ_Factor
    ', outcome_var, ' ~ cprime*MEQ_Factor + b*', stab_var, '
    MEQ1 ~~ MEQ2
    MEQ2 ~~ MEQ10
    MEQ11 ~~ MEQ15
    MEQ17 ~~ MEQ18
    MEQ2 ~~ MEQ8
  ')
  list(without = m0, with = m1)
}

compute_stabilization <- function(data, moderator, models) {
  d <- data[!is.na(data[[moderator]]), ]
  d[[moderator]] <- factor(d[[moderator]])
  groups <- levels(d[[moderator]])
  if (length(groups) < 2) return(NULL)
  
  c_paths <- c(); cprime_paths <- c()
  group_ns <- c(); group_names <- c()
  
  for (g in groups) {
    gd <- d[d[[moderator]] == g, ]
    if (nrow(gd) < 30) next
    fit0 <- tryCatch(sem(models$without, data=gd, std.lv=TRUE, estimator="MLR", missing="fiml"), error=function(e) NULL)
    fit1 <- tryCatch(sem(models$with, data=gd, std.lv=TRUE, estimator="MLR", missing="fiml"), error=function(e) NULL)
    if (!is.null(fit0) && !is.null(fit1)) {
      s0 <- standardizedSolution(fit0); s1 <- standardizedSolution(fit1)
      cp <- s0[s0$label=="c","est.std"][1]; cpp <- s1[s1$label=="cprime","est.std"][1]
      if (!is.na(cp) && !is.na(cpp)) {
        c_paths <- c(c_paths, cp); cprime_paths <- c(cprime_paths, cpp)
        group_ns <- c(group_ns, nrow(gd)); group_names <- c(group_names, g)
      }
    }
  }
  if (length(c_paths) < 2) return(NULL)
  
  wmean <- function(x, w) weighted.mean(x, w)
  wsd   <- function(x, w) sqrt(sum(w * (x - wmean(x, w))^2) / sum(w))
  
  mean_c <- wmean(c_paths, group_ns); mean_cp <- wmean(cprime_paths, group_ns)
  sd_c <- wsd(c_paths, group_ns); sd_cp <- wsd(cprime_paths, group_ns)
  
  eps_mu <- 0.0001
  cv_before <- (sd_c / max(abs(mean_c), eps_mu)) * 100
  cv_after  <- (sd_cp / max(abs(mean_cp), eps_mu)) * 100
  
  eps_log <- 1e-8
  dl <- log(cv_before + eps_log) - log(cv_after + eps_log)
  pct_red <- 100 * (1 - exp(-dl))
  
  I_g <- sapply(seq_along(c_paths), function(i) {
    ifelse(abs(cprime_paths[i] - mean_cp) < abs(c_paths[i] - mean_c), 1, 0)
  })
  ocr <- mean(I_g)
  
  eps_os <- 1e-8
  delta_log_mu    <- log(abs(mean_cp) + eps_os) - log(abs(mean_c) + eps_os)
  delta_log_sigma <- log(sd_c + eps_os) - log(sd_cp + eps_os)
  os <- max(0, delta_log_mu) / (max(0, delta_log_mu) + max(0, delta_log_sigma) + eps_os)
  
  list(moderator=moderator, delta_log=dl, cv_before=cv_before, cv_after=cv_after,
       pct_reduction=pct_red, ocr=ocr, os=os, n_groups=length(c_paths), total_n=sum(group_ns))
}

run_svt <- function(data, candidate_name, outcome_name, models, moderators,
                    mi_per_moderator, verbose = TRUE) {
  label <- paste0(candidate_name, " on MEQ->", outcome_name)
  if (verbose) {
    cat("\nSVT:", label, "\n")
    cat(paste(rep("-", 60), collapse = ""), "\n")
  }
  
  results <- list()
  for (mod in moderators) {
    if (!(mod %in% names(data))) next
    r <- compute_stabilization(data, mod, models)
    if (!is.null(r)) {
      mi <- mi_per_moderator[[mod]]
      r$mi_score <- if (!is.null(mi)) mi$mi_score else NA
      r$is_invariant <- if (!is.null(mi)) mi$is_invariant else NA
      results[[mod]] <- r
      if (verbose) {
        cat(sprintf("  %-25s: MI=%s (%.3f), Dl=%6.3f, CV:%5.1f%%->%5.1f%% (%+.1f%%)\n",
                    mod,
                    ifelse(is.na(r$is_invariant), "???",
                           ifelse(r$is_invariant, "OK ", "BAD")),
                    ifelse(is.na(r$mi_score), 0, r$mi_score),
                    r$delta_log, r$cv_before, r$cv_after, r$pct_reduction))
      }
    }
  }
  
  if (length(results) == 0) {
    return(list(candidate=candidate_name, outcome=outcome_name,
                decision="No valid results", p_boot=NA, p_binom=NA,
                mean_dl=NA, d_cohen=NA, weighted_os=NA,
                n_positive=NA, n_total=NA, summary=NULL, detailed=list()))
  }
  
  delta_vec <- sapply(results, function(x) x$delta_log)
  wt_vec    <- sapply(results, function(x) x$total_n)
  ocr_vec   <- sapply(results, function(x) x$ocr)
  os_vec    <- sapply(results, function(x) x$os)
  mi_scores <- sapply(results, function(x) x$mi_score)
  is_inv    <- sapply(results, function(x) x$is_invariant)
  
  boot_mean <- function(d, idx) {
    w <- wt_vec[idx] / sum(wt_vec[idx])
    sum(d[idx] * w)
  }
  boot_res <- boot(delta_vec, boot_mean, R = 1000)
  mean_dl  <- boot_res$t0
  se_dl    <- sd(boot_res$t)
  z_stat   <- mean_dl / se_dl
  p_boot   <- pnorm(z_stat, lower.tail = FALSE)
  
  ci_obj <- tryCatch(boot.ci(boot_res, type = "perc"), error = function(e) NULL)
  ci_lo  <- if (!is.null(ci_obj)) ci_obj$percent[4] else NA
  ci_hi  <- if (!is.null(ci_obj)) ci_obj$percent[5] else NA
  
  wvar    <- sum(wt_vec * (delta_vec - mean_dl)^2) / sum(wt_vec)
  d_cohen <- if (wvar > 0) mean_dl / sqrt(wvar) else NA
  
  crit2_pass <- (delta_vec > 0) | (ocr_vec >= 0.5)
  S_obs   <- sum(crit2_pass)
  M_total <- length(crit2_pass)
  p_binom <- binom.test(S_obs, M_total, p = 0.5, alternative = "greater")$p.value
  
  weighted_os <- weighted.mean(os_vec, wt_vec)
  
  bootstrap_pass <- (p_boot < 0.05 && mean_dl > 0)
  binom_pass     <- (p_binom < 0.05)
  is_stabilizer  <- bootstrap_pass && binom_pass
  
  if (is_stabilizer) {
    if (weighted_os < 0.3) mech <- "Type A"
    else if (weighted_os > 0.7) mech <- "Type B"
    else mech <- "Type AB"
  } else {
    mech <- "None"
  }
  
  r_mi <- suppressWarnings(cor(mi_scores[!is.na(mi_scores)], delta_vec[!is.na(mi_scores)]))
  n_mi_bad <- sum(!is_inv, na.rm = TRUE)
  n_mi_good <- sum(is_inv, na.rm = TRUE)
  
  if (verbose) {
    cat(sprintf("\n  p_boot=%.4f %s | p_binom=%.4f %s | Dl=%.3f | d=%.3f | %d/%d pos | %s\n",
                p_boot, ifelse(bootstrap_pass, "[PASS]", "[FAIL]"),
                p_binom, ifelse(binom_pass, "[PASS]", "[FAIL]"),
                mean_dl, d_cohen, sum(delta_vec > 0), M_total,
                ifelse(is_stabilizer, paste("STABILIZER -", mech), "NOT A STABILIZER")))
    cat(sprintf("  MI: %d satisfied, %d violated | cor(S_MI, Dl) = %.3f\n",
                n_mi_good, n_mi_bad, ifelse(is.na(r_mi), 0, r_mi)))
  }
  
  summary_df <- data.frame(
    moderator = names(results), delta_log = delta_vec,
    pct_reduction = sapply(results, function(x) x$pct_reduction),
    cv_before = sapply(results, function(x) x$cv_before),
    cv_after = sapply(results, function(x) x$cv_after),
    ocr = ocr_vec, os = os_vec,
    mi_score = mi_scores, mi_satisfied = is_inv,
    n_groups = sapply(results, function(x) x$n_groups),
    total_n = wt_vec, stringsAsFactors = FALSE)
  
  list(candidate = candidate_name, outcome = outcome_name,
       decision = ifelse(is_stabilizer, "Stabilizer", "Not a stabilizer"),
       mechanism = mech, p_boot = p_boot, p_binom = p_binom,
       mean_dl = mean_dl, se_dl = se_dl, d_cohen = d_cohen,
       ci_lower = ci_lo, ci_upper = ci_hi, weighted_os = weighted_os,
       n_positive = sum(delta_vec > 0), n_total = M_total,
       n_mi_satisfied = n_mi_good, n_mi_violated = n_mi_bad,
       cor_mi_dl = r_mi,
       summary = summary_df, detailed = results)
}

candidates <- list(
  list(name = "BRIAN",  var = "BRIAN",       role = "Positive control"),
  list(name = "PSQI",   var = "PSQI_ctrl",   role = "Negative control"),
  list(name = "KIDMED", var = "KIDMED_ctrl", role = "Negative control")
)

focal_paths <- list(
  list(name = "YJIPAQ", var = "YJIPAQ",     note = "Physical activity (primary)"),
  list(name = "PSF12",  var = "PSF12_ctrl",  note = "Physical health (PSQI->PSF12 beta=-0.300)"),
  list(name = "KIDMED", var = "KIDMED_ctrl", note = "Diet quality (BRIAN->KIDMED=-0.174)")
)

all_results <- list()
result_table <- list()

for (fp in focal_paths) {
  cat("\n")
  cat(paste(rep("=", 70), collapse = ""), "\n")
  cat("FOCAL PATH: MEQ ->", fp$name, " [", fp$note, "]\n")
  cat(paste(rep("=", 70), collapse = ""), "\n")
  
  for (cand in candidates) {
    if (cand$var == fp$var) {
      cat(sprintf("\n  Skipping %s: candidate = outcome variable\n", cand$name))
      next
    }
    
    models <- build_sem_models(fp$var, cand$var)
    res <- run_svt(data, cand$name, fp$name, models, moderators, mi_per_moderator)
    
    key <- paste0(cand$name, "_", fp$name)
    all_results[[key]] <- res
    
    result_table[[length(result_table) + 1]] <- data.frame(
      focal_path = paste0("MEQ->", fp$name),
      candidate = cand$name,
      role = cand$role,
      decision = res$decision,
      mechanism = res$mechanism,
      p_boot = res$p_boot,
      p_binom = res$p_binom,
      mean_dl = res$mean_dl,
      d_cohen = res$d_cohen,
      weighted_os = res$weighted_os,
      n_positive = res$n_positive,
      n_total = res$n_total,
      n_mi_satisfied = res$n_mi_satisfied,
      n_mi_violated = res$n_mi_violated,
      cor_mi_dl = res$cor_mi_dl,
      ci_lower = res$ci_lower,
      ci_upper = res$ci_upper,
      stringsAsFactors = FALSE)
  }
}

result_df <- do.call(rbind, result_table)

cat("\n")
cat(paste(rep("=", 100), collapse = ""), "\n")
cat("COMPREHENSIVE COMPARISON TABLE\n")
cat(paste(rep("=", 100), collapse = ""), "\n\n")

cat(sprintf("%-14s | %-8s | %8s | %8s | %8s | %8s | %5s | %5s | %s\n",
            "Path", "Candidat", "p_boot", "p_binom", "Mean Dl", "Cohen d", "Pos", "MI_OK", "Decision"))
cat(paste(rep("-", 100), collapse = ""), "\n")

for (i in seq_len(nrow(result_df))) {
  r <- result_df[i, ]
  cat(sprintf("%-14s | %-8s | %8.4f | %8.4f | %8.3f | %8.3f | %3d/%1d | %2d/%2d | %s\n",
              r$focal_path, r$candidate, r$p_boot, r$p_binom,
              r$mean_dl, r$d_cohen, r$n_positive, r$n_total,
              r$n_mi_satisfied, r$n_mi_satisfied + r$n_mi_violated,
              r$decision))
}

cat("\n\nExpected pattern:\n")
cat("  BRIAN should be the ONLY candidate classified as Stabilizer\n")
cat("  on MEQ->YJIPAQ. PSQI and KIDMED should fail on ALL focal paths,\n")
cat("  even on paths where they are structurally connected (e.g. PSQI->PSF12).\n")
cat("  This confirms: structural connection alone does not make a stabilizer.\n")
cat("  C1 (MI-artifact correlation) is specifically required.\n")

saveRDS(list(
  results = all_results,
  comparison = result_df,
  mi_per_moderator = mi_per_moderator,
  calibrated_weights = calibrated_weights,
  metadata = list(
    dataset = "analizliksonAVE.xlsx", n = nrow(data),
    cfa_model = "5-cov", focal_paths = focal_paths,
    candidates = candidates, date = Sys.time())
), file.path(rds_path, "SVT_negative_controls.rds"))

write.csv(result_df, file.path(output_path, "SVT_negative_controls_comparison.csv"), row.names = FALSE)

per_mod_list <- list()
for (key in names(all_results)) {
  r <- all_results[[key]]
  if (!is.null(r$summary)) {
    s <- r$summary
    s$candidate <- r$candidate
    s$focal_outcome <- r$outcome
    per_mod_list[[length(per_mod_list) + 1]] <- s
  }
}
per_mod_df <- do.call(rbind, per_mod_list)
write.csv(per_mod_df, file.path(output_path, "SVT_negative_controls_per_moderator.csv"), row.names = FALSE)

cat("\nSaved:\n")
cat("  RDS/SVT_negative_controls.rds\n")
cat("  Outputs/SVT_negative_controls_comparison.csv\n")
cat("  Outputs/SVT_negative_controls_per_moderator.csv\n")

cat("\nGenerating figures...\n")

col_brian  <- "#27ae60"
col_psqi   <- "#e74c3c"
col_kidmed <- "#3498db"

p1_data <- result_df %>%
  mutate(candidate = factor(candidate, levels = c("BRIAN","PSQI","KIDMED")),
         focal_path = factor(focal_path))

p1 <- ggplot(p1_data, aes(x = candidate, y = mean_dl, fill = candidate)) +
  geom_col(width = 0.6, alpha = 0.85) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.15, linewidth = 0.7) +
  geom_hline(yintercept = 0, color = "gray30") +
  geom_text(aes(label = sprintf("p=%.3f", p_boot)),
            vjust = ifelse(p1_data$mean_dl >= 0, -1.2, 1.5), size = 2.8) +
  facet_wrap(~focal_path, scales = "free_y") +
  scale_fill_manual(values = c("BRIAN" = col_brian, "PSQI" = col_psqi, "KIDMED" = col_kidmed)) +
  labs(x = "", y = expression(paste("Weighted Mean ", Delta, ell)),
       title = "Stabilization Effect Across Focal Paths and Candidates",
       subtitle = "Error bars = 95% bootstrap CI. Only BRIAN should show significant positive effect.") +
  theme_minimal() +
  theme(plot.title = element_text(size = 13, face = "bold"),
        plot.subtitle = element_text(size = 10, color = "gray40"),
        legend.position = "none", strip.text = element_text(face = "bold"))

p2_data <- data.frame(
  candidate = character(), criterion = character(),
  focal_path = character(), p_value = numeric(),
  pass = logical(), stringsAsFactors = FALSE)

for (i in seq_len(nrow(result_df))) {
  r <- result_df[i, ]
  p2_data <- rbind(p2_data,
                   data.frame(candidate = r$candidate, criterion = "Bootstrap\n(magnitude)",
                              focal_path = r$focal_path, p_value = r$p_boot,
                              pass = (r$p_boot < 0.05 & r$mean_dl > 0), stringsAsFactors = FALSE),
                   data.frame(candidate = r$candidate, criterion = "Binomial\n(consistency)",
                              focal_path = r$focal_path, p_value = r$p_binom,
                              pass = (r$p_binom < 0.05), stringsAsFactors = FALSE))
}

p2 <- ggplot(p2_data, aes(x = candidate, y = criterion, fill = pass)) +
  geom_tile(color = "white", linewidth = 1.5) +
  geom_text(aes(label = sprintf("%.4f", p_value)), size = 2.8, fontface = "bold") +
  facet_wrap(~focal_path) +
  scale_fill_manual(values = c("TRUE" = "#27ae60", "FALSE" = "#e74c3c"),
                    labels = c("FAIL", "PASS"), name = "") +
  labs(x = "", y = "", title = "Dual-Criterion Decision Matrix Across All Paths",
       subtitle = "Both criteria must PASS (green) for stabilizer classification") +
  theme_minimal() +
  theme(plot.title = element_text(size = 13, face = "bold"),
        plot.subtitle = element_text(size = 10, color = "gray40"),
        legend.position = "bottom", panel.grid = element_blank(),
        strip.text = element_text(face = "bold"))

yjipaq_mods <- per_mod_df %>%
  filter(focal_outcome == "YJIPAQ") %>%
  mutate(candidate = factor(candidate, levels = c("BRIAN","PSQI","KIDMED")))

p3 <- ggplot(yjipaq_mods, aes(x = reorder(moderator, delta_log), y = delta_log, fill = candidate)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7, alpha = 0.85) +
  geom_hline(yintercept = 0, color = "gray30") +
  coord_flip() +
  scale_fill_manual(values = c("BRIAN" = col_brian, "PSQI" = col_psqi, "KIDMED" = col_kidmed),
                    name = "Candidate") +
  labs(x = "", y = expression(paste(Delta, ell)),
       title = "Per-Moderator Stabilization on MEQ -> YJIPAQ",
       subtitle = "BRIAN shows consistent positive values; PSQI and KIDMED are mixed/negative") +
  theme_minimal() +
  theme(plot.title = element_text(size = 13, face = "bold"),
        plot.subtitle = element_text(size = 10, color = "gray40"),
        legend.position = "bottom")

p4 <- ggplot(yjipaq_mods, aes(x = candidate, y = pct_reduction, fill = candidate)) +
  geom_boxplot(alpha = 0.7, outlier.size = 1.5) +
  geom_hline(yintercept = 0, color = "gray30", linetype = "dashed") +
  scale_fill_manual(values = c("BRIAN" = col_brian, "PSQI" = col_psqi, "KIDMED" = col_kidmed)) +
  labs(x = "", y = "CV Reduction (%)",
       title = "Distribution of Stabilization (MEQ -> YJIPAQ)",
       subtitle = "BRIAN: consistently positive. Controls: centered around zero.") +
  theme_minimal() +
  theme(plot.title = element_text(size = 13, face = "bold"),
        plot.subtitle = element_text(size = 10, color = "gray40"),
        legend.position = "none")

ggsave(file.path(figure_path, "neg_ctrl_multipath_effect.png"), p1, width = 14, height = 6, dpi = 600)
ggsave(file.path(figure_path, "neg_ctrl_decision_matrix.png"), p2, width = 14, height = 6, dpi = 600)
ggsave(file.path(figure_path, "neg_ctrl_per_moderator.png"), p3, width = 12, height = 7, dpi = 600)
ggsave(file.path(figure_path, "neg_ctrl_boxplot.png"), p4, width = 8, height = 6, dpi = 600)

combined <- grid.arrange(p1, p2, p3, p4, ncol = 2, nrow = 2,
                         top = grid::textGrob("SVT Negative Control Analysis: Multi-Path Comparison",
                                              gp = grid::gpar(fontsize = 16, fontface = "bold")))
ggsave(file.path(figure_path, "neg_ctrl_combined_4panel.png"),
       combined, width = 18, height = 14, dpi = 600)

cat("Figures saved to Figures/\n")
cat("\nDone.\n")