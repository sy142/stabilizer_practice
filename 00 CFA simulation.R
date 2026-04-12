rm(list = ls())
suppressPackageStartupMessages({
  library(readxl)
  library(lavaan)
  library(MASS)
  library(dplyr)
  library(stringr)
  library(boot)
  library(foreach)
  library(doParallel)
  library(ggplot2)
  library(gridExtra)
})

set.seed(9186)

cat("CFA-BASED SVT MONTE CARLO SIMULATION\n")
cat("Full Pipeline: MI Cascade + Stabilization + Interaction Filter\n")
cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

base_path   <- "C:/Users/Salim/OneDrive/Desktop/Structural Equation Modeling Companion"
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

cat("PHASE 0: Extracting Population Parameters\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

data_real <- read_excel(file.path(data_path, "analizliksonAVE.xlsx"))
data_real <- data_real %>% mutate(across(where(is.character), turkce_ascii))
data_real$Z <- data_real$BRO_Toplam

x_items <- c("MEQ1","MEQ2","MEQ8","MEQ9","MEQ10","MEQ11","MEQ15","MEQ17","MEQ18","MEQ19")
x_labels <- paste0("x", 1:10)

for (i in seq_along(x_items)) data_real[[x_labels[i]]] <- data_real[[x_items[i]]]
data_real$Y <- data_real$YJIPAQ

cfa_model <- '
  X_Factor =~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10
  x1 ~~ x2
  x2 ~~ x5
  x6 ~~ x7
  x8 ~~ x9
  x2 ~~ x3
'

fit_pop <- cfa(cfa_model, data = data_real, std.lv = TRUE, estimator = "MLR")
pop_params <- list(
  loadings_std = standardizedSolution(fit_pop) %>% filter(op == "=~") %>% pull(est.std),
  item_labels = x_labels, n_items = 10
)

sem_model_pop <- '
  X_Factor =~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10
  Z ~ a*X_Factor
  Y ~ cprime*X_Factor + b*Z
  x1 ~~ x2
  x2 ~~ x5
  x6 ~~ x7
  x8 ~~ x9
  x2 ~~ x3
'

fit_sem_pop <- sem(sem_model_pop, data = data_real, std.lv = TRUE, estimator = "MLR", missing = "fiml")
std_sem <- standardizedSolution(fit_sem_pop)
pop_params$a_path <- std_sem[std_sem$label == "a", "est.std"]
pop_params$cprime_path <- std_sem[std_sem$label == "cprime", "est.std"]
pop_params$b_path <- std_sem[std_sem$label == "b", "est.std"]

cat("Population parameters:\n")
for (i in seq_along(x_labels)) cat(sprintf("  %s: %.3f\n", x_labels[i], pop_params$loadings_std[i]))
cat(sprintf("  a=%.3f  cprime=%.3f  b=%.3f  c=%.3f\n\n",
            pop_params$a_path, pop_params$cprime_path, pop_params$b_path,
            pop_params$cprime_path + pop_params$a_path * pop_params$b_path))

cat("PHASE 1: Simulation Design\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

N_REPS  <- 500
N_CORES <- 14
BOOT_R  <- 1000
SA <- 0.5
SB <- 0.4

MOD_CONFIGS <- list(
  M7  = list(structure = list(G1=2,G2=2,G3=2,G4=3,G5=3,G6=4,G7=5),
             label = "7mod(3x2,2x3,1x4,1x5)"),
  M10 = list(structure = list(G1=2,G2=2,G3=2,G4=2,G5=3,G6=3,G7=3,G8=4,G9=4,G10=5),
             label = "10mod(4x2,3x3,2x4,1x5)"),
  M12 = list(structure = list(G1=2,G2=2,G3=2,G4=2,G5=2,G6=3,G7=3,G8=3,G9=3,G10=4,G11=4,G12=5),
             label = "12mod(5x2,4x3,2x4,1x5)")
)

SAMPLE_SIZES <- c(500, 1000, 1500, 2000)
MI_LEVELS    <- c(0, 0.2, 0.3, 0.5, 0.7)
SCENARIOS    <- c("H0", "H1_TypeA", "H1_TypeB", "H1_TypeAB", "H1_Moderator")
INT_STRENGTHS <- c(0.05, 0.10, 0.15, 0.20, 0.25)

cat("Scenarios: H0, H1_TypeA, H1_TypeB, H1_TypeAB, H1_Moderator\n")
cat(sprintf("Moderator interaction strengths: %s\n", paste(INT_STRENGTHS, collapse = ", ")))
cat(sprintf("Configs: %s | N: %s | MI: %s\n",
            paste(names(MOD_CONFIGS), collapse = ","),
            paste(SAMPLE_SIZES, collapse = ","),
            paste(MI_LEVELS, collapse = ",")))
cat(sprintf("Reps=%d | Bootstrap=%d | Cores=%d\n", N_REPS, BOOT_R, N_CORES))
cat(sprintf("Phase A conditions: %d (H0+H1_Type*+Mod_025)\n",
            length(MOD_CONFIGS) * length(SAMPLE_SIZES) * length(MI_LEVELS) * length(SCENARIOS)))
cat(sprintf("Phase B conditions: %d (Mod graded 0.05-0.20)\n",
            length(MOD_CONFIGS) * length(SAMPLE_SIZES) * length(MI_LEVELS) * 4))
cat(sprintf("Total SVT runs: %d\n\n",
            length(MOD_CONFIGS) * length(SAMPLE_SIZES) * length(MI_LEVELS) * (5 + 4) * N_REPS))

cat("Decision rule:\n")
cat("  SVT: bootstrap p<0.05 AND mean_delta_log>0 AND (binom_delta p<0.05 OR binom_orient p<0.05)\n")
cat("  Interaction filter (H1_Moderator only): Residual-Centered Product Indicator\n")
cat("    Little, Bovaird & Widaman (2006). gamma p<0.05 => interaction detected => Not Stabilizer\n")
cat("  Final: SVT_pass AND NOT interaction_flag\n\n")

generate_sim_data <- function(n_total, mod_structure, mi_severity, scenario,
                              pop_loadings_std, x_labels, pop_params, seed,
                              sa, sb, int_strength = 0.25) {
  set.seed(seed)
  n_items <- length(x_labels)
  mod_names <- names(mod_structure)
  mod_ngroups <- as.integer(mod_structure)
  n_mods <- length(mod_structure)
  n_per_group <- max(n_total %/% sum(mod_ngroups), 30)
  all_data <- list()
  for (m in seq_len(n_mods)) {
    n_g <- mod_ngroups[m]
    for (g in seq_len(n_g)) {
      perturbation_sd <- mi_severity * 0.10
      if (perturbation_sd > 0) {
        delta <- rnorm(n_items, mean = 0, sd = perturbation_sd)
        direction <- (2 * g - n_g - 1) / max(n_g - 1, 1)
        delta <- abs(delta) * direction * 0.5
      } else {
        delta <- rep(0, n_items)
      }
      lambda_g <- pmin(pmax(pop_loadings_std + delta, 0.20), 0.95)
      eta <- rnorm(n_per_group)
      X <- matrix(NA, n_per_group, n_items)
      for (j in seq_len(n_items)) {
        rv <- max(1 - lambda_g[j]^2, 0.05)
        X[, j] <- lambda_g[j] * eta + rnorm(n_per_group, 0, sqrt(rv))
      }
      colnames(X) <- x_labels
      mi_artifact <- rowMeans(X[, 1:6]) - mean(pop_loadings_std[1:6]) * eta
      if (scenario == "H1_TypeA") {
        a_true <- pop_params$a_path; cprime_true <- pop_params$cprime_path; b_true <- pop_params$b_path
        Z <- a_true * eta + rnorm(n_per_group, 0, sqrt(max(1 - a_true^2, 0.1)))
        b_perturb <- rnorm(1, mean = 0, sd = mi_severity * sa)
        rv <- max(1 - cprime_true^2 - (b_true + b_perturb)^2, 0.1)
        Y <- cprime_true * eta + (b_true + b_perturb) * Z + rnorm(n_per_group, 0, sqrt(rv))
      } else if (scenario == "H1_TypeB") {
        a_true <- pop_params$a_path; cprime_true <- pop_params$cprime_path; b_true <- pop_params$b_path
        Z <- a_true * eta + sb * mi_artifact + rnorm(n_per_group, 0, sqrt(max(1 - a_true^2, 0.1)))
        rv <- max(1 - cprime_true^2 - b_true^2, 0.1)
        Y <- cprime_true * eta + b_true * Z + rnorm(n_per_group, 0, sqrt(rv))
      } else if (scenario == "H1_TypeAB") {
        a_true <- pop_params$a_path; cprime_true <- pop_params$cprime_path; b_true <- pop_params$b_path
        Z <- a_true * eta + sb * mi_artifact + rnorm(n_per_group, 0, sqrt(max(1 - a_true^2, 0.1)))
        b_perturb <- rnorm(1, mean = 0, sd = mi_severity * sa)
        rv <- max(1 - cprime_true^2 - (b_true + b_perturb)^2, 0.1)
        Y <- cprime_true * eta + (b_true + b_perturb) * Z + rnorm(n_per_group, 0, sqrt(rv))
      } else if (scenario == "H1_Moderator") {
        a_true <- pop_params$a_path; cprime_true <- pop_params$cprime_path; b_true <- pop_params$b_path
        Z <- a_true * eta + rnorm(n_per_group, 0, sqrt(max(1 - a_true^2, 0.1)))
        interaction_term <- int_strength * eta * Z
        rv <- max(1 - cprime_true^2 - b_true^2, 0.1)
        Y <- cprime_true * eta + b_true * Z + interaction_term + rnorm(n_per_group, 0, sqrt(rv))
      } else {
        Z <- rnorm(n_per_group)
        c_true <- pop_params$cprime_path + pop_params$a_path * pop_params$b_path
        Y <- c_true * eta + rnorm(n_per_group, 0, sqrt(max(1 - c_true^2, 0.1)))
      }
      gdf <- as.data.frame(X); gdf$Z <- Z; gdf$Y <- Y
      for (mm in seq_len(n_mods)) {
        gdf[[mod_names[mm]]] <- if (mm == m) g else sample(seq_len(mod_ngroups[mm]), n_per_group, replace = TRUE)
      }
      all_data[[length(all_data) + 1]] <- gdf
    }
  }
  sim_data <- do.call(rbind, all_data)
  for (mm in mod_names) sim_data[[mm]] <- factor(sim_data[[mm]])
  sim_data
}

run_svt_full <- function(dat, mod_names, cfa_model, x_labels, boot_r) {
  gd_abs <- function(f1, f2, idx) {
    v1 <- as.numeric(f1[idx]); v2 <- as.numeric(f2[idx])
    if (is.na(v1) || is.na(v2)) NA else abs(v1 - v2)
  }
  gd_s <- function(f1, f2, idx) {
    v1 <- as.numeric(f1[idx]); v2 <- as.numeric(f2[idx])
    if (is.na(v1) || is.na(v2)) NA else v1 - v2
  }
  all_dm <- list(); all_ds <- list(); all_dt <- list()
  all_inv <- c(); mi_det <- list(); valid_mods <- c()
  for (mod in mod_names) {
    if (!(mod %in% names(dat))) next
    d <- dat[!is.na(dat[[mod]]), ]; d[[mod]] <- factor(d[[mod]])
    if (length(levels(d[[mod]])) < 2) next
    fc <- tryCatch(cfa(cfa_model, data=d, group=mod, std.lv=TRUE, estimator="MLR"), error=function(e) NULL)
    if (is.null(fc)) next
    fm <- tryCatch(cfa(cfa_model, data=d, group=mod, group.equal="loadings", std.lv=TRUE, estimator="MLR"), error=function(e) NULL)
    if (is.null(fm)) next
    fs <- tryCatch(cfa(cfa_model, data=d, group=mod, group.equal=c("loadings","intercepts"), std.lv=TRUE, estimator="MLR"), error=function(e) NULL)
    if (is.null(fs)) next
    ft <- tryCatch(cfa(cfa_model, data=d, group=mod, group.equal=c("loadings","intercepts","residuals"), std.lv=TRUE, estimator="MLR"), error=function(e) NULL)
    if (is.null(ft)) next
    fmc <- fitmeasures(fc); fmm <- fitmeasures(fm); fms <- fitmeasures(fs); fmt <- fitmeasures(ft)
    dm <- c(CFI=gd_abs(fmc,fmm,"cfi.scaled"), TLI=gd_abs(fmc,fmm,"tli.scaled"),
            RMSEA=gd_abs(fmc,fmm,"rmsea.scaled"), SRMR=gd_abs(fmc,fmm,"srmr"))
    ds <- c(CFI=gd_abs(fmm,fms,"cfi.scaled"), TLI=gd_abs(fmm,fms,"tli.scaled"),
            RMSEA=gd_abs(fmm,fms,"rmsea.scaled"), SRMR=gd_abs(fmm,fms,"srmr"))
    dt <- c(CFI=gd_abs(fms,fmt,"cfi.scaled"), TLI=gd_abs(fms,fmt,"tli.scaled"),
            RMSEA=gd_abs(fms,fmt,"rmsea.scaled"), SRMR=gd_abs(fms,fmt,"srmr"))
    dcm <- gd_s(fmc,fmm,"cfi.scaled"); dtm <- gd_s(fmc,fmm,"tli.scaled")
    drm <- -gd_s(fmc,fmm,"rmsea.scaled"); dsm <- -gd_s(fmc,fmm,"srmr")
    dcs <- gd_s(fmm,fms,"cfi.scaled"); dts <- gd_s(fmm,fms,"tli.scaled")
    drs <- -gd_s(fmm,fms,"rmsea.scaled"); dss <- -gd_s(fmm,fms,"srmr")
    dct <- gd_s(fms,fmt,"cfi.scaled"); dtt <- gd_s(fms,fmt,"tli.scaled")
    drt <- -gd_s(fms,fmt,"rmsea.scaled"); dst <- -gd_s(fms,fmt,"srmr")
    if (any(is.na(c(dm, ds, dt)))) next
    all_dm[[mod]] <- dm; all_ds[[mod]] <- ds; all_dt[[mod]] <- dt
    tc <- 0.010; tt <- 0.010; tr <- 0.015; ts <- 0.030
    inv <- (dcm<=tc && dtm<=tt && drm<=tr && dsm<=ts &&
              dcs<=tc && dts<=tt && drs<=tr && dss<=ts &&
              dct<=tc && dtt<=tt && drt<=tr && dst<=ts)
    all_inv <- c(all_inv, inv); names(all_inv)[length(all_inv)] <- mod
    mi_det[[mod]] <- list(is_invariant=inv,
                          metric=c(dCFI=dcm, dTLI=dtm, dRMSEA=drm, dSRMR=dsm),
                          scalar=c(dCFI=dcs, dTLI=dts, dRMSEA=drs, dSRMR=dss),
                          strict=c(dCFI=dct, dTLI=dtt, dRMSEA=drt, dSRMR=dst),
                          fit_config=as.numeric(fmc[c("cfi.scaled","tli.scaled","rmsea.scaled","srmr")]),
                          fit_metric=as.numeric(fmm[c("cfi.scaled","tli.scaled","rmsea.scaled","srmr")]),
                          fit_scalar=as.numeric(fms[c("cfi.scaled","tli.scaled","rmsea.scaled","srmr")]),
                          fit_strict=as.numeric(fmt[c("cfi.scaled","tli.scaled","rmsea.scaled","srmr")]))
    valid_mods <- c(valid_mods, mod)
  }
  if (length(all_dm) >= 3) {
    dm_mat <- do.call(rbind, all_dm); ds_mat <- do.call(rbind, all_ds); dt_mat <- do.call(rbind, all_dt)
    dmax <- pmax(dm_mat, ds_mat, dt_mat); cormat <- cor(dmax)
    nm <- cbind(CFI=pmin(dmax[,"CFI"]/0.01,1), TLI=pmin(dmax[,"TLI"]/0.01,1),
                RMSEA=pmin(dmax[,"RMSEA"]/0.015,1), SRMR=pmin(dmax[,"SRMR"]/0.03,1))
    cv <- apply(nm, 2, function(x) sd(x)/(mean(x)+0.001))
    mc <- apply(abs(cormat), 1, function(x) mean(x[x<1])); rp <- 1/(1+mc)
    dp <- numeric(4); names(dp) <- c("CFI","TLI","RMSEA","SRMR")
    for (idx in names(dp)) {
      if (sum(all_inv)>0 && sum(!all_inv)>0) dp[idx] <- max(0, mean(nm[!all_inv,idx]) - mean(nm[all_inv,idx]))
      else dp[idx] <- sd(nm[,idx])
    }
    vs <- 1/(1+cv); rw <- dp * rp * vs; wt <- as.numeric(rw / (sum(rw) + 1e-8))
    smi_per_mod <- numeric(length(valid_mods)); names(smi_per_mod) <- valid_mods
    for (i in seq_along(valid_mods)) smi_per_mod[i] <- sum(wt * nm[i,])
    mean_smi <- mean(smi_per_mod); adaptive_weights <- wt
    names(adaptive_weights) <- c("CFI","TLI","RMSEA","SRMR")
  } else if (length(all_dm) >= 1) {
    wt <- c(0.4, 0.1, 0.3, 0.2); smi_per_mod <- rep(0.5, length(valid_mods))
    names(smi_per_mod) <- valid_mods; mean_smi <- 0.5
    adaptive_weights <- wt; names(adaptive_weights) <- c("CFI","TLI","RMSEA","SRMR")
  } else {
    return(list(
      step1=list(n_valid_mods=0, n_mi_sat=NA, n_mi_viol=NA, mean_smi=NA, adaptive_weights=rep(NA,4)),
      step2=list(mean_dl=NA, sd_dl=NA, n_stab_mods=NA),
      step3=list(p_boot=NA, p_binom_delta=NA, p_binom_orient=NA, boot_pass=NA,
                 binom_delta_pass=NA, binom_orient_pass=NA, decision="Error"),
      summary=list(decision="Error", mean_dl=NA, d_cohen=NA, n_delta_pos=NA,
                   n_orient_good=NA, n_mods_valid=0)))
  }
  model_without <- paste(c(
    paste0("X_Factor =~ ", paste(x_labels, collapse = " + ")),
    "Y ~ c*X_Factor",
    "x1 ~~ x2", "x2 ~~ x5", "x6 ~~ x7", "x8 ~~ x9", "x2 ~~ x3"), collapse = "\n")
  model_with <- paste(c(
    paste0("X_Factor =~ ", paste(x_labels, collapse = " + ")),
    "Z ~ a*X_Factor", "Y ~ cprime*X_Factor + b*Z",
    "x1 ~~ x2", "x2 ~~ x5", "x6 ~~ x7", "x8 ~~ x9", "x2 ~~ x3"), collapse = "\n")
  stab <- list(); stab_detail <- list()
  for (mod in valid_mods) {
    d <- dat[!is.na(dat[[mod]]), ]; d[[mod]] <- factor(d[[mod]])
    grps <- levels(d[[mod]]); if (length(grps) < 2) next
    cp <- c(); cpp <- c(); gn <- c(); gnames <- c()
    for (g in grps) {
      gd <- d[d[[mod]]==g, ]; if (nrow(gd) < 30) next
      f0 <- tryCatch(sem(model_without, data=gd, std.lv=TRUE, estimator="MLR"), error=function(e) NULL)
      f1 <- tryCatch(sem(model_with, data=gd, std.lv=TRUE, estimator="MLR"), error=function(e) NULL)
      if (!is.null(f0) && !is.null(f1)) {
        s0 <- standardizedSolution(f0); s1 <- standardizedSolution(f1)
        c1 <- s0[s0$label=="c","est.std"][1]; c2 <- s1[s1$label=="cprime","est.std"][1]
        if (!is.na(c1) && !is.na(c2)) {
          cp <- c(cp, c1); cpp <- c(cpp, c2); gn <- c(gn, nrow(gd)); gnames <- c(gnames, g)
        }
      }
    }
    if (length(cp) < 2) next
    wmn <- function(x, w) weighted.mean(x, w)
    wsd <- function(x, w) sqrt(sum(w * (x - wmn(x,w))^2) / sum(w))
    mc_val <- wmn(cp, gn); mcp <- wmn(cpp, gn)
    sc <- wsd(cp, gn); scp <- wsd(cpp, gn)
    em <- 0.0001; el <- 1e-8
    cv0 <- (sc / max(abs(mc_val), em)) * 100; cv1 <- (scp / max(abs(mcp), em)) * 100
    dl <- log(cv0 + el) - log(cv1 + el); pct_red <- 100 * (1 - exp(-dl))
    ig <- sapply(seq_along(cp), function(i) ifelse(abs(cpp[i]-mcp) < abs(cp[i]-mc_val), 1, 0))
    ocr <- mean(ig); d_k <- abs(cp - mc_val) - abs(cpp - mcp)
    d_w <- wmn(d_k, gn); d_w_sd <- wsd(d_k, gn)
    stab[[mod]] <- list(delta_log=dl, ocr=ocr, total_n=sum(gn),
                        cv_without=cv0, cv_with=cv1, pct_reduction=pct_red, d_w=d_w, d_w_sd=d_w_sd)
    stab_detail[[mod]] <- list(n_groups=length(cp), group_names=gnames, group_ns=gn,
                               c_paths=cp, cprime_paths=cpp, mean_c=mc_val, mean_cprime=mcp, sd_c=sc, sd_cprime=scp,
                               I_g=ig, d_k=d_k,
                               smi=ifelse(mod %in% names(smi_per_mod), smi_per_mod[mod], NA),
                               is_invariant=ifelse(mod %in% names(all_inv), all_inv[mod], NA))
  }
  dvec <- sapply(stab, function(x) x$delta_log); dvec <- dvec[!is.na(dvec)]
  if (length(dvec) < 2) {
    return(list(
      step1=list(n_valid_mods=length(valid_mods), n_mi_sat=sum(all_inv), n_mi_viol=sum(!all_inv),
                 mean_smi=mean_smi, adaptive_weights=adaptive_weights, mi_detail=mi_det),
      step2=list(mean_dl=NA, sd_dl=NA, n_stab_mods=length(dvec),
                 moderator_results=stab, moderator_detail=stab_detail),
      step3=list(p_boot=NA, p_binom_delta=NA, p_binom_orient=NA, boot_pass=NA,
                 binom_delta_pass=NA, binom_orient_pass=NA, decision="Error"),
      summary=list(decision="Error", mean_dl=NA, d_cohen=NA,
                   n_delta_pos=NA, n_orient_good=NA, n_mods_valid=length(dvec))))
  }
  wts <- sapply(names(dvec), function(m) stab[[m]]$total_n)
  ocrs <- sapply(names(dvec), function(m) stab[[m]]$ocr)
  mdl <- weighted.mean(dvec, wts)
  sdl <- sqrt(sum(wts * (dvec - mdl)^2) / sum(wts))
  bfn <- function(d, idx) { w <- wts[idx]/sum(wts[idx]); sum(d[idx]*w) }
  br <- tryCatch(boot(dvec, bfn, R = boot_r), error = function(e) NULL)
  if (!is.null(br)) {
    se_boot <- sd(br$t)
    p_boot <- if (se_boot > 0) pnorm(br$t0 / se_boot, lower.tail = FALSE) else ifelse(br$t0 > 0, 0, 1)
    ci_boot <- tryCatch(boot.ci(br, type = "perc", conf = 0.95), error = function(e) NULL)
    if (!is.null(ci_boot)) { ci_lower <- ci_boot$percent[4]; ci_upper <- ci_boot$percent[5] }
    else { ci_lower <- NA; ci_upper <- NA }
  } else { p_boot <- NA; se_boot <- NA; ci_lower <- NA; ci_upper <- NA }
  d_cohen <- if (sdl > 0) mdl / sdl else NA
  n_delta_pos <- sum(dvec > 0); n_orient_good <- sum(ocrs >= 0.5); n_total <- length(dvec)
  p_bd <- binom.test(n_delta_pos, n_total, p = 0.5, alternative = "greater")$p.value
  p_bo <- binom.test(n_orient_good, n_total, p = 0.5, alternative = "greater")$p.value
  boot_pass <- (!is.na(p_boot) && p_boot < 0.05 && mdl > 0)
  bd_pass <- (p_bd < 0.05); bo_pass <- (p_bo < 0.05)
  decision <- ifelse(boot_pass && (bd_pass || bo_pass), "Stabilizer", "Not")
  os_values <- ocrs[names(dvec)]; mean_os <- mean(os_values, na.rm = TRUE)
  if (decision == "Stabilizer") {
    if (mean_os < 0.3) mech_class <- "TypeA"
    else if (mean_os > 0.7) mech_class <- "TypeB"
    else mech_class <- "TypeAB"
  } else { mech_class <- "None" }
  list(
    step1=list(n_valid_mods=length(valid_mods), n_mi_sat=sum(all_inv), n_mi_viol=sum(!all_inv),
               mean_smi=mean_smi, smi_per_mod=smi_per_mod, adaptive_weights=adaptive_weights, mi_detail=mi_det),
    step2=list(mean_dl=mdl, sd_dl=sdl, se_boot=se_boot, ci_lower=ci_lower, ci_upper=ci_upper,
               d_cohen=d_cohen, pct_reduction_mean=mean(sapply(stab, function(x) x$pct_reduction)),
               cv_without_mean=mean(sapply(stab, function(x) x$cv_without)),
               cv_with_mean=mean(sapply(stab, function(x) x$cv_with)),
               n_stab_mods=length(dvec), delta_vec=dvec, ocr_vec=ocrs, weight_vec=wts,
               moderator_results=stab, moderator_detail=stab_detail),
    step3=list(p_boot=p_boot, p_binom_delta=p_bd, p_binom_orient=p_bo,
               boot_pass=boot_pass, binom_delta_pass=bd_pass, binom_orient_pass=bo_pass,
               decision=decision, mechanism_class=mech_class, mean_orientation_share=mean_os),
    summary=list(decision=decision, mechanism=mech_class, mean_dl=mdl, d_cohen=d_cohen,
                 p_boot=p_boot, p_binom_delta=p_bd, p_binom_orient=p_bo,
                 n_delta_pos=n_delta_pos, n_orient_good=n_orient_good, n_mods_valid=n_total,
                 mean_smi=mean_smi, n_mi_sat=sum(all_inv), n_mi_viol=sum(!all_inv)))
}

pi_items <- c("x1", "x3", "x5", "x7", "x9")
pi_names <- paste0("pi_", pi_items, "Z")

sem_interaction_model <- paste(c(
  paste0("X_Factor =~ ", paste(x_labels, collapse = " + ")),
  paste0("XZ =~ ", paste(pi_names, collapse = " + ")),
  "Z ~ a*X_Factor",
  "Y ~ cprime*X_Factor + b*Z + gamma*XZ",
  "X_Factor ~~ 0*XZ",
  "x1 ~~ x2", "x2 ~~ x5", "x6 ~~ x7", "x8 ~~ x9", "x2 ~~ x3"
), collapse = "\n")

interaction_test_cfa <- function(dat, x_labels, pi_items, pi_names, sem_interaction_model) {
  dat_pi <- dat
  for (k in seq_along(pi_items)) {
    xi <- pi_items[k]
    raw_prod <- dat_pi[[xi]] * dat_pi$Z
    rc_fit <- lm(raw_prod ~ dat_pi[[xi]] + dat_pi$Z)
    dat_pi[[pi_names[k]]] <- residuals(rc_fit)
  }
  fit_int <- tryCatch(sem(sem_interaction_model, data = dat_pi, std.lv = TRUE, estimator = "MLR"), error = function(e) NULL)
  if (is.null(fit_int)) return(list(flag = FALSE, p = NA, gamma_est = NA, gamma_std = NA))
  if (!tryCatch(lavInspect(fit_int, "converged"), error = function(e) FALSE))
    return(list(flag = FALSE, p = NA, gamma_est = NA, gamma_std = NA))
  pe <- tryCatch(parameterEstimates(fit_int), error = function(e) NULL)
  if (is.null(pe)) return(list(flag = FALSE, p = NA, gamma_est = NA, gamma_std = NA))
  g_row <- pe[pe$label == "gamma", ]
  if (nrow(g_row) == 0 || is.na(g_row$pvalue[1]))
    return(list(flag = FALSE, p = NA, gamma_est = NA, gamma_std = NA))
  ss <- tryCatch(standardizedSolution(fit_int), error = function(e) NULL)
  gamma_std <- NA
  if (!is.null(ss)) { gs <- ss[ss$label == "gamma", ]; if (nrow(gs) > 0) gamma_std <- gs$est.std[1] }
  list(flag = (!is.na(g_row$pvalue[1]) && g_row$pvalue[1] < 0.05),
       p = g_row$pvalue[1], gamma_est = g_row$est[1], gamma_std = gamma_std)
}

cat("PHASE 2A: Main Simulation (H0 + H1_Type* + H1_Moderator gamma=0.25)\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

conditions_a <- expand.grid(
  config = names(MOD_CONFIGS), sample_size = SAMPLE_SIZES,
  mi_severity = MI_LEVELS, scenario = SCENARIOS, stringsAsFactors = FALSE
)

cl <- makeCluster(N_CORES)
registerDoParallel(cl)
clusterExport(cl, c("pop_params", "x_labels", "cfa_model",
                    "generate_sim_data", "run_svt_full",
                    "interaction_test_cfa", "sem_interaction_model",
                    "pi_items", "pi_names",
                    "MOD_CONFIGS", "BOOT_R", "SA", "SB"))
clusterEvalQ(cl, { suppressPackageStartupMessages({
  library(lavaan); library(MASS); library(dplyr); library(boot) }) })

t_start <- Sys.time()
all_sim_results <- list()
summary_rows <- list()

for (cond_idx in seq_len(nrow(conditions_a))) {
  cond <- conditions_a[cond_idx, ]
  cfg <- cond$config; n_sam <- cond$sample_size; mi_sev <- cond$mi_severity; scen <- cond$scenario
  mc <- MOD_CONFIGS[[cfg]]; mod_struct <- mc$structure; mod_names_local <- names(mod_struct)
  if (scen == "H1_TypeA") { sa_use <- SA; sb_use <- 0 }
  else if (scen == "H1_TypeB") { sa_use <- 0; sb_use <- SB }
  else if (scen == "H1_TypeAB") { sa_use <- SA; sb_use <- SB }
  else { sa_use <- 0; sb_use <- 0 }
  is_moderator <- (scen == "H1_Moderator")
  clusterExport(cl, c("mod_struct", "mod_names_local", "n_sam", "mi_sev", "scen",
                      "sa_use", "sb_use", "is_moderator"), envir = environment())
  cond_results <- foreach(rep = 1:N_REPS, .packages = c("lavaan","MASS","dplyr","boot"),
                          .errorhandling = "pass") %dopar% {
                            seed_val <- 9186 + cond_idx * 10000 + rep
                            sim_dat <- generate_sim_data(n_total = n_sam, mod_structure = mod_struct,
                                                         mi_severity = mi_sev, scenario = scen, pop_loadings_std = pop_params$loadings_std,
                                                         x_labels = x_labels, pop_params = pop_params, seed = seed_val,
                                                         sa = sa_use, sb = sb_use, int_strength = 0.25)
                            res <- tryCatch(run_svt_full(sim_dat, mod_names_local, cfa_model, x_labels, BOOT_R),
                                            error = function(e) list(
                                              step1=list(n_valid_mods=NA, n_mi_sat=NA, n_mi_viol=NA, mean_smi=NA),
                                              step2=list(mean_dl=NA, d_cohen=NA),
                                              step3=list(p_boot=NA, p_binom_delta=NA, p_binom_orient=NA, boot_pass=NA,
                                                         binom_delta_pass=NA, binom_orient_pass=NA, decision="Error"),
                                              summary=list(decision="Error", mean_dl=NA, d_cohen=NA, n_delta_pos=NA,
                                                           n_orient_good=NA, n_mods_valid=NA, mean_smi=NA, n_mi_sat=NA, n_mi_viol=NA)))
                            if (is_moderator && res$summary$decision != "Error") {
                              it <- tryCatch(interaction_test_cfa(sim_dat, x_labels, pi_items, pi_names, sem_interaction_model),
                                             error = function(e) list(flag=FALSE, p=NA, gamma_est=NA, gamma_std=NA))
                              int_flag <- ifelse(is.na(it$flag), FALSE, it$flag)
                              svt_pass <- (res$summary$decision == "Stabilizer")
                              final_decision <- ifelse(svt_pass && !int_flag, "Stabilizer", "Not")
                              res$step3$interaction_flag <- int_flag
                              res$step3$p_interaction <- it$p
                              res$step3$gamma_est <- it$gamma_est
                              res$step3$gamma_std <- it$gamma_std
                              res$step3$decision <- final_decision
                              res$summary$decision <- final_decision
                            }
                            res$rep <- rep; res$config <- cfg; res$scenario <- scen
                            res$mi_severity <- mi_sev; res$sample_size <- n_sam
                            if (is_moderator) res$int_strength <- 0.25
                            res
                          }
  all_sim_results <- c(all_sim_results, cond_results)
  decs <- sapply(cond_results, function(x) x$summary$decision)
  n_stab <- sum(decs == "Stabilizer", na.rm = TRUE)
  n_valid <- sum(decs != "Error", na.rm = TRUE)
  n_err <- sum(decs == "Error", na.rm = TRUE)
  rate <- n_stab / max(n_valid, 1)
  dls <- sapply(cond_results, function(x) x$step2$mean_dl)
  bp <- sapply(cond_results, function(x) x$step3$boot_pass)
  bd <- sapply(cond_results, function(x) x$step3$binom_delta_pass)
  bo <- sapply(cond_results, function(x) x$step3$binom_orient_pass)
  ndp <- sapply(cond_results, function(x) x$summary$n_delta_pos)
  nsat <- sapply(cond_results, function(x) x$summary$n_mi_sat)
  nviol <- sapply(cond_results, function(x) x$summary$n_mi_viol)
  smi_vals <- sapply(cond_results, function(x) x$summary$mean_smi)
  dc <- sapply(cond_results, function(x) x$step2$d_cohen)
  metric <- ifelse(scen == "H0", "TypeI", ifelse(scen == "H1_Moderator", "FPR", "Power"))
  cat(sprintf("%-4s | %-12s | %.1f | %5d | %s=%.3f | mdl=%.3f\n",
              cfg, scen, mi_sev, n_sam, metric, rate, mean(dls, na.rm=TRUE)))
  summary_rows[[length(summary_rows)+1]] <- data.frame(
    config=cfg, scenario=scen, mi_severity=mi_sev, sample_size=n_sam,
    int_strength=ifelse(is_moderator, 0.25, NA),
    rate=rate, metric=metric, n_valid=n_valid, n_error=n_err,
    mean_dl=mean(dls, na.rm=TRUE), sd_dl=sd(dls, na.rm=TRUE),
    mean_d_cohen=mean(dc, na.rm=TRUE),
    boot_rate=mean(bp, na.rm=TRUE), bdelta_rate=mean(bd, na.rm=TRUE), borient_rate=mean(bo, na.rm=TRUE),
    mean_n_delta_pos=mean(ndp, na.rm=TRUE),
    mean_n_orient_good=mean(sapply(cond_results, function(x) x$summary$n_orient_good), na.rm=TRUE),
    mean_n_mods_valid=mean(sapply(cond_results, function(x) x$summary$n_mods_valid), na.rm=TRUE),
    mean_mi_sat=mean(nsat, na.rm=TRUE), mean_mi_viol=mean(nviol, na.rm=TRUE),
    mean_smi=mean(smi_vals, na.rm=TRUE), stringsAsFactors=FALSE)
  if (cond_idx %% 10 == 0) {
    cat(sprintf("  >>> %d/%d completed\n", cond_idx, nrow(conditions_a)))
    saveRDS(list(summary=do.call(rbind, summary_rows), n_completed=cond_idx, n_total=nrow(conditions_a)),
            file.path(rds_path, "CFA_full_checkpoint.rds"))
  }
}

cat(sprintf("\nPhase 2A complete: %d conditions\n\n", nrow(conditions_a)))

cat("PHASE 2B: Graded Moderator (gamma = 0.05, 0.10, 0.15, 0.20)\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

INT_GRADED <- c(0.05, 0.10, 0.15, 0.20)
conditions_b <- expand.grid(
  config = names(MOD_CONFIGS), sample_size = SAMPLE_SIZES,
  mi_severity = MI_LEVELS, int_str = INT_GRADED, stringsAsFactors = FALSE
)

clusterExport(cl, c("INT_GRADED"), envir = environment())

for (ci in seq_len(nrow(conditions_b))) {
  mc <- conditions_b[ci, ]
  cfg <- mc$config; n_sam <- mc$sample_size; mi_sev <- mc$mi_severity; int_str <- mc$int_str
  mod_struct <- MOD_CONFIGS[[cfg]]$structure; mod_names_local <- names(mod_struct)
  clusterExport(cl, c("mod_struct", "mod_names_local", "n_sam", "mi_sev", "int_str", "ci"),
                envir = environment())
  cond_results <- foreach(rep = 1:N_REPS, .packages = c("lavaan","MASS","dplyr","boot"),
                          .errorhandling = "pass") %dopar% {
                            seed_val <- 50000 + ci * 10000 + rep
                            sim_dat <- generate_sim_data(n_total = n_sam, mod_structure = mod_struct,
                                                         mi_severity = mi_sev, scenario = "H1_Moderator",
                                                         pop_loadings_std = pop_params$loadings_std,
                                                         x_labels = x_labels, pop_params = pop_params, seed = seed_val,
                                                         sa = 0, sb = 0, int_strength = int_str)
                            res <- tryCatch(run_svt_full(sim_dat, mod_names_local, cfa_model, x_labels, BOOT_R),
                                            error = function(e) list(
                                              step1=list(n_valid_mods=NA, n_mi_sat=NA, n_mi_viol=NA, mean_smi=NA),
                                              step2=list(mean_dl=NA, d_cohen=NA),
                                              step3=list(p_boot=NA, p_binom_delta=NA, p_binom_orient=NA, boot_pass=NA,
                                                         binom_delta_pass=NA, binom_orient_pass=NA, decision="Error"),
                                              summary=list(decision="Error", mean_dl=NA, d_cohen=NA, n_delta_pos=NA,
                                                           n_orient_good=NA, n_mods_valid=NA, mean_smi=NA, n_mi_sat=NA, n_mi_viol=NA)))
                            if (res$summary$decision != "Error") {
                              it <- tryCatch(interaction_test_cfa(sim_dat, x_labels, pi_items, pi_names, sem_interaction_model),
                                             error = function(e) list(flag=FALSE, p=NA, gamma_est=NA, gamma_std=NA))
                              int_flag <- ifelse(is.na(it$flag), FALSE, it$flag)
                              svt_pass <- (res$summary$decision == "Stabilizer")
                              final_decision <- ifelse(svt_pass && !int_flag, "Stabilizer", "Not")
                              res$step3$interaction_flag <- int_flag
                              res$step3$p_interaction <- it$p
                              res$step3$gamma_est <- it$gamma_est
                              res$step3$gamma_std <- it$gamma_std
                              res$step3$decision <- final_decision
                              res$summary$decision <- final_decision
                            }
                            res$rep <- rep; res$config <- cfg; res$scenario <- "H1_Moderator"
                            res$mi_severity <- mi_sev; res$sample_size <- n_sam; res$int_strength <- int_str
                            res
                          }
  all_sim_results <- c(all_sim_results, cond_results)
  decs <- sapply(cond_results, function(x) x$summary$decision)
  n_stab <- sum(decs == "Stabilizer", na.rm = TRUE)
  n_valid <- sum(decs != "Error", na.rm = TRUE)
  n_err <- sum(decs == "Error", na.rm = TRUE)
  rate <- n_stab / max(n_valid, 1)
  dls <- sapply(cond_results, function(x) x$step2$mean_dl)
  dc <- sapply(cond_results, function(x) x$step2$d_cohen)
  bp <- sapply(cond_results, function(x) x$step3$boot_pass)
  bd <- sapply(cond_results, function(x) x$step3$binom_delta_pass)
  bo <- sapply(cond_results, function(x) x$step3$binom_orient_pass)
  ndp <- sapply(cond_results, function(x) x$summary$n_delta_pos)
  nsat <- sapply(cond_results, function(x) x$summary$n_mi_sat)
  nviol <- sapply(cond_results, function(x) x$summary$n_mi_viol)
  smi_vals <- sapply(cond_results, function(x) x$summary$mean_smi)
  cat(sprintf("%-4s | H1_Mod_%.2f | %.1f | %5d | FPR=%.3f | mdl=%.3f\n",
              cfg, int_str, mi_sev, n_sam, rate, mean(dls, na.rm=TRUE)))
  summary_rows[[length(summary_rows)+1]] <- data.frame(
    config=cfg, scenario="H1_Moderator", mi_severity=mi_sev, sample_size=n_sam,
    int_strength=int_str, rate=rate, metric="FPR", n_valid=n_valid, n_error=n_err,
    mean_dl=mean(dls, na.rm=TRUE), sd_dl=sd(dls, na.rm=TRUE),
    mean_d_cohen=mean(dc, na.rm=TRUE),
    boot_rate=mean(bp, na.rm=TRUE), bdelta_rate=mean(bd, na.rm=TRUE), borient_rate=mean(bo, na.rm=TRUE),
    mean_n_delta_pos=mean(ndp, na.rm=TRUE),
    mean_n_orient_good=mean(sapply(cond_results, function(x) x$summary$n_orient_good), na.rm=TRUE),
    mean_n_mods_valid=mean(sapply(cond_results, function(x) x$summary$n_mods_valid), na.rm=TRUE),
    mean_mi_sat=mean(nsat, na.rm=TRUE), mean_mi_viol=mean(nviol, na.rm=TRUE),
    mean_smi=mean(smi_vals, na.rm=TRUE), stringsAsFactors=FALSE)
  if (ci %% 10 == 0) {
    cat(sprintf("  >>> %d/%d graded conditions completed\n", ci, nrow(conditions_b)))
  }
}

stopCluster(cl)
total_hours <- as.numeric(difftime(Sys.time(), t_start, units = "hours"))

cat(paste(rep("=", 70), collapse = ""), "\n")
cat(sprintf("Total computation time: %.1f hours\n\n", total_hours))

cat("PHASE 3: Assembling Results\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

summary_df <- do.call(rbind, summary_rows)

cat("TYPE I ERROR (H0):\n")
h0 <- summary_df[summary_df$scenario == "H0", ]
for (i in seq_len(nrow(h0))) {
  r <- h0[i, ]
  flag <- ifelse(r$rate <= 0.025, "CONSERV", ifelse(r$rate <= 0.075, "BRADLEY", "LIBERAL"))
  cat(sprintf("  %s | N=%4d | MI=%.1f: TypeI=%.3f [%s]\n",
              r$config, r$sample_size, r$mi_severity, r$rate, flag))
}

cat("\nPOWER (H1 stabilization):\n")
h1 <- summary_df[summary_df$scenario %in% c("H1_TypeA","H1_TypeB","H1_TypeAB"), ]
for (scen in c("H1_TypeA","H1_TypeB","H1_TypeAB")) {
  sub <- h1[h1$scenario == scen, ]
  cat(sprintf("  %-10s: Mean=%.3f | Adequate(>=0.80): %d/%d\n",
              scen, mean(sub$rate), sum(sub$rate >= 0.80), nrow(sub)))
}

cat("\nFPR (H1_Moderator - graded interaction):\n")
hmod <- summary_df[summary_df$scenario == "H1_Moderator", ]
for (is_val in sort(unique(hmod$int_strength))) {
  sub <- hmod[hmod$int_strength == is_val, ]
  brad <- sum(sub$rate <= 0.075)
  cat(sprintf("  gamma=%.2f: Mean_FPR=%.3f [%.3f-%.3f] | Bradley: %d/%d\n",
              is_val, mean(sub$rate), min(sub$rate), max(sub$rate), brad, nrow(sub)))
}

saveRDS(list(
  summary = summary_df,
  all_results = all_sim_results,
  pop_params = pop_params,
  design = list(
    n_reps = N_REPS, sample_sizes = SAMPLE_SIZES, mi_levels = MI_LEVELS,
    scenarios = SCENARIOS, interaction_strengths = INT_STRENGTHS,
    mod_configs = MOD_CONFIGS, sa = SA, sb = SB, boot_r = BOOT_R,
    interaction_filter = list(
      method = "Residual-Centered Product Indicator (Little, Bovaird & Widaman, 2006)",
      pi_items = pi_items, constraint = "X_Factor ~~ 0*XZ",
      decision_rule = "SVT_pass AND NOT interaction_flag")),
  total_hours = total_hours,
  date = Sys.time()
), file.path(rds_path, "CFA_full_simulation_results.rds"))

write.csv(summary_df, file.path(output_path, "CFA_full_simulation_summary.csv"), row.names = FALSE)

cat("\nPHASE 4: Generating Figures\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

h0_plot <- summary_df[summary_df$scenario == "H0", ]
hmod_plot <- summary_df[summary_df$scenario == "H1_Moderator", ]
h1_plot <- summary_df[summary_df$scenario %in% c("H1_TypeA","H1_TypeB","H1_TypeAB"), ]

p1 <- ggplot(h0_plot, aes(x = factor(mi_severity), y = rate, fill = factor(sample_size))) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6, alpha = 0.85) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red", linewidth = 0.8) +
  geom_hline(yintercept = c(0.025, 0.075), linetype = "dotted", color = rgb(243,156,18,maxColorValue=255)) +
  facet_wrap(~config, ncol = 1) +
  scale_fill_brewer(palette = "Blues", name = "N") +
  labs(x = "MI Severity", y = "Type I Error Rate",
       title = "Type I Error Control (H0: Z Independent)") +
  theme_minimal() + theme(plot.title = element_text(size = 14, face = "bold"), legend.position = "bottom")

p2 <- ggplot(h1_plot, aes(x = factor(mi_severity), y = rate, fill = factor(sample_size))) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6, alpha = 0.85) +
  geom_hline(yintercept = 0.80, linetype = "dashed", color = rgb(39,174,96,maxColorValue=255)) +
  facet_grid(config ~ scenario) +
  scale_fill_brewer(palette = "Greens", name = "N") +
  labs(x = "MI Severity", y = "Statistical Power",
       title = "Power by Scenario, MI Severity, and Sample Size") +
  theme_minimal() + theme(plot.title = element_text(size = 14, face = "bold"), legend.position = "bottom")

p3 <- ggplot(hmod_plot, aes(x = factor(mi_severity), y = rate, fill = factor(sample_size))) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6, alpha = 0.85) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red", linewidth = 0.8) +
  geom_hline(yintercept = 0.075, linetype = "dotted", color = rgb(243,156,18,maxColorValue=255)) +
  facet_grid(config ~ factor(int_strength, labels = paste0("gamma=", sort(unique(hmod_plot$int_strength))))) +
  scale_fill_brewer(palette = "Oranges", name = "N") +
  labs(x = "MI Severity", y = "False Positive Rate",
       title = "Moderator FPR by Interaction Strength (with CFA Interaction Filter)") +
  theme_minimal() + theme(plot.title = element_text(size = 14, face = "bold"), legend.position = "bottom")

p4 <- ggplot(hmod_plot, aes(x = int_strength, y = rate, color = factor(sample_size))) +
  geom_line(linewidth = 1) + geom_point(size = 2.5) +
  geom_hline(yintercept = c(0.05, 0.075), linetype = c("dashed","dotted"), color = c("red",rgb(243,156,18,maxColorValue=255))) +
  facet_wrap(~config, ncol = 1) +
  scale_color_brewer(palette = "Set1", name = "N") +
  scale_x_continuous(breaks = INT_STRENGTHS) +
  labs(x = "Interaction Strength (gamma)", y = "False Positive Rate",
       title = "SVT Discriminability: FPR as a Function of Interaction Strength") +
  theme_minimal() + theme(plot.title = element_text(size = 14, face = "bold"), legend.position = "bottom")

ggsave(file.path(figure_path, "type1_error.png"), p1, width = 10, height = 8, dpi = 600)
ggsave(file.path(figure_path, "power.png"), p2, width = 14, height = 8, dpi = 600)
ggsave(file.path(figure_path, "moderator_fpr_graded.png"), p3, width = 16, height = 10, dpi = 600)
ggsave(file.path(figure_path, "discriminability_curves.png"), p4, width = 10, height = 8, dpi = 600)

combined <- grid.arrange(p1, p2, p3, p4, ncol = 2, nrow = 2,
                         top = grid::textGrob("CFA-Based SVT Monte Carlo Simulation",
                                              gp = grid::gpar(fontsize = 16, fontface = "bold")))
ggsave(file.path(figure_path, "combined_4panel.png"), combined, width = 22, height = 18, dpi = 600)

cat("Figures saved.\n\n")
cat("SIMULATION COMPLETE\n")
cat(paste(rep("=", 70), collapse = ""), "\n")
cat(sprintf("Total time: %.1f hours\n", total_hours))
cat(sprintf("Total conditions: %d\n", nrow(summary_df)))
cat(sprintf("Total SVT runs: %d\n", length(all_sim_results)))
cat(sprintf("RDS: %s\n", file.path(rds_path, "CFA_full_simulation_results.rds")))
cat(sprintf("CSV: %s\n", file.path(output_path, "CFA_full_simulation_summary.csv")))
cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")