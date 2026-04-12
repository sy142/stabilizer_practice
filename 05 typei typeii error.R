rm(list = ls())
suppressPackageStartupMessages({
  library(readxl)
  library(lavaan)
  library(dplyr)
  library(stringr)
  library(boot)
  library(foreach)
  library(doParallel)
  library(ggplot2)
  library(gridExtra)
})

set.seed(9186)

cat("TYPE I AND TYPE II ERROR ANALYSIS\n")
cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

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

cat("Dataset: analizliksonAVE.xlsx, N =", nrow(data), "\n\n")

N_PERM  <- 1000
N_BOOT  <- 1000
N_CORES <- 14
ALPHA   <- 0.05

cat("Design:\n")
cat("  H0 iterations (permutation): ", N_PERM, "\n")
cat("  H1 iterations (bootstrap):   ", N_BOOT, "\n")
cat("  Cores:                        ", N_CORES, "\n")
cat("  Alpha:                        ", ALPHA, "\n")
cat("  Bradley criterion:            [", ALPHA/2, ",", ALPHA*1.5, "]\n\n")

cfa_model <- '
  MEQ_Factor =~ MEQ1 + MEQ2 + MEQ8 + MEQ9 + MEQ10 + MEQ11 + MEQ15 + MEQ17 + MEQ18 + MEQ19
  MEQ1 ~~ MEQ2
  MEQ2 ~~ MEQ10
  MEQ11 ~~ MEQ15
  MEQ17 ~~ MEQ18
  MEQ2 ~~ MEQ8
'

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

moderators <- c("YasGruplari","Cinsiyet","EgitimDuzeyi","MedeniDurum",
                "GelirDuzeyi","BKISiniflamasi","SigaraKullanma",
                "AlkolKullanma","KronikHastalikVarligi","TibbiBeslenmeTedavisiAlma")

run_svt_core <- function(dat, moderators, cfa_model, model_without, model_with) {
  
  gd_abs <- function(f1, f2, idx) {
    v1 <- as.numeric(f1[idx]); v2 <- as.numeric(f2[idx])
    if (is.na(v1)||is.na(v2)) NA else abs(v1-v2)
  }
  gd_signed <- function(f1, f2, idx) {
    v1 <- as.numeric(f1[idx]); v2 <- as.numeric(f2[idx])
    if (is.na(v1)||is.na(v2)) NA else v1-v2
  }
  
  all_deltas_metric <- list()
  all_deltas_scalar <- list()
  all_invariance <- c()
  mi_details <- list()
  
  for (mod in moderators) {
    if (!(mod %in% names(dat))) next
    d <- dat[!is.na(dat[[mod]]), ]
    d[[mod]] <- factor(d[[mod]])
    if (length(levels(d[[mod]])) < 2) next
    
    fc <- tryCatch(cfa(cfa_model, data=d, group=mod, std.lv=TRUE, estimator="MLR", missing="fiml"), error=function(e) NULL)
    if (is.null(fc)) next
    fm <- tryCatch(cfa(cfa_model, data=d, group=mod, group.equal="loadings", std.lv=TRUE, estimator="MLR", missing="fiml"), error=function(e) NULL)
    if (is.null(fm)) next
    fs <- tryCatch(cfa(cfa_model, data=d, group=mod, group.equal=c("loadings","intercepts"), std.lv=TRUE, estimator="MLR", missing="fiml"), error=function(e) NULL)
    if (is.null(fs)) next
    ft <- tryCatch(cfa(cfa_model, data=d, group=mod, group.equal=c("loadings","intercepts","residuals"), std.lv=TRUE, estimator="MLR", missing="fiml"), error=function(e) NULL)
    if (is.null(ft)) next
    
    fmc <- fitmeasures(fc); fmm <- fitmeasures(fm); fms <- fitmeasures(fs); fmt <- fitmeasures(ft)
    
    dm_abs <- c(CFI=gd_abs(fmc,fmm,"cfi.scaled"), TLI=gd_abs(fmc,fmm,"tli.scaled"),
                RMSEA=gd_abs(fmc,fmm,"rmsea.scaled"), SRMR=gd_abs(fmc,fmm,"srmr"))
    ds_abs <- c(CFI=gd_abs(fmm,fms,"cfi.scaled"), TLI=gd_abs(fmm,fms,"tli.scaled"),
                RMSEA=gd_abs(fmm,fms,"rmsea.scaled"), SRMR=gd_abs(fmm,fms,"srmr"))
    
    dc_m <- gd_signed(fmc,fmm,"cfi.scaled"); dt_m <- gd_signed(fmc,fmm,"tli.scaled")
    dr_m <- -gd_signed(fmc,fmm,"rmsea.scaled"); ds_m <- -gd_signed(fmc,fmm,"srmr")
    dc_s <- gd_signed(fmm,fms,"cfi.scaled"); dt_s <- gd_signed(fmm,fms,"tli.scaled")
    dr_s <- -gd_signed(fmm,fms,"rmsea.scaled"); ds_s <- -gd_signed(fmm,fms,"srmr")
    dc_t <- gd_signed(fms,fmt,"cfi.scaled"); dt_t <- gd_signed(fms,fmt,"tli.scaled")
    dr_t <- -gd_signed(fms,fmt,"rmsea.scaled"); ds_t <- -gd_signed(fms,fmt,"srmr")
    
    if (any(is.na(c(dm_abs, ds_abs, dc_t, dt_t, dr_t, ds_t)))) next
    
    all_deltas_metric[[mod]] <- dm_abs
    all_deltas_scalar[[mod]] <- ds_abs
    
    th_c <- 0.010; th_t <- 0.010; th_r <- 0.015; th_s <- 0.030
    is_inv <- (dc_m<=th_c && dt_m<=th_t && dr_m<=th_r && ds_m<=th_s &&
                 dc_s<=th_c && dt_s<=th_t && dr_s<=th_r && ds_s<=th_s &&
                 dc_t<=th_c && dt_t<=th_t && dr_t<=th_r && ds_t<=th_s)
    all_invariance <- c(all_invariance, is_inv)
    names(all_invariance)[length(all_invariance)] <- mod
    
    mi_details[[mod]] <- list(dc_m=dc_m,dt_m=dt_m,dr_m=dr_m,ds_m=ds_m,
                              dc_s=dc_s,dt_s=dt_s,dr_s=dr_s,ds_s=ds_s,
                              dc_t=dc_t,dt_t=dt_t,dr_t=dr_t,ds_t=ds_t,
                              is_invariant=is_inv)
  }
  
  if (length(all_deltas_metric) >= 3) {
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
    weights <- as.numeric(raw_w/sum(raw_w))
  } else {
    weights <- c(0.40, 0.10, 0.30, 0.20)
  }
  
  mi_results <- list()
  for (mod in names(mi_details)) {
    md <- mi_details[[mod]]
    nc <- min(max(pmax(0,md$dc_m),pmax(0,md$dc_s),pmax(0,md$dc_t))/0.010, 1)
    nt <- min(max(pmax(0,md$dt_m),pmax(0,md$dt_s),pmax(0,md$dt_t))/0.010, 1)
    nr <- min(max(pmax(0,md$dr_m),pmax(0,md$dr_s),pmax(0,md$dr_t))/0.015, 1)
    ns <- min(max(pmax(0,md$ds_m),pmax(0,md$ds_s),pmax(0,md$ds_t))/0.030, 1)
    mi_score <- weights[1]*nc + weights[2]*nt + weights[3]*nr + weights[4]*ns
    mi_results[[mod]] <- list(mi_score=mi_score, is_invariant=md$is_invariant)
  }
  
  stab_results <- list()
  for (mod in moderators) {
    if (!(mod %in% names(dat))) next
    d <- dat[!is.na(dat[[mod]]), ]
    d[[mod]] <- factor(d[[mod]])
    groups <- levels(d[[mod]])
    if (length(groups) < 2) next
    
    c_paths <- c(); cprime_paths <- c()
    group_ns <- c()
    
    for (g in groups) {
      gd <- d[d[[mod]] == g, ]
      if (nrow(gd) < 30) next
      f0 <- tryCatch(sem(model_without, data=gd, std.lv=TRUE, estimator="MLR", missing="fiml"),
                     error=function(e) NULL)
      f1 <- tryCatch(sem(model_with, data=gd, std.lv=TRUE, estimator="MLR", missing="fiml"),
                     error=function(e) NULL)
      if (!is.null(f0) && !is.null(f1)) {
        s0 <- standardizedSolution(f0); s1 <- standardizedSolution(f1)
        cp <- s0[s0$label=="c","est.std"][1]
        cpp <- s1[s1$label=="cprime","est.std"][1]
        if (!is.na(cp) && !is.na(cpp)) {
          c_paths <- c(c_paths, cp); cprime_paths <- c(cprime_paths, cpp)
          group_ns <- c(group_ns, nrow(gd))
        }
      }
    }
    
    if (length(c_paths) < 2) next
    
    wmn <- function(x, w) weighted.mean(x, w)
    wsd <- function(x, w) sqrt(sum(w * (x - wmn(x,w))^2) / sum(w))
    
    mc <- wmn(c_paths, group_ns); mcp <- wmn(cprime_paths, group_ns)
    sc <- wsd(c_paths, group_ns); scp <- wsd(cprime_paths, group_ns)
    em <- 0.0001
    cv0 <- (sc / max(abs(mc), em)) * 100
    cv1 <- (scp / max(abs(mcp), em)) * 100
    el <- 1e-8
    dl <- log(cv0 + el) - log(cv1 + el)
    
    ig <- sapply(seq_along(c_paths), function(i) {
      ifelse(abs(cprime_paths[i] - mcp) < abs(c_paths[i] - mc), 1, 0)
    })
    ocr <- mean(ig)
    
    stab_results[[mod]] <- list(delta_log = dl, ocr = ocr, total_n = sum(group_ns))
  }
  
  dvec <- sapply(stab_results, function(x) x$delta_log)
  dvec <- dvec[!is.na(dvec)]
  if (length(dvec) < 2) return(list(p_boot=NA, p_binom=NA, mean_dl=NA, decision="Error",
                                    n_mi_violated=NA, n_mi_satisfied=NA))
  
  wts <- sapply(names(dvec), function(m) stab_results[[m]]$total_n)
  ocrs <- sapply(names(dvec), function(m) stab_results[[m]]$ocr)
  
  bfn <- function(d, idx) { w <- wts[idx]/sum(wts[idx]); sum(d[idx]*w) }
  br <- tryCatch(boot(dvec, bfn, R = 1000), error = function(e) NULL)
  if (!is.null(br)) {
    mdl <- br$t0; sedl <- sd(br$t)
    zs <- mdl / sedl
    pb <- pnorm(zs, lower.tail = FALSE)
  } else {
    mdl <- weighted.mean(dvec, wts); pb <- NA
  }
  
  c2p <- (dvec > 0) | (ocrs >= 0.5)
  pbn <- binom.test(sum(c2p), length(c2p), p = 0.5, alternative = "greater")$p.value
  
  bpass <- (!is.na(pb) && pb < 0.05 && mdl > 0)
  bnpass <- (pbn < 0.05)
  
  all_inv <- sapply(mi_results, function(x) x$is_invariant)
  
  list(p_boot = pb, p_binom = pbn, mean_dl = mdl,
       decision = ifelse(bpass && bnpass, "Stabilizer", "Not"),
       n_mi_violated = sum(!all_inv), n_mi_satisfied = sum(all_inv))
}

cat("Setting up parallel cluster...\n")
cl <- makeCluster(N_CORES)
registerDoParallel(cl)

clusterExport(cl, c("data", "moderators", "model_without", "model_with",
                    "cfa_model", "run_svt_core"))
clusterEvalQ(cl, {
  suppressPackageStartupMessages({ library(lavaan); library(dplyr); library(boot) })
})

cat("\nPHASE 1: Type I Error (H0 - Permuted BRIAN)\n")
cat(paste(rep("=", 60), collapse = ""), "\n")
cat("Shuffling BRIAN to break true relationship...\n\n")

t1 <- Sys.time()

h0_results <- foreach(iter = 1:N_PERM, .packages = c("lavaan","dplyr","boot"),
                      .errorhandling = "pass") %dopar% {
                        dat_perm <- data
                        dat_perm$BRIAN <- sample(dat_perm$BRIAN)
                        res <- tryCatch(
                          run_svt_core(dat_perm, moderators, cfa_model, model_without, model_with),
                          error = function(e) list(p_boot = NA, p_binom = NA, mean_dl = NA, decision = "Error")
                        )
                        res$iter <- iter
                        res
                      }

elapsed1 <- difftime(Sys.time(), t1, units = "mins")
cat(sprintf("  Completed %d H0 iterations in %.1f minutes\n\n", N_PERM, as.numeric(elapsed1)))

cat("PHASE 2: Power (H1 - Bootstrap Resampled Real Data)\n")
cat(paste(rep("=", 60), collapse = ""), "\n")
cat("Resampling with replacement from real data...\n\n")

t2 <- Sys.time()

h1_results <- foreach(iter = 1:N_BOOT, .packages = c("lavaan","dplyr","boot"),
                      .errorhandling = "pass") %dopar% {
                        idx <- sample(nrow(data), replace = TRUE)
                        dat_boot <- data[idx, ]
                        dat_boot$GozlemNo <- 1:nrow(dat_boot)
                        res <- tryCatch(
                          run_svt_core(dat_boot, moderators, cfa_model, model_without, model_with),
                          error = function(e) list(p_boot = NA, p_binom = NA, mean_dl = NA, decision = "Error")
                        )
                        res$iter <- iter
                        res
                      }

elapsed2 <- difftime(Sys.time(), t2, units = "mins")
cat(sprintf("  Completed %d H1 iterations in %.1f minutes\n\n", N_BOOT, as.numeric(elapsed2)))

stopCluster(cl)

h0_pboot  <- sapply(h0_results, function(x) if(is.null(x$p_boot)) NA else x$p_boot)
h0_pbinom <- sapply(h0_results, function(x) if(is.null(x$p_binom)) NA else x$p_binom)
h0_dl     <- sapply(h0_results, function(x) if(is.null(x$mean_dl)) NA else x$mean_dl)
h0_dec    <- sapply(h0_results, function(x) if(is.null(x$decision)) "Error" else x$decision)

h1_pboot  <- sapply(h1_results, function(x) if(is.null(x$p_boot)) NA else x$p_boot)
h1_pbinom <- sapply(h1_results, function(x) if(is.null(x$p_binom)) NA else x$p_binom)
h1_dl     <- sapply(h1_results, function(x) if(is.null(x$mean_dl)) NA else x$mean_dl)
h1_dec    <- sapply(h1_results, function(x) if(is.null(x$decision)) "Error" else x$decision)

type1_boot  <- mean(h0_pboot < ALPHA, na.rm = TRUE)
type1_binom <- mean(h0_pbinom < ALPHA, na.rm = TRUE)
type1_dual  <- mean(h0_dec == "Stabilizer", na.rm = TRUE)

power_boot  <- mean(h1_pboot < ALPHA & h1_dl > 0, na.rm = TRUE)
power_binom <- mean(h1_pbinom < ALPHA, na.rm = TRUE)
power_dual  <- mean(h1_dec == "Stabilizer", na.rm = TRUE)

n_valid_h0 <- sum(!is.na(h0_pboot))
n_valid_h1 <- sum(!is.na(h1_pboot))

bradley_lower <- ALPHA / 2
bradley_upper <- ALPHA * 1.5
bradley_pass  <- (type1_dual >= bradley_lower) && (type1_dual <= bradley_upper)

cat("RESULTS\n")
cat(paste(rep("=", 60), collapse = ""), "\n\n")

cat("Valid iterations: H0 =", n_valid_h0, "/", N_PERM, "| H1 =", n_valid_h1, "/", N_BOOT, "\n\n")

cat("Type I Error Rates (H0: permuted BRIAN):\n")
cat("  Bootstrap criterion only:   ", sprintf("%.4f", type1_boot), "\n")
cat("  Binomial criterion only:    ", sprintf("%.4f", type1_binom), "\n")
cat("  Dual-criterion (Eq. 48):    ", sprintf("%.4f", type1_dual), "\n\n")

cat("Statistical Power (H1: real data bootstrap):\n")
cat("  Bootstrap criterion only:   ", sprintf("%.4f", power_boot), "\n")
cat("  Binomial criterion only:    ", sprintf("%.4f", power_binom), "\n")
cat("  Dual-criterion (Eq. 48):    ", sprintf("%.4f", power_dual), "\n\n")

cat("Bradley's Liberal Criterion:\n")
cat("  Bounds: [", sprintf("%.3f", bradley_lower), ",", sprintf("%.3f", bradley_upper), "]\n")
cat("  Observed Type I (dual):     ", sprintf("%.4f", type1_dual), "\n")
cat("  Within bounds:              ", ifelse(bradley_pass, "YES", "NO"), "\n\n")

if (bradley_pass) {
  cat("CONCLUSION: SVT dual-criterion maintains nominal Type I error rate.\n")
  cat("  Bradley's liberal criterion SATISFIED.\n")
} else if (type1_dual < bradley_lower) {
  cat("CONCLUSION: SVT is CONSERVATIVE (Type I below lower bound).\n")
  cat("  This is acceptable - reduces false positives at cost of some power.\n")
} else {
  cat("CONCLUSION: SVT Type I error EXCEEDS Bradley upper bound.\n")
  cat("  Investigate further.\n")
}

cat(sprintf("\nPower = %.1f%% with Type I = %.2f%%\n", power_dual*100, type1_dual*100))

mc_results <- list(
  type1_boot = type1_boot, type1_binom = type1_binom, type1_dual = type1_dual,
  power_boot = power_boot, power_binom = power_binom, power_dual = power_dual,
  bradley_lower = bradley_lower, bradley_upper = bradley_upper, bradley_pass = bradley_pass,
  h0_pboot = h0_pboot, h0_pbinom = h0_pbinom, h0_dl = h0_dl, h0_dec = h0_dec,
  h1_pboot = h1_pboot, h1_pbinom = h1_pbinom, h1_dl = h1_dl, h1_dec = h1_dec,
  n_perm = N_PERM, n_boot = N_BOOT, alpha = ALPHA,
  n_valid_h0 = n_valid_h0, n_valid_h1 = n_valid_h1,
  total_time_min = as.numeric(elapsed1) + as.numeric(elapsed2),
  date = Sys.time()
)

saveRDS(mc_results, file.path(rds_path, "TypeI_TypeII_error.rds"))

summary_df <- data.frame(
  Metric = c("Type I (bootstrap only)", "Type I (binomial only)", "Type I (dual-criterion)",
             "Power (bootstrap only)", "Power (binomial only)", "Power (dual-criterion)",
             "Bradley lower", "Bradley upper", "Bradley pass",
             "N permutations", "N bootstrap", "Valid H0", "Valid H1"),
  Value = c(sprintf("%.4f", type1_boot), sprintf("%.4f", type1_binom), sprintf("%.4f", type1_dual),
            sprintf("%.4f", power_boot), sprintf("%.4f", power_binom), sprintf("%.4f", power_dual),
            sprintf("%.3f", bradley_lower), sprintf("%.3f", bradley_upper), as.character(bradley_pass),
            N_PERM, N_BOOT, n_valid_h0, n_valid_h1),
  stringsAsFactors = FALSE
)
write.csv(summary_df, file.path(output_path, "TypeI_TypeII_summary.csv"), row.names = FALSE)

cat("\nSaved:\n")
cat("  RDS/TypeI_TypeII_error.rds\n")
cat("  Outputs/TypeI_TypeII_summary.csv\n")

cat("\nGenerating figures...\n")

df_pvals <- data.frame(
  p_value = c(h0_pboot[!is.na(h0_pboot)], h1_pboot[!is.na(h1_pboot)]),
  Hypothesis = factor(rep(c("Null (H0: permuted)", "Alternative (H1: real)"),
                          c(sum(!is.na(h0_pboot)), sum(!is.na(h1_pboot)))),
                      levels = c("Null (H0: permuted)", "Alternative (H1: real)"))
)

p1 <- ggplot(df_pvals, aes(x = p_value, fill = Hypothesis)) +
  geom_histogram(aes(y = after_stat(density)), binwidth = 0.025, alpha = 0.6,
                 position = "identity", color = "white", linewidth = 0.2) +
  geom_density(aes(color = Hypothesis), alpha = 0, linewidth = 1.2, adjust = 1.5) +
  geom_vline(xintercept = ALPHA, linetype = "dashed", color = "#E74C3C", linewidth = 1) +
  geom_vline(xintercept = bradley_lower, linetype = "dotted", color = "#F39C12", linewidth = 0.8) +
  geom_vline(xintercept = bradley_upper, linetype = "dotted", color = "#F39C12", linewidth = 0.8) +
  scale_fill_manual(values = c("Null (H0: permuted)" = "#3498DB", "Alternative (H1: real)" = "#2ECC71")) +
  scale_color_manual(values = c("Null (H0: permuted)" = "#2874A6", "Alternative (H1: real)" = "#1E8449")) +
  scale_x_continuous(breaks = seq(0, 1, 0.1), limits = c(0, 1)) +
  annotate("text", x = ALPHA, y = Inf, label = paste0("alpha = ", ALPHA),
           hjust = -0.2, vjust = 2, color = "#E74C3C", size = 3.5, fontface = "italic") +
  annotate("text", x = bradley_lower, y = Inf, label = "Bradley\nlower",
           hjust = 1.2, vjust = 3.5, color = "#F39C12", size = 2.8, fontface = "italic") +
  annotate("text", x = bradley_upper, y = Inf, label = "Bradley\nupper",
           hjust = -0.2, vjust = 3.5, color = "#F39C12", size = 2.8, fontface = "italic") +
  labs(title = "Bootstrap p-value Distributions Under H0 and H1",
       subtitle = sprintf("Type I Error = %.3f | Power = %.3f | Bradley: %s",
                          type1_dual, power_dual, ifelse(bradley_pass, "PASS", "FAIL")),
       x = "p-value (bootstrap)", y = "Density") +
  theme_minimal() +
  theme(plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 11),
        legend.position = "bottom")

df_dl <- data.frame(
  mean_dl = c(h0_dl[!is.na(h0_dl)], h1_dl[!is.na(h1_dl)]),
  Hypothesis = factor(rep(c("H0 (permuted)", "H1 (real)"),
                          c(sum(!is.na(h0_dl)), sum(!is.na(h1_dl)))),
                      levels = c("H0 (permuted)", "H1 (real)"))
)

p2 <- ggplot(df_dl, aes(x = mean_dl, fill = Hypothesis)) +
  geom_histogram(aes(y = after_stat(density)), binwidth = 0.1, alpha = 0.6,
                 position = "identity", color = "white", linewidth = 0.2) +
  geom_density(aes(color = Hypothesis), alpha = 0, linewidth = 1.2, adjust = 1.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
  scale_fill_manual(values = c("H0 (permuted)" = "#3498DB", "H1 (real)" = "#2ECC71")) +
  scale_color_manual(values = c("H0 (permuted)" = "#2874A6", "H1 (real)" = "#1E8449")) +
  labs(title = expression(paste("Distribution of Weighted Mean ", Delta, ell, " Under H0 and H1")),
       subtitle = sprintf("H0 mean = %.3f | H1 mean = %.3f | Separation confirms test validity",
                          mean(h0_dl, na.rm=TRUE), mean(h1_dl, na.rm=TRUE)),
       x = expression(paste("Weighted Mean ", Delta, ell)), y = "Density") +
  theme_minimal() +
  theme(plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 11),
        legend.position = "bottom")

decision_data <- data.frame(
  Criterion = factor(rep(c("Bootstrap\nonly", "Binomial\nonly", "Dual-criterion\n(Eq. 48)"), 2),
                     levels = c("Bootstrap\nonly", "Binomial\nonly", "Dual-criterion\n(Eq. 48)")),
  Hypothesis = factor(rep(c("Type I Error", "Power"), each = 3),
                      levels = c("Type I Error", "Power")),
  Rate = c(type1_boot, type1_binom, type1_dual, power_boot, power_binom, power_dual)
)

p3 <- ggplot(decision_data, aes(x = Criterion, y = Rate, fill = Hypothesis)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6, alpha = 0.85) +
  geom_hline(yintercept = ALPHA, linetype = "dashed", color = "red", linewidth = 0.8) +
  geom_hline(yintercept = bradley_lower, linetype = "dotted", color = "#F39C12") +
  geom_hline(yintercept = bradley_upper, linetype = "dotted", color = "#F39C12") +
  geom_text(aes(label = sprintf("%.3f", Rate)), position = position_dodge(width = 0.7),
            vjust = -0.5, size = 3) +
  scale_fill_manual(values = c("Type I Error" = "#e74c3c", "Power" = "#27ae60")) +
  labs(title = "Type I Error and Power by Decision Criterion",
       subtitle = "Orange dotted = Bradley bounds. Red dashed = nominal alpha.",
       x = "", y = "Rate") +
  theme_minimal() +
  theme(plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 11),
        legend.position = "bottom")

alphas <- seq(0.01, 0.20, by = 0.005)
roc_data <- data.frame(
  alpha = alphas,
  type1 = sapply(alphas, function(a) mean(h0_dec == "Stabilizer" | h0_pboot < a, na.rm = TRUE)),
  power = sapply(alphas, function(a) mean(h1_pboot < a & h1_dl > 0 & h1_pbinom < a, na.rm = TRUE))
)

p4 <- ggplot(roc_data, aes(x = alpha)) +
  geom_line(aes(y = type1, color = "Type I Error"), linewidth = 1.2) +
  geom_line(aes(y = power, color = "Power"), linewidth = 1.2) +
  geom_vline(xintercept = ALPHA, linetype = "dashed", color = "gray40") +
  geom_hline(yintercept = 0.80, linetype = "dotted", color = "gray60") +
  scale_color_manual(values = c("Type I Error" = "#e74c3c", "Power" = "#27ae60")) +
  labs(title = "Type I Error and Power Across Alpha Levels",
       subtitle = "Vertical line = alpha = 0.05. Horizontal = 80% power threshold.",
       x = "Significance Level (alpha)", y = "Rate", color = "") +
  theme_minimal() +
  theme(plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 11),
        legend.position = "bottom")

ggsave(file.path(figure_path, "error_pvalue_distributions.png"), p1, width = 10, height = 6, dpi = 600)
ggsave(file.path(figure_path, "error_dl_distributions.png"), p2, width = 10, height = 6, dpi = 600)
ggsave(file.path(figure_path, "error_criterion_comparison.png"), p3, width = 10, height = 6, dpi = 600)
ggsave(file.path(figure_path, "error_power_curve.png"), p4, width = 10, height = 6, dpi = 600)

combined <- grid.arrange(p1, p2, p3, p4, ncol = 2, nrow = 2,
                         top = grid::textGrob("SVT Error Rate Analysis: Empirical Type I Error and Power",
                                              gp = grid::gpar(fontsize = 16, fontface = "bold")))
ggsave(file.path(figure_path, "error_combined_4panel.png"),
       combined, width = 16, height = 12, dpi = 600)

cat("Figures saved to Figures/\n")
cat(sprintf("\nTotal time: %.1f minutes\n", as.numeric(elapsed1) + as.numeric(elapsed2)))
cat("\nDone.\n")
cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")