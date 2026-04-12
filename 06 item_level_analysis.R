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

cat("ITEM-LEVEL STABILIZATION ANALYSIS\n")
cat("Dataset: N =", nrow(data), "\n\n")

all_items <- c("MEQ1","MEQ2","MEQ8","MEQ9","MEQ10","MEQ11","MEQ15","MEQ17","MEQ18","MEQ19")

cov_pairs <- list(c("MEQ1","MEQ2"), c("MEQ2","MEQ10"), c("MEQ11","MEQ15"),
                  c("MEQ17","MEQ18"), c("MEQ2","MEQ8"))

moderators <- c("YasGruplari","Cinsiyet","EgitimDuzeyi","MedeniDurum",
                "GelirDuzeyi","BKISiniflamasi","SigaraKullanma",
                "AlkolKullanma","KronikHastalikVarligi","TibbiBeslenmeTedavisiAlma")

run_svt_with_items <- function(data, items, covs, moderators) {
  cfa_spec <- paste0("MEQ_Factor =~ ", paste(items, collapse = " + "))
  cov_lines <- c()
  for (cp in covs) {
    if (cp[1] %in% items && cp[2] %in% items) {
      cov_lines <- c(cov_lines, paste0(cp[1], " ~~ ", cp[2]))
    }
  }
  cfa_model <- paste(c(cfa_spec, cov_lines), collapse = "\n")

  m_without <- paste(c(cfa_spec, "YJIPAQ ~ c*MEQ_Factor", cov_lines), collapse = "\n")
  m_with <- paste(c(cfa_spec, "BRIAN ~ a*MEQ_Factor",
                     "YJIPAQ ~ cprime*MEQ_Factor + b*BRIAN", cov_lines), collapse = "\n")

  gd_abs <- function(f1, f2, idx) {
    v1 <- as.numeric(f1[idx]); v2 <- as.numeric(f2[idx])
    if (is.na(v1)||is.na(v2)) NA else abs(v1-v2)
  }
  gd_signed <- function(f1, f2, idx) {
    v1 <- as.numeric(f1[idx]); v2 <- as.numeric(f2[idx])
    if (is.na(v1)||is.na(v2)) NA else v1-v2
  }

  all_dm <- list(); all_ds <- list(); all_inv <- c()
  mi_details <- list()

  for (mod in moderators) {
    if (!(mod %in% names(data))) next
    d <- data[!is.na(data[[mod]]), ]
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
    dm <- c(CFI=gd_abs(fmc,fmm,"cfi.scaled"), TLI=gd_abs(fmc,fmm,"tli.scaled"),
            RMSEA=gd_abs(fmc,fmm,"rmsea.scaled"), SRMR=gd_abs(fmc,fmm,"srmr"))
    ds <- c(CFI=gd_abs(fmm,fms,"cfi.scaled"), TLI=gd_abs(fmm,fms,"tli.scaled"),
            RMSEA=gd_abs(fmm,fms,"rmsea.scaled"), SRMR=gd_abs(fmm,fms,"srmr"))
    dc_m <- gd_signed(fmc,fmm,"cfi.scaled"); dt_m <- gd_signed(fmc,fmm,"tli.scaled")
    dr_m <- -gd_signed(fmc,fmm,"rmsea.scaled"); ds_m <- -gd_signed(fmc,fmm,"srmr")
    dc_s <- gd_signed(fmm,fms,"cfi.scaled"); dt_s <- gd_signed(fmm,fms,"tli.scaled")
    dr_s <- -gd_signed(fmm,fms,"rmsea.scaled"); ds_s <- -gd_signed(fmm,fms,"srmr")
    dc_t <- gd_signed(fms,fmt,"cfi.scaled"); dt_t <- gd_signed(fms,fmt,"tli.scaled")
    dr_t <- -gd_signed(fms,fmt,"rmsea.scaled"); ds_t <- -gd_signed(fms,fmt,"srmr")

    if (any(is.na(c(dm, ds, dc_t, dt_t, dr_t, ds_t)))) next
    all_dm[[mod]] <- dm; all_ds[[mod]] <- ds
    th_c <- 0.010; th_t <- 0.010; th_r <- 0.015; th_s <- 0.030
    is_inv <- (dc_m<=th_c && dt_m<=th_t && dr_m<=th_r && ds_m<=th_s &&
               dc_s<=th_c && dt_s<=th_t && dr_s<=th_r && ds_s<=th_s &&
               dc_t<=th_c && dt_t<=th_t && dr_t<=th_r && ds_t<=th_s)
    all_inv <- c(all_inv, is_inv); names(all_inv)[length(all_inv)] <- mod
    mi_details[[mod]] <- list(dc_m=dc_m,dt_m=dt_m,dr_m=dr_m,ds_m=ds_m,
                              dc_s=dc_s,dt_s=dt_s,dr_s=dr_s,ds_s=ds_s,
                              dc_t=dc_t,dt_t=dt_t,dr_t=dr_t,ds_t=ds_t, is_invariant=is_inv)
  }

  if (length(all_dm) >= 3) {
    dm_mat <- do.call(rbind, all_dm); ds_mat <- do.call(rbind, all_ds)
    dmax <- pmax(dm_mat, ds_mat)
    cormat <- cor(dmax)
    nm <- cbind(CFI=pmin(dmax[,"CFI"]/0.010,1), TLI=pmin(dmax[,"TLI"]/0.010,1),
                RMSEA=pmin(dmax[,"RMSEA"]/0.015,1), SRMR=pmin(dmax[,"SRMR"]/0.030,1))
    cv_v <- apply(nm, 2, function(x) sd(x)/(mean(x)+0.001))
    mcor <- apply(abs(cormat), 1, function(x) mean(x[x<1]))
    rp <- 1/(1+mcor)
    dp <- numeric(4); names(dp) <- c("CFI","TLI","RMSEA","SRMR")
    for (idx in names(dp)) {
      if (sum(all_inv)>0 && sum(!all_inv)>0) dp[idx] <- max(0, mean(nm[!all_inv,idx])-mean(nm[all_inv,idx]))
      else dp[idx] <- sd(nm[,idx])
    }
    vs <- 1/(1+cv_v); raw_w <- dp*rp*vs; weights <- as.numeric(raw_w/sum(raw_w))
  } else weights <- c(0.40, 0.10, 0.30, 0.20)

  stab <- list()
  for (mod in moderators) {
    if (!(mod %in% names(data))) next
    d <- data[!is.na(data[[mod]]), ]
    d[[mod]] <- factor(d[[mod]])
    groups <- levels(d[[mod]])
    if (length(groups) < 2) next
    c_p <- c(); cp_p <- c(); g_ns <- c()
    for (g in groups) {
      gd <- d[d[[mod]]==g, ]
      if (nrow(gd) < 30) next
      f0 <- tryCatch(sem(m_without, data=gd, std.lv=TRUE, estimator="MLR", missing="fiml"), error=function(e) NULL)
      f1 <- tryCatch(sem(m_with, data=gd, std.lv=TRUE, estimator="MLR", missing="fiml"), error=function(e) NULL)
      if (!is.null(f0) && !is.null(f1)) {
        s0 <- standardizedSolution(f0); s1 <- standardizedSolution(f1)
        cp <- s0[s0$label=="c","est.std"][1]; cpp <- s1[s1$label=="cprime","est.std"][1]
        if (!is.na(cp) && !is.na(cpp)) { c_p <- c(c_p, cp); cp_p <- c(cp_p, cpp); g_ns <- c(g_ns, nrow(gd)) }
      }
    }
    if (length(c_p) < 2) next
    wmn <- function(x,w) weighted.mean(x,w)
    wsd <- function(x,w) sqrt(sum(w*(x-wmn(x,w))^2)/sum(w))
    mc <- wmn(c_p,g_ns); mcp <- wmn(cp_p,g_ns); sc <- wsd(c_p,g_ns); scp <- wsd(cp_p,g_ns)
    em <- 0.0001; el <- 1e-8
    cv0 <- (sc/max(abs(mc),em))*100; cv1 <- (scp/max(abs(mcp),em))*100
    dl <- log(cv0+el)-log(cv1+el)
    ig <- sapply(seq_along(c_p), function(i) ifelse(abs(cp_p[i]-mcp)<abs(c_p[i]-mc),1,0))
    stab[[mod]] <- list(delta_log=dl, ocr=mean(ig), total_n=sum(g_ns))
  }

  dvec <- sapply(stab, function(x) x$delta_log); dvec <- dvec[!is.na(dvec)]
  if (length(dvec) < 2) return(list(mean_dl=NA, p_boot=NA, p_binom=NA, d_cohen=NA, decision="Error"))

  wts <- sapply(names(dvec), function(m) stab[[m]]$total_n)
  ocrs <- sapply(names(dvec), function(m) stab[[m]]$ocr)
  bfn <- function(d, idx) { w <- wts[idx]/sum(wts[idx]); sum(d[idx]*w) }
  br <- tryCatch(boot(dvec, bfn, R=1000), error=function(e) NULL)
  if (!is.null(br)) { mdl <- br$t0; pb <- pnorm(mdl/sd(br$t), lower.tail=FALSE) }
  else { mdl <- weighted.mean(dvec, wts); pb <- NA }

  wvar <- sum(wts*(dvec-mdl)^2)/sum(wts)
  dc <- if(wvar>0) mdl/sqrt(wvar) else NA
  c2p <- (dvec>0)|(ocrs>=0.5)
  pbn <- binom.test(sum(c2p), length(c2p), p=0.5, alternative="greater")$p.value
  bpass <- (!is.na(pb) && pb<0.05 && mdl>0); bnpass <- (pbn<0.05)

  list(mean_dl=mdl, p_boot=pb, p_binom=pbn, d_cohen=dc,
       decision=ifelse(bpass&&bnpass, "Stabilizer", "Not"),
       n_mi_violated=sum(!all_inv), n_mi_satisfied=sum(all_inv))
}

cat("PART 1: Leave-One-Out Item Analysis\n")
cat(paste(rep("=", 50), collapse = ""), "\n\n")

cat("Full model (10 items) baseline...\n")
full_res <- run_svt_with_items(data, all_items, cov_pairs, moderators)
cat(sprintf("  Full: Dl=%.3f, p_boot=%.4f, p_binom=%.4f, d=%.3f, MI: %d/%d, %s\n\n",
            full_res$mean_dl, full_res$p_boot, full_res$p_binom, full_res$d_cohen,
            full_res$n_mi_satisfied, full_res$n_mi_satisfied+full_res$n_mi_violated,
            full_res$decision))

loo_results <- list()
for (item in all_items) {
  reduced_items <- setdiff(all_items, item)
  cat(sprintf("  Removing %s (9 items)... ", item))
  res <- run_svt_with_items(data, reduced_items, cov_pairs, moderators)
  loo_results[[item]] <- res
  delta_from_full <- full_res$mean_dl - res$mean_dl
  cat(sprintf("Dl=%.3f (change=%+.3f), d=%.3f, %s\n",
              res$mean_dl, delta_from_full, res$d_cohen, res$decision))
}

loo_df <- data.frame(
  item_removed = names(loo_results),
  dl = sapply(loo_results, function(x) x$mean_dl),
  p_boot = sapply(loo_results, function(x) x$p_boot),
  p_binom = sapply(loo_results, function(x) x$p_binom),
  d_cohen = sapply(loo_results, function(x) x$d_cohen),
  decision = sapply(loo_results, function(x) x$decision),
  dl_change = full_res$mean_dl - sapply(loo_results, function(x) x$mean_dl),
  stringsAsFactors = FALSE
)

cat("\nPART 2: Per-Item Loading Variance Across Groups\n")
cat(paste(rep("=", 50), collapse = ""), "\n\n")

cfa_model_full <- paste(c(
  paste0("MEQ_Factor =~ ", paste(all_items, collapse = " + ")),
  sapply(cov_pairs, function(x) paste0(x[1], " ~~ ", x[2]))
), collapse = "\n")

loading_var_results <- list()

for (mod in moderators) {
  if (!(mod %in% names(data))) next
  d <- data[!is.na(data[[mod]]), ]
  d[[mod]] <- factor(d[[mod]])
  groups <- levels(d[[mod]])
  if (length(groups) < 2) next

  group_loadings <- list()
  for (g in groups) {
    gd <- d[d[[mod]] == g, ]
    if (nrow(gd) < 50) next
    fit <- tryCatch(cfa(cfa_model_full, data=gd, std.lv=TRUE, estimator="MLR"),
                    error=function(e) NULL)
    if (!is.null(fit)) {
      std <- standardizedSolution(fit)
      lds <- std[std$op == "=~", c("rhs", "est.std")]
      group_loadings[[g]] <- setNames(lds$est.std, lds$rhs)
    }
  }

  if (length(group_loadings) >= 2) {
    ld_mat <- do.call(rbind, group_loadings)
    item_vars <- apply(ld_mat, 2, var)
    loading_var_results[[mod]] <- item_vars
  }
}

if (length(loading_var_results) > 0) {
  lv_mat <- do.call(rbind, loading_var_results)
  mean_var_per_item <- colMeans(lv_mat)
  cat("Mean loading variance across moderators (higher = more MI violation):\n")
  for (item in names(sort(mean_var_per_item, decreasing = TRUE))) {
    cat(sprintf("  %-6s: %.6f\n", item, mean_var_per_item[item]))
  }
}

cat("\nPART 3: Correlation Between Item-Level MI and Stabilization\n")
cat(paste(rep("=", 50), collapse = ""), "\n\n")

if (length(loading_var_results) > 0) {
  cat("Item | Loading Var | Dl Change (when removed) | Interpretation\n")
  cat(paste(rep("-", 70), collapse = ""), "\n")
  for (item in all_items) {
    lv <- mean_var_per_item[item]
    dc <- loo_df[loo_df$item_removed == item, "dl_change"]
    interp <- if (dc > 0.1) "Stabilization DECREASES when removed (contributes to stabilization)"
              else if (dc < -0.1) "Stabilization INCREASES when removed (hinders stabilization)"
              else "Minimal impact"
    cat(sprintf("  %-6s | %.6f    | %+.3f                  | %s\n", item, lv, dc, interp))
  }
}

write.csv(loo_df, file.path(output_path, "item_level_LOO.csv"), row.names = FALSE)
if (length(loading_var_results) > 0) {
  lv_df <- as.data.frame(lv_mat)
  lv_df$moderator <- rownames(lv_df)
  write.csv(lv_df, file.path(output_path, "item_level_loading_variance.csv"), row.names = FALSE)
}

saveRDS(list(full_result = full_res, loo_results = loo_results, loo_df = loo_df,
             loading_variance = loading_var_results, mean_var_per_item = mean_var_per_item),
        file.path(rds_path, "item_level_analysis.rds"))

cat("\nGenerating figures...\n")

p1 <- ggplot(loo_df, aes(x = reorder(item_removed, dl_change), y = dl_change)) +
  geom_col(aes(fill = dl_change > 0), alpha = 0.85) +
  geom_hline(yintercept = 0, color = "gray30") +
  coord_flip() +
  scale_fill_manual(values = c("TRUE" = "#27ae60", "FALSE" = "#e74c3c"),
                    labels = c("Hinders", "Contributes"), name = "Role") +
  labs(x = "Item Removed", y = expression(paste(Delta, ell, " Change (full - reduced)")),
       title = "Leave-One-Out: Item Contribution to Stabilization",
       subtitle = "Positive = removing item DECREASES stabilization (item contributes)") +
  theme_minimal() + theme(plot.title = element_text(size = 13, face = "bold"),
                           legend.position = "bottom")

if (length(loading_var_results) > 0) {
  mv_df <- data.frame(item = names(mean_var_per_item), variance = mean_var_per_item)
  p2 <- ggplot(mv_df, aes(x = reorder(item, variance), y = variance)) +
    geom_col(fill = "#3498db", alpha = 0.85) + coord_flip() +
    labs(x = "Item", y = "Mean Loading Variance Across Groups",
         title = "Cross-Group Loading Variability per Item",
         subtitle = "Higher = more measurement non-invariance at item level") +
    theme_minimal() + theme(plot.title = element_text(size = 13, face = "bold"))
} else {
  p2 <- ggplot() + theme_void()
}

ggsave(file.path(figure_path, "item_LOO_contribution.png"), p1, width = 10, height = 6, dpi = 600)
ggsave(file.path(figure_path, "item_loading_variance.png"), p2, width = 10, height = 6, dpi = 600)

cat("Saved to Outputs/, RDS/, Figures/\n")
cat("\nDone.\n")
