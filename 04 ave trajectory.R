rm(list = ls())
suppressPackageStartupMessages({
  library(readxl)
  library(lavaan)
  library(dplyr)
  library(stringr)
  library(boot)
  library(foreach)
  library(doParallel)
})

set.seed(9186)

cat("AVE TRAJECTORY ANALYSIS - FULL SVT AT EVERY STEP\n")
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

cat("Loading data...\n")
data_full <- read_excel(file.path(data_path, "analizlik2.xlsx"))
data_full <- data_full %>% mutate(across(where(is.character), turkce_ascii))
data_full$BRIAN <- data_full$BRO_Toplam
cat("  N =", nrow(data_full), "\n")

if (!"YJIPAQ" %in% names(data_full)) {
  if (requireNamespace("bestNormalize", quietly = TRUE)) {
    ipaq_col <- grep("IPAQ|ToplamMET|FA_Toplam", names(data_full), value = TRUE)[1]
    data_full$YJIPAQ <- bestNormalize::yeojohnson(data_full[[ipaq_col]])$x.t
    cat("  YJIPAQ created from", ipaq_col, "\n")
  } else stop("bestNormalize needed for Yeo-Johnson")
}

cat("Loading trajectories...\n")
traj1 <- read_excel(file.path(data_path, "AVE_trajectory_detailedCORRECTED.xlsx"))
traj1 <- traj1[order(traj1$Step), ]
traj1 <- traj1[!is.na(traj1$RemovedGozlemNo) & traj1$RemovedGozlemNo != "", ]
traj1$RemovedGozlemNo <- as.numeric(traj1$RemovedGozlemNo)
cat("  Phase 1:", nrow(traj1), "steps\n")

traj2_file <- file.path(data_path, "AVE_trajectory_MAX.xlsx")
if (file.exists(traj2_file)) {
  traj2 <- read_excel(traj2_file)
  traj2 <- traj2[order(traj2$Iterasyon), ]
  gcol <- grep("GozlemNo|Son_Cikarilan", names(traj2), value = TRUE)[1]
  traj2$RemovedGozlemNo <- as.numeric(traj2[[gcol]])
  traj2 <- traj2[!is.na(traj2$RemovedGozlemNo) & traj2$RemovedGozlemNo != 0, ]
  cat("  Phase 2:", nrow(traj2), "steps\n")
  removal_ids <- c(traj1$RemovedGozlemNo, traj2$RemovedGozlemNo)
  if (any(duplicated(removal_ids))) {
    removal_ids <- unique(removal_ids)
    cat("  Duplicates removed\n")
  }
} else {
  cat("  Phase 2 not found, Phase 1 only\n")
  removal_ids <- traj1$RemovedGozlemNo
}

n_steps <- length(removal_ids)
n_phase1 <- nrow(traj1)
cat("  Total:", n_steps, "steps\n\n")

meq_items <- c("MEQ1","MEQ2","MEQ8","MEQ9","MEQ10","MEQ11","MEQ15","MEQ17","MEQ18","MEQ19")

CFA_MODEL <- '
  MEQ_Factor =~ MEQ1 + MEQ2 + MEQ8 + MEQ9 + MEQ10 + MEQ11 + MEQ15 + MEQ17 + MEQ18 + MEQ19
  MEQ1 ~~ MEQ2
  MEQ2 ~~ MEQ10
  MEQ11 ~~ MEQ15
  MEQ17 ~~ MEQ18
  MEQ2 ~~ MEQ8
'

MODEL_WITHOUT <- '
  MEQ_Factor =~ MEQ1 + MEQ2 + MEQ8 + MEQ9 + MEQ10 + MEQ11 + MEQ15 + MEQ17 + MEQ18 + MEQ19
  YJIPAQ ~ c*MEQ_Factor
  MEQ1 ~~ MEQ2
  MEQ2 ~~ MEQ10
  MEQ11 ~~ MEQ15
  MEQ17 ~~ MEQ18
  MEQ2 ~~ MEQ8
'

MODEL_WITH <- '
  MEQ_Factor =~ MEQ1 + MEQ2 + MEQ8 + MEQ9 + MEQ10 + MEQ11 + MEQ15 + MEQ17 + MEQ18 + MEQ19
  BRIAN ~ a*MEQ_Factor
  YJIPAQ ~ cprime*MEQ_Factor + b*BRIAN
  MEQ1 ~~ MEQ2
  MEQ2 ~~ MEQ10
  MEQ11 ~~ MEQ15
  MEQ17 ~~ MEQ18
  MEQ2 ~~ MEQ8
'

MODEL_MEQB <- '
  MEQ_Factor =~ MEQ1 + MEQ2 + MEQ8 + MEQ9 + MEQ10 + MEQ11 + MEQ15 + MEQ17 + MEQ18 + MEQ19
  BRIAN ~ meqb*MEQ_Factor
  MEQ1 ~~ MEQ2
  MEQ2 ~~ MEQ10
  MEQ11 ~~ MEQ15
  MEQ17 ~~ MEQ18
  MEQ2 ~~ MEQ8
'

MODERATORS <- c("YasGruplari","Cinsiyet","EgitimDuzeyi","MedeniDurum",
                "GelirDuzeyi","BKISiniflamasi","SigaraKullanma",
                "AlkolKullanma","KronikHastalikVarligi","TibbiBeslenmeTedavisiAlma")

run_full_svt <- function(dat, meq_items, cfa_model, model_without, model_with, model_meqb, moderators) {
  
  calc_ave <- function(d) {
    md <- na.omit(d[, meq_items])
    if (nrow(md) < 100) return(NA)
    tryCatch({
      fit <- cfa(cfa_model, data = md, estimator = "MLR")
      if (lavInspect(fit, "converged")) {
        lds <- standardizedSolution(fit)[standardizedSolution(fit)$op == "=~", "est.std"]
        mean(lds^2)
      } else NA
    }, error = function(e) NA)
  }
  
  ave_val <- calc_ave(dat)
  
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
    gd_abs <- function(f1, f2, idx) { v1<-as.numeric(f1[idx]); v2<-as.numeric(f2[idx]); if(is.na(v1)||is.na(v2)) NA else abs(v1-v2) }
    gd_signed <- function(f1, f2, idx) { v1<-as.numeric(f1[idx]); v2<-as.numeric(f2[idx]); if(is.na(v1)||is.na(v2)) NA else v1-v2 }
    
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
    
    mi_details[[mod]] <- list(
      dc_m=dc_m, dt_m=dt_m, dr_m=dr_m, ds_m=ds_m,
      dc_s=dc_s, dt_s=dt_s, dr_s=dr_s, ds_s=ds_s,
      dc_t=dc_t, dt_t=dt_t, dr_t=dr_t, ds_t=ds_t,
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
    nc <- min(max(pmax(0,md$dc_m), pmax(0,md$dc_s), pmax(0,md$dc_t)) / 0.010, 1)
    nt <- min(max(pmax(0,md$dt_m), pmax(0,md$dt_s), pmax(0,md$dt_t)) / 0.010, 1)
    nr <- min(max(pmax(0,md$dr_m), pmax(0,md$dr_s), pmax(0,md$dr_t)) / 0.015, 1)
    ns <- min(max(pmax(0,md$ds_m), pmax(0,md$ds_s), pmax(0,md$ds_t)) / 0.030, 1)
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
    
    c_paths <- c(); cprime_paths <- c(); meqb_paths <- c()
    group_ns <- c(); group_names <- c()
    
    for (g in groups) {
      gd <- d[d[[mod]] == g, ]
      if (nrow(gd) < 30) next
      f0 <- tryCatch(sem(model_without, data=gd, std.lv=TRUE, estimator="MLR", missing="fiml"), error=function(e) NULL)
      f1 <- tryCatch(sem(model_with, data=gd, std.lv=TRUE, estimator="MLR", missing="fiml"), error=function(e) NULL)
      fb <- tryCatch(sem(model_meqb, data=gd, std.lv=TRUE, estimator="MLR", missing="fiml"), error=function(e) NULL)
      if (!is.null(f0) && !is.null(f1) && !is.null(fb)) {
        s0 <- standardizedSolution(f0); s1 <- standardizedSolution(f1); sb <- standardizedSolution(fb)
        cp <- s0[s0$label=="c","est.std"][1]
        cpp <- s1[s1$label=="cprime","est.std"][1]
        mbp <- sb[sb$label=="meqb","est.std"][1]
        if (!is.na(cp) && !is.na(cpp) && !is.na(mbp)) {
          c_paths <- c(c_paths, cp); cprime_paths <- c(cprime_paths, cpp); meqb_paths <- c(meqb_paths, mbp)
          group_ns <- c(group_ns, nrow(gd)); group_names <- c(group_names, g)
        }
      }
    }
    
    if (length(c_paths) < 2) { stab_results[[mod]] <- list(delta_log=NA); next }
    
    wmn <- function(x,w) weighted.mean(x,w)
    wsd <- function(x,w) sqrt(sum(w*(x-wmn(x,w))^2)/sum(w))
    mc <- wmn(c_paths,group_ns); mcp <- wmn(cprime_paths,group_ns); mmb <- wmn(meqb_paths,group_ns)
    sc <- wsd(c_paths,group_ns); scp <- wsd(cprime_paths,group_ns); smb <- wsd(meqb_paths,group_ns)
    em <- 0.0001
    cv0 <- (sc/max(abs(mc),em))*100; cv1 <- (scp/max(abs(mcp),em))*100; cvmb <- (smb/max(abs(mmb),em))*100
    el <- 1e-8
    dl <- log(cv0+el) - log(cv1+el)
    pct <- 100*(1-exp(-dl))
    ig <- sapply(seq_along(c_paths), function(i) ifelse(abs(cprime_paths[i]-mcp)<abs(c_paths[i]-mc),1,0))
    ocr <- mean(ig)
    eo <- 1e-8
    dlm <- log(abs(mcp)+eo) - log(abs(mc)+eo)
    dls <- log(sc+eo) - log(scp+eo)
    os <- max(0,dlm)/(max(0,dlm)+max(0,dls)+eo)
    
    stab_results[[mod]] <- list(delta_log=dl, cv_before=cv0, cv_after=cv1, cv_meqb=cvmb,
                                pct_reduction=pct, ocr=ocr, os=os, n_groups=length(c_paths), total_n=sum(group_ns))
  }
  
  dvec <- sapply(stab_results, function(x) x$delta_log)
  dvec <- dvec[!is.na(dvec)]
  if (length(dvec) == 0) {
    return(list(ave=ave_val, n=nrow(dat), p_boot=NA, p_binom=NA, mean_dl=NA,
                d_cohen=NA, weighted_os=NA, n_positive=NA, n_total=0, decision="No data",
                n_mi_violated=NA, n_mi_satisfied=NA, mi_results=mi_results, stab_results=stab_results))
  }
  
  wts <- sapply(names(dvec), function(m) stab_results[[m]]$total_n)
  ocrs <- sapply(names(dvec), function(m) stab_results[[m]]$ocr)
  oss <- sapply(names(dvec), function(m) stab_results[[m]]$os)
  
  bfn <- function(d, idx) { w <- wts[idx]/sum(wts[idx]); sum(d[idx]*w) }
  br <- tryCatch(boot(dvec, bfn, R=1000), error=function(e) NULL)
  if (!is.null(br)) {
    mdl <- br$t0; sedl <- sd(br$t); zs <- mdl/sedl
    pb <- pnorm(zs, lower.tail=FALSE)
    ci <- tryCatch(boot.ci(br, type="perc")$percent[4:5], error=function(e) c(NA,NA))
  } else {
    mdl <- weighted.mean(dvec, wts); pb <- NA; ci <- c(NA,NA)
  }
  
  wvar <- sum(wts*(dvec-mdl)^2)/sum(wts)
  dc <- if(wvar>0) mdl/sqrt(wvar) else NA
  
  c2p <- (dvec > 0) | (ocrs >= 0.5)
  so <- sum(c2p); mt <- length(c2p)
  pbn <- binom.test(so, mt, p=0.5, alternative="greater")$p.value
  
  wos <- weighted.mean(oss, wts)
  bpass <- (!is.na(pb) && pb < 0.05 && mdl > 0)
  bnpass <- (pbn < 0.05)
  is_stab <- bpass && bnpass
  
  if (is_stab) {
    if (wos < 0.3) mech <- "Type A"
    else if (wos > 0.7) mech <- "Type B"
    else mech <- "Type AB"
  } else mech <- "None"
  
  all_inv <- sapply(mi_results, function(x) x$is_invariant)
  
  list(ave=ave_val, n=nrow(dat), p_boot=pb, p_binom=pbn, mean_dl=mdl,
       d_cohen=dc, ci_lower=ci[1], ci_upper=ci[2], weighted_os=wos,
       n_positive=sum(dvec>0), n_total=mt, decision=ifelse(is_stab,"Stabilizer","Not"),
       mechanism=mech, n_mi_violated=sum(!all_inv), n_mi_satisfied=sum(all_inv),
       mi_results=mi_results, stab_results=stab_results, weights=weights)
}

n_cores <- 14
cat("Parallel setup:", n_cores, "cores\n")
cl <- makeCluster(n_cores)
registerDoParallel(cl)

clusterExport(cl, c("data_full","removal_ids","meq_items",
                    "CFA_MODEL","MODEL_WITHOUT","MODEL_WITH","MODEL_MEQB","MODERATORS","run_full_svt"))
clusterEvalQ(cl, {
  suppressPackageStartupMessages({ library(lavaan); library(dplyr); library(boot) })
})

cat("Starting full trajectory computation...\n")
cat("Steps:", n_steps, "| Cores:", n_cores, "\n\n")
t_start <- Sys.time()

batch_size <- 25
n_batches <- ceiling(n_steps / batch_size)
all_results <- list()

for (batch in 1:n_batches) {
  bs <- (batch-1)*batch_size + 1
  be <- min(batch*batch_size, n_steps)
  steps <- bs:be
  
  cat(sprintf("[Batch %d/%d] Steps %d-%d...\n", batch, n_batches, bs, be))
  
  batch_res <- foreach(step = steps, .packages = c("lavaan","dplyr","boot"),
                       .errorhandling = "pass") %dopar% {
                         ids_rm <- removal_ids[1:step]
                         step_data <- data_full[!data_full$GozlemNo %in% ids_rm, ]
                         res <- tryCatch(
                           run_full_svt(step_data, meq_items, CFA_MODEL, MODEL_WITHOUT, MODEL_WITH, MODEL_MEQB, MODERATORS),
                           error = function(e) list(ave=NA, n=nrow(step_data), p_boot=NA, p_binom=NA,
                                                    mean_dl=NA, d_cohen=NA, weighted_os=NA, n_positive=NA, n_total=0,
                                                    decision="Error", error=conditionMessage(e))
                         )
                         res$step <- step
                         res
                       }
  
  all_results <- c(all_results, batch_res)
  saveRDS(all_results, file.path(rds_path, "AVE_trajectory_checkpoint.rds"))
  
  elapsed <- difftime(Sys.time(), t_start, units="mins")
  pct <- be/n_steps*100
  eta <- as.numeric(elapsed)/(be/n_steps) - as.numeric(elapsed)
  cat(sprintf("  Done: %d/%d (%.0f%%) | %.1f min elapsed | ETA: %.1f min\n", be, n_steps, pct, as.numeric(elapsed), eta))
}

stopCluster(cl)
cat(sprintf("\nTotal time: %.1f hours\n", as.numeric(difftime(Sys.time(), t_start, units="hours"))))

cat("Computing step 0 (full data)...\n")
res0 <- run_full_svt(data_full, meq_items, CFA_MODEL, MODEL_WITHOUT, MODEL_WITH, MODEL_MEQB, MODERATORS)
res0$step <- 0

cat("Assembling results...\n")

rows <- list()
rows[[1]] <- data.frame(step=0, n_removed=0, n_remaining=nrow(data_full), phase=0,
                        ave=res0$ave, p_boot=res0$p_boot, p_binom=res0$p_binom, mean_dl=res0$mean_dl,
                        d_cohen=res0$d_cohen, ci_lower=res0$ci_lower, ci_upper=res0$ci_upper,
                        weighted_os=res0$weighted_os, n_positive=res0$n_positive, n_total=res0$n_total,
                        decision=res0$decision, mechanism=res0$mechanism,
                        n_mi_violated=res0$n_mi_violated, n_mi_satisfied=res0$n_mi_satisfied,
                        stringsAsFactors=FALSE)

safe <- function(x, d=NA) if(is.null(x)) d else x

for (i in seq_along(all_results)) {
  r <- all_results[[i]]
  if (is.null(r$step)) next
  rows[[length(rows)+1]] <- data.frame(
    step=r$step, n_removed=r$step, n_remaining=safe(r$n),
    phase=ifelse(r$step <= n_phase1, 1, 2),
    ave=safe(r$ave), p_boot=safe(r$p_boot), p_binom=safe(r$p_binom),
    mean_dl=safe(r$mean_dl), d_cohen=safe(r$d_cohen),
    ci_lower=safe(r$ci_lower), ci_upper=safe(r$ci_upper),
    weighted_os=safe(r$weighted_os), n_positive=safe(r$n_positive),
    n_total=safe(r$n_total), decision=safe(r$decision,"Error"),
    mechanism=safe(r$mechanism,"None"),
    n_mi_violated=safe(r$n_mi_violated), n_mi_satisfied=safe(r$n_mi_satisfied),
    stringsAsFactors=FALSE)
}

summary_df <- do.call(rbind, rows)
summary_df <- summary_df[order(summary_df$step), ]

mod_rows <- list()
for (i in seq_along(all_results)) {
  r <- all_results[[i]]
  if (is.null(r$stab_results)) next
  for (mn in names(r$stab_results)) {
    mr <- r$stab_results[[mn]]
    mi <- if (!is.null(r$mi_results[[mn]])) r$mi_results[[mn]] else list(mi_score=NA, is_invariant=NA)
    mod_rows[[length(mod_rows)+1]] <- data.frame(
      step=r$step, ave=safe(r$ave), moderator=mn,
      mi_score=safe(mi$mi_score), mi_satisfied=safe(mi$is_invariant),
      delta_log=safe(mr$delta_log), cv_before=safe(mr$cv_before),
      cv_after=safe(mr$cv_after), pct_reduction=safe(mr$pct_reduction),
      ocr=safe(mr$ocr), os=safe(mr$os), stringsAsFactors=FALSE)
  }
}
mod_df <- do.call(rbind, mod_rows)
mod_df <- mod_df[order(mod_df$step), ]

saveRDS(list(summary=summary_df, moderator_detail=mod_df, step0=res0,
             all_results=all_results,
             metadata=list(n_steps=n_steps, n_phase1=n_phase1, n_cores=n_cores,
                           cfa_model="5-cov",
                           total_hours=as.numeric(difftime(Sys.time(),t_start,units="hours")),
                           date=Sys.time())),
        file.path(rds_path, "AVE_trajectory_SVT.rds"))

write.csv(summary_df, file.path(output_path, "AVE_trajectory_summary.csv"), row.names=FALSE)
write.csv(mod_df, file.path(output_path, "AVE_trajectory_moderator_detail.csv"), row.names=FALSE)

cat("\nSaved:\n")
cat("  RDS/AVE_trajectory_SVT.rds\n")
cat("  Outputs/AVE_trajectory_summary.csv\n")
cat("  Outputs/AVE_trajectory_moderator_detail.csv\n")

v <- summary_df[!is.na(summary_df$ave), ]
cat(sprintf("\nAVE:  %.4f -> %.4f\n", min(v$ave,na.rm=T), max(v$ave,na.rm=T)))
cat(sprintf("N:    %d -> %d\n", max(v$n_remaining,na.rm=T), min(v$n_remaining,na.rm=T)))

cat("\nMilestones:\n")
for (ms in c(0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60)) {
  cl <- v[which.min(abs(v$ave - ms)), ]
  cat(sprintf("  AVE~%.2f (N=%d): Dl=%.3f, p_boot=%.4f, p_binom=%.4f, d=%.3f, OS=%.3f, %d/%d pos, %s\n",
              cl$ave, cl$n_remaining, cl$mean_dl, cl$p_boot, cl$p_binom,
              cl$d_cohen, cl$weighted_os, cl$n_positive, cl$n_total, cl$decision))
}

cat("\nGenerating figures...\n")
if (requireNamespace("ggplot2", quietly=TRUE)) {
  library(ggplot2)
  library(gridExtra)
  
  vv <- summary_df[!is.na(summary_df$mean_dl) & !is.na(summary_df$ave), ]
  
  p1 <- ggplot(vv, aes(x=ave, y=mean_dl)) +
    geom_point(aes(color=decision), alpha=0.4, size=0.8) +
    geom_smooth(method="loess", span=0.3, color="#e74c3c", se=TRUE, alpha=0.2) +
    geom_vline(xintercept=0.50, linetype="dashed", color="gray40") +
    scale_color_manual(values=c("Stabilizer"="#27ae60","Not"="#e74c3c","Error"="gray","No data"="gray")) +
    labs(x="AVE", y=expression(paste("Weighted Mean ",Delta,ell)),
         title="Stabilization Magnitude vs Measurement Quality") +
    theme_minimal() + theme(legend.position="bottom")
  
  p2 <- ggplot(vv, aes(x=ave, y=p_boot)) +
    geom_point(alpha=0.3, size=0.6, color="#2c3e50") +
    geom_hline(yintercept=0.05, linetype="dashed", color="red") +
    scale_y_log10() +
    labs(x="AVE", y="p_boot (log scale)", title="Bootstrap p-value Across AVE") +
    theme_minimal()
  
  p3 <- ggplot(vv, aes(x=ave, y=d_cohen)) +
    geom_point(alpha=0.3, size=0.6, color="#8e44ad") +
    geom_smooth(method="loess", span=0.3, color="#8e44ad", se=TRUE, alpha=0.2) +
    geom_hline(yintercept=c(0.2,0.5,0.8), linetype="dotted", color="gray60") +
    labs(x="AVE", y="Cohen's d", title="Effect Size Across AVE") +
    theme_minimal()
  
  p4 <- ggplot(vv, aes(x=ave, y=n_positive/n_total*100)) +
    geom_point(alpha=0.3, size=0.6, color="#27ae60") +
    geom_smooth(method="loess", span=0.3, color="#27ae60", se=TRUE, alpha=0.2) +
    geom_hline(yintercept=50, linetype="dashed", color="gray60") +
    labs(x="AVE", y="% Positive Moderators", title="Directional Consistency") +
    theme_minimal()
  
  ggsave(file.path(figure_path,"traj_Dl_vs_AVE.png"), p1, width=10, height=6, dpi=600)
  ggsave(file.path(figure_path,"traj_pboot_vs_AVE.png"), p2, width=10, height=6, dpi=600)
  ggsave(file.path(figure_path,"traj_d_vs_AVE.png"), p3, width=10, height=6, dpi=600)
  ggsave(file.path(figure_path,"traj_consistency_vs_AVE.png"), p4, width=10, height=6, dpi=600)
  
  comb <- grid.arrange(p1,p2,p3,p4, ncol=2, nrow=2,
                       top=grid::textGrob("SVT Full Trajectory: Stabilization as a Function of Measurement Quality",
                                          gp=grid::gpar(fontsize=14,fontface="bold")))
  ggsave(file.path(figure_path,"traj_combined_4panel.png"), comb, width=16, height=12, dpi=600)
  cat("  Figures saved\n")
}

cat("\nDone.\n")
cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")