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

cat("COMPREHENSIVE ROBUSTNESS ANALYSIS\n")
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

cat("Dataset: N =", nrow(data), "\n\n")

N_ITER  <- 30
N_CORES <- 14

SAMPLE_SIZES <- seq(100, 1700, by = 100)

K_VALUES <- c(3, 5, 7, 10)

IMBALANCE_RATIOS <- c(0.5, 0.7, 0.9)

MISSING_PCTS <- c(0, 0.05, 0.10, 0.20, 0.30)

cat("Configuration:\n")
cat("  Iterations per test:   ", N_ITER, "\n")
cat("  Cores:                 ", N_CORES, "\n")
cat("  Sample sizes:          ", paste(SAMPLE_SIZES, collapse=", "), "\n")
cat("  Moderator counts (K):  ", paste(K_VALUES, collapse=", "), "\n")
cat("  Imbalance ratios:      ", paste(IMBALANCE_RATIOS, collapse=", "), "\n")
cat("  Missing data pcts:     ", paste(MISSING_PCTS*100, collapse=", "), "%\n\n")

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

all_moderators <- c("YasGruplari","Cinsiyet","EgitimDuzeyi","MedeniDurum",
                    "GelirDuzeyi","BKISiniflamasi","SigaraKullanma",
                    "AlkolKullanma","KronikHastalikVarligi","TibbiBeslenmeTedavisiAlma")

run_svt_full <- function(dat, moderators, cfa_model, model_without, model_with) {
  gd_abs <- function(f1,f2,idx) { v1<-as.numeric(f1[idx]);v2<-as.numeric(f2[idx]);if(is.na(v1)||is.na(v2)) NA else abs(v1-v2) }
  gd_s <- function(f1,f2,idx) { v1<-as.numeric(f1[idx]);v2<-as.numeric(f2[idx]);if(is.na(v1)||is.na(v2)) NA else v1-v2 }
  
  all_dm<-list();all_ds<-list();all_inv<-c();mi_det<-list()
  for(mod in moderators){
    if(!(mod%in%names(dat)))next;d<-dat[!is.na(dat[[mod]]),];d[[mod]]<-factor(d[[mod]])
    if(length(levels(d[[mod]]))<2)next
    fc<-tryCatch(cfa(cfa_model,data=d,group=mod,std.lv=T,estimator="MLR",missing="fiml"),error=function(e)NULL);if(is.null(fc))next
    fm<-tryCatch(cfa(cfa_model,data=d,group=mod,group.equal="loadings",std.lv=T,estimator="MLR",missing="fiml"),error=function(e)NULL);if(is.null(fm))next
    fs<-tryCatch(cfa(cfa_model,data=d,group=mod,group.equal=c("loadings","intercepts"),std.lv=T,estimator="MLR",missing="fiml"),error=function(e)NULL);if(is.null(fs))next
    ft<-tryCatch(cfa(cfa_model,data=d,group=mod,group.equal=c("loadings","intercepts","residuals"),std.lv=T,estimator="MLR",missing="fiml"),error=function(e)NULL);if(is.null(ft))next
    fmc<-fitmeasures(fc);fmm<-fitmeasures(fm);fms<-fitmeasures(fs);fmt<-fitmeasures(ft)
    dm<-c(CFI=gd_abs(fmc,fmm,"cfi.scaled"),TLI=gd_abs(fmc,fmm,"tli.scaled"),RMSEA=gd_abs(fmc,fmm,"rmsea.scaled"),SRMR=gd_abs(fmc,fmm,"srmr"))
    ds<-c(CFI=gd_abs(fmm,fms,"cfi.scaled"),TLI=gd_abs(fmm,fms,"tli.scaled"),RMSEA=gd_abs(fmm,fms,"rmsea.scaled"),SRMR=gd_abs(fmm,fms,"srmr"))
    dcm<-gd_s(fmc,fmm,"cfi.scaled");dtm<-gd_s(fmc,fmm,"tli.scaled");drm<- -gd_s(fmc,fmm,"rmsea.scaled");dsm<- -gd_s(fmc,fmm,"srmr")
    dcs<-gd_s(fmm,fms,"cfi.scaled");dts<-gd_s(fmm,fms,"tli.scaled");drs<- -gd_s(fmm,fms,"rmsea.scaled");dss<- -gd_s(fmm,fms,"srmr")
    dct<-gd_s(fms,fmt,"cfi.scaled");dtt<-gd_s(fms,fmt,"tli.scaled");drt<- -gd_s(fms,fmt,"rmsea.scaled");dst<- -gd_s(fms,fmt,"srmr")
    if(any(is.na(c(dm,ds,dct,dtt,drt,dst))))next
    all_dm[[mod]]<-dm;all_ds[[mod]]<-ds
    tc<-0.010;tt<-0.010;tr<-0.015;ts<-0.030
    inv<-(dcm<=tc&&dtm<=tt&&drm<=tr&&dsm<=ts&&dcs<=tc&&dts<=tt&&drs<=tr&&dss<=ts&&dct<=tc&&dtt<=tt&&drt<=tr&&dst<=ts)
    all_inv<-c(all_inv,inv);names(all_inv)[length(all_inv)]<-mod
    mi_det[[mod]]<-list(dcm=dcm,dtm=dtm,drm=drm,dsm=dsm,dcs=dcs,dts=dts,drs=drs,dss=dss,dct=dct,dtt=dtt,drt=drt,dst=dst,is_invariant=inv)
  }
  if(length(all_dm)>=3){
    dm_m<-do.call(rbind,all_dm);ds_m<-do.call(rbind,all_ds);dmax<-pmax(dm_m,ds_m)
    cr<-cor(dmax);nm<-cbind(CFI=pmin(dmax[,"CFI"]/0.01,1),TLI=pmin(dmax[,"TLI"]/0.01,1),RMSEA=pmin(dmax[,"RMSEA"]/0.015,1),SRMR=pmin(dmax[,"SRMR"]/0.03,1))
    cv<-apply(nm,2,function(x)sd(x)/(mean(x)+0.001));mc<-apply(abs(cr),1,function(x)mean(x[x<1]));rp<-1/(1+mc)
    dp<-numeric(4);names(dp)<-c("CFI","TLI","RMSEA","SRMR")
    for(idx in names(dp)){if(sum(all_inv)>0&&sum(!all_inv)>0)dp[idx]<-max(0,mean(nm[!all_inv,idx])-mean(nm[all_inv,idx]))else dp[idx]<-sd(nm[,idx])}
    vs<-1/(1+cv);rw<-dp*rp*vs;wt<-as.numeric(rw/sum(rw))
  }else wt<-c(0.4,0.1,0.3,0.2)
  
  mi_res<-list()
  for(mod in names(mi_det)){md<-mi_det[[mod]]
  nc<-min(max(pmax(0,md$dcm),pmax(0,md$dcs),pmax(0,md$dct))/0.01,1)
  nt<-min(max(pmax(0,md$dtm),pmax(0,md$dts),pmax(0,md$dtt))/0.01,1)
  nr<-min(max(pmax(0,md$drm),pmax(0,md$drs),pmax(0,md$drt))/0.015,1)
  ns<-min(max(pmax(0,md$dsm),pmax(0,md$dss),pmax(0,md$dst))/0.03,1)
  mi_res[[mod]]<-list(mi_score=wt[1]*nc+wt[2]*nt+wt[3]*nr+wt[4]*ns,is_invariant=md$is_invariant)
  }
  
  stab<-list()
  for(mod in moderators){
    if(!(mod%in%names(dat)))next;d<-dat[!is.na(dat[[mod]]),];d[[mod]]<-factor(d[[mod]])
    grps<-levels(d[[mod]]);if(length(grps)<2)next
    cp<-c();cpp<-c();gn<-c()
    for(g in grps){gd<-d[d[[mod]]==g,];if(nrow(gd)<30)next
    f0<-tryCatch(sem(model_without,data=gd,std.lv=T,estimator="MLR",missing="fiml"),error=function(e)NULL)
    f1<-tryCatch(sem(model_with,data=gd,std.lv=T,estimator="MLR",missing="fiml"),error=function(e)NULL)
    if(!is.null(f0)&&!is.null(f1)){s0<-standardizedSolution(f0);s1<-standardizedSolution(f1)
    c1<-s0[s0$label=="c","est.std"][1];c2<-s1[s1$label=="cprime","est.std"][1]
    if(!is.na(c1)&&!is.na(c2)){cp<-c(cp,c1);cpp<-c(cpp,c2);gn<-c(gn,nrow(gd))}}}
    if(length(cp)<2)next
    wmn<-function(x,w)weighted.mean(x,w);wsd<-function(x,w)sqrt(sum(w*(x-wmn(x,w))^2)/sum(w))
    mc<-wmn(cp,gn);mcp<-wmn(cpp,gn);sc<-wsd(cp,gn);scp<-wsd(cpp,gn)
    em<-0.0001;el<-1e-8
    cv0<-(sc/max(abs(mc),em))*100;cv1<-(scp/max(abs(mcp),em))*100
    dl<-log(cv0+el)-log(cv1+el)
    ig<-sapply(seq_along(cp),function(i)ifelse(abs(cpp[i]-mcp)<abs(cp[i]-mc),1,0))
    stab[[mod]]<-list(delta_log=dl,ocr=mean(ig),total_n=sum(gn))
  }
  
  dvec<-sapply(stab,function(x)x$delta_log);dvec<-dvec[!is.na(dvec)]
  if(length(dvec)<2)return(list(mean_dl=NA,p_boot=NA,p_binom=NA,d_cohen=NA,decision="Error",n_mi_sat=NA,n_mi_viol=NA))
  wts<-sapply(names(dvec),function(m)stab[[m]]$total_n)
  ocrs<-sapply(names(dvec),function(m)stab[[m]]$ocr)
  bfn<-function(d,idx){w<-wts[idx]/sum(wts[idx]);sum(d[idx]*w)}
  br<-tryCatch(boot(dvec,bfn,R=1000),error=function(e)NULL)
  if(!is.null(br)){mdl<-br$t0;pb<-pnorm(mdl/sd(br$t),lower.tail=FALSE)}else{mdl<-weighted.mean(dvec,wts);pb<-NA}
  wvar<-sum(wts*(dvec-mdl)^2)/sum(wts);dc<-if(wvar>0)mdl/sqrt(wvar)else NA
  c2p<-(dvec>0)|(ocrs>=0.5);pbn<-binom.test(sum(c2p),length(c2p),p=0.5,alternative="greater")$p.value
  bp<-(!is.na(pb)&&pb<0.05&&mdl>0);bn<-(pbn<0.05)
  ai<-sapply(mi_res,function(x)x$is_invariant)
  list(mean_dl=mdl,p_boot=pb,p_binom=pbn,d_cohen=dc,decision=ifelse(bp&&bn,"Stabilizer","Not"),
       n_positive=sum(dvec>0),n_total=length(dvec),n_mi_sat=sum(ai),n_mi_viol=sum(!ai))
}

cl <- makeCluster(N_CORES)
registerDoParallel(cl)
clusterExport(cl, c("data","all_moderators","cfa_model","model_without","model_with","run_svt_full"))
clusterEvalQ(cl, { suppressPackageStartupMessages({library(lavaan);library(dplyr);library(boot)}) })

t_start <- Sys.time()

cat("TEST 1/6: Sample Size Sensitivity\n")
cat(paste(rep("=", 50), collapse=""), "\n")

ss_results <- list()
for (n in SAMPLE_SIZES) {
  cat(sprintf("  N = %d... ", n))
  res_list <- foreach(i = 1:N_ITER, .packages=c("lavaan","dplyr","boot"), .errorhandling="pass") %dopar% {
    set.seed(9186 + i + n)
    idx <- sample(nrow(data), min(n, nrow(data)), replace = FALSE)
    d_sub <- data[idx, ]
    tryCatch(run_svt_full(d_sub, all_moderators, cfa_model, model_without, model_with),
             error = function(e) list(mean_dl=NA, p_boot=NA, p_binom=NA, d_cohen=NA, decision="Error"))
  }
  dls <- sapply(res_list, function(x) if(is.null(x$mean_dl)) NA else x$mean_dl)
  pbs <- sapply(res_list, function(x) if(is.null(x$p_boot)) NA else x$p_boot)
  decs <- sapply(res_list, function(x) if(is.null(x$decision)) "Error" else x$decision)
  valid <- !is.na(dls)
  ss_results[[as.character(n)]] <- list(n=n, mean_dl=mean(dls[valid]), sd_dl=sd(dls[valid]),
                                        prop_sig=mean(decs=="Stabilizer",na.rm=T), n_valid=sum(valid), n_iter=N_ITER)
  cat(sprintf("Dl=%.3f, Prop_Stab=%.2f (%d/%d valid)\n",
              mean(dls[valid]), mean(decs=="Stabilizer",na.rm=T), sum(valid), N_ITER))
}

cat("\nTEST 2/6: Number of Moderators\n")
cat(paste(rep("=", 50), collapse=""), "\n")

k_results <- list()
for (k in K_VALUES) {
  k_use <- min(k, length(all_moderators))
  cat(sprintf("  K = %d... ", k_use))
  res_list <- foreach(i = 1:N_ITER, .packages=c("lavaan","dplyr","boot"), .errorhandling="pass") %dopar% {
    set.seed(9186 + i*100 + k)
    sel_mods <- sample(all_moderators, k_use, replace = FALSE)
    tryCatch(run_svt_full(data, sel_mods, cfa_model, model_without, model_with),
             error = function(e) list(mean_dl=NA, decision="Error"))
  }
  dls <- sapply(res_list, function(x) if(is.null(x$mean_dl)) NA else x$mean_dl)
  decs <- sapply(res_list, function(x) if(is.null(x$decision)) "Error" else x$decision)
  valid <- !is.na(dls)
  k_results[[as.character(k_use)]] <- list(k=k_use, mean_dl=mean(dls[valid]), sd_dl=sd(dls[valid]),
                                           prop_sig=mean(decs=="Stabilizer",na.rm=T), n_valid=sum(valid))
  cat(sprintf("Dl=%.3f, Prop_Stab=%.2f (%d valid)\n",
              mean(dls[valid]), mean(decs=="Stabilizer",na.rm=T), sum(valid)))
}

cat("\nTEST 3/6: Group Imbalance\n")
cat(paste(rep("=", 50), collapse=""), "\n")

imb_results <- list()
test_mod <- "Cinsiyet"

for (ratio in IMBALANCE_RATIOS) {
  label <- paste0(round(ratio*100), "-", round((1-ratio)*100))
  cat(sprintf("  Ratio %s... ", label))
  res_list <- foreach(i = 1:N_ITER, .packages=c("lavaan","dplyr","boot"), .errorhandling="pass") %dopar% {
    set.seed(9186 + i*200 + round(ratio*100))
    d <- data[!is.na(data[[test_mod]]), ]
    d[[test_mod]] <- factor(d[[test_mod]])
    grps <- levels(d[[test_mod]])
    if (length(grps) < 2) return(list(mean_dl=NA, decision="Error"))
    g1 <- d[d[[test_mod]]==grps[1], ]; g2 <- d[d[[test_mod]]==grps[2], ]
    n_max <- min(nrow(g1), nrow(g2))
    n_total <- min(1000, floor(n_max / max(ratio, 1-ratio)))
    n1 <- round(n_total * ratio); n2 <- n_total - n1
    if (n1 < 30 || n2 < 30) return(list(mean_dl=NA, decision="Error"))
    s1 <- g1[sample(nrow(g1), n1, replace=FALSE), ]
    s2 <- g2[sample(nrow(g2), n2, replace=FALSE), ]
    d_imb <- rbind(s1, s2)
    tryCatch(run_svt_full(d_imb, all_moderators, cfa_model, model_without, model_with),
             error = function(e) list(mean_dl=NA, decision="Error"))
  }
  dls <- sapply(res_list, function(x) if(is.null(x$mean_dl)) NA else x$mean_dl)
  decs <- sapply(res_list, function(x) if(is.null(x$decision)) "Error" else x$decision)
  valid <- !is.na(dls)
  imb_results[[label]] <- list(ratio=label, mean_dl=mean(dls[valid]), sd_dl=sd(dls[valid]),
                               prop_sig=mean(decs=="Stabilizer",na.rm=T), n_valid=sum(valid))
  cat(sprintf("Dl=%.3f (%d valid)\n", mean(dls[valid]), sum(valid)))
}

cat("\nTEST 4/6: Missing Data Tolerance\n")
cat(paste(rep("=", 50), collapse=""), "\n")

meq_items <- c("MEQ1","MEQ2","MEQ8","MEQ9","MEQ10","MEQ11","MEQ15","MEQ17","MEQ18","MEQ19")
miss_results <- list()

for (pct in MISSING_PCTS) {
  cat(sprintf("  Missing %.0f%%... ", pct*100))
  res_list <- foreach(i = 1:N_ITER, .packages=c("lavaan","dplyr","boot"), .errorhandling="pass") %dopar% {
    set.seed(9186 + i + round(pct*1000))
    d_miss <- data
    if (pct > 0) {
      n_miss <- round(nrow(data) * pct)
      rows <- sample(nrow(data), n_miss)
      for (r in rows) {
        nv <- sample(1:2, 1)
        vars <- sample(c(meq_items, "BRIAN", "YJIPAQ"), nv)
        d_miss[r, vars] <- NA
      }
    }
    tryCatch(run_svt_full(d_miss, all_moderators, cfa_model, model_without, model_with),
             error = function(e) list(mean_dl=NA, decision="Error"))
  }
  dls <- sapply(res_list, function(x) if(is.null(x$mean_dl)) NA else x$mean_dl)
  decs <- sapply(res_list, function(x) if(is.null(x$decision)) "Error" else x$decision)
  valid <- !is.na(dls)
  miss_results[[as.character(pct)]] <- list(pct=pct, mean_dl=mean(dls[valid]), sd_dl=sd(dls[valid]),
                                            prop_sig=mean(decs=="Stabilizer",na.rm=T), n_valid=sum(valid))
  cat(sprintf("Dl=%.3f, Prop_Stab=%.2f (%d valid)\n",
              mean(dls[valid]), mean(decs=="Stabilizer",na.rm=T), sum(valid)))
}

cat("\nTEST 5/6: Bootstrap Stability\n")
cat(paste(rep("=", 50), collapse=""), "\n")

boot_results <- foreach(i = 1:N_ITER, .packages=c("lavaan","dplyr","boot"), .errorhandling="pass") %dopar% {
  set.seed(9186 + i*500)
  idx <- sample(nrow(data), nrow(data), replace = TRUE)
  d_boot <- data[idx, ]
  tryCatch(run_svt_full(d_boot, all_moderators, cfa_model, model_without, model_with),
           error = function(e) list(mean_dl=NA, decision="Error"))
}

boot_dls <- sapply(boot_results, function(x) if(is.null(x$mean_dl)) NA else x$mean_dl)
boot_decs <- sapply(boot_results, function(x) if(is.null(x$decision)) "Error" else x$decision)
boot_valid <- !is.na(boot_dls)
cat(sprintf("  Bootstrap: mean Dl=%.3f (SD=%.3f), Prop_Stab=%.2f (%d valid)\n",
            mean(boot_dls[boot_valid]), sd(boot_dls[boot_valid]),
            mean(boot_decs=="Stabilizer",na.rm=T), sum(boot_valid)))

if (sum(boot_valid) > 2) {
  boot_ci <- quantile(boot_dls[boot_valid], c(0.025, 0.975))
  cat(sprintf("  95%% CI: [%.3f, %.3f]\n", boot_ci[1], boot_ci[2]))
  cat(sprintf("  CV: %.1f%%\n", sd(boot_dls[boot_valid])/mean(boot_dls[boot_valid])*100))
}

cat("\nTEST 6/6: Outlier Sensitivity\n")
cat(paste(rep("=", 50), collapse=""), "\n")

brian_z <- abs(scale(data$BRIAN[!is.na(data$BRIAN)]))
yjipaq_z <- abs(scale(data$YJIPAQ[!is.na(data$YJIPAQ)]))
outlier_mask <- rep(FALSE, nrow(data))
if (!is.na(data$BRIAN[1])) outlier_mask <- outlier_mask | (abs(scale(data$BRIAN)) > 3)
if (!is.na(data$YJIPAQ[1])) outlier_mask <- outlier_mask | (abs(scale(data$YJIPAQ)) > 3)
outlier_mask[is.na(outlier_mask)] <- FALSE

n_outliers <- sum(outlier_mask)
cat(sprintf("  Identified %d outliers (|z| > 3)\n", n_outliers))

outlier_result <- NULL
if (n_outliers > 0 && n_outliers < nrow(data) * 0.1) {
  data_clean <- data[!outlier_mask, ]
  cat(sprintf("  Testing without outliers (N = %d)... ", nrow(data_clean)))
  outlier_result <- tryCatch(
    run_svt_full(data_clean, all_moderators, cfa_model, model_without, model_with),
    error = function(e) NULL)
  if (!is.null(outlier_result)) {
    cat(sprintf("Dl=%.3f, d=%.3f, %s\n", outlier_result$mean_dl, outlier_result$d_cohen, outlier_result$decision))
  } else cat("Failed\n")
}

stopCluster(cl)
total_time <- as.numeric(difftime(Sys.time(), t_start, units = "hours"))
cat(sprintf("\nTotal computation time: %.1f hours\n", total_time))

robustness_results <- list(
  sample_size = ss_results,
  n_moderators = k_results,
  group_imbalance = imb_results,
  missing_data = miss_results,
  bootstrap = list(dls = boot_dls[boot_valid], decs = boot_decs[boot_valid],
                   mean_dl = mean(boot_dls[boot_valid]), sd_dl = sd(boot_dls[boot_valid]),
                   prop_sig = mean(boot_decs[boot_valid]=="Stabilizer")),
  outliers = list(n_outliers = n_outliers, result = outlier_result),
  metadata = list(n_iter = N_ITER, n_cores = N_CORES, n_original = nrow(data),
                  total_hours = total_time, date = Sys.time())
)

saveRDS(robustness_results, file.path(rds_path, "robustness_analysis.rds"))

ss_df <- do.call(rbind, lapply(ss_results, function(x) data.frame(N=x$n, Mean_Dl=x$mean_dl,
                                                                  SD_Dl=x$sd_dl, Prop_Stabilizer=x$prop_sig, N_Valid=x$n_valid, stringsAsFactors=FALSE)))
k_df <- do.call(rbind, lapply(k_results, function(x) data.frame(K=x$k, Mean_Dl=x$mean_dl,
                                                                SD_Dl=x$sd_dl, Prop_Stabilizer=x$prop_sig, N_Valid=x$n_valid, stringsAsFactors=FALSE)))
imb_df <- do.call(rbind, lapply(imb_results, function(x) data.frame(Ratio=x$ratio, Mean_Dl=x$mean_dl,
                                                                    SD_Dl=x$sd_dl, Prop_Stabilizer=x$prop_sig, N_Valid=x$n_valid, stringsAsFactors=FALSE)))
miss_df <- do.call(rbind, lapply(miss_results, function(x) data.frame(Missing_Pct=x$pct*100,
                                                                      Mean_Dl=x$mean_dl, SD_Dl=x$sd_dl, Prop_Stabilizer=x$prop_sig, N_Valid=x$n_valid, stringsAsFactors=FALSE)))

write.csv(ss_df, file.path(output_path, "robustness_sample_size.csv"), row.names=FALSE)
write.csv(k_df, file.path(output_path, "robustness_moderator_count.csv"), row.names=FALSE)
write.csv(imb_df, file.path(output_path, "robustness_imbalance.csv"), row.names=FALSE)
write.csv(miss_df, file.path(output_path, "robustness_missing_data.csv"), row.names=FALSE)

cat("\nSaved:\n")
cat("  RDS/robustness_analysis.rds\n")
cat("  Outputs/robustness_*.csv\n")

cat("\nGenerating figures...\n")

p1 <- ggplot(ss_df, aes(x=N)) +
  geom_ribbon(aes(ymin=Mean_Dl-1.96*SD_Dl, ymax=Mean_Dl+1.96*SD_Dl), fill="#3498db", alpha=0.2) +
  geom_line(aes(y=Mean_Dl), color="#2c3e50", linewidth=1.2) +
  geom_point(aes(y=Mean_Dl), color="#e74c3c", size=3) +
  geom_hline(yintercept=0, linetype="dashed", color="gray40") +
  scale_x_continuous(breaks=seq(0,1700,200)) +
  labs(x="Sample Size", y=expression(paste("Mean ",Delta,ell)),
       title="Sample Size Sensitivity", subtitle="Shaded = 95% CI") +
  theme_minimal() + theme(plot.title=element_text(size=13,face="bold"))

p2 <- ggplot(ss_df, aes(x=N, y=Prop_Stabilizer)) +
  geom_area(fill="#27ae60", alpha=0.3) +
  geom_line(color="#27ae60", linewidth=1.2) +
  geom_point(color="#27ae60", size=3) +
  geom_hline(yintercept=0.80, linetype="dotted", color="gray60") +
  scale_x_continuous(breaks=seq(0,1700,200)) +
  scale_y_continuous(labels=scales::percent, limits=c(0,1)) +
  labs(x="Sample Size", y="Proportion Stabilizer Detected",
       title="Detection Rate Across Sample Sizes") +
  theme_minimal() + theme(plot.title=element_text(size=13,face="bold"))

p3 <- ggplot(k_df, aes(x=factor(K), y=Mean_Dl)) +
  geom_col(aes(fill=Mean_Dl), width=0.6, alpha=0.85, show.legend=FALSE) +
  geom_errorbar(aes(ymin=Mean_Dl-SD_Dl, ymax=Mean_Dl+SD_Dl), width=0.2) +
  geom_hline(yintercept=0, color="gray30") +
  scale_fill_gradient(low="#3498db", high="#e74c3c") +
  labs(x="Number of Moderators (K)", y=expression(paste("Mean ",Delta,ell)),
       title="Effect of Moderator Count") +
  theme_minimal() + theme(plot.title=element_text(size=13,face="bold"))

p4 <- ggplot(miss_df, aes(x=Missing_Pct, y=Mean_Dl)) +
  geom_ribbon(aes(ymin=Mean_Dl-SD_Dl, ymax=Mean_Dl+SD_Dl), fill="#1abc9c", alpha=0.2) +
  geom_line(color="#16a085", linewidth=1.2) +
  geom_point(color="#16a085", size=3) +
  labs(x="Missing Data (%)", y=expression(paste("Mean ",Delta,ell)),
       title="Missing Data Tolerance") +
  theme_minimal() + theme(plot.title=element_text(size=13,face="bold"))

if (sum(boot_valid) > 2) {
  bdf <- data.frame(dl = boot_dls[boot_valid])
  p5 <- ggplot(bdf, aes(x=dl)) +
    geom_histogram(aes(y=after_stat(density)), bins=20, fill="#3498db", alpha=0.7, color="white") +
    geom_density(color="#2c3e50", linewidth=1) +
    geom_vline(xintercept=mean(bdf$dl), linetype="dashed", color="#e74c3c", linewidth=1) +
    labs(x=expression(paste(Delta,ell)), y="Density",
         title="Bootstrap Distribution", subtitle=sprintf("Mean=%.3f, 95%% CI=[%.3f, %.3f]",
                                                          mean(bdf$dl), boot_ci[1], boot_ci[2])) +
    theme_minimal() + theme(plot.title=element_text(size=13,face="bold"))
} else {
  p5 <- ggplot() + theme_void() + labs(title="Bootstrap: insufficient data")
}

imb_df$Ratio <- factor(imb_df$Ratio, levels=c("50-50","70-30","90-10"))
p6 <- ggplot(imb_df, aes(x=Ratio, y=Mean_Dl)) +
  geom_col(fill="#8e44ad", alpha=0.85, width=0.5) +
  geom_errorbar(aes(ymin=Mean_Dl-SD_Dl, ymax=Mean_Dl+SD_Dl), width=0.15) +
  labs(x="Group Ratio", y=expression(paste("Mean ",Delta,ell)),
       title="Group Imbalance Effects") +
  theme_minimal() + theme(plot.title=element_text(size=13,face="bold"))

ggsave(file.path(figure_path,"rob_sample_size_dl.png"), p1, width=10, height=6, dpi=600)
ggsave(file.path(figure_path,"rob_sample_size_detection.png"), p2, width=10, height=6, dpi=600)
ggsave(file.path(figure_path,"rob_moderator_count.png"), p3, width=8, height=6, dpi=600)
ggsave(file.path(figure_path,"rob_missing_data.png"), p4, width=10, height=6, dpi=600)
ggsave(file.path(figure_path,"rob_bootstrap_dist.png"), p5, width=8, height=6, dpi=600)
ggsave(file.path(figure_path,"rob_imbalance.png"), p6, width=8, height=6, dpi=600)

comb <- grid.arrange(p1,p2,p3,p4,p5,p6, ncol=2, nrow=3,
                     top=grid::textGrob("SVT Comprehensive Robustness Analysis",
                                        gp=grid::gpar(fontsize=16,fontface="bold")))
ggsave(file.path(figure_path,"rob_combined_6panel.png"), comb, width=16, height=18, dpi=600)

cat("Figures saved\n\nDone.\n")
cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")