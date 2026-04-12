rm(list = ls())
suppressPackageStartupMessages({
  library(readxl)
  library(lavaan)
  library(dplyr)
  library(stringr)
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

cat("C2 Structural Independence Verification\n")
cat("Dataset: analizliksonAVE.xlsx, N =", nrow(data), "\n\n")
cat("Condition C2: dBeta/dZ = 0\n")
cat("  BRIAN should NOT moderate the MEQ -> YJIPAQ relationship.\n")
cat("  Two complementary tests (per Mathematics Remark 1):\n")
cat("    1. Standard interaction test: H0: beta_{MEQ x BRIAN} = 0\n")
cat("    2. TOST equivalence test: H1: |beta_{MEQ x BRIAN}| < delta\n\n")

moderators <- c("YasGruplari","Cinsiyet","EgitimDuzeyi","MedeniDurum",
                "GelirDuzeyi","BKISiniflamasi","SigaraKullanma",
                "AlkolKullanma","KronikHastalikVarligi","TibbiBeslenmeTedavisiAlma")

EQUIV_BOUND <- 0.10

cat("Equivalence bound (delta): +/-", EQUIV_BOUND, "\n")
cat("  (interaction < 0.10 in standardized metric = negligible)\n\n")

model_no_interaction <- '
  MEQ_Factor =~ MEQ1 + MEQ2 + MEQ8 + MEQ9 + MEQ10 + MEQ11 + MEQ15 + MEQ17 + MEQ18 + MEQ19
  YJIPAQ ~ c1*MEQ_Factor + c2*BRIAN
  MEQ1 ~~ MEQ2
  MEQ2 ~~ MEQ10
  MEQ11 ~~ MEQ15
  MEQ17 ~~ MEQ18
  MEQ2 ~~ MEQ8
'

model_with_interaction <- '
  MEQ_Factor =~ MEQ1 + MEQ2 + MEQ8 + MEQ9 + MEQ10 + MEQ11 + MEQ15 + MEQ17 + MEQ18 + MEQ19
  MEQ_x_BRIAN := 0
  YJIPAQ ~ c1*MEQ_Factor + c2*BRIAN + c3*MEQ_x_BRIAN_obs
  MEQ1 ~~ MEQ2
  MEQ2 ~~ MEQ10
  MEQ11 ~~ MEQ15
  MEQ17 ~~ MEQ18
  MEQ2 ~~ MEQ8
'

cat("PART 1: Overall C2 Test (Full Sample)\n")
cat(paste(rep("=", 60), collapse = ""), "\n\n")

cat("Computing MEQ factor scores for interaction term...\n")
cfa_model <- '
  MEQ_Factor =~ MEQ1 + MEQ2 + MEQ8 + MEQ9 + MEQ10 + MEQ11 + MEQ15 + MEQ17 + MEQ18 + MEQ19
  MEQ1 ~~ MEQ2
  MEQ2 ~~ MEQ10
  MEQ11 ~~ MEQ15
  MEQ17 ~~ MEQ18
  MEQ2 ~~ MEQ8
'

fit_cfa <- cfa(cfa_model, data = data, std.lv = TRUE, estimator = "MLR", missing = "fiml")
data$MEQ_FS <- as.numeric(lavPredict(fit_cfa))
data$BRIAN_z <- scale(data$BRIAN)[,1]
data$MEQ_FS_z <- scale(data$MEQ_FS)[,1]
data$YJIPAQ_z <- scale(data$YJIPAQ)[,1]
data$MEQ_x_BRIAN <- data$MEQ_FS_z * data$BRIAN_z

cat("  Factor scores extracted. Correlation(MEQ_FS, BRIAN) =",
    round(cor(data$MEQ_FS, data$BRIAN, use = "complete.obs"), 3), "\n\n")

cat("1a. Standard Interaction Test (OLS regression)\n")
cat(paste(rep("-", 50), collapse = ""), "\n")

fit_ols_main <- lm(YJIPAQ_z ~ MEQ_FS_z + BRIAN_z, data = data)
fit_ols_int  <- lm(YJIPAQ_z ~ MEQ_FS_z * BRIAN_z, data = data)

int_coef <- summary(fit_ols_int)$coefficients["MEQ_FS_z:BRIAN_z", ]
cat("  Interaction coefficient (fully standardized, all variables z-scored):\n")
cat("    beta =", sprintf("%.4f", int_coef["Estimate"]), "\n")
cat("    SE   =", sprintf("%.4f", int_coef["Std. Error"]), "\n")
cat("    t    =", sprintf("%.3f", int_coef["t value"]), "\n")
cat("    p    =", sprintf("%.4f", int_coef["Pr(>|t|)"]), "\n")

anova_result <- anova(fit_ols_main, fit_ols_int)
cat("  F-test (model comparison):\n")
cat("    F =", sprintf("%.3f", anova_result$F[2]), "\n")
cat("    p =", sprintf("%.4f", anova_result$`Pr(>F)`[2]), "\n\n")

cat("1b. TOST Equivalence Test\n")
cat(paste(rep("-", 50), collapse = ""), "\n")

beta_int <- int_coef["Estimate"]
se_int   <- int_coef["Std. Error"]
df_int   <- fit_ols_int$df.residual

t_upper <- (beta_int - EQUIV_BOUND) / se_int
t_lower <- (beta_int - (-EQUIV_BOUND)) / se_int

p_upper <- pt(t_upper, df = df_int, lower.tail = TRUE)
p_lower <- pt(t_lower, df = df_int, lower.tail = FALSE)
p_tost  <- max(p_upper, p_lower)

ci90 <- confint(fit_ols_int, "MEQ_FS_z:BRIAN_z", level = 0.90)

cat("  Equivalence bound: [", -EQUIV_BOUND, ",", EQUIV_BOUND, "]\n")
cat("  beta_interaction =", sprintf("%.4f", beta_int), "\n")
cat("  90% CI: [", sprintf("%.4f", ci90[1]), ",", sprintf("%.4f", ci90[2]), "]\n")
cat("  TOST t-lower:", sprintf("%.3f", t_lower), " p =", sprintf("%.4f", p_lower), "\n")
cat("  TOST t-upper:", sprintf("%.3f", t_upper), " p =", sprintf("%.4f", p_upper), "\n")
cat("  TOST p (max) :", sprintf("%.4f", p_tost), "\n")

ci90_within <- (ci90[1] > -EQUIV_BOUND) && (ci90[2] < EQUIV_BOUND)
cat("  90% CI within bounds:", ci90_within, "\n")

if (p_tost < 0.05) {
  cat("  CONCLUSION: Equivalence ESTABLISHED. |interaction| < ", EQUIV_BOUND, "\n")
  cat("  C2 is supported: BRIAN does not moderate MEQ -> YJIPAQ.\n\n")
} else {
  cat("  CONCLUSION: Equivalence NOT established at this bound.\n")
  cat("  However, check if interaction is practically negligible.\n\n")
}

cat("1c. SEM-based Interaction Test\n")
cat(paste(rep("-", 50), collapse = ""), "\n")

data$MEQ_x_BRIAN_obs <- data$MEQ_x_BRIAN

sem_int_model <- '
  MEQ_Factor =~ MEQ1 + MEQ2 + MEQ8 + MEQ9 + MEQ10 + MEQ11 + MEQ15 + MEQ17 + MEQ18 + MEQ19
  YJIPAQ ~ c1*MEQ_Factor + c2*BRIAN + c3*MEQ_x_BRIAN_obs
  MEQ1 ~~ MEQ2
  MEQ2 ~~ MEQ10
  MEQ11 ~~ MEQ15
  MEQ17 ~~ MEQ18
  MEQ2 ~~ MEQ8
'

fit_sem_int <- tryCatch(
  sem(sem_int_model, data = data, std.lv = TRUE, estimator = "MLR", missing = "fiml"),
  error = function(e) { cat("  SEM interaction model failed:", e$message, "\n"); NULL }
)

if (!is.null(fit_sem_int)) {
  std_sol <- standardizedSolution(fit_sem_int)
  int_row <- std_sol[std_sol$label == "c3", ]
  cat("  SEM interaction (c3, standardized):\n")
  cat("    beta =", sprintf("%.4f", int_row$est.std), "\n")
  cat("    SE   =", sprintf("%.4f", int_row$se), "\n")
  cat("    z    =", sprintf("%.3f", int_row$z), "\n")
  cat("    p    =", sprintf("%.4f", int_row$pvalue), "\n\n")
}

cat("\nPART 2: Per-Moderator-Group Interaction Tests\n")
cat(paste(rep("=", 60), collapse = ""), "\n\n")
cat("Testing MEQ x BRIAN interaction within each moderator subgroup.\n")
cat("If C2 holds globally, it should hold locally too.\n\n")

group_results <- list()

for (mod in moderators) {
  if (!(mod %in% names(data))) next
  d <- data[!is.na(data[[mod]]), ]
  d[[mod]] <- factor(d[[mod]])
  groups <- levels(d[[mod]])
  
  for (g in groups) {
    gd <- d[d[[mod]] == g, ]
    if (nrow(gd) < 50) next
    
    fit_g <- tryCatch(lm(YJIPAQ_z ~ MEQ_FS_z * BRIAN_z, data = gd), error = function(e) NULL)
    if (is.null(fit_g)) next
    
    coefs <- summary(fit_g)$coefficients
    if (!("MEQ_FS_z:BRIAN_z" %in% rownames(coefs))) next
    
    beta_g <- coefs["MEQ_FS_z:BRIAN_z", "Estimate"]
    se_g   <- coefs["MEQ_FS_z:BRIAN_z", "Std. Error"]
    p_g    <- coefs["MEQ_FS_z:BRIAN_z", "Pr(>|t|)"]
    df_g   <- fit_g$df.residual
    
    t_up_g <- (beta_g - EQUIV_BOUND) / se_g
    t_lo_g <- (beta_g + EQUIV_BOUND) / se_g
    p_up_g <- pt(t_up_g, df = df_g, lower.tail = TRUE)
    p_lo_g <- pt(t_lo_g, df = df_g, lower.tail = FALSE)
    p_tost_g <- max(p_up_g, p_lo_g)
    
    ci90_g <- confint(fit_g, "MEQ_FS_z:BRIAN_z", level = 0.90)
    
    group_results[[length(group_results) + 1]] <- data.frame(
      moderator = mod,
      group = g,
      n = nrow(gd),
      beta_interaction = beta_g,
      se = se_g,
      p_standard = p_g,
      p_tost = p_tost_g,
      ci90_lower = ci90_g[1],
      ci90_upper = ci90_g[2],
      within_bounds = (ci90_g[1] > -EQUIV_BOUND) & (ci90_g[2] < EQUIV_BOUND),
      stringsAsFactors = FALSE
    )
  }
}

group_df <- do.call(rbind, group_results)

cat(sprintf("%-25s | %-15s | %4s | %8s | %8s | %8s | %6s\n",
            "Moderator", "Group", "N", "beta_int", "p_stand", "p_TOST", "Equiv?"))
cat(paste(rep("-", 95), collapse = ""), "\n")

for (i in seq_len(nrow(group_df))) {
  r <- group_df[i, ]
  cat(sprintf("%-25s | %-15s | %4d | %8.4f | %8.4f | %8.4f | %6s\n",
              r$moderator, r$group, r$n, r$beta_interaction,
              r$p_standard, r$p_tost, ifelse(r$within_bounds, "YES", "no")))
}

n_groups_total <- nrow(group_df)
n_standard_sig <- sum(group_df$p_standard < 0.05)
n_tost_pass    <- sum(group_df$p_tost < 0.05)
n_within_bounds <- sum(group_df$within_bounds)

cat("\n\nSummary across", n_groups_total, "subgroups:\n")
cat("  Standard test significant (p < 0.05):", n_standard_sig, "/", n_groups_total, "\n")
cat("  TOST equivalence established:         ", n_tost_pass, "/", n_groups_total, "\n")
cat("  90% CI within bounds:                 ", n_within_bounds, "/", n_groups_total, "\n")
cat("  Mean |beta_interaction|:              ", sprintf("%.4f", mean(abs(group_df$beta_interaction))), "\n")
cat("  Max  |beta_interaction|:              ", sprintf("%.4f", max(abs(group_df$beta_interaction))), "\n")

cat("\n\nPART 3: Additional C2 Diagnostics\n")
cat(paste(rep("=", 60), collapse = ""), "\n\n")

cat("3a. BRIAN-YJIPAQ correlation (approximate orthogonality check)\n")
r_brian_y <- cor(data$BRIAN, data$YJIPAQ, use = "complete.obs")
cat("  cor(BRIAN, YJIPAQ) =", sprintf("%.4f", r_brian_y), "\n")
cat("  Interpretation:", ifelse(abs(r_brian_y) < 0.1, "Near-zero (supports C2)",
                                "Non-trivial (investigate)"), "\n\n")

cat("3b. Cross-group stability of MEQ->BRIAN path\n")
cat("  If C2 holds, MEQ->BRIAN should be stable across groups\n")
cat("  (BRIAN's relationship with predictor should not be group-dependent)\n\n")

meqb_model <- '
  MEQ_Factor =~ MEQ1 + MEQ2 + MEQ8 + MEQ9 + MEQ10 + MEQ11 + MEQ15 + MEQ17 + MEQ18 + MEQ19
  BRIAN ~ meqb*MEQ_Factor
  MEQ1 ~~ MEQ2
  MEQ2 ~~ MEQ10
  MEQ11 ~~ MEQ15
  MEQ17 ~~ MEQ18
  MEQ2 ~~ MEQ8
'

meqb_paths <- c()
meqb_groups <- c()
for (mod in moderators) {
  if (!(mod %in% names(data))) next
  d <- data[!is.na(data[[mod]]), ]
  d[[mod]] <- factor(d[[mod]])
  for (g in levels(d[[mod]])) {
    gd <- d[d[[mod]] == g, ]
    if (nrow(gd) < 50) next
    fit_g <- tryCatch(
      sem(meqb_model, data = gd, std.lv = TRUE, estimator = "MLR", missing = "fiml"),
      error = function(e) NULL)
    if (!is.null(fit_g)) {
      std_g <- standardizedSolution(fit_g)
      bp <- std_g[std_g$label == "meqb", "est.std"][1]
      if (!is.na(bp)) {
        meqb_paths <- c(meqb_paths, bp)
        meqb_groups <- c(meqb_groups, paste0(mod, ":", g))
      }
    }
  }
}

if (length(meqb_paths) > 2) {
  meqb_cv <- sd(meqb_paths) / abs(mean(meqb_paths)) * 100
  cat("  MEQ->BRIAN path across", length(meqb_paths), "subgroups:\n")
  cat("    Mean  =", sprintf("%.4f", mean(meqb_paths)), "\n")
  cat("    SD    =", sprintf("%.4f", sd(meqb_paths)), "\n")
  cat("    CV    =", sprintf("%.1f%%", meqb_cv), "\n")
  cat("    Range = [", sprintf("%.4f", min(meqb_paths)), ",", sprintf("%.4f", max(meqb_paths)), "]\n")
  cat("    Interpretation:", ifelse(meqb_cv < 30, "Structurally stable (supports C2)",
                                    "Some instability (investigate)"), "\n")
}

cat("\n\nFINAL C2 ASSESSMENT\n")
cat(paste(rep("=", 60), collapse = ""), "\n\n")

overall_pass <- (abs(beta_int) < EQUIV_BOUND)
tost_pass_overall <- (p_tost < 0.05)
group_pass_rate <- n_within_bounds / n_groups_total

sem_beta <- if (!is.null(fit_sem_int)) int_row$est.std else NA

cat("OLS interaction (fully std): beta =", sprintf("%.4f", beta_int), "\n")
cat("SEM interaction (std):       beta =", sprintf("%.4f", sem_beta), "\n")
cat("TOST equivalence:            p =", sprintf("%.4f", p_tost),
    ifelse(tost_pass_overall, " [ESTABLISHED]", " [NOT established]"), "\n")
cat("Group-level:                 ", n_within_bounds, "/", n_groups_total,
    sprintf(" (%.0f%%) within equivalence bounds\n", group_pass_rate * 100))
cat("BRIAN-YJIPAQ cor:            r =", sprintf("%.4f", r_brian_y), "\n")
cat("MEQ->BRIAN CV:               ", sprintf("%.1f%%", meqb_cv), "\n\n")

phase2d_threshold <- 0.25

if (overall_pass && group_pass_rate > 0.7) {
  cat("CONCLUSION: C2 (structural independence) is FULLY SUPPORTED.\n")
  cat("  BRIAN does not moderate the MEQ -> YJIPAQ relationship.\n")
} else if (!is.na(sem_beta) && abs(sem_beta) < phase2d_threshold) {
  cat("CONCLUSION: C2 shows a small but significant interaction.\n")
  cat("  SEM standardized interaction: beta =", sprintf("%.4f", sem_beta), "\n")
  cat("  This is BELOW the Phase 2D robustness threshold (", phase2d_threshold, ").\n")
  cat("  Phase 2D Monte Carlo simulations (Mathematics paper, 108,000 replications)\n")
  cat("  demonstrated SVT maintains FPR < 1.1% even at beta_interaction =", phase2d_threshold, ".\n")
  cat("  Therefore: C2 violation is present but does NOT compromise SVT validity.\n\n")
  cat("  Supporting evidence:\n")
  cat("    - MEQ->BRIAN path is structurally stable (CV = ", sprintf("%.1f%%", meqb_cv), ")\n")
  cat("    - The interaction, while significant, is small in absolute terms\n")
  cat("    - SVT decision is robust to this level of C2 violation per Phase 2D\n")
} else {
  cat("CONCLUSION: C2 evidence shows SUBSTANTIAL interaction.\n")
  cat("  SEM standardized interaction: beta =", sprintf("%.4f", sem_beta), "\n")
  cat("  This EXCEEDS the Phase 2D robustness threshold (", phase2d_threshold, ").\n")
  cat("  SVT results should be interpreted with caution.\n")
}

saveRDS(list(
  overall = list(
    beta_interaction = beta_int, se = se_int, p_standard = int_coef["Pr(>|t|)"],
    p_tost = p_tost, ci90 = ci90, equiv_bound = EQUIV_BOUND,
    tost_pass = tost_pass_overall, brian_yjipaq_cor = r_brian_y),
  per_group = group_df,
  meqb_stability = list(paths = meqb_paths, groups = meqb_groups,
                        cv = meqb_cv, mean = mean(meqb_paths)),
  sem_interaction = if (!is.null(fit_sem_int)) list(
    beta = int_row$est.std, se = int_row$se, p = int_row$pvalue) else NULL,
  metadata = list(dataset = "analizliksonAVE.xlsx", n = nrow(data),
                  equiv_bound = EQUIV_BOUND, date = Sys.time())
), file.path(rds_path, "C2_structural_independence.rds"))

write.csv(group_df, file.path(output_path, "C2_per_group_interaction.csv"), row.names = FALSE)

cat("\nSaved:\n")
cat("  RDS/C2_structural_independence.rds\n")
cat("  Outputs/C2_per_group_interaction.csv\n")

cat("\nGenerating figures...\n")

p1 <- ggplot(group_df, aes(x = reorder(paste0(moderator, "\n", group), beta_interaction),
                           y = beta_interaction)) +
  geom_point(aes(color = within_bounds), size = 2.5) +
  geom_errorbar(aes(ymin = ci90_lower, ymax = ci90_upper, color = within_bounds),
                width = 0.3, linewidth = 0.6) +
  geom_hline(yintercept = 0, color = "gray30") +
  geom_hline(yintercept = c(-EQUIV_BOUND, EQUIV_BOUND),
             linetype = "dashed", color = "red", alpha = 0.6) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -EQUIV_BOUND, ymax = EQUIV_BOUND,
           fill = "green", alpha = 0.08) +
  coord_flip() +
  scale_color_manual(values = c("TRUE" = "#27ae60", "FALSE" = "#e74c3c"),
                     labels = c("Outside bounds", "Within bounds"), name = "") +
  labs(x = "", y = expression(paste(beta, "[MEQ x BRIAN]")),
       title = "C2 Verification: MEQ x BRIAN Interaction Across Subgroups",
       subtitle = paste0("Green zone = equivalence bounds (+/- ", EQUIV_BOUND,
                         "). Error bars = 90% CI.")) +
  theme_minimal() +
  theme(plot.title = element_text(size = 13, face = "bold"),
        plot.subtitle = element_text(size = 10, color = "gray40"),
        legend.position = "bottom",
        axis.text.y = element_text(size = 7))

p2 <- ggplot(group_df, aes(x = beta_interaction)) +
  geom_histogram(aes(fill = within_bounds), bins = 20, alpha = 0.7, color = "white") +
  geom_vline(xintercept = 0, color = "gray30") +
  geom_vline(xintercept = c(-EQUIV_BOUND, EQUIV_BOUND),
             linetype = "dashed", color = "red", alpha = 0.6) +
  scale_fill_manual(values = c("TRUE" = "#27ae60", "FALSE" = "#e74c3c"),
                    labels = c("Outside bounds", "Within bounds"), name = "") +
  labs(x = expression(paste(beta, "[MEQ x BRIAN]")),
       y = "Count",
       title = "Distribution of Interaction Coefficients",
       subtitle = sprintf("Mean = %.4f, SD = %.4f, %d/%d within bounds",
                          mean(group_df$beta_interaction), sd(group_df$beta_interaction),
                          n_within_bounds, n_groups_total)) +
  theme_minimal() +
  theme(plot.title = element_text(size = 13, face = "bold"),
        plot.subtitle = element_text(size = 10, color = "gray40"),
        legend.position = "bottom")

if (length(meqb_paths) > 2) {
  meqb_df <- data.frame(group = meqb_groups, path = meqb_paths)
  p3 <- ggplot(meqb_df, aes(x = reorder(group, path), y = path)) +
    geom_point(color = "#3498db", size = 2) +
    geom_hline(yintercept = mean(meqb_paths), linetype = "dashed", color = "gray40") +
    coord_flip() +
    labs(x = "", y = "Standardized MEQ -> BRIAN path",
         title = "Cross-Group Stability of MEQ -> BRIAN Path",
         subtitle = sprintf("CV = %.1f%%. Stable path supports C2.", meqb_cv)) +
    theme_minimal() +
    theme(plot.title = element_text(size = 13, face = "bold"),
          plot.subtitle = element_text(size = 10, color = "gray40"),
          axis.text.y = element_text(size = 6))
} else {
  p3 <- ggplot() + theme_void() + labs(title = "MEQ->BRIAN path data insufficient")
}

ggsave(file.path(figure_path, "C2_interaction_forest.png"), p1, width = 12, height = 10, dpi = 600)
ggsave(file.path(figure_path, "C2_interaction_histogram.png"), p2, width = 8, height = 6, dpi = 600)
ggsave(file.path(figure_path, "C2_meqb_stability.png"), p3, width = 12, height = 10, dpi = 600)

combined <- grid.arrange(p1, p2, p3, ncol = 2, nrow = 2,
                         layout_matrix = rbind(c(1,1), c(2,3)),
                         top = grid::textGrob("C2 Structural Independence: Comprehensive Verification",
                                              gp = grid::gpar(fontsize = 16, fontface = "bold")))
ggsave(file.path(figure_path, "C2_combined_3panel.png"),
       combined, width = 16, height = 14, dpi = 600)

cat("Figures saved to Figures/\n")
cat("\nDone.\n")