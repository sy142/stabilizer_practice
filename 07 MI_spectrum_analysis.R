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

cat("MI-STABILIZATION SPECTRUM ANALYSIS\n")
cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

svt_file <- file.path(rds_path, "SVT_results_adaptive.rds")
if (file.exists(svt_file)) {
  cat("Loading pre-computed SVT results from 01...\n")
  svt_results <- readRDS(svt_file)
  summary_df <- svt_results$summary
  cat("  Loaded", nrow(summary_df), "moderator results\n\n")
} else {
  stop("Run 01_SVT_Empirical.R first to generate SVT_results_adaptive.rds")
}

cat("Per-moderator MI Score vs Stabilization:\n")
cat(sprintf("%-25s | %8s | %8s | %8s | %6s\n", "Moderator", "S_MI", "Dl", "Pct_Red", "MI_OK"))
cat(paste(rep("-", 80), collapse = ""), "\n")

for (i in seq_len(nrow(summary_df))) {
  r <- summary_df[i, ]
  cat(sprintf("%-25s | %8.3f | %8.3f | %+7.1f%% | %6s\n",
              r$moderator, r$mi_score, r$delta_log, r$pct_reduction,
              ifelse(r$mi_satisfied, "OK", "BAD")))
}

summary_df$mi_category <- cut(summary_df$mi_score,
                               breaks = c(-Inf, 0.3, 0.7, Inf),
                               labels = c("Low MI Violation", "Moderate MI Violation", "Severe MI Violation"))

cat("\n\nCategory Distribution:\n")
print(table(summary_df$mi_category))

cat("\nMean Stabilization by MI Category:\n")
cat_stats <- summary_df %>%
  group_by(mi_category) %>%
  summarise(
    n = n(),
    mean_dl = mean(delta_log, na.rm = TRUE),
    sd_dl = sd(delta_log, na.rm = TRUE),
    mean_pct = mean(pct_reduction, na.rm = TRUE),
    mean_ocr = mean(orientation_consistency, na.rm = TRUE),
    mean_os = mean(orientation_share, na.rm = TRUE),
    .groups = "drop"
  )

for (i in seq_len(nrow(cat_stats))) {
  cs <- cat_stats[i, ]
  cat(sprintf("  %s (n=%d): Mean Dl=%.3f (SD=%.3f), Mean %%Red=%.1f%%, OCR=%.2f, OS=%.3f\n",
              cs$mi_category, cs$n, cs$mean_dl, cs$sd_dl, cs$mean_pct, cs$mean_ocr, cs$mean_os))
}

cat("\nLinear and Quadratic Models:\n")
cat(paste(rep("=", 50), collapse = ""), "\n\n")

linear_model <- lm(delta_log ~ mi_score, data = summary_df)
quad_model <- lm(delta_log ~ mi_score + I(mi_score^2), data = summary_df)

cat("Linear model: Dl ~ S_MI\n")
cat(sprintf("  R-squared: %.4f\n", summary(linear_model)$r.squared))
cat(sprintf("  beta = %.3f, p = %.4f\n",
            coef(linear_model)[2], summary(linear_model)$coefficients[2, 4]))

cat("\nQuadratic model: Dl ~ S_MI + S_MI^2\n")
cat(sprintf("  R-squared: %.4f\n", summary(quad_model)$r.squared))
cat(sprintf("  beta_1 = %.3f, p = %.4f\n",
            coef(quad_model)[2], summary(quad_model)$coefficients[2, 4]))
cat(sprintf("  beta_2 = %.3f, p = %.4f\n",
            coef(quad_model)[3], summary(quad_model)$coefficients[3, 4]))

anova_result <- anova(linear_model, quad_model)
cat(sprintf("\n  Model comparison F = %.3f, p = %.4f\n",
            anova_result$F[2], anova_result$`Pr(>F)`[2]))

a_coef <- coef(quad_model)[3]
b_coef <- coef(quad_model)[2]
if (!is.na(a_coef) && a_coef != 0) {
  optimal_mi <- -b_coef / (2 * a_coef)
  if (optimal_mi > 0 && optimal_mi < 1) {
    predicted_max <- predict(quad_model, data.frame(mi_score = optimal_mi))
    cat(sprintf("\n  Optimal S_MI for max stabilization: %.3f\n", optimal_mi))
    cat(sprintf("  Predicted max Dl at optimum: %.3f\n", predicted_max))
    if (optimal_mi >= 0.3 && optimal_mi <= 0.7) {
      cat("  ** Optimal point in MODERATE violation range (confirms theory) **\n")
    }
  }
}

cat("\nStatistical Tests:\n")
cat(paste(rep("=", 50), collapse = ""), "\n\n")

if (nrow(cat_stats) >= 2) {
  kw <- kruskal.test(delta_log ~ mi_category, data = summary_df)
  cat(sprintf("Kruskal-Wallis: chi-sq = %.3f, p = %.4f\n", kw$statistic, kw$p.value))

  cats <- unique(summary_df$mi_category)
  if (length(cats) >= 2) {
    cat("\nPairwise Wilcoxon:\n")
    for (i in 1:(length(cats)-1)) {
      for (j in (i+1):length(cats)) {
        s1 <- summary_df$delta_log[summary_df$mi_category == cats[i]]
        s2 <- summary_df$delta_log[summary_df$mi_category == cats[j]]
        if (length(s1) > 0 && length(s2) > 0) {
          wt <- wilcox.test(s1, s2)
          cat(sprintf("  %s vs %s: p = %.4f %s\n",
                      cats[i], cats[j], wt$p.value, ifelse(wt$p.value < 0.05, "*", "")))
        }
      }
    }
  }
}

cat("\nMI-Complementarity Diagnostics:\n")
cat(paste(rep("=", 50), collapse = ""), "\n\n")

r_mi_dl <- cor(summary_df$mi_score, summary_df$delta_log, use = "complete.obs")
cat(sprintf("cor(S_MI, Dl) = %.3f\n", r_mi_dl))

if (sum(summary_df$mi_satisfied == FALSE, na.rm = TRUE) > 0) {
  dl_mi_bad <- summary_df$delta_log[!summary_df$mi_satisfied]
  dl_mi_good <- summary_df$delta_log[summary_df$mi_satisfied]
  cat(sprintf("Mean Dl (MI violated):  %.3f\n", mean(dl_mi_bad)))
  cat(sprintf("Mean Dl (MI satisfied): %.3f\n", mean(dl_mi_good)))

  if (length(dl_mi_bad) > 1 && length(dl_mi_good) > 1) {
    wt <- wilcox.test(dl_mi_bad, dl_mi_good)
    cat(sprintf("Wilcoxon test: p = %.4f\n", wt$p.value))
  }
}

cat("\nDestabilization Analysis (Severe MI):\n")
severe <- summary_df[summary_df$mi_score > 0.7, ]
if (nrow(severe) > 0) {
  n_destab <- sum(severe$delta_log < 0)
  cat(sprintf("  Severe violations: %d moderators\n", nrow(severe)))
  cat(sprintf("  Showing destabilization (Dl < 0): %d (%.1f%%)\n",
              n_destab, n_destab/nrow(severe)*100))
  cat(sprintf("  Mean Dl in severe: %.3f\n", mean(severe$delta_log)))
} else {
  cat("  No moderators with S_MI > 0.7\n")
}

write.csv(summary_df, file.path(output_path, "MI_spectrum_analysis.csv"), row.names = FALSE)
write.csv(cat_stats, file.path(output_path, "MI_spectrum_categories.csv"), row.names = FALSE)

saveRDS(list(summary = summary_df, category_stats = cat_stats,
             linear_model = linear_model, quadratic_model = quad_model,
             cor_mi_dl = r_mi_dl),
        file.path(rds_path, "MI_spectrum_analysis.rds"))

cat("\nGenerating figures...\n")

p1 <- ggplot(summary_df, aes(x = mi_score, y = delta_log)) +
  geom_point(aes(color = mi_category), size = 4, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, color = "blue", linetype = "dashed", alpha = 0.2) +
  geom_smooth(method = "lm", formula = y ~ x + I(x^2), se = TRUE, color = "red", alpha = 0.2) +
  geom_vline(xintercept = c(0.3, 0.7), linetype = "dotted", alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "solid", alpha = 0.3) +
  scale_color_manual(values = c("Low MI Violation" = "#27ae60",
                                "Moderate MI Violation" = "#f39c12",
                                "Severe MI Violation" = "#e74c3c")) +
  geom_text(aes(label = moderator), size = 2.5, vjust = -1, alpha = 0.7) +
  labs(title = "MI Violation Severity vs Stabilization Effect",
       subtitle = "Blue: linear | Red: quadratic. Dashed lines at S_MI = 0.3, 0.7",
       x = expression(paste("MI Score (", S[MI], ")")),
       y = expression(paste(Delta, ell, " (stabilization)")),
       color = "MI Category") +
  theme_minimal() +
  theme(plot.title = element_text(size = 14, face = "bold"),
        legend.position = "bottom")

p2 <- ggplot(summary_df, aes(x = mi_category, y = delta_log, fill = mi_category)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 2) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  scale_fill_manual(values = c("Low MI Violation" = "#27ae60",
                               "Moderate MI Violation" = "#f39c12",
                               "Severe MI Violation" = "#e74c3c")) +
  labs(title = "Stabilization Distribution by MI Category",
       x = "MI Violation Category",
       y = expression(paste(Delta, ell))) +
  theme_minimal() +
  theme(plot.title = element_text(size = 14, face = "bold"), legend.position = "none")

p3 <- ggplot(summary_df, aes(x = mi_score, y = pct_reduction)) +
  geom_point(aes(color = mi_category), size = 4, alpha = 0.8) +
  geom_smooth(method = "loess", se = TRUE, color = "gray40", alpha = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_manual(values = c("Low MI Violation" = "#27ae60",
                                "Moderate MI Violation" = "#f39c12",
                                "Severe MI Violation" = "#e74c3c")) +
  labs(title = "MI Score vs CV Reduction",
       x = expression(paste("MI Score (", S[MI], ")")),
       y = "CV Reduction (%)", color = "MI Category") +
  theme_minimal() +
  theme(plot.title = element_text(size = 14, face = "bold"), legend.position = "bottom")

ggsave(file.path(figure_path, "MI_spectrum_scatter.png"), p1, width = 10, height = 7, dpi = 600)
ggsave(file.path(figure_path, "MI_spectrum_boxplot.png"), p2, width = 8, height = 6, dpi = 600)
ggsave(file.path(figure_path, "MI_spectrum_cv_reduction.png"), p3, width = 10, height = 7, dpi = 600)

combined <- grid.arrange(p1, p2, p3, ncol = 2, nrow = 2,
  layout_matrix = rbind(c(1,1), c(2,3)),
  top = grid::textGrob("MI-Stabilization Spectrum Analysis",
                       gp = grid::gpar(fontsize = 16, fontface = "bold")))
ggsave(file.path(figure_path, "MI_spectrum_combined.png"), combined, width = 14, height = 12, dpi = 600)

cat("Saved to Outputs/, RDS/, Figures/\n")
cat("\nDone.\n")
