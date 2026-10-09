#!/usr/bin/env Rscript
# Year-to-year variability in Rao's Q (NIRv and all bands) across NEON sites
# Usage: Rscript Code/DataAnalysis/RaoQYearVariability.R [csv] [data_dir] [fig_dir]
#   csv      default "./Data/spectral_diversity_by_year.csv"
#   data_dir default "./Data"    (raoq_year_tests_R.csv, raoq_site_temporal_cv_R.csv)
#   fig_dir  default "./Figures" (all PNGs)

required_pkgs <- c("lme4", "patchwork", "readr", "dplyr", "tidyr", "ggplot2")
missing_pkgs <- required_pkgs[!required_pkgs %in% rownames(installed.packages())]
if (length(missing_pkgs) > 0) {
  install.packages(missing_pkgs, repos = "https://cloud.r-project.org", quiet = TRUE)
}

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr)
  library(ggplot2); library(patchwork); library(lme4)
})

args     <- commandArgs(trailingOnly = TRUE)
csv      <- if (length(args) >= 1) args[1] else "./Data/spectral_diversity_by_year.csv"
data_dir <- if (length(args) >= 2) args[2] else "./Data"
fig_dir  <- if (length(args) >= 3) args[3] else "./Figures"
n_perm <- 5000
set.seed(42)

dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir,  showWarnings = FALSE, recursive = TRUE)

if (!file.exists(csv)) stop("Input csv not found: ", csv)

metrics <- c(raoq_nirv = "Rao's Q (NIRv)", raoq_allbands = "Rao's Q (all bands)")
cols    <- c(raoq_nirv = "#2a9d8f", raoq_allbands = "#e76f51")

dat <- read_csv(csv, show_col_types = FALSE)

# ---- input validation --------------------------------------------------------
required_cols <- c("neon_site", "year", "raoq_nirv", "raoq_allbands")
missing_cols <- setdiff(required_cols, names(dat))
if (length(missing_cols) > 0) {
  stop("Missing required column(s): ", paste(missing_cols, collapse = ", "),
       "\nColumns present in file: ", paste(names(dat), collapse = ", "))
}

dat <- dat %>% mutate(year = as.integer(year))

status_col <- intersect(c("status", "ok", "error"), names(dat))
if (length(status_col) > 0) {
  status_col <- status_col[1]
  cat("\n=== Status column found:", status_col, "===\n")
  print(table(dat[[status_col]], useNA = "ifany"))
  if (status_col %in% c("status", "ok")) {
    success_vals <- c("success", "ok", "TRUE", TRUE, 1)
    dat <- dat %>% filter(.data[[status_col]] %in% success_vals)
  } else {
    dat <- dat %>% filter(is.na(.data[[status_col]]) | .data[[status_col]] == "")
  }
  cat("Filtered to successful rows using column '", status_col, "'\n", sep = "")
} else {
  cat("\nNo status/success column found (checked: status, ok, error). Proceeding without filtering.\n")
}

dups <- dat %>% count(neon_site, year) %>% filter(n > 1)
if (nrow(dups) > 0) {
  offending <- dat %>% semi_join(dups, by = c("neon_site", "year")) %>%
    arrange(neon_site, year)
  print(offending)
  stop("Duplicate neon_site x year rows found (ANOVA/LMM assume one value per site-year). ",
       "See offending rows printed above.")
}

cat("\n=== Post-filter summary ===\n")
cat("n site-years:", nrow(dat), "\n")
cat("n sites:", n_distinct(dat$neon_site), "\n")
cat("year range:", paste(range(dat$year, na.rm = TRUE), collapse = " - "), "\n\n")

fmt_p <- function(p) ifelse(p < 0.001, "p < 0.001", sprintf("p = %.3f", p))

# ---- F statistic for year effect, adjusted for site (fast, via QR) -------------
# Precompute design matrices/ranks/dfs once per metric; only y changes across permutations.
year_F_prep <- function(d) {
  Xr <- model.matrix(~ site, d)
  Xf <- model.matrix(~ site + factor(year), d)
  qr_r <- qr(Xr); qr_f <- qr(Xf)
  rank_r <- qr_r$rank; rank_f <- qr_f$rank
  df1 <- rank_f - rank_r
  df2 <- nrow(d) - rank_f
  list(Xr = Xr, Xf = Xf, df1 = df1, df2 = df2)
}

year_F_fast <- function(y, prep) {
  fr <- lm.fit(prep$Xr, y); ff <- lm.fit(prep$Xf, y)
  rss_r <- sum(fr$residuals^2); rss_f <- sum(ff$residuals^2)
  ((rss_r - rss_f) / prep$df1) / (rss_f / prep$df2)
}

safe_lmer <- function(formula, data, metric_label, model_label) {
  fit <- withCallingHandlers(
    lmer(formula, data = data, REML = FALSE),
    warning = function(w) {
      message(sprintf("[%s / %s] lmer warning: %s", metric_label, model_label, conditionMessage(w)))
      invokeRestart("muffleWarning")
    }
  )
  fit
}

analyse <- function(metric) {
  d <- dat %>% select(site = neon_site, year, y = all_of(metric)) %>%
    filter(!is.na(y)) %>% group_by(site) %>% filter(n() >= 2) %>% ungroup() %>%
    mutate(site = factor(site))

  # 1. two-way ANOVA (site, then year)
  an    <- anova(lm(y ~ site + factor(year), d))
  ss    <- an[["Sum Sq"]]; names(ss) <- rownames(an)
  F_yr  <- an["factor(year)", "F value"]

  # 2. permutation test: shuffle y within site
  prep  <- year_F_prep(d)
  F_obs <- year_F_fast(d$y, prep)
  idx   <- split(seq_len(nrow(d)), d$site)
  F_prm <- replicate(n_perm, {
    yp <- d$y
    for (i in idx) yp[i] <- d$y[i][sample.int(length(i))]
    year_F_fast(yp, prep)
  })
  p_perm <- (sum(F_prm >= F_obs) + 1) / (n_perm + 1)

  # 3. mixed model LRT (random intercept for site)
  m0 <- safe_lmer(y ~ 1 + (1 | site), d, metric, "null model (m0)")
  m1 <- safe_lmer(y ~ factor(year) + (1 | site), d, metric, "year model (m1)")
  lrt <- anova(m0, m1)
  vc  <- as.data.frame(VarCorr(m0))
  icc <- vc$vcov[1] / sum(vc$vcov)

  # site-adjusted year means + 95% CI
  d <- d %>% group_by(site) %>% mutate(y_c = y - mean(y) + mean(d$y)) %>% ungroup()
  ym <- d %>% group_by(year) %>%
    summarise(m = mean(y_c), s = sd(y_c), n = n(), .groups = "drop") %>%
    mutate(ci = ifelse(n > 1, qt(0.975, n - 1) * s / sqrt(n), NA_real_))

  # per-site temporal CV (>= 3 years)
  g <- d %>% group_by(site) %>%
    summarise(mean = mean(y), sd = sd(y), n = n(), .groups = "drop") %>%
    filter(n >= 3) %>% mutate(cv_pct = 100 * sd / mean)
  site_means <- d %>% group_by(site) %>% summarise(m = mean(y), .groups = "drop")

  res <- tibble(
    metric = metric, n_obs = nrow(d), n_sites = nlevels(d$site), n_years = n_distinct(d$year),
    F_year = F_yr, df1 = an["factor(year)", "Df"], df2 = an["Residuals", "Df"],
    p_anova = an["factor(year)", "Pr(>F)"], p_perm = p_perm,
    LRT = lrt$Chisq[2], p_lmm = lrt$`Pr(>Chisq)`[2], ICC_site = icc,
    pct_var_site = 100 * ss["site"] / sum(ss),
    pct_var_year = 100 * ss["factor(year)"] / sum(ss),
    pct_var_resid = 100 * ss["Residuals"] / sum(ss),
    median_within_site_cv = median(g$cv_pct),
    among_site_cv = 100 * sd(site_means$m) / mean(d$y))
  list(d = d, res = res, ym = ym, g = g)
}

out <- tryCatch(
  { o <- lapply(names(metrics), analyse); names(o) <- names(metrics); o },
  error = function(e) {
    message("Analysis step failed: ", conditionMessage(e))
    quit(status = 1)
  }
)
res <- bind_rows(lapply(out, `[[`, "res"))
cvs <- bind_rows(lapply(names(out), function(k) mutate(out[[k]]$g, metric = k)))
write_csv(res, file.path(data_dir, "raoq_year_tests_R.csv"))
write_csv(cvs, file.path(data_dir, "raoq_site_temporal_cv_R.csv"))

cat("\n=== Year-effect tests (site-adjusted) ===\n")
print(as.data.frame(t(res %>% mutate(across(where(is.numeric), ~ round(.x, 4))))),
      right = FALSE)
cat("\n=== Top 10 most variable sites (temporal CV, %) ===\n")
print(cvs %>% group_by(metric) %>% slice_max(cv_pct, n = 10) %>%
        select(metric, site, n_years = n, cv_pct) %>% mutate(cv_pct = round(cv_pct, 1)),
      n = 20)

# ---- plots ------------------------------------------------------------------------
theme_set(theme_classic(base_size = 11))

ggplot2_ge_35 <- packageVersion("ggplot2") >= "3.5.0"

label_annotation <- function(..., label.size = 0.2) {
  if (ggplot2_ge_35) {
    annotate("label", ..., linewidth = label.size)
  } else {
    annotate("label", ..., label.size = label.size)
  }
}

cv_legend_theme <- if (ggplot2_ge_35) {
  theme(legend.position = "inside", legend.position.inside = c(0.85, 0.85))
} else {
  theme(legend.position = c(0.85, 0.85))
}

traj_plot <- function(k) {
  o <- out[[k]]; r <- o$res
  lab <- sprintf("Year effect (site-adjusted)\nANOVA F(%d,%d) = %.2f, %s\nPermutation %s\nMixed model LRT %s",
                 r$df1, r$df2, r$F_year, fmt_p(r$p_anova), fmt_p(r$p_perm), fmt_p(r$p_lmm))
  ggplot() +
    geom_line(data = o$d, aes(year, y, group = site), colour = "grey75", linewidth = 0.3) +
    geom_point(data = o$d, aes(year, y), colour = "grey55", size = 0.8) +
    geom_errorbar(data = o$ym, aes(year, ymin = m - ci, ymax = m + ci),
                  colour = cols[k], width = 0.25, linewidth = 0.8, na.rm = TRUE) +
    geom_line(data = o$ym, aes(year, m), colour = cols[k], linewidth = 1.1) +
    geom_point(data = o$ym, aes(year, m), colour = cols[k], size = 3) +
    label_annotation(x = -Inf, y = Inf, label = lab, hjust = -0.05, vjust = 1.1,
             size = 2.6, fill = "white") +
    scale_x_continuous(breaks = sort(unique(o$d$year))) +
    scale_y_continuous(expand = expansion(mult = c(0.03, 0.32))) +
    labs(x = "Year", y = metrics[k], title = metrics[k],
         caption = "Grey = individual sites; colour = site-adjusted year mean +/- 95% CI") +
    theme(plot.title = element_text(face = "bold", colour = cols[k]),
          plot.caption = element_text(size = 7, colour = "grey40"),
          axis.text.x = element_text(size = 8))
}

for (k in names(metrics))
  ggsave(file.path(fig_dir, paste0(k, "_by_year_R.png")), traj_plot(k),
         width = 7, height = 4.8, dpi = 200)

pA <- traj_plot("raoq_nirv"); pB <- traj_plot("raoq_allbands")

vp <- res %>% select(metric, Site = pct_var_site, Year = pct_var_year, Residual = pct_var_resid) %>%
  pivot_longer(-metric, names_to = "component", values_to = "pct") %>%
  mutate(metric = factor(metrics[metric], levels = metrics),
         component = factor(component, levels = c("Residual", "Year", "Site")))
pC <- ggplot(vp, aes(metric, pct, fill = component)) +
  geom_col(colour = "white") +
  geom_text(data = filter(vp, pct > 4), aes(label = sprintf("%.0f%%", pct)),
            position = position_stack(vjust = 0.5), size = 3,
            colour = ifelse(filter(vp, pct > 4)$component == "Site", "white", "black")) +
  scale_fill_manual(values = c(Site = "#264653", Year = "#e9c46a", Residual = "#d9d9d9")) +
  scale_x_discrete(labels = c("NIRv", "All bands")) +
  labs(x = NULL, y = "% of total variance", fill = NULL, title = "Variance partition") +
  theme(plot.title = element_text(face = "bold"), legend.position = "bottom")

ord <- cvs %>% group_by(site) %>% summarise(m = mean(cv_pct), .groups = "drop") %>%
  arrange(desc(m)) %>% pull(site)
pD <- cvs %>% mutate(site = factor(site, levels = ord), metric = factor(metrics[metric], levels = metrics)) %>%
  ggplot(aes(site, cv_pct, fill = metric)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  scale_fill_manual(values = setNames(cols, metrics)) +
  labs(x = NULL, y = "Temporal CV across years (%)", fill = NULL,
       title = "Year-to-year variability within each site (sites with >= 3 years)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 6.5),
        plot.title = element_text(face = "bold")) +
  cv_legend_theme

verdict <- function(p) ifelse(p < 0.05, "significant", ifelse(p < 0.10, "marginal", "not significant"))
txt <- paste(c("Summary", "", unlist(lapply(names(metrics), function(k) {
  r <- filter(res, metric == k)
  c(metrics[k],
    sprintf("  Year effect: %s", verdict(r$p_perm)),
    sprintf("  perm. %s; LMM %s", fmt_p(r$p_perm), fmt_p(r$p_lmm)),
    sprintf("  Median within-site CV: %.0f%%", r$median_within_site_cv),
    sprintf("  Among-site CV: %.0f%%", r$among_site_cv), "")
})), sprintf("n = %d site-years, %d sites", res$n_obs[1], res$n_sites[1])), collapse = "\n")
pE <- ggplot() + annotate("text", x = 0, y = 1, label = txt, hjust = 0, vjust = 1,
                          family = "mono", size = 3) +
  xlim(0, 1) + ylim(0, 1) + theme_void()

summary_fig <- ((pA | pB | pC) / (pD | pE)) +
  plot_layout(heights = c(1.1, 1)) + plot_layout(widths = c(1, 1, 1)) +
  plot_annotation(
    title = "Year-to-year variation in spectral diversity (Rao's Q) across NEON sites",
    caption = "One value per site-year, so tests ask whether years differ after accounting for site.",
    tag_levels = "a",
    theme = theme(plot.title = element_text(face = "bold", size = 13)))
ggsave(file.path(fig_dir, "raoq_summary_figure_R.png"), summary_fig,
       width = 13, height = 9, dpi = 200)
cat("\nSaved figures to", fig_dir, "and CSVs to", data_dir, "\n")
