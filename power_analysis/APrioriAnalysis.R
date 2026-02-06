# ===============================================================
# Power planner for my RM project
# What this script does:
#   - BETWEEN-sets: compare 3 independent cohorts on a single scalar (e.g., the 2×2 interaction contrast)
#   - WITHIN-set: plan a 2-level within effect inside one cohort (e.g., “true RM”, Visibility, Path, or the 2×2)
#   - MIXED: check whether a within effect (df1=1) differs by Set (i.e., Set × Within)
#
# Notes to future-me:
# - All effect sizes are Cohen’s f for the *exact* hypothesis I care about (already collapsed as I’ll analyze).
# - I don’t know if sphericity holds, so I’ll show both ε = 1.00 and a conservative ε = 0.75.
# - If I plan multiple tests, I’ll use Bonferroni in planning via m_family; Holm/Hochberg later in analysis.
#
# Dependency: ggplot2 (auto-installed if missing)
# ===============================================================

# ---- 0) Setup: keep it lightweight ----
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
suppressPackageStartupMessages(library(ggplot2))

# ---- 1) User parameters (edit these, but keep the structure) ----

## Global settings (same as before)
alpha_nominal <- 0.05        # familywise alpha before multiplicity
target_power  <- 0.80        # desired power level
m_family      <- 1           # Bonferroni family size (1 = single primary test)
stopifnot(m_family >= 1)     # sanity check: don’t accidentally divide by 0
alpha_per_test <- alpha_nominal / m_family

## Design: I’m planning three separate cohorts (e.g., Up / Down / Horizontal)
k_sets <- 3

## Effect sizes (Cohen's f) I’m using right now (pilot-informed, adjust when I have better estimates)
# A) BETWEEN-sets on a scalar per person (e.g., subject-level 2×2 interaction contrast)
f_between_sets_on_contrast <- 0.33

# B) WITHIN a single set (df1=1): could be “true RM”, Visibility, Path, or the 2×2 within interaction
f_within_single_set <- 0.40  # from pilot Geom×Vis; good medium–large anchor

# C) MIXED: Set (between, 3) × Within (df1=1): moderation size I expect across sets
f_mixed_set_by_within <- 0.33

## Plot ranges for the x-axes (per-cohort n). Solvers expand internally; this is just for visuals.
n_plot_min <- 6
n_plot_max <- 60

## Sphericity sensitivity: show both “assume sphericity” and “more conservative” cases
sensitivity_eps <- c(1.00, 0.75)

# ---- 2) Helper functions (same logic as before; just clearer comments) ----
eta2p_to_f <- function(eta2p) sqrt(eta2p / (1 - eta2p))
f_to_eta2p <- function(f)     f^2 / (1 + f^2)

# ----- BETWEEN-SUBJECTS (balanced one-way ANOVA; k groups) -----
# Compare a single scalar per participant across sets.
# Noncentrality: lambda = N_total * f^2 = n_per * k * f^2
power_between_oneway <- function(n_per, k, f, alpha = 0.05) {
  if (n_per < 2) return(NA_real_)
  df1 <- k - 1
  df2 <- k * (n_per - 1)
  lambda <- n_per * k * f^2
  Fcrit  <- qf(1 - alpha, df1, df2)
  1 - pf(Fcrit, df1, df2, ncp = lambda)
}

# Solve minimal integer n_per to hit target power (auto-expands upper bound if needed)
solve_n_per_between <- function(k, f, target = 0.80, alpha = 0.05, cap = 10000) {
  g <- function(n) power_between_oneway(n, k, f, alpha) - target
  lo <- 2; hi <- 12
  while (is.na(g(hi)) || g(hi) < 0) {
    lo <- hi; hi <- min(cap, ceiling(hi * 1.5))
    if (hi >= cap && g(hi) < 0) stop("Target power not reached (between). Increase f or cap.")
  }
  ceiling(uniroot(g, lower = lo + 1e-6, upper = hi, tol = 1e-6)$root)
}

# ----- WITHIN-SUBJECTS (single set, 2-level effect: df1 = 1) -----
# Under sphericity: df2 = (N - 1); GG: df1' = eps * 1; df2' = eps * (N - 1)
# Noncentrality: lambda = N * f^2
power_within_twolevel <- function(N, f, alpha = 0.05, eps = 1.0) {
  if (N <= 2) return(NA_real_)
  df1p <- eps * 1
  df2p <- eps * (N - 1)
  lambda <- N * f^2
  Fcrit  <- qf(1 - alpha, df1p, df2p)
  1 - pf(Fcrit, df1p, df2p, ncp = lambda)
}

# Solve per-cohort N for the within effect (treat per-cohort n as the within-subjects N)
solve_n_per_within <- function(f, target = 0.80, alpha = 0.05, eps = 1.0,
                               lower = 4, upper_start = 40, cap = 5000) {
  g <- function(N) power_within_twolevel(N, f, alpha, eps) - target
  lo <- lower; hi <- upper_start
  while (is.na(g(hi)) || g(hi) < 0) {
    lo <- hi; hi <- min(cap, ceiling(hi * 1.5))
    if (hi >= cap && g(hi) < 0) stop("Target power not reached (within). Increase f or cap.")
  }
  ceiling(uniroot(g, lower = lo + 1e-6, upper = hi, tol = 1e-6)$root)
}

# ----- MIXED: Set (between, k) × Within (df1=1) -----
# Numerator df: df1 = (k - 1) * 1 = k - 1
# Denominator df (GG): df2' = eps * (N_total - k)
# Noncentrality: lambda = N_total * df1 * f^2
power_mixed_set_within <- function(n_per, k, f, alpha = 0.05, eps = 1.0) {
  if (n_per < 2) return(NA_real_)
  Ntot <- n_per * k
  df1  <- (k - 1)
  df1p <- eps * df1
  df2p <- eps * (Ntot - k)
  lambda <- Ntot * df1 * f^2
  Fcrit  <- qf(1 - alpha, df1p, df2p)
  1 - pf(Fcrit, df1p, df2p, ncp = lambda)
}

# Solve per-cohort n for the mixed Set × Within test
solve_n_per_mixed <- function(k, f, target = 0.80, alpha = 0.05, eps = 1.0, cap = 10000) {
  g <- function(n) power_mixed_set_within(n, k, f, alpha, eps) - target
  lo <- 2; hi <- 12
  while (is.na(g(hi)) || g(hi) < 0) {
    lo <- hi; hi <- min(cap, ceiling(hi * 1.5))
    if (hi >= cap && g(hi) < 0) stop("Target power not reached (mixed). Increase f or cap.")
  }
  ceiling(uniroot(g, lower = lo + 1e-6, upper = hi, tol = 1e-6)$root)
}

# ----- Tiny plotting helper: same visuals, tidier title subtitles -----
plot_power_curve <- function(n_seq, power_seq, n_mark, title, subtitle, xlab = "Per-cohort sample size (n)") {
  df <- data.frame(n = n_seq, Power = power_seq)
  ggplot(df, aes(x = n, y = Power)) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 1.6) +
    geom_hline(yintercept = 0.80, linetype = "dashed") +
    geom_hline(yintercept = 0.90, linetype = "dashed") +
    geom_vline(xintercept = n_mark, linetype = "dotted") +
    annotate("text", x = max(n_seq) - 2, y = 0.80, label = "80% power", vjust = -0.6) +
    annotate("text", x = max(n_seq) - 2, y = 0.90, label = "90% power", vjust = -0.6) +
    annotate("text", x = n_mark, y = 0.55, label = paste0("n = ", n_mark), hjust = -0.05) +
    labs(title = title, subtitle = subtitle, x = xlab, y = "Power (1 - β)") +
    scale_y_continuous(limits = c(0.4, 1.0), breaks = seq(0.4, 1.0, 0.1)) +
    scale_x_continuous(breaks = seq(min(n_seq), max(n_seq), by = 5)) +
    theme_minimal(base_size = 14)
}

# ---- 3) Run the three planning branches, plus ε-sensitivity where it matters ----

## A) BETWEEN-sets (k=3) on the subject-level 2×2 interaction contrast
n_per_between <- solve_n_per_between(k_sets, f_between_sets_on_contrast, target_power, alpha_per_test)
p_between     <- power_between_oneway(n_per_between, k_sets, f_between_sets_on_contrast, alpha_per_test)

cat("\n=== BETWEEN (k=3) on interaction contrast ===\n")
cat(sprintf("f = %.3f (η_p^2 ≈ %.3f), α%s = %.3f\n",
            f_between_sets_on_contrast, f_to_eta2p(f_between_sets_on_contrast),
            if (m_family > 1) "'" else "", if (m_family > 1) alpha_per_test else alpha_nominal))
cat(sprintf("n per cohort = %d (Total N = %d), achieved power = %.3f\n",
            n_per_between, n_per_between * k_sets, p_between))

n_seq <- seq(min(n_plot_min, n_per_between - 20), max(n_plot_max, n_per_between + 40), by = 1)
pow_seq <- vapply(n_seq, function(n) power_between_oneway(n, k_sets, f_between_sets_on_contrast, alpha_per_test), numeric(1))
print(plot_power_curve(
  n_seq, pow_seq, n_per_between,
  title = "Power: BETWEEN-sets on 2×2 Interaction Contrast (1-way ANOVA)",
  subtitle = sprintf("f = %.2f (η_p^2 ≈ %.2f), k = %d, α%s = %.3f, target = %.2f",
                     f_between_sets_on_contrast, f_to_eta2p(f_between_sets_on_contrast),
                     k_sets, if (m_family > 1) "'" else "", if (m_family > 1) alpha_per_test else alpha_nominal,
                     target_power)
))

## B) WITHIN a single set (df1=1): I’ll show both ε = 1.00 and 0.75
cat("\n=== WITHIN a single set (df1=1) ===\n")
for (eps_within in sensitivity_eps) {
  n_per_within <- solve_n_per_within(f_within_single_set, target_power, alpha_per_test, eps_within)
  p_within     <- power_within_twolevel(n_per_within, f_within_single_set, alpha_per_test, eps_within)
  
  cat(sprintf("\n[ε = %.2f] f = %.3f (η_p^2 ≈ %.3f), α%s = %.3f\n",
              eps_within, f_within_single_set, f_to_eta2p(f_within_single_set),
              if (m_family > 1) "'" else "", if (m_family > 1) alpha_per_test else alpha_nominal))
  cat(sprintf("per-set N = %d, achieved power = %.3f\n", n_per_within, p_within))
  
  n_seq <- seq(min(n_plot_min, n_per_within - 20), max(n_plot_max, n_per_within + 40), by = 1)
  pow_seq <- vapply(n_seq, function(n) power_within_twolevel(n, f_within_single_set, alpha_per_test, eps_within), numeric(1))
  print(plot_power_curve(
    n_seq, pow_seq, n_per_within,
    title = "Power: WITHIN a Single Set (2-level effect, df1=1)",
    subtitle = sprintf("f = %.2f (η_p^2 ≈ %.2f), α%s = %.3f, ε = %.2f, target = %.2f",
                       f_within_single_set, f_to_eta2p(f_within_single_set),
                       if (m_family > 1) "'" else "", if (m_family > 1) alpha_per_test else alpha_nominal,
                       eps_within, target_power),
    xlab = "Per-cohort N (within-subjects in that set)"
  ))
}

## C) MIXED: Set × Within (df1 = k-1): again, show both ε values
cat("\n=== MIXED: Set × Within (df1 = k-1) ===\n")
for (eps_mixed in sensitivity_eps) {
  n_per_mixed <- solve_n_per_mixed(k_sets, f_mixed_set_by_within, target_power, alpha_per_test, eps_mixed)
  p_mixed     <- power_mixed_set_within(n_per_mixed, k_sets, f_mixed_set_by_within, alpha_per_test, eps_mixed)
  
  cat(sprintf("\n[ε = %.2f] f = %.3f (η_p^2 ≈ %.3f), α%s = %.3f\n",
              eps_mixed, f_mixed_set_by_within, f_to_eta2p(f_mixed_set_by_within),
              if (m_family > 1) "'" else "", if (m_family > 1) alpha_per_test else alpha_nominal))
  cat(sprintf("n per cohort = %d (Total N = %d), achieved power = %.3f\n",
              n_per_mixed, n_per_mixed * k_sets, p_mixed))
  
  n_seq <- seq(min(n_plot_min, n_per_mixed - 20), max(n_plot_max, n_per_mixed + 40), by = 1)
  pow_seq <- vapply(n_seq, function(n) power_mixed_set_within(n, k_sets, f_mixed_set_by_within, alpha_per_test, eps_mixed), numeric(1))
  print(plot_power_curve(
    n_seq, pow_seq, n_per_mixed,
    title = "Power: MIXED (Set × Within)",
    subtitle = sprintf("f = %.2f (η_p^2 ≈ %.2f), k = %d, α%s = %.3f, ε = %.2f, target = %.2f",
                       f_mixed_set_by_within, f_to_eta2p(f_mixed_set_by_within),
                       k_sets, if (m_family > 1) "'" else "", if (m_family > 1) alpha_per_test else alpha_nominal,
                       eps_mixed, target_power)
  ))
}

cat("\nNOTES:\n")
cat("• I print results for ε = 1.00 (assume sphericity) and ε = 0.75 (more conservative GG).\n")
cat("• If I collapse direction, I should estimate f on the *collapsed* metric and plug that here.\n")
cat("• BETWEEN-sets on the subject-level 2×2 interaction contrast equals the Set×A×B interaction algebraically,\n")
cat("  but this contrast route avoids messy mixed-model df details in planning.\n")
cat("• Using Holm/Hochberg in analysis is usually better; Bonferroni here is just a safe planning stand-in via m_family.\n\n")
