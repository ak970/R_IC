# Title: 02_Multivariable_Panel_Modeling.R
# Scope: Dynamic Annual State (Perennial vs Intermittent), Continuous Climate/Hydro Drivers
# Features: Panel Fixed-Effects (Within Estimator), Repeated Measures Correlation (rmcorr), Lagged Cross-Correlation
# It takes the panel dataset created by Script 1.
# You can paste or merge your actual annual climate/hydrology columns into panel_df.
# It runs the within-catchment Fixed Effects regression (plm), Repeated Measures correlations (rmcorr), and lagged responses.

required_pkgs <- c("dplyr", "tidyr", "ggplot2", "cowplot", "plm", "rmcorr")
new_pkgs <- required_pkgs[!(required_pkgs %in% installed.packages()[, "Package"])]
if (length(new_pkgs) > 0) install.packages(new_pkgs)

library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(plm)
library(rmcorr)

###################################
# 1. SETUP & DATA INGESTION
###################################
out_dir <- "F:/Peninsular India/Connectivity/IC/Plots/Multivariable_Panel"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# 1. Load the 31-year extracted IC panel data (created by Script 1)
panel_ic_path <- "F:/Peninsular India/Connectivity/IC/Plots/Temporal_1991_2021/Panel_IC_1991_2021_Extracted.csv"
panel_df      <- read.csv(panel_ic_path, stringsAsFactors = FALSE)
panel_df$gauge_id <- sprintf("%05d", as.numeric(panel_df$gauge_id))

# -------------------------------------------------------------------------
# 2. ANNUAL DYNAMIC COVARIATES SETUP (FILL IN YOUR DATA HERE)
# -------------------------------------------------------------------------
# In your project, you will load or merge your table containing:
#   - gauge_id
#   - year
#   - annual_flow_state : Binary dynamic variable: 0 = Perennial, 1 = Intermittent
#   - annual_precip     : Continuous dynamic variable (mm)
#   - annual_temp       : Continuous dynamic variable (°C)
#   - annual_aridity    : Continuous dynamic variable (-)
#   - annual_runoff     : Continuous dynamic variable (mm)

# --- PLACEHOLDER SYNTHESIS (Replace this block with your actual merged CSV) ---
cat("Synthesizing dynamic covariates structure (replace with your dynamic CSV)...\n")
set.seed(42)
panel_df$annual_flow_state <- ifelse(panel_df$Regime == "1", 
                                     rbinom(nrow(panel_df), 1, 0.75), 
                                     rbinom(nrow(panel_df), 1, 0.20))
panel_df$annual_precip     <- rnorm(nrow(panel_df), mean = 1100, sd = 250)
panel_df$annual_temp       <- rnorm(nrow(panel_df), mean = 26.5, sd = 1.5)
panel_df$annual_aridity    <- rnorm(nrow(panel_df), mean = 0.65, sd = 0.15)
# -----------------------------------------------------------------------------

# Ensure proper factoring
panel_df$annual_flow_state <- factor(panel_df$annual_flow_state, levels = c(0, 1), 
                                     labels = c("Perennial", "Intermittent"))


###################################
# 2. DYNAMIC BINARY STATE ANALYSIS (Perennial vs. Intermittent)
###################################
cat("\n--- Running Dynamic Annual Flow State Analysis ---\n")

# Compare catchments that experienced BOTH states during the 31 years
state_comparison <- panel_df %>%
  group_by(gauge_id, annual_flow_state) %>%
  summarize(mean_IC_state = mean(median_IC), .groups = "drop") %>%
  pivot_wider(names_from = annual_flow_state, values_from = mean_IC_state) %>%
  filter(!is.na(Perennial) & !is.na(Intermittent))

paired_state_test <- wilcox.test(state_comparison$Intermittent, state_comparison$Perennial, paired = TRUE)
mean_state_diff   <- mean(state_comparison$Intermittent - state_comparison$Perennial)

cat(sprintf("Paired Wilcoxon Test (Intermittent vs Perennial Years): p = %.4e %s\n", 
            paired_state_test$p.value, ifelse(paired_state_test$p.value < 0.05, "[SIGNIFICANT]", "[NOT SIG.]")))
cat(sprintf("Mean Connectivity Shift during Intermittent Years: %+.4f IC units\n", mean_state_diff))


###################################
# 3. PANEL FIXED-EFFECTS ("WITHIN") REGRESSION
###################################
cat("\n--- Running Panel Fixed-Effects (Within-Catchment) Regression ---\n")
# The Within-Estimator eliminates all static variables (elevation, geology, area)
# and evaluates strictly how annual climate anomalies perturb connectivity.

pdata <- pdata.frame(panel_df, index = c("gauge_id", "year"))

fe_model <- plm(
  median_IC ~ annual_flow_state + annual_precip + annual_temp + annual_aridity,
  data  = pdata,
  model = "within"
)

cat("Panel Fixed-Effects Model Summary:\n")
fe_summary <- summary(fe_model)
print(fe_summary)
write.csv(as.data.frame(fe_summary$coefficients), file.path(out_dir, "Panel_FixedEffects_Coefficients.csv"))


###################################
# 4. REPEATED MEASURES CORRELATION (rmcorr)
###################################
cat("\n--- Running Repeated Measures Correlation (rmcorr) ---\n")
# Evaluates common within-basin linear association between continuous climate drivers and IC

rmc_precip <- rmcorr::rmcorr(participant = gauge_id, measure1 = annual_precip, measure2 = median_IC, dataset = panel_df)
rmc_temp   <- rmcorr::rmcorr(participant = gauge_id, measure1 = annual_temp,   measure2 = median_IC, dataset = panel_df)
rmc_arid   <- rmcorr::rmcorr(participant = gauge_id, measure1 = annual_aridity,measure2 = median_IC, dataset = panel_df)

cat(sprintf("rmcorr (IC vs Precip) : r = %+.3f, p = %.3e\n", rmc_precip$r, rmc_precip$p))
cat(sprintf("rmcorr (IC vs Temp)   : r = %+.3f, p = %.3e\n", rmc_temp$r,   rmc_temp$p))
cat(sprintf("rmcorr (IC vs Aridity): r = %+.3f, p = %.3e\n", rmc_arid$r,   rmc_arid$p))


###################################
# 5. LAGGED CROSS-CORRELATION (Lags: t, t-1, t-2)
###################################
cat("\n--- Computing Lagged Cross-Correlations (t-2 to t+2) ---\n")

panel_lagged <- panel_df %>%
  arrange(gauge_id, year) %>%
  group_by(gauge_id) %>%
  mutate(
    precip_lag1 = lag(annual_precip, 1),
    precip_lag2 = lag(annual_precip, 2)
  ) %>%
  ungroup() %>%
  filter(!is.na(precip_lag2))

r_lag0 <- cor(panel_lagged$annual_precip, panel_lagged$median_IC, use = "complete.obs")
r_lag1 <- cor(panel_lagged$precip_lag1,   panel_lagged$median_IC, use = "complete.obs")
r_lag2 <- cor(panel_lagged$precip_lag2,   panel_lagged$median_IC, use = "complete.obs")

cat(sprintf("Precipitation Lags Correlation with IC:\n  Lag 0 (Current Year) : r = %+.3f\n  Lag 1 (Prior Year)   : r = %+.3f\n  Lag 2 (Two Years Ago): r = %+.3f\n", 
            r_lag0, r_lag1, r_lag2))


###################################
# 6. VISUALIZATION PIPELINE
###################################
cat("\nGenerating Multi-Variable Figures...\n")

# FIGURE 1: Dynamic Flow State Impact (Dumbbell Slope & Boxplot)
p1_box <- ggplot(panel_df, aes(x = annual_flow_state, y = median_IC, fill = annual_flow_state)) +
  geom_boxplot(alpha = 0.65, outlier.shape = NA, width = 0.45) +
  geom_jitter(width = 0.1, alpha = 0.25, size = 1.2, color = "gray20") +
  scale_fill_manual(values = c("Perennial" = "#2b5c8f", "Intermittent" = "#d95f02")) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.3,
           label = sprintf("Within-Catchment Paired Wilcoxon:\np = %.3e\nMean Delta = %+.3f IC units", 
                           paired_state_test$p.value, mean_state_diff),
           size = 4, fontface = "bold.italic") +
  labs(
    title = "(a) Annual Connectivity by Active Flow State",
    x = "Annual Hydrologic State",
    y = "Catchment Median IC [-]"
  ) +
  theme_bw(base_size = 13) +
  theme(legend.position = "none", panel.grid.minor = element_blank())

# FIGURE 2: Repeated Measures Correlation Plots (Precipitation & Temperature)
p2_precip <- ggplot(panel_df, aes(x = annual_precip, y = median_IC)) +
  geom_point(alpha = 0.2, color = "gray40", size = 1.2) +
  geom_smooth(method = "lm", color = "#2b5c8f", fill = "#2b5c8f", alpha = 0.25, size = 1.2) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.3,
           label = sprintf("rmcorr: r = %+.3f\np = %.2e", rmc_precip$r, rmc_precip$p),
           size = 4, fontface = "bold.italic") +
  labs(title = "(b) IC Sensitivity to Annual Precipitation", x = "Annual Precipitation [mm]", y = "Median IC [-]") +
  theme_bw(base_size = 13) + theme(panel.grid.minor = element_blank())

p2_temp <- ggplot(panel_df, aes(x = annual_temp, y = median_IC)) +
  geom_point(alpha = 0.2, color = "gray40", size = 1.2) +
  geom_smooth(method = "lm", color = "#d95f02", fill = "#d95f02", alpha = 0.25, size = 1.2) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.3,
           label = sprintf("rmcorr: r = %+.3f\np = %.2e", rmc_temp$r, rmc_temp$p),
           size = 4, fontface = "bold.italic") +
  labs(title = "(c) IC Sensitivity to Annual Temperature", x = "Annual Mean Temperature [°C]", y = "Median IC [-]") +
  theme_bw(base_size = 13) + theme(panel.grid.minor = element_blank())

# FIGURE 3: Lagged Cross-Correlation Barplot
lag_df <- data.frame(
  Lag = factor(c("Lag 0 (t)", "Lag 1 (t-1)", "Lag 2 (t-2)"), 
               levels = c("Lag 0 (t)", "Lag 1 (t-1)", "Lag 2 (t-2)")),
  Correlation = c(r_lag0, r_lag1, r_lag2)
)

p3_lag <- ggplot(lag_df, aes(x = Lag, y = Correlation, fill = Correlation > 0)) +
  geom_bar(stat = "identity", width = 0.5, alpha = 0.85) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_fill_manual(values = c("TRUE" = "#2b5c8f", "FALSE" = "#d95f02")) +
  labs(title = "(d) Temporal Memory: Lagged Climate Correlation", x = "Time Lag", y = "Correlation with IC (r)") +
  theme_bw(base_size = 13) + theme(legend.position = "none", panel.grid.minor = element_blank())

fig2_synthesis <- plot_grid(plot_grid(p1_box, p3_lag, ncol = 2), 
                            plot_grid(p2_precip, p2_temp, ncol = 2), 
                            nrow = 2)

ggsave(file.path(out_dir, "Figure_Multivariable_Panel_Synthesis.png"), 
       plot = fig2_synthesis, width = 12, height = 9, dpi = 300)

cat(sprintf("\nScript 2 execution completed. Results saved to: %s\n", out_dir))