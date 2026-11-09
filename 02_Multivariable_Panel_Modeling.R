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

# ###################################
# # 1. SETUP & DATA INGESTION
# ###################################
# out_dir <- "F:/Peninsular India/Connectivity/IC/Plots/Multivariable_Panel"
# if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
# 
# # 1. Load the 31-year extracted IC panel data (created by Script 1)
# panel_ic_path <- "F:/Peninsular India/Connectivity/IC/Plots/Temporal_1991_2021/Panel_IC_1991_2021_Extracted.csv"
# # ==============================================================================
# # UPDATED DATA INGESTION & COVARIATE SPECIFICATION (SCRIPT 2)
# # ==============================================================================
# 
# # ------------------------------------------------------------------------------
# # 1. MERGE ANNUAL DYNAMIC COVARIATES
# # ------------------------------------------------------------------------------
# # Assuming your dynamic CSV has columns: gauge_id, year, flow_state, precip, temp, pet
# # If reading an external annual CSV:
# # dynamic_df <- read.csv("F:/.../annual_climate_and_flow.csv", stringsAsFactors = FALSE)
# # dynamic_df$gauge_id <- sprintf("%05d", as.numeric(dynamic_df$gauge_id))
# # panel_df <- panel_df %>% left_join(dynamic_df, by = c("gauge_id", "year"))
# 
# # Rename/assign your specific column names:
# # (Replace the column assignments below with your actual merged column names)
# panel_df <- panel_df %>%
#   mutate(
#     # Binary state: 0 = Perennial, 1 = Intermittent (can contain NA)
#     annual_flow_state = factor(Regime, levels = c(0, 1), 
#                                labels = c("Perennial", "Intermittent")),
#     
#     # Continuous annual climate variables
#     annual_precip = as.numeric(annual_precip), # mm
#     annual_temp   = as.numeric(annual_temp),   # °C
#     annual_pet    = as.numeric(annual_pet),    # mm
#     
#     # Dynamic annual Aridity Index: P / PET
#     annual_ai     = annual_precip / annual_pet
#   )




# ==============================================================================
# SECTION 1: SETUP & MERGING SEPARATE ANNUAL DATASET (SCRIPT 2)
# ==============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(plm)
library(rmcorr)

out_dir <- "F:/Peninsular India/Connectivity/IC/Plots/Multivariable_Panel"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ------------------------------------------------------------------------------
# 1. FILE PATHS SETUP
# ------------------------------------------------------------------------------
# Output from Script 1 (contains gauge_id, year, median_IC, and static geology/climate)
panel_ic_path   <- "F:/Peninsular India/Connectivity/IC/Plots/Temporal_1991_2021/Panel_IC_1991_2021_Extracted.csv"

# YOUR SEPARATE ANNUAL DATASET CSV (contains year-by-year climate and flow records)
annual_csv_path <- "D:/GitHub/IntermittencyTrends/results/Mine/Master_Dataset_Updated2.csv"

# Column names inside your annual CSV:
col_ann_gid    <- "gauge_id"           # Column for gauge ID in annual CSV
col_ann_year   <- "hydro_year"         # Column for Year in annual CSV
col_ann_flow   <- "Flow_State"         # Column for binary flow state (0 = Perennial, 1 = Intermittent)
col_ann_precip <- "p"                  # Column for annual rainfall (mm)
col_ann_temp   <- "tmax"        # Column for annual temperature (°C)
col_ann_pet    <- "pet"                # Column for annual PET (mm)
col_soil_mois  <- "sm_lvl1(kg/m2)"     # Column for soil moisture
col_wind       <- "wind_u(m/s)"        # column for wind speed
col_srad       <- "srad_sw(w/m2)"      # column for short-wave radiation

# ------------------------------------------------------------------------------
# 2. LOAD & STANDARDIZE THE SEPARATE ANNUAL CSV
# ------------------------------------------------------------------------------
cat("Loading separate annual time-series dataset...\n")
annual_raw <- read.csv(annual_csv_path, stringsAsFactors = FALSE, check.names = FALSE)

# Format gauge_id with leading zeroes to match the 5-digit IC dataset
annual_raw$gauge_id_clean <- sprintf("%05d", as.numeric(trimws(as.character(annual_raw[[col_ann_gid]]))))
annual_raw$year_clean     <- as.integer(trimws(as.character(annual_raw[[col_ann_year]])))

# Select and rename your relevant annual variables
annual_clean <- annual_raw %>%
  dplyr::select(
    gauge_id          = gauge_id_clean,
    year              = year_clean,
    annual_flow_state = all_of(col_ann_flow),
    annual_precip     = all_of(col_ann_precip),
    annual_temp       = all_of(col_ann_temp),
    annual_pet        = all_of(col_ann_pet),
    annual_soil_mois  = all_of(col_soil_mois),
    annual_wind       = all_of(col_wind),
    annual_srad       = all_of(col_srad)
  )

# ==============================================================================
# SECTIONS 3–6: MULTI-VARIABLE MODELING, LAGGED CORRELATION & MASTER PLOTTING
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. DEFINE CONTINUOUS VARIABLES AND DISPLAY LABELS
# ------------------------------------------------------------------------------
# Add as many variables as exist in panel_df!
var_definitions <- c(
  "annual_precip"        = "Annual Precipitation [mm]",
  "annual_temp"          = "Annual Mean Temperature [°C]",
  "annual_pet"           = "Annual Potential Evapotranspiration [mm]",
  "annual_soil_moisture" = "Soil Moisture Layer 1 [kg/m²]",
  "annual_ai"            = "Aridity Index (P/PET) [-]"
)

# Filter to variables that exist in your data frame
active_vars <- var_definitions[names(var_definitions) %in% names(panel_df)]
cat(sprintf("Active variables identified for modeling (%d):\n", length(active_vars)))
print(names(active_vars))


# ------------------------------------------------------------------------------
# 2. DYNAMIC PANEL FIXED-EFFECTS ("WITHIN") REGRESSION
# ------------------------------------------------------------------------------
cat("\n--- Running Dynamic Panel Fixed-Effects Regression ---\n")

complete_vars <- c("gauge_id", "year", "median_IC", "annual_flow_state", names(active_vars))
panel_complete <- panel_df %>%
  dplyr::select(all_of(complete_vars)) %>%
  filter(complete.cases(.))

pdata <- pdata.frame(panel_complete, index = c("gauge_id", "year"))

fe_formula <- as.formula(paste("median_IC ~ annual_flow_state +", paste(names(active_vars), collapse = " + ")))
cat("Model Formula:", deparse(fe_formula), "\n\n")

fe_model <- plm(fe_formula, data = pdata, model = "within")
fe_summary <- summary(fe_model)
print(fe_summary)

write.csv(as.data.frame(fe_summary$coefficients), 
          file.path(out_dir, "Panel_FixedEffects_All_Variables.csv"))


# ------------------------------------------------------------------------------
# 3. REPEATED-MEASURES CORRELATION & SCATTER PLOTS
# ------------------------------------------------------------------------------
cat("\n--- Generating Automated Plots and rmcorr for All Variables ---\n")

plot_list <- list()
rmc_summary <- data.frame(
  Variable = character(),
  Label    = character(),
  r_rmcorr = numeric(),
  p_value  = numeric(),
  stringsAsFactors = FALSE
)

plot_colors <- c("#2b5c8f", "#d95f02", "#1b9e77", "#7570b3", "#e7298a", "#e6ab02", "#a6761d")

for (idx in seq_along(active_vars)) {
  var_col   <- names(active_vars)[idx]
  var_label <- active_vars[idx]
  cur_color <- plot_colors[((idx - 1) %% length(plot_colors)) + 1]
  
  # Within-catchment correlation using rmcorr
  sub_df <- panel_df %>% 
    dplyr::select(gauge_id, median_IC, x_val = all_of(var_col)) %>%
    filter(!is.na(x_val) & !is.na(median_IC))
  
  rmc_res <- rmcorr::rmcorr(participant = gauge_id, measure1 = x_val, measure2 = median_IC, dataset = sub_df)
  
  rmc_summary <- rbind(rmc_summary, data.frame(
    Variable = var_col,
    Label    = var_label,
    r_rmcorr = rmc_res$r,
    p_value  = rmc_res$p,
    stringsAsFactors = FALSE
  ))
  
  # Panel tag starts from (c) since (a) is Flow State and (b) is Lagged Correlation
  tag_letter <- letters[idx + 2]
  
  p_cur <- ggplot(sub_df, aes(x = x_val, y = median_IC)) +
    geom_point(alpha = 0.15, color = "gray30", size = 1.1) +
    geom_smooth(method = "lm", color = cur_color, fill = cur_color, alpha = 0.22, size = 1.2) +
    annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.3,
             label = sprintf("rmcorr: r = %+.3f\np = %.2e", rmc_res$r, rmc_res$p),
             size = 3.8, fontface = "bold.italic", color = "gray20") +
    labs(
      title = paste0("(", tag_letter, ") Sensitivity to ", sub(" \\[.*", "", var_label)),
      x     = var_label,
      y     = "Median IC [-]"
    ) +
    theme_bw(base_size = 12) +
    theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"))
  
  plot_list[[var_col]] <- p_cur
}

write.csv(rmc_summary, file.path(out_dir, "rmcorr_All_Variables_Summary.csv"), row.names = FALSE)


# ------------------------------------------------------------------------------
# 4. ANNUAL FLOW STATE PLOT — Panel (a)
# ------------------------------------------------------------------------------
p1_box <- ggplot(panel_df %>% filter(!is.na(annual_flow_state)), 
                 aes(x = annual_flow_state, y = median_IC, fill = annual_flow_state)) +
  geom_boxplot(alpha = 0.65, outlier.shape = NA, width = 0.45) +
  geom_jitter(width = 0.08, alpha = 0.2, size = 1.1, color = "gray20") +
  scale_fill_manual(values = c("Perennial" = "#2b5c8f", "Intermittent" = "#d95f02")) +
  labs(
    title = "(a) IC by Hydrologic State",
    x     = "Annual Flow Regime",
    y     = "Median IC [-]"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none", panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"))


# ------------------------------------------------------------------------------
# 5. LAGGED CLIMATE CROSS-CORRELATION — Panel (b)
# ------------------------------------------------------------------------------
cat("\n--- Computing Lagged Climate Cross-Correlations (t, t-1, t-2) ---\n")

panel_lagged <- panel_df %>%
  arrange(gauge_id, year) %>%
  group_by(gauge_id) %>%
  mutate(
    precip_lag1 = lag(annual_precip, 1),
    precip_lag2 = lag(annual_precip, 2)
  ) %>%
  ungroup() %>%
  filter(!is.na(precip_lag2) & !is.na(median_IC))

# Compute correlation for Lag 0, Lag 1, and Lag 2
r_lag0 <- cor(panel_lagged$annual_precip, panel_lagged$median_IC, use = "complete.obs")
r_lag1 <- cor(panel_lagged$precip_lag1,   panel_lagged$median_IC, use = "complete.obs")
r_lag2 <- cor(panel_lagged$precip_lag2,   panel_lagged$median_IC, use = "complete.obs")

cat(sprintf("Precipitation Lag Correlations:\n  Lag 0 (t)   : r = %+.4f\n  Lag 1 (t-1) : r = %+.4f\n  Lag 2 (t-2) : r = %+.4f\n", 
            r_lag0, r_lag1, r_lag2))

lag_df <- data.frame(
  Lag = factor(c("Lag 0 (t)", "Lag 1 (t-1)", "Lag 2 (t-2)"), 
               levels = c("Lag 0 (t)", "Lag 1 (t-1)", "Lag 2 (t-2)")),
  Correlation = c(r_lag0, r_lag1, r_lag2)
)

# Barplot colored dynamically (Blue if positive, Orange if negative)
p_lag <- ggplot(lag_df, aes(x = Lag, y = Correlation, fill = Correlation > 0)) +
  geom_bar(stat = "identity", width = 0.5, alpha = 0.85) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_fill_manual(values = c("TRUE" = "#2b5c8f", "FALSE" = "#d95f02")) +
  labs(
    title = "(b) Temporal Memory: Lagged Climate Correlation",
    x     = "Time Lag",
    y     = "Correlation with IC (r)"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none", panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"))

# Save standalone lagged plot
ggsave(file.path(out_dir, "Figure_Lagged_Climate_Correlation.png"), 
       plot = p_lag, width = 6.5, height = 5, dpi = 300)


# ------------------------------------------------------------------------------
# 6. ASSEMBLE MASTER MULTI-PANEL SYNTHESIS FIGURE
# ------------------------------------------------------------------------------
cat("\nCombining all panels into master synthesis figure...\n")

# Top row: (a) Flow State + (b) Lagged Climate Correlation
# Subsequent rows: (c), (d)... continuous scatter plots
all_panels <- c(list(p1_box, p_lag), plot_list)

n_plots  <- length(all_panels)
num_rows <- ceiling(n_plots / 2)

fig_master <- cowplot::plot_grid(plotlist = all_panels, ncol = 2, align = "hv")
fig_height <- max(8, num_rows * 3.8)

out_fig_master <- file.path(out_dir, "Figure_Multivariable_Synthesis_All_Vars.png")
ggsave(out_fig_master, plot = fig_master, width = 12, height = fig_height, dpi = 300)

cat(sprintf("\nMaster synthesis figure saved to:\n  -> %s\n", out_fig_master))