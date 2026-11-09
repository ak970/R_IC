# Title: 01_Temporal_Trends_and_Zones.R
# Scope: 31-Year (1991-2021) IC Trajectories, Breakpoints, and Mixed-Effects ANOVA
# Features: Multi-core parallel extraction, Mann-Kendall, Sen's Slope, Pettitt Test, LMM
# It scans your folders across your multi-core CPU and generates Panel_IC_1991_2021_Extracted.csv.
# It performs the regional trajectories, Sen's slope rates of change, Pettitt breakpoint detection, and Type-III Mixed-Effects ANOVA across Geology, Climate, and Regime.
######### If IC files are located in nested subdirectories, replace the following line in section 3: fpath <- file.path(ic_dir, paste0("IC_", yval, "_", gid, ".tif")) 
# with
# fpath <- file.path(ic_dir, as.character(yval), paste0("IC_", yval, "_", gid, ".tif"))

# Packages
required_pkgs <- c("terra", "raster", "parallel", "foreach", "doParallel", 
                   "dplyr", "tidyr", "ggplot2", "cowplot", "trend", "lme4", "lmerTest")
new_pkgs <- required_pkgs[!(required_pkgs %in% installed.packages()[, "Package"])]
if (length(new_pkgs) > 0) install.packages(new_pkgs)

library(terra)
library(raster)
library(parallel)
library(foreach)
library(doParallel)
library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(trend)
library(lme4)
library(lmerTest)

###################################
# 1. SETUP & DIRECTORIES
###################################
# Directory where all IC_year_gaugeID.tif files reside
ic_dir   <- "F:/Peninsular India/Connectivity/IC/Output_Parallel2" 
csv_path <- "F:/Peninsular India/CAMELS_IND/CAMELS_IND_Catchments_Streamflow_Sufficient3/attributes_csv/camels_ind_geol.csv"
out_dir  <- "F:/Peninsular India/Connectivity/IC/Plots/Temporal_1991_2021"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

col_gauge_id <- "gauge_id"
col_regime   <- "Regime"
col_geol     <- "geol_class_1st"
col_climate  <- "Climate"

years <- 1991:2021

# Geology Grouping
geo_grouping <- c(
  'Acid Plutonic Rocks'             = 'Igneous',    
  'Basic Volcanic Rocks'            = 'Igneous',
  'Intermediate Volcanic Rocks'     = 'Igneous',
  'Metamorphic Rocks'               = 'Metamorphic',
  'Mixed Sedimentary Rocks'         = 'Sedimentary',
  'Siliciclastic Sedimentary Rocks' = 'Sedimentary',
  'Unconsolidated Sediments'        = 'Unconsolidated'
)

# Number of parallel cores for data extraction and trend testing
n_cores <- min(32, max(1, parallel::detectCores() - 2))


###################################
# 2. LOAD METADATA
###################################
cat("Loading and preparing metadata...\n")
meta_df <- read.csv(csv_path, stringsAsFactors = FALSE)
meta_df$gauge_id_clean <- sprintf("%05d", as.numeric(trimws(as.character(meta_df[[col_gauge_id]]))))
meta_df$Regime         <- as.factor(trimws(as.character(meta_df[[col_regime]])))
meta_df$Climate        <- as.factor(trimws(as.character(meta_df[[col_climate]])))
meta_df$raw_geol       <- trimws(as.character(meta_df[[col_geol]]))
meta_df$Geology_Class  <- factor(geo_grouping[meta_df$raw_geol], 
                                 levels = c('Igneous', 'Metamorphic', 'Sedimentary', 'Unconsolidated'))

# Build expected panel grid
panel_tasks <- expand.grid(
  gauge_id = unique(meta_df$gauge_id_clean),
  year     = years,
  stringsAsFactors = FALSE
)

cat(sprintf("Scanning for %d potential catchment-year GeoTIFFs across %d cores...\n", 
            nrow(panel_tasks), n_cores))


###################################
# 3. PARALLEL RASTER EXTRACTION
###################################
cl <- parallel::makeCluster(n_cores)
doParallel::registerDoParallel(cl)

clusterExport(cl, varlist = c("ic_dir", "panel_tasks"))
clusterEvalQ(cl, { library(terra) })

# Extract annual summary stats for each basin
panel_extracted <- foreach(
  i = 1:nrow(panel_tasks), 
  .combine = "rbind", 
  .packages = "terra",
  .errorhandling = "pass"
) %dopar% {
  gid  <- panel_tasks$gauge_id[i]
  yval <- panel_tasks$year[i]
  
  fpath <- file.path(ic_dir, paste0("IC_", yval, "_", gid, ".tif"))
  
  if (file.exists(fpath)) {
    r <- terra::rast(fpath)
    vals <- terra::values(r, mat = FALSE)
    clean <- vals[!is.na(vals) & !is.nan(vals) & !is.infinite(vals)]
    
    if (length(clean) > 50) {
      return(data.frame(
        gauge_id    = gid,
        year        = yval,
        median_IC   = median(clean),
        mean_IC     = mean(clean),
        q25_IC      = quantile(clean, 0.25),
        q75_IC      = quantile(clean, 0.75),
        valid_cells = length(clean),
        stringsAsFactors = FALSE
      ))
    }
  }
  return(NULL)
}

parallel::stopCluster(cl)

# Merge with static attributes
panel_df <- panel_extracted %>%
  left_join(meta_df %>% dplyr::select(gauge_id_clean, Regime, Climate, Geology_Class), 
            by = c("gauge_id" = "gauge_id_clean")) %>%
  filter(!is.na(Geology_Class) & !is.na(Climate) & !is.na(Regime))

write.csv(panel_df, file.path(out_dir, "Panel_IC_1991_2021_Extracted.csv"), row.names = FALSE)
cat(sprintf("Extracted %d valid catchment-year observations across 31 years.\n", nrow(panel_df)))


###################################
# 4. MANN-KENDALL & PETTITT BREAKPOINT TESTS
###################################
cat("\nRunning parallel Mann-Kendall and Pettitt tests for each catchment...\n")

unique_basins <- unique(panel_df$gauge_id)

cl_trend <- parallel::makeCluster(n_cores)
doParallel::registerDoParallel(cl_trend)
clusterExport(cl_trend, varlist = c("panel_df", "unique_basins"))
clusterEvalQ(cl_trend, { library(trend); library(dplyr) })

trend_results <- foreach(
  b = seq_along(unique_basins), 
  .combine = "rbind", 
  .packages = "trend"
) %dopar% {
  gid <- unique_basins[b]
  sub_data <- panel_df %>% filter(gauge_id == gid) %>% arrange(year)
  
  # Require at least 20 years of valid observations
  if (nrow(sub_data) >= 20) {
    ts_ic <- ts(sub_data$median_IC, start = min(sub_data$year))
    
    # 1. Mann-Kendall Trend Test
    mk  <- trend::mk.test(ts_ic)
    # 2. Sen's Slope
    sen <- trend::sens.slope(ts_ic)
    # 3. Pettitt Breakpoint Test
    pet <- trend::pettitt.test(ts_ic)
    
    break_year <- sub_data$year[pet$estimate]
    
    return(data.frame(
      gauge_id     = gid,
      MK_tau       = mk$estimates["tau"],
      MK_pvalue    = mk$p.value,
      Sens_slope   = sen$estimates,
      Pettitt_year = break_year,
      Pettitt_pval = pet$p.value,
      stringsAsFactors = FALSE
    ))
  }
  return(NULL)
}

parallel::stopCluster(cl_trend)

trend_df <- trend_results %>%
  left_join(meta_df %>% dplyr::select(gauge_id_clean, Regime, Climate, Geology_Class), 
            by = c("gauge_id" = "gauge_id_clean"))

write.csv(trend_df, file.path(out_dir, "Catchment_Trend_and_Breakpoints_Summary.csv"), row.names = FALSE)


###################################
# 5. LINEAR MIXED-EFFECTS ANOVA (LMM)
###################################
cat("\n--- Running Linear Mixed-Effects Model (LMM / Repeated Measures ANOVA) ---\n")

# Center year for numerical stability
panel_df$year_c <- panel_df$year - 1991

# Model: Main effects of Geology, Climate, Regime, Time + Geology x Time interaction
# Random intercept per catchment: (1 | gauge_id)
lmm_model <- lmer(median_IC ~ year_c * Geology_Class + Climate + Regime + (1 | gauge_id), 
                  data = panel_df)

cat("Type-III Analysis of Variance (Satterthwaite's method):\n")
lmm_anova <- anova(lmm_model, type = 3)
print(lmm_anova)
write.csv(as.data.frame(lmm_anova), file.path(out_dir, "LMM_TypeIII_ANOVA_Table.csv"))


###################################
# 6. VISUALIZATION PIPELINE
###################################
cat("\nGenerating publication figures...\n")

# FIGURE 1: 31-Year Regional Trajectory with 95% Confidence Interval
regional_annual <- panel_df %>%
  group_by(year) %>%
  summarize(
    mean_IC  = mean(median_IC),
    se_IC    = sd(median_IC) / sqrt(n()),
    ci_lower = mean_IC - 1.96 * se_IC,
    ci_upper = mean_IC + 1.96 * se_IC,
    .groups  = "drop"
  )

p1 <- ggplot(regional_annual, aes(x = year, y = mean_IC)) +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), fill = "#2b5c8f", alpha = 0.25) +
  geom_line(color = "#1a365d", size = 1.3) +
  geom_point(color = "#1a365d", size = 2) +
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed", color = "red", size = 0.8) +
  labs(
    title = "(a) Regional 31-Year Sediment Connectivity Trajectory (1991–2021)",
    subtitle = "Solid curve: Mean regional IC | Shaded: 95% CI | Red dashed: Linear trajectory",
    x = "Year",
    y = "Catchment Median IC [-]"
  ) +
  theme_bw(base_size = 14) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"))

# FIGURE 2: Temporal Evolution by Geological Class (4 Ribbons)
geol_annual <- panel_df %>%
  group_by(Geology_Class, year) %>%
  summarize(
    mean_IC  = mean(median_IC),
    se_IC    = sd(median_IC) / sqrt(n()),
    ci_lower = mean_IC - 1.96 * se_IC,
    ci_upper = mean_IC + 1.96 * se_IC,
    .groups  = "drop"
  )

geol_palette <- c("Igneous" = "#e41a1c", "Metamorphic" = "#984ea3", 
                  "Sedimentary" = "#4daf4a", "Unconsolidated" = "#ff7f00")

p2 <- ggplot(geol_annual, aes(x = year, y = mean_IC, color = Geology_Class, fill = Geology_Class)) +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.18, color = NA) +
  geom_line(size = 1.2) +
  scale_color_manual(values = geol_palette, name = "Geology") +
  scale_fill_manual(values = geol_palette, name = "Geology") +
  labs(
    title = "(b) Long-Term Connectivity Evolution Stratified by Lithology",
    subtitle = "Solid curves: Mean annual IC per lithological group | Shaded bands: 95% CIs",
    x = "Year",
    y = "Catchment Median IC [-]"
  ) +
  theme_bw(base_size = 14) +
  theme(legend.position = "right", panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"))

# FIGURE 3: Distribution of Sen's Slopes & Pettitt Breakpoint Years
p3a <- ggplot(trend_df, aes(x = Geology_Class, y = Sens_slope * 10, fill = Geology_Class)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray30") +
  geom_boxplot(alpha = 0.7, outlier.shape = NA, width = 0.5) +
  geom_jitter(width = 0.15, alpha = 0.6, size = 1.8) +
  scale_fill_manual(values = geol_palette) +
  labs(
    title = "(c) Rate of Change (Decadal Sen's Slope)",
    x = "Geology Class",
    y = expression("Sen's Slope [ " * Delta * "IC / decade ]")
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none", panel.grid.minor = element_blank())

p3b <- ggplot(trend_df %>% filter(Pettitt_pval < 0.05), aes(x = Pettitt_year)) +
  geom_histogram(binwidth = 2, fill = "#2b5c8f", color = "white", alpha = 0.8) +
  labs(
    title = "(d) Statistically Significant Tipping Years (Pettitt Test p < 0.05)",
    x = "Detected Breakpoint Year",
    y = "Number of Catchments"
  ) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank())

fig_synthesis <- plot_grid(p1, p2, plot_grid(p3a, p3b, ncol = 2), ncol = 1, rel_heights = c(1, 1, 0.9))
ggsave(file.path(out_dir, "Figure_31Year_Temporal_Geology_Synthesis.png"), 
       plot = fig_synthesis, width = 12, height = 15, dpi = 300)

cat(sprintf("\nScript 1 execution completed. Results saved to: %s\n", out_dir))