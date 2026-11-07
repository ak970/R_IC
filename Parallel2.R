# ==============================================================================
# SEQUENTIAL RESCUE SCRIPT FOR THE 15 FAILED CATCHMENTS (MULTI-YEAR & CLEAN NAMING)
# ==============================================================================

library(terra)
library(raster)
library(compiler)

# 1. Resolve TauDEM PROJ environment conflict
if (dir.exists(system.file("proj", package = "terra"))) {
  Sys.setenv(PROJ_LIB = system.file("proj", package = "terra"))
} else {
  Sys.unsetenv("PROJ_LIB")
}
rasterOptions(tolerance = 0.05)

# Directories
wd           <- "F:/Peninsular India/Connectivity/IC/Output_Parallel"
dem_dir      <- "F:/Peninsular India/Connectivity/IC/Input/DEM"
# Parent directory containing multiple subdirectories with NDVI files
w_factor_dir <- "F:/Peninsular India/Connectivity/IC/Input/NDVI"

if (!dir.exists(wd)) dir.create(wd, recursive = TRUE)
setwd(wd)

# Target coordinate system & resolution
target_crs    <- "+proj=aea +lat_1=7 +lat_2=32 +lat_0=20 +lon_0=80 +x_0=2000000 +y_0=2000000 +datum=WGS84 +units=m +no_defs"
target_res    <- 100
mw            <- 5
mpi_cores     <- 4  # Safe to use 4 cores for TauDEM since we process 1 basin at a time

# The 15 failed gauge IDs (formatted to 5 digits)
failed_ids_raw <- c(3043, 3073, 3075, 4060, 4062, 8038, 12004, 12020, 12025, 12035, 12046, 13007, 13016, 13017, 14008)
failed_ids     <- sprintf("%05d", failed_ids_raw)

# ------------------------------------------------------------------------------
# Scan parent directory recursively and match NDVI files for the failed gauges
# ------------------------------------------------------------------------------
cat("Scanning parent directory and subdirectories recursively for NDVI files...\n")
all_ndvi <- list.files(
  w_factor_dir, 
  pattern = "NDVI_[0-9]{4}_[0-9]+\\.(tif|TIF)$", 
  recursive = TRUE, 
  full.names = TRUE
)

rescue_tasks <- data.frame(
  gauge_id      = character(),
  year          = character(),
  dem_path      = character(),
  w_factor_path = character(),
  stringsAsFactors = FALSE
)

for (nf in all_ndvi) {
  fname    <- basename(nf)
  year_val <- sub(".*NDVI_([0-9]{4})_([0-9]+)\\..*", "\\1", fname, ignore.case = TRUE)
  raw_gid  <- sub(".*NDVI_([0-9]{4})_([0-9]+)\\..*", "\\2", fname, ignore.case = TRUE)
  gid      <- sprintf("%05d", as.numeric(raw_gid))
  
  if (gid %in% failed_ids) {
    dem_p <- file.path(dem_dir, paste0("DEM_", gid, ".tif"))
    if (file.exists(dem_p)) {
      rescue_tasks <- rbind(rescue_tasks, data.frame(
        gauge_id      = gid,
        year          = year_val,
        dem_path      = dem_p,
        w_factor_path = nf,
        stringsAsFactors = FALSE
      ))
    }
  }
}

rescue_tasks <- rescue_tasks[order(rescue_tasks$year, rescue_tasks$gauge_id), ]
cat(sprintf("Found %d rescue task(s) for the %d failed gauge IDs.\n", nrow(rescue_tasks), length(failed_ids)))


# ------------------------------------------------------------------------------
# Byte-compiled flow-length engine
# ------------------------------------------------------------------------------
flow_length_engine <- function(dir8_m, W_FL_m, Wgt_m, target_idx, res) {
  posrow <- target_idx[, "row"] + 1
  poscol <- target_idx[, "col"] + 1
  
  r1 <- posrow;     c1 <- poscol - 1; idx1 <- which(dir8_m[cbind(r1, c1)] == 1); pos1row_0 <- r1[idx1]; pos1col_0 <- c1[idx1]; W_FL_m[cbind(pos1row_0, pos1col_0)] <- 0
  r2 <- posrow + 1; c2 <- poscol - 1; idx2 <- which(dir8_m[cbind(r2, c2)] == 2); pos2row_0 <- r2[idx2]; pos2col_0 <- c2[idx2]; W_FL_m[cbind(pos2row_0, pos2col_0)] <- 0
  r3 <- posrow + 1; c3 <- poscol;     idx3 <- which(dir8_m[cbind(r3, c3)] == 3); pos3row_0 <- r3[idx3]; pos3col_0 <- c3[idx3]; W_FL_m[cbind(pos3row_0, pos3col_0)] <- 0
  r4 <- posrow + 1; c4 <- poscol + 1; idx4 <- which(dir8_m[cbind(r4, c4)] == 4); pos4row_0 <- r4[idx4]; pos4col_0 <- c4[idx4]; W_FL_m[cbind(pos4row_0, pos4col_0)] <- 0
  r5 <- posrow;     c5 <- poscol + 1; idx5 <- which(dir8_m[cbind(r5, c5)] == 5); pos5row_0 <- r5[idx5]; pos5col_0 <- c5[idx5]; W_FL_m[cbind(pos5row_0, pos5col_0)] <- 0
  r6 <- posrow - 1; c6 <- poscol + 1; idx6 <- which(dir8_m[cbind(r6, c6)] == 6); pos6row_0 <- r6[idx6]; pos6col_0 <- c6[idx6]; W_FL_m[cbind(pos6row_0, pos6col_0)] <- 0
  r7 <- posrow - 1; c7 <- poscol;     idx7 <- which(dir8_m[cbind(r7, c7)] == 7); pos7row_0 <- r7[idx7]; pos7col_0 <- c7[idx7]; W_FL_m[cbind(pos7row_0, pos7col_0)] <- 0
  r8 <- posrow - 1; c8 <- poscol - 1; idx8 <- which(dir8_m[cbind(r8, c8)] == 8); pos8row_0 <- r8[idx8]; pos8col_0 <- c8[idx8]; W_FL_m[cbind(pos8row_0, pos8col_0)] <- 0
  
  count <- 1
  while ((length(pos1col_0) | length(pos2col_0) | length(pos3col_0) | length(pos4col_0) | 
          length(pos5col_0) | length(pos6col_0) | length(pos7col_0) | length(pos8col_0)) > 0) {
    
    pos1row_0_bis <- pos1row_0; pos1col_0_bis <- pos1col_0 + 1
    pos2row_0_bis <- pos2row_0 - 1; pos2col_0_bis <- pos2col_0 + 1
    pos3row_0_bis <- pos3row_0 - 1; pos3col_0_bis <- pos3col_0
    pos4row_0_bis <- pos4row_0 - 1; pos4col_0_bis <- pos4col_0 - 1
    pos5row_0_bis <- pos5row_0;     pos5col_0_bis <- pos5col_0 - 1
    pos6row_0_bis <- pos6row_0 + 1; pos6col_0_bis <- pos6col_0 - 1
    pos7row_0_bis <- pos7row_0 + 1; pos7col_0_bis <- pos7col_0
    pos8row_0_bis <- pos8row_0 + 1; pos8col_0_bis <- pos8col_0 + 1
    
    if (count == 1) {
      W_FL_m[cbind(pos1row_0, pos1col_0)] <- 0; W_FL_m[cbind(pos2row_0, pos2col_0)] <- 0
      W_FL_m[cbind(pos3row_0, pos3col_0)] <- 0; W_FL_m[cbind(pos4row_0, pos4col_0)] <- 0
      W_FL_m[cbind(pos5row_0, pos5col_0)] <- 0; W_FL_m[cbind(pos6row_0, pos6col_0)] <- 0
      W_FL_m[cbind(pos7row_0, pos7col_0)] <- 0; W_FL_m[cbind(pos8row_0, pos8col_0)] <- 0
    } else {
      W_FL_m[cbind(pos1row_0, pos1col_0)] <- W_FL_m[cbind(pos1row_0_bis, pos1col_0_bis)] + (res * ((Wgt_m[cbind(pos1row_0_bis, pos1col_0_bis)] + Wgt_m[cbind(pos1row_0, pos1col_0)]) / 2))
      W_FL_m[cbind(pos2row_0, pos2col_0)] <- W_FL_m[cbind(pos2row_0_bis, pos2col_0_bis)] + (res * sqrt(2) * ((Wgt_m[cbind(pos2row_0_bis, pos2col_0_bis)] + Wgt_m[cbind(pos2row_0, pos2col_0)]) / 2))
      W_FL_m[cbind(pos3row_0, pos3col_0)] <- W_FL_m[cbind(pos3row_0_bis, pos3col_0_bis)] + (res * ((Wgt_m[cbind(pos3row_0_bis, pos3col_0_bis)] + Wgt_m[cbind(pos3row_0, pos3col_0)]) / 2))
      W_FL_m[cbind(pos4row_0, pos4col_0)] <- W_FL_m[cbind(pos4row_0_bis, pos4col_0_bis)] + (res * sqrt(2) * ((Wgt_m[cbind(pos4row_0_bis, pos4col_0_bis)] + Wgt_m[cbind(pos4row_0, pos4col_0)]) / 2))
      W_FL_m[cbind(pos5row_0, pos5col_0)] <- W_FL_m[cbind(pos5row_0_bis, pos5col_0_bis)] + (res * ((Wgt_m[cbind(pos5row_0_bis, pos5col_0_bis)] + Wgt_m[cbind(pos5row_0, pos5col_0)]) / 2))
      W_FL_m[cbind(pos6row_0, pos6col_0)] <- W_FL_m[cbind(pos6row_0_bis, pos6col_0_bis)] + (res * sqrt(2) * ((Wgt_m[cbind(pos6row_0_bis, pos6col_0_bis)] + Wgt_m[cbind(pos6row_0, pos6col_0)]) / 2))
      W_FL_m[cbind(pos7row_0, pos7col_0)] <- W_FL_m[cbind(pos7row_0_bis, pos7col_0_bis)] + (res * ((Wgt_m[cbind(pos7row_0_bis, pos7col_0_bis)] + Wgt_m[cbind(pos7row_0, pos7col_0)]) / 2))
      W_FL_m[cbind(pos8row_0, pos8col_0)] <- W_FL_m[cbind(pos8row_0_bis, pos8col_0_bis)] + (res * sqrt(2) * ((Wgt_m[cbind(pos8row_0_bis, pos8col_0_bis)] + Wgt_m[cbind(pos8row_0, pos8col_0)]) / 2))
    }
    
    pos0row_bis <- c(pos1row_0, pos2row_0, pos3row_0, pos4row_0, pos5row_0, pos6row_0, pos7row_0, pos8row_0)
    pos0col_bis <- c(pos1col_0, pos2col_0, pos3col_0, pos4col_0, pos5col_0, pos6col_0, pos7col_0, pos8col_0)
    
    r1 <- pos0row_bis;     c1 <- pos0col_bis - 1; idx1 <- which(dir8_m[cbind(r1, c1)] == 1); pos1row_0 <- r1[idx1]; pos1col_0 <- c1[idx1]
    r2 <- pos0row_bis + 1; c2 <- pos0col_bis - 1; idx2 <- which(dir8_m[cbind(r2, c2)] == 2); pos2row_0 <- r2[idx2]; pos2col_0 <- c2[idx2]
    r3 <- pos0row_bis + 1; c3 <- pos0col_bis;     idx3 <- which(dir8_m[cbind(r3, c3)] == 3); pos3row_0 <- r3[idx3]; pos3col_0 <- c3[idx3]
    r4 <- pos0row_bis + 1; c4 <- pos0col_bis + 1; idx4 <- which(dir8_m[cbind(r4, c4)] == 4); pos4row_0 <- r4[idx4]; pos4col_0 <- c4[idx4]
    r5 <- pos0row_bis;     c5 <- pos0col_bis + 1; idx5 <- which(dir8_m[cbind(r5, c5)] == 5); pos5row_0 <- r5[idx5]; pos5col_0 <- c5[idx5]
    r6 <- pos0row_bis - 1; c6 <- pos0col_bis + 1; idx6 <- which(dir8_m[cbind(r6, c6)] == 6); pos6row_0 <- r6[idx6]; pos6col_0 <- c6[idx6]
    r7 <- pos0row_bis - 1; c7 <- pos0col_bis;     idx7 <- which(dir8_m[cbind(r7, c7)] == 7); pos7row_0 <- r7[idx7]; pos7col_0 <- c7[idx7]
    r8 <- pos0row_bis - 1; c8 <- pos0col_bis - 1; idx8 <- which(dir8_m[cbind(r8, c8)] == 8); pos8row_0 <- r8[idx8]; pos8col_0 <- c8[idx8]
    
    count <- count + 1
  }
  return(list(W_FL_m = W_FL_m, count = count))
}

flow_length_engine_cmp <- compiler::cmpfun(flow_length_engine)

# ------------------------------------------------------------------------------
# Process each failed catchment sequentially with dedicated RAM
# ------------------------------------------------------------------------------
rescue_results <- list()

for (k in 1:nrow(rescue_tasks)) {
  gid      <- rescue_tasks$gauge_id[k]
  year_str <- rescue_tasks$year[k]
  dem_p    <- rescue_tasks$dem_path[k]
  wf_p     <- rescue_tasks$w_factor_path[k]
  
  cat(sprintf("\n[%d/%d] Processing Gauge %s | Year %s...\n", k, nrow(rescue_tasks), gid, year_str))
  
  # Sandbox folder per catchment and year
  temp_dir <- file.path(wd, paste0("temp_rescue_", gid, "_", year_str))
  if (!dir.exists(temp_dir)) dir.create(temp_dir, recursive = TRUE)
  setwd(temp_dir)
  
  t0 <- Sys.time()
  
  tryCatch({
    # 1. Reproject DEM & PitRemove
    dem_raw    <- terra::rast(dem_p)
    dem_proj_t <- terra::project(dem_raw, target_crs, res = target_res, method = "bilinear")
    writeRaster(raster(dem_proj_t), "dem_raw_proj.tif", overwrite = TRUE)
    
    system(paste("mpiexec -n", mpi_cores, "PitRemove -z dem_raw_proj.tif -fel dem_pitfilled.tif"))
    z_fill  <- raster("dem_pitfilled.tif")
    crs_out <- crs(z_fill)
    
    # 2. Align NDVI & compute C-factor
    ndvi_raw  <- terra::rast(wf_p)
    wf_proj_t <- terra::project(ndvi_raw, terra::rast("dem_pitfilled.tif"), method = "bilinear")
    W_factor  <- raster(wf_proj_t)
    
    vals <- getValues(W_factor)
    vals_clean <- vals[!is.na(vals) & !is.nan(vals)]
    ndvi_min <- min(vals_clean); ndvi_max <- max(vals_clean)
    
    f_g <- if (ndvi_max == ndvi_min) W_factor * 0 else ((W_factor - ndvi_min) / (ndvi_max - ndvi_min)) * 100
    f_clamped <- clamp(f_g, lower = 0.095, upper = 78.3)
    C_factor  <- 0.6508 - (0.3436 * log10(f_clamped))
    C_factor[f_g <= 0]     <- 1.0
    C_factor[f_g > 78.3]   <- 0.001
    C_factor[C_factor > 1] <- 1.0
    C_factor[C_factor < 0.001] <- 0.001
    writeRaster(C_factor, "W_factor.tif", overwrite = TRUE)
    W_factor <- C_factor
    
    # 3. D8 Flow directions & Accumulation
    system(paste("mpiexec -n", mpi_cores, "D8Flowdir -p D8_flow_dir.tif -sd8 slope_D8.tif -fel dem_pitfilled.tif"))
    flow_dir_D8 <- raster("D8_flow_dir.tif")
    slope_D8    <- raster("slope_D8.tif")
    
    system(paste("mpiexec -n", mpi_cores, "AreaD8 -p D8_flow_dir.tif -ad8 accD8.tif -nc"))
    accD8 <- raster("accD8.tif")
    
    # 4. Target Catchment Outlet (Max Accumulation)
    max_cell   <- which.max(getValues(accD8))
    snapped_xy <- xyFromCell(accD8, max_cell)
    max_acc    <- max(getValues(accD8), na.rm = TRUE)
    
    tgt_snapped      <- terra::vect(matrix(snapped_xy, ncol = 2), crs = target_crs)
    tgt_snapped$burn <- 1
    target_rst       <- raster(terra::rasterize(tgt_snapped, terra::rast(z_fill), field = "burn", background = 2))
    target_rst       <- reclassify(target_rst, c(0.9, 1.1, NA, 1.9, 2.1, 1))
    
    flow_dir_D8[is.na(flow_dir_D8)] <- -9999
    flow_dir_D8 <- flow_dir_D8 * target_rst
    
    slope_D8_r_m <- reclassify(slope_D8, c(-2, 0.005, 0.005, 1, 10^50, 1))
    writeRaster(slope_D8_r_m, "slope_D8_rec_m.tif", overwrite = TRUE)
    
    inv_CS <- 1 / (W_factor * slope_D8_r_m)
    
    # 5. Matrix Setup & Flow-Length Engine
    dir8_raw_m <- as.matrix(flow_dir_D8)
    target_idx <- which(is.na(dir8_raw_m), arr.ind = TRUE)
    
    coladd <- vector(mode = "numeric", length = nrow(dir8_raw_m)); coladd[1:length(coladd)] <- NA
    dir8_m <- cbind(coladd, dir8_raw_m, coladd)
    W_FL_m <- cbind(coladd, as.matrix(inv_CS), coladd)
    
    rowadd <- vector(mode = "numeric", length = ncol(dir8_m)); rowadd[1:length(rowadd)] <- NA
    dir8_m <- rbind(rowadd, dir8_m, rowadd)
    W_FL_m <- rbind(rowadd, W_FL_m, rowadd)
    
    loop_res <- flow_length_engine_cmp(dir8_m, W_FL_m, W_FL_m, target_idx, target_res)
    W_FL_m   <- loop_res$W_FL_m
    count    <- loop_res$count
    
    W_FL_m <- W_FL_m[-c(1, nrow(W_FL_m)), -c(1, ncol(W_FL_m))]
    ext    <- as.matrix(extent(inv_CS))
    W_FL_fin <- raster(W_FL_m, xmn = ext[1, 1], xmx = ext[1, 2], ymn = ext[2, 1], ymx = ext[2, 2], crs = crs_out)
    
    dif_WFL    <- W_FL_fin - inv_CS
    dif_WFL[dif_WFL == 0] <- NA
    mask_cells <- dif_WFL * 0 + 1
    
    Ddn <- ((W_FL_fin == 0) * inv_CS) + W_FL_fin
    crs(Ddn) <- crs_out
    
    # 6. Upslope Component (Dup)
    system(paste("mpiexec -n", mpi_cores, "DinfFlowdir -ang ang.tif -slp slp.tif -fel dem_pitfilled.tif"))
    ang <- raster("ang.tif") * target_rst
    writeRaster(ang, "ang.tif", overwrite = TRUE)
    
    system(paste("mpiexec -n", mpi_cores, "AreaDinf -ang ang.tif -sca accS.tif -wg slope_D8_rec_m.tif -nc"))
    accS <- raster("accS.tif")
    
    system(paste("mpiexec -n", mpi_cores, "AreaDinf -ang ang.tif -sca sca.tif -nc"))
    sca <- raster("sca.tif")
    
    system(paste("mpiexec -n", mpi_cores, "AreaDinf -ang ang.tif -sca accW.tif -wg W_factor.tif -nc"))
    accW <- raster("accW.tif") * target_res
    
    constant_res <- z_fill * 0 + target_res
    Acc_fin      <- sca / constant_res
    smean        <- (accS + slope_D8_r_m) / Acc_fin
    cmean        <- (accW + W_factor) / Acc_fin
    Dup          <- cmean * smean * (Acc_fin * constant_res^2)^(1/2)
    crs(Dup)     <- crs_out
    
    # 7. Final IC
    IC <- log10(Dup / Ddn) * mask_cells
    crs(IC) <- crs_out
    
    # Write final files strictly as IC_year_gaugeID
    writeRaster(IC, file.path(wd, paste0("IC_", year_str, "_", gid, ".tif")), overwrite = TRUE)
    writeRaster(Ddn, file.path(wd, paste0("Ddn_", year_str, "_", gid, ".tif")), overwrite = TRUE)
    writeRaster(Dup, file.path(wd, paste0("Dup_", year_str, "_", gid, ".tif")), overwrite = TRUE)
    
    pal <- colorRampPalette(c("blue", "lightyellow1", "red"))
    png(file.path(wd, paste0("IC_", year_str, "_", gid, ".png")), width = 800, height = 800)
    plot(IC, col = pal(100), main = paste("IC - Gauge", gid, "(", year_str, ")"))
    dev.off()
    
    clean <- getValues(IC)
    clean <- clean[!is.na(clean)]
    
    setwd(wd)
    unlink(temp_dir, recursive = TRUE)
    
    mins <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
    cat(sprintf("  -> SUCCESS: %d valid cells in %.2f minutes (Iterations: %d)\n", length(clean), mins, count))
    
    rescue_results[[k]] <- data.frame(
      gauge_id    = gid,
      year        = year_str,
      drainage    = as.integer(max_acc),
      valid_cells = length(clean),
      max_IC      = max(clean, na.rm = TRUE),
      min_IC      = min(clean, na.rm = TRUE),
      median_IC   = median(clean, na.rm = TRUE),
      SD_IC       = sd(clean, na.rm = TRUE),
      runtime_min = round(mins, 2),
      status      = "Success"
    )
    
  }, error = function(e) {
    setwd(wd)
    unlink(temp_dir, recursive = TRUE)
    cat(sprintf("  -> ERROR on Gauge %s (%s): %s\n", gid, year_str, e$message))
    rescue_results[[k]] <- data.frame(
      gauge_id    = gid,
      year        = year_str,
      drainage    = NA,
      valid_cells = 0,
      max_IC      = NA,
      min_IC      = NA,
      median_IC   = NA,
      SD_IC       = NA,
      runtime_min = 0,
      status      = paste("Error:", e$message)
    )
  })
  
  gc()
}

# ------------------------------------------------------------------------------
# Update the master CSV table
# ------------------------------------------------------------------------------
rescue_df <- dplyr::bind_rows(rescue_results)
master_csv <- file.path(wd, "stat_IC_batch_parallel.csv")

if (file.exists(master_csv) && nrow(rescue_df) > 0) {
  old_df <- read.csv(master_csv, stringsAsFactors = FALSE)
  old_df$gauge_id <- sprintf("%05d", as.numeric(old_df$gauge_id))
  
  if ("year" %in% names(old_df)) {
    old_df$year    <- as.character(old_df$year)
    rescue_df$year <- as.character(rescue_df$year)
    keys_old       <- paste(old_df$gauge_id, old_df$year, sep = "_")
    keys_rescue    <- paste(rescue_df$gauge_id, rescue_df$year, sep = "_")
    updated_df     <- rbind(old_df[!keys_old %in% keys_rescue, ], rescue_df)
  } else {
    updated_df <- rbind(old_df[!old_df$gauge_id %in% rescue_df$gauge_id, ], rescue_df)
  }
  
  updated_df <- updated_df[order(updated_df$gauge_id), ]
  write.csv(updated_df, master_csv, row.names = FALSE)
  
  cat("\n==================================================\n")
  cat(sprintf("Master summary updated: %d entries recorded in %s\n", 
              nrow(updated_df), basename(master_csv)))
  cat("==================================================\n")
}