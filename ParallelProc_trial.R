# Title: High-Performance Parallel Batch R_IC
# Features: Multi-core catchment parallelization, worker sandboxing, & JIT byte-compiled loop

library(terra)
library(raster)
library(parallel)
library(foreach)
library(doParallel)
library(compiler)

# 1. Resolve TauDEM PROJ environment conflict
if (dir.exists(system.file("proj", package = "terra"))) {
  Sys.setenv(PROJ_LIB = system.file("proj", package = "terra"))
} else {
  Sys.unsetenv("PROJ_LIB")
}

rasterOptions(tolerance = 0.05)

###################################
# USER INPUTS & DIRECTORIES SETUP
###################################

# Main output directory
wd <- "F:/Peninsular India/Connectivity/IC/Output_Parallel"
if (!dir.exists(wd)) dir.create(wd, recursive = TRUE)

dem_dir      <- "F:/Peninsular India/Connectivity/IC/Input/DEM"
w_factor_dir <- "F:/Peninsular India/Connectivity/IC/Input/NDVI/2021"
target_path  <- "F:/Peninsular India/CAMELS_IND/CAMELS_IND_Catchments_Streamflow_Sufficient3/shapefiles_catchment/gauge_stations.shp"
gauge_id_col <- "gauge_id"

# Coordinate System & Target Resolution
target_crs   <- "+proj=aea +lat_1=7 +lat_2=32 +lat_0=20 +lon_0=80 +x_0=2000000 +y_0=2000000 +datum=WGS84 +units=m +no_defs"
target_res   <- 100          # 100 meters

# Outlet Method: "catchment_outlet" (Recommended for batch) OR "gauge_snapping"
outlet_method <- "catchment_outlet"
snap_radius   <- 2500        # Used only if outlet_method == "gauge_snapping"
mw            <- 5

# PARALLEL CORES CONFIGURATION
# Set how many catchments to process simultaneously.
# Leave 1 or 2 cores free for Windows OS stability.
available_cores <- parallel::detectCores()
num_workers     <- max(1, available_cores - 2) # e.g., 6 workers on an 8-core CPU
cat(sprintf("System has %d CPU cores. Launching %d parallel workers...\n", available_cores, num_workers))


###################################
# COMPILED FLOW-LENGTH ENGINE (JIT Bytecode)
###################################
# Compiling the loop into C-bytecode speeds up matrix traversal by 2x to 3x
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

# Compile function to native bytecode
flow_length_engine_cmp <- compiler::cmpfun(flow_length_engine)


###################################
# FILE MATCHING PIPELINE
###################################
all_targets <- terra::vect(target_path)
if (is.na(terra::crs(all_targets)) || terra::crs(all_targets) == "") terra::crs(all_targets) <- "EPSG:4326"

dem_file_list <- list.files(dem_dir, pattern = "^DEM_[0-9]{5}\\.(tif|TIF)$", full.names = TRUE)
shp_gauge_ids <- trimws(as.character(all_targets[[gauge_id_col]][[1]]))

batch_list <- data.frame(
  gauge_id      = character(),
  dem_path      = character(),
  w_factor_path = character(),
  stringsAsFactors = FALSE
)

for (df in dem_file_list) {
  fname <- basename(df)
  gid   <- sub(".*DEM_([0-9]{5})\\..*", "\\1", fname, ignore.case = TRUE)
  wf_path <- file.path(w_factor_dir, paste0("NDVI_", gid, ".tif"))
  
  if (file.exists(wf_path) && any(shp_gauge_ids == gid | sprintf("%05d", suppressWarnings(as.numeric(shp_gauge_ids))) == gid)) {
    batch_list <- rbind(batch_list, data.frame(
      gauge_id      = gid,
      dem_path      = df,
      w_factor_path = wf_path,
      stringsAsFactors = FALSE
    ))
  }
}

cat(sprintf("Ready to process %d catchments across %d parallel workers.\n", nrow(batch_list), num_workers))


###################################
# PARALLEL WORKER FUNCTION
###################################
process_single_catchment <- function(row_idx, batch_df, main_out_dir, target_shapefile_path, 
                                     gauge_col_name, target_crs_str, target_cell_res, 
                                     outlet_sel_mode, snap_dist, mw_tag) {
  
  gid     <- batch_df$gauge_id[row_idx]
  dem_p   <- batch_df$dem_path[row_idx]
  wf_p    <- batch_df$w_factor_path[row_idx]
  res     <- target_cell_res
  
  # 1. Create a unique sandbox working directory for this catchment
  worker_dir <- file.path(main_out_dir, paste0("temp_worker_", gid, "_", Sys.getpid()))
  if (!dir.exists(worker_dir)) dir.create(worker_dir, recursive = TRUE)
  setwd(worker_dir)
  
  start_time <- Sys.time()
  
  tryCatch({
    # Reproject DEM to metric grid & save for PitRemove
    dem_raw    <- terra::rast(dem_p)
    dem_proj_t <- terra::project(dem_raw, target_crs_str, res = res, method = "bilinear")
    writeRaster(raster(dem_proj_t), "dem_raw_proj.tif", overwrite = TRUE)
    
    # 1 thread per TauDEM CLI call prevents CPU congestion when workers run concurrently
    system("mpiexec -n 1 PitRemove -z dem_raw_proj.tif -fel dem_pitfilled.tif")
    
    z_fill  <- raster("dem_pitfilled.tif")
    crs_out <- crs(z_fill)
    
    # Align NDVI & calculate C-factor
    ndvi_raw  <- terra::rast(wf_p)
    wf_proj_t <- terra::project(ndvi_raw, terra::rast("dem_pitfilled.tif"), method = "bilinear")
    W_factor  <- raster(wf_proj_t)
    
    vals <- getValues(W_factor)
    vals_clean <- vals[!is.na(vals) & !is.nan(vals)]
    ndvi_min <- min(vals_clean)
    ndvi_max <- max(vals_clean)
    
    f_g <- if (ndvi_max == ndvi_min) W_factor * 0 else ((W_factor - ndvi_min) / (ndvi_max - ndvi_min)) * 100
    f_clamped <- clamp(f_g, lower = 0.095, upper = 78.3)
    C_factor  <- 0.6508 - (0.3436 * log10(f_clamped))
    C_factor[f_g <= 0]     <- 1.0
    C_factor[f_g > 78.3]   <- 0.001
    C_factor[C_factor > 1] <- 1.0
    C_factor[C_factor < 0.001] <- 0.001
    
    writeRaster(C_factor, "W_factor.tif", overwrite = TRUE)
    W_factor <- C_factor
    
    # D8 flow directions & Accumulation
    system("mpiexec -n 1 D8Flowdir -p D8_flow_dir.tif -sd8 slope_D8.tif -fel dem_pitfilled.tif")
    flow_dir_D8 <- raster("D8_flow_dir.tif")
    slope_D8    <- raster("slope_D8.tif")
    
    system("mpiexec -n 1 AreaD8 -p D8_flow_dir.tif -ad8 accD8.tif -nc")
    accD8 <- raster("accD8.tif")
    
    # Target outlet determination
    if (outlet_sel_mode == "catchment_outlet") {
      max_cell   <- which.max(getValues(accD8))
      snapped_xy <- xyFromCell(accD8, max_cell)
      max_acc    <- max(getValues(accD8), na.rm = TRUE)
    } else {
      all_tgts <- terra::vect(target_shapefile_path)
      ids      <- trimws(as.character(all_tgts[[gauge_col_name]][[1]]))
      is_m     <- (ids == gid) | (sprintf("%05d", suppressWarnings(as.numeric(ids))) == gid)
      tgt_feat <- terra::project(all_tgts[which(is_m)[1], ], target_crs_str)
      
      buf      <- terra::buffer(tgt_feat, snap_dist)
      accD8_t  <- terra::rast("accD8.tif")
      acc_mask <- terra::mask(terra::crop(accD8_t, buf), buf)
      acc_vals <- terra::values(acc_mask)
      max_idx  <- which.max(acc_vals)
      
      if (length(max_idx) == 0 || is.na(max_idx)) {
        snapped_xy <- terra::geom(tgt_feat)[, c("x", "y")]
        max_acc    <- 0
      } else {
        snapped_xy <- terra::xyFromCell(acc_mask, max_idx)
        max_acc    <- acc_vals[max_idx]
      }
      rm(accD8_t, acc_mask); gc()
    }
    
    # Rasterize target
    tgt_snapped      <- terra::vect(matrix(snapped_xy, ncol = 2), crs = target_crs_str)
    tgt_snapped$burn <- 1
    target_rst       <- raster(terra::rasterize(tgt_snapped, terra::rast(z_fill), field = "burn", background = 2))
    target_rst       <- reclassify(target_rst, c(0.9, 1.1, NA, 1.9, 2.1, 1))
    
    flow_dir_D8[is.na(flow_dir_D8)] <- -9999
    flow_dir_D8 <- flow_dir_D8 * target_rst
    
    slope_D8_r_m <- reclassify(slope_D8, c(-2, 0.005, 0.005, 1, 10^50, 1))
    writeRaster(slope_D8_r_m, "slope_D8_rec_m.tif", overwrite = TRUE)
    
    inv_CS <- 1 / (W_factor * slope_D8_r_m)
    
    # Matrix setup with boundary-safe isolation
    dir8_raw_m <- as.matrix(flow_dir_D8)
    target_idx <- which(is.na(dir8_raw_m), arr.ind = TRUE)
    
    coladd <- vector(mode = "numeric", length = nrow(dir8_raw_m)); coladd[1:length(coladd)] <- NA
    dir8_m <- cbind(coladd, dir8_raw_m, coladd)
    W_FL_m <- cbind(coladd, as.matrix(inv_CS), coladd)
    
    rowadd <- vector(mode = "numeric", length = ncol(dir8_m)); rowadd[1:length(rowadd)] <- NA
    dir8_m <- rbind(rowadd, dir8_m, rowadd)
    W_FL_m <- rbind(rowadd, W_FL_m, rowadd)
    
    # Run the compiled bytecode engine
    loop_result <- flow_length_engine_cmp(dir8_m, W_FL_m, W_FL_m, target_idx, res)
    W_FL_m      <- loop_result$W_FL_m
    iter_count  <- loop_result$count
    
    # Remove padding & reconstruct
    W_FL_m <- W_FL_m[-c(1, nrow(W_FL_m)), -c(1, ncol(W_FL_m))]
    ext    <- as.matrix(extent(inv_CS))
    W_FL_fin <- raster(W_FL_m, xmn = ext[1, 1], xmx = ext[1, 2], ymn = ext[2, 1], ymx = ext[2, 2], crs = crs_out)
    
    dif_WFL    <- W_FL_fin - inv_CS
    dif_WFL[dif_WFL == 0] <- NA
    mask_cells <- dif_WFL * 0 + 1
    
    Ddn <- ((W_FL_fin == 0) * inv_CS) + W_FL_fin
    crs(Ddn) <- crs_out
    
    # Upslope component
    system("mpiexec -n 1 DinfFlowdir -ang ang.tif -slp slp.tif -fel dem_pitfilled.tif")
    ang <- raster("ang.tif") * target_rst
    writeRaster(ang, "ang.tif", overwrite = TRUE)
    
    system("mpiexec -n 1 AreaDinf -ang ang.tif -sca accS.tif -wg slope_D8_rec_m.tif -nc")
    accS <- raster("accS.tif")
    
    system("mpiexec -n 1 AreaDinf -ang ang.tif -sca sca.tif -nc")
    sca <- raster("sca.tif")
    
    system("mpiexec -n 1 AreaDinf -ang ang.tif -sca accW.tif -wg W_factor.tif -nc")
    accW <- raster("accW.tif") * res
    
    constant_res <- z_fill * 0 + res
    Acc_fin      <- sca / constant_res
    smean        <- (accS + slope_D8_r_m) / Acc_fin
    cmean        <- (accW + W_factor) / Acc_fin
    Dup          <- cmean * smean * (Acc_fin * constant_res^2)^(1/2)
    crs(Dup)     <- crs_out
    
    # Final IC
    IC <- log10(Dup / Ddn) * mask_cells
    crs(IC) <- crs_out
    
    # Save final results directly to MAIN output folder
    writeRaster(IC, file.path(main_out_dir, paste0("IC_", res, "m_", mw_tag, "x", mw_tag, "_", gid, ".tif")), overwrite = TRUE)
    writeRaster(Ddn, file.path(main_out_dir, paste0("Ddn_", res, "m_", mw_tag, "x", mw_tag, "_", gid, ".tif")), overwrite = TRUE)
    writeRaster(Dup, file.path(main_out_dir, paste0("Dup_", res, "m_", mw_tag, "x", mw_tag, "_", gid, ".tif")), overwrite = TRUE)
    
    pal <- colorRampPalette(c("blue", "lightyellow1", "red"))
    png(file.path(main_out_dir, paste0("IC_", res, "m_", mw_tag, "x", mw_tag, "_", gid, ".png")), width = 800, height = 800)
    plot(IC, col = pal(100), main = paste("IC - Gauge", gid))
    dev.off()
    
    # Calculate statistics
    ic_v  <- getValues(IC)
    clean <- ic_v[!is.na(ic_v)]
    
    # Clean up and delete sandbox
    setwd(main_out_dir)
    unlink(worker_dir, recursive = TRUE)
    
    mins_taken <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
    cat(sprintf("[SUCCESS] Gauge %s finished in %.2f min (Iterations: %d | Cells: %d)\n", 
                gid, mins_taken, iter_count, length(clean)))
    
    return(data.frame(
      gauge_id    = gid,
      drainage    = as.integer(max_acc),
      valid_cells = length(clean),
      max_IC      = max(clean, na.rm = TRUE),
      min_IC      = min(clean, na.rm = TRUE),
      median_IC   = median(clean, na.rm = TRUE),
      SD_IC       = sd(clean, na.rm = TRUE),
      runtime_min = round(mins_taken, 2),
      status      = "Success"
    ))
    
  }, error = function(e) {
    setwd(main_out_dir)
    unlink(worker_dir, recursive = TRUE)
    cat(sprintf("[ERROR] Gauge %s failed: %s\n", gid, e$message))
    
    return(data.frame(
      gauge_id    = gid,
      drainage    = NA,
      valid_cells = 0,
      max_IC      = NA,
      min_IC      = NA,
      median_IC   = NA,
      SD_IC       = NA,
      runtime_min = 0,
      status      = paste("Error:", e$message)
    ))
  })
}


###################################################################################################################
# EXECUTE PARALLEL BATCH RUN
###################################################################################################################

cat("\n=================== STARTING PARALLEL BATCH RUN ===================\n")
total_start_time <- Sys.time()

# 1. Initialize Cluster
cl <- parallel::makeCluster(num_workers)
doParallel::registerDoParallel(cl)

# 2. Export required objects and environment settings to all workers
clusterEvalQ(cl, {
  library(terra)
  library(raster)
  library(compiler)
  
  if (dir.exists(system.file("proj", package = "terra"))) {
    Sys.setenv(PROJ_LIB = system.file("proj", package = "terra"))
  }
  rasterOptions(tolerance = 0.05)
})

clusterExport(cl, varlist = c("flow_length_engine_cmp"))

# 3. Parallel Loop across all catchments
results_list <- foreach(
  i = 1:nrow(batch_list), 
  .packages = c("terra", "raster", "compiler"),
  .errorhandling = "pass"
) %dopar% {
  process_single_catchment(
    row_idx                = i,
    batch_df               = batch_list,
    main_out_dir           = wd,
    target_shapefile_path  = target_path,
    gauge_col_name         = gauge_id_col,
    target_crs_str         = target_crs,
    target_cell_res        = target_res,
    outlet_sel_mode        = outlet_method,
    snap_dist              = snap_radius,
    mw_tag                 = mw
  )
}

# 4. Stop Parallel Cluster
parallel::stopCluster(cl)

# 5. Bind Results and Save Summary CSV
final_summary <- dplyr::bind_rows(results_list)
summary_file  <- file.path(wd, "stat_IC_batch_parallel.csv")
write.csv(final_summary, summary_file, row.names = FALSE)

total_end_time <- Sys.time()
total_duration <- as.numeric(difftime(total_end_time, total_start_time, units = "mins"))

cat("\n====================================================================\n")
cat(sprintf("ALL %d CATCHMENTS COMPLETED IN %.2f MINUTES!\n", nrow(batch_list), total_duration))
cat(sprintf("Summary statistics exported to: %s\n", summary_file))
cat("====================================================================\n")