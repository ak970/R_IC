#### South Asia Albers Equal Area Conic (ESRI:102028)

# ==============================================================================
# SCRIPT 2: STANDALONE DEBUG SCRIPT (INDIA ALBERS EQUAL-AREA PROJECTION)
# ==============================================================================

# 1. Resolve TauDEM PROJ environment conflict
if (dir.exists(system.file("proj", package = "terra"))) {
  Sys.setenv(PROJ_LIB = system.file("proj", package = "terra"))
} else {
  Sys.unsetenv("PROJ_LIB")
}

library(terra)
library(raster)

rasterOptions(tolerance = 0.05)

# Set working / output directory
debug_dir <- "F:/Peninsular India/Connectivity/IC/Debug"
if (!dir.exists(debug_dir)) dir.create(debug_dir, recursive = TRUE)
setwd(debug_dir)

# File Paths & Parameters
gauge_id     <- "04039"
year         <- "2021"
dem_dir      <- "F:/Peninsular India/Connectivity/IC/Input/DEM"
# Parent directory containing multiple subdirectories with NDVI files
w_factor_dir <- "F:/Peninsular India/Connectivity/IC/Input/NDVI"

dem_path     <- file.path(dem_dir, paste0("DEM_", gauge_id, ".tif"))

# Recursively locate NDVI_xxxx_gaugeID.tif in parent directory & subdirectories
ndvi_pattern <- paste0("NDVI_", year, "_0*", as.numeric(gauge_id), "\\.(tif|TIF)$")
ndvi_matches <- list.files(w_factor_dir, pattern = ndvi_pattern, recursive = TRUE, full.names = TRUE)

if (length(ndvi_matches) == 0) {
  stop(sprintf("No NDVI file matching year %s and gauge %s found in %s", year, gauge_id, w_factor_dir))
}
ndvi_path    <- ndvi_matches[1]

# India / South Asia Albers Equal Area Conic (ESRI:102028)
target_crs   <- "+proj=aea +lat_1=7 +lat_2=32 +lat_0=20 +lon_0=80 +x_0=2000000 +y_0=2000000 +datum=WGS84 +units=m +no_defs"
target_res   <- 100 # 100 meters
mpi_cores    <- 4

cat(sprintf("\n=================== DEBUG RUN: GAUGE %s (%s) [ALBERS EQUAL AREA] ===================\n\n", gauge_id, year))

# ------------------------------------------------------------------------------
# STEP 1: Load Raw Inputs & Reproject DEM to India Albers
# ------------------------------------------------------------------------------
cat("1. Loading raw DEM & NDVI, reprojecting to India Albers Equal-Area (100m)...\n")
dem_raw  <- terra::rast(dem_path)
ndvi_raw <- terra::rast(ndvi_path)

dem_proj_t <- terra::project(dem_raw, target_crs, res = target_res, method = "bilinear")
writeRaster(raster(dem_proj_t), "dem_raw_proj.tif", overwrite = TRUE)

# ------------------------------------------------------------------------------
# STEP 2: Pit-Filling Upfront (TauDEM PitRemove)
# ------------------------------------------------------------------------------
cat("2. Running TauDEM PitRemove to resolve all surface depressions...\n")
system(paste("mpiexec -n", mpi_cores, "PitRemove -z dem_raw_proj.tif -fel dem_pitfilled.tif"))

z_fill  <- raster("dem_pitfilled.tif")
res     <- res(z_fill)[1]
crs_out <- crs(z_fill) # Dynamically inherits India Albers

# ------------------------------------------------------------------------------
# STEP 3: Align NDVI & Compute C-factor (Cai et al., 2000)
# ------------------------------------------------------------------------------
cat("3. Computing C-factor from NDVI...\n")
wf_proj_t <- terra::project(ndvi_raw, terra::rast("dem_pitfilled.tif"), method = "bilinear")
W_factor  <- raster(wf_proj_t)

vals <- getValues(W_factor)
vals_clean <- vals[!is.na(vals) & !is.nan(vals)]
if (length(vals_clean) == 0) stop("Error: NDVI raster is empty (all NA) after reprojection!")

ndvi_min <- min(vals_clean)
ndvi_max <- max(vals_clean)

if (ndvi_max == ndvi_min) {
  f_g <- W_factor * 0
} else {
  f_g <- ((W_factor - ndvi_min) / (ndvi_max - ndvi_min)) * 100
}

f_clamped <- clamp(f_g, lower = 0.095, upper = 78.3)
C_factor  <- 0.6508 - (0.3436 * log10(f_clamped))

C_factor[f_g <= 0]     <- 1.0
C_factor[f_g > 78.3]   <- 0.001
C_factor[C_factor > 1] <- 1.0
C_factor[C_factor < 0.001] <- 0.001

W_factor <- C_factor
writeRaster(W_factor, "W_factor.tif", overwrite = TRUE)

# ------------------------------------------------------------------------------
# STEP 4: Flow Directions & True Catchment Outlet Identification
# ------------------------------------------------------------------------------
cat("4. Computing D8 Flow Directions & Locating Outlet...\n")
system(paste("mpiexec -n", mpi_cores, "D8Flowdir -p D8_flow_dir.tif -sd8 slope_D8.tif -fel dem_pitfilled.tif"))
flow_dir_D8 <- raster("D8_flow_dir.tif")
slope_D8    <- raster("slope_D8.tif")

system(paste("mpiexec -n", mpi_cores, "AreaD8 -p D8_flow_dir.tif -ad8 accD8.tif -nc"))
accD8 <- raster("accD8.tif")

# Automatic catchment outlet determination via max accumulation
max_cell   <- which.max(getValues(accD8))
snapped_xy <- xyFromCell(accD8, max_cell)
max_accum  <- max(getValues(accD8), na.rm = TRUE)

cat(sprintf("  Catchment Outlet located at X=%.1f, Y=%.1f (Total Drainage: %d cells)\n", 
            snapped_xy[1], snapped_xy[2], as.integer(max_accum)))

target_snapped      <- terra::vect(matrix(snapped_xy, ncol = 2), crs = target_crs)
target_snapped$burn <- 1

target_rst_t <- terra::rasterize(target_snapped, terra::rast(z_fill), field = "burn", background = 2)
target_rst   <- raster(target_rst_t)
target_rst   <- reclassify(target_rst, c(0.9, 1.1, NA, 1.9, 2.1, 1))

# ------------------------------------------------------------------------------
# STEP 5: Impedance & Corrected 8-Direction Downslope Traversal
# ------------------------------------------------------------------------------
cat("5. Tracing Downslope Flow Length from Outlet...\n")
flow_dir_D8[is.na(flow_dir_D8)] <- -9999
flow_dir_D8 <- flow_dir_D8 * target_rst

slope_D8_r_m <- reclassify(slope_D8, c(-2, 0.005, 0.005, 1, 10^50, 1))
writeRaster(slope_D8_r_m, "slope_D8_rec_m.tif", overwrite = TRUE)

inv_CS <- 1 / (W_factor * slope_D8_r_m)
writeRaster(inv_CS, paste0("inv_CS_", year, "_", gauge_id, ".tif"), overwrite = TRUE)

# Isolate target index before padding matrix to prevent boundary artifacts
dir8_raw_m <- as.matrix(flow_dir_D8)
target_idx <- which(is.na(dir8_raw_m), arr.ind = TRUE)

coladd <- vector(mode = "numeric", length = nrow(dir8_raw_m)); coladd[1:length(coladd)] <- NA
dir8_m <- cbind(coladd, dir8_raw_m, coladd)
W_FL_m <- cbind(coladd, as.matrix(inv_CS), coladd)

rowadd <- vector(mode = "numeric", length = ncol(dir8_m)); rowadd[1:length(rowadd)] <- NA
dir8_m <- rbind(rowadd, dir8_m, rowadd)
W_FL_m <- rbind(rowadd, W_FL_m, rowadd)

Wgt_m <- W_FL_m 

# Sinks to trace from: ONLY the true target cell
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
  
  # Clean 8-direction offset search
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

cat(sprintf("  Flow length loop finished in %d iterations!\n", count))

# Remove boundary padding
W_FL_m <- W_FL_m[-1, ]; W_FL_m <- W_FL_m[-nrow(W_FL_m), ]
dir8_m <- dir8_m[-1, ]; dir8_m <- dir8_m[-nrow(dir8_m), ]
W_FL_m <- W_FL_m[, -1]; W_FL_m <- W_FL_m[, -ncol(W_FL_m)] 
dir8_m <- dir8_m[, -1]; dir8_m <- dir8_m[, -ncol(dir8_m)]

ext <- as.matrix(extent(inv_CS))
W_FL_fin <- raster(W_FL_m, xmn = ext[1, 1], xmx = ext[1, 2], ymn = ext[2, 1], ymx = ext[2, 2], crs = crs_out)
writeRaster(W_FL_fin, "Weighted_Flow_length.tif", overwrite = TRUE)

dif_WFL <- W_FL_fin - inv_CS
dif_WFL[dif_WFL == 0] <- NA
mask_cells <- dif_WFL * 0 + 1
writeRaster(mask_cells, "mask_cells.tif", overwrite = TRUE)

Ddn <- ((W_FL_fin == 0) * inv_CS) + W_FL_fin
crs(Ddn) <- crs_out
writeRaster(Ddn, paste0("Ddn_", year, "_", gauge_id, ".tif"), overwrite = TRUE)

rm(Wgt_m, W_FL_m, dir8_m, dif_WFL, inv_CS, flow_dir_D8); gc()

# ------------------------------------------------------------------------------
# STEP 6: Upslope Component (Dup)
# ------------------------------------------------------------------------------
cat("6. Computing Upslope Component (Dup)...\n")
system(paste("mpiexec -n", mpi_cores, "DinfFlowdir -ang ang.tif -slp slp.tif -fel dem_pitfilled.tif"))
ang <- raster("ang.tif") * target_rst
writeRaster(ang, "ang.tif", overwrite = TRUE)
gc()

system(paste("mpiexec -n", mpi_cores, "AreaDinf -ang ang.tif -sca accS.tif -wg slope_D8_rec_m.tif -nc"))
accS <- raster("accS.tif")

system(paste("mpiexec -n", mpi_cores, "AreaDinf -ang ang.tif -sca sca.tif -nc"))
sca <- raster("sca.tif")

system(paste("mpiexec -n", mpi_cores, "AreaDinf -ang ang.tif -sca accW.tif -wg W_factor.tif -nc"))
accW <- raster("accW.tif") * res

constant_res <- z_fill * 0 + res
Acc_fin      <- sca / constant_res
smean        <- (accS + slope_D8_r_m) / Acc_fin
cmean        <- (accW + W_factor) / Acc_fin

Dup <- cmean * smean * (Acc_fin * constant_res^2)^(1/2)
crs(Dup) <- crs_out
writeRaster(Dup, paste0("Dup_", year, "_", gauge_id, ".tif"), overwrite = TRUE)

rm(constant_res, ang, Acc_fin, smean, accW, cmean); gc()

# ------------------------------------------------------------------------------
# STEP 7: Final Index of Connectivity (IC) Named IC_year_gaugeID
# ------------------------------------------------------------------------------
cat("7. Calculating Final Index of Connectivity with mask_cells...\n")
ratio <- Dup / Ddn
IC    <- log10(ratio)
IC    <- IC * mask_cells
crs(IC) <- crs_out

# Preview Plot strictly named IC_year_gaugeID.png
pal <- colorRampPalette(c("blue", "lightyellow1", "red"))
png(paste0("IC_", year, "_", gauge_id, ".png"), width = 800, height = 800)
plot(IC, col = pal(100), main = paste("Index of Connectivity (IC) - Catchment", gauge_id, "(", year, ")"))
dev.off()

# Final Raster strictly named IC_year_gaugeID.tif
writeRaster(IC, paste0("IC_", year, "_", gauge_id, ".tif"), overwrite = TRUE)

cat("\n=================== FINAL STATUS ===================\n")
valid_mask  <- sum(!is.na(getValues(mask_cells)))
valid_ic    <- sum(!is.na(getValues(IC)))
total_cells <- length(getValues(IC))

cat(sprintf("Valid Cells in mask_cells : %d (%.1f%% of DEM)\n", valid_mask, (valid_mask / total_cells) * 100))
cat(sprintf("Valid Cells in Final IC   : %d (%.1f%% of DEM)\n", valid_ic, (valid_ic / total_cells) * 100))
cat("IC Range                  : [", min(getValues(IC), na.rm = TRUE), "to", max(getValues(IC), na.rm = TRUE), "]\n")
cat("Debug run successfully completed!\n")