## max. code width ============================================================

load_and_time <- function(
    file_path,
    reader_fun,
    ...,
    verbose = TRUE
) {
  
  # Check file exists
  if (!file.exists(file_path)) {
    stop("File does not exist: ", file_path)
  }
  
  # File size in MB
  file_size_mb <- file.info(file_path)$size / (1024^2)
  
  # Timing
  start_time <- Sys.time()
  
  data <- reader_fun(file_path, ...)
  
  end_time <- Sys.time()
  
  elapsed_sec <- as.numeric(
    difftime(end_time, start_time, units = "secs")
  )
  
  result <- list(
    data = data,
    file = file_path,
    size_mb = round(file_size_mb, 3),
    load_time_sec = round(elapsed_sec, 3)
  )
  
  if (verbose) {
    cat("File loaded:", basename(file_path), "\n")
    cat("- object size (MB):", result$size_mb, "\n")
    cat("- load time (sec):", result$load_time_sec, "\n")
  }
  
  return(result)
}

## Experimental description ===================================================

## Title: Arc Institute in vivo pertubation screens
## Author: DKatzlberger
## Date: 21.01.2026

## Libraries =================================================================

# Scource R scripts
source("R/utils.R")
source("R/JB22_pseudobulk_utils.R")

# Activate renv + load libs.
.script_startup(
  c(
    "Seurat"
  )
)

## Directories ================================================================

# Shared directory (location: vscratch)
JB22_pb_dir <- "/vscratch/wes/arc/share"

# Define paths
res_dir <- file.path("results", "JB22_Trutschnig")
fig_dir <- file.path(res_dir, "fig")
tab_dir <- file.path(res_dir, "tab")

# Create output dirs
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

## Datasets ===================================================================

i <- readRDS(file.path(JB22_pb_dir, "JB22_Spleen/rna_integrated.rds"))

i
str(i)
