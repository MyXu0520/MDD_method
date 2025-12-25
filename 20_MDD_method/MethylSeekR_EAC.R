###############################   Compute EAC PMDs  ###################################
suppressPackageStartupMessages({
  library(MethylSeekR)
  library(GenomicRanges)
  library(data.table)
})

# Prevent certain functions from requiring X11 for png()
options(bitmapType = "cairo")

# ========== Set Paths ==========
input_dir  <- "/data3/xumy_PMD/data"
output_dir <- file.path(input_dir, "EAC_GEJ")

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Use only these 12 samples:
sample_names <- c(
  paste0("EAC_", c(1:4,6)),
  paste0("GEJ_", 1:7)
)
# Corresponding input filenames: EAC_1methyCov.bed ... GEJ_7methyCov.bed
files <- file.path(input_dir, paste0(sample_names, "methyCov.bed"))

# Check if all files exist
missing <- !file.exists(files)
if (any(missing)) {
  stop("The following files were not found. Please check file names or paths:\n",
       paste(files[missing], collapse = "\n"))
}

# ========== Set Parameters (modify as needed) ==========
minCov    <- 5      # Minimum read coverage for CpG sites
num.cores <- 4      # Number of cores for HMM (set to 1 if no parallel processing)

# ========== Helper Function: Load Data and Build GRanges Object ==========
# Current file format: V1=chr, V2=pos, V3=beta, V4=meth, V5=unmeth
load_methyl_gr <- function(f, minCov = 5) {
  message("Loading file: ", f)
  dt <- fread(f)
  
  if (ncol(dt) < 5) {
    stop("File has < 5 columns. Check methyCov.bed format: ", f)
  }
  setnames(dt, 1:5, c("chr", "pos", "beta", "meth", "unmeth"))
  
  dt[, T := meth + unmeth]  # Total reads
  dt[, M := meth]           # Methylated reads
  
  dt <- dt[T >= minCov]     # Filter by minimum coverage
  if (nrow(dt) == 0) {
    stop("No CpG sites remain after filtering. Check minCov or raw data: ", f)
  }
  
  gr <- GRanges(
    seqnames = dt$chr,
    ranges   = IRanges(start = dt$pos, end = dt$pos),
    strand   = "*",
    T        = as.integer(dt$T),  # Total coverage
    M        = as.integer(dt$M)   # Methylated counts
  )
  
  sort(gr)
}

# ========== Main Loop: Run segmentPMDs for Each Sample ==========
for (i in seq_along(sample_names)) {
  sample_name <- sample_names[i]
  f           <- files[i]
  
  message("===== Processing sample: ", sample_name, " =====")
  
  meth.gr <- load_methyl_gr(f, minCov = minCov)
  
  # Estimate seqLengths from the current data
  seqLengths_vec <- tapply(end(meth.gr), as.character(seqnames(meth.gr)), max)
  seqLengths     <- as.numeric(seqLengths_vec)
  names(seqLengths) <- names(seqLengths_vec)
  
  # Select a chromosome for HMM training
  chr.sel <- names(seqLengths)[1]
  message("Chromosome selected for HMM training: ", chr.sel)
  
  # ⭐ Key modification: Force plotting to PDF instead of screen
  pdf_file <- file.path(output_dir, paste0(sample_name, "_PMD_HMM.pdf"))
  
  PMDsegments.gr <- segmentPMDs(
    m           = meth.gr,
    chr.sel     = chr.sel,
    pdfFilename = pdf_file,   # Not NULL anymore
    seqLengths  = seqLengths,
    num.cores   = num.cores
  )
  
  # Keep only PMD regions (type == "PMD")
  if (!"type" %in% colnames(mcols(PMDsegments.gr))) {
    stop("No 'type' column in PMDsegments.gr. Use str(PMDsegments.gr) to check actual column names.")
  }
  PMD_only <- PMDsegments.gr[mcols(PMDsegments.gr)$type == "PMD"]
  
  if (length(PMD_only) == 0) {
    warning("Sample ", sample_name, " has no detected PMD regions.")
  }
  
  # Output BED format: chr, start(0-based), end, nCG
  pmd_dt <- data.table(
    chr   = as.character(seqnames(PMD_only)),
    start = start(PMD_only) - 1L,   # Convert to 0-based for BED
    end   = end(PMD_only),          # BED uses closed end
    nCG   = mcols(PMD_only)$nCG
  )
  
  out_bed <- file.path(output_dir, paste0(sample_name, "_PMD.bed"))
  fwrite(pmd_dt, out_bed, sep = "\t", col.names = TRUE)
  message("PMD results saved to: ", out_bed)
}

message("PMD calculation completed for all 12 samples.")


###############################   Compute Core EAC PMDs  ###################################
suppressPackageStartupMessages({
  library(GenomicRanges)
  library(data.table)
})

# ========== Set Paths and Sample Names ==========
pmd_dir   <- "/data3/xumy_PMD/data/EAC_GEJ"
out_file  <- file.path(pmd_dir, "core_PMD_all12.bed")

sample_names <- c(
  paste0("EAC_", c(1:4,6)),
  paste0("GEJ_", 1:7)
)

pmd_files <- file.path(pmd_dir, paste0(sample_names, "_PMD.bed"))

# Check if all PMD files exist
missing <- !file.exists(pmd_files)
if (any(missing)) {
  stop("The following PMD files were not found. Check paths or names:\n",
       paste(pmd_files[missing], collapse = "\n"))
}

# ========== Read BED Files and Convert to GRanges ==========
# PMD.bed format: chr, start(0-based), end, nCG
read_pmd_gr <- function(f) {
  message("Loading PMD file: ", f)
  dt <- fread(f)
  if (ncol(dt) < 3) {
    stop("PMD BED file should have at least 3 columns: chr, start, end: ", f)
  }
  setnames(dt, 1:3, c("chr", "start", "end"))
  
  GRanges(
    seqnames = dt$chr,
    ranges   = IRanges(
      start = dt$start + 1L,   # Convert from 0-based BED to 1-based GRanges
      end   = dt$end
    ),
    strand   = "*"
  )
}

gr_list <- lapply(pmd_files, read_pmd_gr)

# ========== Compute Intersection of PMDs Across All 12 Samples (Core PMDs) ==========
message("Calculating intersection of PMDs across 12 samples (core PMDs)...")

core_pmd <- Reduce(
  f = function(x, y) GenomicRanges::intersect(x, y, ignore.strand = TRUE),
  x = gr_list
)

# Optional: Reduce to merge overlapping intervals (usually not needed after intersection)
core_pmd <- reduce(core_pmd, ignore.strand = TRUE)

if (length(core_pmd) == 0L) {
  warning("Note: No common PMD regions found across all 12 samples. core_pmd is empty.")
}

# ========== Output Core PMDs in BED Format ==========
core_dt <- data.table(
  chr   = as.character(seqnames(core_pmd)),
  start = start(core_pmd) - 1L,   # Convert back to 0-based BED
  end   = end(core_pmd)
)

fwrite(core_dt, out_file, sep = "\t", col.names = TRUE)
message("Core PMDs saved to: ", out_file)