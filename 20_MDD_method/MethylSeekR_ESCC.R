############################################################
# Calculate PMDs for ESCC (Squamous Cell Carcinoma) Samples
############################################################
suppressPackageStartupMessages({
  library(MethylSeekR)    # For PMD segmentation
  library(GenomicRanges)  # For genomic interval operations
  library(data.table)     # For fast file I/O
})

# Prevent X11 dependency for png() calls
options(bitmapType = "cairo")

# ========== Path Configuration ==========
input_dir  <- "/data3/xumy_PMD/data"
output_dir <- file.path(input_dir, "ESCC")

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Sample list (22 samples total)
sample_names <- c(
  paste0("ESCC_", c(1:17,19:22))
)
# Corresponding input files: ESCC_1methyCov.bed ... ESCC_22methyCov.bed
files <- file.path(input_dir, paste0(sample_names, "methyCov.bed"))

# Verify all input files exist
missing <- !file.exists(files)
if (any(missing)) {
  stop("The following files are missing. Please check file names and paths:\n",
       paste(files[missing], collapse = "\n"))
}

# ========== Adjustable Parameters ==========
minCov    <- 5      # Minimum read coverage required per CpG
num.cores <- 4      # Number of cores for HMM (set to 1 if no parallelization)

# ========== Helper Function: Load Data and Create GRanges ==========
# Expected file format: V1=chromosome, V2=position, V3=beta, V4=methylated, V5=unmethylated
load_methyl_gr <- function(f, minCov = 5) {
  message("Loading file: ", f)
  dt <- fread(f)
  
  # Validate column count
  if (ncol(dt) < 5) {
    stop("File has fewer than 5 columns. Check format of: ", f)
  }
  setnames(dt, 1:5, c("chr", "pos", "beta", "meth", "unmeth"))
  
  # Calculate total coverage and methylated counts
  dt[, T := meth + unmeth]  # Total coverage
  dt[, M := meth]           # Methylated counts
  
  # Filter by minimum coverage
  dt <- dt[T >= minCov]
  if (nrow(dt) == 0) {
    stop("No CpGs remain after filtering. Check minCov or input data: ", f)
  }
  
  # Create GRanges object
  gr <- GRanges(
    seqnames = dt$chr,
    ranges   = IRanges(start = dt$pos, end = dt$pos),
    strand   = "*",
    T        = as.integer(dt$T),  # Total coverage
    M        = as.integer(dt$M)  # Methylated count
  )
  
  sort(gr)  # Return sorted genomic ranges
}

# ========== Main Loop: PMD Segmentation for Each Sample ==========
for (i in seq_along(sample_names)) {
  sample_name <- sample_names[i]
  f           <- files[i]
  
  message("===== Processing sample: ", sample_name, " =====")
  
  # Load methylation data
  meth.gr <- load_methyl_gr(f, minCov = minCov)
  
  # Estimate chromosome lengths from data
  seqLengths_vec <- tapply(end(meth.gr), as.character(seqnames(meth.gr)), max)
  seqLengths     <- as.numeric(seqLengths_vec)
  names(seqLengths) <- names(seqLengths_vec)
  
  # Select first chromosome for HMM training
  chr.sel <- names(seqLengths)[1]
  message("Chromosome selected for HMM training: ", chr.sel)
  
  # Run PMD segmentation (plots will be saved to PDF)
  pdf_file <- file.path(output_dir, paste0(sample_name, "_PMD_HMM.pdf"))
  
  PMDsegments.gr <- segmentPMDs(
    m           = meth.gr,
    chr.sel     = chr.sel,
    pdfFilename = pdf_file,   # Save plots to PDF instead of screen
    seqLengths  = seqLengths,
    num.cores   = num.cores
  )
  
  # Extract only PMD regions (excluding non-PMD segments)
  if (!"type" %in% colnames(mcols(PMDsegments.gr))) {
    stop("PMDsegments.gr lacks 'type' column. Check structure with str(PMDsegments.gr).")
  }
  PMD_only <- PMDsegments.gr[mcols(PMDsegments.gr)$type == "PMD"]
  
  if (length(PMD_only) == 0) {
    warning("Sample ", sample_name, " has no detected PMD regions.")
  }
  
  # Format as BED: chromosome, 0-based start, end, CpG count
  pmd_dt <- data.table(
    chr   = as.character(seqnames(PMD_only)),
    start = start(PMD_only) - 1L,   # Convert to 0-based for BED
    end   = end(PMD_only),          # BED uses 1-based end
    nCG   = mcols(PMD_only)$nCG     # Number of CpGs in segment
  )
  
  # Save to file
  out_bed <- file.path(output_dir, paste0(sample_name, "_PMD.bed"))
  fwrite(pmd_dt, out_bed, sep = "\t", col.names = TRUE)
  message("PMD results saved to: ", out_bed)
}

message("PMD calculation completed for all 22 samples.")



############################################################
# Calculate Core PMD (Intersection Across All Samples)
############################################################
suppressPackageStartupMessages({
  library(GenomicRanges)  # For genomic interval operations
  library(data.table)     # For fast file I/O
})

# ========== Paths and Sample Names ==========
pmd_dir   <- "/data3/xumy_PMD/data/ESCC"
out_file  <- file.path(pmd_dir, "core_PMD_all22.bed")

sample_names <- c(
  paste0("ESCC_", c(1:17,19:22))
)

pmd_files <- file.path(pmd_dir, paste0(sample_names, "_PMD.bed"))

# Verify PMD files exist
missing <- !file.exists(pmd_files)
if (any(missing)) {
  stop("The following PMD files are missing:\n",
       paste(pmd_files[missing], collapse = "\n"))
}

# ========== Helper Function: Read PMD BED -> GRanges ==========
# PMD BED format: chr, 0-based start, end, nCG
read_pmd_gr <- function(f) {
  message("Loading PMD file: ", f)
  dt <- fread(f)
  if (ncol(dt) < 3) {
    stop("PMD BED must have at least 3 columns: chr, start, end: ", f)
  }
  setnames(dt, 1:3, c("chr", "start", "end"))
  
  # Convert 0-based BED to 1-based GRanges
  GRanges(
    seqnames = dt$chr,
    ranges   = IRanges(
      start = dt$start + 1L,   # Convert to 1-based
      end   = dt$end
    ),
    strand   = "*"
  )
}

# Load all PMD sets
gr_list <- lapply(pmd_files, read_pmd_gr)

# ========== Calculate Intersection Across All 22 Samples ==========
message("Calculating intersection of PMDs across 22 samples (core PMD)...")

# Compute sequential intersection: intersect(..., intersect(intersect(s1, s2), s3), ...)
core_pmd <- Reduce(
  f = function(x, y) GenomicRanges::intersect(x, y, ignore.strand = TRUE),
  x = gr_list
)

# Merge adjacent/overlapping intervals (optional, usually clean after intersection)
core_pmd <- reduce(core_pmd, ignore.strand = TRUE)

if (length(core_pmd) == 0L) {
  warning("Note: No common PMD regions found across 22 samples. Core PMD is empty.")
}

# ========== Output Core PMD as BED ==========
core_dt <- data.table(
  chr   = as.character(seqnames(core_pmd)),
  start = start(core_pmd) - 1L,   # Convert back to 0-based BED
  end   = end(core_pmd)
)

fwrite(core_dt, out_file, sep = "\t", col.names = TRUE)
message("Core PMD saved to: ", out_file)