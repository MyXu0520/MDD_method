suppressPackageStartupMessages({
  library(MethylSeekR)
  library(GenomicRanges)
  library(rtracklayer)
  library(data.table)
  library(BSgenome.Hsapiens.UCSC.hg38)
})

# =============================================
# 1. Paths and samples
# =============================================
data_dir <- "/data3/xumy_PMD/data"
output_dir <- "/data3/xumy_PMD/MethylSeekR_out"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# EAC (GEJ)
eac_files <- list.files(data_dir, pattern="GEJ", full=TRUE)
# ESCC
escc_files <- list.files(data_dir, pattern="ESCC", full=TRUE)

cat("EAC sample count:", length(eac_files), "\n")
cat("ESCC sample count:", length(escc_files), "\n")

# hg38 sequence lengths
sLengths <- seqlengths(BSgenome.Hsapiens.UCSC.hg38)

# =============================================
# 2. Prepare a function to read BED files
# Input format must be: chr, start, T, M
# T = mC + uC, M = mC
# =============================================
read_to_methylome <- function(file) {
  cat("Reading:", basename(file), "\n")
  dt <- fread(file, header = FALSE)
  colnames(dt)[1:5] <- c("chr","start","meth","mC","uC")
  
  dt <- dt[!is.na(mC) & !is.na(uC)]
  dt$T <- dt$mC + dt$uC
  dt$M <- dt$mC
  
  dt$end <- dt$start
  
  GRanges(
    seqnames = dt$chr,
    ranges   = IRanges(dt$start, dt$end),
    T = dt$T,
    M = dt$M
  )
}

# =============================================
# 3. Use MethylSeekR's segmentPMDs
# =============================================
run_methylseekr_pmd <- function(gr, sample_id, out_dir) {
  
  all_pmds <- GRanges()
  
  for (chr in unique(seqnames(gr))) {
    cat("  Analyzing chromosome:", chr, "\n")
    chr_gr <- gr[seqnames(gr)==chr]
    
    if (length(chr_gr) < 500) next  # Too few to train HMM
    
    pdf_file <- file.path(out_dir, paste0(sample_id, "_", chr, "_HMM.pdf"))
    
    seg <- segmentPMDs(
      m = chr_gr,
      chr.sel = as.character(chr),
      pdfFilename = pdf_file,
      seqLengths = sLengths
    )
    
    all_pmds <- c(all_pmds, seg[seg$type == "PMD"])
  }
  
  if (length(all_pmds) > 0)
    all_pmds <- reduce(all_pmds)
  
  # Save results
  savePMDSegments(
    PMDs = all_pmds,
    GRangesFilename = file.path(out_dir, paste0(sample_id, "_PMDs.rds")),
    TableFilename = file.path(out_dir, paste0(sample_id, "_PMDs.tab"))
  )
  
  return(all_pmds)
}

# =============================================
# 4. Process a group
# =============================================
process_group <- function(files, group_name){
  cat("\n====== Processing group:", group_name, "======\n")
  
  group_dir <- file.path(output_dir, group_name)
  dir.create(group_dir, showWarnings = FALSE)
  
  group_pmds <- GRanges()
  
  for (f in files){
    sid <- tools::file_path_sans_ext(basename(f))  
    gr <- read_to_methylome(f)
    pmd <- run_methylseekr_pmd(gr, sid, group_dir)
    group_pmds <- c(group_pmds, pmd)
  }
  
  if (length(group_pmds) > 0)
    group_pmds <- reduce(group_pmds)
  
  # Save overall PMDs for this group
  savePMDSegments(
    PMDs = group_pmds,
    GRangesFilename = file.path(output_dir, paste0(group_name, "_group_PMDS.rds")),
    TableFilename = file.path(output_dir, paste0(group_name, "_group_PMDS.tab"))
  )
  
  return(group_pmds)
}

# =============================================
# 5. Execute analysis
# =============================================
cat("Starting MethylSeekR PMD analysis...\n")

EAC_pmd <- process_group(eac_files, "EAC")
ESCC_pmd <- process_group(escc_files, "ESCC")

# =============================================
# 6. Find shared/specific PMDs
# =============================================
shared_pmd <- intersect(EAC_pmd, ESCC_pmd)
EAC_specific <- setdiff(EAC_pmd, ESCC_pmd)
ESCC_specific <- setdiff(ESCC_pmd, EAC_pmd)

# Save results
savePMDSegments(shared_pmd,
                file.path(output_dir, "Shared_PMDs.rds"),
                file.path(output_dir, "Shared_PMDs.tab"))

savePMDSegments(EAC_specific,
                file.path(output_dir, "EAC_specific_PMDs.rds"),
                file.path(output_dir, "EAC_specific_PMDs.tab"))

savePMDSegments(ESCC_specific,
                file.path(output_dir, "ESCC_specific_PMDs.rds"),
                file.path(output_dir, "ESCC_specific_PMDs.tab"))

cat("Analysis completed!\nResults directory:", output_dir, "\n")