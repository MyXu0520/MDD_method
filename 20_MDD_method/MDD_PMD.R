# Install PMDfinder
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("PMDfinder")

# Use PMDfinder
library(PMDfinder)


if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("MethylSeekR")

library(MethylSeekR)
library(GenomicRanges)
library(rtracklayer)


if (!requireNamespace("BiocManager", quietly = TRUE)){install.packages("BiocManager")}
BiocManager::install("MethylSeekR")
BiocManager::install("valr")
if(!requireNamespace("devtools", quietly = TRUE)) {install.packages("devtools")}
devtools::install_github("ran-ran/zoocat")
devtools::install("MMSeekR.data")
devtools::install("MMSeekR") 
library("MMSeekR")
?runMultiModel


# Load required packages
library(MethylSeekR)
library(GenomicRanges)
library(rtracklayer)
library(dplyr)
library(ggplot2)

# Set working directory
setwd("D:/DWB/bedfile")

# 1. Read BED files and prepare data
read_bed_files <- function() {
  # Get all BED files
  bed_files <- list.files(pattern = "\\.bed$")
  
  # Classify samples based on filenames
  eac_tumor <- bed_files[grepl("^EAC|^GEJ", bed_files) & !grepl("Nonmalignant", bed_files)]
  eac_normal <- bed_files[grepl("GEJ_Nonmalignant", bed_files)]
  escc_tumor <- bed_files[grepl("^ESCC", bed_files) & !grepl("Nonmalignant", bed_files)]
  escc_normal <- bed_files[grepl("ESCC_Nonmalignant", bed_files)]
  
  cat("Found samples:\n")
  cat("EAC tumor:", length(eac_tumor), "samples:", paste(eac_tumor, collapse = ", "), "\n")
  cat("EAC normal:", length(eac_normal), "samples:", paste(eac_normal, collapse = ", "), "\n")
  cat("ESCC tumor:", length(escc_tumor), "samples:", paste(escc_tumor, collapse = ", "), "\n")
  cat("ESCC normal:", length(escc_normal), "samples:", paste(escc_normal, collapse = ", "), "\n")
  
  return(list(
    eac_tumor = eac_tumor,
    eac_normal = eac_normal,
    escc_tumor = escc_tumor,
    escc_normal = escc_normal
  ))
}

# 2. Convert BED files to MethylSeekR format
bed_to_methylseekr_format <- function(bed_file) {
  cat("Processing file:", bed_file, "\n")
  
  # Read BED file
  bed_data <- read.table(bed_file, header = FALSE, sep = "\t", 
                         col.names = c("chr", "pos", "meth_rate", "methylated", "unmethylated"))
  
  # Filter low coverage sites (minimum 5 reads)
  bed_data <- bed_data[bed_data$methylated + bed_data$unmethylated >= 5, ]
  
  # Create GRanges object
  gr <- GRanges(
    seqnames = bed_data$chr,
    ranges = IRanges(start = bed_data$pos, end = bed_data$pos),
    methylated = bed_data$methylated,
    total = bed_data$methylated + bed_data$unmethylated
  )
  
  return(gr)
}

# 3. Use MethylSeekR to identify PMDs
identify_pmds_with_methylseekr <- function(meth_gr, sample_name, output_dir = "methylseekr_results") {
  dir.create(output_dir, showWarnings = FALSE)
  
  cat("Identifying PMDs for sample", sample_name, "...\n")
  
  tryCatch({
    # Calculate methylation level
    meth_gr$meth <- meth_gr$methylated / meth_gr$total
    
    # Set hg38 chromosome lengths
    seqLengths <- c(
      chr1 = 248956422, chr2 = 242193529, chr3 = 198295559, chr4 = 190214555,
      chr5 = 181538259, chr6 = 170805979, chr7 = 159345973, chr8 = 145138636,
      chr9 = 138394717, chr10 = 133797422, chr11 = 135086622, chr12 = 133275309,
      chr13 = 114364328, chr14 = 107043718, chr15 = 101991189, chr16 = 90338345,
      chr17 = 83257441, chr18 = 80373285, chr19 = 58617616, chr20 = 64444167,
      chr21 = 46709983, chr22 = 50818468, chrX = 156040895, chrY = 57227415
    )
    
    # Keep only chromosomes with length information
    available_chrs <- intersect(names(seqLengths), as.character(unique(seqnames(meth_gr))))
    meth_gr <- meth_gr[as.character(seqnames(meth_gr)) %in% available_chrs]
    seqLengths <- seqLengths[available_chrs]
    
    # Set MethylSeekR parameters
    FDR.cutoff <- 5  # FDR cutoff for PMD segmentation
    m.cutoff <- 0.5  # methylation cutoff for PMD boundaries
    
    # Identify PMDs
    pmds <- segmentPMDs(m = meth_gr, 
                        seqLengths = seqLengths,
                        num.cores = 1,  # Use single core to avoid memory issues
                        FDR.cutoff = FDR.cutoff,
                        m.cutoff = m.cutoff)
    
    # Save results
    pmds_file <- file.path(output_dir, paste0(sample_name, "_pmds.bed"))
    export.bed(pmds, pmds_file)
    
    cat("Successfully identified", length(pmds), "PMDs for", sample_name, "\n")
    
    return(pmds)
    
  }, error = function(e) {
    cat("Error processing sample", sample_name, ":", e$message, "\n")
    return(GRanges())  # Return empty GRanges object
  })
}

# 4. Process all samples
process_all_samples <- function(sample_files) {
  all_pmds <- list()
  
  # Process EAC tumor samples
  for (file in sample_files$eac_tumor) {
    sample_name <- gsub("\\.bed$", "", file)
    meth_gr <- bed_to_methylseekr_format(file)
    pmds <- identify_pmds_with_methylseekr(meth_gr, sample_name)
    all_pmds[[sample_name]] <- list(pmds = pmds, type = "EAC_tumor")
  }
  
  # Process EAC normal samples
  for (file in sample_files$eac_normal) {
    sample_name <- gsub("\\.bed$", "", file)
    meth_gr <- bed_to_methylseekr_format(file)
    pmds <- identify_pmds_with_methylseekr(meth_gr, sample_name)
    all_pmds[[sample_name]] <- list(pmds = pmds, type = "EAC_normal")
  }
  
  # Process ESCC tumor samples
  for (file in sample_files$escc_tumor) {
    sample_name <- gsub("\\.bed$", "", file)
    meth_gr <- bed_to_methylseekr_format(file)
    pmds <- identify_pmds_with_methylseekr(meth_gr, sample_name)
    all_pmds[[sample_name]] <- list(pmds = pmds, type = "ESCC_tumor")
  }
  
  # Process ESCC normal samples
  for (file in sample_files$escc_normal) {
    sample_name <- gsub("\\.bed$", "", file)
    meth_gr <- bed_to_methylseekr_format(file)
    pmds <- identify_pmds_with_methylseekr(meth_gr, sample_name)
    all_pmds[[sample_name]] <- list(pmds = pmds, type = "ESCC_normal")
  }
  
  return(all_pmds)
}

# 5. Identify shared and specific PMDs
identify_shared_specific_pmds <- function(all_pmds) {
  # Extract PMDs for each group
  eac_tumor_pmds <- do.call(c, lapply(all_pmds[grep("EAC_tumor", sapply(all_pmds, function(x) x$type))], function(x) x$pmds))
  escc_tumor_pmds <- do.call(c, lapply(all_pmds[grep("ESCC_tumor", sapply(all_pmds, function(x) x$type))], function(x) x$pmds))
  
  # Ensure no duplicate regions
  eac_tumor_pmds <- reduce(eac_tumor_pmds)
  escc_tumor_pmds <- reduce(escc_tumor_pmds)
  
  # Calculate overlaps
  eac_only <- setdiff(eac_tumor_pmds, escc_tumor_pmds)
  escc_only <- setdiff(escc_tumor_pmds, eac_tumor_pmds)
  shared <- intersect(eac_tumor_pmds, escc_tumor_pmds)
  
  cat("Identification results:\n")
  cat("EAC-specific PMDs:", length(eac_only), "\n")
  cat("ESCC-specific PMDs:", length(escc_only), "\n")
  cat("Shared PMDs:", length(shared), "\n")
  
  return(list(
    eac_specific = eac_only,
    escc_specific = escc_only,
    shared = shared
  ))
}

# 6. Save results
save_pmd_results <- function(pmd_results, output_dir = "pmd_results") {
  dir.create(output_dir, showWarnings = FALSE)
  
  # Save EAC-specific PMDs
  if (length(pmd_results$eac_specific) > 0) {
    export.bed(pmd_results$eac_specific, file.path(output_dir, "EAC_specific_PMDs.bed"))
  }
  
  # Save ESCC-specific PMDs
  if (length(pmd_results$escc_specific) > 0) {
    export.bed(pmd_results$escc_specific, file.path(output_dir, "ESCC_specific_PMDs.bed"))
  }
  
  # Save shared PMDs
  if (length(pmd_results$shared) > 0) {
    export.bed(pmd_results$shared, file.path(output_dir, "Shared_PMDs.bed"))
  }
  
  # Create statistical report
  stats <- data.frame(
    PMD_Type = c("EAC_Specific", "ESCC_Specific", "Shared"),
    Count = c(length(pmd_results$eac_specific), 
              length(pmd_results$escc_specific), 
              length(pmd_results$shared)),
    Total_Size_bp = c(sum(width(pmd_results$eac_specific)),
                      sum(width(pmd_results$escc_specific)),
                      sum(width(pmd_results$shared)))
  )
  
  write.csv(stats, file.path(output_dir, "PMD_statistics.csv"), row.names = FALSE)
  
  # Create visualization
  p <- ggplot(stats, aes(x = PMD_Type, y = Count, fill = PMD_Type)) +
    geom_bar(stat = "identity") +
    theme_minimal() +
    labs(title = "PMD Distribution by Type", x = "PMD Type", y = "Count")
  
  ggsave(file.path(output_dir, "PMD_distribution.png"), p, width = 8, height = 6)
  
  cat("Results saved to", output_dir, "directory\n")
}

# 7. Main function
main <- function() {
  cat("Starting PMD analysis...\n")
  
  # Read files
  sample_files <- read_bed_files()
  
  # Process all samples
  all_pmds <- process_all_samples(sample_files)
  
  # Identify shared and specific PMDs
  pmd_results <- identify_shared_specific_pmds(all_pmds)
  
  # Save results
  save_pmd_results(pmd_results)
  
  cat("Analysis completed!\n")
  return(pmd_results)
}

# Run main function
results <- main()

# If memory issues occur, process large samples in batches
# Alternative batch processing:

process_large_samples <- function(sample_files, batch_size = 3) {
  all_pmds <- list()
  
  # Process EAC tumor samples in batches
  eac_batches <- split(sample_files$eac_tumor, 
                       ceiling(seq_along(sample_files$eac_tumor)/batch_size))
  
  for (i in seq_along(eac_batches)) {
    cat("Processing EAC tumor batch", i, "/", length(eac_batches), "\n")
    for (file in eac_batches[[i]]) {
      sample_name <- gsub("\\.bed$", "", file)
      meth_gr <- bed_to_methylseekr_format(file)
      pmds <- identify_pmds_with_methylseekr(meth_gr, sample_name)
      all_pmds[[sample_name]] <- list(pmds = pmds, type = "EAC_tumor")
    }
    gc()  # Garbage collection to free memory
  }
  
  # Similarly process other groups...
  # Add ESCC tumor, normal sample processing here
  
  return(all_pmds)
}

# If main function fails due to memory issues, try:
# results <- process_large_samples(read_bed_files())
# pmd_results <- identify_shared_specific_pmds(results)
# save_pmd_results(pmd_results)