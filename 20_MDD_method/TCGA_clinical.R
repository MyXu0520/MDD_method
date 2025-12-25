# =============================================================================
# Clinical Variable Impact Analysis Scheme
# =============================================================================

# Load necessary packages
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(broom)
library(rstatix)
library(patchwork)

# Set seed for reproducibility
set.seed(123)

# =============================================================================
# 1. Data Loading and ID Conversion
# =============================================================================

cat("Starting to load real data...\n")

# Load MDD probe list
load("D:/PMD_result/4/cg_TCGA_select.RData")
cat("MDD probe list loaded\n")

# Load methylation data
TCGA_ESCA <- read.table("E:/33cancer/33cancers/knn_data/ESCA.methylation450.knn.txt", 
                        header = TRUE, row.names = 1)
cat("Methylation data loaded, dimensions:", dim(TCGA_ESCA), "\n")

# Load clinical data
clinical_EAC <- read_tsv("D:/PMD_result/4/EAC/clinical.cart.2024-11-20/clinical.tsv", 
                         show_col_types = FALSE)
clinical_ESCC <- read_tsv("D:/PMD_result/4/ESCC/clinical.cart.2024-11-20/clinical.tsv", 
                          show_col_types = FALSE)

cat("Clinical data loaded:\n")
cat("EAC clinical samples:", nrow(clinical_EAC), "\n")
cat("ESCC clinical samples:", nrow(clinical_ESCC), "\n")

# =============================================================================
# 2. ID Format Conversion Functions
# =============================================================================

# Function: Convert methylation sample IDs to TCGA standard format
convert_methyl_id_to_tcga <- function(methyl_ids) {
  # Input format: "TCGA.LN.A8I0.01A"
  # Output format: "TCGA-LN-A8I0"
  tcga_ids <- gsub("\\.", "-", methyl_ids)  # Replace dots with hyphens
  tcga_ids <- substr(tcga_ids, 1, 12)  # Take first 12 characters
  return(tcga_ids)
}

# =============================================================================
# 3. Calculate MDD Methylation Features
# =============================================================================

calculate_MDD_methylation <- function(methyl_data, cpg_list) {
  MDD_cpgs <- unique(unlist(cpg_list))
  common_cpgs <- intersect(rownames(methyl_data), MDD_cpgs)
  
  if (length(common_cpgs) > 0) {
    MDD_avg <- colMeans(methyl_data[common_cpgs, ], na.rm = TRUE)
    return(MDD_avg)
  }
  return(NULL)
}

# Filter primary tumor samples (ending with 01A)
primary_samples <- grep("01A$", colnames(TCGA_ESCA), value = TRUE)
cat("Primary tumor samples:", length(primary_samples), "\n")

# Calculate average methylation for three MDD regions
MDD_features <- data.frame(
  sample_id = primary_samples,
  patient_id = convert_methyl_id_to_tcga(primary_samples),
  shared_MDD = calculate_MDD_methylation(TCGA_ESCA[, primary_samples], shared_list),
  EAC_MDD = calculate_MDD_methylation(TCGA_ESCA[, primary_samples], EAC_spec_list),
  ESCC_MDD = calculate_MDD_methylation(TCGA_ESCA[, primary_samples], ESCC_spec_list)
)

cat("MDD features calculated, samples:", nrow(MDD_features), "\n")
cat("Example patient IDs in MDD data:", head(MDD_features$patient_id, 5), "\n")

# =============================================================================
# 4. Process Clinical Data and Merge (Fix Subtype Assignment)
# =============================================================================

# Add subtype labels to clinical data (before merging)
clinical_EAC$subtype <- "EAC"
clinical_ESCC$subtype <- "ESCC"

# Merge clinical data
clinical_all <- bind_rows(clinical_EAC, clinical_ESCC)

cat("\nMerged clinical data rows:", nrow(clinical_all), "\n")
cat("Subtype distribution in merged clinical data:\n")
print(table(clinical_all$subtype))

# Clean clinical data (keep one row per patient with subtype information)
clean_clinical_data <- function(clinical_df) {
  # First determine subtype for each patient (take first non-missing value)
  subtype_info <- clinical_df %>%
    group_by(case_submitter_id) %>%
    summarise(
      subtype = first(na.omit(subtype)),
      .groups = 'drop'
    )
  
  # Then deduplicate to get clinical information
  clinical_clean <- clinical_df %>%
    distinct(case_submitter_id, .keep_all = TRUE) %>%
    select(
      patient_id = case_submitter_id,
      age_at_index,
      gender,
      race,
      ethnicity,
      vital_status,
      days_to_death,
      days_to_last_follow_up,
      ajcc_pathologic_stage,
      primary_diagnosis
    ) %>%
    mutate(
      # Age
      age = as.numeric(age_at_index),
      age_group = case_when(
        age < 50 ~ "<50",
        age >= 50 & age < 60 ~ "50-59",
        age >= 60 & age < 70 ~ "60-69",
        age >= 70 ~ "≥70",
        TRUE ~ "Unknown"
      ),
      
      # Standardize gender
      gender = case_when(
        tolower(gender) == "male" ~ "Male",
        tolower(gender) == "female" ~ "Female",
        TRUE ~ NA_character_
      ),
      
      # Simplify tumor stage (adjust based on your actual data)
      stage = case_when(
        grepl("Stage I[A]?$", ajcc_pathologic_stage, ignore.case = TRUE) ~ "I",
        grepl("Stage II[A|B]?$", ajcc_pathologic_stage, ignore.case = TRUE) ~ "II",
        grepl("Stage III[A|B|C]?$", ajcc_pathologic_stage, ignore.case = TRUE) ~ "III",
        grepl("Stage IV[A]?$", ajcc_pathologic_stage, ignore.case = TRUE) ~ "IV",
        ajcc_pathologic_stage == "'--" ~ "Unknown",
        TRUE ~ "Unknown"
      )
    )
  
  # Merge subtype information
  clinical_clean <- clinical_clean %>%
    left_join(subtype_info, by = c("patient_id" = "case_submitter_id"))
  
  return(clinical_clean)
}

# Clean clinical data
clinical_clean <- clean_clinical_data(clinical_all)

cat("\nClinical data cleaned:\n")
cat("Patients:", nrow(clinical_clean), "\n")
cat("EAC patients:", sum(clinical_clean$subtype == "EAC", na.rm = TRUE), "\n")
cat("ESCC patients:", sum(clinical_clean$subtype == "ESCC", na.rm = TRUE), "\n")

# Merge MDD data and clinical data
combined_data <- MDD_features %>%
  left_join(clinical_clean, by = c("patient_id" = "patient_id")) %>%
  filter(
    !is.na(shared_MDD) & !is.na(EAC_MDD) & !is.na(ESCC_MDD) &
      !is.na(age) & !is.na(gender) & !is.na(stage) & !is.na(subtype)
  )

cat("\nData merged, total samples:", nrow(combined_data), "\n")
cat("EAC samples after merging:", sum(combined_data$subtype == "EAC", na.rm = TRUE), "\n")
cat("ESCC samples after merging:", sum(combined_data$subtype == "ESCC", na.rm = TRUE), "\n")
cat("Age range:", range(combined_data$age, na.rm = TRUE), "\n")
cat("Gender distribution: Male =", sum(combined_data$gender == "Male", na.rm = TRUE), 
    ", Female =", sum(combined_data$gender == "Female", na.rm = TRUE), "\n")
cat("Stage distribution:\n")
print(table(combined_data$stage))

# =============================================================================
# 5. Analysis Functions (Optimized for Compact A4 Layout)
# =============================================================================

# Create unified theme function
create_compact_theme <- function(base_size = 9) {
  theme_bw(base_size = base_size) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = base_size + 1, margin = margin(b = 3)),
      axis.title = element_text(size = base_size),
      axis.text = element_text(size = base_size - 1),
      strip.text = element_text(size = base_size, face = "bold", margin = margin(2, 0, 2, 0)),
      legend.position = "none",
      panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      plot.margin = margin(3, 3, 3, 3),
      panel.spacing = unit(0.5, "lines"),
      strip.background = element_rect(fill = "grey95")
    )
}

# 5.1 Age Effect Analysis
analyze_age_effect <- function(data) {
  analysis_data <- data %>% filter(!is.na(age))
  
  if (nrow(analysis_data) < 10) return(NULL)
  
  plot_data <- analysis_data %>%
    select(age, shared_MDD, EAC_MDD, ESCC_MDD) %>%
    pivot_longer(cols = -age, names_to = "MDD_type", values_to = "methylation") %>%
    mutate(
      MDD_type = factor(MDD_type,
                        levels = c("shared_MDD", "EAC_MDD", "ESCC_MDD"),
                        labels = c("Shared MDD", "EAC-specific MDD", "ESCC-specific MDD"))
    )
  
  # Calculate statistics
  stat_data <- plot_data %>%
    group_by(MDD_type) %>%
    summarise(
      cor = cor(age, methylation, use = "complete.obs"),
      p_value = cor.test(age, methylation)$p.value,
      .groups = 'drop'
    ) %>%
    mutate(
      label = ifelse(p_value < 0.001, "***",
                     ifelse(p_value < 0.01, "**",
                            ifelse(p_value < 0.05, "*", "ns"))),
      x_pos = min(plot_data$age) + 0.7 * diff(range(plot_data$age)),
      y_pos = max(plot_data$methylation) - 0.05 * diff(range(plot_data$methylation))
    )
  
  p <- ggplot(plot_data, aes(x = age, y = methylation)) +
    geom_point(alpha = 0.5, color = "steelblue", size = 1.5) +
    geom_smooth(method = "lm", formula = y ~ x, color = "red", se = FALSE, linewidth = 0.8) +
    facet_wrap(~ MDD_type, scales = "free_y", ncol = 3) +
    geom_text(data = stat_data, 
              aes(x = x_pos, y = y_pos, 
                  label = sprintf("r = %.3f\n%s", cor, label)),
              size = 2.8, hjust = 0, inherit.aes = FALSE) +
    labs(x = "Age (years)", y = "Methylation Level",
         title = "A. Age Effects on MDD Methylation Levels") +
    create_compact_theme(base_size = 9)
  
  return(p)
}

# 5.2 Gender Effect Analysis
analyze_gender_effect <- function(data) {
  gender_data <- data %>% filter(gender %in% c("Male", "Female"))
  if (nrow(gender_data) < 10) return(NULL)
  
  plot_data <- gender_data %>%
    select(gender, shared_MDD, EAC_MDD, ESCC_MDD) %>%
    pivot_longer(cols = -gender, names_to = "MDD_type", values_to = "methylation") %>%
    mutate(
      MDD_type = factor(MDD_type,
                        levels = c("shared_MDD", "EAC_MDD", "ESCC_MDD"),
                        labels = c("Shared MDD", "EAC-specific MDD", "ESCC-specific MDD")),
      gender = factor(gender, levels = c("Male", "Female"))
    )
  
  p <- ggplot(plot_data, aes(x = gender, y = methylation, fill = gender)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA, width = 0.6, linewidth = 0.4) +
    geom_jitter(width = 0.15, alpha = 0.5, size = 1) +
    facet_wrap(~ MDD_type, scales = "free_y", ncol = 3) +
    stat_compare_means(
      method = "t.test",
      comparisons = list(c("Male", "Female")),
      label = "p.signif",
      symnum.args = list(
        cutpoints = c(0, 0.001, 0.01, 0.05, 1),
        symbols = c("***", "**", "*", "ns")
      ),
      vjust = 0.5,
      size = 3.5,
      bracket.size = 0.3,
      tip.length = 0.02
    ) +
    scale_fill_manual(values = c("Male" = "#1f77b4", "Female" = "#ff7f0e")) +
    labs(x = "Gender", y = "Methylation Level",
         title = "B. Gender Effects on MDD Methylation Levels") +
    create_compact_theme(base_size = 9)
  
  return(p)
}

# 5.3 Tumor Stage Effect Analysis (Complete Significance Annotation)
analyze_stage_effect <- function(data) {
  stage_levels <- c("I", "II", "III", "IV")
  stage_data <- data %>% 
    filter(stage %in% stage_levels) %>%
    mutate(stage = factor(stage, levels = stage_levels))
  
  if (nrow(stage_data) < 10) return(NULL)
  
  plot_data <- stage_data %>%
    select(stage, shared_MDD, EAC_MDD, ESCC_MDD) %>%
    pivot_longer(cols = -stage, names_to = "MDD_type", values_to = "methylation") %>%
    mutate(
      MDD_type = factor(MDD_type,
                        levels = c("shared_MDD", "EAC_MDD", "ESCC_MDD"),
                        labels = c("Shared MDD", "EAC-specific MDD", "ESCC-specific MDD"))
    )
  
  # Perform ANOVA and pairwise comparisons for each MDD type
  plots_list <- list()
  
  for (mdd_label in c("Shared MDD", "EAC-specific MDD", "ESCC-specific MDD")) {
    mdd_data <- plot_data %>% filter(MDD_type == mdd_label)
    
    # Check sample count per stage
    stage_counts <- table(mdd_data$stage)
    valid_stages <- names(stage_counts[stage_counts >= 2])
    
    # Create base plot
    p_single <- ggplot(mdd_data, aes(x = stage, y = methylation, fill = stage)) +
      geom_boxplot(alpha = 0.7, outlier.shape = NA, width = 0.6, linewidth = 0.4) +
      geom_jitter(width = 0.15, alpha = 0.5, size = 1) +
      scale_fill_brewer(palette = "Purples") +
      labs(x = "Tumor Stage", y = "Methylation Level") +
      theme_bw(base_size = 9) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 10),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 7),
        legend.position = "none",
        panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
        panel.grid.minor = element_blank(),
        plot.margin = margin(3, 3, 3, 3),
        panel.spacing = unit(0.5, "lines")
      )
    
    # If at least two groups with at least 2 samples each, perform pairwise comparisons
    if (length(valid_stages) >= 2) {
      # Prepare pairwise comparisons
      comparisons <- combn(valid_stages, 2, simplify = FALSE)
      
      # Perform t-tests with multiple testing correction
      tryCatch({
        stat_test <- mdd_data %>%
          filter(stage %in% valid_stages) %>%
          group_by(MDD_type) %>%
          t_test(methylation ~ stage) %>%
          adjust_pvalue(method = "BH") %>%
          add_significance("p.adj") %>%
          add_xy_position(x = "stage", step.increase = 0.1)
        
        if (nrow(stat_test) > 0) {
          p_single <- p_single +
            stat_pvalue_manual(
              stat_test,
              label = "p.adj.signif",
              tip.length = 0.01,
              size = 3,
              bracket.size = 0.3
            )
        }
      }, error = function(e) {
        cat(paste0("  ", mdd_label, " pairwise comparison error: ", e$message, "\n"))
      })
    }
    
    plots_list[[mdd_label]] <- p_single
  }
  
  # Combine three subplots into one figure
  if (length(plots_list) == 3) {
    # Use patchwork to combine three plots
    combined_plot <- plots_list[[1]] + plots_list[[2]] + plots_list[[3]] +
      plot_layout(ncol = 3) +
      plot_annotation(title = "C. Tumor Stage Effects on MDD Methylation Levels",
                      theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 10)))
  } else {
    # Fallback: Use facet_wrap
    combined_plot <- ggplot(plot_data, aes(x = stage, y = methylation, fill = stage)) +
      geom_boxplot(alpha = 0.7, outlier.shape = NA) +
      geom_jitter(width = 0.2, alpha = 0.5, size = 1) +
      facet_wrap(~ MDD_type, scales = "free_y", ncol = 3) +
      scale_fill_brewer(palette = "Purples") +
      labs(x = "Tumor Stage", y = "Methylation Level",
           title = "C. Tumor Stage Effects on MDD Methylation Levels") +
      create_compact_theme(base_size = 9)
  }
  
  return(combined_plot)
}

# 5.4 Subtype Effect Analysis
analyze_subtype_effect <- function(data) {
  subtype_data <- data %>% 
    filter(subtype %in% c("EAC", "ESCC")) %>%
    mutate(subtype = factor(subtype, levels = c("EAC", "ESCC")))
  
  if (nrow(subtype_data) < 10) return(NULL)
  
  plot_data <- subtype_data %>%
    select(subtype, shared_MDD, EAC_MDD, ESCC_MDD) %>%
    pivot_longer(cols = -subtype, names_to = "MDD_type", values_to = "methylation") %>%
    mutate(
      MDD_type = factor(MDD_type,
                        levels = c("shared_MDD", "EAC_MDD", "ESCC_MDD"),
                        labels = c("Shared MDD", "EAC-specific MDD", "ESCC-specific MDD"))
    )
  
  p <- ggplot(plot_data, aes(x = subtype, y = methylation, fill = subtype)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA, width = 0.6, linewidth = 0.4) +
    geom_jitter(width = 0.15, alpha = 0.5, size = 1) +
    facet_wrap(~ MDD_type, scales = "free_y", ncol = 3) +
    stat_compare_means(
      method = "t.test",
      comparisons = list(c("EAC", "ESCC")),
      label = "p.signif",
      symnum.args = list(
        cutpoints = c(0, 0.001, 0.01, 0.05, 1),
        symbols = c("***", "**", "*", "ns")
      ),
      vjust = 0.5,
      size = 3.5,
      bracket.size = 0.3,
      tip.length = 0.02
    ) +
    scale_fill_manual(values = c("EAC" = "#F5C748", "ESCC" = "#07888A")) +
    labs(x = "Cancer Subtype", y = "Methylation Level",
         title = "D. Subtype Differences in MDD Methylation Levels") +
    create_compact_theme(base_size = 9)
  
  return(p)
}

# =============================================================================
# 6. Execute All Analyses and Create Compact A4 Vertical Layout
# =============================================================================

# Create output directory
output_dir <- "D:/PMD_result/revised_result"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

cat("\nStarting all analyses...\n")

# Execute analyses
all_plots <- list()

age_plot <- analyze_age_effect(combined_data)
if (!is.null(age_plot)) all_plots$age <- age_plot

gender_plot <- analyze_gender_effect(combined_data)
if (!is.null(gender_plot)) all_plots$gender <- gender_plot

stage_plot <- analyze_stage_effect(combined_data)
if (!is.null(stage_plot)) all_plots$stage <- stage_plot

subtype_plot <- analyze_subtype_effect(combined_data)
if (!is.null(subtype_plot)) all_plots$subtype <- subtype_plot

# Check which plots were successfully generated
cat("\nSuccessfully generated plots:\n")
for (plot_name in c("age", "gender", "stage", "subtype")) {
  if (!is.null(all_plots[[plot_name]])) {
    cat(paste0(plot_name, ": Successful\n"))
  } else {
    cat(paste0(plot_name, ": Failed (insufficient data)\n"))
  }
}

# =============================================================================
# 7. Create Compact A4 Vertical Layout (No Main Title)
# =============================================================================

cat("\nCreating compact A4 vertical layout...\n")

if (length(all_plots) > 0) {
  # A4 dimensions: 210×297mm or 8.27×11.69 inches
  a4_width <- 8.27  # inches
  a4_height <- 11.69  # inches
  
  # Adjust height allocation per plot based on number of available plots
  n_plots <- length(all_plots)
  
  # Create combined figure (no main title)
  combined_figure <- NULL
  
  # Add plots in order
  plot_order <- c("age", "gender", "stage", "subtype")
  for (plot_name in plot_order) {
    if (!is.null(all_plots[[plot_name]])) {
      if (is.null(combined_figure)) {
        combined_figure <- all_plots[[plot_name]]
      } else {
        combined_figure <- combined_figure / all_plots[[plot_name]]
      }
    }
  }
  
  # Add small caption (bottom, minimal space)
  combined_figure <- combined_figure +
    plot_annotation(
      caption = paste("Generated on", format(Sys.Date(), "%Y-%m-%d"),
                      "| *** p<0.001, ** p<0.01, * p<0.05, ns not significant"),
      theme = theme(
        plot.caption = element_text(size = 8, hjust = 1, color = "grey50", 
                                    margin = margin(t = 5))
      )
    ) &
    theme(
      plot.margin = margin(5, 5, 5, 5)
    )
  
  # Calculate appropriate height (compact layout)
  if (n_plots == 1) {
    final_height <- 4
  } else if (n_plots == 2) {
    final_height <- 6
  } else if (n_plots == 3) {
    final_height <- 8.5
  } else {
    # 4 plots: ~2.5 inches each, plus margins
    final_height <- 10.5
  }
  
  # Save compact A4 vertical layout
  tryCatch({
    ggsave(
      file.path(output_dir, "A4_compact_clinical_effects_on_MDD.pdf"),
      combined_figure,
      width = a4_width,
      height = final_height,
      device = "pdf"
    )
    cat("Compact A4 vertical layout saved:", file.path(output_dir, "A4_compact_clinical_effects_on_MDD.pdf"), "\n")
    cat("Dimensions: ", a4_width, "×", final_height, "inches\n")
  }, error = function(e) {
    cat("Error saving PDF:", e$message, "\n")
    # Try standard A4 height
    ggsave(
      file.path(output_dir, "A4_compact_clinical_effects_on_MDD.pdf"),
      combined_figure,
      width = a4_width,
      height = a4_height,
      device = "pdf"
    )
  })
  
  # Also save PNG version (high resolution)
  tryCatch({
    ggsave(
      file.path(output_dir, "A4_compact_clinical_effects_on_MDD.png"),
      combined_figure,
      width = a4_width,
      height = final_height,
      dpi = 600,  # Higher resolution
      bg = "white"
    )
    cat("PNG version saved (600 dpi)\n")
  }, error = function(e) {
    cat("Error saving PNG:", e$message, "\n")
    # Try 300 dpi
    ggsave(
      file.path(output_dir, "A4_compact_clinical_effects_on_MDD.png"),
      combined_figure,
      width = a4_width,
      height = final_height,
      dpi = 300,
      bg = "white"
    )
  })
  
  # Save each plot individually for inspection
  for (plot_name in names(all_plots)) {
    tryCatch({
      ggsave(
        file.path(output_dir, paste0(plot_name, "_effect_individual.pdf")),
        all_plots[[plot_name]],
        width = 10,
        height = 3.5,  # More compact
        device = "pdf"
      )
    }, error = function(e) {
      cat(paste0("Error saving ", plot_name, " individual plot:", e$message, "\n"))
    })
  }
  
} else {
  cat("No analysis plots available\n")
}

# =============================================================================
# 8. Create Data Summary Plot (Optional)
# =============================================================================

create_data_summary_plot <- function(data) {
  # Create simple data summary
  summary_stats <- data.frame(
    Category = c("Total Samples", "EAC", "ESCC", "Male", "Female", 
                 "Stage I", "Stage II", "Stage III", "Stage IV"),
    Count = c(
      nrow(data),
      sum(data$subtype == "EAC"),
      sum(data$subtype == "ESCC"),
      sum(data$gender == "Male"),
      sum(data$gender == "Female"),
      sum(data$stage == "I"),
      sum(data$stage == "II"),
      sum(data$stage == "III"),
      sum(data$stage == "IV")
    )
  )
  
  # Add colors
  summary_stats <- summary_stats %>%
    mutate(
      Group = case_when(
        Category %in% c("Total Samples") ~ "Overall",
        Category %in% c("EAC", "ESCC") ~ "Subtype",
        Category %in% c("Male", "Female") ~ "Gender",
        Category %in% c("Stage I", "Stage II", "Stage III", "Stage IV") ~ "Stage",
        TRUE ~ "Other"
      ),
      Group = factor(Group, levels = c("Overall", "Subtype", "Gender", "Stage"))
    )
  
  p <- ggplot(summary_stats, aes(x = reorder(Category, -Count), y = Count, fill = Group)) +
    geom_bar(stat = "identity", alpha = 0.8) +
    geom_text(aes(label = Count), vjust = -0.5, size = 2.5) +
    scale_fill_brewer(palette = "Set2") +
    labs(x = "", y = "Count", 
         title = "Dataset Summary",
         subtitle = paste("Total:", nrow(data), "samples")) +
    theme_bw(base_size = 9) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 10),
      plot.subtitle = element_text(hjust = 0.5, size = 9),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = 8),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      plot.margin = margin(5, 5, 5, 20)
    )
  
  return(p)
}

# Generate and save data summary plot
summary_plot <- create_data_summary_plot(combined_data)
tryCatch({
  ggsave(
    file.path(output_dir, "data_summary_plot.pdf"),
    summary_plot,
    width = 6,
    height = 4,
    device = "pdf"
  )
  cat("Data summary plot saved\n")
}, error = function(e) {
  cat("Error saving data summary plot:", e$message, "\n")
})

# =============================================================================
# 9. Save Data Files
# =============================================================================

# Save merged data
tryCatch({
  write.csv(combined_data,
            file.path(output_dir, "combined_MDD_clinical_data.csv"),
            row.names = FALSE)
  cat("Merged data saved\n")
}, error = function(e) {
  cat("Error saving merged data:", e$message, "\n")
})

# Save data summary
data_summary <- data.frame(
  Metric = c("Total Samples", "EAC Samples", "ESCC Samples", 
             "Male", "Female", "Mean Age ± SD", "Age Range",
             "Stage I", "Stage II", "Stage III", "Stage IV"),
  Value = c(
    nrow(combined_data),
    sum(combined_data$subtype == "EAC"),
    sum(combined_data$subtype == "ESCC"),
    sum(combined_data$gender == "Male"),
    sum(combined_data$gender == "Female"),
    paste0(round(mean(combined_data$age, na.rm = TRUE), 1), 
           " ± ", round(sd(combined_data$age, na.rm = TRUE), 1)),
    paste(round(min(combined_data$age, na.rm = TRUE), 1), 
          "-", round(max(combined_data$age, na.rm = TRUE), 1)),
    sum(combined_data$stage == "I"),
    sum(combined_data$stage == "II"),
    sum(combined_data$stage == "III"),
    sum(combined_data$stage == "IV")
  )
)

write.csv(data_summary, 
          file.path(output_dir, "data_summary.csv"),
          row.names = FALSE)

cat("Data summary saved\n")

# =============================================================================
# 10. Generate Brief Run Report
# =============================================================================

sink(file.path(output_dir, "run_report.txt"))
cat("Compact A4 Vertical Layout Generation Report\n")
cat("==============================\n\n")
cat("Generation date:", format(Sys.Date(), "%Y-%m-%d"), "\n")
cat("Generation time:", format(Sys.time(), "%H:%M:%S"), "\n\n")
cat("Data overview:\n")
cat("Total samples:", nrow(combined_data), "\n")
cat("EAC samples:", sum(combined_data$subtype == "EAC"), "\n")
cat("ESCC samples:", sum(combined_data$subtype == "ESCC"), "\n\n")
cat("Generated plots:\n")
for (plot_name in c("age", "gender", "stage", "subtype")) {
  if (!is.null(all_plots[[plot_name]])) {
    cat(paste0("- ", plot_name, " effect plot: Successful\n"))
  } else {
    cat(paste0("- ", plot_name, " effect plot: Failed\n"))
  }
}
cat("\nOutput files:\n")
cat("1. A4_compact_clinical_effects_on_MDD.pdf: Compact A4 vertical layout\n")
cat("2. A4_compact_clinical_effects_on_MDD.png: PNG version (600 dpi)\n")
cat("3. data_summary_plot.pdf: Data summary plot\n")
cat("4. combined_MDD_clinical_data.csv: Merged data\n")
cat("5. data_summary.csv: Data summary\n")
cat("6. run_report.txt: This report\n")
sink()

cat("Run report generated\n")

# =============================================================================
# 11. Completion
# =============================================================================

cat(paste0("\n", strrep("=", 60), "\n"))
cat("Analysis completed!\n")
cat(strrep("=", 60), "\n\n")

cat("Main output file:\n")
cat("Compact A4 vertical layout: ", file.path(output_dir, "A4_compact_clinical_effects_on_MDD.pdf"), "\n")
cat("Dimensions: ", a4_width, "×", final_height, "inches (compact layout)\n\n")

cat("Plot features:\n")
cat("1. No main title, compact A4 layout\n")
cat("2. All plots vertically arranged\n")
cat("3. Complete pairwise significance annotations in plot C\n")
cat("4. Unified compact theme and style\n")
cat("5. Brief caption at bottom\n\n")

cat("Note: If significance annotations in plot C are incomplete, possible reasons:\n")
cat("1. Insufficient samples in some stages\n")
cat("2. Non-significant differences between groups\n")
cat("3. Statistical tests cannot be performed\n")
cat("Check console output for details\n")