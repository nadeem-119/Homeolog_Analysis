# Load required libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(viridis)
library(ggrepel)

# TPM mapping
cat("TPM1 = H, TPM2 = St1, TPM3 = St2 for all tissues.\n")
cat("TPM1, TPM2, TPM3 are pre-averaged values across replicates for each subgenome.\n")

# Data Loading and Preprocessing
Root <- read.delim("Root_homologous_expression.tsv", header = TRUE, sep = "\t")
Stem <- read.delim("Stem_homologous_expression.tsv", header = TRUE, sep = "\t")
Leaf <- read.delim("Leaf_homologous_expression.tsv", header = TRUE, sep = "\t")
Flower <- read.delim("Flower_homologous_expression.tsv", header = TRUE, sep = "\t")

H_gene_list <- Root$chrH
St1_gene_list <- Root$chrSt1
St2_gene_list <- Root$chrSt2

check_gene_consistency <- function(tissue_data, tissue_name) {
  if (!all(tissue_data$chrH == H_gene_list)) warning(sprintf("chrH in %s does not match Root", tissue_name))
  if (!all(tissue_data$chrSt1 == St1_gene_list)) warning(sprintf("chrSt1 in %s does not match Root", tissue_name))
  if (!all(tissue_data$chrSt2 == St2_gene_list)) warning(sprintf("chrSt2 in %s does not match Root", tissue_name))
}

check_gene_consistency(Stem, "Stem")
check_gene_consistency(Leaf, "Leaf")
check_gene_consistency(Flower, "Flower")

prepare_expression <- function(data, tissue_name) {
  required_cols <- c("TPM1", "TPM2", "TPM3")
  if (!all(required_cols %in% colnames(data))) {
    stop(sprintf("Missing required TPM columns for %s. Found columns: %s", 
                 tissue_name, paste(colnames(data), collapse = ", ")))
  }
  data.frame(H_expr = data$TPM1, St1_expr = data$TPM2, St2_expr = data$TPM3)
}

expression_data <- list(
  Root = prepare_expression(Root, "Root"),
  Stem = prepare_expression(Stem, "Stem"),
  Leaf = prepare_expression(Leaf, "Leaf"),
  Flower = prepare_expression(Flower, "Flower")
)

valid_rows <- Reduce("&", lapply(expression_data, function(df) complete.cases(df$H_expr, df$St1_expr, df$St2_expr)))
expression_data <- lapply(expression_data, function(df) df[valid_rows, ])
H_gene_list <- H_gene_list[valid_rows]
St1_gene_list <- St1_gene_list[valid_rows]
St2_gene_list <- St2_gene_list[valid_rows]

# Analysis Parameters
pseudocount <- min(unlist(lapply(expression_data, function(df) unlist(df[df > 0])))) / 10

# Helper Functions
calculate_homoeolog_ratios <- function(expression_data, pseudocount) {
  ratios <- list()
  for (tissue in names(expression_data)) {
    df <- expression_data[[tissue]]
    h_expr <- pmax(df$H_expr, pseudocount)
    st1_expr <- pmax(df$St1_expr, pseudocount)
    st2_expr <- pmax(df$St2_expr, pseudocount)
    ratios[[tissue]] <- data.frame(
      GeneID = H_gene_list,
      H_vs_St1 = log2(h_expr / st1_expr),
      H_vs_St2 = log2(h_expr / st2_expr),
      St1_vs_St2 = log2(st1_expr / st2_expr),
      H_expr = df$H_expr,
      St1_expr = df$St1_expr,
      St2_expr = df$St2_expr,
      row.names = NULL
    )
  }
  return(ratios)
}

# Dominance Classification
classify_homoeolog_dominance <- function(ratios, pseudocount, tolerance = 0.05) {
  dominance_list <- list()
  
  # Define ideal proportions for each category
  ideal_proportions <- list(
    Balanced = c(H = 0.33, St1 = 0.33, St2 = 0.33),
    H_dominant = c(H = 1, St1 = 0, St2 = 0),
    St1_dominant = c(H = 0, St1 = 1, St2 = 0),
    St2_dominant = c(H = 0, St1 = 0, St2 = 1),
    H_suppressed = c(H = 0, St1 = 0.5, St2 = 0.5),
    St1_suppressed = c(H = 0.5, St1 = 0, St2 = 0.5),
    St2_suppressed = c(H = 0.5, St1 = 0.5, St2 = 0)
  )
  
  for (tissue in names(ratios)) {
    ratio_df <- ratios[[tissue]]
    dominance <- vector("character", nrow(ratio_df))
    
    for (i in 1:nrow(ratio_df)) {
      h_expr <- pmax(ratio_df$H_expr[i], pseudocount)
      st1_expr <- pmax(ratio_df$St1_expr[i], pseudocount)
      st2_expr <- pmax(ratio_df$St2_expr[i], pseudocount)
      
      # Calculate total expression and proportions
      total_expr <- h_expr + st1_expr + st2_expr
      h_prop <- h_expr / total_expr
      st1_prop <- st1_expr / total_expr
      st2_prop <- st2_expr / total_expr
      observed_props <- c(H = h_prop, St1 = st1_prop, St2 = st2_prop)
      
      # Check each category
      classified <- FALSE
      for (category in names(ideal_proportions)) {
        ideal <- ideal_proportions[[category]]
        if (all(abs(observed_props - ideal) <= tolerance)) {
          dominance[i] <- category
          classified <- TRUE
          break
        }
      }
      
      # If no category matches within tolerance, default to "Balanced"
      if (!classified) {
        dominance[i] <- "Balanced"
      }
    }
    
    ratio_df$Dominance <- dominance
    dominance_list[[tissue]] <- ratio_df
  }
  return(dominance_list)
}

# Combined Histogram
generate_combined_histogram <- function(ratios_data, dominance_classification, subA, subB, pseudocount) {
  plot_data <- bind_rows(lapply(names(ratios_data), function(tissue) {
    tissue_ratios <- ratios_data[[tissue]]
    subA_expr <- tissue_ratios[[paste0(subA, "_expr")]] + pseudocount
    subB_expr <- tissue_ratios[[paste0(subB, "_expr")]] + pseudocount
    log2_ratio <- log2(subA_expr / subB_expr)
    data.frame(
      GeneID = tissue_ratios$GeneID,
      Tissue = tissue,
      Log2_Ratio = log2_ratio,
      SubA_expr = tissue_ratios[[paste0(subA, "_expr")]],
      SubB_expr = tissue_ratios[[paste0(subB, "_expr")]]
    )
  }))
  
  dominance_data <- bind_rows(dominance_classification, .id = "Tissue") %>%
    select(GeneID, Tissue, Dominance)
  
  plot_data <- plot_data %>%
    left_join(dominance_data, by = c("GeneID", "Tissue")) %>%
    filter(!is.na(Log2_Ratio) & !is.infinite(Log2_Ratio) & Log2_Ratio >= -15 & Log2_Ratio <= 15)
  
  plot_data <- plot_data %>%
    mutate(
      Bias = case_when(
        Log2_Ratio > 1 ~ paste0(subA, "-biased"),
        Log2_Ratio < -1 ~ paste0(subB, "-biased"),
        TRUE ~ "Unbiased"
      )
    )
  
  max_counts <- plot_data %>%
    group_by(Tissue) %>%
    summarise(Max_Count = max(hist(Log2_Ratio, breaks = 100, plot = FALSE)$counts))
  
  stats <- plot_data %>%
    group_by(Tissue) %>%
    summarise(
      SD = sd(Log2_Ratio, na.rm = TRUE),
      SubA_biased = sum(Log2_Ratio > 1) / n() * 100,
      SubB_biased = sum(Log2_Ratio < -1) / n() * 100,
      Unbiased = sum(Log2_Ratio >= -1 & Log2_Ratio <= 1) / n() * 100,
      Weighted_Avg_Log2 = sum(Log2_Ratio * (SubA_expr + SubB_expr), na.rm = TRUE) / sum(SubA_expr + SubB_expr, na.rm = TRUE)
    ) %>%
    left_join(max_counts, by = "Tissue") %>%
    mutate(
      Label = sprintf("SD: %.2f\n%s-biased: %.1f%%\n%s-biased: %.1f%%\nUnbiased: %.1f%%",
                      SD, subA, SubA_biased, subB, SubB_biased, Unbiased),
      Y_Position = Max_Count * 0.9
    )
  
  cat(sprintf("\nSummary for Combined Histogram (%s vs %s):\n", subA, subB))
  print(stats)
  
  x_label <- substitute(log[2] ~ (a/b), list(a = subA, b = subB))
  
  p <- ggplot(plot_data, aes(x = Log2_Ratio, fill = Bias)) +
    geom_histogram(bins = 100, aes(color = after_stat(fill)), linewidth = 0.2, alpha = 0.7) +
    scale_fill_manual(values = setNames(
      c("#1f77b4", "#ff7f0e", "#808080"),
      c(paste0(subA, "-biased"), paste0(subB, "-biased"), "Unbiased")
    )) +
    scale_color_manual(values = setNames(
      c("#1f77b4", "#ff7f0e", "#808080"),
      c(paste0(subA, "-biased"), paste0(subB, "-biased"), "Unbiased")
    )) +
    guides(color = "none") +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "red", linewidth = 1) +
    geom_text(data = stats, aes(x = Inf, y = Y_Position, label = Label), 
              hjust = 1, vjust = 1, size = 3.5, color = "black", inherit.aes = FALSE) +
    facet_wrap(~ Tissue, ncol = 2) +
    scale_x_continuous(breaks = seq(-15, 15, by = 5), limits = c(-15, 15)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.3))) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      axis.line = element_line(color = "black"),
      axis.text = element_text(color = "black", size = 10, face = "bold"),
      axis.title = element_text(size = 12, face = "bold"),
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      strip.text = element_text(size = 12, face = "bold"),
      panel.spacing = unit(1, "lines"),
      legend.position = "right",
      legend.title = element_text(size = 10, face = "bold"),
      legend.text = element_text(size = 8)
    ) +
    labs(title = paste(subA, "vs", subB), x = x_label, y = "Count", fill = "Bias")
  
  return(p)
}

# Averaged Histogram
generate_avg_histogram <- function(ratios_data, dominance_classification, subA, subB, pseudocount) {
  avg_data <- data.frame(
    GeneID = H_gene_list,
    Log2_Ratio = rowMeans(sapply(ratios_data, function(tissue) {
      subA_expr <- tissue[[paste0(subA, "_expr")]] + pseudocount
      subB_expr <- tissue[[paste0(subB, "_expr")]] + pseudocount
      log2(subA_expr / subB_expr)
    }), na.rm = TRUE),
    SubA_expr = rowMeans(sapply(ratios_data, function(tissue) {
      tissue[[paste0(subA, "_expr")]]
    }), na.rm = TRUE),
    SubB_expr = rowMeans(sapply(ratios_data, function(tissue) {
      tissue[[paste0(subB, "_expr")]]
    }), na.rm = TRUE)
  )
  
  dominance_across_tissues <- bind_rows(dominance_classification, .id = "Tissue") %>%
    select(GeneID, Tissue, Dominance) %>%
    pivot_wider(names_from = Tissue, values_from = Dominance) %>%
    rowwise() %>%
    mutate(
      Dominance = names(sort(table(c(Root, Stem, Leaf, Flower)), decreasing = TRUE))[1]
    ) %>%
    ungroup() %>%
    select(GeneID, Dominance)
  
  avg_data <- avg_data %>%
    left_join(dominance_across_tissues, by = "GeneID") %>%
    filter(!is.na(Log2_Ratio) & !is.infinite(Log2_Ratio) & Log2_Ratio >= -15 & Log2_Ratio <= 15)
  
  avg_data <- avg_data %>%
    mutate(
      Bias = case_when(
        Log2_Ratio > 1 ~ paste0(subA, "-biased"),
        Log2_Ratio < -1 ~ paste0(subB, "-biased"),
        TRUE ~ "Unbiased"
      )
    )
  
  hist_data <- hist(avg_data$Log2_Ratio, breaks = 100, plot = FALSE)
  max_count <- max(hist_data$counts)
  
  stats <- avg_data %>%
    summarise(
      SD = sd(Log2_Ratio, na.rm = TRUE),
      SubA_biased = sum(Log2_Ratio > 1) / n() * 100,
      SubB_biased = sum(Log2_Ratio < -1) / n() * 100,
      Unbiased = sum(Log2_Ratio >= -1 & Log2_Ratio <= 1) / n() * 100,
      Weighted_Avg_Log2 = sum(Log2_Ratio * (SubA_expr + SubB_expr), na.rm = TRUE) / sum(SubA_expr + SubB_expr, na.rm = TRUE)
    ) %>%
    mutate(
      Label = sprintf("SD: %.2f\n%s-biased: %.1f%%\n%s-biased: %.1f%%\nUnbiased: %.1f%%",
                      SD, subA, SubA_biased, subB, SubB_biased, Unbiased),
      Y_Position = max_count * 0.9
    )
  
  cat(sprintf("\nSummary for Averaged Histogram (%s vs %s):\n", subA, subB))
  print(stats)
  
  x_label <- substitute(log[2] ~ (a/b), list(a = subA, b = subB))
  
  p <- ggplot(avg_data, aes(x = Log2_Ratio, fill = Bias)) +
    geom_histogram(bins = 100, aes(color = after_stat(fill)), linewidth = 0.2, alpha = 0.7) +
    scale_fill_manual(values = setNames(
      c("#1f77b4", "#ff7f0e", "#808080"),
      c(paste0(subA, "-biased"), paste0(subB, "-biased"), "Unbiased")
    )) +
    scale_color_manual(values = setNames(
      c("#1f77b4", "#ff7f0e", "#808080"),
      c(paste0(subA, "-biased"), paste0(subB, "-biased"), "Unbiased")
    )) +
    guides(color = "none") +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "red", linewidth = 1) +
    geom_text(data = stats, aes(x = Inf, y = Y_Position, label = Label), 
              hjust = 1, vjust = 1, size = 3.5, color = "black", inherit.aes = FALSE) +
    scale_x_continuous(breaks = seq(-15, 15, by = 5), limits = c(-15, 15)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.3))) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      axis.line = element_line(color = "black"),
      axis.text = element_text(color = "black", size = 10, face = "bold"),
      axis.title = element_text(size = 12, face = "bold"),
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      strip.text = element_text(size = 12, face = "bold"),
      panel.spacing = unit(1, "lines"),
      legend.position = "right",
      legend.title = element_text(size = 10, face = "bold"),
      legend.text = element_text(size = 8)
    ) +
    labs(title = paste(subA, "vs", subB), x = x_label, y = "Count", fill = "Bias")
  
  return(p)
}

# Analysis and Visualization
homoeolog_ratios <- calculate_homoeolog_ratios(expression_data, pseudocount)
dominance_classification <- classify_homoeolog_dominance(homoeolog_ratios, pseudocount)

# 1. Combined Histograms (Per Tissue)
comp_list <- list(c("H", "St1"), c("H", "St2"), c("St1", "St2"))
combined_histogram_plots <- lapply(comp_list, function(comp) {
  p <- generate_combined_histogram(homoeolog_ratios, dominance_classification, comp[1], comp[2], pseudocount)
  ggsave(paste0("histogram_", comp[1], "_vs_", comp[2], "_across_tissues.jpg"), p, width = 10, height = 6, dpi = 600)
  ggsave(paste0("histogram_", comp[1], "_vs_", comp[2], "_across_tissues.pdf"), p, width = 10, height = 6)
  return(p)
})

# 2. Averaged Histograms (Across All Tissues)
avg_histogram_plots <- lapply(comp_list, function(comp) {
  p <- generate_avg_histogram(homoeolog_ratios, dominance_classification, comp[1], comp[2], pseudocount)
  ggsave(paste0("histogram_avg_", comp[1], "_vs_", comp[2], "_all_tissues.jpg"), p, width = 8, height = 6, dpi = 600)
  ggsave(paste0("histogram_avg_", comp[1], "_vs_", comp[2], "_all_tissues.pdf"), p, width = 8, height = 6)
  return(p)
})

# 3. Combined Figure for Averaged Histograms (Panels a, b, c)
avg_histogram_plots[[1]] <- avg_histogram_plots[[1]] + labs(tag = "a")
avg_histogram_plots[[2]] <- avg_histogram_plots[[2]] + labs(tag = "b")
avg_histogram_plots[[3]] <- avg_histogram_plots[[3]] + labs(tag = "c")

combined_avg_plot <- avg_histogram_plots[[1]] / avg_histogram_plots[[2]] / avg_histogram_plots[[3]]
ggsave("histogram_avg_all_comparisons_combined.jpg", combined_avg_plot, width = 8, height = 10, dpi = 600)
ggsave("histogram_avg_all_comparisons_combined.pdf", combined_avg_plot, width = 8, height = 10)

# Session Info
sessionInfo()
