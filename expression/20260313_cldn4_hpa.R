## CLDN4 HPA Normal vs Lung Cancer Distribution
#by e-bruli
#20260310

# !! HPAanalyze v1.24.0 required (bundles HPA v22.0 data) !!

#Package Version : 1.24.0
#HPA Version: 22.0
# If not installed, run:
# url  <- "https://bioconductor.org/packages/3.20/bioc/src/contrib/HPAanalyze_1.24.0.tar.gz"
# dest <- tempfile(fileext = ".tar.gz")
# download.file(url, dest)
# install.packages(dest, repos = NULL, type = "source")


# ══════════════════════════════════════════════════════════════════
# CLDN4 HPA IHC Staining Plot
# Usage: Rscript cldn4_HPA_IHC.R /path/to/output/dir
# ══════════════════════════════════════════════════════════════════

# ── Install if needed ─────────────────────────────────────────────
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#if (!requireNamespace("HPAanalyze", quietly = TRUE))
#  BiocManager::install("HPAanalyze")

library(HPAanalyze)
packageVersion("HPAanalyze")
# ── Get output directory from command line argument ───────────────
args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  stop("Please provide an output directory path.
    Usage: Rscript cldn4_HPA_IHC.R /path/to/output/dir")
}

OUT_DIR <- args[1]

if (!dir.exists(OUT_DIR)) {
  dir.create(OUT_DIR, recursive = TRUE)
  cat("Created output directory:", OUT_DIR, "\n")
}
cat("Output will be saved to:", OUT_DIR, "\n")

# ── Load built-in HPA histology data ─────────────────────────────
data("hpa_histology_data")
cat("HPA data loaded. Version:", hpa_histology_data$metadata$HPAversion, "\n")

# ── Filter CLDN4 ─────────────────────────────────────────────────
cldn4_cancer <- hpa_histology_data$pathology[
  hpa_histology_data$pathology$gene == "CLDN4" &
    hpa_histology_data$pathology$cancer == "lung cancer", ]

cldn4_normal <- hpa_histology_data$normal_tissue[
  hpa_histology_data$normal_tissue$gene == "CLDN4" &
    hpa_histology_data$normal_tissue$tissue == "lung", ]

cat("\nCLDN4 Lung Cancer IHC:\n"); print(cldn4_cancer)
cat("\nCLDN4 Normal Lung IHC:\n"); print(cldn4_normal)

# ── Prepare cancer data ───────────────────────────────────────────
cancer_long <- data.frame(
  group = "Lung Cancer",
  level = c("High", "Medium", "Low", "Not detected"),
  count = c(cldn4_cancer$high, cldn4_cancer$medium,
            cldn4_cancer$low,  cldn4_cancer$not_detected)
)

# ── Prepare normal data ───────────────────────────────────────────
normal_counts <- table(cldn4_normal$level)
normal_long <- data.frame(
  group = "Normal Lung",
  level = names(normal_counts),
  count = as.integer(normal_counts)
)

# Add missing levels with count 0
all_levels <- c("High", "Medium", "Low", "Not detected")
missing <- setdiff(all_levels, normal_long$level)
if (length(missing) > 0) {
  normal_long <- rbind(normal_long,
                       data.frame(group = "Normal Lung", level = missing, count = 0))
}

# ── Combine & order ───────────────────────────────────────────────
df_plot <- rbind(cancer_long, normal_long)
df_plot$level <- factor(df_plot$level, levels = all_levels)
df_plot$group <- factor(df_plot$group, levels = c("Normal Lung", "Lung Cancer"))

# ── Calculate percentages ─────────────────────────────────────────
df_plot$pct <- ave(df_plot$count, df_plot$group,
                   FUN = function(x) x / sum(x) * 100)

# ── Color scheme options — uncomment the one you want ─────────────

# Option 1: Clinical Red Scale
# level_colors <- c(
#     "High"         = "#8B0000",  # dark red
#     "Medium"       = "#E07B54",  # salmon
#     "Low"          = "#F5CBA7",  # peach
#     "Not detected" = "#D3D3D3"   # grey
# )

# Option 2: Viridis-inspired (colorblind friendly)
# level_colors <- c(
#     "High"         = "#440154",  # dark purple
#     "Medium"       = "#31688E",  # blue
#     "Low"          = "#90D743",  # green
#     "Not detected" = "#ECECEC"   # light grey
# )

# Option 3: Okabe-Ito (Nature recommended, colorblind friendly)
level_colors <- c(
     "High"         = "#E69F00",  # orange
     "Medium"       = "#56B4E9",  # sky blue
     "Low"          = "#009E73",  # green
     "Not detected" = "#ECECEC"   # light grey
 )

# Option 4: Paul Tol Sequential (colorblind friendly) ← ACTIVE
#level_colors <- c(
 # "High"         = "#364B9A",  # dark blue
  #"Medium"       = "#6EA6CD",  # medium blue
  #"Low"          = "#C2E4EF",  # light blue
  #"Not detected" = "#ECECEC"   # grey
#)

# Option 5: Warm Yellow→Red (sequential, common in IHC papers)
# level_colors <- c(
#     "High"         = "#B03A2E",  # dark red
#     "Medium"       = "#E67E22",  # orange
#     "Low"          = "#F9E79F",  # yellow
#     "Not detected" = "#ECECEC"   # grey
# )

# Option 6: Greyscale (for journals without color)
# level_colors <- c(
#     "High"         = "#1C1C1C",  # near black
#     "Medium"       = "#666666",  # dark grey
#     "Low"          = "#BBBBBB",  # light grey
#     "Not detected" = "#ECECEC"   # near white
# )

# ── Sample counts per group (for annotation) ─────────────────────
sample_counts <- aggregate(count ~ group, data = df_plot, FUN = sum)
# Create x-axis labels with n= underneath
sample_counts$label <- paste0(sample_counts$group,
                              "\n(n=", sample_counts$count, ")")
# Match order to factor levels
group_labels <- setNames(sample_counts$label, sample_counts$group)

# ── Plot ──────────────────────────────────────────────────────────
p <- ggplot2::ggplot(df_plot,
                     ggplot2::aes(x = group, y = pct, fill = level)) +
  ggplot2::geom_bar(stat = "identity", width = 0.6,
                    color = "white", linewidth = 0.3) +
  ggplot2::geom_text(
    ggplot2::aes(label = ifelse(pct > 5, paste0(round(pct), "%"), "")),
    position = ggplot2::position_stack(vjust = 0.5),
    size = 3.5, color = "white", fontface = "bold") +
  ggplot2::scale_fill_manual(values = level_colors,
                             name = "IHC Staining Level") +
  ggplot2::scale_x_discrete(labels = group_labels) +  # n= in x labels
  ggplot2::scale_y_continuous(expand = c(0, 0),
                              labels = function(x) paste0(x, "%")) +
  ggplot2::labs(
    title    = "CLDN4 Protein Expression (IHC)",
    subtitle = "Human Protein Atlas v22 — Normal Lung vs. Lung Cancer",
    x        = NULL,
    y        = "Percentage of samples (%)"
  ) +
  ggplot2::theme_classic(base_size = 13) +
  ggplot2::theme(
    plot.title         = ggplot2::element_text(face = "bold", size = 14),
    plot.subtitle      = ggplot2::element_text(color = "grey50", size = 10),
    legend.position    = "right",
    axis.text.x        = ggplot2::element_text(size = 12, face = "bold"),
    panel.grid.major.y = ggplot2::element_line(color = "grey90")
  )

print(p)

# ── Save ──────────────────────────────────────────────────────────
ggplot2::ggsave(file.path(OUT_DIR, "cldn4_HPA_IHC.pdf"),
                p, width = 6, height = 5, dpi = 300)
ggplot2::ggsave(file.path(OUT_DIR, "cldn4_HPA_IHC.png"),
                p, width = 6, height = 5, dpi = 300)

cat("Saved: cldn4_HPA_IHC.pdf and cldn4_HPA_IHC.png\n")


