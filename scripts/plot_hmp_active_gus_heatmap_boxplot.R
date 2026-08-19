#!/usr/bin/env Rscript
# Publication-style figures for HMP active-GUS abundance:
#   (1) filled-tile TPM heatmap, unclustered (fixed biological order)
#   (2) filled-tile TPM heatmap, hierarchically clustered (row + column dendrograms)
#   (3) boxplot of genome-wide TPM per active GUS, one point per person
#
# Input : results/hmp_active_gus/report/active_gus_sample_by_reference.tsv
# Output: results/hmp_active_gus/report/figures/05_* 06_* 07_*  (png + pdf)

Sys.setenv(OPENBLAS_NUM_THREADS = "1", OMP_NUM_THREADS = "1")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
  library(ggplot2)
  library(scales)
  library(forcats)
  library(patchwork)
  library(ggbeeswarm)
})

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args) >= 1) args[[1]] else getwd()
project_root <- normalizePath(project_root, mustWork = TRUE)

in_path  <- file.path(project_root, "results/hmp_active_gus/report/active_gus_sample_by_reference.tsv")
fig_dir  <- file.path(project_root, "results/hmp_active_gus/report/figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

dat <- read_tsv(in_path, show_col_types = FALSE) %>%
  mutate(
    sample          = as.character(sample),
    gus_label       = as.character(gus_label),
    tpm_genome_wide = as.numeric(tpm_genome_wide),
    present         = as.integer(present)
  )

value_col <- "tpm_genome_wide"   # the "TPM" shown throughout the report

# ---- Wide matrix: rows = active GUS, columns = people -----------------------
mat <- dat %>%
  select(gus_label, sample, all_of(value_col)) %>%
  pivot_wider(names_from = sample, values_from = all_of(value_col), values_fill = 0) %>%
  column_to_rownames("gus_label") %>%
  as.matrix()

# ---- Orders -----------------------------------------------------------------
# Biological (unclustered) row order: by clade then source, matching the panel.
row_meta <- dat %>%
  distinct(gus_label, clade, source) %>%
  arrange(clade, source, gus_label)
bio_row_order <- row_meta$gus_label

# Unclustered column order: people with more active GUS first, then higher total TPM.
col_stats <- dat %>%
  group_by(sample) %>%
  summarise(n_active = sum(present), total_tpm = sum(.data[[value_col]]), .groups = "drop") %>%
  arrange(desc(n_active), desc(total_tpm), sample)
bio_col_order <- col_stats$sample

# ---- Hierarchical clustering (on log1p TPM) ---------------------------------
lmat <- log1p(mat)
# rows: cluster all 11 references
row_hc <- hclust(dist(lmat, method = "euclidean"), method = "ward.D2")
# columns: drop all-zero people for a meaningful dendrogram, append them at the end
col_nonzero <- colSums(mat) > 0
if (sum(col_nonzero) >= 2) {
  col_hc <- hclust(dist(t(lmat[, col_nonzero, drop = FALSE]), method = "euclidean"),
                   method = "ward.D2")
  clust_col_order <- c(colnames(mat)[col_nonzero][col_hc$order],
                       colnames(mat)[!col_nonzero])
} else {
  col_hc <- NULL
  clust_col_order <- colnames(mat)
}
clust_row_order <- rownames(mat)[row_hc$order]

# ---- Dendrogram -> segment coordinates --------------------------------------
# Returns segments in (pos, ht) space where `pos` is the leaf axis (1..n in the
# dendrogram leaf order) and `ht` is merge height.
dend_segments <- function(hc) {
  n <- length(hc$order)
  leaf_pos <- numeric(n)
  leaf_pos[hc$order] <- seq_len(n)
  node_x <- numeric(nrow(hc$merge))
  get <- function(idx) if (idx < 0) list(x = leaf_pos[-idx], h = 0)
                       else          list(x = node_x[idx],    h = hc$height[idx])
  segs <- vector("list", nrow(hc$merge) * 3L)
  k <- 0L
  for (i in seq_len(nrow(hc$merge))) {
    L <- get(hc$merge[i, 1]); R <- get(hc$merge[i, 2]); h <- hc$height[i]
    segs[[k <- k + 1L]] <- data.frame(pos = L$x, posend = L$x, ht = L$h, htend = h)
    segs[[k <- k + 1L]] <- data.frame(pos = R$x, posend = R$x, ht = R$h, htend = h)
    segs[[k <- k + 1L]] <- data.frame(pos = L$x, posend = R$x, ht = h,   htend = h)
    node_x[i] <- (L$x + R$x) / 2
  }
  do.call(rbind, segs)
}

# ---- Shared fill scale ------------------------------------------------------
# Absent (TPM == 0) -> light grey; present -> plasma gradient on a log scale.
tpm_max  <- max(mat)
fill_scale <- scale_fill_viridis_c(
  option    = "C",
  name      = "Genome-wide TPM",
  trans     = "log10",
  limits    = c(1, tpm_max),
  breaks    = c(1, 2, 5, 10, 20, 50),
  labels    = c("1", "2", "5", "10", "20", "50"),
  na.value  = "grey90",
  guide     = guide_colorbar(barheight = unit(9, "lines"), barwidth = unit(0.7, "lines"))
)

# Long data frame keyed by integer positions so tiles align with dendrograms.
make_long <- function(row_order, col_order) {
  dat %>%
    mutate(
      gy = match(gus_label, row_order),
      cx = match(sample, col_order),
      fillv = ifelse(.data[[value_col]] > 0, .data[[value_col]], NA_real_)
    ) %>%
    filter(!is.na(gy), !is.na(cx))
}

heatmap_base_theme <- theme_minimal(base_size = 12) +
  theme(
    panel.grid      = element_blank(),
    axis.title      = element_blank(),
    axis.text.x     = element_blank(),
    axis.ticks.x    = element_blank(),
    axis.text.y     = element_text(size = 10, color = "#222222"),
    plot.title      = element_text(face = "bold"),
    plot.subtitle   = element_text(color = "#555555", size = 9.5),
    legend.title    = element_text(size = 10),
    plot.margin     = margin(4, 6, 4, 6)
  )

n_samples <- ncol(mat)
n_refs    <- nrow(mat)

# =============================================================================
# (1) Unclustered heatmap
# =============================================================================
long_bio <- make_long(bio_row_order, bio_col_order)

p_hm_nocl <- ggplot(long_bio, aes(x = cx, y = gy, fill = fillv)) +
  geom_tile(color = "white", linewidth = 0.15) +
  fill_scale +
  scale_x_continuous(expand = c(0, 0), limits = c(0.5, n_samples + 0.5)) +
  scale_y_continuous(
    expand = c(0, 0), limits = c(0.5, n_refs + 0.5),
    breaks = seq_len(n_refs), labels = bio_row_order
  ) +
  heatmap_base_theme

# =============================================================================
# (2) Clustered heatmap with dendrograms
# =============================================================================
long_cl <- make_long(clust_row_order, clust_col_order)

p_hm_core <- ggplot(long_cl, aes(x = cx, y = gy, fill = fillv)) +
  geom_tile(color = "white", linewidth = 0.15) +
  fill_scale +
  scale_x_continuous(expand = c(0, 0), limits = c(0.5, n_samples + 0.5)) +
  scale_y_continuous(
    expand = c(0, 0), limits = c(0.5, n_refs + 0.5),
    breaks = seq_len(n_refs), labels = clust_row_order, position = "right"
  ) +
  heatmap_base_theme +
  theme(plot.margin = margin(2, 6, 4, 2))

# Column dendrogram (top). Positions 1..k map onto the clustered non-zero block.
blank_axes <- theme_void() + theme(plot.margin = margin(0, 0, 0, 0))

if (!is.null(col_hc)) {
  col_seg <- dend_segments(col_hc)  # pos in 1..k over the non-zero people
  p_col_dendro <- ggplot(col_seg) +
    geom_segment(aes(x = pos, xend = posend, y = ht, yend = htend),
                 linewidth = 0.3, color = "#333333", lineend = "round") +
    scale_x_continuous(expand = c(0, 0), limits = c(0.5, n_samples + 0.5)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
    blank_axes
} else {
  p_col_dendro <- ggplot() + blank_axes
}

# Row dendrogram (left). Height on x (reversed so root points left), leaves on y.
row_seg <- dend_segments(row_hc)
p_row_dendro <- ggplot(row_seg) +
  geom_segment(aes(x = ht, xend = htend, y = pos, yend = posend),
               linewidth = 0.3, color = "#333333", lineend = "round") +
  scale_x_reverse(expand = expansion(mult = c(0.02, 0))) +
  scale_y_continuous(expand = c(0, 0), limits = c(0.5, n_refs + 0.5)) +
  blank_axes

# 2x2 grid: A = spacer, B = column dendrogram, C = row dendrogram, D = heatmap.
p_clustered <- patchwork::wrap_plots(
  A = plot_spacer(),
  B = p_col_dendro,
  C = p_row_dendro,
  D = p_hm_core,
  design  = "AB\nCD",
  widths  = c(0.13, 1),
  heights = c(0.16, 1),
  guides  = "collect"
)

# =============================================================================
# (3) Boxplot: TPM per active GUS, one point per person (present only)
# =============================================================================
present_dat <- dat %>%
  filter(.data[[value_col]] > 0) %>%
  mutate(tpm = .data[[value_col]])

box_order <- present_dat %>%
  group_by(gus_label) %>%
  summarise(n = n(), med = median(tpm), .groups = "drop") %>%
  arrange(desc(n), desc(med))

# Keep every reference on the axis (even those never detected).
all_labels <- row_meta$gus_label
present_dat <- present_dat %>%
  mutate(gus_label = factor(gus_label, levels = box_order$gus_label))

n_label <- dat %>%
  group_by(gus_label) %>%
  summarise(n_present = sum(present), .groups = "drop")

box_df <- present_dat
box_axis_levels <- box_order$gus_label
# Detected counts for annotation, aligned to plotted order.
count_annot <- n_label %>%
  filter(gus_label %in% box_axis_levels) %>%
  mutate(gus_label = factor(gus_label, levels = box_axis_levels))

p_box <- ggplot(box_df, aes(x = gus_label, y = tpm)) +
  geom_boxplot(width = 0.6, outlier.shape = NA, fill = "#DCE6EC",
               color = "#2F6F8F", linewidth = 0.5, alpha = 0.9) +
  ggbeeswarm::geom_quasirandom(width = 0.18, size = 1.9, shape = 21,
                               fill = "#E15759", color = "white", stroke = 0.3,
                               alpha = 0.9) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.06))) +
  labs(x = NULL, y = "Genome-wide TPM") +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1, size = 10),
    plot.margin = margin(8, 12, 8, 8)
  )

# =============================================================================
# Save
# =============================================================================
f_nocl <- file.path(fig_dir, "05_abundance_heatmap_tpm_nocluster.png")
f_cl   <- file.path(fig_dir, "06_abundance_heatmap_tpm_clustered.png")
f_box  <- file.path(fig_dir, "07_tpm_boxplot_by_active_gus.png")

ggsave(f_nocl, p_hm_nocl,     width = 15, height = 5.2, dpi = 200, bg = "white")
ggsave(f_cl,   p_clustered,   width = 15, height = 5.6, dpi = 200, bg = "white")
ggsave(f_box,  p_box,         width = 10, height = 6.0, dpi = 200, bg = "white")

ggsave(sub("\\.png$", ".pdf", f_nocl), p_hm_nocl,   width = 15, height = 5.2, bg = "white")
ggsave(sub("\\.png$", ".pdf", f_cl),   p_clustered, width = 15, height = 5.6, bg = "white")
ggsave(sub("\\.png$", ".pdf", f_box),  p_box,       width = 10, height = 6.0, bg = "white")

cat("Wrote figures:\n  ", f_nocl, "\n  ", f_cl, "\n  ", f_box, "\n")
