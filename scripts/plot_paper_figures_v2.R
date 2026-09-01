#!/usr/bin/env Rscript
# Paper figures v2: restyle of the original ggplot2 panels.
#   - Arial, larger type
#   - no plot titles / subtitles on data panels
#   - writes a NEW folder; does not overwrite the originals
#
# Data  : results/hmp_assembly_gus/report/*.tsv
# Output: results/hmp_assembly_gus/report/figures_v2/
#
# Based on scripts/report_hmp_assembly_gus.R and
#          scripts/plot_hmp_active_gus_heatmap_boxplot.R

Sys.setenv(OPENBLAS_NUM_THREADS = "1", OMP_NUM_THREADS = "1")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
  library(ggplot2)
  library(scales)
  library(systemfonts)
  library(grid)
  library(ggtext)
})

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args) >= 1) args[[1]] else getwd()
project_root <- normalizePath(project_root, mustWork = TRUE)

in_dir  <- file.path(project_root, "results/hmp_assembly_gus/report")
fig_dir <- file.path(in_dir, "figures_v2")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

font_dir <- file.path(project_root, "scripts/fonts")
if (!"Arial" %in% systemfonts::system_fonts()$family) {
  register_font(
    "Arial",
    plain = file.path(font_dir, "Arial.ttf"),
    bold  = file.path(font_dir, "Arial_Bold.ttf")
  )
}

# Embed TrueType in PDF so Arial stays editable in Illustrator.
if (capabilities("cairo")) {
  pdf_dev <- cairo_pdf
} else {
  pdf_dev <- grDevices::pdf
}

save_panel <- function(plot, stem, width, height) {
  png_path <- file.path(fig_dir, paste0(stem, ".png"))
  pdf_path <- file.path(fig_dir, paste0(stem, ".pdf"))
  ggsave(
    png_path, plot,
    width = width, height = height, dpi = 300, bg = "white",
    device = ragg::agg_png
  )
  ggsave(
    pdf_path, plot,
    width = width, height = height, bg = "white",
    device = pdf_dev, family = "Arial"
  )
  cat("  ", basename(png_path), "  ", basename(pdf_path), "\n", sep = "")
}

# -----------------------------------------------------------------------------
# Data (same objects the original report used to draw 01 / 04 / 05)
# -----------------------------------------------------------------------------
grid <- read_tsv(file.path(in_dir, "active_gus_sample_by_reference.tsv"), show_col_types = FALSE) %>%
  mutate(
    sample          = as.character(sample),
    gus_label       = as.character(gus_label),
    clade           = as.character(clade),
    source          = as.character(source),
    tpm_genome_wide = as.numeric(tpm_genome_wide),
    present         = as.integer(present)
  )

count_dist <- read_tsv(file.path(in_dir, "active_gus_count_distribution.tsv"), show_col_types = FALSE) %>%
  mutate(n_active_f = factor(n_active, levels = 0:11))

sample_order_df <- grid %>%
  group_by(sample) %>%
  summarise(
    n_active  = sum(present),
    total_tpm = sum(tpm_genome_wide),
    .groups   = "drop"
  ) %>%
  arrange(desc(n_active), desc(total_tpm), sample) %>%
  mutate(sample_f = factor(sample, levels = sample))

row_meta <- grid %>%
  distinct(gus_label, clade, source) %>%
  arrange(clade, source, gus_label)

gus_label_levels <- row_meta$gus_label
grid <- grid %>%
  mutate(
    gus_label = factor(gus_label, levels = gus_label_levels),
    sample_f  = factor(sample, levels = levels(sample_order_df$sample_f))
  )

rel <- grid %>%
  group_by(sample) %>%
  mutate(
    total_tpm = sum(tpm_genome_wide),
    pct = if_else(total_tpm > 0, 100 * tpm_genome_wide / total_tpm, 0)
  ) %>%
  mutate(
    pct = if_else(total_tpm > 0 & sum(pct) > 0, 100 * pct / sum(pct), 0)
  ) %>%
  ungroup()

distinct_cols <- c(
  "#4E79A7", "#F28E2B", "#E15759", "#76B7B2", "#59A14F",
  "#EDC948", "#B07AA1", "#FF9DA7", "#9C755F", "#BAB0AC", "#86BCB6"
)
gus_colors <- setNames(distinct_cols[seq_len(length(gus_label_levels))], gus_label_levels)

n_people <- n_distinct(grid$sample)
n_refs   <- length(gus_label_levels)
n_zero   <- sum(sample_order_df$n_active == 0)
n_any    <- sum(sample_order_df$n_active > 0)
n_max    <- max(sample_order_df$n_active)

# =============================================================================
# 01  Count distribution  (title/subtitle removed)
# =============================================================================
p_count <- ggplot(count_dist, aes(x = n_active_f, y = n_people)) +
  geom_col(fill = "#2F6F8F", width = 0.72) +
  geom_text(
    aes(label = n_people),
    vjust = -0.35, size = 5, color = "#222222", family = "Arial"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.14))) +
  labs(
    x = "Number of distinct active GUS detected",
    y = "Number of people"
  ) +
  theme_classic(base_size = 16, base_family = "Arial") +
  theme(
    plot.title    = element_blank(),
    plot.subtitle = element_blank(),
    axis.title    = element_text(size = 16, color = "black"),
    axis.text     = element_text(size = 14, color = "black"),
    axis.line     = element_line(linewidth = 0.5),
    axis.ticks    = element_line(linewidth = 0.5)
  )

# =============================================================================
# 04  Stacked relative abundance  (title/subtitle removed)
# =============================================================================
bdorei_label <- "B.dorei (NL)"
stack_samples <- sample_order_df %>%
  filter(n_active > 0) %>%
  select(sample)

stack_order <- rel %>%
  semi_join(stack_samples, by = "sample") %>%
  mutate(bdorei_pct = if_else(as.character(gus_label) == bdorei_label, pct, 0)) %>%
  group_by(sample) %>%
  summarise(bdorei_pct = sum(bdorei_pct), .groups = "drop") %>%
  arrange(bdorei_pct, sample)

stack_all <- rel %>%
  semi_join(stack_samples, by = "sample") %>%
  mutate(
    gus_label = factor(as.character(gus_label), levels = gus_label_levels),
    sample_f  = factor(sample, levels = stack_order$sample)
  )

p_stack <- ggplot(stack_all, aes(x = sample_f, y = pct, fill = gus_label)) +
  geom_col(width = 0.92, color = NA, position = "stack") +
  scale_fill_manual(values = gus_colors, name = "Active GUS", drop = FALSE) +
  scale_y_continuous(breaks = seq(0, 100, 25), expand = c(0, 0)) +
  coord_cartesian(ylim = c(0, 100), clip = "off") +
  labs(x = NULL, y = "Relative abundance (%)") +
  theme_classic(base_size = 16, base_family = "Arial") +
  theme(
    plot.title    = element_blank(),
    plot.subtitle = element_blank(),
    axis.title.y  = element_text(size = 15, color = "black"),
    axis.text.y   = element_text(size = 13, color = "black"),
    axis.text.x   = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 7, color = "black"),
    axis.ticks.x  = element_line(linewidth = 0.4),
    axis.line     = element_line(linewidth = 0.5),
    legend.position = "right",
    legend.title  = element_text(size = 14),
    legend.text   = element_text(size = 12),
    legend.key.size = unit(0.42, "cm"),
    plot.margin    = margin(8, 8, 12, 8)
  )

# =============================================================================
# 05  Unclustered TPM heatmap  (same layout as plot_hmp_active_gus_heatmap_boxplot.R)
# =============================================================================
value_col <- "tpm_genome_wide"
mat <- grid %>%
  select(gus_label, sample, all_of(value_col)) %>%
  pivot_wider(names_from = sample, values_from = all_of(value_col), values_fill = 0) %>%
  column_to_rownames("gus_label") %>%
  as.matrix()

bio_row_order <- row_meta$gus_label
bio_col_order <- sample_order_df$sample
n_samples <- ncol(mat)
n_refs_hm <- nrow(mat)
tpm_max <- max(mat)

fill_scale <- scale_fill_viridis_c(
  option   = "C",
  name     = "Genome-wide TPM",
  trans    = "log10",
  limits   = c(1, tpm_max),
  breaks   = c(1, 2, 5, 10, 20, 50),
  labels   = c("1", "2", "5", "10", "20", "50"),
  na.value = "grey90",
  guide    = guide_colorbar(
    barheight   = unit(11, "lines"),
    barwidth    = unit(0.85, "lines"),
    title.theme = element_text(family = "Arial", size = 14),
    label.theme = element_text(family = "Arial", size = 12)
  )
)

long_bio <- grid %>%
  mutate(
    gy    = match(as.character(gus_label), bio_row_order),
    cx    = match(sample, bio_col_order),
    fillv = ifelse(.data[[value_col]] > 0, .data[[value_col]], NA_real_)
  ) %>%
  filter(!is.na(gy), !is.na(cx))

p_hm_nocl <- ggplot(long_bio, aes(x = cx, y = gy, fill = fillv)) +
  geom_tile(color = "white", linewidth = 0.15) +
  fill_scale +
  scale_x_continuous(expand = c(0, 0), limits = c(0.5, n_samples + 0.5)) +
  scale_y_continuous(
    expand = c(0, 0), limits = c(0.5, n_refs_hm + 0.5),
    breaks = seq_len(n_refs_hm), labels = bio_row_order
  ) +
  theme_minimal(base_size = 16, base_family = "Arial") +
  theme(
    panel.grid       = element_blank(),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background  = element_rect(fill = "white", colour = NA),
    axis.title       = element_blank(),
    axis.text.x   = element_blank(),
    axis.ticks.x  = element_blank(),
    axis.text.y   = element_text(size = 13, color = "#222222", family = "Arial"),
    legend.title  = element_text(size = 14, family = "Arial"),
    legend.text   = element_text(size = 12, family = "Arial"),
    plot.title    = element_blank(),
    plot.subtitle = element_blank(),
    plot.margin   = margin(6, 10, 6, 6)
  )

# =============================================================================
# 00  Analysis scheme  (assembly pipeline + active GUS only)
# =============================================================================
# Coordinate system in ggplot (x 0–12, y 0–10). Boxes are drawn with grid
# roundrects so the schematic stays in R / ggplot2 rather than a second stack.
round_box <- function(xmin, xmax, ymin, ymax, fill, colour, lwd = 1.2) {
  annotation_custom(
    grob = roundrectGrob(
      r  = unit(0.12, "npc"),
      gp = gpar(fill = fill, col = colour, lwd = lwd)
    ),
    xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax
  )
}

label_box <- function(x, y, title, subtitle, colour, title_size = 4.6, sub_size = 3.7) {
  list(
    annotate(
      "text", x = x, y = y + 0.18, label = title,
      family = "Arial", fontface = "bold", size = title_size, colour = colour
    ),
    annotate(
      "text", x = x, y = y - 0.28, label = subtitle,
      family = "Arial", size = sub_size, colour = colour, lineheight = 0.95
    )
  )
}

down_arrow <- function(x, y1, y2) {
  annotate(
    "segment", x = x, xend = x, y = y1, yend = y2,
    arrow = arrow(length = unit(0.16, "cm"), type = "closed"),
    colour = "#2B3A4A", linewidth = 0.7
  )
}

diag_arrow <- function(x1, y1, x2, y2) {
  annotate(
    "segment", x = x1, xend = x2, y = y1, yend = y2,
    arrow = arrow(length = unit(0.16, "cm"), type = "closed"),
    colour = "#2B3A4A", linewidth = 0.7
  )
}

navy <- "#1F4E79"
ink  <- "#1C2833"
teal_fill <- "#E8F2F5"
teal_edge <- "#2F6F8F"

p_scheme <- ggplot() +
  # 1 input
  round_box(3.35, 8.65, 8.45, 9.70, fill = navy, colour = navy) +
  label_box(
    6, 9.18,
    "HMP stool shotgun metagenomes",
    "139 HMASM samples  ·  paired trimmed FASTQ",
    "white", title_size = 4.7, sub_size = 3.8
  ) +
  down_arrow(6, 8.45, 7.75) +
  # 2 assembly
  round_box(3.35, 8.65, 6.45, 7.75, fill = "white", colour = navy) +
  label_box(
    6, 7.20,
    "MEGAHIT de novo assembly",
    "contigs ≥ 1,000 bp",
    ink, title_size = 4.7, sub_size = 3.8
  ) +
  diag_arrow(4.55, 6.45, 2.70, 5.55) +
  diag_arrow(7.45, 6.45, 9.30, 5.55) +
  annotate("text", x = 2.70, y = 5.68, label = "proteins",
           family = "Arial", fontface = "italic", size = 3.3, colour = "#4A5560") +
  annotate("text", x = 9.30, y = 5.68, label = "reads",
           family = "Arial", fontface = "italic", size = 3.3, colour = "#4A5560") +
  # 3 prodigal
  round_box(0.45, 4.95, 4.15, 5.50, fill = teal_fill, colour = teal_edge) +
  label_box(
    2.70, 4.92,
    "Gene prediction",
    "Prodigal  ·  meta mode",
    ink, title_size = 4.5, sub_size = 3.7
  ) +
  # 4 mapping
  round_box(7.05, 11.55, 4.15, 5.50, fill = teal_fill, colour = teal_edge) +
  label_box(
    9.30, 4.92,
    "Read mapping and counting",
    "Bowtie2: reads → same assembly\nfeatureCounts on predicted CDS",
    ink, title_size = 4.4, sub_size = 3.5
  ) +
  down_arrow(2.70, 4.15, 3.40) +
  # 5 diamond
  round_box(0.45, 4.95, 1.85, 3.40, fill = "white", colour = navy) +
  label_box(
    2.70, 2.72,
    "Active GUS homology search",
    "DIAMOND blastp vs 11 active GUS\n≥ 95% identity, ≥ 90% coverage",
    ink, title_size = 4.4, sub_size = 3.5
  ) +
  # merge to analysis
  diag_arrow(2.70, 1.85, 4.55, 1.35) +
  annotate(
    "segment", x = 9.30, xend = 9.30, y = 4.15, yend = 1.70,
    colour = "#2B3A4A", linewidth = 0.7
  ) +
  diag_arrow(9.30, 1.70, 7.45, 1.35) +
  # 6 analysis
  round_box(3.35, 8.65, 0.15, 1.35, fill = navy, colour = navy) +
  label_box(
    6, 0.82,
    "Active GUS abundance analysis",
    "best hit per gene  ·  genome-wide TPM\npresence, abundance, and composition",
    "white", title_size = 4.5, sub_size = 3.5
  ) +
  # step numbers
  geom_point(
    data = data.frame(
      x = c(3.55, 3.55, 0.68, 7.28, 0.68, 3.55),
      y = c(9.52, 7.57, 5.32, 5.32, 3.22, 1.17)
    ),
    aes(x, y),
    inherit.aes = FALSE, size = 6.2, colour = navy, fill = navy, shape = 21, stroke = 0.6
  ) +
  geom_text(
    data = data.frame(
      x = c(3.55, 3.55, 0.68, 7.28, 0.68, 3.55),
      y = c(9.52, 7.57, 5.32, 5.32, 3.22, 1.17),
      lab = as.character(1:6)
    ),
    aes(x, y, label = lab),
    inherit.aes = FALSE, family = "Arial", fontface = "bold", size = 3.3, colour = "white"
  ) +
  coord_cartesian(xlim = c(0, 12), ylim = c(0, 10), expand = FALSE) +
  theme_void(base_family = "Arial") +
  theme(plot.margin = margin(8, 8, 8, 8), plot.background = element_rect(fill = "white", colour = NA))

# =============================================================================
# Save
# =============================================================================
cat("Output directory:\n  ", fig_dir, "\n", sep = "")
save_panel(p_count,    "01_active_gus_count_distribution", 8.6, 5.4)
save_panel(p_stack,    "04_stacked_relative_abundance",    18, 7.2)
save_panel(p_hm_nocl,  "05_abundance_heatmap_tpm_nocluster", 15.2, 5.8)
save_panel(p_scheme,   "00_analysis_scheme",             11.4, 8.6)

# Figure legends (English, paper style)
prev <- grid %>%
  group_by(gus_label) %>%
  summarise(n = sum(present), .groups = "drop") %>%
  arrange(desc(n), gus_label)
prev_txt <- paste0(prev$gus_label, " (n = ", prev$n, ")", collapse = "; ")
dist_txt <- count_dist %>%
  filter(n_people > 0) %>%
  mutate(bit = sprintf("%s (%d people, %.1f%%)", n_active, n_people, pct_people)) %>%
  pull(bit) %>%
  paste(collapse = "; ")

legends <- paste0(
  "# Figure legends (paper figures v2, assembly pipeline)\n\n",
  "Source tables: `results/hmp_assembly_gus/report/`  \n",
  "Figures: `results/hmp_assembly_gus/report/figures_v2/`  \n",
  "Drawn with ggplot2. Font: Arial. Plot titles/subtitles removed from data panels.\n\n",
  "---\n\n",
  "**Figure. Distribution of active GUS richness in the HMP stool cohort.**  \n",
  "Number of distinct active β-glucuronidase (GUS) enzymes detected per person ",
  "in ", n_people, " HMP HMASM stool shotgun metagenomes (MEGAHIT assembly; ",
  "identity ≥ 95%, coverage ≥ 90%). Presence was called when a predicted protein ",
  "matched one of ", n_refs, " experimentally supported active GUS references. ",
  "Numbers above bars are headcounts. Detected richness ranged from 0 to ",
  n_max, " enzymes per person (", dist_txt, "). ", n_zero,
  " people had no active GUS and ", n_any, " people had at least one.\n\n",
  "---\n\n",
  "**Figure. Genome-wide abundance of 11 active GUS enzymes across ", n_people, " people.**  \n",
  "Heatmap of genome-wide transcripts per million (TPM) for ", n_refs,
  " active GUS references (rows) in ", n_people, " HMP stool metagenomes (columns). ",
  "Columns are ordered by decreasing number of distinct active GUS detected, ",
  "then by total genome-wide TPM. Rows follow clade then source. Color shows ",
  "log10 genome-wide TPM for detected genes (plasma scale); light grey ",
  "indicates absence (TPM = 0). Prevalence: ", prev_txt, ".\n\n",
  "---\n\n",
  "**Figure. Relative composition of active GUS among carriers.**  \n",
  "Stacked bars show the relative abundance of the ", n_refs,
  " active GUS references within each person who carried at least one enzyme ",
  "(n = ", n_any, "). For each person, genome-wide TPM values of the ", n_refs,
  " references were scaled to 100%. Each bar is one sample (SRS ID). ",
  "Individuals are ordered by the relative abundance of *B. dorei* (NL) ",
  "from low to high (people without *B. dorei* appear first). Colors denote ",
  "the source organism and GUS clade (FMN, L1, L2, NL, hybrid, mL1).\n\n",
  "---\n\n",
  "**Figure. Analysis scheme for HMP active GUS quantification.**  \n",
  "HMP HMASM stool shotgun metagenomes (n = 139; paired trimmed FASTQ) were ",
  "assembled with MEGAHIT (contigs ≥ 1,000 bp). Predicted proteins (Prodigal, ",
  "meta mode) were searched with DIAMOND blastp against 11 active GUS references ",
  "(identity ≥ 95%, query and subject coverage ≥ 90%). The same reads were mapped ",
  "back to each sample’s assembly (Bowtie2) and counted on predicted CDS ",
  "(featureCounts). Best hits were combined with CDS counts to obtain genome-wide ",
  "TPM, which was used for presence, abundance, and composition analyses of active GUS.\n"
)
writeLines(legends, file.path(fig_dir, "figure_legends.md"), useBytes = TRUE)
cat("  figure_legends.md\nDone.\n")
