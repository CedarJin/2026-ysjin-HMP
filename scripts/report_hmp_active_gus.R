#!/usr/bin/env Rscript
# HMP active GUS presence / abundance report for 139 stool WGS samples.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(scales)
  library(stringr)
  library(forcats)
  library(RColorBrewer)
  library(viridisLite)
  library(base64enc)
})

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args) >= 1) args[[1]] else getwd()
project_root <- normalizePath(project_root, mustWork = TRUE)

sample_list_path <- file.path(project_root, "data/HMP_GUS/manifests/HMASM_WGS_stool_SRS_samples.txt")
meta_path <- file.path(project_root, "results/reference/active_gus_metadata.tsv")
abundance_path <- file.path(project_root, "results/hmp_active_gus/summary/active_gus_abundance_all_samples.tsv")
out_dir <- file.path(project_root, "results/hmp_active_gus/report")
fig_dir <- file.path(out_dir, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

samples <- read_lines(sample_list_path) %>% str_trim() %>% .[. != ""]
meta <- read_tsv(meta_path, show_col_types = FALSE) %>%
  mutate(
    gus_label = sprintf("%s (%s)", source, clade),
    reference_id = as.character(reference_id)
  )

# Stable row order: clade then source.
meta <- meta %>%
  arrange(clade, source, reference_id) %>%
  mutate(gus_label = factor(gus_label, levels = gus_label))

abundance <- read_tsv(abundance_path, show_col_types = FALSE) %>%
  mutate(
    sample = as.character(sample),
    sseqid = as.character(sseqid),
    tpm_genome_wide = as.numeric(tpm_genome_wide),
    tpm_active_gus_db = as.numeric(tpm_active_gus_db),
    raw_count = as.numeric(raw_count)
  )

# One value per sample × reference (sum if multiple gene hits map to same ref).
hit_summary <- abundance %>%
  group_by(sample, sseqid) %>%
  summarise(
    raw_count = sum(raw_count, na.rm = TRUE),
    tpm_genome_wide = sum(tpm_genome_wide, na.rm = TRUE),
    tpm_active_gus_db = sum(tpm_active_gus_db, na.rm = TRUE),
    .groups = "drop"
  )

grid <- tidyr::expand_grid(
  sample = samples,
  reference_id = meta$reference_id
) %>%
  left_join(meta %>% select(reference_id, gus_label, clade, source), by = "reference_id") %>%
  left_join(hit_summary, by = c("sample" = "sample", "reference_id" = "sseqid")) %>%
  mutate(
    raw_count = replace_na(raw_count, 0),
    tpm_genome_wide = replace_na(tpm_genome_wide, 0),
    tpm_active_gus_db = replace_na(tpm_active_gus_db, 0),
    present = as.integer(raw_count > 0 | tpm_genome_wide > 0)
  )

# Sample ordering: more active GUS first, then total genome-wide TPM.
sample_order_df <- grid %>%
  group_by(sample) %>%
  summarise(
    n_active = sum(present),
    total_tpm = sum(tpm_genome_wide),
    .groups = "drop"
  ) %>%
  arrange(desc(n_active), desc(total_tpm), sample) %>%
  mutate(sample_f = factor(sample, levels = sample))

grid <- grid %>%
  mutate(sample_f = factor(sample, levels = levels(sample_order_df$sample_f)))

# Relative abundance among the 11 active GUS references (sums to 100% when any present).
rel <- grid %>%
  group_by(sample) %>%
  mutate(
    total_tpm = sum(tpm_genome_wide),
    pct = dplyr::if_else(
      total_tpm > 0,
      100 * tpm_genome_wide / total_tpm,
      0
    )
  ) %>%
  mutate(
    # Guard against floating-point drift so stacked bars stay within [0, 100].
    pct = dplyr::if_else(
      total_tpm > 0 & sum(pct) > 0,
      100 * pct / sum(pct),
      0
    )
  ) %>%
  ungroup() %>%
  mutate(sample_f = factor(sample, levels = levels(sample_order_df$sample_f)))

# ---- Count distribution ----------------------------------------------------
count_dist <- sample_order_df %>%
  count(n_active, name = "n_people") %>%
  complete(n_active = 0:nrow(meta), fill = list(n_people = 0)) %>%
  mutate(
    pct_people = 100 * n_people / length(samples),
    n_active_f = factor(n_active, levels = 0:nrow(meta))
  )

write_tsv(count_dist, file.path(out_dir, "active_gus_count_distribution.tsv"))
write_tsv(
  sample_order_df %>% select(sample, n_active, total_tpm),
  file.path(out_dir, "active_gus_per_sample.tsv")
)
write_tsv(
  grid %>%
    select(sample, reference_id, gus_label, clade, source, present, tpm_genome_wide, tpm_active_gus_db, raw_count),
  file.path(out_dir, "active_gus_sample_by_reference.tsv")
)

p_count <- ggplot(count_dist, aes(x = n_active_f, y = n_people)) +
  geom_col(fill = "#2F6F8F", width = 0.72) +
  geom_text(aes(label = n_people), vjust = -0.35, size = 3.4, color = "#222222") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(
    title = "Number of people with N active GUS",
    subtitle = sprintf("HMP stool WGS cohort (n = %d); reference panel = %d active GUS", length(samples), nrow(meta)),
    x = "Number of distinct active GUS detected",
    y = "Number of people"
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(color = "#555555")
  )

# ---- Shared color / theme helpers -----------------------------------------
gus_colors <- setNames(
  colorRampPalette(brewer.pal(12, "Set3"))(nrow(meta)),
  levels(meta$gus_label)
)
# Make labels more distinguishable: override with a curated distinct set.
distinct_cols <- c(
  "#4E79A7", "#F28E2B", "#E15759", "#76B7B2", "#59A14F",
  "#EDC948", "#B07AA1", "#FF9DA7", "#9C755F", "#BAB0AC", "#86BCB6"
)
gus_colors <- setNames(distinct_cols[seq_len(nrow(meta))], levels(meta$gus_label))

heatmap_theme <- theme_minimal(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 5.5),
    axis.text.y = element_text(size = 9),
    legend.position = "top",
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(color = "#555555"),
    plot.margin = margin(8, 12, 8, 8)
  )

# ---- Binary presence heatmap (circular cells) ------------------------------
p_bin <- ggplot(grid, aes(x = sample_f, y = gus_label)) +
  geom_point(
    aes(color = factor(present, levels = c(0, 1), labels = c("Absent", "Present"))),
    size = 2.1, shape = 16
  ) +
  scale_color_manual(values = c("Absent" = "#D9D9D9", "Present" = "#1F4E79"), name = NULL) +
  labs(
    title = "Active GUS presence (binary)",
    subtitle = "Columns = people (ordered by # active GUS); rows = active GUS reference; filled = present",
    x = NULL,
    y = NULL
  ) +
  heatmap_theme

# ---- TPM abundance heatmap (circular cells) --------------------------------
tpm_nonzero <- grid$tpm_genome_wide[grid$tpm_genome_wide > 0]
tpm_max <- if (length(tpm_nonzero)) quantile(tpm_nonzero, 0.98) else 1
grid_plot <- grid %>%
  mutate(tpm_plot = pmin(tpm_genome_wide, tpm_max))

p_tpm <- ggplot(grid_plot, aes(x = sample_f, y = gus_label)) +
  geom_point(
    aes(color = tpm_plot, alpha = present == 1),
    size = 2.1, shape = 16
  ) +
  scale_color_viridis_c(
    option = "C",
    name = "Genome-wide TPM",
    limits = c(0, tpm_max),
    breaks = pretty_breaks(4)
  ) +
  scale_alpha_manual(values = c("FALSE" = 0.18, "TRUE" = 1), guide = "none") +
  labs(
    title = "Active GUS abundance (genome-wide TPM)",
    subtitle = sprintf(
      "Color = genome-wide TPM (color scale capped at 98th percentile of non-zero values = %.2f)",
      tpm_max
    ),
    x = NULL,
    y = NULL
  ) +
  heatmap_theme

# ---- Stacked relative abundance --------------------------------------------
# Drop people with 0 active GUS; order x-axis by B.dorei (NL) relative % low → high.
bdorei_label <- "B.dorei (NL)"
stack_samples <- sample_order_df %>%
  filter(n_active > 0) %>%
  select(sample)

stack_order <- rel %>%
  semi_join(stack_samples, by = "sample") %>%
  mutate(
    bdorei_pct = if_else(as.character(gus_label) == bdorei_label, pct, 0)
  ) %>%
  group_by(sample) %>%
  summarise(bdorei_pct = sum(bdorei_pct), .groups = "drop") %>%
  arrange(bdorei_pct, sample)

stack_all <- rel %>%
  semi_join(stack_samples, by = "sample") %>%
  mutate(
    gus_label = factor(as.character(gus_label), levels = levels(meta$gus_label)),
    sample_f = factor(sample, levels = stack_order$sample)
  )

p_stack <- ggplot(stack_all, aes(x = sample_f, y = pct, fill = gus_label)) +
  geom_col(width = 0.92, color = NA, position = "stack") +
  scale_fill_manual(values = gus_colors, name = "Active GUS", drop = FALSE) +
  scale_y_continuous(breaks = seq(0, 100, 25), expand = c(0, 0)) +
  coord_cartesian(ylim = c(0, 100), clip = "off") +
  labs(
    title = "Relative abundance among detected active GUS",
    subtitle = sprintf(
      "People with >=1 active GUS (n = %d); x ordered by %s relative %% low -> high; y = %% among the 11 active GUS (=100%%)",
      nrow(stack_order), bdorei_label
    ),
    x = NULL,
    y = "Relative abundance (%)"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 5.5),
    legend.position = "right",
    legend.text = element_text(size = 9),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(color = "#555555", size = 9.5),
    plot.margin = margin(8, 8, 8, 8)
  )

# ---- Save figures ----------------------------------------------------------
fig_count <- file.path(fig_dir, "01_active_gus_count_distribution.png")
fig_bin <- file.path(fig_dir, "02_presence_heatmap_binary.png")
fig_tpm <- file.path(fig_dir, "03_abundance_heatmap_tpm.png")
fig_stack <- file.path(fig_dir, "04_stacked_relative_abundance.png")

ggsave(fig_count, p_count, width = 8.5, height = 5.2, dpi = 160, bg = "white")
ggsave(fig_bin, p_bin, width = 18, height = 6.8, dpi = 160, bg = "white")
ggsave(fig_tpm, p_tpm, width = 18, height = 6.8, dpi = 160, bg = "white")
ggsave(fig_stack, p_stack, width = 18, height = 7.2, dpi = 160, bg = "white")

# Also PDF copies for publication-quality viewing.
ggsave(sub("\\.png$", ".pdf", fig_count), p_count, width = 8.5, height = 5.2, bg = "white")
ggsave(sub("\\.png$", ".pdf", fig_bin), p_bin, width = 18, height = 6.8, bg = "white")
ggsave(sub("\\.png$", ".pdf", fig_tpm), p_tpm, width = 18, height = 6.8, bg = "white")
ggsave(sub("\\.png$", ".pdf", fig_stack), p_stack, width = 18, height = 7.2, bg = "white")

# ---- HTML report -----------------------------------------------------------
escape <- function(x) {
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x
}

img_tag <- function(path, alt) {
  b64 <- base64encode(path)
  sprintf(
    '<img src="data:image/png;base64,%s" alt="%s" style="max-width:100%%;height:auto;border:1px solid #e5e5e5;border-radius:4px;"/>',
    b64, escape(alt)
  )
}

count_rows <- apply(count_dist, 1, function(row) {
  sprintf(
    "<tr><td>%s</td><td>%s</td><td>%.1f%%</td></tr>",
    escape(row[["n_active"]]),
    escape(row[["n_people"]]),
    as.numeric(row[["pct_people"]])
  )
})

ref_rows <- apply(meta, 1, function(row) {
  sprintf(
    "<tr><td>%s</td><td>%s</td><td>%s</td><td>%s</td></tr>",
    escape(row[["reference_id"]]),
    escape(row[["clade"]]),
    escape(row[["source"]]),
    escape(row[["gus_label"]])
  )
})

n_with_any <- sum(sample_order_df$n_active > 0)
n_with_zero <- sum(sample_order_df$n_active == 0)
refs_detected <- grid %>%
  group_by(reference_id, gus_label) %>%
  summarise(n_people = sum(present), .groups = "drop") %>%
  arrange(desc(n_people), gus_label)

ref_prev_rows <- apply(refs_detected, 1, function(row) {
  sprintf(
    "<tr><td>%s</td><td>%s</td><td>%s</td><td>%.1f%%</td></tr>",
    escape(row[["gus_label"]]),
    escape(row[["reference_id"]]),
    escape(row[["n_people"]]),
    100 * as.numeric(row[["n_people"]]) / length(samples)
  )
})

html <- paste0(
  "<!DOCTYPE html>\n<html lang=\"en\">\n<head>\n<meta charset=\"utf-8\"/>\n",
  "<title>HMP Active GUS Report</title>\n",
  "<style>\n",
  "body{font-family:Helvetica,Arial,sans-serif;margin:32px;color:#222;line-height:1.45;background:#fafafa;}\n",
  "main{max-width:1200px;margin:0 auto;background:#fff;padding:28px 32px;border:1px solid #e8e8e8;border-radius:8px;}\n",
  "h1{font-size:1.7rem;margin:0 0 0.4rem 0;} h2{margin-top:2rem;border-bottom:1px solid #eee;padding-bottom:0.35rem;}\n",
  "p,li{font-size:0.98rem;} .muted{color:#666;} table{border-collapse:collapse;margin:0.8rem 0 1.2rem;font-size:0.92rem;}\n",
  "th,td{border:1px solid #ddd;padding:6px 10px;text-align:left;} th{background:#f3f6f8;} .fig{margin:1rem 0 1.6rem;}\n",
  ".kpi{display:flex;gap:12px;flex-wrap:wrap;margin:1rem 0;} .kpi div{background:#f3f6f8;padding:10px 14px;border-radius:6px;min-width:140px;}\n",
  ".kpi strong{display:block;font-size:1.35rem;}\n",
  "</style>\n</head>\n<body><main>\n",
  "<h1>HMP Active GUS Report</h1>\n",
  "<p class=\"muted\">Generated ", format(Sys.time(), "%Y-%m-%d %H:%M %Z"),
  " &middot; cohort = HMASM stool WGS SRS samples</p>\n",
  "<div class=\"kpi\">",
  sprintf("<div><strong>%d</strong>people (samples)</div>", length(samples)),
  sprintf("<div><strong>%d</strong>active GUS in reference panel</div>", nrow(meta)),
  sprintf("<div><strong>%d</strong>people with &ge;1 active GUS</div>", n_with_any),
  sprintf("<div><strong>%d</strong>people with 0 active GUS</div>", n_with_zero),
  "</div>\n",
  "<h2>1. How many active GUS per person?</h2>\n",
  "<p>Presence is defined as a best-hit gene mapping to an active GUS reference ",
  "(identity / coverage / e-value thresholds from the pipeline) with non-zero counts.</p>\n",
  "<table><thead><tr><th># active GUS</th><th># people</th><th>% of cohort</th></tr></thead><tbody>\n",
  paste(count_rows, collapse = "\n"),
  "\n</tbody></table>\n",
  "<div class=\"fig\">", img_tag(fig_count, "Count distribution"), "</div>\n",
  "<h2>2. Prevalence of each active GUS</h2>\n",
  "<table><thead><tr><th>Active GUS</th><th>Reference ID</th><th># people</th><th>% of cohort</th></tr></thead><tbody>\n",
  paste(ref_prev_rows, collapse = "\n"),
  "\n</tbody></table>\n",
  "<h2>3. Binary presence heatmap</h2>\n",
  "<p>Each column is one person; each row is one active GUS. Circles are colored by presence/absence.</p>\n",
  "<div class=\"fig\">", img_tag(fig_bin, "Binary presence heatmap"), "</div>\n",
  "<h2>4. TPM abundance heatmap</h2>\n",
  "<p>Same layout as above; circle color encodes genome-wide TPM abundance.</p>\n",
  "<div class=\"fig\">", img_tag(fig_tpm, "TPM abundance heatmap"), "</div>\n",
  "<h2>5. Stacked relative abundance</h2>\n",
  "<p>Each bar is one person with ≥1 active GUS (n = ", nrow(stack_order),
  "). Colors are active GUS. Y-axis is relative abundance (%) among the ",
  nrow(meta), " active GUS references (denominator = sum of their genome-wide TPM = 100%). ",
  "X-axis is ordered by <strong>B.dorei (NL)</strong> relative % from low to high ",
  "(people without B.dorei appear first).</p>\n",
  "<div class=\"fig\">", img_tag(fig_stack, "Stacked relative abundance"), "</div>\n",
  "<h2>Reference panel</h2>\n",
  "<table><thead><tr><th>Reference ID</th><th>Clade</th><th>Source</th><th>Plot label</th></tr></thead><tbody>\n",
  paste(ref_rows, collapse = "\n"),
  "\n</tbody></table>\n",
  "<h2>Methods notes</h2>\n<ul>\n",
  "<li>Input matrix rebuilt from <code>active_gus_abundance_all_samples.tsv</code> plus the full ",
  nrow(meta), "-member active GUS metadata, so references never observed in the cohort still appear as empty rows.</li>\n",
  "<li>Heatmaps: samples ordered by decreasing number of distinct active GUS, then by total genome-wide TPM.</li>\n",
  "<li>Stacked plot: excludes people with 0 active GUS; x ordered by B.dorei (NL) relative % low → high.</li>\n",
  "<li>Relative stacked percentages use genome-wide TPM summed within each person across the active GUS panel.</li>\n",
  "<li>Figure files are also saved under <code>results/hmp_active_gus/report/figures/</code>.</li>\n",
  "</ul>\n",
  "</main></body></html>\n"
)

html_path <- file.path(out_dir, "hmp_active_gus_report.html")
writeLines(html, html_path, useBytes = TRUE)

# Convenience copy under reports/
reports_dir <- file.path(project_root, "reports")
dir.create(reports_dir, showWarnings = FALSE, recursive = TRUE)
file.copy(html_path, file.path(reports_dir, "hmp_active_gus_report.html"), overwrite = TRUE)

cat("Wrote:\n")
cat(" ", html_path, "\n")
cat(" ", file.path(reports_dir, "hmp_active_gus_report.html"), "\n")
cat(" Figures in ", fig_dir, "\n", sep = "")
print(as.data.frame(count_dist %>% select(n_active, n_people, pct_people)))
