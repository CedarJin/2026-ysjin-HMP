# HMP stool WGS GUS abundance

This repository quantifies β-glucuronidase (GUS) genes in 139 HMP stool shotgun metagenomes (HMASM WGS). HMP-preprocessed reads are assembled, genes are predicted, reads are mapped back to each sample’s own assembly, predicted proteins are searched against paper GUS reference sets, and abundances are reported as **genome-wide TPM**.

The working analysis uses identity ≥ 95% and query/subject coverage ≥ 90% (`pident95_cov90`). GUS databases and abundance scripts are shared with [`2026-ysjin-gus`](https://github.com/CedarJin/2026-ysjin-gus).

## Dataset

| Item | Value |
| --- | --- |
| Source | HMP HMASM stool WGS |
| Samples | 139 SRS IDs in `data/HMP_GUS/manifests/HMASM_WGS_stool_SRS_samples.txt` (same cohort as Pollet et al., 2017) |
| Reads | HMP preprocessed paired FASTQ (`denovo_duplicates_marked.trimmed`) |

Sample `SRS023176` uses a repaired read pair; the public HMASM tarball has broken pairing.

## GUS reference sets

| Set | Source | Typical size |
| --- | --- | --- |
| **Total GUS** | 279 GUS sequences from Pollet et al., 2017. 129 FMN-binding GUS sequences from Simpson et al., 2022 | 408 |
| **Candidate GUS** | Table S3 (18 experimentally discussed proteins) from the preprint (Zhang et al., 2026) | 18 |
| **Active GUS** | First 11 rows of Table S3 from the preprint (Zhang et al., 2026) | 11 |

Active GUS IDs (Table S3 order):

`NP_416134.1` (*E. coli*), `WP_003467686.1` (*C. perfringens*), `6MVG_A` (*R. gnavus*), `WP_022501306` (*Eubacterium*), `WP_000966715` (*S. agalactiae*), `WP_012740861` (*E. eligens*), `3CMG_A` (*B. fragilis*), `6NCY_A` (*F. saccharivorans*), `WP_270534885` (*B. mass*), `WP_007841259` (*B. dorei*), `WP_217761416` (*P. copri*).

## Pipeline

```text
HMP HMASM trimmed FASTQ
  → MEGAHIT assembly (--min-contig-len 1000)
  → Prodigal (meta) on final.contigs.fa
  → Bowtie2: same reads → same assembly
  → featureCounts on Prodigal CDS
  → DIAMOND blastp vs total, candidate, and active GUS
  → filter hits → best hit per predicted gene
  → genome-wide TPM and GUS-set TPM
```

- Workflow: `workflow/hmp_assembly_gus.smk`
- Config: `config/hmp_assembly_gus.yaml`
- Submit all samples: `sbatch scripts/submit_hmp_assembly_gus_array.sbatch`
- Review array: `bash scripts/review_hmp_assembly_gus_array.sh`

| Profile | Identity | Query cov. | Subject cov. | Output directory |
| --- | --- | --- | --- | --- |
| total / candidate | ≥ 40% | ≥ 50% | ≥ 50% | `results/hmp_assembly_gus/summary/` |
| **active (`pident95_cov90`)** | **≥ 95%** | **≥ 90%** | **≥ 90%** | **`results/hmp_assembly_gus/summary/pident95_cov90/`** |

## Key outputs

Per sample, under `results/hmp_assembly_gus/`:

- `summary/{SRS}.min1000.total_gus_abundance.tsv`
- `summary/{SRS}.min1000.candidate_gus_abundance.tsv`
- `summary/pident95_cov90/abundance/{SRS}.min1000.active_gus_abundance.tsv`

Each row is one predicted gene (best GUS hit) plus `raw_count`, `cds_length_bp`, `rpk`, `cpm`, `tpm_gus_db`, and `tpm_genome_wide`.

Cohort report: `results/hmp_assembly_gus/report/` (`scripts/report_hmp_assembly_gus.R`).

## How to run

```bash
bash scripts/setup_hmp_assembly_gus_refs.sh
sbatch scripts/submit_hmp_assembly_gus_array.sbatch
bash scripts/review_hmp_assembly_gus_array.sh
Rscript scripts/report_hmp_assembly_gus.R .
Rscript scripts/plot_hmp_active_gus_heatmap_boxplot.R .
```

Pilot the first five samples:

```bash
sbatch --array=1-5%5 scripts/submit_hmp_assembly_gus_array.sbatch
```

The abundance workflow expects HMP trimmed pairs at the paths in `config/hmp_assembly_gus.yaml` (`hmp_reads_dir`; repaired reads for `SRS023176`). GUS `.faa` files and `summarize_gus_hits.py` / `calculate_gus_abundance.py` come from the companion `2026-ysjin-gus` checkout (`gus_root`).

## Repository layout

```text
config/       sample list and hmp_assembly_gus.yaml
workflow/     assembly + GUS abundance Snakemake
scripts/      sbatch, setup, review, and R report
data/         HMASM packages and manifests (gitignored)
results/      assemblies, searches, abundance, report (gitignored)
```

## References

Pollet RM, D'Agostino EH, Walton WG, Xu Y, Little MS, Biernat KA, Pellock SJ, Patterson LM, Creekmore BC, Isenberg HN, Bahethi RR, Bhatt AP, Liu J, Gharaibeh RZ, Redinbo MR. An Atlas of β-Glucuronidases in the Human Intestinal Microbiome. Structure. 2017 Jul 5;25(7):967-977.e5. doi: 10.1016/j.str.2017.05.003.

Simpson, J. B., Sekela, J. J., Graboski, A. L., Borlandelli, V. B., Bivins, M. M., Barker, N. K., … Redinbo, Matthew. R. (2022). Metagenomics combined with activity-based proteomics point to gut bacterial enzymes that reactivate mycophenolate. Gut Microbes, 14(1). https://doi.org/10.1080/19490976.2022.2107289

Guodong Zhang, Nan Jing, Guodong Cao et al. Gut Microbiota Mediate the Metabolism of Colonic Prostaglandins, 20 February 2026, PREPRINT (Version 1) available at Research Square [https://doi.org/10.21203/rs.3.rs-8856024/v1]
