# HMP active GUS and RPKM quantification

## Active GUS

This analysis quantifies active β-glucuronidase (GUS) genes in 139 HMP stool shotgun metagenomes.

**Active GUS** is the first 11 rows of Table S3 from the preprint (Zhang et al., 2026).

Active GUS IDs (Table S3 order):

`NP_416134.1` (*E. coli*), `WP_003467686.1` (*C. perfringens*), `6MVG_A` (*R. gnavus*), `WP_022501306` (*Eubacterium*), `WP_000966715` (*S. agalactiae*), `WP_012740861` (*E. eligens*), `3CMG_A` (*B. fragilis*), `6NCY_A` (*F. saccharivorans*), `WP_270534885` (*B. mass*), `WP_007841259` (*B. dorei*), `WP_217761416` (*P. copri*).

Predicted genes are assigned to their best active GUS reference hit using identity ≥ 95% and query/subject coverage ≥ 90%.

## RPKM quantification

Gene abundance is normalized by gene length and the total counts assigned to all genes in the sample:

```text
RPKM = raw_count × 10^9 / (gene_length_bp × library_size)
```

`library_size` is the sum of all gene counts in the sample's featureCounts output, not just active GUS counts. RPKM values for genes assigned to the same active GUS reference are summed within each sample.

Outputs under `results/hmp_active_gus/summary/`:

- `active_gus_rpkm_genome_wide_matrix.tsv`: active GUS reference × sample RPKM matrix.
- `active_gus_abundance_all_samples.tsv`: gene-level abundance table including `rpkm_genome_wide`.

Recalculate and validate RPKM from existing per-sample abundance tables and featureCounts outputs:

```bash
mkdir -p logs
sbatch scripts/submit_hmp_active_gus_rpkm.sbatch
```

The recalculation preserves the sample and reference order of the existing abundance matrix and checks saved counts against featureCounts outputs.
