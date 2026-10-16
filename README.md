# Active GUS abundance in HMP stool metagenomes

This project quantifies genes matching the **active β-glucuronidase (GUS) reference set** in 139 Human Microbiome Project (HMP) stool shotgun metagenomes. Abundance is reported as **RPKM**, normalized by gene length and the total number of reads mapped to each sample-specific assembly.

The working analysis requires **≥95% amino-acid identity and ≥90% query and reference coverage** (`pident95_cov90`). A match identifies a homolog of an experimentally active reference protein.

## Samples and input data

- The cohort comprises the 139 HMP stool metagenomic samples, with SRS sample IDs listed in `data/HMP_GUS/manifests/HMASM_WGS_stool_SRS_samples.txt`.
- Paired-end reads come from HMP HMASM stool WGS packages, using the HMP-preprocessed `denovo_duplicates_marked.trimmed` FASTQ pairs.
- Coding sequences are predicted from sample-specific assemblies generated from these reads.
- `SRS023176` requires repaired or rebuilt read pairs because its public HMASM package has broken pairing. The local configuration routes this sample to the prepared replacement reads.


## Active GUS references

The active set comprises the first 11 entries of Table S3 in Zhang et al. (2026), *Gut Microbiota Mediate the Metabolism of Colonic Prostaglandins*. Reference identifiers are retained as supplied in the analysis, following the companion `2026-ysjin-gus` project.

| Reference ID | Source label |
| --- | --- |
| `NP_416134.1` | E. coli |
| `WP_003467686.1` | C. perfringens |
| `6MVG_A` | R. gnavus |
| `WP_022501306` | Eubacterium |
| `WP_000966715` | S. agalactiae |
| `WP_012740861` | E. eligens |
| `3CMG_A` | B. fragilis |
| `6NCY_A` | F. saccharivorans |
| `WP_270534885` | B. mass |
| `WP_007841259` | B. dorei |
| `WP_217761416` | P. copri |

Reference sequences and annotations are supplied by the companion `2026-ysjin-gus` checkout through the `active_faa` and `active_metadata` paths in `config/hmp_assembly_gus.yaml`.

## Analysis workflow

```text
HMP-preprocessed paired-end reads (139 stool metagenomes)
  → MEGAHIT v1.2.9 assembly (retain contigs ≥1,000 bp)
  → Prodigal v2.6.3 coding-sequence prediction (metagenomic mode)
  → Bowtie 2 v2.5.5 mapping of reads back to sample-specific assemblies
  → SAMtools v1.23.1 sorting of alignments
  → featureCounts v2.1.1 assignment of reads to coding sequences
  → DIAMOND BLASTp v2.1.11 screening of predicted proteins against active GUS references
  → filter hits (identity ≥95%, query and subject coverage ≥90%, E-value ≤1e-5)
  → retain the highest-bit-score qualifying match per gene
  → Abundance calculation (RPKM)
```

The assembly, mapping, coding-sequence counting, and active GUS screening workflow is defined in `workflow/hmp_assembly_gus.smk`, with settings in `config/hmp_assembly_gus.yaml`.


## Repository layout

```text
config/       Sample selection and assembly workflow configuration
workflow/     Assembly, mapping, coding-sequence counting, and GUS screening
scripts/      Reference setup, job submission, hit filtering, and reporting
data/         Local HMP read packages and sample manifests
results/      Local assemblies, predicted genes, alignments, counts, and abundance tables
logs/         Slurm job logs
```

Large input and result files are stored locally and are excluded from Git. 
