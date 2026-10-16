# Active GUS abundance in HMP stool metagenomes

This project quantifies genes matching the **active β-glucuronidase (GUS) reference set** in 139 Human Microbiome Project (HMP) stool shotgun metagenomes. Abundance is reported as **RPKM**, normalized by gene length and the total number of reads mapped to each sample-specific assembly.

The working analysis requires **≥95% amino-acid identity and ≥90% query and reference coverage** (`pident95_cov90`). A match identifies a homolog of an experimentally active reference protein.

## Samples and input data

- The cohort comprises the 139 HMP stool metagenomic samples used to establish the GUS Atlas, with SRS sample IDs listed in `data/HMP_GUS/manifests/HMASM_WGS_stool_SRS_samples.txt`.
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
  → calculate gene-level RPKM using total reads mapped to the sample-specific assembly
  → sum RPKM by active GUS reference within each sample
```

HMP-preprocessed paired-end reads were assembled separately for each sample using MEGAHIT, retaining contigs of at least 1,000 bp. Prodigal predicted coding sequences in metagenomic mode. Reads were mapped back to the corresponding sample-specific assembly using Bowtie 2, sorted using SAMtools, and assigned to coding sequences using featureCounts.

Predicted protein sequences were screened against the active GUS references using DIAMOND BLASTp. Qualifying hits required ≥95% amino-acid identity, ≥90% query and subject coverage, and an E-value ≤10⁻⁵. The highest-bit-score qualifying match was retained per gene, giving each gene one active GUS reference assignment.

The assembly, mapping, coding-sequence counting, and active GUS screening workflow is defined in `workflow/hmp_assembly_gus.smk`, with settings in `config/hmp_assembly_gus.yaml`.

## RPKM quantification

For GUS gene *g* in sample *s*:

```text
RPKM(g, s) = C(g, s) × 10^9 / (L(g) × N(s))
```

| Symbol | Definition |
| --- | --- |
| `C(g, s)` | Number of reads assigned to the GUS gene by featureCounts |
| `L(g)` | Coding-sequence length of the GUS gene in base pairs |
| `N(s)` | Total number of reads mapped to the corresponding sample-specific assembly |

The denominator uses all reads mapped to the sample-specific assembly, including mapped reads that are not assigned to a coding sequence. It must therefore be obtained from the assembly read mappings rather than calculated by summing gene counts.

Gene-level RPKM values are summed for genes assigned to the same active GUS reference within each sample to obtain reference-level abundance.

## Analysis files

Assembly-based analysis files are stored under `results/hmp_assembly_gus/`:

| Relative path | Contents |
| --- | --- |
| `mapping/{SRS}.min1000.sorted.bam` | Reads mapped to the sample-specific assembly, sorted by coordinate |
| `prodigal/{SRS}.min1000.gff` | Predicted coding-sequence coordinates |
| `prodigal/{SRS}.min1000.faa` | Predicted protein sequences used for screening |
| `counts/{SRS}.min1000.featurecounts.txt` | Reads assigned to predicted coding sequences and their lengths |
| `summary/pident95_cov90/hits/{SRS}.min1000.active_gus_best_hits.tsv` | Highest-bit-score qualifying active GUS match per gene |
| `summary/pident95_cov90/abundance/{SRS}.min1000.active_gus_abundance.tsv` | Gene-level active GUS abundance output from the assembly workflow |

The mapped-read total for RPKM normalization comes from the corresponding assembly mapping. Gene counts and coding-sequence lengths come from the featureCounts output.

## Repository layout

```text
config/       Sample selection and assembly workflow configuration
workflow/     Assembly, mapping, coding-sequence counting, and GUS screening
scripts/      Reference setup, job submission, hit filtering, and reporting
data/         Local HMP read packages and sample manifests
results/      Local assemblies, predicted genes, alignments, counts, and abundance tables
logs/         Slurm job logs
```

Large input and result files are stored locally and are excluded from Git. The assembly workflow uses reference files and shared abundance scripts from the companion `2026-ysjin-gus` checkout, configured through `gus_root` and the reference paths in `config/hmp_assembly_gus.yaml`.
