import os

from snakemake.utils import min_version

min_version("7.0")

configfile: "config/hmp_active_gus.yaml"

TOOL_BIN_DIR = config.get("tool_bin_dir", "")
if TOOL_BIN_DIR:
    shell.prefix(f"set -euo pipefail; export PATH={TOOL_BIN_DIR}:$PATH; ")
else:
    shell.prefix("set -euo pipefail; ")

RESULTS = config["results_dir"]
SAMPLE_LIST = config["sample_list"]
HMGI_ARCHIVE = config["hmgI_gene_archive"]
READ_DIR = config["hmasm_read_dir"]
REF_FAA = config["reference_faa"]
REF_META = config["reference_metadata"]
THREADS = int(config.get("threads", 8))
THREADS_BOWTIE2 = int(config.get("threads_bowtie2", THREADS))
THREADS_FEATURECOUNTS = int(config.get("threads_featurecounts", THREADS))
THREADS_DIAMOND = int(config.get("threads_diamond", THREADS))
READ_PAIR_LIMIT = int(config.get("read_pair_limit", 0))
REPAIR_SAMPLES = set(config.get("repair_samples", []))
SRA_REBUILD = config.get("sra_rebuild", {})
SRA_SAMPLE = SRA_REBUILD.get("sample", "")
SRA_RAW_DIR = SRA_REBUILD.get("raw_dir", "data/HMP_GUS/SRA")
SRA_RUN_COUNTS = SRA_REBUILD.get("runs", {})
SRA_RUNS = list(SRA_RUN_COUNTS)
SRA_RUN_PATTERN = "|".join(SRA_RUNS) if SRA_RUNS else r"a^"


def load_samples(path):
    with open(path) as handle:
        return [line.strip() for line in handle if line.strip()]


SAMPLES = load_samples(SAMPLE_LIST)


def read1_for_mapping(sample):
    if sample in REPAIR_SAMPLES:
        return f"{RESULTS}/repaired_reads/{sample}/{sample}.denovo_duplicates_marked.trimmed.repaired.1.fastq"
    return f"{RESULTS}/reads/{sample}/{sample}.denovo_duplicates_marked.trimmed.1.fastq"


def read2_for_mapping(sample):
    if sample in REPAIR_SAMPLES:
        return f"{RESULTS}/repaired_reads/{sample}/{sample}.denovo_duplicates_marked.trimmed.repaired.2.fastq"
    return f"{RESULTS}/reads/{sample}/{sample}.denovo_duplicates_marked.trimmed.2.fastq"


def sra_read(run, mate):
    return f"{SRA_RAW_DIR}/{run}/{run}_{mate}.fastq"


rule all:
    input:
        f"{RESULTS}/summary/active_gus_abundance_all_samples.tsv",
        f"{RESULTS}/summary/active_gus_tpm_genome_wide_matrix.tsv",
        f"{RESULTS}/summary/active_gus_tpm_within_active_db_matrix.tsv",


rule diamond_makedb_active:
    input:
        REF_FAA,
    output:
        f"{RESULTS}/reference/active_gus.dmnd",
    params:
        db_prefix=f"{RESULTS}/reference/active_gus",
    shell:
        """
        mkdir -p {RESULTS}/reference
        diamond makedb --in {input} --db {params.db_prefix}
        """


rule extract_hmgi_gene_fna:
    input:
        HMGI_ARCHIVE,
    output:
        fna=f"{RESULTS}/genes/{{sample}}.gene-seqs.fna",
    params:
        member=lambda wc: f"./{wc.sample}-gene-seqs.fsa",
    shell:
        """
        mkdir -p {RESULTS}/genes
        tar -xOzf {input} {params.member} > {output.fna}
        """


rule prepare_gene_features:
    input:
        fna=f"{RESULTS}/genes/{{sample}}.gene-seqs.fna",
    output:
        faa=f"{RESULTS}/genes/{{sample}}.gene-seqs.faa",
        gff=f"{RESULTS}/genes/{{sample}}.gene-seqs.gff",
        saf=f"{RESULTS}/genes/{{sample}}.gene-seqs.saf",
    shell:
        """
        python3 scripts/prepare_hmgi_gene_features.py \
            --fna {input.fna} --faa {output.faa} --gff {output.gff} --saf {output.saf}
        """


rule extract_reads:
    input:
        tarball=lambda wc: os.path.join(READ_DIR, f"{wc.sample}.tar.bz2"),
    output:
        r1=f"{RESULTS}/reads/{{sample}}/{{sample}}.denovo_duplicates_marked.trimmed.1.fastq",
        r2=f"{RESULTS}/reads/{{sample}}/{{sample}}.denovo_duplicates_marked.trimmed.2.fastq",
    params:
        sample_dir=lambda wc: wc.sample,
        read_pair_limit=READ_PAIR_LIMIT,
    shell:
        """
        mkdir -p {RESULTS}/reads
        python3 scripts/extract_hmp_read_pair.py \
            --tarball {input.tarball} \
            --sample {wildcards.sample} \
            --r1 {output.r1} \
            --r2 {output.r2} \
            --read-pair-limit {params.read_pair_limit}
        """


rule sra_fasterq_dump:
    output:
        r1=f"{SRA_RAW_DIR}/{{run}}/{{run}}_1.fastq",
        r2=f"{SRA_RAW_DIR}/{{run}}/{{run}}_2.fastq",
    threads: THREADS
    wildcard_constraints:
        run=SRA_RUN_PATTERN,
    params:
        outdir=lambda wc: f"{SRA_RAW_DIR}/{wc.run}",
    shell:
        r"""
        mkdir -p {params.outdir}
        tmpdir="${{TMPDIR:-/tmp}}/{wildcards.run}_${{SLURM_JOB_ID:-local}}"
        mkdir -p "$tmpdir"
        fasterq-dump {wildcards.run} \
            --split-files --use-name \
            --seq-defline '@$sn/$ri' --qual-defline '+$sn/$ri' \
            --threads {threads} --temp "$tmpdir" --outdir {params.outdir} --progress
        rm -rf "$tmpdir"
        """


rule validate_sra_run_fastq:
    input:
        r1=f"{SRA_RAW_DIR}/{{run}}/{{run}}_1.fastq",
        r2=f"{SRA_RAW_DIR}/{{run}}/{{run}}_2.fastq",
    output:
        report=f"{RESULTS}/sra_rebuild/{SRA_SAMPLE}/{{run}}.fastq_validation.tsv",
    wildcard_constraints:
        run=SRA_RUN_PATTERN,
    params:
        expected=lambda wc: SRA_RUN_COUNTS[wc.run],
        target=SRA_REBUILD.get("target_read_id", ""),
    shell:
        """
        python3 scripts/validate_hmp_fastq_pairs.py \
            --r1 {input.r1} --r2 {input.r2} --report {output.report} \
            --expected-pairs {params.expected} --track-read-id {params.target}
        """


rule build_hmp_sequence_dedup_trim:
    input:
        "scripts/hmp_sequence_dedup_trim.cpp",
    output:
        executable=f"{RESULTS}/tools/hmp_sequence_dedup_trim",
    shell:
        """
        mkdir -p {RESULTS}/tools
        g++ -O3 -std=c++17 -fopenmp {input} -o {output.executable}
        """


rule sra_hmp_preprocess:
    input:
        executable=f"{RESULTS}/tools/hmp_sequence_dedup_trim",
        r1=[sra_read(run, 1) for run in SRA_RUNS],
        r2=[sra_read(run, 2) for run in SRA_RUNS],
        validation=expand(
            f"{RESULTS}/sra_rebuild/{SRA_SAMPLE}/{{run}}.fastq_validation.tsv",
            run=SRA_RUNS,
        ),
    output:
        r1=f"{RESULTS}/repaired_reads/{SRA_SAMPLE}/{SRA_SAMPLE}.denovo_duplicates_marked.trimmed.repaired.1.fastq",
        r2=f"{RESULTS}/repaired_reads/{SRA_SAMPLE}/{SRA_SAMPLE}.denovo_duplicates_marked.trimmed.repaired.2.fastq",
        singleton=f"{RESULTS}/repaired_reads/{SRA_SAMPLE}/{SRA_SAMPLE}.denovo_duplicates_marked.trimmed.singleton.fastq",
        report=f"{RESULTS}/repaired_reads/{SRA_SAMPLE}/{SRA_SAMPLE}.repair_report.tsv",
        validation=f"{RESULTS}/sra_rebuild/{SRA_SAMPLE}/{SRA_SAMPLE}.processed_fastq_validation.tsv",
    threads: THREADS
    params:
        r1_args=" ".join(f"--r1 {sra_read(run, 1)}" for run in SRA_RUNS),
        r2_args=" ".join(f"--r2 {sra_read(run, 2)}" for run in SRA_RUNS),
        seed=SRA_REBUILD.get("dedup_min_identical_bases", 5),
        diff=SRA_REBUILD.get("dedup_max_diff_rate", 0.03),
        meanq=SRA_REBUILD.get("dedup_min_mean_quality", 20),
        trimq=SRA_REBUILD.get("trim_quality", 3),
        minlen=SRA_REBUILD.get("minimum_trimmed_length", 60),
    shell:
        """
        mkdir -p {RESULTS}/repaired_reads/{SRA_SAMPLE} {RESULTS}/sra_rebuild/{SRA_SAMPLE}
        {input.executable} {params.r1_args} {params.r2_args} \
            --out-r1 {output.r1} --out-r2 {output.r2} \
            --singleton {output.singleton} --report {output.report} \
            --threads {threads} --min-identical-bases {params.seed} \
            --max-diff-rate {params.diff} --min-mean-quality {params.meanq} \
            --trim-quality {params.trimq} --min-length {params.minlen}
        python3 scripts/validate_hmp_fastq_pairs.py \
            --r1 {output.r1} --r2 {output.r2} --report {output.validation}
        """


ruleorder: sra_hmp_preprocess > repair_reads


rule repair_reads:
    input:
        r1=f"{RESULTS}/reads/{{sample}}/{{sample}}.denovo_duplicates_marked.trimmed.1.fastq",
        r2=f"{RESULTS}/reads/{{sample}}/{{sample}}.denovo_duplicates_marked.trimmed.2.fastq",
    output:
        r1=f"{RESULTS}/repaired_reads/{{sample}}/{{sample}}.denovo_duplicates_marked.trimmed.repaired.1.fastq",
        r2=f"{RESULTS}/repaired_reads/{{sample}}/{{sample}}.denovo_duplicates_marked.trimmed.repaired.2.fastq",
        report=f"{RESULTS}/repaired_reads/{{sample}}/{{sample}}.repair_report.tsv",
    shell:
        """
        python3 scripts/repair_hmp_fastq_pairs.py \
            --r1 {input.r1} --r2 {input.r2} \
            --out-r1 {output.r1} --out-r2 {output.r2} \
            --report {output.report}
        """


rule bowtie2_index_genes:
    input:
        fna=f"{RESULTS}/genes/{{sample}}.gene-seqs.fna",
    output:
        idx1=temp(f"{RESULTS}/mapping/{{sample}}.1.bt2"),
        idx2=temp(f"{RESULTS}/mapping/{{sample}}.2.bt2"),
        idx3=temp(f"{RESULTS}/mapping/{{sample}}.3.bt2"),
        idx4=temp(f"{RESULTS}/mapping/{{sample}}.4.bt2"),
        idxrev1=temp(f"{RESULTS}/mapping/{{sample}}.rev.1.bt2"),
        idxrev2=temp(f"{RESULTS}/mapping/{{sample}}.rev.2.bt2"),
    threads: THREADS_BOWTIE2
    shell:
        """
        mkdir -p {RESULTS}/mapping
        bowtie2-build {input.fna} {RESULTS}/mapping/{wildcards.sample} > {RESULTS}/mapping/{wildcards.sample}.bowtie2-build.log 2>&1
        """


rule bowtie2_align_genes:
    input:
        idx=rules.bowtie2_index_genes.output,
        r1=lambda wc: read1_for_mapping(wc.sample),
        r2=lambda wc: read2_for_mapping(wc.sample),
    output:
        bam=f"{RESULTS}/mapping/{{sample}}.sorted.bam",
        bai=f"{RESULTS}/mapping/{{sample}}.sorted.bam.bai",
    threads: THREADS_BOWTIE2
    shell:
        """
        bowtie2 -x {RESULTS}/mapping/{wildcards.sample} \
            -1 {input.r1} -2 {input.r2} -p {threads} 2> {RESULTS}/mapping/{wildcards.sample}.bowtie2.log \
          | samtools sort -@ {threads} -o {output.bam}
        samtools index {output.bam}
        """


rule featurecounts_genes:
    input:
        bam=f"{RESULTS}/mapping/{{sample}}.sorted.bam",
        saf=f"{RESULTS}/genes/{{sample}}.gene-seqs.saf",
    output:
        counts=f"{RESULTS}/counts/{{sample}}.featurecounts.txt",
    threads: THREADS_FEATURECOUNTS
    shell:
        """
        mkdir -p {RESULTS}/counts
        python3 scripts/count_hmp_gene_alignments.py \
            --bam {input.bam} --saf {input.saf} --output {output.counts}
        """


rule diamond_blastp_active:
    input:
        faa=f"{RESULTS}/genes/{{sample}}.gene-seqs.faa",
        db=f"{RESULTS}/reference/active_gus.dmnd",
    output:
        hits=f"{RESULTS}/search/{{sample}}.active_gus.tsv",
    params:
        db_prefix=f"{RESULTS}/reference/active_gus",
        evalue=config.get("max_evalue", 1e-5),
        max_target=config.get("max_target_seqs", 25),
    threads: THREADS_DIAMOND
    shell:
        """
        mkdir -p {RESULTS}/search
        diamond blastp --query {input.faa} --db {params.db_prefix} --out {output.hits} \
            --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen \
            --evalue {params.evalue} --max-target-seqs {params.max_target} --threads {threads}
        """


rule summarize_active_hits:
    input:
        hits=f"{RESULTS}/search/{{sample}}.active_gus.tsv",
        metadata=REF_META,
    output:
        best=f"{RESULTS}/summary/{{sample}}.active_gus_best_hits.tsv",
        clade=f"{RESULTS}/summary/{{sample}}.active_gus_clade_abundance.tsv",
    params:
        min_identity=config.get("min_identity", 95.0),
        min_query_cov=config.get("min_query_cov", 0.9),
        min_subject_cov=config.get("min_subject_cov", 0.9),
        max_evalue=config.get("max_evalue", 1e-5),
    shell:
        """
        mkdir -p {RESULTS}/summary
        python3 scripts/summarize_hmp_active_gus_hits.py \
            --hits {input.hits} --metadata {input.metadata} \
            --best-hits {output.best} --clade-summary {output.clade} \
            --min-identity {params.min_identity} \
            --min-query-cov {params.min_query_cov} \
            --min-subject-cov {params.min_subject_cov} \
            --max-evalue {params.max_evalue}
        """


rule abundance_active:
    input:
        hits=f"{RESULTS}/summary/{{sample}}.active_gus_best_hits.tsv",
        counts=f"{RESULTS}/counts/{{sample}}.featurecounts.txt",
        faa=f"{RESULTS}/genes/{{sample}}.gene-seqs.faa",
        gff=f"{RESULTS}/genes/{{sample}}.gene-seqs.gff",
    output:
        f"{RESULTS}/summary/{{sample}}.active_gus_abundance.tsv",
    shell:
        """
        python3 scripts/calculate_hmp_active_gus_abundance.py \
            --hits {input.hits} --counts {input.counts} --faa {input.faa} --gff {input.gff} \
            --output {output}
        """


rule merge_active_abundance:
    input:
        expand(f"{RESULTS}/summary/{{sample}}.active_gus_abundance.tsv", sample=SAMPLES),
    output:
        table=f"{RESULTS}/summary/active_gus_abundance_all_samples.tsv",
        genome_matrix=f"{RESULTS}/summary/active_gus_tpm_genome_wide_matrix.tsv",
        active_matrix=f"{RESULTS}/summary/active_gus_tpm_within_active_db_matrix.tsv",
    params:
        samples=",".join(SAMPLES),
    shell:
        """
        python3 scripts/merge_hmp_active_gus_abundance.py \
            --samples {params.samples} \
            --output-table {output.table} \
            --genome-matrix {output.genome_matrix} \
            --active-matrix {output.active_matrix} \
            {input}
        """
