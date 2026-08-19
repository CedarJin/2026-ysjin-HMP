# HMP per-sample assembly + GUS abundance workflow
# (aligned with 2026-ysjin-gus: MEGAHIT -> bowtie2 -> Prodigal -> featureCounts -> DIAMOND)
#
# Designed for Slurm job arrays: one sample per task via --config sample=SRS....
#
# From project root:
#   snakemake -s workflow/hmp_assembly_gus.smk --configfile config/hmp_assembly_gus.yaml \
#     --config sample=SRS011061 -c 8 -n

import os
from pathlib import Path

from snakemake.utils import min_version

min_version("7.0")

configfile: "config/hmp_assembly_gus.yaml"

shell.prefix("set -euo pipefail; ")

RESULTS = config["results_dir"]
READ_STAGE = config["read_stage_dir"]
ASSEMBLY_DIR = config["assembly_dir"]
CONTIGS_BASENAME = config["contigs_basename"]
MINLEN = int(config["min_contig_length"])
HMP_READS = config["hmp_reads_dir"]
HMP_REPAIRED = config["hmp_repaired_reads_dir"]
REPAIR_SAMPLES = set(config.get("repair_samples", []))
GUS_SCRIPTS = config["gus_scripts_dir"]
REF_DIR = Path(RESULTS) / "reference"
TOTAL_FAA = config["total_faa"]
TOTAL_META = config["total_metadata"]
CAND_FAA = config["candidate_faa"]
CAND_META = config["candidate_metadata"]
ACTIVE_FAA = config["active_faa"]
ACTIVE_META = config["active_metadata"]
DIAMOND = config.get("diamond", {"evalue": 1e-5, "max_target_seqs": 25})
SUMMARIZE = config.get(
    "summarize",
    {"min_identity": 40.0, "min_query_cov": 0.5, "min_subject_cov": 0.5, "max_evalue": 1e-5},
)
STRICT = config.get(
    "strict_active",
    {"min_identity": 95.0, "min_query_cov": 0.9, "min_subject_cov": 0.9, "max_evalue": 1e-5},
)

THREADS_MEGAHIT = int(config.get("threads_megahit", config.get("threads", 8)))
THREADS_BOWTIE2 = int(config.get("threads_bowtie2", config.get("threads", 8)))
THREADS_FEATURECOUNTS = int(config.get("threads_featurecounts", config.get("threads", 8)))
THREADS_PRODIGAL = int(config.get("threads_prodigal", config.get("threads", 8)))
THREADS_DIAMOND = int(config.get("threads_diamond", config.get("threads", 8)))

SAMPLE = config.get("sample")
if not SAMPLE:
    raise WorkflowError("Pass a single sample via --config sample=SRS....")


def source_r1(sample):
    if sample in REPAIR_SAMPLES:
        return f"{HMP_REPAIRED}/{sample}/{sample}.denovo_duplicates_marked.trimmed.repaired.1.fastq"
    return f"{HMP_READS}/{sample}/{sample}.denovo_duplicates_marked.trimmed.1.fastq"


def source_r2(sample):
    if sample in REPAIR_SAMPLES:
        return f"{HMP_REPAIRED}/{sample}/{sample}.denovo_duplicates_marked.trimmed.repaired.2.fastq"
    return f"{HMP_READS}/{sample}/{sample}.denovo_duplicates_marked.trimmed.2.fastq"


def contigs(sample):
    return os.path.join(ASSEMBLY_DIR, sample, CONTIGS_BASENAME)


wildcard_constraints:
    sample=r"SRS\d+"


rule all:
    input:
        f"{RESULTS}/summary/{SAMPLE}.min{MINLEN}.total_gus_abundance.tsv",
        f"{RESULTS}/summary/{SAMPLE}.min{MINLEN}.candidate_gus_abundance.tsv",
        f"{RESULTS}/summary/pident95_cov90/abundance/{SAMPLE}.min{MINLEN}.active_gus_abundance.tsv",
        f"{RESULTS}/summary/pident95_cov90/hits/{SAMPLE}.min{MINLEN}.active_gus_best_hits.tsv",


rule stage_reads:
    input:
        r1=lambda wc: source_r1(wc.sample),
        r2=lambda wc: source_r2(wc.sample),
    output:
        r1=f"{READ_STAGE}/{{sample}}/{{sample}}_1.fastq",
        r2=f"{READ_STAGE}/{{sample}}/{{sample}}_2.fastq",
    shell:
        """
        mkdir -p {READ_STAGE}/{wildcards.sample}
        ln -sfr {input.r1} {output.r1}
        ln -sfr {input.r2} {output.r2}
        """


rule megahit:
    input:
        r1=f"{READ_STAGE}/{{sample}}/{{sample}}_1.fastq",
        r2=f"{READ_STAGE}/{{sample}}/{{sample}}_2.fastq",
    output:
        contigs=f"{ASSEMBLY_DIR}/{{sample}}/{CONTIGS_BASENAME}",
    log:
        f"{RESULTS}/logs/megahit/{{sample}}.log",
    threads: THREADS_MEGAHIT
    params:
        odir=lambda wc: os.path.join(ASSEMBLY_DIR, wc.sample),
        mlen=MINLEN,
    shell:
        r"""
        mkdir -p "$(dirname {log})"
        rm -rf "{params.odir}"
        megahit -1 {input.r1} -2 {input.r2} -o "{params.odir}" \
            --min-contig-len {params.mlen} -t {threads} > {log} 2>&1
        """


rule diamond_makedb_total:
    input:
        TOTAL_FAA,
    output:
        f"{RESULTS}/reference/gus_279_plus_fmn.dmnd",
    params:
        prefix=f"{RESULTS}/reference/gus_279_plus_fmn",
    shell:
        """
        mkdir -p {RESULTS}/reference
        diamond makedb --in {input} --db {params.prefix}
        """


rule diamond_makedb_candidate:
    input:
        CAND_FAA,
    output:
        f"{RESULTS}/reference/candidate_gus.dmnd",
    params:
        prefix=f"{RESULTS}/reference/candidate_gus",
    shell:
        """
        mkdir -p {RESULTS}/reference
        diamond makedb --in {input} --db {params.prefix}
        """


rule diamond_makedb_active:
    input:
        ACTIVE_FAA,
    output:
        f"{RESULTS}/reference/active_gus.dmnd",
    params:
        prefix=f"{RESULTS}/reference/active_gus",
    shell:
        """
        mkdir -p {RESULTS}/reference
        diamond makedb --in {input} --db {params.prefix}
        """


rule bowtie2_index:
    input:
        contigs=lambda wc: contigs(wc.sample),
    output:
        idx1=temp(f"{RESULTS}/mapping/{{sample}}.min{MINLEN}.1.bt2"),
        idx2=temp(f"{RESULTS}/mapping/{{sample}}.min{MINLEN}.2.bt2"),
        idx3=temp(f"{RESULTS}/mapping/{{sample}}.min{MINLEN}.3.bt2"),
        idx4=temp(f"{RESULTS}/mapping/{{sample}}.min{MINLEN}.4.bt2"),
        idxrev1=temp(f"{RESULTS}/mapping/{{sample}}.min{MINLEN}.rev.1.bt2"),
        idxrev2=temp(f"{RESULTS}/mapping/{{sample}}.min{MINLEN}.rev.2.bt2"),
    threads: THREADS_BOWTIE2
    shell:
        """
        mkdir -p {RESULTS}/mapping
        bowtie2-build {input.contigs} {RESULTS}/mapping/{wildcards.sample}.min{MINLEN} > /dev/null 2>&1
        """


rule bowtie2_align:
    input:
        idx=rules.bowtie2_index.output,
        r1=f"{READ_STAGE}/{{sample}}/{{sample}}_1.fastq",
        r2=f"{READ_STAGE}/{{sample}}/{{sample}}_2.fastq",
    output:
        bam=f"{RESULTS}/mapping/{{sample}}.min{MINLEN}.sorted.bam",
        bai=f"{RESULTS}/mapping/{{sample}}.min{MINLEN}.sorted.bam.bai",
    threads: THREADS_BOWTIE2
    shell:
        """
        bowtie2 -x {RESULTS}/mapping/{wildcards.sample}.min{MINLEN} \
            -1 {input.r1} -2 {input.r2} -p {threads} \
            2> {RESULTS}/mapping/{wildcards.sample}.min{MINLEN}.bowtie2.log \
          | samtools sort -@ {threads} -o {output.bam}
        samtools index {output.bam}
        """


rule prodigal:
    input:
        contigs=lambda wc: contigs(wc.sample),
    output:
        faa=f"{RESULTS}/prodigal/{{sample}}.min{MINLEN}.faa",
        ffn=f"{RESULTS}/prodigal/{{sample}}.min{MINLEN}.ffn",
        gff=f"{RESULTS}/prodigal/{{sample}}.min{MINLEN}.gff",
    threads: THREADS_PRODIGAL
    shell:
        """
        mkdir -p {RESULTS}/prodigal
        prodigal -i {input.contigs} -p meta -q \
            -a {output.faa} -d {output.ffn} -o {output.gff} -f gff
        """


rule featurecounts:
    input:
        bam=f"{RESULTS}/mapping/{{sample}}.min{MINLEN}.sorted.bam",
        gff=f"{RESULTS}/prodigal/{{sample}}.min{MINLEN}.gff",
    output:
        counts=f"{RESULTS}/counts/{{sample}}.min{MINLEN}.featurecounts.txt",
    threads: THREADS_FEATURECOUNTS
    shell:
        """
        mkdir -p {RESULTS}/counts
        featureCounts -T {threads} -p -B -C -t CDS -g ID -a {input.gff} \
            -o {output.counts} {input.bam}
        """


rule diamond_blastp_total:
    input:
        faa=f"{RESULTS}/prodigal/{{sample}}.min{MINLEN}.faa",
        db=f"{RESULTS}/reference/gus_279_plus_fmn.dmnd",
    output:
        hits=f"{RESULTS}/search/{{sample}}.min{MINLEN}.total_gus.tsv",
    params:
        db_prefix=f"{RESULTS}/reference/gus_279_plus_fmn",
        evalue=DIAMOND["evalue"],
        max_target=DIAMOND["max_target_seqs"],
    threads: THREADS_DIAMOND
    shell:
        """
        mkdir -p {RESULTS}/search
        diamond blastp --query {input.faa} --db {params.db_prefix} --out {output.hits} \
            --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen \
            --evalue {params.evalue} --max-target-seqs {params.max_target} --threads {threads}
        """


rule diamond_blastp_candidate:
    input:
        faa=f"{RESULTS}/prodigal/{{sample}}.min{MINLEN}.faa",
        db=f"{RESULTS}/reference/candidate_gus.dmnd",
    output:
        hits=f"{RESULTS}/search/{{sample}}.min{MINLEN}.candidate_gus.tsv",
    params:
        db_prefix=f"{RESULTS}/reference/candidate_gus",
        evalue=DIAMOND["evalue"],
        max_target=DIAMOND["max_target_seqs"],
    threads: THREADS_DIAMOND
    shell:
        """
        mkdir -p {RESULTS}/search
        diamond blastp --query {input.faa} --db {params.db_prefix} --out {output.hits} \
            --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen \
            --evalue {params.evalue} --max-target-seqs {params.max_target} --threads {threads}
        """


rule diamond_blastp_active:
    input:
        faa=f"{RESULTS}/prodigal/{{sample}}.min{MINLEN}.faa",
        db=f"{RESULTS}/reference/active_gus.dmnd",
    output:
        hits=f"{RESULTS}/search/{{sample}}.min{MINLEN}.active_gus.tsv",
    params:
        db_prefix=f"{RESULTS}/reference/active_gus",
        evalue=DIAMOND["evalue"],
        max_target=DIAMOND["max_target_seqs"],
    threads: THREADS_DIAMOND
    shell:
        """
        mkdir -p {RESULTS}/search
        diamond blastp --query {input.faa} --db {params.db_prefix} --out {output.hits} \
            --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen \
            --evalue {params.evalue} --max-target-seqs {params.max_target} --threads {threads}
        """


rule summarize_total_hits:
    input:
        hits=f"{RESULTS}/search/{{sample}}.min{MINLEN}.total_gus.tsv",
        contigs=lambda wc: contigs(wc.sample),
        metadata=TOTAL_META,
    output:
        best=temp(f"{RESULTS}/summary/{{sample}}.min{MINLEN}.total_gus_best_hits.tsv"),
        clade=temp(f"{RESULTS}/summary/{{sample}}.min{MINLEN}.total_gus_clade_abundance.tsv"),
    params:
        min_identity=SUMMARIZE["min_identity"],
        min_query_cov=SUMMARIZE["min_query_cov"],
        min_subject_cov=SUMMARIZE["min_subject_cov"],
        max_evalue=SUMMARIZE["max_evalue"],
    shell:
        """
        mkdir -p {RESULTS}/summary
        python3 {GUS_SCRIPTS}/summarize_gus_hits.py \
            --hits {input.hits} --contigs {input.contigs} --metadata {input.metadata} \
            --best-hits {output.best} --clade-summary {output.clade} \
            --min-identity {params.min_identity} \
            --min-query-cov {params.min_query_cov} \
            --min-subject-cov {params.min_subject_cov} \
            --max-evalue {params.max_evalue}
        """


rule summarize_candidate_hits:
    input:
        hits=f"{RESULTS}/search/{{sample}}.min{MINLEN}.candidate_gus.tsv",
        contigs=lambda wc: contigs(wc.sample),
        metadata=CAND_META,
    output:
        best=temp(f"{RESULTS}/summary/{{sample}}.min{MINLEN}.candidate_gus_best_hits.tsv"),
        clade=temp(f"{RESULTS}/summary/{{sample}}.min{MINLEN}.candidate_gus_clade_abundance.tsv"),
    params:
        min_identity=SUMMARIZE["min_identity"],
        min_query_cov=SUMMARIZE["min_query_cov"],
        min_subject_cov=SUMMARIZE["min_subject_cov"],
        max_evalue=SUMMARIZE["max_evalue"],
    shell:
        """
        mkdir -p {RESULTS}/summary
        python3 {GUS_SCRIPTS}/summarize_gus_hits.py \
            --hits {input.hits} --contigs {input.contigs} --metadata {input.metadata} \
            --best-hits {output.best} --clade-summary {output.clade} \
            --min-identity {params.min_identity} \
            --min-query-cov {params.min_query_cov} \
            --min-subject-cov {params.min_subject_cov} \
            --max-evalue {params.max_evalue}
        """


rule summarize_active_hits_strict:
    input:
        hits=f"{RESULTS}/search/{{sample}}.min{MINLEN}.active_gus.tsv",
        contigs=lambda wc: contigs(wc.sample),
        metadata=ACTIVE_META,
    output:
        best=f"{RESULTS}/summary/pident95_cov90/hits/{{sample}}.min{MINLEN}.active_gus_best_hits.tsv",
        clade=f"{RESULTS}/summary/pident95_cov90/hits/{{sample}}.min{MINLEN}.active_gus_clade_abundance.tsv",
    params:
        min_identity=STRICT["min_identity"],
        min_query_cov=STRICT["min_query_cov"],
        min_subject_cov=STRICT["min_subject_cov"],
        max_evalue=STRICT["max_evalue"],
    shell:
        """
        mkdir -p {RESULTS}/summary/pident95_cov90/hits
        python3 {GUS_SCRIPTS}/summarize_gus_hits.py \
            --hits {input.hits} --contigs {input.contigs} --metadata {input.metadata} \
            --best-hits {output.best} --clade-summary {output.clade} \
            --min-identity {params.min_identity} \
            --min-query-cov {params.min_query_cov} \
            --min-subject-cov {params.min_subject_cov} \
            --max-evalue {params.max_evalue}
        """


rule abundance_total:
    input:
        hits=f"{RESULTS}/summary/{{sample}}.min{MINLEN}.total_gus_best_hits.tsv",
        counts=f"{RESULTS}/counts/{{sample}}.min{MINLEN}.featurecounts.txt",
        faa=f"{RESULTS}/prodigal/{{sample}}.min{MINLEN}.faa",
        gff=f"{RESULTS}/prodigal/{{sample}}.min{MINLEN}.gff",
    output:
        f"{RESULTS}/summary/{{sample}}.min{MINLEN}.total_gus_abundance.tsv",
    shell:
        """
        python3 {GUS_SCRIPTS}/calculate_gus_abundance.py \
            --hits {input.hits} --counts {input.counts} --faa {input.faa} --gff {input.gff} \
            --output {output}
        """


rule abundance_candidate:
    input:
        hits=f"{RESULTS}/summary/{{sample}}.min{MINLEN}.candidate_gus_best_hits.tsv",
        counts=f"{RESULTS}/counts/{{sample}}.min{MINLEN}.featurecounts.txt",
        faa=f"{RESULTS}/prodigal/{{sample}}.min{MINLEN}.faa",
        gff=f"{RESULTS}/prodigal/{{sample}}.min{MINLEN}.gff",
    output:
        f"{RESULTS}/summary/{{sample}}.min{MINLEN}.candidate_gus_abundance.tsv",
    shell:
        """
        python3 {GUS_SCRIPTS}/calculate_gus_abundance.py \
            --hits {input.hits} --counts {input.counts} --faa {input.faa} --gff {input.gff} \
            --output {output}
        """


rule abundance_active_strict:
    input:
        hits=f"{RESULTS}/summary/pident95_cov90/hits/{{sample}}.min{MINLEN}.active_gus_best_hits.tsv",
        counts=f"{RESULTS}/counts/{{sample}}.min{MINLEN}.featurecounts.txt",
        faa=f"{RESULTS}/prodigal/{{sample}}.min{MINLEN}.faa",
        gff=f"{RESULTS}/prodigal/{{sample}}.min{MINLEN}.gff",
    output:
        f"{RESULTS}/summary/pident95_cov90/abundance/{{sample}}.min{MINLEN}.active_gus_abundance.tsv",
    shell:
        """
        mkdir -p {RESULTS}/summary/pident95_cov90/abundance
        python3 {GUS_SCRIPTS}/calculate_gus_abundance.py \
            --hits {input.hits} --counts {input.counts} --faa {input.faa} --gff {input.gff} \
            --output {output}
        """
