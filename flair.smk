#!/usr/bin/env snakemake

import oncopipe as op
import os
import pandas

##### SNAKEMAKE COMMAND
# snakemake -p --use-conda --scheduler greedy --jobs 5000 --latency-wait 120 --keep-going --default-resources mem_mb=380000 disk_mb=2000  --cluster-sync "srun -p upgrade -n 1 -N 1 -J {rule} --mem {resources.mem_mb} --cpus-per-task {threads}" -s last.dnarrange.smk all -np
# snakemake -p --use-conda --scheduler greedy --jobs 5000 --latency-wait 120 --keep-going --default-resources mem_mb=2000 disk_mb=2000  --cluster-sync "srun -p upgrade -n 1 -N 1 -J {rule} --mem {resources.mem_mb} --cpus-per-task {threads}" -s last.dnarrange.smk all -np

##### SETUP VARIABLES

##### CONFIG FILES

##### SETUP SAMPLES
SAMPLES = pandas.read_csv('/projects/rmorin_scratch/ONT_scratch/samples.last.tsv', sep = "\t")
CONTROLS = op.filter_samples(SAMPLES, pathology = ["BL", "DLBCL"])
MCLS = op.filter_samples(SAMPLES, pathology = "MCL")

#print(SAMPLES)
rule symlink_fastq_bam:
    input:
        fq = "/projects/rmorin/projects/gambl-repos/gambl-kcoyle/data/promethION_fastq/genome/{sample_id}.fastq.gz",
        bam = ""
    output:
        fq = "01-inputs/fastq/{sample_id}.fq.gz",
        bam = "01-inputs/fastq/{sample_id}.bam"
    run:
        op.absolute_symlink(input.fq, output.fq)
        op.absolute_symlink(input.fq, output.fq)



rule bam_to_bed:
    input:
        bam = ""
    output:
        bed = ""
    conda:
    shell:
        op.as_one_line("""
            bam2Bed12 -i ../01-inputs/01-15563T.hg38.bam > 01-15563T.hg38.bed
        """)

rule flair_correct:
    input:
        bam = ""
    output:
        bed = ""
    conda:
    shell:
        op.as_one_line("""
            bam2Bed12 -i ../01-inputs/01-15563T.hg38.bam > 01-15563T.hg38.bed
        """)

rule flair_collapse:
    input:
        bam = ""
    output:
        bed = ""
    conda:
    shell:
        op.as_one_line("""
            bam2Bed12 -i ../01-inputs/01-15563T.hg38.bam > 01-15563T.hg38.bed
        """)

rule flair_quantify:
    input:
        bam = ""
    output:
        bed = ""
    conda:
    shell:
        op.as_one_line("""
            bam2Bed12 -i ../01-inputs/01-15563T.hg38.bam > 01-15563T.hg38.bed
        """)



rule all:
    input:
        expand(rules.last_alignment.output.maf,
        zip,
        sample_id=SAMPLES['sample_id']),
        expand(rules.get_all_plot_ids.output.txt,
        zip,
        sample_id=MCLS['sample_id']),
        expand(rules.parse_maf.output.alignment_df,
        zip,
        sample_id=MCLS['sample_id'])
