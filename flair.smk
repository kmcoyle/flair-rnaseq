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
SAMPLES = pandas.read_csv('/projects/rmorin_scratch/ONT_scratch/results/flair/MCL_promethion_seq.tsv', sep = "\t")
RNAS = op.filter_samples(SAMPLES, seq_type = "promethION_mrna")

#print(SAMPLES)

## import files; bam must be mimimap2 with -ax splice
rule symlink_fastq_bam:
    input:
        fq = "/projects/rmorin/projects/gambl-repos/gambl-kcoyle/data/promethION_fastq/mrna/{sample_id}.fastq.gz",
        bam = "/projects/rmorin/projects/gambl-repos/gambl-kcoyle/data/promethION_bams/mrna/{sample_id}.hg38.bam"
    output:
        fq = "01-inputs/fastq/{sample_id}.fq.gz",
        bam = "01-inputs/bam/{sample_id}.bam"
    run:
        op.absolute_symlink(input.fq, output.fq)
        op.absolute_symlink(input.fq, output.fq)

## If your input sequences are Oxford nanopore reads, please use Pychopper before running Flair.
## Pychopper v2 is a tool to identify, orient and trim full-length Nanopore cDNA reads. The tool is also able to rescue fused reads.

rule pychopper:
    input:
        fq = 
    output:
        report = 
        unclassified_fq = 
        rescued_fq = 
        fq = 
    conda:
    shell:
        op.as_one_line("""
            pychopper -r report.pdf -u unclassified.fq -w rescued.fq input.fq full_length_output.fq
        """)


rule bam_to_bed:
    input:
        bam = str(rules.symlink_fastq_bam.output.bam)
    output:
        bed = "02-bed12/{sample_id}.bed"
    conda: "/projects/rmorin_scratch/conda_environments/flair/flair.yaml"
    shell:
        op.as_one_line("""
            bam2Bed12 -i {input.bam} > {output.bed}
        """)


## how to make this output to specific directory? just called "flair_..." atm.
rule flair_correct:
    input:
        bed = str(rules.bam_to_bed.output.bed),
        genome_gtf = "/projects/rmorin_scratch/ONT_scratch/results/ref/gencode.annotation.grch38.gtf".
        genome_fa = "/projects/rmorin_scratch/ONT_scratch/results/ref/GRCh38_no_alt.fa"
    output:
        corrected_bed = "03-flair-correct/{sample_id}.bed",
        inconsistent_bed = "03-flair-correct/{sample_id}.bed",
        cannot_verify_bed = "03-flair-correct/{sample_id}.bed"
    conda: "/projects/rmorin_scratch/conda_environments/flair/flair.yaml"
    shell:
        op.as_one_line("""
            flair correct -q {input.bed} -f {input.genome_gtf} -g {input.genome_fa}
        """)

## This needs to include all samples
## flair-corrected read bed files should be concatenated prior to running flair-collapse. 

rule concat_bed:
    input:
        all_beds = expand(str(rules.last_alignment.output.maf), 
        zip,
        sample_id=CONTROLS['sample_id'])
    output:
        single_bed = 

## In addition, all raw read fastq/fasta files should either be specified after --reads with space/comma separators 
## or concatenated into a single file.
rule flair_collapse:
    input:
        bed = "",
        fq = "",
        genome_fa = "",
        genome_gtf = ""
    output:
        bed = ""
        isoforms.bed

isoforms.gtf

isoforms.fa
    conda: "/projects/rmorin_scratch/conda_environments/flair/flair.yaml"
    shell:
        op.as_one_line("""
            flair collapse -g ../../ref/GRCh38_no_alt.fa -q ../03-flair-correct/flair_all_corrected.bed -r ../01-inputs/01-15563T.fastq.gz --gtf ../../ref/gencode.annotation.grch38.gtf
            """)


rule flair_quantify:
    input:
        bam = ""
    output:
        bed = ""
    conda: "/projects/rmorin_scratch/conda_environments/flair/flair.yaml"
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
