#!/usr/bin/env snakemake

import oncopipe as op
import os
import pandas

##### SNAKEMAKE COMMAND
# snakemake -p --use-conda --scheduler greedy --jobs 5000 --latency-wait 120 --keep-going --default-resources mem_mb=380000 disk_mb=2000  --cluster-sync "srun -p upgrade -n 1 -N 1 -J {rule} --mem {resources.mem_mb} --cpus-per-task {threads}" -s last.dnarrange.smk all -np
# snakemake -p --use-conda --scheduler greedy --jobs 5000 --latency-wait 120 --keep-going --default-resources mem_mb=2000 disk_mb=2000  --cluster-sync "srun -p upgrade -n 1 -N 1 -J {rule} --mem {resources.mem_mb} --cpus-per-task {threads}" -s flair.smk all -np

##### SETUP VARIABLES

##### CONFIG FILES

##### SETUP SAMPLES
SAMPLES = pandas.read_csv('/projects/rmorin_scratch/ONT_scratch/results/flair/MCL_promethion_seq.tsv', sep = "\t")
SAMPLES = op.filter_samples(SAMPLES, seq_type = "promethION_mrna")

#print(SAMPLES)

## import files; bam must be mimimap2 with -ax splice
rule symlink_fastq_bam:
    input:
        fq = "/projects/rmorin/projects/gambl-repos/gambl-kcoyle/data/promethION_fastq/mrna/{sample_id}.fastq.gz",
        bam = "/projects/rmorin/projects/gambl-repos/gambl-kcoyle/data/promethION_bams/mrna/{sample_id}.hg38.bam"
    output:
        fq = "00-inputs/fastq/{sample_id}.fq.gz",
        bam = "00-inputs/bam/{sample_id}.bam"
    run:
        op.absolute_symlink(input.fq, output.fq),
        op.absolute_symlink(input.bam, output.bam)

## If your input sequences are Oxford nanopore reads, please use Pychopper before running Flair.
## Pychopper v2 is a tool to identify, orient and trim full-length Nanopore cDNA reads. The tool is also able to rescue fused reads.
# we are using 2.7.9 and gz capabilities added in v2.6.0

#try not specifying an output file. That should make pychopper write output to stdout so you can stream directly into gzip.

rule pychopper:
    input:
        fq = str(rules.symlink_fastq_bam.output.fq)
    output:
        report = "01-pychopper/reports/{sample_id}.pdf",
        unclassified_fq = "01-pychopper/fastq/{sample_id}_unclassified.fq",
        rescued_fq = "01-pychopper/fastq/{sample_id}_rescued.fq",
        fq = "01-pychopper/fastq/{sample_id}_full_length.fq.gz"
    params:
        threads = 24
    conda: "/projects/rmorin_scratch/ONT_scratch/results/flair/pychopper.yaml"
    shell:
        op.as_one_line("""
            pychopper -r {output.report} -u {output.unclassified_fq} -w {output.rescued_fq} -t {params.threads} {input.fq} | gzip > {output.fq}
        """)

 ## damn if we use pychopper I have to run minimap again

rule flair_align:
    input:
        genome_fa = "/projects/rmorin_scratch/ONT_scratch/results/ref/GRCh38_no_alt.fa",
        fq = "/projects/rmorin_scratch/ONT_scratch/results/flair/01-pychopper/fastq/{sample_id}_full_length.fq.gz"
    output:
        bam = "02-bed12/{sample_id}/flair.aligned.bam",
        bai = "02-bed12/{sample_id}/flair.aligned.bam.bai",
        bed = "02-bed12/{sample_id}/flair.aligned.bed",
    conda: "/projects/rmorin_scratch/ONT_scratch/results/flair/flair.yaml"
    params:
        threads = 36
    shell:
        op.as_one_line("""
            cd 02-bed12/{wildcards.sample_id} &&
            flair align -g {input.genome_fa} -r {input.fq} --threads {params.threads}
        """)
## not sure where this is being output???

## how to make this output to specific directory? just called "flair_..." atm.
rule flair_correct:
    input:
        bed = str(rules.flair_align.output.bed),
        genome_gtf = "/projects/rmorin_scratch/ONT_scratch/results/ref/gencode.annotation.grch38.gtf",
        genome_fa = "/projects/rmorin_scratch/ONT_scratch/results/ref/GRCh38_no_alt.fa"
    output:
        corrected_bed = "03-flair-correct/{sample_id}/all_corrected.bed",
        inconsistent_bed = "03-flair-correct/{sample_id}/all_inconsistent.bed",
        cannot_verify_bed = "03-flair-correct/{sample_id}/cannot_verify.bed"
    conda: "/projects/rmorin_scratch/ONT_scratch/results/flair/flair.yaml"
    shell:
        op.as_one_line("""
            cd 03-flair-correct/{wildcards.sample_id} &&
            flair correct -q {input.bed} -f {input.genome_gtf} -g {input.genome_fa}
        """)

## This needs to include all samples
## flair-corrected read bed files should be concatenated prior to running flair-collapse. 

rule concat_bed:
    input:
        all_beds = expand(str(rules.flair_correct.output.corrected_bed), 
        zip,
        sample_id=SAMPLES['sample_id'])
    output:
        single_bed = "bed.bed"
    conda: "/projects/rmorin_scratch/ONT_scratch/results/flair/bedtools.yaml"
    shell:
        op.as_one_line("""
            cat {input.all_beds} | sort -k 1,1 -k2,2n | bedtools merge > {single_bed}
        """)

## In addition, all raw read fastq/fasta files should either be specified after --reads with space/comma separators 
## or concatenated into a single file.
rule flair_collapse:
    input:
        bed = str(rules.concat_bed.output.bed),,
        all_reads = expand(str(rules.pychopper.output.fq), 
        zip,
        sample_id=SAMPLES['sample_id']),
        genome_fa = "../../ref/GRCh38_no_alt.fa",
        genome_gtf = "../../ref/gencode.annotation.grch38.gtf"
    output:
        bed = "isoforms.bed",
        gtf = "isoforms.gtf",
        fa = "isoforms.fa"
    conda: "/projects/rmorin_scratch/ONT_scratch/results/flair/flair.yaml"
    shell:
        op.as_one_line("""
            flair collapse -g {input.genome_fa} -q {input.bed} -r {input.fq} --gtf {input.genome_gtf}
            """)


rule flair_quantify:
    input:
        bam = ""
    output:
        bed = ""
    conda: "/projects/rmorin_scratch/ONT_scratch/results/flair/flair.yaml"
    shell:
        op.as_one_line("""
            bam2Bed12 -i ../01-inputs/01-15563T.hg38.bam > 01-15563T.hg38.bed
        """)



rule all:
    input:
        expand(rules.flair_align.output.bed,
        zip,
        sample_id=SAMPLES['sample_id'])
