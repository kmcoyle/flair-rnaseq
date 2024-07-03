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
        fq = "/projects/rmorin/projects/gambl-repos/gambl-kcoyle/data/promethION_fastq/mrna/{sample_id}.fastq.gz"
    output:
        fq = "00-inputs/fastq/{sample_id}.fq.gz",
        bam = "00-inputs/bam/{sample_id}.bam"
    run:
        op.absolute_symlink(input.fq, output.fq)
        op.absolute_symlink(input.fq, output.fq)

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
    conda: "/projects/rmorin_scratch/conda_environments/flair/pychopper.yaml"
    shell:
        op.as_one_line("""
            pychopper -r {output.pdf} -u {output.unclassified_fq} -w {output.rescued_fq} {input.fq} | gzip > {output.fq}
        """)

 ## damn if we use pychopper I have to run minimap again

rule flair_align:
    input:
        genome_fa = "/projects/rmorin_scratch/ONT_scratch/results/ref/GRCh38_no_alt.fa",
        fq = str(rules.pychopper.output.fq)
    output:
        bam = "02-bed12/{sample_id}/flair.aligned.bam",
        bai = "02-bed12/{sample_id}/flair.aligned.bam.bai",
        bed = "02-bed12/{sample_id}/flair.aligned.bed",
    conda: "/projects/rmorin_scratch/conda_environments/flair/flair.yaml"
    params:
        threads = 24
    shell:
        op.as_one_line("""
            flair align -g {input.genome_fa} -r {input.fq} --output 02-bed12/{sample_id}/ --threads {params.threads}
        """)


## how to make this output to specific directory? just called "flair_..." atm.
rule flair_correct:
    input:
        bed = str(rules.flair_align.output.bed),
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
    conda:
    shell:
        op.as_one_line("""
            
        """)

## In addition, all raw read fastq/fasta files should either be specified after --reads with space/comma separators 
## or concatenated into a single file.
rule flair_collapse:
    input:
        bed = "../03-flair-correct/flair_all_corrected.bed",
        fq = "",
        genome_fa = "../../ref/GRCh38_no_alt.fa",
        genome_gtf = "../../ref/gencode.annotation.grch38.gtf"
    output:
        bed = "isoforms.bed"
        gtf = "isoforms.gtf"
        fa = "isoforms.fa"
    conda: "/projects/rmorin_scratch/conda_environments/flair/flair.yaml"
    shell:
        op.as_one_line("""
            flair collapse -g {input.genome_fa} -q {input.bed} -r {input.fq} --gtf {input.genome_gtf}
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
