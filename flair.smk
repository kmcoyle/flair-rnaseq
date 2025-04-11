#!/usr/bin/env snakemake

import oncopipe as op
import os
import pandas as pd
import inspect

##### SNAKEMAKE COMMAND
# snakemake -p --use-conda --scheduler greedy --jobs 5000 --latency-wait 120 --keep-going --default-resources mem_mb=380000 disk_mb=2000  --cluster-sync "srun -p upgrade -n 1 -N 1 -J {rule} --mem {resources.mem_mb} --cpus-per-task {threads}" -s flair.smk all -np

##### SETUP VARIABLES

##### CONFIG FILES

##### SETUP SAMPLES
SAMPLES = pd.read_csv('/projects/rmorin_scratch/ONT_scratch/results/flair/MCL_promethion_seq.tsv', sep = "\t")
SAMPLES = op.filter_samples(SAMPLES, seq_type = "promethION_mrna")


# localrules:
#     flair_get_chrs

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


# checkpoint flair_get_chrs:
#     input:
#         chrs = "/projects/rmorin/projects/gambl-repos/gambl-kcoyle/ref/lcr-modules-references/genomes/hg38/genome_fasta/main_chromosomes.txt"
#     output:
#         chrs = "00-inputs/chroms/hg38.csv"
#     run:
#         # obtain list of main chromosomes
#         main_chrs = pd.read_csv(input.chrs, comment='#', sep='\t', header=None)
#         main_chrs = main_chrs.iloc[:, 0]#.astype(str).unique().tolist()
#         #main_chrs = pd.DataFrame(main_chrs)
#         # write out the file with all chromosomes
#         main_chrs.to_csv(output.chrs, index=False, header=False)

def flair_get_chr_fa(wildcards):
    chrs = "/projects/rmorin/projects/gambl-repos/gambl-kcoyle/ref/lcr-modules-references/genomes/hg38/genome_fasta/main_chromosomes.txt"
    with open(chrs) as file:
        chrs = file.read().rstrip("\n").split("\n")
    fa = expand(
        "02-bed12/" + "{{sample_id}}/chromosomes/{chrom}_reads.fq",
        chrom = chrs
    )
    return(fa)


def flair_get_chr_bed(wildcards):
    chrs = "/projects/rmorin/projects/gambl-repos/gambl-kcoyle/ref/lcr-modules-references/genomes/hg38/genome_fasta/main_chromosomes.txt"
    with open(chrs) as file:
        chrs = file.read().rstrip("\n").split("\n")
    bed = expand(
        "02-bed12/" + "{{sample_id}}/chromosomes/{chrom}_reads.bed",
        chrom = chrs
    )
    return(bed)



rule flair_split_align:
    input:
        bam = str(rules.flair_align.output.bam)
    output:
        fa = "02-bed12/" + "{sample_id}/chromosomes/{chrom}_reads.fa",
        fq = "02-bed12/" + "{sample_id}/chromosomes/{chrom}_reads.fq",
        bam = "02-bed12/" + "{sample_id}/chromosomes/{chrom}.bam",
        bai = "02-bed12/" + "{sample_id}/chromosomes/{chrom}.bam.bai",
        bed = "02-bed12/" + "{sample_id}/chromosomes/{chrom}_reads.bed",
        done = "02-bed12/" + "{sample_id}/chromosomes/{chrom}_done.txt",
    conda: "/projects/rmorin_scratch/ONT_scratch/results/flair/flair.yaml"
    params:
    shell:
        op.as_one_line("""
                samtools view -b {input.bam} {wildcards.chrom} > {output.bam}
                samtools index {output.bam} {output.bai}

                samtools fasta {output.bam} > {output.fa}
                gzip {output.fa}

                bam2Bed12 -i {output.bam} > {output.bed}
                touch {output.done}
        """)

def flair_get_chr_bed_correct(wildcards):
    chrs = "/projects/rmorin/projects/gambl-repos/gambl-kcoyle/ref/lcr-modules-references/genomes/hg38/genome_fasta/main_chromosomes.txt"
    with open(chrs) as file:
        chrs = file.read().rstrip("\n").split("\n")
    bed_corr = expand(
        "03-flair-correct/" +"{{sample_id}}/{chrom}",
        chrom = chrs
    )
    return(bed_corr)

rule flair_correct:
    input:
        bed = str(rules.flair_split_align.output.bed),
        genome_gtf = "/projects/rmorin_scratch/ONT_scratch/results/ref/gencode.annotation.grch38.gtf",
        genome_fa = "/projects/rmorin_scratch/ONT_scratch/results/ref/GRCh38_no_alt.fa"
    output:
        outdir = directory("03-flair-correct/{sample_id}/{chrom}"),
        bed = "03-flair-correct/{sample_id}/{chrom}_all_corrected.bed",
        done = "03-flair-correct/" + "{sample_id}/chromosomes/{chrom}_done.txt"
    conda: "/projects/rmorin_scratch/ONT_scratch/results/flair/flair.yaml"
    shell:
        op.as_one_line("""
            flair correct -q {input.bed} \
                -f {input.genome_gtf} \
                -g {input.genome_fa} \
                --output {output.outdir}
        """)

## In addition, all raw read fastq/fasta files should either be specified after --reads with space/comma separators 
## or concatenated into a single file.
## suggested flair collapse -g genome.fa --gtf gene_annotations.gtf -q reads.flair_all_corrected.bed -r reads.fastq
## --stringent --check_splice --generate_map --annotation_reliant generate
rule flair_collapse:
    input:
        bed = str(rules.flair_correct.output.bed),
        all_reads = str(rules.flair_split_align.output.fq),
        genome_fa = "/projects/rmorin_scratch/ONT_scratch/results/ref/GRCh38_no_alt.fa",
        genome_gtf = "/projects/rmorin_scratch/ONT_scratch/results/ref/gencode.annotation.grch38.gtf"
    output:
        bed = "04-flair-collapse/{sample_id}/{chrom}/isoforms.bed",
        gtf = "04-flair-collapse/{sample_id}/{chrom}/isoforms.gtf",
        fa = "04-flair-collapse/{sample_id}/{chrom}/flair.collapse.isoforms.fa"
    threads: 12
    log:
        stdout = "04-flair-collapse/{sample_id}/{chrom}/stdout.txt",
        stderr = "04-flair-collapse/{sample_id}/{chrom}/stderr.txt"
    params:
        out_dir = "04-flair-collapse/{sample_id}/{chrom}/"
    conda: "/projects/rmorin_scratch/ONT_scratch/results/flair/flair.yaml"
    shell:
        op.as_one_line("""
            flair collapse -g {input.genome_fa} -q {input.bed} -r {input.all_reads} --gtf {input.genome_gtf} \
            --threads {threads} --output {params.out_dir} \
            --stringent --check_splice --generate_map --annotation_reliant generate \
            > {log.stdout} 2> {log.stderr}
        """)


# rule flair_quantify:
#     input:
#         isoforms = str(rules.flair_collapse.output.fa),
#         manifest = "/projects/rmorin_scratch/ONT_scratch/results/flair/manifest.tsv"
#     output:
#         complete = "05-flair-quantify/complete"
#     threads: 12
#     params:
#         out_dir = "05-flair-quantify/"
#     conda: "/projects/rmorin_scratch/ONT_scratch/results/flair/flair.yaml"
#     shell:
#         op.as_one_line("""
#             flair quantify -r {input.manifest} -i {input.isoforms} \
#             --threads {threads} --output {params.out_dir} --sample_id_only && touch {output.complete}
#         """)



rule all:
    input:
        expand(
            [
                str(rules.flair_collapse.output.bed)
            ],
            #product(),  # Run expand() with zip(), not product()
            sample_id=SAMPLES['sample_id'],
            chrom=['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX'])




## This needs to include all samples
## flair-corrected read bed files should be concatenated prior to running flair-collapse. 

# rule concat_bed:
#     input:
#         all_beds = expand(str(rules.flair_correct.output.corrected_bed), 
#         zip,
#         sample_id=SAMPLES['sample_id'])
#     output:
#         single_bed = "04-flair-collapse/bed.bed"
#     conda: "/projects/rmorin_scratch/ONT_scratch/results/flair/bedtools.yaml"
#     shell:
#         op.as_one_line("""
#             cat {input.all_beds} | sort -k 1,1 -k2,2n | bedtools merge > {output.single_bed}
#         """)