# List all FASTQ files based on pattern
FASTQ_FILES = glob_wildcards("data/{sample}.fastq.gz").sample
# OUTPUT_FILES = glob_wildcards("data_output/{output}.fastq.gz").output

# print(expand("data_output_multiqc/multiqc_report.html", output=OUTPUT_FILES))
# print(output)

rule all:
    input:
        # expand("data_output/{sample}_fastqc.html", sample=FASTQ_FILES)
        "data_output_multiqc/multiqc_report.html",
        expand("data_output_bbduk_trimed/{sample}_bbduk.fastq", sample=FASTQ_FILES),
        # expand("data_output_fastqc_trimed/{sample}_fastqc_trimed.html", sample=FASTQ_FILES) 
        expand("data_output_fastqc_trimed/{sample}_bbduk_fastqc.html", sample=FASTQ_FILES) 

rule fastqc:
    input:
        "data/{sample}.fastq.gz"
    output:
        "data_output/{sample}_fastqc.html", 
    shell:
        "fastqc -o data_output {input}"

rule multiqc:
    input:
        expand("data_output/{sample}_fastqc.html", sample=FASTQ_FILES)
    output:
        "data_output_multiqc/multiqc_report.html", 
    shell:
        "multiqc data_output -o data_output_multiqc" 

rule bbduk_trim_barcodes:
    input:
        "data/{sample}.fastq.gz"
    output:
        "data_output_bbduk_trimed/{sample}_bbduk.fastq", 
    shell:
        "bbduk.sh in={input} out={output} ref=data/adapters.fa ktrim=r k=23 mink=11 hdist=1 tpe tbo qtrim=r trimq=10"

rule fastqc_trimed:
    input:
        "data_output_bbduk_trimed/{sample}_bbduk.fastq"
    output:
        "data_output_fastqc_trimed/{sample}_bbduk_fastqc.html"
    shell:
        "fastqc -o data_output_fastqc_trimed {input}"




# configfile: "config.yaml"
# SAMPLES = ["A", "B", "C"]
# rule all:
#     input:
#         "plots/quals.svg"
# rule download_data:
#     output:
#         "snakemake-tutorial-data-v5.4.5.tar.gz"
#     shell:
#         "curl -L -o {output} https://github.com/snakemake/snakemake-tutorial-data/archive/v5.4.5.tar.gz"
# rule extract_data:
#     input:
#         "snakemake-tutorial-data-v5.4.5.tar.gz"
#     output:
#         "data/genome.fa"
#     shell:
#         """
#         mkdir -p data
#         tar -xzf {input} -C data --strip-components 2
#         """
# rule bwa_map:
#     input:
#         "data/genome.fa",
#         "data/samples/{sample}.fastq"
#     output:
#         "mapped_reads/{sample}.bam"
#     shell:
#         "bwa mem {input} | samtools view -Sb - > {output}"
# rule samtools_sort:
#     input:
#         "mapped_reads/{sample}.bam"
#     output:
#         "sorted_reads/{sample}.bam"
#     shell:
#         "samtools sort -T sorted_reads/{wildcards.sample} "
#         "-O bam {input} > {output}"
# rule samtools_index:
#     input:
#         "sorted_reads/{sample}.bam"
#     output:
#         "sorted_reads/{sample}.bam.bai"
#     shell:
#         "samtools index {input}"
# rule bcftools_call:
#     input:
#         fa="data/genome.fa",
#         bam=expand("sorted_reads/{sample}.bam", sample=config["samples"]),
#         bai=expand("sorted_reads/{sample}.bam.bai", sample=config["samples"])
#     output:
#         "calls/all.vcf"
#     shell:
#         "bcftools mpileup -f {input.fa} {input.bam} | "
#         "bcftools call -mv - > {output}"
# rule plot_quals:
#     input:
#         "calls/all.vcf"
#     output:
#         "plots/quals.svg"
#     script:
#         "scripts/plot-quals.py"
