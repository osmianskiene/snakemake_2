# List all FASTQ files based on pattern
FASTQ_FILES = glob_wildcards("data/{sample}.fastq.gz").sample
TRIMMED_FILES = glob_wildcards("data/{trimmed}_R1_001.fastq.gz").trimmed

# print(expand("data_output_multiqc/multiqc_report.html", output=OUTPUT_FILES))
print(TRIMMED_FILES)

rule all:
    input:
        # expand("data_output/{sample}_fastqc.html", sample=FASTQ_FILES)
        "data_output_multiqc/multiqc_report.html",
        expand("data_output_bbduk_trimmed/{sample}_bbduk.fastq", sample=FASTQ_FILES),
        # expand("data_output_fastqc_trimmed/{sample}_fastqc_trimmed.html", sample=FASTQ_FILES) 
        expand("data_output_fastqc_trimmed/{sample}_bbduk_fastqc.html", sample=FASTQ_FILES),
        "data_output_multiqc_trimmed/multiqc_report.html",
        "play_data_ref_annot/Genome",
        expand("data_output_star_mapped/{trimmed}_Aligned.sortedByCoord.out.bam", trimmed=TRIMMED_FILES), 
        expand("data_output_samtools_indexed/{trimmed}.bam.bai", trimmed=TRIMMED_FILES) 


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
        "data_output_bbduk_trimmed/{sample}_bbduk.fastq", 
    shell:
        "bbduk.sh in={input} out={output} ref=data/adapters.fa ktrim=r k=23 mink=11 hdist=1 tpe tbo qtrim=r trimq=10"

rule fastqc_trimmed:
    input:
        "data_output_bbduk_trimmed/{sample}_bbduk.fastq"
    output:
        "data_output_fastqc_trimmed/{sample}_bbduk_fastqc.html"
    shell:
        "fastqc -o data_output_fastqc_trimmed {input}"

rule multiqc_trimmed:
    input:
        expand("data_output_fastqc_trimmed/{sample}_bbduk_fastqc.html", sample=FASTQ_FILES)
    output:
        "data_output_multiqc_trimmed/multiqc_report.html", 
    shell:
        "multiqc data_output_fastqc_trimmed -o data_output_multiqc_trimmed"
 

# How to change output directory?
# Must be run --cores 8
rule star_index:
    input:
        "play_data_ref_annot/chr19_20Mb.fa"
    output:
        "play_data_ref_annot/Genome",
    shell:
        """
        star/STAR \
            --runThreadN 8 \
            --runMode genomeGenerate \
            --genomeDir play_data_ref_annot \
            --genomeFastaFiles {input} \
            --sjdbGTFfile play_data_ref_annot/chr19_20Mb.gtf
        """

# Must be run --cores 1
# Mapped paired, total 8 and not 16. Is it ok?
rule star_mapping:
    message: 
        "MONIKA MONIKA MONIKA MONIKA MONIKA MONIKA MONIKA MONIKA",
    input:
        read1="data_output_bbduk_trimmed/{trimmed}_R1_001_bbduk.fastq",
        read2="data_output_bbduk_trimmed/{trimmed}_R2_001_bbduk.fastq"
    output:
        "data_output_star_mapped/{trimmed}_Aligned.sortedByCoord.out.bam"
    # threads: 1
    shell:
        """
        star/STAR \
            --runThreadN 8 \
            --genomeDir play_data_ref_annot \
            --readFilesIn {input.read1} {input.read2} \
            --outSAMtype BAM SortedByCoordinate \
            --outFileNamePrefix data_output_star_mapped/{wildcards.trimmed}_
        """
rule samtools_index:
    input:
        # "sorted_reads/{sample}.bam"
        "data_output_star_mapped/{trimmed}_Aligned.sortedByCoord.out.bam"
    output:
        "data_output_samtools_indexed/{trimmed}.bam.bai"
    shell:
        "samtools index {input} -o {output}"





