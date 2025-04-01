# List all FASTQ files based on pattern
FASTQ_FILES = glob_wildcards("data/{sample}.fastq.gz").sample
TRIMMED_FILES = glob_wildcards("data/{trimmed}_R1_001.fastq.gz").trimmed
TRIMMED_FILES = glob_wildcards("data/{trimmed}_R1_001.fastq.gz").trimmed
COLLIBRI_FILES = glob_wildcards("data/Collibri_{collibri}_R1_001.fastq.gz").collibri
KAPA_FILES = glob_wildcards("data/KAPA_{kapa}_R1_001.fastq.gz").kapa

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
        expand("data_output_samtools_indexed/{trimmed}.bam.bai", trimmed=TRIMMED_FILES),
        expand("data_output_featureCount/{trimmed}_s1.txt", trimmed=TRIMMED_FILES),
        expand("data_output_featureCount/{trimmed}_s2.txt", trimmed=TRIMMED_FILES),
        expand("data_output_selected/{trimmed}.txt", trimmed=TRIMMED_FILES),
        "data_output_matrix/Collibri.txt",
        "data_output_matrix/KAPA.txt"         

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

rule featureCounts s1:
    input:
        "data_output_star_mapped/{trimmed}_Aligned.sortedByCoord.out.bam"
    output:
        "data_output_featureCount/{trimmed}_s1.txt"    
    shell:
        "featureCounts -p -t exon -g gene_id -O -T 8 -a play_data_ref_annot/chr19_20Mb.gtf -o {output} {input} -s 1" 

rule featureCounts s2:
    input:
        "data_output_star_mapped/{trimmed}_Aligned.sortedByCoord.out.bam"
    output:
        "data_output_featureCount/{trimmed}_s2.txt"    
    shell:
        "featureCounts -p -t exon -g gene_id -O -T 8 -a play_data_ref_annot/chr19_20Mb.gtf -o {output} {input} -s 2" 

rule select_s1_or_s2:
    input:
        "data_output_featureCount/{trimmed}_s1.txt",
        "data_output_featureCount/{trimmed}_s2.txt"
    output:
        "data_output_selected/{trimmed}.txt"    
    shell:
        "python3 select_best.py {input} {output}" 

rule collibri_matrix:
    input:
        expand("data_output_selected/Collibri_{collibri}.txt", collibri=COLLIBRI_FILES)
    output:
        "data_output_matrix/Collibri.txt"    
    shell:
        "python3 count_matrix.py {input} {output}" 

rule kapa_matrix:
    input:
        expand("data_output_selected/KAPA_{kapa}.txt", kapa=KAPA_FILES)
    output:
        "data_output_matrix/KAPA.txt"    
    shell:
        "python3 count_matrix.py {input} {output}" 



# KAPA_mRNA_HyperPrep_-UHRR-KAPA-100_ng_total_RNA-2_S7_L001_s2.bam 
# Assigned	79965

# Prainkti viena is dvieju ir vienam meginiui visus failu duomenis ir padarytu countu eilute. Y - stulpeliai meginys arba failas, o eilute - genai. Selection rule, kad pasirinktu viena is s1 s2





