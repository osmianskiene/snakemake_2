```
conda activate
conda deactivate
conda activate snakemake-tutorial
cd snakemake_2
conda list

conda list fastqc
conda install fastqc
fastqc --version
fastqc --help
$ fastqc -o data_output data/Collibri_standard_protocol-HBR-Collibri-100_ng-2_S1_L001_R1_001.fastq.gz

/home/mn/snakemake_2/data/Collibri_standard_protocol-HBR-Collibri-100_ng-2_S1_L001_R1_001.fastq.gz

snakemake --cores 3

conda install multiqc
multiqc  --version
multiqc  --help

cd data_output
multiqc .
multiqc data_output -o data_output_multiqc
cd ..
snakemake --cores 3
 snakemake multiqc  --cores 3

conda install bbmap
bbduk.sh --version
snakemake --cores 3

snakemake -R --dag samtools_sort
snakemake -nR samtools_sort
snakemake -F
snakemake --dag | dot -Tsvg > dag_temp.svg
snakemake --help > snakemake_help.txt
snakemake plots/quals.svg
tar --wildcards -xf snakemake-tutorial-data.tar.gz --strip 1 "*/data" "*/environment.yaml"

touch data/samples/{A,B,C}.fastq

pip install WeasyPrint

echo $PATH
```
1. Vykdyk sia komanda tik jeigu kazkas paprase output filo ir jo nera.
2. Tada, jei vykdai isitikink, kad yra input failas
3.  Jei nera - surask kita komanda, kuri generuoja ir paleisk ja.
4.  Jei yra - tai vykdyk shell komandą.
