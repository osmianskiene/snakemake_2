```
conda activate snakemake-tutorial
cd snakemake_2



conda install bbmap

conda list

conda activate
conda deactivate

cd ..

snakemake --cores 3
snakemake multiqc  --cores 3
snakemake -R --dag samtools_sort
snakemake -nR samtools_sort
snakemake -F
snakemake --dag | dot -Tsvg > dag_temp.svg
snakemake --help > snakemake_help.txt
snakemake plots/quals.svg
tar --wildcards -xf snakemake-tutorial-data.tar.gz --strip 1 "*/data" "*/environment.yaml"

touch data/samples/{A,B,C}.fastq

multiqc .
multiqc  --version

bbduk.sh --version




pip install WeasyPrint

echo $PATH
```