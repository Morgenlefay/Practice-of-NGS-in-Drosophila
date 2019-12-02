### Raw data download
```bash
conda activate Download
cat Fly/SRR_Acc_List.txt | while read id; do prefetch $id -O Fly/Raw & done
cp -r Fly/Raw/*/*.sra Fly/Data
ls Fly/Data/* | while read id; do fastq-dump --gzip --split-3 -O Fly/Data/ ${id} & done
rm Fly/Data/*.sra
conda deactivate
```
### FastQC and multiQC of raw data
```bash
fastqc=/home/morgenlefay/miniconda2/pkgs/fastqc-0.11.8-1/bin/fastqc
ls Fly/Data/*fastq.gz |  while read id; do 
$fastqc -f fastq -o Fly/FastQC_1/./ ${id} &
done
multiqc Fly/FastQC_1 -o Fly/FastQC_1
```

### Build Index
```bash
hisat2-build Fly/Drosophila.fa Fly/Hisat2_fly/fly
bowtie2-build Fly/Drosophila.fa Fly/Bowtie2_fly/fly
```


