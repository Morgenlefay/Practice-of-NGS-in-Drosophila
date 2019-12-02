# Practice of NGS in _Drosophila_
Global changes of H3K27me3 domains and Polycomb group protein distribution in the absence of recruiters Spps or Pho
(https://www.ncbi.nlm.nih.gov/pubmed/29432187)
## Preprocess
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
conda activate QC
ls Fly/Data/*fastq.gz |  while read id; do fastqc -f fastq -o Fly/FastQC_1/./ ${id} & done
multiqc Fly/FastQC_1 -o Fly/FastQC_1
conda deactivate
```
### Build Index
```bash
hisat2-build Fly/Drosophila.fa Fly/Hisat2_fly/fly
bowtie2-build Fly/Drosophila.fa Fly/Bowtie2_fly/fly
```
## RNA-Seq



## ChIP-Seq
