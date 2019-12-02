# Practice of NGS in _Drosophila_
Global changes of H3K27me3 domains and Polycomb group protein distribution in the absence of recruiters Spps or Pho
(https://www.ncbi.nlm.nih.gov/pubmed/29432187)
## Preprocess
#### Raw data download
```bash
conda activate Download
cat Fly/SRR_Acc_List.txt | while read id; do prefetch $id -O Fly/Raw & done
cp -r Fly/Raw/*/*.sra Fly/Data
ls Fly/Data/* | while read id; do fastq-dump --gzip --split-3 -O Fly/Data/ ${id} & done
rm Fly/Data/*.sra
conda deactivate
```
#### FastQC and multiQC of raw data
```bash
conda activate QC
ls Fly/Data/*fastq.gz |  while read id; do fastqc -f fastq -o Fly/FastQC_1/./ ${id} & done
multiqc Fly/FastQC_1 -o Fly/FastQC_1
conda deactivate
```
#### Build Index
```bash
hisat2-build Fly/Drosophila.fa Fly/Hisat2_fly/fly
bowtie2-build Fly/Drosophila.fa Fly/Bowtie2_fly/fly
```
## RNA-Seq
#### Hisat2
```bash
hisat2 -t -p 8 -x Fly/Hisat2_fly/fly -1 Fly/Data/SRR5907469_1.fastq.gz -2 Fly/Data/SRR5907469_2.fastq.gz | samtools sort -O bam -o Fly/RNA_seq/PhoKO_1.bam
hisat2 -t -p 8 -x Fly/Hisat2_fly/fly -1 Fly/Data/SRR5907470_1.fastq.gz -2 Fly/Data/SRR5907470_2.fastq.gz | samtools sort -O bam -o Fly/RNA_seq/PhoKO_2.bam
hisat2 -t -p 8 -x Fly/Hisat2_fly/fly -1 Fly/Data/SRR5907471_1.fastq.gz -2 Fly/Data/SRR5907471_2.fastq.gz | samtools sort -O bam -o Fly/RNA_seq/PhoKO_3.bam
hisat2 -t -p 8 -x Fly/Hisat2_fly/fly -1 Fly/Data/SRR5907472_1.fastq.gz -2 Fly/Data/SRR5907472_2.fastq.gz | samtools sort -O bam -o Fly/RNA_seq/SppsKO_1.bam
hisat2 -t -p 8 -x Fly/Hisat2_fly/fly -1 Fly/Data/SRR5907473_1.fastq.gz -2 Fly/Data/SRR5907473_2.fastq.gz | samtools sort -O bam -o Fly/RNA_seq/SppsKO_2.bam
hisat2 -t -p 8 -x Fly/Hisat2_fly/fly -1 Fly/Data/SRR5907474_1.fastq.gz -2 Fly/Data/SRR5907474_2.fastq.gz | samtools sort -O bam -o Fly/RNA_seq/SppsKO_3.bam
hisat2 -t -p 8 -x Fly/Hisat2_fly/fly -1 Fly/Data/SRR5907475_1.fastq.gz -2 Fly/Data/SRR5907475_2.fastq.gz | samtools sort -O bam -o Fly/RNA_seq/SppsKO_4.bam
hisat2 -t -p 8 -x Fly/Hisat2_fly/fly -1 Fly/Data/SRR5907476_1.fastq.gz -2 Fly/Data/SRR5907476_2.fastq.gz | samtools sort -O bam -o Fly/RNA_seq/WT_1.bam
hisat2 -t -p 8 -x Fly/Hisat2_fly/fly -1 Fly/Data/SRR5907477_1.fastq.gz -2 Fly/Data/SRR5907477_2.fastq.gz | samtools sort -O bam -o Fly/RNA_seq/WT_2.bam
hisat2 -t -p 8 -x Fly/Hisat2_fly/fly -1 Fly/Data/SRR5907478_1.fastq.gz -2 Fly/Data/SRR5907478_2.fastq.gz | samtools sort -O bam -o Fly/RNA_seq/WT_3.bam
```
#### FeatureCounts
```bash
featureCounts -T 8 -p -t exon -g gene_name -a Fly/Drosophila.gtf -o PhoKO_1.txt Fly/RNA_seq/PhoKO_1.bam
featureCounts -T 8 -p -t exon -g gene_name -a Fly/Drosophila.gtf -o PhoKO_2.txt Fly/RNA_seq/PhoKO_2.bam
featureCounts -T 8 -p -t exon -g gene_name -a Fly/Drosophila.gtf -o PhoKO_3.txt Fly/RNA_seq/PhoKO_3.bam
featureCounts -T 8 -p -t exon -g gene_name -a Fly/Drosophila.gtf -o SppsKO_1.txt Fly/RNA_seq/SppsKO_1.bam
featureCounts -T 8 -p -t exon -g gene_name -a Fly/Drosophila.gtf -o SppsKO_2.txt Fly/RNA_seq/SppsKO_2.bam
featureCounts -T 8 -p -t exon -g gene_name -a Fly/Drosophila.gtf -o SppsKO_3.txt Fly/RNA_seq/SppsKO_3.bam
featureCounts -T 8 -p -t exon -g gene_name -a Fly/Drosophila.gtf -o SppsKO_4.txt Fly/RNA_seq/SppsKO_4.bam
featureCounts -T 8 -p -t exon -g gene_name -a Fly/Drosophila.gtf -o WT_1.txt Fly/RNA_seq/WT_1.bam
featureCounts -T 8 -p -t exon -g gene_name -a Fly/Drosophila.gtf -o WT_2.txt Fly/RNA_seq/WT_2.bam
featureCounts -T 8 -p -t exon -g gene_name -a Fly/Drosophila.gtf -o WT_3.txt Fly/RNA_seq/WT_3.bam
cut -f 1,7 PhoKO_1.txt |grep -v '^#' >Fly/Count/PhoKO_1.txt
cut -f 1,7 PhoKO_2.txt |grep -v '^#' >Fly/Count/PhoKO_2.txt
cut -f 1,7 PhoKO_3.txt |grep -v '^#' >Fly/Count/PhoKO_3.txt
cut -f 1,7 SppsKO_1.txt |grep -v '^#' >Fly/Count/SppsKO_1.txt
cut -f 1,7 SppsKO_2.txt |grep -v '^#' >Fly/Count/SppsKO_2.txt
cut -f 1,7 SppsKO_3.txt |grep -v '^#' >Fly/Count/SppsKO_3.txt
cut -f 1,7 SppsKO_4.txt |grep -v '^#' >Fly/Count/SppsKO_4.txt
cut -f 1,7 WT_1.txt |grep -v '^#' >Fly/Count/WT_1.txt
cut -f 1,7 WT_2.txt |grep -v '^#' >Fly/Count/WT_2.txt
cut -f 1,7 WT_3.txt |grep -v '^#' >Fly/Count/WT_3.txt
paste Fly/Count/WT_1.txt Fly/Count/WT_2.txt Fly/Count/WT_3.txt Fly/Count/PhoKO_1.txt Fly/Count/PhoKO_2.txt Fly/Count/PhoKO_3.txt | awk '{printf $1"\t";for(i=2;i<=NF;i=i+2)printf $i"\t";print $i}' > PhoKO.txt
paste Fly/Count/WT_1.txt Fly/Count/WT_2.txt Fly/Count/WT_3.txt Fly/Count/SppsKO_1.txt Fly/Count/SppsKO_2.txt Fly/Count/SppsKO_3.txt Fly/Count/SppsKO_4.txt | awk '{printf $1"\t";for(i=2;i<=NF;i=i+2)printf $i"\t";print $i}' > SppsKO.txt
paste Fly/Count/WT_1.txt Fly/Count/WT_2.txt Fly/Count/WT_3.txt Fly/Count/SppsKO_1.txt Fly/Count/SppsKO_2.txt Fly/Count/SppsKO_3.txt Fly/Count/SppsKO_4.txt Fly/Count/PhoKO_1.txt Fly/Count/PhoKO_2.txt Fly/Count/PhoKO_3.txt | awk '{printf $1"\t";for(i=2;i<=NF;i=i+2)printf $i"\t";print $i}' > Matrix.txt
```

## ChIP-Seq
