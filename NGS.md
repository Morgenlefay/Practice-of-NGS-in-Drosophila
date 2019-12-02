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
#### Bam2Bigwig
```bash
ls Fly/RNA_seq/*.bam |  while read id; do samtools index $id & done
ls Fly/RNA_seq/*.bam |  while read id; do bamCoverage -b $id -o $id.bw --normalizeUsing RPKM -p 1 & done
```
#### Differential gene expression analysis
##### Pho KO
```r
###--------------------------
rm(list=ls())
options(stringsAsFactors = F)
library(DESeq2)

###构建表达矩阵
Mydata <- read.table("PhoKO.txt", header=TRUE)
row.names(Mydata) <- Mydata$Geneid
CountData <- Mydata[ ,-1]
condition <- factor(c("WT","WT","WT","PhoKO","PhoKO","PhoKO"))
colData <- data.frame(row.names=colnames(CountData), condition)

###DEseq标准化dds
dds <- DESeqDataSetFromMatrix(countData=CountData, colData=colData, design=~condition)
dds <- DESeq(dds)
res <- results(dds)
mcols(res, use.names = TRUE)
summary(res)
plotMA(res, ylim = c(-6,6))
topGene <- rownames(res)[which.min(res$padj)]
with(res[topGene, ], {
  points(baseMean, log2FoldChange, col="dodgerblue", cex=2, lwd=2)
  text(baseMean, log2FoldChange, topGene, pos=2, col="dodgerblue")
})

###提取差异分析结果
res <- res[order(res$padj),]
diff_gene_deseq2 <-subset(res,padj < 0.05 & (log2FoldChange > 1 | log2FoldChange < -1))
diff_gene_deseq2 <- row.names(diff_gene_deseq2)
resdata <-  merge(as.data.frame(res),as.data.frame(counts(dds,normalize=TRUE)),by="row.names",sort=FALSE)
write.csv(resdata,file= "PhoKO_DEG.csv",row.names = F)
subset(res,padj < 0.05) -> diff
subset(diff,log2FoldChange < -1) -> down
subset(diff,log2FoldChange > 1) -> up
as.data.frame(down) -> down_gene
as.data.frame(up) -> up_gene
write.csv(up_gene, file="PhoKO_Up.csv",row.names = T)
write.csv(down_gene, file="PhoKO_Down.csv",row.names = T)
```

## ChIP-Seq
#### Bowtie2
```bash
bowtie2 -t -p 8 -x Fly/Bowtie2_fly/fly -U Fly/Data/SRR5907429.fastq.gz | samtools sort -O bam -o Fly/ChIP_seq/Ez_WT_1.bam
bowtie2 -t -p 8 -x Fly/Bowtie2_fly/fly -U Fly/Data/SRR5907430.fastq.gz | samtools sort -O bam -o Fly/ChIP_seq/Ez_WT_2.bam
bowtie2 -t -p 8 -x Fly/Bowtie2_fly/fly -U Fly/Data/SRR5907431.fastq.gz | samtools sort -O bam -o Fly/ChIP_seq/Ez_WT_3.bam
bowtie2 -t -p 8 -x Fly/Bowtie2_fly/fly -U Fly/Data/SRR5907432.fastq.gz | samtools sort -O bam -o Fly/ChIP_seq/Pc_WT_1.bam
bowtie2 -t -p 8 -x Fly/Bowtie2_fly/fly -U Fly/Data/SRR5907433.fastq.gz | samtools sort -O bam -o Fly/ChIP_seq/Pc_WT_2.bam
bowtie2 -t -p 8 -x Fly/Bowtie2_fly/fly -U Fly/Data/SRR5907434.fastq.gz | samtools sort -O bam -o Fly/ChIP_seq/Ph_WT_1.bam
bowtie2 -t -p 8 -x Fly/Bowtie2_fly/fly -U Fly/Data/SRR6490544.fastq.gz | samtools sort -O bam -o Fly/ChIP_seq/Ph_WT_2.bam
bowtie2 -t -p 8 -x Fly/Bowtie2_fly/fly -U Fly/Data/SRR5907436.fastq.gz | samtools sort -O bam -o Fly/ChIP_seq/Pho_WT_1.bam
bowtie2 -t -p 8 -x Fly/Bowtie2_fly/fly -U Fly/Data/SRR5907437.fastq.gz | samtools sort -O bam -o Fly/ChIP_seq/Pho_WT_2.bam
bowtie2 -t -p 8 -x Fly/Bowtie2_fly/fly -U Fly/Data/SRR5907438.fastq.gz | samtools sort -O bam -o Fly/ChIP_seq/Pho_WT_3.bam
bowtie2 -t -p 8 -x Fly/Bowtie2_fly/fly -U Fly/Data/SRR5907439.fastq.gz | samtools sort -O bam -o Fly/ChIP_seq/Psc_WT_1.bam
bowtie2 -t -p 8 -x Fly/Bowtie2_fly/fly -U Fly/Data/SRR5907440.fastq.gz | samtools sort -O bam -o Fly/ChIP_seq/Psc_WT_2.bam
bowtie2 -t -p 8 -x Fly/Bowtie2_fly/fly -U Fly/Data/SRR5907441.fastq.gz | samtools sort -O bam -o Fly/ChIP_seq/Spps_WT_1.bam
bowtie2 -t -p 8 -x Fly/Bowtie2_fly/fly -U Fly/Data/SRR5907442.fastq.gz | samtools sort -O bam -o Fly/ChIP_seq/Spps_WT_2.bam
bowtie2 -t -p 8 -x Fly/Bowtie2_fly/fly -U Fly/Data/SRR5907443.fastq.gz | samtools sort -O bam -o Fly/ChIP_seq/Input_WT_1.bam
bowtie2 -t -p 8 -x Fly/Bowtie2_fly/fly -U Fly/Data/SRR5907444.fastq.gz | samtools sort -O bam -o Fly/ChIP_seq/Input_WT_2.bam
bowtie2 -t -p 8 -x Fly/Bowtie2_fly/fly -U Fly/Data/SRR5907445.fastq.gz | samtools sort -O bam -o Fly/ChIP_seq/H3K27_SppsKO_1.bam
bowtie2 -t -p 8 -x Fly/Bowtie2_fly/fly -U Fly/Data/SRR5907446.fastq.gz | samtools sort -O bam -o Fly/ChIP_seq/H3K27_SppsKO_2.bam
bowtie2 -t -p 8 -x Fly/Bowtie2_fly/fly -U Fly/Data/SRR5907447.fastq.gz | samtools sort -O bam -o Fly/ChIP_seq/H3K27_SppsKO_3.bam
bowtie2 -t -p 8 -x Fly/Bowtie2_fly/fly -U Fly/Data/SRR5907448.fastq.gz | samtools sort -O bam -o Fly/ChIP_seq/Ez_SppsKO_1.bam
bowtie2 -t -p 8 -x Fly/Bowtie2_fly/fly -U Fly/Data/SRR5907449.fastq.gz | samtools sort -O bam -o Fly/ChIP_seq/Ez_SppsKO_2.bam
bowtie2 -t -p 8 -x Fly/Bowtie2_fly/fly -U Fly/Data/SRR5907450.fastq.gz | samtools sort -O bam -o Fly/ChIP_seq/Ez_SppsKO_3.bam
bowtie2 -t -p 8 -x Fly/Bowtie2_fly/fly -U Fly/Data/SRR5907451.fastq.gz | samtools sort -O bam -o Fly/ChIP_seq/Ph_SppsKO_1.bam
bowtie2 -t -p 8 -x Fly/Bowtie2_fly/fly -U Fly/Data/SRR5907452.fastq.gz | samtools sort -O bam -o Fly/ChIP_seq/Ph_SppsKO_2.bam
bowtie2 -t -p 8 -x Fly/Bowtie2_fly/fly -U Fly/Data/SRR5907423.fastq.gz | samtools sort -O bam -o Fly/ChIP_seq/Ph_SppsKO_3.bam
bowtie2 -t -p 8 -x Fly/Bowtie2_fly/fly -U Fly/Data/SRR5907454.fastq.gz | samtools sort -O bam -o Fly/ChIP_seq/Pho_SppsKO_1.bam
bowtie2 -t -p 8 -x Fly/Bowtie2_fly/fly -U Fly/Data/SRR5907455.fastq.gz | samtools sort -O bam -o Fly/ChIP_seq/Pho_SppsKO_2.bam
bowtie2 -t -p 8 -x Fly/Bowtie2_fly/fly -U Fly/Data/SRR5907456.fastq.gz | samtools sort -O bam -o Fly/ChIP_seq/Psc_SppsKO_1.bam
bowtie2 -t -p 8 -x Fly/Bowtie2_fly/fly -U Fly/Data/SRR5907457.fastq.gz | samtools sort -O bam -o Fly/ChIP_seq/Psc_SppsKO_2.bam
bowtie2 -t -p 8 -x Fly/Bowtie2_fly/fly -U Fly/Data/SRR5907458.fastq.gz | samtools sort -O bam -o Fly/ChIP_seq/Psc_SppsKO_3.bam
bowtie2 -t -p 8 -x Fly/Bowtie2_fly/fly -U Fly/Data/SRR5907459.fastq.gz | samtools sort -O bam -o Fly/ChIP_seq/Input_SppsKO_1.bam
bowtie2 -t -p 8 -x Fly/Bowtie2_fly/fly -U Fly/Data/SRR5907460.fastq.gz | samtools sort -O bam -o Fly/ChIP_seq/Input_SppsKO_2.bam
bowtie2 -t -p 8 -x Fly/Bowtie2_fly/fly -U Fly/Data/SRR5907461.fastq.gz | samtools sort -O bam -o Fly/ChIP_seq/H3K27_PhoKO_1.bam
bowtie2 -t -p 8 -x Fly/Bowtie2_fly/fly -U Fly/Data/SRR5907462.fastq.gz | samtools sort -O bam -o Fly/ChIP_seq/H3K27_PhoKO_2.bam
bowtie2 -t -p 8 -x Fly/Bowtie2_fly/fly -U Fly/Data/SRR5907463.fastq.gz | samtools sort -O bam -o Fly/ChIP_seq/H3K27_PhoKO_3.bam
bowtie2 -t -p 8 -x Fly/Bowtie2_fly/fly -U Fly/Data/SRR5907464.fastq.gz | samtools sort -O bam -o Fly/ChIP_seq/Ph_PhoKO_1.bam
bowtie2 -t -p 8 -x Fly/Bowtie2_fly/fly -U Fly/Data/SRR5907465.fastq.gz | samtools sort -O bam -o Fly/ChIP_seq/Ph_PhoKO_2.bam
bowtie2 -t -p 8 -x Fly/Bowtie2_fly/fly -U Fly/Data/SRR5907466.fastq.gz | samtools sort -O bam -o Fly/ChIP_seq/Ph_PhoKO_3.bam
bowtie2 -t -p 8 -x Fly/Bowtie2_fly/fly -U Fly/Data/SRR5907467.fastq.gz | samtools sort -O bam -o Fly/ChIP_seq/Input_PhoKO_1.bam
bowtie2 -t -p 8 -x Fly/Bowtie2_fly/fly -U Fly/Data/SRR5907468.fastq.gz | samtools sort -O bam -o Fly/ChIP_seq/Input_PhoKO_2.bam
bowtie2 -t -p 8 -x Fly/Bowtie2_fly/fly -U Fly/Data/SRR9967698.fastq.gz | samtools sort -O bam -o Fly/ChIP_seq/Input_PhoKO_2.bam
```
#### Deduplicates
```bash

```


