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
bowtie2-build Fly/dm6.fa Fly/Bowtie2_fly/fly
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

### shrinkage estimators
resNorm <- lfcShrink(dds, coef=2, type="normal")
par(mfrow=c(1,1))
xlim <- c(1,1e5); ylim <- c(-3,3)
plotMA(resNorm, xlim=xlim, ylim=ylim, main="normal")

### 表达矩阵Normalization
rld <- rlogTransformation(dds)
exprSet=assay(rld)
par(cex = 0.7)
n.sample=ncol(exprSet)
if(n.sample>40) par(cex = 0.5)
cols <- rainbow(n.sample*1.2)
par(mfrow=c(2,2))
boxplot(exprSet, col = cols,main="expression value",las=2)
hist(exprSet)

###提取差异分析结果
res <- res[order(res$padj),]
diff_gene_deseq2 <-subset(res,padj < 0.05)
diff_gene_deseq2 <- row.names(diff_gene_deseq2)
resdata <-  merge(as.data.frame(res),as.data.frame(counts(dds,normalize=TRUE)),by="row.names",sort=FALSE)
write.csv(resdata,file= "PhoKO_DEG.csv",row.names = F)
subset(res,padj < 0.05) -> diff
subset(diff,log2FoldChange < 0) -> down
subset(diff,log2FoldChange > 0) -> up
as.data.frame(down) -> down_gene
as.data.frame(up) -> up_gene
write.csv(up_gene, file="PhoKO_Up.csv",row.names = T)
write.csv(down_gene, file="PhoKO_Down.csv",row.names = T)
```
##### Spps KO
```bash
###--------------------------
rm(list=ls())
options(stringsAsFactors = F)
library(DESeq2)

###构建表达矩阵
Mydata <- read.table("SppsKO.txt", header=TRUE)
row.names(Mydata) <- Mydata$Geneid
CountData <- Mydata[ ,-1]
condition <- factor(c("WT","WT","WT","SppsKO","SppsKO","SppsKO","SppsKO"))
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

### shrinkage estimators
resNorm <- lfcShrink(dds, coef=2, type="normal")
par(mfrow=c(1,1))
xlim <- c(1,1e5); ylim <- c(-3,3)
plotMA(resNorm, xlim=xlim, ylim=ylim, main="normal")

### 表达矩阵Normalization
rld <- rlogTransformation(dds)
exprSet=assay(rld)
par(cex = 0.7)
n.sample=ncol(exprSet)
if(n.sample>40) par(cex = 0.5)
cols <- rainbow(n.sample*1.2)
par(mfrow=c(2,2))
boxplot(exprSet, col = cols,main="expression value",las=2)
hist(exprSet)

###提取差异分析结果
res <- res[order(res$padj),]
diff_gene_deseq2 <-subset(res,padj < 0.05)
diff_gene_deseq2 <- row.names(diff_gene_deseq2)
resdata <-  merge(as.data.frame(res),as.data.frame(counts(dds,normalize=TRUE)),by="row.names",sort=FALSE)
write.csv(resdata,file= "SppsKO_DEG.csv",row.names = F)
subset(res,padj < 0.05) -> diff
subset(diff,log2FoldChange < 0) -> down
subset(diff,log2FoldChange > 0) -> up
as.data.frame(down) -> down_gene
as.data.frame(up) -> up_gene
write.csv(up_gene, file="SppsKO_Up.csv",row.names = T)
write.csv(down_gene, file="SppsKO_Down.csv",row.names = T)
```
#### Correlation
##### Expression
```r
###--------------------------
rm(list=ls())
options(stringsAsFactors = F)
library(pheatmap)
Data <- read.table("Matrix.txt", header=TRUE)
row.names(Data) <- Data$Geneid
Data <- Data[ ,-1]
library(stringr)
ac=data.frame(group=str_split(colnames(Data),'_',simplify = T)[,1])
rownames(ac)=colnames(Data)
pheatmap::pheatmap(cor(log(Data+1)),annotation_col = ac)
cg=Data[,colnames(Data)[grepl('_1',colnames(Data))]]
library(ggpubr)
pairs(log(cg+1))
```
##### log2FC
```r
###--------------------------
rm(list=ls())
options(stringsAsFactors = F)
library(ggpubr)
library(plyr)
PhoKO_DEG <- read.csv("PhoKO_DEG.csv",header = TRUE)
SppsKO_DEG <-read.csv("SppsKO_DEG.csv",header = TRUE)
PhoKO_DEG <- PhoKO_DEG[,1:3]
SppsKO_DEG <- SppsKO_DEG[,1:3]
PhoKO_DEG <- PhoKO_DEG[,-2]
SppsKO_DEG <- SppsKO_DEG[,-2]
PhoKO_DEG <- PhoKO_DEG[order(PhoKO_DEG[,1]),]
SppsKO_DEG <- SppsKO_DEG[order(SppsKO_DEG[,1]),]
df=data.frame(SppsKO=SppsKO_DEG$log2FoldChange,
              PhoKO=PhoKO_DEG$log2FoldChange)
sp <- ggscatter(df, "SppsKO","PhoKO",
                add = "reg.line",  # Add regressin line
                add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                conf.int = TRUE # Add confidence interval
                )
# Add correlation coefficient
sp + stat_cor(method = "pearson", label.x = 3, label.y = 30)
```
#### Expression
```r
###--------------------------
rm(list = ls())
options(stringsAsFactors = F)
library(ggpubr)
library(stringr)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
Data <- read.table("Matrix.txt", header=TRUE)
## 定义Barplot主题
my_bar_theme <- theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 16), 
                      # 因为 x 轴标签要旋转 90°
                      axis.text.y = element_text(size = 16),
                      axis.title.y = element_text(size = 18,
                                                  face = "bold",),
                      legend.position = "none") +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5)) # 基因名居中


## 提取 `Pho` 基因对应的行以及表达量信息
cg <- Data[Data[,1] == 'pho',2:11]
dat <- data.frame(gene = as.numeric(cg),
                  sample = str_split(names(cg),'\\.',simplify = T)[,1],
                  group = str_split(names(cg),'_',simplify = T)[,1]
)
labels = c(paste0("WT", "_",1:3), paste0("SppsKO", "_", 1:4), paste0("PhoKO", "_", 1:3))

# 使用ggbarplot()函数进行绘图
ggbarplot(dat,x = 'sample', y = 'gene', 
          color = 'black', fill = "group", 
          size = 1) +
  # 使用 Takecolor 软件进行对准样图取色
  scale_fill_manual(values = c(WT = "#9C4B25", 
                               SppsKO = "#DE2C1C", 
                               PhoKO = "#43A542")) +
  scale_color_manual(values = "black") +
  scale_x_discrete(limits = labels) +
  labs(y = "Normalized read count",
       x = "",
       title = "Pho") +
  my_bar_theme +
  scale_y_continuous(expand = c(0, 0)) + 
  geom_text(aes(y = gene * 1.1, label = "")) 

## 封装画图函数
my_barplot <- function(gene){
  cg <- Data[Data[,1] == gene, 2:11] # 提取候选的表达量对应的行和列
  
  dat <- data.frame(gene = as.numeric(cg),
                    sample = str_split(names(cg),'\\.',simplify = T)[,1], # 这里将样品后面的 bam文件后缀去掉
                    group = str_split(names(cg),'_',simplify = T)[,1]
  )
  
  # x 轴标签的顺序，这里是按照原图的顺序来的
  labels = c(paste0("WT", "_",1:3), paste0("SppsKO", "_", 1:4), paste0("PhoKO", "_", 1:3))
  
  ggbarplot(dat, x = 'sample', y = 'gene', 
            color = 'black', fill = "group", 
            size = 1) +
    scale_fill_manual(values = c(WT = "#9C4B25", 
                                 SppsKO = "#DE2C1C", 
                                 PhoKO = "#43A542")) +
    scale_color_manual(values = "black") +
    scale_x_discrete(limits = labels) +
    labs(y = "Normalized read count",
         x = "",
         title = gene) +
    scale_y_continuous(expand = c(0, 0)) + 
    geom_text(aes(y = gene * 1.1, label = "")) +
    my_bar_theme 
}


## 多个基因组图
Spps <- Data[, 1][3045]
Pho <- Data[, 1][17671]
Spps_plot <- my_barplot("Spps")
Pho_plot <- my_barplot("Pho")
gene_exp_plot <- plot_grid(Pho_plot, Spps_plot, 
                           labels = c("A", "B"))
ggsave("gene_exp_plot.pdf", gene_exp_plot, width = 10, height = 5) 

## 批量多个基因组图
gene_list <- Data[, 1][1:10]
test <- lapply(gene_list, my_barplot)
gene <- plot_grid(plotlist = test, ncol = 5)
ggsave("gene.pdf", gene, width = 20, height = 5* length(test) %/% 5)
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
bowtie2 -t -p 8 -x Fly/Bowtie2_fly/fly -U Fly/Data/SRR5907453.fastq.gz | samtools sort -O bam -o Fly/ChIP_seq/Ph_SppsKO_3.bam
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
bowtie2 -t -p 8 -x Fly/Bowtie2_fly/fly -U Fly/Data/SRR9967698.fastq.gz | samtools sort -O bam -o Fly/ChIP_seq/Input_PhoKO_3.bam
```
#### Deduplicates
```bash
ls Fly/ChIP_seq/*.bam |  while read id; do sambamba markdup -r $id $id.dedup.bam & done
```
#### Merge
```bash
samtools merge Fly/ChIP_seq/Merge/Ez_SppsKO.bam Fly/ChIP_seq/Ez_SppsKO_1.bam.dedup.bam  Fly/ChIP_seq/Ez_SppsKO_2.bam.dedup.bam Fly/ChIP_seq/Ez_SppsKO_3.bam.dedup.bam
samtools merge Fly/ChIP_seq/Merge/Ez_WT.bam Fly/ChIP_seq/Ez_WT_1.bam.dedup.bam Fly/ChIP_seq/Ez_WT_2.bam.dedup.bam Fly/ChIP_seq/Ez_WT_3.bam.dedup.bam
samtools merge Fly/ChIP_seq/Merge/H3K27_PhoKO.bam Fly/ChIP_seq/H3K27_PhoKO_1.bam.dedup.bam Fly/ChIP_seq/H3K27_PhoKO_2.bam.dedup.bam Fly/ChIP_seq/H3K27_PhoKO_3.bam.dedup.bam
samtools merge Fly/ChIP_seq/Merge/H3K27_SppsKO.bam Fly/ChIP_seq/H3K27_SppsKO_1.bam.dedup.bam Fly/ChIP_seq/H3K27_SppsKO_2.bam.dedup.bam Fly/ChIP_seq/H3K27_SppsKO_3.bam.dedup.bam
samtools merge Fly/ChIP_seq/Merge/Input_PhoKO.bam Fly/ChIP_seq/Input_PhoKO_1.bam.dedup.bam Fly/ChIP_seq/Input_PhoKO_2.bam.dedup.bam Fly/ChIP_seq/Input_PhoKO_3.bam.dedup.bam
samtools merge Fly/ChIP_seq/Merge/Input_SppsKO.bam Fly/ChIP_seq/Input_SppsKO_1.bam.dedup.bam Fly/ChIP_seq/Input_SppsKO_2.bam.dedup.bam
samtools merge Fly/ChIP_seq/Merge/Input_WT.bam Fly/ChIP_seq/Input_WT_1.bam.dedup.bam Fly/ChIP_seq/Input_WT_2.bam.dedup.bam
samtools merge Fly/ChIP_seq/Merge/Pc_WT.bam Fly/ChIP_seq/Pc_WT_1.bam.dedup.bam Fly/ChIP_seq/Pc_WT_2.bam.dedup.bam
samtools merge Fly/ChIP_seq/Merge/Ph_PhoKO.bam Fly/ChIP_seq/Ph_PhoKO_1.bam.dedup.bam Fly/ChIP_seq/Ph_PhoKO_2.bam.dedup.bam Fly/ChIP_seq/Ph_PhoKO_3.bam.dedup.bam
samtools merge Fly/ChIP_seq/Merge/Ph_SppsKO.bam Fly/ChIP_seq/Ph_SppsKO_1.bam.dedup.bam Fly/ChIP_seq/Ph_SppsKO_2.bam.dedup.bam Fly/ChIP_seq/Ph_SppsKO_3.bam.dedup.bam
samtools merge Fly/ChIP_seq/Merge/Ph_WT.bam Fly/ChIP_seq/Ph_WT_1.bam.dedup.bam Fly/ChIP_seq/Ph_WT_2.bam.dedup.bam
samtools merge Fly/ChIP_seq/Merge/Pho_SppsKO.bam Fly/ChIP_seq/Pho_SppsKO_1.bam.dedup.bam Fly/ChIP_seq/Pho_SppsKO_2.bam.dedup.bam
samtools merge Fly/ChIP_seq/Merge/Pho_WT.bam Fly/ChIP_seq/Pho_WT_1.bam.dedup.bam Fly/ChIP_seq/Pho_WT_2.bam.dedup.bam Fly/ChIP_seq/Pho_WT_3.bam.dedup.bam
samtools merge Fly/ChIP_seq/Merge/Psc_SppsKO.bam Fly/ChIP_seq/Psc_SppsKO_1.bam.dedup.bam Fly/ChIP_seq/Psc_SppsKO_2.bam.dedup.bam Fly/ChIP_seq/Psc_SppsKO_3.bam.dedup.bam
samtools merge Fly/ChIP_seq/Merge/Psc_WT.bam Fly/ChIP_seq/Psc_WT_1.bam.dedup.bam Fly/ChIP_seq/Psc_WT_2.bam.dedup.bam
samtools merge Fly/ChIP_seq/Merge/Spps_WT.bam Fly/ChIP_seq/Spps_WT_1.bam.dedup.bam Fly/ChIP_seq/Spps_WT_2.bam.dedup.bam
```
#### Peak calling
```bash
macs2 callpeak -t Fly/ChIP_seq/Ez_WT.bam -c Fly/ChIP_seq/Input_WT.bam -f BAM -g dm -n Ez_WT -q 0.01 --outdir Peak
macs2 callpeak -t Fly/ChIP_seq/Pc_WT.bam -c Fly/ChIP_seq/Input_WT.bam -f BAM -g dm -n Pc_WT -q 0.01 --outdir Peak
macs2 callpeak -t Fly/ChIP_seq/Ph_WT.bam -c Fly/ChIP_seq/Input_WT.bam -f BAM -g dm -n Ph_WT -q 0.01 --outdir Peak
macs2 callpeak -t Fly/ChIP_seq/Pho_WT.bam -c Fly/ChIP_seq/Input_WT.bam -f BAM -g dm -n Pho_WT -q 0.01 --outdir Peak
macs2 callpeak -t Fly/ChIP_seq/Psc_WT.bam -c Fly/ChIP_seq/Input_WT.bam -f BAM -g dm -n Psc_WT -q 0.01 --outdir Peak
macs2 callpeak -t Fly/ChIP_seq/Spps_WT.bam -c Fly/ChIP_seq/Input_WT.bam -f BAM -g dm -n Spps_WT -q 0.01 --outdir Peak
```
#### Peak annotation
##### Annotation
```r
rm(list=ls())
library("ChIPseeker")
library("org.Dm.eg.db")
library("TxDb.Dmelanogaster.UCSC.dm6.ensGene")
txdb <- TxDb.Dmelanogaster.UCSC.dm6.ensGene
library("clusterProfiler")
anno_bed <- function(bedPeaksfile){peak<-readPeakFile(bedPeaksfile)
  #keepChr= !grepl("Het",seqlevels(peak))
  #table(keepChr)
  #seqlevels(peak, pruning.mode="coarse") <- seqlevels(peak)[keepChr]
  peakAnno <- annotatePeak(peak, tssRegion=c(-3000, 3000), 
                          TxDb=txdb, annoDb="org.Dm.eg.db") 
  peakAnno_df=as.data.frame(peakAnno)
  write.table(peakAnno_df,sep = "\t",row.names = F,file = paste0(basename(strsplit(bedPeaksfile, '\\.')[[1]][1]),".xls"))
  return(peakAnno_df)
}
tmp=lapply(list.files(path = 'Peak/',pattern = 'WT',full.names = T),anno_bed)
```
##### Peak distribution
```r
rm(list=ls())
library("ChIPseeker")
library("org.Dm.eg.db")
library("TxDb.Dmelanogaster.UCSC.dm6.ensGene")
txdb <- TxDb.Dmelanogaster.UCSC.dm6.ensGene
library("clusterProfiler")
anno_bed <- function(bedPeaksfile){peak<-readPeakFile(bedPeaksfile)
#keepChr= !grepl("Het",seqlevels(peak))
#table(keepChr)
#seqlevels(peak, pruning.mode="coarse") <- seqlevels(peak)[keepChr]
peakAnno <- annotatePeak(peak, tssRegion=c(-3000, 3000), 
                         TxDb=txdb, annoDb="org.Dm.eg.db") 
peakAnno_df=as.data.frame(peakAnno)
sampleName=basename(strsplit(bedPeaksfile, '\\.')[[1]][1])
return(peakAnno_df)
}
tmp=lapply(list.files(path = 'Peak/',pattern = 'WT',full.names = T),anno_bed)

df=lapply(tmp, function(x){
  #table(x$annotation)
  numb1=length(grep("Promoter", x$annotation))
  numb2=length(grep("5' UTR", x$annotation))
  numb3=length(grep("Exon", x$annotation))
  numb4=length(grep("Intron", x$annotation))
  numb5=length(grep("3' UTR", x$annotation))
  numb6=length(grep("Intergenic", x$annotation))
  return(c(numb1,numb2,numb3,numb4,numb5,numb6))
})
df=do.call(rbind,df)
colnames(df)=c("Promoter","5' UTR","Exon","Intron","3' UTR","Intergenic")
rownames(df)=unlist(lapply(list.files(path = 'Peak/',pattern = 'WT',full.names = T), function(x)
{sampleName=basename(strsplit(x, '\\.')[[1]][1])}))

library(ggpubr)
library(reshape2)
df2=melt(apply(df,1,function(x) x/sum(x)))
colnames(df2)=c("Location","Sample","Fraction")

ggbarplot(df2, "Sample", "Fraction",
          fill = "Location", color = "Location", palette = "lancet")
```
##### Venn
```r
rm(list=ls())
library("ChIPseeker")
library("org.Dm.eg.db")
library("TxDb.Dmelanogaster.UCSC.dm6.ensGene")
txdb <- TxDb.Dmelanogaster.UCSC.dm6.ensGene
library("clusterProfiler")
library("ChIPpeakAnno")
library("stringr")
BedFiles=list.files(path = 'Peak/',pattern = '*bed$',full.names = T)
fs=lapply(BedFiles,readPeakFile)
ol <- findOverlapsOfPeaks(fs[[1]],fs[[2]],fs[[3]])
makeVennDiagram(ol,
                NameOfPeaks = str_split(BedFiles,'_',simplify = T),
                TxDb = txdb)
```
##### Visualization
```r
rm(list=ls())
library("ChIPseeker")
library("org.Dm.eg.db")
library("TxDb.Dmelanogaster.UCSC.dm6.ensGene")
txdb <- TxDb.Dmelanogaster.UCSC.dm6.ensGene
library("clusterProfiler")
setwd("Peak")
Ez<-readPeakFile("Ez_WT_peaks.narrowPeak")
Pc<-readPeakFile("Pc_WT_peaks.narrowPeak")
Ph<-readPeakFile("Ph_WT_peaks.narrowPeak")
Pho<-readPeakFile("Pho_WT_peaks.narrowPeak")
Psc<-readPeakFile("Psc_WT_peaks.narrowPeak")
Spps<-readPeakFile("Spps_WT_peaks.narrowPeak")
peaks <- list(Ez=Ez,Pc=Pc,Ph=Ph,Pho=Pho,Psc=Psc,Spps=Spps)
promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
tagMatrixList <- lapply(peaks, getTagMatrix, windows=promoter)
plotAvgProf(tagMatrixList, xlim=c(-3000, 3000))
tagHeatmap(tagMatrixList, xlim=c(-3000, 3000), color=NULL)
```






















