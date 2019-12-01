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

```

### Build Index
```bash
hisat2-build Fly/Drosophila.fa Fly/Hisat2_fly/fly
bowtie2-build Fly/Drosophila.fa Fly/Bowtie2_fly/fly
```


