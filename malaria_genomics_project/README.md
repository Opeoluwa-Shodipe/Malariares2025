# 🧬 Comparative Genomic Analysis of *Plasmodium falciparum* Field Isolates to Identify Genetic Variants Associated with Antimalarial Drug Resistance

---

## 🚀 Project Overview

### 🎯 Main Aim
This project aims to identify and analyze genomic variants—such as SNPs, indels, and CNVs—in *Plasmodium falciparum* field isolates from various regions and associate them with known or novel antimalarial drug resistance markers.

### 🌍 Why It Matters
- Antimalarial drug resistance is a major global health concern.
- Genomic surveillance supports evidence-based policy decisions for treatment.
- This pipeline builds skills in WGS analysis, variant calling, and interpretation in a real-world malaria research context.

---

## 📌 Workflow Summary

### 🛠 Tools & Pipelines
| Process                | Tools Used                                       |
|------------------------|--------------------------------------------------|
| Quality Control        | FastQC, Trimmomatic                              |
| Read Mapping           | BWA, Bowtie2                                     |
| Variant Calling        | GATK, FreeBayes                                  |
| Annotation             | SnpEff, VEP                                      |
| Population Analysis    | Plink, bcftools, R (adegenet, ggplot2)           |
| Visualization          | IGV, ggplot2, matplotlib, seaborn                |

---

## 📚 Step-by-Step Roadmap

### ✅ PHASE 1: Project Design
- Define research questions on resistance-associated variants.
- Select datasets from NCBI SRA https://www.ncbi.nlm.nih.gov/sra/ (Nigeria).
- Choose reference genome: *P. falciparum* 3D7 (PlasmoDB- https://plasmodb.org/).

---

### ✅ PHASE 2: Data Acquisition
```bash
# Download reads from SRA
prefetch SRR1234567
fastq-dump --split-files SRR1234567

# Download reference genome (FASTA + GFF)
wget https://plasmodb.org/common/downloads/release-64/Pf3D7/fasta/data/PlasmoDB-64_Pf3D7_Genome.fasta
wget https://plasmodb.org/common/downloads/release-64/Pf3D7/gff/data/PlasmoDB-64_Pf3D7.gff.
---

```
### ✅ PHASE 3: Read Mapping & Variant Calling
```...
# Run FastQC
fastqc *.fastq

# Trimming
trimmomatic PE sample_R1.fastq sample_R2.fastq \
sample_R1_paired.fq sample_R1_unpaired.fq \
sample_R2_paired.fq sample_R2_unpaired.fq \
ILLUMINACLIP:adapters.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:50
---

```
### ✅ PHASE 4: Read Mapping & Variant Calling
```...
# Index reference genome
bwa index Pf3D7.fasta

# Align reads
bwa mem Pf3D7.fasta sample_R1_paired.fq sample_R2_paired.fq > sample.sam

# Convert to BAM, sort, and index
samtools view -bS sample.sam | samtools sort -o sample.sorted.bam
samtools index sample.sorted.bam

# Variant calling
gatk HaplotypeCaller -R Pf3D7.fasta -I sample.sorted.bam -O sample.g.vcf -ERC GVCF

# Joint genotyping
gatk GenotypeGVCFs -R Pf3D7.fasta -V gendb://your_database -O cohort.vcf


