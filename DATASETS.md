# Datasets

## Overview

We selected ten FASTQ datasets to ensure a comprehensive evaluation across diverse genomic features, sequencing technologies, and file scales. The collection spans multiple research domains including human genomics, environmental metagenomics, transgenic monitoring, and microbial genomics, with file sizes ranging from several hundred megabytes to over 15 GB.

## Dataset Selection Rationale

Our dataset selection primarily includes a subset of the benchmark test set proposed by Numanagić and Bonfield et al., specifically the high-repeat, high-GC *S. cerevisiae* dataset (SRR554369), which is recognized by the ISO/IEC 23092 (MPEG-G) committee for compression evaluation. To validate performance on larger scales and diverse libraries, we incorporated:

- **Human genomic data** (SRR1210085\_1) sequenced on Illumina HiSeq 2000, representing a large-scale whole-genome sequencing scenario (~15 GB).
- **Environmental metagenomes** (DRR057887\_1), presenting mixed-species complexity that challenges compression algorithms lacking a single reference genome.
- **Transgenic monitoring datasets** (SRR14139158 and SRR14626645), each with paired-end reads (\_1 and \_2), covering plant genomic and environmental monitoring contexts on Illumina HiSeq 4000.
- **Complementary library strategies for *B. thuringiensis* HS18-1**, featuring both a short-insert paired-end library (SRR2093871\_1) and a large-insert mate-pair library (SRR2093872\_1), enabling evaluation across different library preparation methods.
- **Microbial genomics** (SRR29245815), a *B. subtilis* GL-4 isolate from bamboo rat gut, adding diversity in organism type and sequencing depth.

This selection enables a robust evaluation of algorithm scalability and efficiency across heterogeneous sequencing data.

## Dataset Details

| Dataset | Organism / Project Context | Technology | Size (bytes) | SRA Link |
|---|---|---|---|---|
| DRR057887\_1 | Metagenome (Environmental sample) | DDBJ/Illumina | 570,579,485 | [DRR057887](https://www.ncbi.nlm.nih.gov/sra/?term=DRR057887) |
| SRR1210085\_1 | *Homo sapiens* (Genomic DNA) | Illumina HiSeq 2000 | 15,126,475,212 | [SRR1210085](https://www.ncbi.nlm.nih.gov/sra/?term=SRR1210085) |
| SRR14139158\_1 | Transgenic monitoring (*S. lycopersicum*) | Illumina HiSeq 4000 | 10,859,255,381 | [SRR14139158](https://www.ncbi.nlm.nih.gov/sra/?term=SRR14139158) |
| SRR14139158\_2 | Transgenic monitoring (*S. lycopersicum*) | Illumina HiSeq 4000 | 10,859,255,381 | [SRR14139158](https://www.ncbi.nlm.nih.gov/sra/?term=SRR14139158) |
| SRR14626645\_1 | Environmental monitoring | Illumina HiSeq 4000 | 8,845,364,612 | [SRR14626645](https://www.ncbi.nlm.nih.gov/sra/?term=SRR14626645) |
| SRR14626645\_2 | Environmental monitoring | Illumina HiSeq 4000 | 8,845,364,612 | [SRR14626645](https://www.ncbi.nlm.nih.gov/sra/?term=SRR14626645) |
| SRR2093871\_1 | *B. thuringiensis* HS18-1 (PE library) | Illumina HiSeq 2000 | 1,817,852,165 | [SRR2093871](https://www.ncbi.nlm.nih.gov/sra/?term=SRR2093871) |
| SRR2093872\_1 | *B. thuringiensis* HS18-1 (Mate-Pair) | Illumina HiSeq 2000 | 1,907,047,029 | [SRR2093872](https://www.ncbi.nlm.nih.gov/sra/?term=SRR2093872) |
| SRR29245815 | *B. subtilis* GL-4 (Bamboo rat gut) | Illumina HiSeq 4000 | 3,649,758,240 | [SRR29245815](https://www.ncbi.nlm.nih.gov/sra/?term=SRR29245815) |
| SRR554369 | *S. cerevisiae* (MPEG-G Benchmark) | Illumina GA II | 788,017,922 | [SRR554369](https://www.ncbi.nlm.nih.gov/sra/?term=SRR554369) |

## Download

All datasets are publicly available from the [NCBI Sequence Read Archive (SRA)](https://www.ncbi.nlm.nih.gov/sra/). Paired-end reads (\_1 / \_2) originate from the same accession. You can download them using the [SRA Toolkit](https://github.com/ncbi/sra-tools):

```bash
# Example: download SRR554369
prefetch SRR554369
fastq-dump --split-files SRR554369
```
