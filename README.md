# MethMap

**MethMap** is a Python tool for methylation analysis using **5-base sequencing chemistry**, in which **methylated cytosines (5mC) are selectively converted to thymine** while unmethylated cytosines remain unchanged. This is the **inverse** of traditional bisulfite sequencing (BS-seq).

MethMap supports both **directional** (Y-adaptor) and **non-directional** (Tn5-based) library preparations.

---

## Table of Contents

- [Chemistry](#chemistry)
- [How It Works](#how-it-works)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Library Types](#library-types)
- [Commands](#commands)
- [SNV Discrimination](#snv-discrimination)
- [Output Files](#output-files)
- [Full Parameter Reference](#full-parameter-reference)
- [Batch Processing (SLURM)](#batch-processing-slurm)

---

## Chemistry

| Base | Traditional BS-seq | 5-base chemistry (MethMap) |
|------|--------------------|---------------------------|
| Methylated C (5mC) | C → C (protected) | **C → T (converted)** |
| Unmethylated C | C → T (converted) | **C → C (unchanged)** |

This chemistry enables direct detection of methylation without harsh bisulfite treatment and supports built-in **SNV discrimination** using paired-strand evidence.

---

## How It Works

MethMap aligns reads directly to CT- and GA-converted reference genomes:

\`\`\`
Reference genome
       │
       ├── C → T  →  CT genome  +  bowtie2 index
       └── G → A  →  GA genome  +  bowtie2 index

Input reads
       │
       ├── C → T  →  align to CT genome  (OT / CTOT strands)
       └── G → A  →  align to GA genome  (OB / CTOB strands)
                │
                └── Best alignment selected (by AS score)
                            │
                            ▼
                   Methylation extraction
                   ref=C, read=T  →  METHYLATED
                   ref=C, read=C  →  UNMETHYLATED
\`\`\`

### Strand assignment

| Genome | Orientation | Strand |
|--------|-------------|--------|
| CT | Forward | OT (Original Top) |
| CT | Reverse | CTOT (Complementary to OT) |
| GA | Forward | CTOB (Complementary to OB) |
| GA | Reverse | OB (Original Bottom) |

---

## Installation

### Requirements

- Python ≥ 3.7
- [Bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/) ≥ 2.3.0
- [samtools](https://www.htslib.org/) ≥ 1.10 *(for deduplication; ≥ 1.13 recommended for Tn5 libraries)*

### Install from source

\`\`\`bash
git clone https://github.com/YOUR_USERNAME/methmap.git
cd methmap
pip install .

# Verify
methmap -h
\`\`\`

### Install from release tarball

\`\`\`bash
pip install methmap-1.3.0.tar.gz
\`\`\`

---

## Quick Start

### Full pipeline (one command)

\`\`\`bash
methmap run \
    --genome        /path/to/genome.fa \
    --index         /path/to/methmap_index \
    --reads         sample_R1.fastq.gz \
    --reads2        sample_R2.fastq.gz \
    --output        ./results/MySample \
    --sample        MySample \
    --bowtie2       /path/to/bowtie2 \
    --bowtie2-build /path/to/bowtie2-build \
    --samtools      /path/to/samtools \
    --threads       8 \
    --remove-duplicates \
    --snv-filter \
    --cytosine-report \
    --allc
\`\`\`

### Step by step

\`\`\`bash
# 1. Build genome index (once per reference)
methmap prepare \
    --genome        /path/to/genome.fa \
    --output        /path/to/methmap_index \
    --bowtie2-build /path/to/bowtie2-build

# 2. Align reads
methmap align \
    --index         /path/to/methmap_index \
    --reads         sample_R1.fastq.gz \
    --reads2        sample_R2.fastq.gz \
    --output        ./results/MySample \
    --sample        MySample \
    --bowtie2       /path/to/bowtie2 \
    --threads       8

# 3. Remove duplicates
methmap dedup \
    --sam           ./results/MySample/MySample_aligned.sam \
    --output        ./results/MySample \
    --sample        MySample \
    --samtools      /path/to/samtools \
    --threads       8

# 4. Extract methylation
methmap extract \
    --genome        /path/to/genome.fa \
    --sam           ./results/MySample/MySample_dedup.sam \
    --output        ./results/MySample \
    --sample        MySample \
    --min-coverage  5 \
    --snv-filter \
    --cytosine-report \
    --allc
\`\`\`

---

## Library Types

MethMap supports two library preparation strategies via \`--library-type\`:

### \`directional\` (default) — Y-adaptor ligation

Standard 5-base kit using Y-adaptors. R1 always reads the forward strand, R2 always reads the reverse strand. Only OT and OB strands are present.

\`\`\`bash
methmap align --library-type directional ...
methmap dedup --library-type directional ...
\`\`\`

Deduplication uses paired-end mode: \`fixmate\` + \`markdup\` (template-based).

### \`non_directional\` — Tn5 transposase

Tn5 cuts randomly and ligates adaptors to both ends, producing a non-directional library with all four strands (OT/OB/CTOT/CTOB).

\`\`\`bash
methmap align --library-type non_directional ...
methmap dedup --library-type non_directional ...
\`\`\`

- Alignment: R1 and R2 aligned independently as single-end reads, then merged
- Deduplication: \`samtools markdup -m s\` (sequenced-strand mode, no fixmate); falls back to positional dedup for samtools < 1.13

> **Note:** \`methmap extract\` does not require \`--library-type\` — methylation calling logic is identical for all four strands.

---

## Commands

### \`methmap prepare\`

Builds CT and GA converted genome indices. Run once per reference genome.

\`\`\`bash
methmap prepare \
    --genome        /path/to/genome.fa \
    --output        /path/to/index_dir \
    --bowtie2-build /path/to/bowtie2-build
\`\`\`

Index directory structure:

\`\`\`
index_dir/
├── CT_conversion/
│   ├── genome_ct.fa
│   └── genome.*.bt2l
├── GA_conversion/
│   ├── genome_ga.fa
│   └── genome.*.bt2l
└── genome_info.txt
\`\`\`

---

### \`methmap align\`

Aligns reads to both CT and GA genomes and selects the best alignment per read.

\`\`\`bash
# Directional (default)
methmap align \
    --index         /path/to/index_dir \
    --reads         R1.fastq.gz \
    --reads2        R2.fastq.gz \
    --output        ./results/MySample \
    --sample        MySample \
    --library-type  directional \
    --bowtie2       /path/to/bowtie2 \
    --threads       8

# Non-directional (Tn5)
methmap align \
    --index         /path/to/index_dir \
    --reads         R1.fastq.gz \
    --reads2        R2.fastq.gz \
    --output        ./results/MySample \
    --sample        MySample \
    --library-type  non_directional \
    --bowtie2       /path/to/bowtie2 \
    --threads       8
\`\`\`

Custom SAM tags added by MethMap:

| Tag | Values | Meaning |
|-----|--------|---------|
| \`XG\` | \`CT\` / \`GA\` | Source genome used for alignment |
| \`XR\` | \`OT\` / \`OB\` / \`CTOT\` / \`CTOB\` | Strand identity |

---

### \`methmap dedup\`

Removes PCR duplicates using \`samtools markdup\`.

\`\`\`bash
methmap dedup \
    --sam           ./results/MySample/MySample_aligned.sam \
    --output        ./results/MySample \
    --sample        MySample \
    --library-type  non_directional \
    --samtools      /path/to/samtools \
    --threads       8
\`\`\`

Output files:
- \`MySample_dedup.sam\` — deduplicated alignments
- \`MySample_dedup.bam\` — sorted, indexed BAM
- \`MySample_dedup_metrics.txt\` — formatted duplication report

Example metrics:
\`\`\`
==============================================================
  MethMap Deduplication Report
==============================================================
  Sample       : MySample
  Library type : non_directional
  Tool         : samtools markdup
==============================================================
  Metric                                   Count  Fraction
  ----------------------------------------------------------
  Total reads (input)                  4,728,397
  Unique reads (retained)                 76,187    1.61%
  Duplicate reads (removed)            4,652,210   98.39%
    of which: duplicate single         4,652,210
    of which: duplicate pair                   0
  ----------------------------------------------------------
  Duplicate ratio                        0.9839
  Duplication rate                       98.39%
==============================================================
\`\`\`

---

### \`methmap extract\`

Calls methylation from aligned SAM/BAM with multiple output formats.

\`\`\`bash
methmap extract \
    --genome        /path/to/genome.fa \
    --sam           ./results/MySample/MySample_dedup.sam \
    --output        ./results/MySample \
    --sample        MySample \
    --min-coverage  5 \
    --snv-filter \
    --snv-ratio     0.2 \
    --snv-report \
    --bedgraph \
    --cytosine-report \
    --allc \
    --vcf \
    --threads       8
\`\`\`

---

## SNV Discrimination

At a CpG site (C on + strand at position P, G on − strand at position P+1):

\`\`\`
                    + strand (pos=P)    − strand (pos=P+1)
Unmethylated        read=C             read=G
Symmetric meth      read=T             read=A   (5mC→T, RC(T)=A)
Hemi-methylation    read=T             read=G   (only + strand methylated)
Homozygous SNV      read=T             read=T   (G→A mutation, RC(A)=T)
\`\`\`

The discriminating signal is **high \`minus_T_frac\`** at the paired G position:
- Methylation → minus strand A accumulates (counted normally)
- SNV → minus strand T accumulates (ignored at ref=G, so only T count rises)

Decision rule:
\`\`\`
SNV if:  plus_T_frac > snv_ratio  AND  minus_T_frac > snv_ratio
\`\`\`

---

## Output Files

| File | Flag | Description |
|------|------|-------------|
| \`{sample}.bismark.cov.gz\` | default | Bismark-compatible methylation coverage |
| \`{sample}.bedGraph\` | \`--bedgraph\` | UCSC browser BedGraph |
| \`{sample}.CX_report.txt\` | \`--cytosine-report\` | Full cytosine report (CpG/CHG/CHH) |
| \`{sample}.SNV_report.txt\` | \`--snv-filter --snv-report\` | Per-site SNV evidence table |
| \`{sample}.snv.vcf\` | \`--snv-filter --vcf\` | Standard VCF of SNV sites |
| \`{sample}.allc.tsv.gz\` | \`--allc\` | methylpy/ALLCools allc format (gzip) |
| \`{sample}_dedup_metrics.txt\` | after \`dedup\` | Duplication statistics |

---

## Full Parameter Reference

### \`methmap prepare\`

| Parameter | Default | Description |
|-----------|---------|-------------|
| \`--genome / -g\` | required | Reference genome FASTA |
| \`--output / -o\` | required | Output index directory |
| \`--bowtie2-build\` | \`bowtie2-build\` | Path to bowtie2-build |

### \`methmap align\`

| Parameter | Default | Description |
|-----------|---------|-------------|
| \`--index / -x\` | required | MethMap index directory |
| \`--reads / -1\` | required | Input FASTQ (SE or R1) |
| \`--reads2 / -2\` | — | R2 FASTQ (paired-end only) |
| \`--output / -o\` | required | Output directory |
| \`--sample / -s\` | required | Sample name |
| \`--library-type\` | \`directional\` | \`directional\` or \`non_directional\` (Tn5) |
| \`--bowtie2\` | \`bowtie2\` | Path to bowtie2 |
| \`--threads / -p\` | \`1\` | Number of threads |
| \`--score-min\` | \`L,-0.6,-0.6\` | Bowtie2 score threshold |
| \`--no-discordant\` | False | Discard discordant pairs |
| \`--no-mixed\` | False | Discard mixed alignments |

### \`methmap dedup\`

| Parameter | Default | Description |
|-----------|---------|-------------|
| \`--sam\` | required | Aligned SAM file |
| \`--output / -o\` | — | Output directory |
| \`--sample / -s\` | — | Sample name |
| \`--library-type\` | \`directional\` | \`directional\` (PE) or \`non_directional\` (SE) |
| \`--samtools\` | \`samtools\` | Path to samtools |
| \`--threads / -p\` | \`1\` | Number of threads |

### \`methmap extract\`

| Parameter | Default | Description |
|-----------|---------|-------------|
| \`--genome / -g\` | required | Reference genome FASTA |
| \`--sam\` | required | Aligned SAM or BAM file |
| \`--output / -o\` | required | Output directory |
| \`--sample / -s\` | required | Sample name |
| \`--min-mapq\` | \`20\` | Minimum mapping quality |
| \`--min-base-qual\` | \`20\` | Minimum base quality |
| \`--min-coverage\` | \`1\` | Minimum site coverage for output |
| \`--context\` | all | Contexts: \`CpG\` \`CHG\` \`CHH\` |
| \`--no-bismark-cov\` | False | Skip Bismark coverage output |
| \`--bedgraph\` | False | Write BedGraph |
| \`--cytosine-report\` | False | Write full CX cytosine report |
| \`--allc\` | False | Write \`.allc.tsv.gz\` (methylpy/ALLCools) |
| \`--snv-filter\` | False | Enable SNV discrimination |
| \`--snv-ratio\` | \`0.2\` | minus_T_frac threshold for SNV calling |
| \`--snv-report\` | False | Write SNV evidence report |
| \`--vcf\` | False | Write VCF of SNV sites |
| \`--threads / -p\` | \`1\` | Accepted for pipeline compatibility |

### \`methmap run\`

Accepts all parameters from \`prepare\`, \`align\`, \`dedup\`, and \`extract\`, plus:

| Parameter | Default | Description |
|-----------|---------|-------------|
| \`--remove-duplicates\` | False | Run deduplication after alignment |

---

## Batch Processing (SLURM)

### Non-directional library (Tn5)

\`\`\`bash
#!/bin/bash
#SBATCH --job-name=methmap_tn5
#SBATCH --cpus-per-task=20
#SBATCH --mem=64G
#SBATCH --time=24:00:00

set -euo pipefail

GENOME=/path/to/genome.fa
INDEX=/path/to/methmap_index
BOWTIE2=/path/to/bowtie2
SAMTOOLS=/path/to/samtools
THREADS=20

for R1 in *_R1_001.fastq.gz; do
    i=${R1/_R1_001.fastq.gz/}

    methmap align \
        --index        ${INDEX} \
        --reads        ${i}_R1_001.fastq.gz \
        --reads2       ${i}_R2_001.fastq.gz \
        --output       ./${i} \
        --sample       ${i} \
        --library-type non_directional \
        --bowtie2      ${BOWTIE2} \
        --threads      ${THREADS}

    methmap dedup \
        --sam          ./${i}/${i}_aligned.sam \
        --output       ./${i} \
        --sample       ${i} \
        --library-type non_directional \
        --samtools     ${SAMTOOLS} \
        --threads      ${THREADS}

    methmap extract \
        --genome       ${GENOME} \
        --sam          ./${i}/${i}_dedup.sam \
        --output       ./${i} \
        --sample       ${i} \
        --cytosine-report \
        --allc \
        --snv-filter \
        --snv-report \
        --threads      ${THREADS}
done
\`\`\`

---

## Comparison with Bismark

| Feature | Bismark (BS-seq) | MethMap (5-base) |
|---------|-----------------|------------------|
| Converted base | Unmethylated C → T | **Methylated C → T** |
| Unchanged base | Methylated C | **Unmethylated C** |
| read=T at C position | Unmethylated | **Methylated** |
| Library types | Directional | **Directional + Tn5/non-directional** |
| SNV discrimination | Not built-in | **Built-in (paired-strand)** |
| Output formats | Bismark | Bismark + BedGraph + CX + allc + VCF |
| Language | Perl | **Python ≥ 3.7** |
| Aligner | Bowtie2 / HISAT2 | Bowtie2 |
