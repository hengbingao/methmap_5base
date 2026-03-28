# MethMap

**MethMap** is a Python tool for methylation analysis using **5-base sequencing chemistry** (Illumina 5-Base Solution), in which **methylated cytosines (5mC) are selectively converted to thymine** while unmethylated cytosines remain unchanged. This is the **inverse** of traditional bisulfite sequencing (BS-seq).

MethMap supports both **directional** (Y-adaptor ligation) and **non-directional** (Tn5 transposase) library preparations, and includes built-in SNV discrimination based on paired-strand evidence.

---

## Table of Contents

- [Chemistry](#chemistry)
- [How It Works](#how-it-works)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Library Types](#library-types)
- [Commands](#commands)
- [SNV Discrimination](#snv-discrimination)
- [Tn5 End-Bias Correction](#tn5-end-bias-correction)
- [Output Files](#output-files)
- [Full Parameter Reference](#full-parameter-reference)
- [Batch Processing (SLURM)](#batch-processing-slurm)

---

## Chemistry

| Base | Traditional BS-seq | 5-base chemistry (MethMap) |
|------|--------------------|---------------------------|
| Methylated C (5mC) | C → C (protected) | **C → T (converted)** |
| Unmethylated C | C → T (converted) | **C → C (unchanged)** |

This chemistry enables direct detection of methylation without harsh bisulfite treatment, maintains high library complexity (only methylated C is converted), and enables simultaneous SNV calling.

---

## How It Works

MethMap aligns reads directly to CT- and GA-converted reference genomes:

```
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
```

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
- [bgzip](https://www.htslib.org/) *(for `.allc.gz` output; falls back to gzip if unavailable)*

### Install from source

```bash
git clone https://github.com/YOUR_USERNAME/methmap.git
cd methmap
pip install .

# Verify
methmap -h
```

### Install from release tarball

```bash
pip install methmap-1.3.5.tar.gz
```

---

## Quick Start

### Full pipeline (one command)

```bash
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
```

### Step by step

```bash
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
    --trim-r1       0 \
    --trim-r2       9 \
    --snv-filter \
    --cytosine-report \
    --allc
```

---

## Library Types

MethMap supports two library preparation strategies via `--library-type`:

### `directional` (default) — Y-adaptor ligation

Standard 5-base kit with Y-adaptors. R1 always reads the forward strand, R2 always reads the reverse strand. Only OT and OB strands are present.

```bash
methmap align --library-type directional ...
methmap dedup --library-type directional ...
```

Deduplication uses paired-end mode: `fixmate` + `markdup` (template-based).

### `non_directional` — Tn5 transposase

Tn5 cuts randomly and ligates adaptors to both ends, producing all four strands (OT/OB/CTOT/CTOB).

```bash
methmap align --library-type non_directional ...
methmap dedup --library-type non_directional ...
```

- Alignment: R1 and R2 are aligned independently as single-end reads, then merged
- Deduplication: `samtools markdup -m s` (sequenced-strand mode, no fixmate); falls back to positional dedup for samtools < 1.13

> **Note:** `methmap extract` does not require `--library-type` — methylation calling logic is identical for all four strands.

---

## Commands

### `methmap prepare`

Builds CT and GA converted genome indices. Run once per reference genome.

```bash
methmap prepare \
    --genome        /path/to/genome.fa \
    --output        /path/to/index_dir \
    --bowtie2-build /path/to/bowtie2-build
```

Index directory structure:

```
index_dir/
├── CT_conversion/
│   ├── genome_ct.fa
│   └── genome.*.bt2l
├── GA_conversion/
│   ├── genome_ga.fa
│   └── genome.*.bt2l
└── genome_info.txt
```

---

### `methmap align`

Aligns reads to both CT and GA genomes and selects the best alignment per read.

```bash
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
```

Custom SAM tags added by MethMap:

| Tag | Values | Meaning |
|-----|--------|---------|
| `XG` | `CT` / `GA` | Source genome used for alignment |
| `XR` | `OT` / `OB` / `CTOT` / `CTOB` | Strand identity |

---

### `methmap dedup`

Removes PCR duplicates using `samtools markdup`.

```bash
methmap dedup \
    --sam           ./results/MySample/MySample_aligned.sam \
    --output        ./results/MySample \
    --sample        MySample \
    --library-type  directional \
    --samtools      /path/to/samtools \
    --threads       8
```

Output files:
- `MySample_dedup.sam` — deduplicated alignments
- `MySample_dedup.bam` — sorted, indexed BAM
- `MySample_dedup_metrics.txt` — formatted duplication report

Example metrics:
```
==============================================================
  MethMap Deduplication Report
==============================================================
  Sample       : MySample
  Library type : directional
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
```

---

### `methmap extract`

Calls methylation from aligned SAM/BAM. Produces multiple output formats and always writes a methylation summary report.

```bash
methmap extract \
    --genome        /path/to/genome.fa \
    --sam           ./results/MySample/MySample_dedup.sam \
    --output        ./results/MySample \
    --sample        MySample \
    --min-coverage  5 \
    --trim-r1       0 \
    --trim-r2       9 \
    --snv-filter \
    --snv-ratio     0.2 \
    --snv-report \
    --bedgraph \
    --cytosine-report \
    --allc \
    --vcf \
    --threads       15
```

---

## SNV Discrimination

Per the Illumina 5-base solution (Figure 4), C→T SNVs are distinguished from methylation using forward-read signal at the adjacent G position:

```
                    + strand (pos=P)    Fwd read at pos=P+1   Interpretation
Unmethylated C      read=C             read=G                 C intact; G intact
5mC (methylated)    read=T             read=G                 5mC→T; G NOT mutated
C>T SNV             read=T             read=A                 C→T mutation; G→A on complement
```

**Key:** only forward reads (OT) at pos=P+1 are used for SNV calling. Reverse reads (OB) at pos=P+1 carry bottom-strand methylation signal and must not be confused with SNV signal.

Decision rule:
```
SNV if:  plus_T_frac > snv_ratio  AND  fwd_A_frac > snv_ratio
         (high T on + strand)          (G→A on complementary strand)
```

---

## Tn5 End-Bias Correction

Tn5 transposase introduces a 9 bp insertion bias at both ends of each fragment. The affected bases should be excluded from methylation calling.

**Which end to trim depends on your library construction:**

| Adapter end | Affected read | Parameter |
|-------------|--------------|-----------|
| P5 end bias | R1 (forward, FLAG 0x10=0) | `--trim-r1 N` |
| P7 end bias | R2 (reverse, FLAG 0x10=1) | `--trim-r2 N` |

Both can be set independently. The default is 0 (no trimming).

```bash
# Only P7 end has bias (most common for standard Tn5 CUT&Tag-style libraries)
methmap extract \
    --trim-r1 0 \
    --trim-r2 9 \
    ...

# Both ends have bias
methmap extract \
    --trim-r1 9 \
    --trim-r2 9 \
    ...
```

**How it works internally:**

```
R1 (forward read, FLAG 0x10=0):
  CIGAR walks low→high coordinate
  R1 5' end (P5) = start of CIGAR = aln_idx=0
  --trim-r1 N  →  skip positions where aln_idx < N

R2 (reverse read, FLAG 0x10=1):
  SAM stores RC of read; CIGAR still walks low→high
  R2 5' end (P7) = original read start = RC'd = end of CIGAR
  --trim-r2 N  →  skip positions where aln_idx >= (total_len - N)
```

---

## Output Files

| File | Flag | Description |
|------|------|-------------|
| `{sample}.bismark.cov.gz` | default | Bismark-compatible methylation coverage |
| `{sample}.bedGraph` | `--bedgraph` | UCSC browser BedGraph |
| `{sample}.CX_report.txt` | `--cytosine-report` | Full cytosine report (CpG/CHG/CHH) |
| `{sample}.SNV_report.txt` | `--snv-filter --snv-report` | Per-site SNV evidence table |
| `{sample}.snv.vcf` | `--snv-filter --vcf` | Standard VCF of SNV sites |
| `{sample}.allc.gz` | `--allc` | bgzip allc format (methylpy/ALLCools, tabix-indexable) |
| `{sample}.methylation_report.txt` | always | Bismark-style methylation summary report |
| `{sample}_dedup_metrics.txt` | after `dedup` | Duplication statistics |

### Methylation report format

Always generated after `methmap extract`. Example:

```
==================================================================
  MethMap Methylation Report
==================================================================
  Sample         : MySample
  Min coverage   : 1
  SNV exclusion  : disabled
==================================================================
  OVERALL (all cytosine contexts)
------------------------------------------------------------------
  Total cytosine sites covered   :       24,188
  Total methylated reads (M)     :        8,312
  Total unmethylated reads (U)   :       71,243
  Overall methylation level      :       10.45%
==================================================================
  PER-CONTEXT BREAKDOWN
------------------------------------------------------------------
  Context        Sites    Methylated  Unmethylated     Meth%    Coverage
  CpG            3,842         7,891        10,234    43.54%         4.7x
  CHG            9,823           218        30,112     0.72%         3.1x
  CHH           10,523           203        30,897     0.65%         3.0x
  --------------------------------------------------------------
  CpG + strand   methylation: 44.95%  (4,100 / 9,121)
  CpG - strand   methylation: 42.10%  (3,791 / 9,004)
==================================================================
  PER-CHROMOSOME CpG METHYLATION
------------------------------------------------------------------
  Chromosome                    Sites    Methylated  Unmethylated     Meth%
  lambda                        2,841           312         8,234     3.65%
  pUC19                         1,001         7,579         2,000    79.12%
==================================================================
```

### allc format (methylpy/ALLCools compatible)

Compressed with bgzip, tabix-indexable. Can be indexed after generation:

```bash
tabix -s 1 -b 2 -e 2 MySample.allc.gz
```

Format (tab-separated, no header):
```
chr1  10497  +  CGT  12  14  1
chr1  10498  -  CGC   8  10  1
```
`chromosome  pos(1-based)  strand  trinucleotide  mc_count  total  methylated(0/1)`

### Bismark coverage format

```
chr1    10497    10497    85.71    12    2
```
`chromosome  start(1-based)  end  methylation%  methylated_count  unmethylated_count`

### CX report format

```
chr1    10497    +    12    2    CpG    CGT
chr1    10525    +     3    2    CHG    CTG
```
`chromosome  pos(1-based)  strand  methylated  unmethylated  context  trinucleotide`

### SNV report format

```
# chr  pos  plus_T  plus_C  plus_T_frac  fwd_G  fwd_A  fwd_A_frac  rev_A  rev_G  verdict
chr1   100   190     10      0.950        185     5      0.026       3      180    methylation
chr1   200   185     15      0.925          2   178      0.989       2       5     SNV
```

---

## Full Parameter Reference

### `methmap prepare`

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--genome / -g` | required | Reference genome FASTA |
| `--output / -o` | required | Output index directory |
| `--bowtie2-build` | `bowtie2-build` | Path to bowtie2-build |

### `methmap align`

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--index / -x` | required | MethMap index directory |
| `--reads / -1` | required | Input FASTQ (SE or R1) |
| `--reads2 / -2` | — | R2 FASTQ (paired-end only) |
| `--output / -o` | required | Output directory |
| `--sample / -s` | required | Sample name |
| `--library-type` | `directional` | `directional` or `non_directional` (Tn5) |
| `--bowtie2` | `bowtie2` | Path to bowtie2 |
| `--threads / -p` | `1` | Number of threads |
| `--score-min` | `L,-0.6,-0.6` | Bowtie2 score threshold |
| `--no-discordant` | False | Discard discordant pairs |
| `--no-mixed` | False | Discard mixed alignments |

### `methmap dedup`

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--sam` | required | Aligned SAM file |
| `--output / -o` | — | Output directory |
| `--sample / -s` | — | Sample name |
| `--library-type` | `directional` | `directional` (PE dedup) or `non_directional` (SE dedup) |
| `--samtools` | `samtools` | Path to samtools |
| `--threads / -p` | `1` | Number of threads |

### `methmap extract`

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--genome / -g` | required | Reference genome FASTA |
| `--sam` | required | Aligned SAM or BAM file |
| `--output / -o` | required | Output directory |
| `--sample / -s` | required | Sample name |
| `--min-mapq` | `20` | Minimum mapping quality |
| `--min-base-qual` | `20` | Minimum base quality |
| `--min-coverage` | `1` | Minimum site coverage for output |
| `--context` | all | Contexts: `CpG` `CHG` `CHH` |
| `--trim-r1` | `0` | Skip first N bp of R1 (forward) reads — P5-end Tn5 bias correction |
| `--trim-r2` | `0` | Skip first N bp of R2 (reverse) reads — P7-end Tn5 bias correction |
| `--no-bismark-cov` | False | Skip Bismark coverage output |
| `--bedgraph` | False | Write BedGraph |
| `--cytosine-report` | False | Write full CX cytosine report |
| `--allc` | False | Write `.allc.gz` (bgzip, methylpy/ALLCools) |
| `--snv-filter` | False | Enable SNV discrimination |
| `--snv-ratio` | `0.2` | fwd_A_frac threshold for SNV calling |
| `--snv-report` | False | Write SNV evidence report |
| `--vcf` | False | Write VCF of SNV sites |
| `--threads / -p` | `1` | Accepted for pipeline compatibility |

### `methmap run`

Accepts all parameters from `prepare`, `align`, `dedup`, and `extract`, plus:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--remove-duplicates` | False | Run deduplication after alignment |

---

## Batch Processing (SLURM)

### Directional library (Y-adaptor, standard 5-base kit)

```bash
#!/bin/bash --login
#SBATCH --job-name=methmap
#SBATCH --partition=ll
#SBATCH --mem-per-cpu=20G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=15
#SBATCH --time=24:00:00
#SBATCH --export=NONE
#SBATCH --mail-user=your.email@university.edu
#SBATCH --mail-type=BEGIN,END

source activate
conda activate power

set -euo pipefail

WORKDIR=/path/to/workdir
GENOME=/path/to/genome.fa
INDEX=/path/to/methmap_index
BOWTIE2=/path/to/bowtie2
SAMTOOLS=/path/to/samtools
THREADS=15

cd ${WORKDIR}

for i in $(ls *_R1_001.fastq.gz | sed 's/_R1_001.fastq.gz//g' | cut -f1 | sort | uniq)
do
    fastp \
        -i ${i}_R1_001.fastq.gz \
        -I ${i}_R2_001.fastq.gz \
        -o ${i}.out.R1.fq.gz \
        -O ${i}.out.R2.fq.gz \
        --thread ${THREADS} \
        --html ${i}.fastp.html \
        --json ${i}.fastp.json

    methmap align \
        --index        ${INDEX} \
        --reads        ${i}.out.R1.fq.gz \
        --reads2       ${i}.out.R2.fq.gz \
        --output       ./${i} \
        --sample       ${i} \
        --library-type directional \
        --bowtie2      ${BOWTIE2} \
        --threads      ${THREADS}

    methmap dedup \
        --sam          ./${i}/${i}_aligned.sam \
        --output       ./${i} \
        --sample       ${i} \
        --library-type directional \
        --samtools     ${SAMTOOLS} \
        --threads      ${THREADS}

    methmap extract \
        --genome       ${GENOME} \
        --sam          ./${i}/${i}_dedup.sam \
        --output       ./${i} \
        --sample       ${i} \
        --trim-r1      0 \
        --trim-r2      9 \
        --cytosine-report \
        --allc \
        --snv-filter \
        --snv-report \
        --threads      ${THREADS}
done
```

### Non-directional library (Tn5-based CUT&Tag-style)

```bash
#!/bin/bash --login
#SBATCH --job-name=methmap_tn5
#SBATCH --cpus-per-task=15
#SBATCH --mem-per-cpu=20G
#SBATCH --time=24:00:00
#SBATCH --export=NONE

source activate
conda activate power

set -euo pipefail

WORKDIR=/path/to/workdir
GENOME=/path/to/genome.fa
INDEX=/path/to/methmap_index
BOWTIE2=/path/to/bowtie2
SAMTOOLS=/path/to/samtools
THREADS=15

cd ${WORKDIR}

for i in $(ls *_R1_001.fastq.gz | sed 's/_R1_001.fastq.gz//g' | cut -f1 | sort | uniq)
do
    fastp \
        -i ${i}_R1_001.fastq.gz \
        -I ${i}_R2_001.fastq.gz \
        -o ${i}.out.R1.fq.gz \
        -O ${i}.out.R2.fq.gz \
        --thread ${THREADS} \
        --html ${i}.fastp.html \
        --json ${i}.fastp.json

    methmap align \
        --index        ${INDEX} \
        --reads        ${i}.out.R1.fq.gz \
        --reads2       ${i}.out.R2.fq.gz \
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
        --trim-r1      9 \
        --trim-r2      9 \
        --cytosine-report \
        --allc \
        --snv-filter \
        --snv-report \
        --threads      ${THREADS}
done
```

---

## Comparison with Bismark

| Feature | Bismark (BS-seq) | MethMap (5-base) |
|---------|-----------------|------------------|
| Converted base | Unmethylated C → T | **Methylated C → T** |
| Unchanged base | Methylated C | **Unmethylated C** |
| read=T at C position | Unmethylated | **Methylated** |
| Library types | Directional | **Directional + Tn5/non-directional** |
| SNV discrimination | Not built-in | **Built-in (paired-strand, per Illumina Fig. 4)** |
| Tn5 end-bias correction | Not built-in | **Built-in (`--trim-r1` / `--trim-r2`)** |
| Output formats | Bismark | Bismark + BedGraph + CX + allc (bgzip) + VCF |
| Methylation report | Basic | **Full report (per-context, per-chromosome, strand-specific)** |
| Language | Perl | **Python ≥ 3.7** |
| Aligner | Bowtie2 / HISAT2 | Bowtie2 |
