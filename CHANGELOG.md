# Changelog




## [1.3.4] - 2026-03

### Changed
- **`write_allc` output**: file renamed from `.allc.tsv.gz` to `.allc.gz` (ALLCools convention)
- **BGZF compression**: allc output is now sorted by chr/pos and compressed with `bgzip` (BGZF format), enabling downstream `tabix` indexing
  - Fallback to standard gzip if `bgzip` is not available (warning logged)
  - New `--bgzip` parameter in `extract` and `run` subcommands (default: `bgzip`)
  - Sort is performed in-memory (equivalent to `sort -k1,1 -k2,2n | bgzip`)
## [1.3.3] - 2026-03

### Added
- **`write_methylation_report()`**: Bismark-style methylation summary report, always written automatically after `methmap extract`
  - Output file: `{sample}.methylation_report.txt`
  - Sections: overall methylation level, per-context breakdown (CpG/CHG/CHH with site count / methylated reads / unmethylated reads / meth% / average coverage), CpG + vs - strand comparison, per-chromosome CpG methylation, SNV filter summary (if enabled)
## [1.3.2] - 2026-03

### Fixed
- **SNV discrimination logic completely rewritten** based on Illumina 5-base Figure 4:
  - Previous logic confused OB-read signal (bottom-strand 5mC) with SNV signal at pos=P+1
  - Correct logic: use ONLY forward reads (OT) at pos=P+1 to judge whether the complementary G is intact (methylation) or mutated to A (C>T SNV)
  - `base_counts` replaced by `plus_base_counts` (forward reads) and `minus_base_counts` (reverse reads), stored separately to prevent signal mixing
  - `get_snv_flags`: now uses `plus_base_counts[pos+1]` A/(A+G) exclusively
  - `write_snv_report`: columns updated to `fwd_G / fwd_A / fwd_A_frac` (SNV signal) and `rev_A / rev_G` (bottom-strand methylation, for reference only)
- **Methylation extraction logic unchanged** — OB reads at ref=G: read=A → bottom-strand 5mC; read=G → bottom-strand umC (was always correct)
## [1.3.1] - 2026-03

### Fixed
- **Trinucleotide context off-by-one error** in `write_cytosine_report` and `write_vcf`:
  - Old (wrong): `ref_seq[pos-1:pos+2]` — started one base before the C, causing first character to be the upstream base instead of C itself (e.g. `TCC`, `GCC` instead of `CCT`, `CCC`)
  - New (correct): `ref_seq[pos:pos+3]` — correctly starts at the C position
- **Minus-strand trinucleotide missing reverse complement**: minus-strand entries now correctly report the bottom-strand trinucleotide as `RC(ref_seq[pos-2:pos+1])` instead of the raw forward-strand sequence

## [1.3.0] - 2026-03

### Added
- **Tn5 / non-directional library support** (`--library-type non_directional` in `align` and `dedup`)
  - R1 and R2 aligned independently as single-end reads to capture all four strands (OT/OB/CTOT/CTOB)
  - Deduplication uses `samtools markdup -m s` (sequenced-strand mode) instead of fixmate+markdup
  - Automatic fallback for samtools < 1.13 (positional dedup with warning)
- **`--vcf` output**: standard VCF file of SNV-flagged CpG sites, compatible with IGV/GATK/AnnoVar
- **`--allc` output**: gzip-compressed allc format (`.allc.tsv.gz`), compatible with methylpy and ALLCools
- **Deduplication metrics report** (`*_dedup_metrics.txt`): human-readable table with duplicate ratio, handles samtools `DULPICATE` typo in output parsing
- `--threads` accepted by `extract` subcommand (for pipeline script compatibility)

### Fixed
- **SNV discrimination logic corrected**: homozygous C→T SNV signature is high `minus_T_frac` (not low `minus_A_frac`). G→A mutation on minus strand → RC(A)=T, which extract.py ignores (not A or G), so T accumulates at ref=G positions
- `argparse` conflict: `--vcf` and `--allc` were registered twice in `extract` subparser
- `cmd_run` pipeline now passes `library_type` to dedup step
- `cmd_run` extract args namespace now includes `vcf` and `allc` fields

### Changed
- SNV report columns updated: `minus_A_frac` → `minus_T_frac` as the discriminating metric
- SNV report header updated with correct methylation vs SNV interpretation table

## [1.2.0] - 2026-02

### Added
- **VCF output** for SNV sites (`write_vcf`)
- **allc format output** for methylpy/ALLCools compatibility (`write_allc`)
- **Library type parameter** (`--library-type directional|non_directional`) in `align` subcommand

### Fixed
- Strand bug in methylation extraction for minus-strand reads

## [1.1.0] - 2026-01

### Added
- **SNV discrimination** at CpG sites using paired-strand evidence (`--snv-filter`, `--snv-ratio`, `--snv-report`)
- BedGraph output (`--bedgraph`)
- Full cytosine report (`--cytosine-report`)

## [1.0.0] - 2025-12

### Initial release
- 5-base chemistry alignment pipeline (5mC → T, unmethylated C → C)
- Bowtie2-based alignment to CT and GA converted genomes
- Bismark-compatible coverage output
- PCR duplicate removal via samtools markdup
- 15 unit tests
