"""
MethMap - Read Alignment Module

Illumina 5-base chemistry:
  5mC (methylated C)    -> T  (enzymatically converted)
  unmC (unmethylated C) -> C  (unchanged, ~95% of all C in human genome)

CORRECT ALIGNMENT STRATEGY (matching DRAGEN 5-base pipeline):
  Align ORIGINAL reads directly to the ORIGINAL (unconverted) reference genome.
  No in-silico conversion needed.

  After alignment, methylation is identified from mismatches:
    + strand:  ref=C, read=T  -> METHYLATED    (5mC was converted to T)
               ref=C, read=C  -> UNMETHYLATED  (C was preserved)
    - strand:  ref=G, read=A  -> METHYLATED    (complement of 5mC on - strand)
               ref=G, read=G  -> UNMETHYLATED  (complement of unmethylated C)

WHY THE PREVIOUS CT/GA DUAL-GENOME STRATEGY WAS WRONG:
  The old bisulfite strategy (e.g. Bismark) converts reads C->T before alignment
  because bisulfite converts UNMETHYLATED C->T, making almost all C disappear
  from reads. This requires converted reference genomes to allow alignment.

  In 5-base chemistry the opposite is true: METHYLATED C->T, and the vast
  majority of C (~95% in human) remain as C in the reads. The reads are
  almost identical to the original genome sequence and can be aligned directly.

  When we incorrectly applied CT/GA in-silico conversion to 5-base reads:
    - All remaining unmethylated C (95%) were converted to T in-silico
    - These T's were indistinguishable from true methylated T's
    - Everything looked methylated -> lambda showed ~70% instead of ~0%

STRAND HANDLING:
  5-base uses a directional library protocol (same as Bismark OT/OB).
  R1 reads the top strand (OT), R2 reads the bottom strand (OB).
  We tag each read with XG:Z:CT (+ strand evidence) or XG:Z:GA (- strand)
  based on read orientation, so extract.py can call methylation correctly.

  For paired-end:
    R1 forward (FLAG & 0x10 == 0): top strand    -> XG=CT, strand=+
    R1 reverse (FLAG & 0x10 != 0): bottom strand -> XG=GA, strand=-
    R2 forward (FLAG & 0x80 & ~0x10):            -> XG=GA, strand=-
    R2 reverse (FLAG & 0x80 &  0x10):            -> XG=CT, strand=+

  Simplified: strand = '+' if read is NOT reverse-complement, '-' otherwise.
  We record this as XG tag for compatibility with extract.py.
"""


# ------------------------------------------------------------------ #
#  Library type constants                                              #
# ------------------------------------------------------------------ #

DIRECTIONAL     = 'directional'      # OT + OB only  (default, Lister protocol)
NON_DIRECTIONAL = 'non_directional'  # OT + OB + CTOT + CTOB

import gzip
import logging
import subprocess
import shutil
from pathlib import Path

logger = logging.getLogger("MethMap.Align")


# ------------------------------------------------------------------ #
#  bowtie2 direct alignment                                            #
# ------------------------------------------------------------------ #

def run_bowtie2_direct(index, reads1, reads2=None, output_sam=None,
                       bowtie2_path='bowtie2', threads=1,
                       extra_args=None, report_file=None):
    """
    Align reads DIRECTLY to the original (unconverted) reference genome.

    Key differences from bisulfite alignment:
      - No --score-min relaxation needed (reads match reference well)
      - No --no-mixed / --no-discordant needed
      - Standard bowtie2 parameters work well

    Returns (return_code, stderr_string).
    """
    cmd = [
        bowtie2_path,
        '-x', str(index),
        '-p', str(threads),
        '--no-unal',
    ]

    if reads2:
        cmd.extend(['-1', str(reads1), '-2', str(reads2)])
    else:
        cmd.extend(['-U', str(reads1)])

    if output_sam:
        cmd.extend(['-S', str(output_sam)])

    if extra_args:
        cmd.extend(extra_args)

    logger.info("Running: {}".format(' '.join(cmd)))

    result = subprocess.run(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True,
    )

    if result.returncode != 0:
        logger.error("bowtie2 error:\n{}".format(result.stderr))
    else:
        logger.info("bowtie2 stats:\n{}".format(result.stderr))

    if report_file and result.stderr:
        with open(str(report_file), 'w') as f:
            f.write(result.stderr)

    return result.returncode, result.stderr


# ------------------------------------------------------------------ #
#  Strand tagging                                                      #
# ------------------------------------------------------------------ #

def tag_sam_strands(sam_in, sam_out, report_file=None):
    """
    Add XG and XR strand tags to each aligned read in the SAM file.

    For 5-base directional libraries:
      Forward read (FLAG & 0x10 == 0): reads + strand -> XG:Z:CT, XR:Z:OT
      Reverse read (FLAG & 0x10 != 0): reads - strand -> XG:Z:GA, XR:Z:OB

    These tags are used by extract.py to determine which base to examine
    (C on + strand, G on - strand) and what constitutes methylation (T or A).

    The XG tag convention matches what extract.py expects:
      XG=CT: read came from + strand; look at C in ref, T in read = methylated
      XG=GA: read came from - strand; look at G in ref, A in read = methylated
    """
    n_fwd = n_rev = 0

    with open(sam_in) as fin, open(sam_out, 'w') as fout:
        for line in fin:
            if line.startswith('@'):
                fout.write(line)
                continue

            parts = line.rstrip('\n').split('\t')
            if len(parts) < 11:
                fout.write(line)
                continue

            flag   = int(parts[1])
            is_rev = bool(flag & 0x10)

            if is_rev:
                # Read aligned to reverse strand: evidence for - strand methylation
                xg, xr = 'GA', 'OB'
                n_rev += 1
            else:
                # Read aligned to forward strand: evidence for + strand methylation
                xg, xr = 'CT', 'OT'
                n_fwd += 1

            parts.extend(['XG:Z:{}'.format(xg), 'XR:Z:{}'.format(xr)])
            fout.write('\t'.join(parts) + '\n')

    total = n_fwd + n_rev
    logger.info("Strand tagging: {:,} total reads "
                "({:,} forward/+strand, {:,} reverse/-strand)".format(
                    total, n_fwd, n_rev))

    if report_file:
        with open(str(report_file), 'a') as f:
            f.write("\n--- Strand Assignment ---\n")
            f.write("Forward (+ strand, XG=CT): {:,} ({:.1f}%)\n".format(
                n_fwd, 100.0 * n_fwd / total if total else 0))
            f.write("Reverse (- strand, XG=GA): {:,} ({:.1f}%)\n".format(
                n_rev, 100.0 * n_rev / total if total else 0))
            f.write("Total aligned: {:,}\n".format(total))

    return total


# ------------------------------------------------------------------ #
#  Single-end alignment pipeline                                       #
# ------------------------------------------------------------------ #

def align_single_end(reads_fastq, genome_index_dir, output_dir,
                     sample_name, bowtie2_path='bowtie2',
                     threads=1, extra_bowtie2_args=None):
    """
    Single-end 5-base alignment pipeline.

    Steps:
      1. Align original reads directly to original reference genome
      2. Tag each read with XG/XR strand information

    Returns path to final tagged SAM file.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    tmp_dir = output_dir / 'tmp_{}'.format(sample_name)
    tmp_dir.mkdir(exist_ok=True)

    genome_index_dir = Path(genome_index_dir)

    # For 5-base, we use a single index of the original genome
    # If user prepared with methmap prepare, use the original genome index
    # Check for direct index first, then fall back to CT_conversion (legacy)
    direct_index = genome_index_dir / 'genome'
    if not (genome_index_dir / 'genome.1.bt2').exists() and \
       not (genome_index_dir / 'genome.1.bt2l').exists():
        # Legacy: try CT_conversion path (will warn user)
        direct_index = genome_index_dir / 'CT_conversion' / 'genome'
        logger.warning(
            "No direct genome index found at {}. "
            "Using CT_conversion index as fallback. "
            "For correct 5-base analysis, rebuild the index with "
            "'methmap prepare' (it now builds a direct genome index).".format(
                genome_index_dir))

    logger.info("=== MethMap Single-End Alignment: {} ===".format(sample_name))
    logger.info("Strategy: direct alignment to original genome (5-base mode)")

    report = output_dir / '{}_alignment_report.txt'.format(sample_name)
    raw_sam = tmp_dir / '{}_raw.sam'.format(sample_name)

    # Step 1: direct alignment
    rc, _ = run_bowtie2_direct(
        direct_index, reads_fastq,
        output_sam=raw_sam,
        bowtie2_path=bowtie2_path,
        threads=threads,
        extra_args=extra_bowtie2_args,
        report_file=report,
    )
    if rc != 0:
        raise RuntimeError("Bowtie2 direct alignment failed")

    # Step 2: add strand tags
    final_sam = output_dir / '{}_aligned.sam'.format(sample_name)
    tag_sam_strands(str(raw_sam), str(final_sam), report_file=report)

    shutil.rmtree(str(tmp_dir))
    logger.info("Alignment report: {}".format(report))
    logger.info("Alignment complete: {}".format(final_sam))
    return str(final_sam)


# ------------------------------------------------------------------ #
#  Paired-end alignment pipeline                                       #
# ------------------------------------------------------------------ #

def align_paired_end(reads1_fastq, reads2_fastq, genome_index_dir,
                     output_dir, sample_name, bowtie2_path='bowtie2',
                     threads=1, extra_bowtie2_args=None,
                     library_type=None):
    """
    Paired-end 5-base alignment pipeline.

    library_type:
      'directional' (default):
        Standard Lister-protocol library. Only OT and OB strands exist.
        R1 always forward-aligns, R2 always reverse-aligns.
        Aligned as proper pairs: bowtie2 -1 R1 -2 R2

      'non_directional':
        All four strands present: OT, OB, CTOT, CTOB.
        R1 or R2 may align in either direction.
        Strategy: align R1 and R2 independently as single-end reads,
        then merge the two SAM files. This captures all four strand types
        without discarding discordant pairs.

    Steps (both modes):
      1. bowtie2 alignment
      2. Tag each aligned read with XG/XR strand info
    """
    import shutil as _shutil

    if library_type is None:
        library_type = DIRECTIONAL

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    tmp_dir = output_dir / 'tmp_{}'.format(sample_name)
    tmp_dir.mkdir(exist_ok=True)

    genome_index_dir = Path(genome_index_dir)
    direct_index = genome_index_dir / 'genome'
    if not (genome_index_dir / 'genome.1.bt2').exists() and \
       not (genome_index_dir / 'genome.1.bt2l').exists():
        direct_index = genome_index_dir / 'CT_conversion' / 'genome'
        logger.warning(
            "No direct genome index found at {}. "
            "Using CT_conversion index as fallback. "
            "Rebuild with 'methmap prepare' for correct 5-base analysis.".format(
                genome_index_dir))

    logger.info("=== MethMap Paired-End Alignment: {} ===".format(sample_name))
    logger.info("Library type: {}".format(library_type))
    logger.info("Strategy: direct alignment to original genome (5-base mode)")

    report  = output_dir / '{}_alignment_report.txt'.format(sample_name)
    raw_sam = tmp_dir    / '{}_raw.sam'.format(sample_name)

    if library_type == NON_DIRECTIONAL:
        # ── Non-directional: align R1 and R2 independently ──────────────
        # This captures all four strand types (OT, OB, CTOT, CTOB).
        # Each read is treated as single-end so bowtie2 does not enforce
        # concordant-pair orientation constraints.
        logger.info("Non-directional mode: aligning R1 and R2 independently")

        raw_r1 = tmp_dir / '{}_r1.sam'.format(sample_name)
        raw_r2 = tmp_dir / '{}_r2.sam'.format(sample_name)

        rc1, _ = run_bowtie2_direct(
            direct_index, reads1_fastq,
            output_sam=raw_r1,
            bowtie2_path=bowtie2_path,
            threads=threads,
            extra_args=extra_bowtie2_args,
            report_file=report,
        )
        if rc1 != 0:
            raise RuntimeError("Bowtie2 failed on R1")

        rc2, _ = run_bowtie2_direct(
            direct_index, reads2_fastq,
            output_sam=raw_r2,
            bowtie2_path=bowtie2_path,
            threads=threads,
            extra_args=extra_bowtie2_args,
            report_file=report,
        )
        if rc2 != 0:
            raise RuntimeError("Bowtie2 failed on R2")

        # Merge: keep header from R1, append all non-header lines from R2
        with open(str(raw_sam), 'w') as fout:
            with open(str(raw_r1)) as f1:
                for line in f1:
                    fout.write(line)
            with open(str(raw_r2)) as f2:
                for line in f2:
                    if not line.startswith('@'):
                        fout.write(line)

        logger.info("Merged R1 + R2 SAM files")

    else:
        # ── Directional: standard paired-end alignment ───────────────────
        logger.info("Directional mode: aligning as proper pairs")
        rc, _ = run_bowtie2_direct(
            direct_index, reads1_fastq, reads2_fastq,
            output_sam=raw_sam,
            bowtie2_path=bowtie2_path,
            threads=threads,
            extra_args=extra_bowtie2_args,
            report_file=report,
        )
        if rc != 0:
            raise RuntimeError("Bowtie2 direct PE alignment failed")

    # Tag strand information (same logic for both library types)
    final_sam = output_dir / '{}_aligned.sam'.format(sample_name)
    tag_sam_strands(str(raw_sam), str(final_sam), report_file=report)

    _shutil.rmtree(str(tmp_dir))
    logger.info("Alignment report: {}".format(report))
    logger.info("Alignment complete: {}".format(final_sam))
    return str(final_sam)


def remove_duplicates(input_sam, output_dir=None, sample_name=None,
                      samtools_path='samtools', threads=1,
                      library_type=None):
    """
    Remove PCR duplicates using samtools markdup.

    Two modes depending on library_type:

    directional (paired-end, default):
      Pipeline: name-sort -> fixmate -m -> coord-sort -> markdup -r
      Uses template-based dedup (-m t), requires proper paired-end reads.
      R1 and R2 are treated as a unit; duplicate = same fragment endpoints.

    non_directional (single-end, Tn5):
      Pipeline: coord-sort -> markdup -r -m s
      Skips fixmate (reads are single-end, no mate info available).
      Uses sequenced-strand mode (-m s): duplicate = same 5' position + strand.
      This is the correct strategy for Tn5-based libraries where R1 and R2
      were aligned independently.
    """
    if library_type is None:
        library_type = DIRECTIONAL

    input_sam = Path(input_sam)
    if output_dir is None:
        output_dir = input_sam.parent
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    if sample_name is None:
        sample_name = input_sam.stem

    tmp_dir     = output_dir / 'tmp_dedup_{}'.format(sample_name)
    tmp_dir.mkdir(exist_ok=True)
    sorted_bam  = tmp_dir    / '{}_sorted.bam'.format(sample_name)
    dedup_bam   = output_dir / '{}_dedup.bam'.format(sample_name)
    dedup_sam   = output_dir / '{}_dedup.sam'.format(sample_name)
    metrics_txt = output_dir / '{}_dedup_metrics.txt'.format(sample_name)

    def _run(cmd, step):
        logger.info("[dedup] {}: {}".format(step, ' '.join(cmd)))
        ret = subprocess.call(cmd)
        if ret != 0:
            raise RuntimeError(
                "samtools failed at step '{}' (exit code {})".format(step, ret))

    logger.info("[dedup] Library type: {}".format(library_type))

    if library_type == NON_DIRECTIONAL:
        # ── Single-end mode (Tn5 / non-directional) ─────────────────────
        # Reads were aligned independently (no mate info), so fixmate is
        # skipped. markdup -m s deduplicates by 5' position + strand.
        # -m s requires samtools >= 1.13; fall back to positional dedup if older.
        _run([samtools_path, 'sort', '-@', str(threads),
              '-o', str(sorted_bam), str(input_sam)], 'sort by coordinate')
        _run([samtools_path, 'index', str(sorted_bam)], 'index')

        # Check samtools version for -m s support (requires >= 1.13)
        import subprocess as _sp2
        _ver_out = _sp2.check_output([samtools_path, '--version'],
                                     stderr=_sp2.STDOUT,
                                     universal_newlines=True).splitlines()[0]
        # Parse "samtools 1.12" or "samtools 1.13+htslib-1.13"
        import re as _re
        _ver_match = _re.search(r'samtools (\d+)\.(\d+)', _ver_out)
        _use_m_s = False
        if _ver_match:
            _major, _minor = int(_ver_match.group(1)), int(_ver_match.group(2))
            _use_m_s = (_major > 1) or (_major == 1 and _minor >= 13)
        logger.info("[dedup] samtools version: {} | -m s support: {}".format(
            _ver_out.strip(), _use_m_s))

        if _use_m_s:
            markdup_cmd = [samtools_path, 'markdup', '-r', '-s',
                           '-m', 's',      # sequenced-strand mode (>= 1.13)
                           '-@', str(threads),
                           str(sorted_bam), str(dedup_bam)]
        else:
            # Fallback for samtools < 1.13:
            # markdup without -m s still works for SE reads using positional dedup.
            # Less accurate (no strand awareness) but functional.
            logger.warning("[dedup] samtools < 1.13: -m s not available, "
                           "using positional dedup (upgrade samtools for best results)")
            markdup_cmd = [samtools_path, 'markdup', '-r', '-s',
                           '-@', str(threads),
                           str(sorted_bam), str(dedup_bam)]
    else:
        # ── Paired-end mode (directional / standard) ─────────────────────
        # Full pipeline: name-sort -> fixmate (fills mate tags) ->
        # coord-sort -> markdup (template-based dedup)
        namesort_bam = tmp_dir / '{}_namesort.bam'.format(sample_name)
        fixmate_bam  = tmp_dir / '{}_fixmate.bam'.format(sample_name)

        _run([samtools_path, 'sort', '-n', '-@', str(threads),
              '-o', str(namesort_bam), str(input_sam)], 'sort by name')
        _run([samtools_path, 'fixmate', '-m', '-@', str(threads),
              str(namesort_bam), str(fixmate_bam)], 'fixmate')
        _run([samtools_path, 'sort', '-@', str(threads),
              '-o', str(sorted_bam), str(fixmate_bam)], 'sort by coordinate')
        _run([samtools_path, 'index', str(sorted_bam)], 'index')

        markdup_cmd = [samtools_path, 'markdup', '-r', '-s',
                       '-@', str(threads),
                       str(sorted_bam), str(dedup_bam)]

    logger.info("[dedup] markdup: {}".format(' '.join(markdup_cmd)))
    proc = subprocess.Popen(markdup_cmd, stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE, universal_newlines=True)
    _, stderr_out = proc.communicate()
    if proc.returncode != 0:
        raise RuntimeError("samtools failed at step 'markdup' (exit code {})".format(
            proc.returncode))
    # Always write a clean metrics file (parse from raw stats or zeros)
    # samtools markdup -s outputs stats to stderr in the format:
    #   READ 4262484 WRITTEN 75882 EXCLUDED 0 EXAMINED 4262484
    #   PAIRED 0 SINGLE 4262484
    #   DULPICATE PAIR 0 DUPLICATE SINGLE 4186602   (note: DULPICATE is a samtools typo)
    #   DUPLICATE TOTAL 4186602
    # All on one or multiple lines, space-separated key-value pairs.
    raw_stats = {}
    raw_text   = stderr_out.strip() if stderr_out.strip() else ""
    # Tokenize: split all whitespace, treat alternating tokens as key(s) + value
    # Strategy: scan for known integer-preceded keywords
    import re as _re2
    # Match patterns like "READ 4262484" or "DUPLICATE SINGLE 4186602"
    for m in _re2.finditer(
            r'(DULPICATE PAIR|DUPLICATE PAIR|DUPLICATE SINGLE|DUPLICATE TOTAL'
            r'|READ|WRITTEN|EXCLUDED|EXAMINED|PAIRED|SINGLE)'
            r'\s+(\d+)',
            raw_text):
        key = m.group(1).replace('DULPICATE', 'DUPLICATE')  # fix samtools typo
        raw_stats[key] = int(m.group(2))

    total      = raw_stats.get('READ',             0)
    written    = raw_stats.get('WRITTEN',           0)
    dup_single = raw_stats.get('DUPLICATE SINGLE',  0)
    dup_pair   = raw_stats.get('DUPLICATE PAIR',    0)
    dup_total  = raw_stats.get('DUPLICATE TOTAL',   dup_single + dup_pair)
    examined   = raw_stats.get('EXAMINED',          total)
    unique     = written
    dup_ratio  = dup_total / total if total > 0 else 0.0
    uniq_ratio = unique    / total if total > 0 else 0.0

    sep = "=" * 62
    report_lines = [
        sep,
        "  MethMap Deduplication Report",
        sep,
        "  Sample       : {}".format(sample_name),
        "  Library type : {}".format(library_type),
        "  Tool         : samtools markdup",
        sep,
        "  {:<32s}  {:>12s}  {:>8s}".format("Metric", "Count", "Fraction"),
        "  " + "-" * 58,
        "  {:<32s}  {:>12,}".format(
            "Total reads (input)",           total),
        "  {:<32s}  {:>12,}  {:>7.2%}".format(
            "Unique reads (retained)",       unique,     uniq_ratio),
        "  {:<32s}  {:>12,}  {:>7.2%}".format(
            "Duplicate reads (removed)",     dup_total,  dup_ratio),
        "  {:<32s}  {:>12,}".format(
            "  of which: duplicate single",  dup_single),
        "  {:<32s}  {:>12,}".format(
            "  of which: duplicate pair",    dup_pair),
        "  " + "-" * 58,
        "  {:<32s}  {:>11.4f}".format(
            "Duplicate ratio",               dup_ratio),
        "  {:<32s}  {:>11.2%}".format(
            "Duplication rate",              dup_ratio),
        sep,
    ]

    if raw_text:
        report_lines += [
            "  Raw samtools output:",
            "  " + raw_text.replace("\n", "\n  "),
            sep,
        ]

    with open(str(metrics_txt), 'w') as mf:
        mf.write("\n".join(report_lines) + "\n")

    # Log the summary
    logger.info("[dedup] " + "=" * 50)
    logger.info("[dedup] Sample:             {}".format(sample_name))
    logger.info("[dedup] Total reads:        {:,}".format(total))
    logger.info("[dedup] Unique (retained):  {:,}  ({:.2%})".format(unique, uniq_ratio))
    logger.info("[dedup] Duplicates removed: {:,}  ({:.2%})".format(dup_total, dup_ratio))
    logger.info("[dedup] Duplicate ratio:    {:.4f}  ({:.2%})".format(dup_ratio, dup_ratio))
    logger.info("[dedup] " + "=" * 50)

    _run([samtools_path, 'index', str(dedup_bam)], 'index dedup')
    _run([samtools_path, 'view', '-h', '-@', str(threads),
          '-o', str(dedup_sam), str(dedup_bam)], 'bam to sam')

    shutil.rmtree(str(tmp_dir))
    logger.info("Deduplication complete: {}".format(dedup_sam))
    return str(dedup_sam)


# ------------------------------------------------------------------ #
#  Legacy helpers (kept for test compatibility)                        #
# ------------------------------------------------------------------ #

def convert_reads_CT(input_fastq, output_fastq):
    """C->T conversion. Not used in 5-base pipeline, kept for compatibility."""
    import gzip as _gz
    opener_in  = _gz.open if str(input_fastq).endswith('.gz') else open
    opener_out = _gz.open if str(output_fastq).endswith('.gz') else open
    mode_out   = 'wt' if str(output_fastq).endswith('.gz') else 'w'
    count = 0
    with opener_in(input_fastq, 'rt') as fin, opener_out(output_fastq, mode_out) as fout:
        while True:
            header = fin.readline().strip()
            if not header:
                break
            seq  = fin.readline().strip()
            plus = fin.readline().strip()
            qual = fin.readline().strip()
            fout.write("{}\n{}\n{}\n{}\n".format(
                header, seq.replace('C', 'T'), plus, qual))
            count += 1
    return count


def convert_reads_GA(input_fastq, output_fastq):
    """G->A conversion. Not used in 5-base pipeline, kept for compatibility."""
    import gzip as _gz
    opener_in  = _gz.open if str(input_fastq).endswith('.gz') else open
    opener_out = _gz.open if str(output_fastq).endswith('.gz') else open
    mode_out   = 'wt' if str(output_fastq).endswith('.gz') else 'w'
    count = 0
    with opener_in(input_fastq, 'rt') as fin, opener_out(output_fastq, mode_out) as fout:
        while True:
            header = fin.readline().strip()
            if not header:
                break
            seq  = fin.readline().strip()
            plus = fin.readline().strip()
            qual = fin.readline().strip()
            fout.write("{}\n{}\n{}\n{}\n".format(
                header, seq.replace('G', 'A'), plus, qual))
            count += 1
    return count
