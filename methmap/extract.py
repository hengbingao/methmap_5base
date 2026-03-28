"""
MethMap - Methylation Extraction Module

New chemistry calling logic:
  + strand:  ref=C, read=T -> METHYLATED   (5mC converted to T)
             ref=C, read=C -> UNMETHYLATED  (C unchanged)
  - strand:  ref=G, read=A -> METHYLATED   (complement of T)
             ref=G, read=G -> UNMETHYLATED  (complement of C)

SNV discrimination (using paired-strand evidence):
  A CpG site has both a C on + strand and a G on - strand.
  For a true SNV (C->T mutation), BOTH strands show the mutation:
    + strand: read=T  (looks like methylation)
    - strand: read=A  (G->A, the complement of T)
  For methylation:
    + strand: read=T  (5mC converted)
    - strand: read=G  (G unchanged, because the C on - strand is NOT methylated)
              OR read=T if symmetric methylation (both strands methylated)

  SNV filter logic (CpG only):
    At a CpG position, collect evidence from BOTH strands:
      + strand C position: T-count (looks meth), C-count (unmeth)
      - strand G position: G-count (unmeth complement), A-count (SNV signal)

    If A-count on - strand > SNV_strand_ratio * total - strand coverage:
      -> likely SNV, flag this site
    Else:
      -> methylation call is trusted

Context:
  CpG : C followed by G
  CHG : C followed by [ACT] followed by G
  CHH : C followed by [ACT] followed by [ACT]
"""

import os
import re
import logging
from pathlib import Path
from collections import defaultdict

logger = logging.getLogger("MethMap.Extract")


def get_complement(base):
    return {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}.get(base.upper(), 'N')


def get_context(ref_seq, pos, strand):
    """Return CpG / CHG / CHH / unknown for a cytosine at (pos, strand)."""
    n = len(ref_seq)

    if strand == '+':
        if pos >= n or ref_seq[pos] != 'C':
            return 'unknown'
        n1 = ref_seq[pos + 1] if pos + 1 < n else 'N'
        n2 = ref_seq[pos + 2] if pos + 2 < n else 'N'
        if n1 == 'G':
            return 'CpG'
        if n2 == 'G':
            return 'CHG'
        return 'CHH'

    else:   # minus strand: G in ref_seq corresponds to C on bottom strand
        if pos >= n or ref_seq[pos] != 'G':
            return 'unknown'
        p1 = ref_seq[pos - 1] if pos - 1 >= 0 else 'N'
        p2 = ref_seq[pos - 2] if pos - 2 >= 0 else 'N'
        c1 = get_complement(p1)
        c2 = get_complement(p2)
        if c1 == 'G':
            return 'CpG'
        if c2 == 'G':
            return 'CHG'
        return 'CHH'


def parse_cigar(cigar_str):
    """Parse CIGAR string -> list of (length, op)."""
    return [(int(l), op)
            for l, op in re.findall(r'(\d+)([MIDNSHP=X])', cigar_str)]


def cigar_to_alignment(cigar_str, read_seq, ref_start, ref_seq):
    """
    Walk CIGAR, yield (ref_pos, read_base, ref_base) for every aligned column.
    Insertions -> ref_pos = None; deletions -> read_base = '-'.
    """
    aligned  = []
    read_pos = 0
    ref_pos  = ref_start

    for length, op in parse_cigar(cigar_str):
        if op in ('M', '=', 'X'):
            for _ in range(length):
                rb = read_seq[read_pos] if read_pos < len(read_seq) else 'N'
                fb = ref_seq[ref_pos].upper() if ref_pos < len(ref_seq) else 'N'
                aligned.append((ref_pos, rb, fb))
                read_pos += 1
                ref_pos  += 1
        elif op == 'I':
            for _ in range(length):
                rb = read_seq[read_pos] if read_pos < len(read_seq) else 'N'
                aligned.append((None, rb, '-'))
                read_pos += 1
        elif op == 'D':
            for _ in range(length):
                fb = ref_seq[ref_pos].upper() if ref_pos < len(ref_seq) else 'N'
                aligned.append((ref_pos, '-', fb))
                ref_pos += 1
        elif op == 'N':
            ref_pos += length
        elif op == 'S':
            read_pos += length

    return aligned


class MethylationExtractor:
    """Extract methylation calls from a MethMap-aligned SAM file."""

    def __init__(self, genome_fasta, min_mapq=20, min_base_qual=20,
                 snv_filter=False, snv_ratio=0.2,
                 trim_r1=0, trim_r2=0):
        """
        Args:
            genome_fasta:   Reference genome FASTA
            min_mapq:       Minimum mapping quality (default 20)
            min_base_qual:  Minimum base quality (default 20)
            snv_filter:     Enable SNV discrimination at CpG sites (default False)
            snv_ratio:      If fraction of A-calls on minus strand > this threshold,
                            flag site as SNV (default 0.2, i.e. 20%)
            trim_r1:        Ignore the first N bp from the 5' end of R1 reads
                            (forward/+ strand reads, FLAG 0x10=0).
                            Use this to exclude Tn5 insertion bias on the P5 end.
                            Corresponds to skipping the first N aligned positions
                            in the CIGAR walk (low-coordinate end).
            trim_r2:        Ignore the first N bp from the 5' end of R2 reads
                            (reverse/- strand reads, FLAG 0x10=1).
                            Use this to exclude Tn5 insertion bias on the P7 end.
                            Because R2 is reverse-complemented in the SAM, the R2
                            5' end maps to the HIGH-coordinate end of the alignment;
                            trim_r2 skips the last N aligned positions in the CIGAR
                            walk.
        """
        self.min_mapq      = min_mapq
        self.min_base_qual = min_base_qual
        self.snv_filter    = snv_filter
        self.snv_ratio     = snv_ratio
        self.trim_r1       = int(trim_r1)
        self.trim_r2       = int(trim_r2)

        logger.info("Loading reference genome: {}".format(genome_fasta))
        from methmap.genome_prepare import read_fasta
        self.genome = read_fasta(genome_fasta)
        logger.info("Loaded {} chromosomes".format(len(self.genome)))

        # Primary methylation counts (strand-aware)
        # {chrom: {pos: {'M': int, 'U': int, 'context': str, 'strand': str}}}
        self.meth_calls = defaultdict(
            lambda: defaultdict(
                lambda: {'M': 0, 'U': 0, 'context': None, 'strand': None}
            )
        )

        # SNV evidence: per-strand base counts at each C/G position
        # Keyed by strand so OT and OB signals are never mixed.
        #
        # plus_base_counts[chrom][pos]: bases seen by FORWARD reads (OT/CTOB)
        #   at ref=C positions and the adjacent ref=G (pos+1)
        #   → pos+1 forward read signal distinguishes 5mC from C>T SNV:
        #       ref=G, read=G  → complementary G intact  → methylation
        #       ref=G, read=A  → complementary G mutated → C>T SNV
        #
        # minus_base_counts[chrom][pos]: bases seen by REVERSE reads (OB/CTOT)
        #   at ref=G positions (= bottom-strand C positions)
        #   → used for bottom-strand methylation calls only, NOT for SNV
        #
        # {chrom: {pos: {'A':0,'T':0,'C':0,'G':0}}}
        self.plus_base_counts  = defaultdict(
            lambda: defaultdict(lambda: {'A': 0, 'T': 0, 'C': 0, 'G': 0})
        )
        self.minus_base_counts = defaultdict(
            lambda: defaultdict(lambda: {'A': 0, 'T': 0, 'C': 0, 'G': 0})
        )

        # NonConverted sample base counts at CpG C positions
        # {chrom: {pos: {'T':0, 'C':0, 'other':0}}}
        # T in NonConverted = true SNV (no chemistry, so T must be a real mutation)
        # C in NonConverted = normal (unmethylated or methylated, but not SNV)
        self.nc_base_counts = defaultdict(
            lambda: defaultdict(
                lambda: {'T': 0, 'C': 0, 'other': 0}
            )
        )

    def process_sam(self, sam_path, context_filter=None, samtools_path='samtools'):
        """
        Parse aligned SAM or BAM and fill self.meth_calls (and self.base_counts if snv_filter).

        Args:
            sam_path:       Path to SAM or BAM. BAM detected by .bam extension.
            context_filter: list of contexts to keep e.g. ['CpG'], or None for all.
            samtools_path:  Path to samtools (used for BAM decoding).

        Returns:
            Number of reads processed
        """
        import subprocess as _sp

        sam_path = str(sam_path)
        is_bam = sam_path.endswith('.bam')

        if is_bam:
            logger.info("BAM input detected, reading via samtools view: {}".format(sam_path))
            proc = _sp.Popen(
                [samtools_path, 'view', '-h', sam_path],
                stdout=_sp.PIPE,
                universal_newlines=True,
            )
            line_iter = proc.stdout
        else:
            proc = None
            line_iter = open(sam_path)

        n_reads = n_calls = n_skipped = 0

        try:
            for line in line_iter:
                if line.startswith('@'):
                    continue

                parts = line.strip().split('\t')
                if len(parts) < 11:
                    continue

                flag = int(parts[1])
                if flag & 0x4:          # unmapped
                    continue

                chrom = parts[2]
                pos   = int(parts[3]) - 1   # 0-based
                mapq  = int(parts[4])
                cigar = parts[5]
                seq   = parts[9].upper()
                qual  = parts[10]

                if mapq < self.min_mapq:
                    n_skipped += 1
                    continue
                if chrom not in self.genome:
                    n_skipped += 1
                    continue
                if cigar == '*':
                    continue

                ref_seq = self.genome[chrom]

                # Determine strand directly from FLAG bit 0x10
                # 5-base direct alignment:
                #   FLAG 0x10 == 0  -> forward read -> + strand
                #                      ref=C, read=T = methylated
                #   FLAG 0x10 != 0  -> reverse read -> - strand
                #                      ref=G, read=A = methylated
                is_rev = bool(flag & 0x10)
                strand = '-' if is_rev else '+'

                try:
                    alignment = cigar_to_alignment(cigar, seq, pos, ref_seq)
                except Exception as e:
                    logger.debug("CIGAR error: {}".format(e))
                    n_skipped += 1
                    continue

                # ── Tn5 end-bias trimming ──────────────────────────────────
                # trim_r1: skip first N aligned positions for forward reads
                #          (R1, FLAG 0x10=0; 5' end = low-coordinate end)
                # trim_r2: skip last  N aligned positions for reverse reads
                #          (R2, FLAG 0x10=1; 5' end = high-coordinate end
                #           because SAM stores RC sequence, so 5' maps to end)
                aln_len = len(alignment) if (self.trim_r2 > 0 and is_rev) else 0

                for aln_idx, (ref_p, read_b, ref_b) in enumerate(alignment):
                    if ref_p is None or read_b == '-':
                        continue

                    # Trim bias region
                    if not is_rev and self.trim_r1 > 0:
                        if aln_idx < self.trim_r1:
                            continue
                    if is_rev and self.trim_r2 > 0:
                        if aln_idx >= aln_len - self.trim_r2:
                            continue

                    # Base quality filter
                    if qual != '*' and aln_idx < len(qual):
                        if ord(qual[aln_idx]) - 33 < self.min_base_qual:
                            continue

                    # New-chemistry methylation calling
                    meth_call = None

                    if strand == '+' and ref_b == 'C':
                        if read_b == 'T':
                            meth_call = 'M'     # methylated: 5mC -> T
                        elif read_b == 'C':
                            meth_call = 'U'     # unmethylated: C -> C

                    elif strand == '-' and ref_b == 'G':
                        if read_b == 'A':
                            meth_call = 'M'     # methylated (complement)
                        elif read_b == 'G':
                            meth_call = 'U'     # unmethylated (complement)

                    if meth_call is None:
                        continue

                    context = get_context(ref_seq, ref_p, strand)
                    if context_filter and context not in context_filter:
                        continue

                    # Record methylation call
                    call = self.meth_calls[chrom][ref_p]
                    call[meth_call] += 1
                    call['context']  = context
                    call['strand']   = strand
                    n_calls += 1

                    # Record per-strand base counts for SNV detection
                    # We collect ALL positions (C and G) on each strand separately
                    # so that get_snv_flags can isolate the forward-read signal
                    # at the G adjacent to each CpG C, which is the true SNV indicator.
                    if self.snv_filter:
                        if strand == '+':
                            bc = self.plus_base_counts[chrom][ref_p]
                        else:
                            bc = self.minus_base_counts[chrom][ref_p]
                        if read_b in bc:
                            bc[read_b] += 1

                n_reads += 1
                if n_reads % 100000 == 0:
                    logger.info("Processed {:,} reads, {:,} calls...".format(
                        n_reads, n_calls))

        finally:
            if proc is None:
                line_iter.close()
            else:
                proc.stdout.close()
                proc.wait()

        logger.info("Done: {:,} reads, {:,} calls, {:,} skipped".format(
            n_reads, n_calls, n_skipped))
        return n_reads

    def process_nonconverted_sam(self, sam_path, samtools_path='samtools'):
        """
        Parse a NonConverted (no chemistry) SAM file to collect base counts
        at CpG C positions on the + strand.

        In NonConverted samples, no chemical conversion is applied:
          - 5mC stays as C  (read=C, same as unmethylated C)
          - SNV (C->T)  stays as T  (read=T, the only way to see T here)

        Therefore:
          read=T at ref=C position in NonConverted -> TRUE SNV
          read=C at ref=C position in NonConverted -> normal (methylated or not)

        This is the gold standard for distinguishing:
          Symmetric methylation: 5base shows T+T, NonConverted shows C+C
          Homozygous SNV:        5base shows T+T, NonConverted shows T+T

        Results stored in self.nc_base_counts[chrom][pos].
        """
        import subprocess as _sp

        logger.info("Processing NonConverted SAM: {}".format(sam_path))
        sam_path = str(sam_path)
        is_bam   = sam_path.endswith('.bam')

        if is_bam:
            proc = _sp.Popen([samtools_path, 'view', '-h', sam_path],
                             stdout=_sp.PIPE, stderr=_sp.PIPE,
                             universal_newlines=True)
            line_iter = proc.stdout
        else:
            proc      = None
            line_iter = open(sam_path)

        n_reads = n_calls = 0

        try:
            for line in line_iter:
                if line.startswith('@'):
                    continue

                parts = line.rstrip('\n').split('\t')
                if len(parts) < 11:
                    continue

                flag = int(parts[1])
                if flag & 0x4:    # unmapped
                    continue
                if flag & 0x100:  # secondary
                    continue
                if flag & 0x800:  # supplementary
                    continue

                mapq = int(parts[4])
                if mapq < self.min_mapq:
                    continue

                chrom  = parts[2]
                pos    = int(parts[3]) - 1   # 0-based
                cigar  = parts[5]
                seq    = parts[9]
                qual   = parts[10]
                is_rev = bool(flag & 0x10)

                if chrom not in self.genome:
                    continue

                ref_seq = self.genome[chrom]

                try:
                    alignment = cigar_to_alignment(cigar, seq, pos, ref_seq)
                except Exception:
                    continue

                for aln_idx, (ref_p, read_b, ref_b) in enumerate(alignment):
                    if ref_p is None or read_b == '-':
                        continue

                    # Base quality filter
                    if qual != '*' and aln_idx < len(qual):
                        if ord(qual[aln_idx]) - 33 < self.min_base_qual:
                            continue

                    # Only look at + strand C positions (CpG context)
                    # In NonConverted: no chemistry, so T at ref=C = SNV
                    if is_rev:
                        continue   # only use forward reads for simplicity

                    if ref_b != 'C':
                        continue

                    # Check CpG context
                    if ref_p + 1 >= len(ref_seq) or ref_seq[ref_p + 1] != 'G':
                        continue

                    bc = self.nc_base_counts[chrom][ref_p]
                    if read_b == 'T':
                        bc['T'] += 1     # SNV signal
                    elif read_b == 'C':
                        bc['C'] += 1     # normal
                    else:
                        bc['other'] += 1

                    n_calls += 1

                n_reads += 1

        finally:
            if proc is None:
                line_iter.close()
            else:
                proc.stdout.close()
                proc.wait()

        logger.info("NonConverted: {:,} reads processed, {:,} CpG C positions recorded".format(
            n_reads, n_calls))
        return n_reads

    def get_snv_flags(self):
        """
        Identify homozygous C->T SNV sites at CpG positions.

        ═══════════════════════════════════════════════════════════
        Illumina 5-base chemistry — SNV vs methylation logic
        (per Illumina 5-base solution, Figure 4)

        At CpG pos=P (ref=C on + strand, ref=G on - strand at pos=P+1):

        The key discriminator is the base seen by FORWARD reads (OT) at
        pos=P+1 (the G adjacent to the CpG C):

          Situation          + strand pos=P    Fwd read at pos=P+1  Explanation
          ─────────────────  ──────────────    ──────────────────   ─────────────────
          Unmethylated C     read=C            read=G               C intact; G intact
          5mC (methylated)   read=T            read=G               5mC→T; G NOT mutated
          C>T SNV            read=T            read=A               C→T mutation; G→A on
                                                                    complementary strand

        CRITICAL: The G at pos=P+1 in a methylation event is NOT mutated —
        it stays G. Only in a true C>T SNV does the complementary strand
        carry an A at that position (because the DNA itself is mutated C→T,
        so the complement G→A as well).

        Note: OB reads at pos=P+1 carry bottom-strand methylation signal
        (5mC→T→RC=A), which looks like A and must NOT be confused with
        the SNV signal. We therefore use ONLY forward-read (plus_base_counts)
        signal at pos=P+1 to judge SNV.

        SNV decision rule:
          1. Forward reads at pos=P:   T/(T+C) > snv_ratio   (looks methylated)
          2. Forward reads at pos=P+1: A/(A+G) > snv_ratio   (G→A complement = SNV)

        Both conditions must hold simultaneously.
        ═══════════════════════════════════════════════════════════

        Returns:
            set of (chrom, pos) tuples that are likely homozygous C->T SNVs
        """
        snv_sites = set()

        for chrom, chrom_data in self.meth_calls.items():
            ref_seq = self.genome.get(chrom, '')

            for pos, call in chrom_data.items():
                if call.get('context') != 'CpG':
                    continue
                if call.get('strand') != '+':
                    continue

                # ── Condition 1: forward reads at pos=P look like T ──────────
                plus_T = call['M']   # OT reads: ref=C, read=T
                plus_C = call['U']   # OT reads: ref=C, read=C
                plus_total = plus_T + plus_C
                if plus_total == 0:
                    continue
                plus_t_frac = plus_T / plus_total
                if plus_t_frac <= self.snv_ratio:
                    continue   # not even T-dominant, skip

                # ── Condition 2: forward reads at pos=P+1 show A not G ───────
                # Use ONLY plus_base_counts (OT reads) at pos=P+1.
                # OT reads align forward; at ref=G (pos=P+1):
                #   read=G → complementary G intact → methylation
                #   read=A → complementary G mutated to A → C>T SNV ← this is the signal
                # DO NOT use minus_base_counts here: OB reads at pos=P+1
                # show A for bottom-strand 5mC, not for SNV.
                paired_pos = pos + 1
                if paired_pos >= len(ref_seq) or ref_seq[paired_pos] != 'G':
                    continue

                fwd_at_g = self.plus_base_counts.get(chrom, {}).get(paired_pos, None)
                if fwd_at_g is None:
                    continue

                fwd_A = fwd_at_g['A']   # SNV signal: complementary G→A
                fwd_G = fwd_at_g['G']   # methylation: G intact
                fwd_AG_total = fwd_A + fwd_G
                if fwd_AG_total == 0:
                    continue

                # A/(A+G) from forward reads at the paired G position
                fwd_a_frac = fwd_A / fwd_AG_total

                if fwd_a_frac > self.snv_ratio:
                    snv_sites.add((chrom, pos))
                    logger.debug(
                        "SNV flagged: {}:{} plus_T_frac={:.3f} fwd_A_frac={:.3f} "
                        "(plus_T={} fwd_A={} fwd_G={})".format(
                            chrom, pos + 1, plus_t_frac, fwd_a_frac,
                            plus_T, fwd_A, fwd_G))

        logger.info("SNV filter: {:,} CpG sites flagged as likely homozygous C->T SNVs".format(
            len(snv_sites)))
        return snv_sites

    # ------------------------------------------------------------------ #
    #  Output writers                                                      #
    # ------------------------------------------------------------------ #

    def write_bismark_coverage(self, output_path, min_coverage=1, exclude_snv=False):
        """
        Bismark-style coverage file.
        Format: chr  start(1-based)  end  methylation%  count_meth  count_unmeth

        Args:
            exclude_snv: If True and snv_filter was enabled, skip SNV-flagged sites
        """
        snv_sites = self.get_snv_flags() if (exclude_snv and self.snv_filter) else set()
        n = n_snv = 0

        with open(output_path, 'w') as f:
            for chrom in sorted(self.meth_calls):
                for pos in sorted(self.meth_calls[chrom]):
                    if (chrom, pos) in snv_sites:
                        n_snv += 1
                        continue
                    c = self.meth_calls[chrom][pos]
                    m, u = c['M'], c['U']
                    if m + u < min_coverage:
                        continue
                    pct = round(100.0 * m / (m + u), 2)
                    f.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(
                        chrom, pos + 1, pos + 1, pct, m, u))
                    n += 1

        if exclude_snv and self.snv_filter:
            logger.info("SNV filter excluded {:,} sites".format(n_snv))
        logger.info("Written {:,} sites to {}".format(n, output_path))
        return n

    def write_bedgraph(self, output_path, min_coverage=1, exclude_snv=False):
        """BedGraph methylation file."""
        snv_sites = self.get_snv_flags() if (exclude_snv and self.snv_filter) else set()
        n = 0

        with open(output_path, 'w') as f:
            f.write("track type=bedGraph\n")
            for chrom in sorted(self.meth_calls):
                for pos in sorted(self.meth_calls[chrom]):
                    if (chrom, pos) in snv_sites:
                        continue
                    c = self.meth_calls[chrom][pos]
                    m, u = c['M'], c['U']
                    if m + u < min_coverage:
                        continue
                    f.write("{}\t{}\t{}\t{}\n".format(
                        chrom, pos, pos + 1,
                        round(100.0 * m / (m + u), 4)))
                    n += 1

        logger.info("Written {:,} sites to BedGraph {}".format(n, output_path))
        return n

    def write_cytosine_report(self, output_path, min_coverage=0, exclude_snv=False):
        """
        Full cytosine report.
        Format: chr  pos(1-based)  strand  count_meth  count_unmeth  context  trinucleotide
        """
        snv_sites = self.get_snv_flags() if (exclude_snv and self.snv_filter) else set()
        n = 0

        with open(output_path, 'w') as f:
            f.write("# MethMap Cytosine Report\n")
            f.write("# chr\tpos(1-based)\tstrand\t"
                    "count_methylated\tcount_unmethylated\tcontext\ttrinucleotide\n")
            for chrom in sorted(self.meth_calls):
                ref_seq = self.genome.get(chrom, '')
                for pos in sorted(self.meth_calls[chrom]):
                    if (chrom, pos) in snv_sites:
                        continue
                    c = self.meth_calls[chrom][pos]
                    m, u = c['M'], c['U']
                    if m + u < min_coverage:
                        continue
                    ctx    = c.get('context', 'unknown')
                    strand = c.get('strand',  '+')
                    # trinucleotide context centred on the cytosine
                    # + strand: C is at pos, next two bases are pos+1, pos+2
                    # - strand: C is the RC of G at pos; the two upstream ref
                    #           bases (pos-1, pos-2) become the downstream context
                    #           after RC, so trinuc = RC(ref[pos-2:pos+1])
                    if strand == '+':
                        if pos + 2 < len(ref_seq):
                            trinuc = ref_seq[pos:pos + 3].upper()
                        else:
                            trinuc = ref_seq[pos:].upper().ljust(3, 'N')
                    else:  # minus strand
                        if pos - 2 >= 0:
                            fwd = ref_seq[pos - 2:pos + 1].upper()
                        else:
                            fwd = ('N' * (2 - pos) + ref_seq[0:pos + 1]).upper()
                        # reverse complement to get the bottom-strand trinuc
                        _comp = str.maketrans('ACGTN', 'TGCAN')
                        trinuc = fwd.translate(_comp)[::-1]
                    f.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                        chrom, pos + 1, strand, m, u, ctx, trinuc))
                    n += 1

        logger.info("Written {:,} sites to cytosine report {}".format(n, output_path))
        return n

    def write_snv_report(self, output_path, min_coverage=1):
        """
        Write SNV evidence report for all CpG sites.

        Columns (per Illumina 5-base Figure 4 logic):
          plus_T / plus_C    forward-read T and C counts at pos=P (the CpG C)
          plus_T_frac        T/(T+C) — high = C is T (either 5mC or SNV)
          fwd_G / fwd_A      forward-read G and A counts at pos=P+1 (the adjacent G)
          fwd_A_frac         A/(A+G) from forward reads at pos=P+1:
                               low  → G intact → 5mC methylation
                               high → G→A mutation → C>T SNV
          rev_A / rev_G      reverse-read A and G at pos=P+1 (bottom-strand meth signal)
          verdict            methylation / SNV
        """
        if not self.snv_filter:
            logger.warning("SNV filter was not enabled. No SNV report written.")
            return 0

        snv_sites = self.get_snv_flags()
        n = 0

        with open(output_path, 'w') as f:
            f.write("# MethMap SNV Report (Illumina 5-base chemistry)\n")
            f.write("# Logic (per Illumina Figure 4):\n")
            f.write("#   5mC:  plus_T_frac HIGH + fwd_A_frac LOW  "
                    "(G intact on complementary strand)\n")
            f.write("#   SNV:  plus_T_frac HIGH + fwd_A_frac HIGH "
                    "(G->A mutation on complementary strand)\n")
            f.write("# IMPORTANT: fwd_A_frac uses ONLY forward reads at pos+1.\n")
            f.write("# rev_A at pos+1 reflects bottom-strand 5mC, not SNV.\n")
            f.write("# snv_ratio threshold: {}\n".format(self.snv_ratio))
            f.write("# chr\tpos(1-based)\t"
                    "plus_T\tplus_C\tplus_T_frac\t"
                    "fwd_G\tfwd_A\tfwd_A_frac\t"
                    "rev_A\trev_G\t"
                    "verdict\n")

            for chrom in sorted(self.meth_calls):
                for pos in sorted(self.meth_calls[chrom]):
                    call = self.meth_calls[chrom][pos]
                    if call.get('context') != 'CpG':
                        continue
                    if call.get('strand') != '+':
                        continue

                    plus_T = call['M']
                    plus_C = call['U']
                    plus_total = plus_T + plus_C
                    if plus_total < min_coverage:
                        continue

                    plus_t_frac = plus_T / plus_total
                    paired_pos  = pos + 1

                    # Forward reads at pos+1: SNV indicator
                    fwd = self.plus_base_counts.get(chrom, {}).get(paired_pos, {})
                    fwd_G = fwd.get('G', 0)
                    fwd_A = fwd.get('A', 0)
                    fwd_ag = fwd_A + fwd_G
                    fwd_a_frac = fwd_A / fwd_ag if fwd_ag > 0 else 0.0

                    # Reverse reads at pos+1: bottom-strand methylation signal
                    rev = self.minus_base_counts.get(chrom, {}).get(paired_pos, {})
                    rev_A = rev.get('A', 0)   # bottom-strand 5mC -> T -> RC=A
                    rev_G = rev.get('G', 0)   # bottom-strand umC -> C -> RC=G

                    verdict = 'SNV' if (chrom, pos) in snv_sites else 'methylation'

                    f.write("{}\t{}\t{}\t{}\t{:.3f}\t{}\t{}\t{:.3f}\t{}\t{}\t{}\n".format(
                        chrom, pos + 1,
                        plus_T, plus_C, plus_t_frac,
                        fwd_G, fwd_A, fwd_a_frac,
                        rev_A, rev_G,
                        verdict))
                    n += 1

        logger.info("Written {:,} CpG sites to SNV report: {}".format(n, output_path))
        return n

    def print_summary(self):
        """Print per-context methylation statistics."""
        stats = defaultdict(lambda: {'M': 0, 'U': 0, 'sites': 0})
        for chrom_data in self.meth_calls.values():
            for call in chrom_data.values():
                ctx = call.get('context', 'unknown')
                stats[ctx]['M']     += call['M']
                stats[ctx]['U']     += call['U']
                stats[ctx]['sites'] += 1

        print("\n=== MethMap Methylation Summary ===")
        print("{:<10} {:>10} {:>8} {:>12} {:>14}".format(
            "Context", "Sites", "Meth%", "Methylated", "Unmethylated"))
        print("-" * 58)
        for ctx in ['CpG', 'CHG', 'CHH']:
            s = stats.get(ctx, {'M': 0, 'U': 0, 'sites': 0})
            tot = s['M'] + s['U']
            pct = "{:.1f}%".format(100.0 * s['M'] / tot) if tot else "N/A"
            print("{:<10} {:>10,} {:>8} {:>12,} {:>14,}".format(
                ctx, s['sites'], pct, s['M'], s['U']))
        print("=" * 58)

        if self.snv_filter:
            snv_count = len(self.get_snv_flags())
            print("SNV-flagged CpG sites: {:,}".format(snv_count))
            print("=" * 58)


    def write_methylation_report(self, output_path, sample_name='sample',
                                  min_coverage=1, exclude_snv=False):
        """
        Write a Bismark-style methylation summary report.

        Sections:
          1. Overall statistics (all cytosines)
          2. Per-context breakdown: CpG / CHG / CHH
          3. Per-chromosome CpG methylation
          4. SNV summary (if snv_filter enabled)
        """
        snv_sites = self.get_snv_flags() if (exclude_snv and self.snv_filter) else set()

        # ── Accumulate stats ──────────────────────────────────────────────────
        # global + per-context counts: M (methylated reads), U (unmethylated reads), sites
        total   = {'M': 0, 'U': 0, 'sites': 0}
        ctx_stats = {ctx: {'M': 0, 'U': 0, 'sites': 0}
                     for ctx in ['CpG', 'CHG', 'CHH']}

        # per-chrom CpG stats
        chrom_cpg = {}

        for chrom in sorted(self.meth_calls):
            cc = {'M': 0, 'U': 0, 'sites': 0}
            for pos, call in self.meth_calls[chrom].items():
                if (chrom, pos) in snv_sites:
                    continue
                m   = call['M']
                u   = call['U']
                cov = m + u
                if cov < min_coverage:
                    continue
                ctx = call.get('context', 'unknown')

                total['M']     += m
                total['U']     += u
                total['sites'] += 1

                if ctx in ctx_stats:
                    ctx_stats[ctx]['M']     += m
                    ctx_stats[ctx]['U']     += u
                    ctx_stats[ctx]['sites'] += 1

                if ctx == 'CpG':
                    cc['M']     += m
                    cc['U']     += u
                    cc['sites'] += 1

            if cc['sites'] > 0:
                chrom_cpg[chrom] = cc

        def _pct(m, u):
            tot = m + u
            return 100.0 * m / tot if tot > 0 else 0.0

        def _fmt(m, u):
            return "{:.2f}%  ({:,} / {:,})".format(_pct(m, u), m, m + u)

        # ── Build report text ─────────────────────────────────────────────────
        sep  = "=" * 66
        sep2 = "-" * 66
        lines = []

        lines += [
            sep,
            "  MethMap Methylation Report",
            sep,
            "  Sample         : {}".format(sample_name),
            "  Min coverage   : {}".format(min_coverage),
            "  SNV exclusion  : {}".format("enabled" if exclude_snv and self.snv_filter
                                            else "disabled"),
            sep,
        ]

        # ── Section 1: Overall ───────────────────────────────────────────────
        lines += [
            "  OVERALL (all cytosine contexts)",
            sep2,
            "  Total cytosine sites covered   : {:>12,}".format(total['sites']),
            "  Total methylated reads (M)     : {:>12,}".format(total['M']),
            "  Total unmethylated reads (U)   : {:>12,}".format(total['U']),
            "  Overall methylation level      : {:>11.2f}%".format(
                _pct(total['M'], total['U'])),
            sep,
        ]

        # ── Section 2: Per-context ───────────────────────────────────────────
        lines += [
            "  PER-CONTEXT BREAKDOWN",
            sep2,
            "  {:<8}  {:>10}  {:>12}  {:>12}  {:>8}  {:>10}".format(
                "Context", "Sites", "Methylated", "Unmethylated",
                "Meth%", "Coverage"),
            "  " + "-" * 62,
        ]

        for ctx in ['CpG', 'CHG', 'CHH']:
            s = ctx_stats[ctx]
            tot_reads = s['M'] + s['U']
            avg_cov   = tot_reads / s['sites'] if s['sites'] > 0 else 0.0
            lines.append(
                "  {:<8}  {:>10,}  {:>12,}  {:>12,}  {:>7.2f}%  {:>10.1f}x".format(
                    ctx, s['sites'], s['M'], s['U'],
                    _pct(s['M'], s['U']), avg_cov))

        # Also report CpG separately as both strands combined
        cpg_plus  = {'M': 0, 'U': 0}
        cpg_minus = {'M': 0, 'U': 0}
        for chrom_data in self.meth_calls.values():
            for pos, call in chrom_data.items():
                if call.get('context') != 'CpG':
                    continue
                m, u = call['M'], call['U']
                if m + u < min_coverage:
                    continue
                if call.get('strand') == '+':
                    cpg_plus['M'] += m; cpg_plus['U'] += u
                else:
                    cpg_minus['M'] += m; cpg_minus['U'] += u

        lines += [
            "  " + "-" * 62,
            "  CpG + strand   methylation: {}".format(
                _fmt(cpg_plus['M'],  cpg_plus['U'])),
            "  CpG - strand   methylation: {}".format(
                _fmt(cpg_minus['M'], cpg_minus['U'])),
            sep,
        ]

        # ── Section 3: Per-chromosome CpG ───────────────────────────────────
        if chrom_cpg:
            lines += [
                "  PER-CHROMOSOME CpG METHYLATION",
                sep2,
                "  {:<25}  {:>8}  {:>12}  {:>12}  {:>8}".format(
                    "Chromosome", "Sites", "Methylated", "Unmethylated", "Meth%"),
                "  " + "-" * 62,
            ]
            for chrom, cc in sorted(chrom_cpg.items()):
                lines.append(
                    "  {:<25}  {:>8,}  {:>12,}  {:>12,}  {:>7.2f}%".format(
                        chrom, cc['sites'], cc['M'], cc['U'],
                        _pct(cc['M'], cc['U'])))
            lines.append(sep)

        # ── Section 4: SNV summary ───────────────────────────────────────────
        if self.snv_filter:
            n_snv = len(self.get_snv_flags()) if not snv_sites else len(snv_sites)
            cpg_total_sites = ctx_stats['CpG']['sites']
            lines += [
                "  SNV FILTER SUMMARY",
                sep2,
                "  CpG sites examined             : {:>12,}".format(cpg_total_sites),
                "  SNV-flagged sites (excluded)   : {:>12,}".format(n_snv),
                "  SNV ratio threshold            : {:>11.2f}".format(self.snv_ratio),
                sep,
            ]

        with open(str(output_path), 'w') as f:
            f.write('\n'.join(lines) + '\n')

        logger.info("Written methylation report: {}".format(output_path))
        return output_path

    def write_vcf(self, output_path, sample_name="SAMPLE", min_coverage=10):
        """
        Write SNV-flagged CpG sites as a standard VCF file.

        Only sites flagged by get_snv_flags() are written (C->T mutation at CpG).
        Genotype inferred from T/(T+C) allele fraction on the plus strand:
          >= 0.8  -> homozygous  1/1
          0.2-0.8 -> heterozygous 0/1
          < 0.2   -> skip (not a real SNV)

        INFO:   DP, VAF, CONTEXT (trinucleotide), MINUS_A, MINUS_DP
        FORMAT: GT:DP:AD:VAF
        """
        if not self.snv_filter:
            logger.warning("SNV filter was not enabled. No VCF written.")
            return 0

        import datetime
        snv_sites = self.get_snv_flags()
        n = 0

        with open(output_path, 'w') as f:
            # VCF header
            for line in [
                "##fileformat=VCFv4.2",
                "##fileDate=" + datetime.date.today().strftime("%Y%m%d"),
                "##source=MethMap",
                "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total read depth\">",
                "##INFO=<ID=VAF,Number=1,Type=Float,Description=\"Variant allele fraction\">",
                "##INFO=<ID=CONTEXT,Number=1,Type=String,Description=\"Trinucleotide context\">",
                "##INFO=<ID=MINUS_A,Number=1,Type=Integer,Description=\"A-count on minus strand at paired G\">",
                "##INFO=<ID=MINUS_DP,Number=1,Type=Integer,Description=\"Total depth on minus strand at paired G\">",
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
                "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read depth\">",
                "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths REF,ALT\">",
                "##FORMAT=<ID=VAF,Number=1,Type=Float,Description=\"Variant allele fraction\">",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + sample_name,
            ]:
                f.write(line + "\n")

            for chrom in sorted(self.meth_calls):
                ref_seq = self.genome.get(chrom, '')
                for pos in sorted(self.meth_calls[chrom]):
                    if (chrom, pos) not in snv_sites:
                        continue
                    call = self.meth_calls[chrom][pos]
                    if call.get('context') != 'CpG':
                        continue
                    if call.get('strand') != '+':
                        continue

                    t_count = call['M']   # read=T (could be SNV or methylation)
                    c_count = call['U']   # read=C (ref allele)
                    depth   = t_count + c_count
                    if depth < min_coverage:
                        continue

                    vaf = t_count / depth
                    if vaf >= 0.8:
                        gt = "1/1"
                    elif vaf >= 0.2:
                        gt = "0/1"
                    else:
                        continue

                    # Minus strand evidence (paired G at pos+1)
                    paired_pos = pos + 1
                    bc       = self.base_counts.get(chrom, {}).get(paired_pos, {})
                    minus_a  = bc.get('A', 0)
                    minus_dp = bc.get('A', 0) + bc.get('G', 0)

                    if pos + 2 < len(ref_seq):
                        trinuc = ref_seq[pos:pos + 3].upper()
                    else:
                        trinuc = ref_seq[pos:].upper().ljust(3, 'N')

                    info = "DP={};VAF={:.4f};CONTEXT={};MINUS_A={};MINUS_DP={}".format(
                        depth, vaf, trinuc, minus_a, minus_dp)
                    samp = "{}:{}:{},{}: {:.4f}".format(
                        gt, depth, c_count, t_count, vaf)

                    f.write("{}\t{}\t.\tC\tT\t.\tPASS\t{}\tGT:DP:AD:VAF\t{}\n".format(
                        chrom, pos + 1, info, samp))
                    n += 1

        logger.info("Written {:,} SNV sites to VCF: {}".format(n, output_path))
        return n

    def write_allc(self, output_path, min_coverage=1, exclude_snv=False):
        """
        Write methylation in allc format (methylpy / ALLCools compatible).

        Format (tab-separated, no header):
          chr  pos(1-based)  strand  trinucleotide  mc_count  total  methylated(0/1)

          trinucleotide: 3-base context always starting with C
            + strand: ref_seq[pos:pos+3]            e.g. CGT, CAT, CCG
            - strand: RC(ref_seq[pos-2:pos+1])      e.g. CAT, CTG
          methylated: 1 if mc_count/total >= 0.9, else 0

        Writes plain text (sorted by chr/pos). Caller is responsible for
        bgzip compression (done by _bgzip() in __main__.py).
        Compatible with: methylpy, ALLCools, allcools_utilities
        """
        _comp = str.maketrans('ACGTN', 'TGCAN')

        def _trinuc(ref_seq, pos, strand):
            if strand == '+':
                raw = ref_seq[pos:pos + 3] if pos + 2 < len(ref_seq) \
                    else ref_seq[pos:].ljust(3, 'N')
                return raw.upper()
            else:
                if pos - 2 >= 0:
                    fwd = ref_seq[pos - 2:pos + 1].upper()
                else:
                    fwd = ('N' * (2 - pos) + ref_seq[0:pos + 1]).upper()
                return fwd.translate(_comp)[::-1]

        snv_sites = self.get_snv_flags() if (exclude_snv and self.snv_filter) else set()
        n = 0

        # Collect and sort rows before writing (required for tabix indexing)
        rows = []
        for chrom in self.meth_calls:
            ref_seq = self.genome.get(chrom, '')
            for pos in self.meth_calls[chrom]:
                if (chrom, pos) in snv_sites:
                    continue
                c = self.meth_calls[chrom][pos]
                m, u = c['M'], c['U']
                if m + u < min_coverage:
                    continue
                strand   = c.get('strand', '+')
                trinuc   = _trinuc(ref_seq, pos, strand)
                meth_bin = 1 if (m / (m + u) >= 0.9) else 0
                rows.append((chrom, pos,
                             "{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                                 chrom, pos + 1, strand,
                                 trinuc, m, m + u, meth_bin)))
                n += 1

        rows.sort(key=lambda r: (r[0], r[1]))

        with open(str(output_path), 'w') as f:
            for _, _, line in rows:
                f.write(line)

        logger.info("Written {:,} sites to allc: {}".format(n, output_path))
        return n
