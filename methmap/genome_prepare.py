"""
MethMap - Genome Preparation Module

For Illumina 5-base chemistry, we align reads DIRECTLY to the original
reference genome (no in-silico conversion needed).

This module builds a standard bowtie2 index from the original genome FASTA.
The CT/GA dual-genome strategy used in bisulfite sequencing is NOT needed here.
"""

import gzip
import logging
import subprocess
from pathlib import Path

logger = logging.getLogger("MethMap.GenomePrepare")


def read_fasta(fasta_path):
    """Read a FASTA file, return dict {chrom: sequence}."""
    sequences = {}
    current_chrom = None
    current_seq = []
    opener = gzip.open if str(fasta_path).endswith('.gz') else open
    with opener(fasta_path, 'rt') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_chrom is not None:
                    sequences[current_chrom] = ''.join(current_seq)
                current_chrom = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line.upper())
        if current_chrom is not None:
            sequences[current_chrom] = ''.join(current_seq)
    return sequences


def write_fasta(sequences, output_path, line_width=60):
    """Write dict {chrom: sequence} to FASTA file."""
    with open(output_path, 'w') as f:
        for chrom, seq in sequences.items():
            f.write('>{}\n'.format(chrom))
            for i in range(0, len(seq), line_width):
                f.write(seq[i:i + line_width] + '\n')


def reverse_complement(seq):
    comp = str.maketrans('ACGTNacgtn', 'TGCANtgcan')
    return seq.translate(comp)[::-1]


# Keep these for backward compatibility / testing
def convert_genome_CT(sequences):
    return {chrom: seq.replace('C', 'T') for chrom, seq in sequences.items()}


def convert_genome_GA(sequences):
    return {chrom: seq.replace('G', 'A') for chrom, seq in sequences.items()}


def _run_bowtie2_build(bowtie2_build_path, fasta, index_prefix):
    """Run bowtie2-build on the given FASTA."""
    cmd = [
        str(bowtie2_build_path),
        '--large-index',
        str(fasta),
        str(index_prefix),
    ]
    logger.info("Running: {}".format(' '.join(cmd)))
    ret = subprocess.call(cmd)
    if ret != 0:
        raise RuntimeError(
            "bowtie2-build failed (exit code {}) for: {}".format(ret, fasta))


def prepare_genome(genome_fasta, output_dir, bowtie2_build_path='bowtie2-build'):
    """
    Prepare bowtie2 index for MethMap 5-base alignment.

    For 5-base chemistry we build ONE standard index of the ORIGINAL genome.
    No CT/GA conversion is needed because reads are aligned directly.

    Output structure:
        output_dir/
            genome.*.bt2l          <- bowtie2 index of original genome
            genome_info.txt        <- chromosome sizes

    Args:
        genome_fasta:       Path to original reference FASTA (e.g. hg38.fa)
        output_dir:         Output directory for the index
        bowtie2_build_path: Path to bowtie2-build executable
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info("Reading genome from {}".format(genome_fasta))
    sequences = read_fasta(genome_fasta)
    total_size = sum(len(s) for s in sequences.values())
    logger.info("Loaded {} chromosomes, total {:,} bp".format(
        len(sequences), total_size))

    # Write genome info
    with open(str(output_dir / 'genome_info.txt'), 'w') as f:
        f.write("# MethMap genome index (5-base direct alignment)\n")
        f.write("# Source: {}\n".format(genome_fasta))
        f.write("# Chromosome\tSize\n")
        for chrom, seq in sequences.items():
            f.write("{}\t{}\n".format(chrom, len(seq)))

    # Build single bowtie2 index from original genome
    logger.info("Building bowtie2 index from original genome...")
    index_prefix = output_dir / 'genome'
    _run_bowtie2_build(bowtie2_build_path, genome_fasta, index_prefix)
    logger.info("Genome index built: {}".format(output_dir))

    return str(output_dir)
