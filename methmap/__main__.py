#!/usr/bin/env python3
"""
MethMap - Methylation Alignment Tool for Illumina 5-Base Sequencing

Chemistry: Methylated C (5mC) -> T  |  Unmethylated C -> C (unchanged)

Alignment strategy: Direct alignment to original genome (no in-silico conversion).
Methylation identified from C->T mismatches in aligned reads.

Commands:
  prepare   Build genome bowtie2 index (run once per reference)
  align     Align reads directly to genome
  dedup     Remove PCR duplicates
  extract   Extract methylation calls
  run       Full pipeline in one step

Examples:
  methmap prepare --genome hg38.fa --output /ref/methmap_index

  methmap align   --index /ref/methmap_index --reads R1.fastq.gz --reads2 R2.fastq.gz \\
                  --output ./out --sample MySample --threads 8

  methmap dedup   --sam ./out/MySample_aligned.sam --output ./out --sample MySample

  methmap extract --genome hg38.fa --sam ./out/MySample_dedup.sam \\
                  --output ./out --sample MySample \\
                  --bedgraph --cytosine-report --snv-filter --snv-report

  methmap run     --genome hg38.fa --index /ref/methmap_index \\
                  --reads R1.fastq.gz --reads2 R2.fastq.gz \\
                  --output ./out --sample MySample --threads 8 \\
                  --remove-duplicates --snv-filter --snv-report
"""

import argparse
import logging
import os
import subprocess
import sys
import gzip
import shutil
from pathlib import Path


def _bgzip(src_path, samtools_path='samtools'):
    """
    Compress a plain-text file with bgzip, producing a tabix-compatible .gz.
    Tries: bgzip → samtools bgzip → Python gzip (fallback).
    The source file is replaced by the .gz version.
    Returns the path to the compressed file.
    """
    src     = str(src_path)
    gz_path = src + '.gz'

    for cmd in ['bgzip', ['samtools', 'bgzip']]:
        args = ([cmd, '-f', src] if isinstance(cmd, str)
                else cmd + ['-f', src])
        try:
            ret = subprocess.call(args,
                stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            if ret == 0:
                return gz_path
        except FileNotFoundError:
            continue

    # Fallback to Python gzip (not tabix-compatible but functional)
    with open(src, 'rb') as fi, gzip.open(gz_path, 'wb') as fo:
        shutil.copyfileobj(fi, fo)
    os.remove(src)
    logger.warning("bgzip not found — used Python gzip (output not tabix-indexable): {}".format(gz_path))
    return gz_path

logging.basicConfig(
    level=logging.INFO,
    format='[%(asctime)s] %(name)s %(levelname)s: %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S',
)
logger = logging.getLogger("MethMap")


# ------------------------------------------------------------------ #
#  Sub-commands                                                        #
# ------------------------------------------------------------------ #

def cmd_prepare(args):
    from methmap.genome_prepare import prepare_genome
    logger.info("=== MethMap Genome Preparation ===")
    prepare_genome(
        genome_fasta=args.genome,
        output_dir=args.output,
        bowtie2_build_path=args.bowtie2_build,
    )
    logger.info("Genome index complete: {}".format(args.output))


def cmd_align(args):
    from methmap.align import align_single_end, align_paired_end

    extra = []
    if getattr(args, 'score_min', None):
        extra.extend(['--score-min', args.score_min])
    if getattr(args, 'no_discordant', False):
        extra.append('--no-discordant')
    if getattr(args, 'no_mixed', False):
        extra.append('--no-mixed')

    library_type = getattr(args, 'library_type', 'directional') or 'directional'

    kw = dict(
        genome_index_dir=args.index,
        output_dir=args.output,
        sample_name=args.sample,
        bowtie2_path=args.bowtie2,
        threads=args.threads,
        extra_bowtie2_args=extra or None,
    )

    if getattr(args, 'reads2', None):
        logger.info("Mode: Paired-end ({})".format(library_type))
        sam = align_paired_end(args.reads, args.reads2,
                               library_type=library_type, **kw)
    else:
        logger.info("Mode: Single-end")
        sam = align_single_end(args.reads, **kw)

    logger.info("Alignment complete: {}".format(sam))
    return sam


def cmd_dedup(args):
    from methmap.align import remove_duplicates
    library_type = getattr(args, 'library_type', 'directional') or 'directional'
    logger.info("=== MethMap Deduplication ({}) ===".format(library_type))
    dedup_sam = remove_duplicates(
        input_sam=args.sam,
        output_dir=getattr(args, 'output', None) or None,
        sample_name=getattr(args, 'sample', None) or None,
        samtools_path=args.samtools,
        threads=args.threads,
        library_type=library_type,
    )
    logger.info("Deduplication complete: {}".format(dedup_sam))
    return dedup_sam


def cmd_extract(args):
    from methmap.extract import MethylationExtractor

    logger.info("=== MethMap Methylation Extraction ===")

    snv_filter = getattr(args, 'snv_filter', False)
    snv_ratio  = getattr(args, 'snv_ratio', 0.2)

    if snv_filter:
        logger.info("SNV discrimination enabled (threshold: {})".format(snv_ratio))

    extractor = MethylationExtractor(
        genome_fasta=args.genome,
        min_mapq=args.min_mapq,
        min_base_qual=args.min_base_qual,
        snv_filter=snv_filter,
        snv_ratio=snv_ratio,
        trim_r1=getattr(args, 'trim_r1', 0) or 0,
        trim_r2=getattr(args, 'trim_r2', 0) or 0,
    )

    context_filter = getattr(args, 'context', None) or None
    if context_filter:
        logger.info("Context filter: {}".format(context_filter))

    extractor.process_sam(
        args.sam,
        context_filter=context_filter,
        samtools_path=getattr(args, 'samtools', 'samtools'),
    )
    extractor.print_summary()

    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)
    s   = args.sample
    cov = getattr(args, 'min_coverage', 1)
    exclude_snv = snv_filter

    if not getattr(args, 'no_bismark_cov', False):
        tmp = out / '{}.bismark.cov'.format(s)
        extractor.write_bismark_coverage(str(tmp), min_coverage=cov,
                                         exclude_snv=exclude_snv)
        gz = _bgzip(tmp)
        logger.info("Bismark coverage: {}".format(gz))

    if getattr(args, 'bedgraph', False):
        bg = out / '{}.bedGraph'.format(s)
        extractor.write_bedgraph(str(bg), min_coverage=cov, exclude_snv=exclude_snv)
        _bgzip(bg)

    if getattr(args, 'cytosine_report', False):
        rpt = out / '{}.CX_report.txt'.format(s)
        extractor.write_cytosine_report(str(rpt), min_coverage=0,
                                        exclude_snv=exclude_snv)
        _bgzip(rpt)

    if getattr(args, 'snv_report', False) and snv_filter:
        snv_rpt = out / '{}.SNV_report.txt'.format(s)
        extractor.write_snv_report(str(snv_rpt), min_coverage=cov)

    if getattr(args, 'vcf', False) and snv_filter:
        vcf_out = out / '{}.snv.vcf'.format(s)
        extractor.write_vcf(str(vcf_out), sample_name=s,
                            min_coverage=getattr(args, 'min_coverage', 10))
        logger.info("VCF: {}".format(vcf_out))

    if getattr(args, 'allc', False):
        # Write plain-text allc then bgzip → {sample}.allc.gz
        # (ALLCools/methylpy convention; tabix-indexable)
        allc_tmp = out / '{}.allc'.format(s)
        extractor.write_allc(str(allc_tmp), min_coverage=cov,
                             exclude_snv=exclude_snv)
        allc_out = _bgzip(allc_tmp)   # produces {sample}.allc.gz
        logger.info("allc (bgzip): {}".format(allc_out))

    # Methylation summary report (always written)
    report_out = out / '{}.methylation_report.txt'.format(s)
    extractor.write_methylation_report(
        str(report_out),
        sample_name=s,
        min_coverage=cov,
        exclude_snv=exclude_snv,
    )
    logger.info("Methylation report: {}".format(report_out))

    # Also print summary to console
    extractor.print_summary()
    logger.info("Extraction complete!")


def cmd_run(args):
    logger.info("=== MethMap Full Pipeline ===")

    # Build index if needed
    idx = Path(args.index)
    if not (idx / 'genome.1.bt2').exists() and \
       not (idx / 'genome.1.bt2l').exists():
        logger.info("Index not found, building...")
        cmd_prepare(args)

    # Align
    sam = cmd_align(args)

    # Dedup (optional)
    if getattr(args, 'remove_duplicates', False):
        from methmap.align import remove_duplicates
        sam = remove_duplicates(
            input_sam=sam,
            output_dir=args.output,
            sample_name=args.sample,
            samtools_path=getattr(args, 'samtools', 'samtools'),
            threads=args.threads,
            library_type=getattr(args, 'library_type', 'directional'),
        )

    # Extract
    extract_args = argparse.Namespace(
        genome=args.genome,
        sam=sam,
        output=args.output,
        sample=args.sample,
        min_mapq=getattr(args, 'min_mapq', 20),
        min_base_qual=getattr(args, 'min_base_qual', 20),
        min_coverage=getattr(args, 'min_coverage', 1),
        context=getattr(args, 'context', None),
        no_bismark_cov=getattr(args, 'no_bismark_cov', False),
        bedgraph=getattr(args, 'bedgraph', False),
        cytosine_report=getattr(args, 'cytosine_report', False),
        snv_filter=getattr(args, 'snv_filter', False),
        snv_ratio=getattr(args, 'snv_ratio', 0.2),
        snv_report=getattr(args, 'snv_report', False),
        vcf=getattr(args, 'vcf', False),
        allc=getattr(args, 'allc', False),
        samtools=getattr(args, 'samtools', 'samtools'),
        bgzip=getattr(args, 'bgzip', 'bgzip'),
        trim_r1=getattr(args, 'trim_r1', 0) or 0,
        trim_r2=getattr(args, 'trim_r2', 0) or 0,
    )
    cmd_extract(extract_args)
    logger.info("=== Pipeline Complete! Output: {} ===".format(args.output))


# ------------------------------------------------------------------ #
#  Argument parser                                                     #
# ------------------------------------------------------------------ #

def _add_snv_args(p):
    p.add_argument('--snv-filter',  action='store_true',
                   help='Enable SNV discrimination at CpG sites')
    p.add_argument('--snv-ratio',   type=float, default=0.2,
                   help='Minus-strand A-fraction threshold for SNV calling (default: 0.2)')
    p.add_argument('--snv-report',  action='store_true',
                   help='Write SNV report file (requires --snv-filter)')
    p.add_argument('--vcf',         action='store_true',
                   help='Write VCF file of SNV sites (requires --snv-filter)')
    p.add_argument('--allc',        action='store_true',
                   help='Write .allc.gz (bgzip, tabix-indexable, methylpy/ALLCools compatible)')


def main():
    parser = argparse.ArgumentParser(
        prog='methmap',
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    sub = parser.add_subparsers(dest='command', help='Command')
    sub.required = True

    # ---- prepare ----
    p = sub.add_parser('prepare', help='Build genome bowtie2 index')
    p.add_argument('--genome',        '-g', required=True,
                   help='Reference genome FASTA (original, unconverted)')
    p.add_argument('--output',        '-o', required=True,
                   help='Output index directory')
    p.add_argument('--bowtie2-build',       default='bowtie2-build',
                   help='Path to bowtie2-build (default: bowtie2-build)')
    p.set_defaults(func=cmd_prepare)

    # ---- align ----
    p = sub.add_parser('align', help='Align reads directly to genome')
    p.add_argument('--index',   '-x', required=True,
                   help='MethMap index directory (from methmap prepare)')
    p.add_argument('--reads',   '-1', required=True,
                   help='Input FASTQ (SE or R1)')
    p.add_argument('--reads2',  '-2', default=None,
                   help='R2 FASTQ (paired-end only)')
    p.add_argument('--output',  '-o', required=True,
                   help='Output directory')
    p.add_argument('--sample',  '-s', required=True,
                   help='Sample name')
    p.add_argument('--bowtie2',       default='bowtie2',
                   help='Path to bowtie2 (default: bowtie2)')
    p.add_argument('--threads', '-p', type=int, default=1)
    p.add_argument('--score-min',     default=None,
                   help='bowtie2 --score-min (default: bowtie2 default)')
    p.add_argument('--no-discordant', action='store_true')
    p.add_argument('--no-mixed',      action='store_true')
    p.add_argument('--library-type',  default='directional',
                   choices=['directional', 'non_directional'],
                   help='Library type: directional (default) or non_directional')
    p.set_defaults(func=cmd_align)

    # ---- dedup ----
    p = sub.add_parser('dedup', help='Remove PCR duplicates')
    p.add_argument('--sam',           required=True,
                   help='Aligned SAM file')
    p.add_argument('--output', '-o',  default=None,
                   help='Output directory (default: same as SAM)')
    p.add_argument('--sample', '-s',  default=None,
                   help='Sample name (default: SAM filename stem)')
    p.add_argument('--samtools',      default='samtools',
                   help='Path to samtools (default: samtools)')
    p.add_argument('--threads', '-p', type=int, default=1)
    p.add_argument('--library-type',  default='directional',
                   choices=['directional', 'non_directional'],
                   help='directional: PE dedup via fixmate+markdup; '
                        'non_directional: SE dedup for Tn5 (markdup -m s, no fixmate)')
    p.set_defaults(func=cmd_dedup)

    # ---- extract ----
    p = sub.add_parser('extract', help='Extract methylation from SAM/BAM')
    p.add_argument('--genome',        '-g', required=True)
    p.add_argument('--sam',                 required=True,
                   help='Aligned SAM or BAM (.bam auto-detected)')
    p.add_argument('--output',        '-o', required=True)
    p.add_argument('--sample',        '-s', required=True)
    p.add_argument('--samtools',            default='samtools',
                   help='Path to samtools (required for BAM input)')
    p.add_argument('--min-mapq',            type=int, default=20)
    p.add_argument('--min-base-qual',       type=int, default=20)
    p.add_argument('--min-coverage',        type=int, default=1)
    p.add_argument('--context',             nargs='+',
                   choices=['CpG', 'CHG', 'CHH'])
    p.add_argument('--no-bismark-cov',      action='store_true')
    p.add_argument('--bedgraph',            action='store_true')
    p.add_argument('--cytosine-report',     action='store_true')
    p.add_argument('--trim-r1',             type=int, default=0, metavar='N',
                   help='Skip first N bp of R1 (forward) reads. P5-end Tn5 bias (default: 0)')
    p.add_argument('--trim-r2',             type=int, default=0, metavar='N',
                   help='Skip first N bp of R2 (reverse) reads. P7-end Tn5 bias (default: 0)')
    p.add_argument('--threads',    '-p',    type=int, default=1,
                   help='Threads (currently unused in extract, accepted for pipeline compatibility)')
    p.add_argument('--bgzip',               default='bgzip',
                   help='Path to bgzip for BGZF-compressed allc output (default: bgzip)')
    _add_snv_args(p)   # adds --snv-filter --snv-ratio --snv-report --vcf --allc
    p.set_defaults(func=cmd_extract)

    # ---- run ----
    p = sub.add_parser('run', help='Run full pipeline')
    p.add_argument('--genome',        '-g', required=True)
    p.add_argument('--index',         '-x', required=True)
    p.add_argument('--reads',         '-1', required=True)
    p.add_argument('--reads2',        '-2', default=None)
    p.add_argument('--output',        '-o', required=True)
    p.add_argument('--sample',        '-s', required=True)
    p.add_argument('--bowtie2',             default='bowtie2')
    p.add_argument('--bowtie2-build',       default='bowtie2-build')
    p.add_argument('--samtools',            default='samtools')
    p.add_argument('--bgzip',               default='bgzip',
                   help='Path to bgzip for BGZF allc output')
    p.add_argument('--threads',       '-p', type=int, default=1)
    p.add_argument('--remove-duplicates',   action='store_true')
    p.add_argument('--min-mapq',            type=int, default=20)
    p.add_argument('--min-base-qual',       type=int, default=20)
    p.add_argument('--min-coverage',        type=int, default=1)
    p.add_argument('--context',             nargs='+',
                   choices=['CpG', 'CHG', 'CHH'])
    p.add_argument('--no-bismark-cov',      action='store_true')
    p.add_argument('--bedgraph',            action='store_true')
    p.add_argument('--cytosine-report',     action='store_true')
    p.add_argument('--score-min',           default=None)
    p.add_argument('--no-discordant',       action='store_true')
    p.add_argument('--no-mixed',            action='store_true')
    p.add_argument('--library-type',        default='directional',
                   choices=['directional', 'non_directional'],
                   help='Library type: directional (default) or non_directional')
    p.add_argument('--trim-r1',             type=int, default=0, metavar='N',
                   help='Skip first N bp of R1 reads (P5-end Tn5 bias, default: 0)')
    p.add_argument('--trim-r2',             type=int, default=0, metavar='N',
                   help='Skip first N bp of R2 reads (P7-end Tn5 bias, default: 0)')
    _add_snv_args(p)
    p.set_defaults(func=cmd_run)

    args = parser.parse_args()
    try:
        args.func(args)
    except KeyboardInterrupt:
        logger.info("Interrupted")
        sys.exit(1)
    except Exception as e:
        logger.error(str(e))
        raise


if __name__ == '__main__':
    main()
