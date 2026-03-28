#!/usr/bin/env python3
"""MethMap unit tests (no bowtie2 required)"""

import sys
import os
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from methmap.genome_prepare import (
    convert_genome_CT, convert_genome_GA,
    read_fasta, write_fasta, reverse_complement,
)
from methmap.align import convert_reads_CT, convert_reads_GA
from methmap.extract import get_context, parse_cigar, cigar_to_alignment


class TestGenomeConversion(unittest.TestCase):

    def setUp(self):
        self.seqs = {
            'chr1': 'ACGCTACGCTAGCGT',
            'chr2': 'CCCCGGGGTTTTAAAA',
        }

    def test_ct_no_c(self):
        ct = convert_genome_CT(self.seqs)
        for chrom, seq in ct.items():
            self.assertNotIn('C', seq)

    def test_ga_no_g(self):
        ga = convert_genome_GA(self.seqs)
        for chrom, seq in ga.items():
            self.assertNotIn('G', seq)

    def test_fasta_roundtrip(self):
        with tempfile.NamedTemporaryFile(suffix='.fa', delete=False) as f:
            tmp = f.name
        write_fasta(self.seqs, tmp)
        loaded = read_fasta(tmp)
        os.unlink(tmp)
        self.assertEqual(self.seqs, loaded)

    def test_reverse_complement(self):
        self.assertEqual(reverse_complement('ACGT'), 'ACGT')
        self.assertEqual(reverse_complement('AACCGG'), 'CCGGTT')


class TestReadConversion(unittest.TestCase):

    def _write(self, reads, path):
        with open(path, 'w') as f:
            for i, s in enumerate(reads):
                f.write('@r{}\n{}\n+\n{}\n'.format(i, s, 'I' * len(s)))

    def _read(self, path):
        seqs = []
        with open(path) as f:
            while True:
                if not f.readline():
                    break
                seqs.append(f.readline().strip())
                f.readline(); f.readline()
        return seqs

    def test_ct(self):
        with tempfile.NamedTemporaryFile(suffix='.fq', delete=False) as f:
            inp = f.name
        with tempfile.NamedTemporaryFile(suffix='.fq', delete=False) as f:
            out = f.name
        self._write(['ACGCT', 'CCCTGGA'], inp)
        convert_reads_CT(inp, out)
        result = self._read(out)
        os.unlink(inp); os.unlink(out)
        self.assertEqual(result[0], 'ATGTT')
        self.assertEqual(result[1], 'TTTTGGA')

    def test_ga(self):
        with tempfile.NamedTemporaryFile(suffix='.fq', delete=False) as f:
            inp = f.name
        with tempfile.NamedTemporaryFile(suffix='.fq', delete=False) as f:
            out = f.name
        self._write(['ACGCT', 'GGGCCC'], inp)
        convert_reads_GA(inp, out)
        result = self._read(out)
        os.unlink(inp); os.unlink(out)
        self.assertEqual(result[0], 'ACACT')
        self.assertEqual(result[1], 'AAACCC')


class TestContext(unittest.TestCase):

    def test_cpg_plus(self):
        self.assertEqual(get_context('ACGT', 1, '+'), 'CpG')

    def test_chg_plus(self):
        self.assertEqual(get_context('ACTG', 1, '+'), 'CHG')

    def test_chh_plus(self):
        self.assertEqual(get_context('ACTA', 1, '+'), 'CHH')

    def test_cpg_minus(self):
        self.assertEqual(get_context('ACGT', 2, '-'), 'CpG')


class TestCigar(unittest.TestCase):

    def test_parse(self):
        self.assertEqual(parse_cigar('5M2I3M'), [(5,'M'),(2,'I'),(3,'M')])

    def test_alignment(self):
        ref  = 'ACGTACGT'
        read = 'ACTTACGT'
        aln  = cigar_to_alignment('8M', read, 0, ref)
        self.assertEqual(len(aln), 8)
        rp, rb, fb = aln[2]
        self.assertEqual(fb, 'G')
        self.assertEqual(rb, 'T')


class TestMethylationLogic(unittest.TestCase):

    def test_meth_plus(self):
        """ref=C read=T -> METHYLATED in new chemistry"""
        aln = cigar_to_alignment('4M', 'ATGT', 0, 'ACGT')
        _, rb, fb = aln[1]
        self.assertEqual(fb, 'C')
        self.assertEqual(rb, 'T')
        self.assertTrue(fb == 'C' and rb == 'T')   # methylated

    def test_unmeth_plus(self):
        """ref=C read=C -> UNMETHYLATED in new chemistry"""
        aln = cigar_to_alignment('4M', 'ACGT', 0, 'ACGT')
        _, rb, fb = aln[1]
        self.assertEqual(fb, 'C')
        self.assertEqual(rb, 'C')
        self.assertTrue(fb == 'C' and rb == 'C')   # unmethylated

    def test_meth_minus(self):
        """ref=G read=A -> METHYLATED on minus strand"""
        aln = cigar_to_alignment('4M', 'ACAT', 0, 'ACGT')
        _, rb, fb = aln[2]
        self.assertEqual(fb, 'G')
        self.assertEqual(rb, 'A')
        self.assertTrue(fb == 'G' and rb == 'A')   # methylated


if __name__ == '__main__':
    unittest.main(verbosity=2)
