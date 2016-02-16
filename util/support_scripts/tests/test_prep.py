import unittest
import shutil
import os
import subprocess
from Bio import SeqIO

# single, leftright
# fasta, fastq
# one file, two files
# forward, reverse
# clear, gzip, bzip


class TestTrinityPrepFlag(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        try:
            os.remove('coverage.log')
        except:
            pass

    def tearDown(self):
        shutil.rmtree('trinity_out_dir', True)
        #pass

    def test_fastq(self):
        self.trinity("left1.fq", "fq")
        self.assertEquals(30575, self.count_seqs(), "Unexpected sequence count")

    def test_fastq_gz(self):
        self.trinity("left1.fq.gz", "fq")
        self.assertEquals(30575, self.count_seqs(), "Unexpected sequence count")

    def test_fastq_bz2(self):
        self.trinity("left1.fq.bz2", "fq")
        self.assertEquals(30575, self.count_seqs(), "Unexpected sequence count")

    def test_fastq_multiple_files_single(self):
        self.trinity("left1.fq,left1.fq.gz", "fq")
        self.assertEquals(61150, self.count_seqs(), "Unexpected sequence count")

    def test_fastq_multiple_files_single_bz2(self):
        self.trinity("left1.fq.bz2,left1.fq.gz", "fq")
        self.assertEquals(61150, self.count_seqs(), "Unexpected sequence count")

    def test_fastq_multiple_files_single_reverse(self):
        self.trinity("left1.fq,left1.fq.gz", "fq", True)
        self.assertEquals(61150, self.count_seqs(), "Unexpected sequence count")

    def test_fasta(self):
        self.trinity("left1.fa")
        self.assertEquals(30575, self.count_seqs(), "Unexpected sequence count")

    def test_fasta_gz(self):
        self.trinity("left1.fa.gz")
        self.assertEquals(30575, self.count_seqs(), "Unexpected sequence count")

    def test_fasta_multiple_files_single(self):
        self.trinity("left1.fa,left1.fa.gz")
        self.assertEquals(61150, self.count_seqs(), "Unexpected sequence count")

    def test_fasta_multiple_files_single_reverse(self):
        self.trinity("left1.fa,left1.fa.gz", reverse=True)
        self.assertEquals(61150, self.count_seqs(), "Unexpected sequence count")

    def test_paired_fastq(self):
        self.trinity("left1.fq", "fq", morefiles="right1.fq")
        self.assertEquals(61150, self.count_seqs(), "Unexpected sequence count")

    def test_paired_fastq_gz(self):
        self.trinity("left1.fq.gz", "fq", morefiles="right1.fq.gz")
        self.assertEquals(61150, self.count_seqs(), "Unexpected sequence count")

    def test_fastq_multiple_files_paired(self):
        self.trinity("left1.fq,left1.fq.gz", "fq", morefiles="right1.fq,right1.fq.gz")
        self.assertEquals(122300, self.count_seqs(), "Unexpected sequence count")

    def test_fastq_multiple_files_paired_reverse(self):
        self.trinity("left1.fq,left1.fq.gz", "fq", reverse=True, morefiles="right1.fq,right1.fq.gz")
        self.assertEquals(122300, self.count_seqs(), "Unexpected sequence count")

    def test_fasta_paired(self):
        self.trinity("left1.fa", morefiles="right1.fa")
        self.assertEquals(61150, self.count_seqs(), "Unexpected sequence count")

    def test_paired_sequences_have_1_or_2_extension(self):
        self.trinity("sra_test.fq", morefiles="sra_test2.fq", seqtype='fq')
        self.assertEquals(0, self.count_bad_endings(), "Found sequences with bad endings")

    def test_fasta_gz_paired(self):
        self.trinity("left1.fa.gz", morefiles="right1.fa.gz")
        self.assertEquals(61150, self.count_seqs(), "Unexpected sequence count")

    def test_fasta_multiple_files_paired(self):
        self.trinity("left1.fa,left1.fa.gz", morefiles="right1.fa,right1.fa.gz")
        self.assertEquals(61150, self.count_seqs(), "Unexpected sequence count")

    def test_fasta_multiple_files_paired(self):
        self.trinity("left1.fa,left1.fa.gz", morefiles="right1.fa,right1.fa.gz", reverse=True)
        self.assertEquals(122300, self.count_seqs(), "Unexpected sequence count")

    def trinity(self, files, seqtype='fa', reverse=False, morefiles=None):
        if morefiles:
            tpl = "Trinity --left %s --right %s --prep --seqType %s --max_memory 2G --no_version_check"
            cmdline = tpl % (files, morefiles, seqtype)
        else:
            tpl = "Trinity --single %s --prep --seqType %s --max_memory 2G --no_version_check"
            cmdline = tpl % (files, seqtype)
        if reverse:
          cmdline += " --SS_lib_type " + ('RF' if morefiles else 'R')
        print "Command line:", cmdline
        with open("coverage.log", 'a') as file_out:
            subprocess.call(cmdline,shell=True, stdout=file_out)

    def count_seqs(self):
        f = "trinity_out_dir/single.fa"
        if os.path.isfile(f):
          handle = open(f, "rU")
        else:
          handle = open("trinity_out_dir/both.fa", "rU")

        seq_count = len([x for x in SeqIO.parse(handle, "fasta")])
        handle.close()
        return seq_count

    def count_bad_endings(self):
        f = "trinity_out_dir/single.fa"
        if os.path.isfile(f):
          handle = open(f, "rU")
        else:
          handle = open("trinity_out_dir/both.fa", "rU")

        seq_count = len(list(x for x in SeqIO.parse(handle, "fasta") if not (x.id.endswith('/1') or x.id.endswith('/2'))))
        handle.close()
        return seq_count

