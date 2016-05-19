import subprocess
from Bio import SeqIO
import unittest
import shutil
import os
import time
import filecmp

# Prereqs:
# module load bowtie/0.12.8
# module load java
# module load samtools
# Trinity
# Copy the .gz files in sample_data/test_Trinity_Assembly to current directory
# Run using nosetests
MEM_FLAG = "--max_memory 2G"
TEMP_FILES = ['both.fa', 'inchworm.K25.L25.fa', 'jellyfish.kmers.fa']


class TestTrinity(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        try:
            os.remove('coverage.log')
        except:
            pass

    def tearDown(self):
        shutil.rmtree('trinity_out_dir', True)

    def test_sample_data_seq_count(self):
        print "When assembling the sample data, the number of sequences assembled should be between 75 and 100"
        self.trinity(
            "Trinity --seqType fq %s --left reads.left.fq.gz,reads2.left.fq.gz --right reads.right.fq.gz,reads2.right.fq.gz --SS_lib_type RF --CPU 4 --no_cleanup" % MEM_FLAG)
        handle = open("trinity_out_dir/Trinity.fasta", "rU")
        seq_count = len([x for x in SeqIO.parse(handle, "fasta")])
        handle.close()
        self.assertTrue(85 <= seq_count <= 110, msg='Found %s sequences' % seq_count)

    def test_sample_data_trimmed_and_normalized(self):
        print "When assembling the sample data with the --trimmomatic --normalize_reads flags, the number of sequences assembled should be between 75 and 85"
        self.trinity(
            "Trinity --seqType fq %s --left reads.left.fq.gz,reads2.left.fq.gz --right reads.right.fq.gz,reads2.right.fq.gz --SS_lib_type RF --CPU 4 --trimmomatic --normalize_reads --no_cleanup" % MEM_FLAG)
        handle = open("trinity_out_dir/Trinity.fasta", "rU")
        seq_count = len([x for x in SeqIO.parse(handle, "fasta")])
        handle.close()
        self.assertTrue(85 <= seq_count <= 100, msg='Found %s sequences' % seq_count)

    def test_no_cleanup_leaves_temp_files(self):
        print "The --no_cleanup flag should ensure that the output directory is left behind"
        self.trinity(
            "Trinity --seqType fq %s --left reads.left.fq.gz,reads2.left.fq.gz --right reads.right.fq.gz,reads2.right.fq.gz --SS_lib_type RF --CPU 4 --no_cleanup" % MEM_FLAG)
        for f in TEMP_FILES:
            self.assertTrue(os.path.exists("trinity_out_dir/%s" % f), msg="%s not found with no_cleanup" % f)

    def test_cleanup_removes_temp_files(self):
        print "The --full_cleanup flag should ensure that the output directory is gone but the output file remains"
        self.trinity(
            "Trinity --seqType fq %s --left reads.left.fq.gz,reads2.left.fq.gz --right reads.right.fq.gz,reads2.right.fq.gz --SS_lib_type RF --CPU 4 --full_cleanup" % MEM_FLAG)
        time.sleep(5) # Make sure the system has time to recognize the directory is gone
        self.assertFalse(os.path.exists("trinity_out_dir"), msg="Did full_cleanup but trinity_out_dir exists")
        self.assertTrue(os.path.isfile("trinity_out_dir.Trinity.fasta"),
                            msg="Did full_cleanup but output file not created")

    def test_single_end_with_rf_lib_type_error(self):
        print "Single reads with an SS_lib_type of RF should result in an error"
        try:
            subprocess.call("Trinity --seqType fq --single reads.left.fq --SS_lib_type RF", shell=True)
        except subprocess.CalledProcessError as e:
            self.assertTrue("Error, with --single reads, the --SS_lib_type can be 'F' or 'R' only." in e.output)

    def test_single_end_with_fq(self):
        print "Single reads with FQ file should succeed"
        self.trinity("Trinity %s --seqType fq --single reads.left.fq --SS_lib_type F" % MEM_FLAG)

    def test_no_run_chrysalis(self):
        print "The --no_run_chrysalis flag should result in no chrysalis subdirectory in the output directory"
        self.trinity("Trinity %s --seqType fq --single reads.left.fq --SS_lib_type F --no_run_chrysalis" % MEM_FLAG)
        self.assertEquals(0, len(os.listdir('trinity_out_dir/chrysalis')))

    def test_no_run_inchworm(self):
        print "The --no_run_inchworm flag should result in no inchworm.finished file"
        self.trinity("Trinity %s --seqType fq --single reads.left.fq --SS_lib_type F --no_run_inchworm" % MEM_FLAG)
        self.assertFalse(os.path.isfile("trinity_out_dir/inchworm.K25.L25.fa.finished"),
                            msg="Inchworm appears to have run although no_run_inchworm was specified")
        self.assertTrue(os.path.isfile("trinity_out_dir/jellyfish.kmers.fa"),
                            msg="jellyfish.kmers.fa was not created")

    def test_no_bowtie(self):
        print "The --no_bowtie flag should result in no bowtie.nameSorted.bam file"
        self.trinity("Trinity %s --seqType fq --single reads.left.fq --SS_lib_type F --no_bowtie" % MEM_FLAG)
        self.assertFalse(os.path.isfile("trinity_out_dir/bowtie.nameSorted.bam"),
                            msg="Bowtie appears to have run although no_bowtie was specified")

    def test_no_distributed_trinity_exec(self):
        print "The --no_distributed_trinity_exec flag should run Jellyfish but not create an output file"
        self.trinity("Trinity %s --seqType fq --single reads.left.fq --SS_lib_type F --no_distributed_trinity_exec" % MEM_FLAG)
        self.assertTrue(os.path.isfile("trinity_out_dir/inchworm.K25.L25.fa.finished"),
                            msg="Inchworm did not appear to run with no_distributed_trinity_exec flag")
        self.assertTrue(os.path.isfile("trinity_out_dir/jellyfish.kmers.fa.histo"),
                            msg="Jellyfish did not appear to run with no_distributed_trinity_exec flag")
        self.assertFalse(os.path.isfile("trinity_out_dir/Trinity.fasta"),
                            msg="Trinity.fasta created with no_distributed_trinity_exec")

    def test_single_end_with_fa_and_reverse(self):
        print "The --no_distributed_trinity_exec flag should run Jellyfish but not create an output file"
        self.fq2fa()
        self.trinity("Trinity %s --seqType fa --single reads.fa --SS_lib_type R" % MEM_FLAG)

    def test_output_correctly_changes_dir(self):
        print "The --output flag should change the output directory"
        shutil.rmtree('trinity_test', True)
        self.trinity("Trinity %s --seqType fq --single reads.left.fq --SS_lib_type F --output trinity_test" % MEM_FLAG)
        self.assertTrue(os.path.exists("trinity_test"), msg="Changed output directory but it was not created")
        shutil.rmtree('trinity_test', True)

    def test_scaffold_iworm_contigs(self):
        print "scaffold_iworm_contigs works as expected"
	os.environ['LD_LIBRARY_PATH'] = os.environ['LD_LIBRARY_PATH'] + ':../src/trinity-plugins/htslib'
        exe = "../src/trinity-plugins/scaffold_iworm_contigs/scaffold_iworm_contigs"
        bamfile = "iworm.bowtie.nameSorted.bam"
        ifile = "inchworm.K25.L25.fa"
        f = subprocess.check_output([exe, bamfile, ifile]).split('\n')
        expected_result = [['a340;25', '339', 'a9;40', '8', '41'],
            ['a719;8', '718', 'a832;15', '831', '33'],
            ['a1;43', '0', 'a346;23', '345', '31'],
            ['a346;23', '345', 'a9;40', '8', '26'],
            ['a339;142', '338', 'a37;14', '36', '25'],
            ['a3;61', '2', 'a432;9', '431', '23'],
            ['a345;34', '344', 'a40;12', '39', '21'],
            ['a354;96', '353', 'a368;13', '367', '18'],
            ['a689;4', '688', 'a774;6', '773', '13']]

        actual_result = [line.split('\t') for line in f if line]
        self.assertEquals(expected_result, actual_result[0:9])
        actual_lengths = [int(s[4].strip()) for s in actual_result]
        expected_order = list(reversed(sorted(actual_lengths)))
        self.assertEquals(expected_order, actual_lengths)

    def test_Inchworm_handles_compressed_files(self):
        print "A compressed single file should be handlred correctly by Inchworm"
        self.trinity('Trinity %s --seqType fq --single reads.left.fq.gz --SS_lib_type F --no_run_chrysalis' % MEM_FLAG);
        num_lines = sum(1 for line in open('trinity_out_dir/inchworm.K25.L25.fa'))
        self.assertTrue(2850 <= num_lines <= 3100, msg='Found %s lines' % num_lines)

    def test_QuantifyGraph_works(self):
        exe = "../src/Chrysalis/QuantifyGraph"
        args = " -g Cbin0/c0.graph.tmp -i Cbin0/c0.reads.tmp -o c0.graph.out -max_reads 200000 -k 24"
        result = subprocess.call(exe + args,shell=True)
        self.assertEquals(0, result, "QuantifyGraph failed")
        self.assertTrue(filecmp.cmp("c0.graph.out", "Cbin0/c0.graph.out.master", shallow=False), "Unexpected result from QuantifyGraph")

### information tests
    def test_cite(self):
        print "Cite flag should return citing information"
        expected = '\n\n* Trinity:\nFull-length transcriptome assembly from RNA-Seq data without a reference genome.\nGrabherr MG, Haas BJ, Yassour M, Levin JZ, Thompson DA, Amit I, Adiconis X, Fan L,\nRaychowdhury R, Zeng Q, Chen Z, Mauceli E, Hacohen N, Gnirke A, Rhind N, di Palma F,\nBirren BW, Nusbaum C, Lindblad-Toh K, Friedman N, Regev A.\nNature Biotechnology 29, 644\xe2\x80\x93652 (2011)\nPaper: http://www.nature.com/nbt/journal/v29/n7/full/nbt.1883.html\nCode:  http://trinityrnaseq.sf.net\n\n\n'
        cite = subprocess.check_output(["Trinity", "--cite"])
        self.assertEqual(expected, cite)

    def test_version(self):
        print "Version flag should return version information"
        try:
            subprocess.check_output(["Trinity", "--version"])
            self.fail("Version returned 0 errorcode!")
        except subprocess.CalledProcessError as e:
            self.assertTrue('Trinity version: __TRINITY_VERSION_TAG__' in e.output)
            self.assertTrue('using Trinity devel version. Note, latest production release is: v2.1.1' in e.output)


    def test_show_full_usage_info(self):
        print "show_full_usage_info flag has several option sections"
        try:
            subprocess.check_output(["Trinity", "--show_full_usage_info"])
        except subprocess.CalledProcessError as e:
            self.assertTrue("Inchworm and K-mer counting-related options" in e.output)
            self.assertTrue("Chrysalis-related options" in e.output)
            self.assertTrue("Butterfly-related options" in e.output)
            self.assertTrue("Quality Trimming Options" in e.output)
            self.assertTrue("In silico Read Normalization Options" in e.output)

### Invalid command line tests
    def test_no_JM_specified_error(self):
        print "max_memory flag is required"
        error = self.get_error("Trinity --seqType fq --single reads.left.fq --SS_lib_type F")
        self.assertTrue("Error, must specify max memory for jellyfish to use, eg.  --max_memory 10G" in error)

    def test_invalid_option_error(self):
        print "Invalid options result in an error"
        error = self.get_error("Trinity --squidward")
        self.assertTrue("ERROR, don't recognize parameter: --squidward" in error)

    def test_set_no_cleanup_and_full_cleanup_error(self):
        print "Setting no_cleanup and full_cleanup together results in an error"
        error = self.get_error("Trinity --no_cleanup --full_cleanup")
        self.assertTrue("cannot set --no_cleanup and --full_cleanup as they contradict" in error)


### Helper methods
    def trinity(self, cmdline):
        print "Command line:", cmdline
        with open("coverage.log", 'a') as file_out:
            file_out.write("COMMAND: %s\n" % cmdline) 
            file_out.flush()
            subprocess.call(cmdline,shell=True, stdout=file_out)
            file_out.write("TEST COMPLETE\n") 

    def get_error(self, cmd):
        try:
            subprocess.check_output(cmd.split(' '))
        except subprocess.CalledProcessError as e:
            return e.output

    def fq2fa(self):
        handle = open("reads.left.fq", "rU")
        records = [x for x in SeqIO.parse(handle, "fastq")]
        handle.close()
        SeqIO.write(records, "reads.fa", "fasta")

