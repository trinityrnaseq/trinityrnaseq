from Bio import SeqIO
import unittest
import os

class TestTrinitySampleData(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.sampledata_dir = os.environ["TRINITY_SAMPLEDATA"]

    def test_genome_guided(self):
        seq_count = self.count_sequences('test_GenomeGuidedTrinity', 'test_GG_use_bam_trinity_outdir',
                      'Trinity-GG.fasta')
        self.assertTrue(50 <= seq_count <= 60, msg='Found %s sequences' % seq_count)

    def test_genome_guided_with_jaccard_clipping(self):
        seq_count = self.count_sequences('test_GenomeGuidedTrinity', 'test_Schizo_trinityGG_jaccard_RF_outdir',
                            'Trinity-GG.fasta')
        self.assertTrue(60 <= seq_count <= 80, msg='Found %s sequences' % seq_count)

    def test_paired_end_normalization(self):
        seq_count = self.count_sequences('test_InSilicoReadNormalization',
                          'reads.left.fq.gz.normalized_K25_C5_pctSD200.fq')
        self.assertTrue(35 <= seq_count <= 50, msg='Found %s sequences' % seq_count)
        seq_count = self.count_sequences('test_InSilicoReadNormalization',
                          'reads.right.fq.gz.normalized_K25_C5_pctSD200.fq')
        self.assertTrue(30 <= seq_count <= 40, msg='Found %s sequences' % seq_count)
        seq_count = self.count_sequences('test_InSilicoReadNormalization',
                          'reads.single.fq.normalized_K25_C5_pctSD200.fq')
        self.assertTrue(50 <= seq_count <= 65, msg='Found %s sequences' % seq_count)

    def test_trinity_assembly(self):
        seq_count = self.count_sequences('test_Trinity_Assembly', 'trinity_out_dir', 'Trinity.fasta')
        self.assertTrue(100 <= seq_count <= 120, msg='Found %s sequences' % seq_count)

    def test_DE_analysis_EdgeR(self):
        check_file = os.path.join(self.sampledata_dir, 'test_DE_analysis', 'edgeR_outdir', 'numDE_feature_counts.P0.001_C2.matrix')
        self.assertTrue(os.path.isfile(check_file))

    def test_align_and_estimate_abundance(self):
        cats = ['PAIRED', 'SINGLE']
        for cat in cats:
            check_file = os.path.join(self.sampledata_dir, 'test_align_and_estimate_abundance',\
               '%s_END_ABUNDANCE_ESTIMATION' % cat, 'RSEM-gene.counts.matrix')
            self.assertTrue(os.path.isfile(check_file))

    def test_full_edgeR_pipeline(self):
        check_file = os.path.join(self.sampledata_dir, 'test_full_edgeR_pipeline', 'read_content_analysis', 'read_content_analysis.nameSorted.bam')
        self.assertTrue(os.path.isfile(check_file))


### Helper methods
    def count_sequences(self, *paths):
        gg = os.path.join(self.sampledata_dir, *paths)
        handle = open(gg, "rU")
        seq_count = len([x for x in SeqIO.parse(handle, "fasta")])
        handle.close()
        return seq_count

