__author__ = 'befulton'

import unittest
import make_data_files

def generate(lines):
    for line in lines:
        yield line

class CollectlTests(unittest.TestCase):

    def test_build_datasets_returns_last_line_time_as_end_time(self):
        lines = ["20140919 13:10:49\n",
                 "20140919 11:10:49\n",
                 "20140919 11:27:49\n"]
        _, _, _, _, end_time = make_data_files.build_datasets(generate(lines))
        self.assertEquals(['20140919', '11:27:49'], end_time)

    def test_build_datasets_adds_inchworm_to_line_dict(self):
        lines = ["20140919 11:08:49 44009 user 20 44008 0 R 10975748 0 10837764 10957692 92 128 3416 9 0.00 4.99 99 "
                 "01:33:35 0 0 0 0 0 0 0 0 1 /N/test/Inchworm/bin/inchworm --kmers jellyfish.kmers.fa --run_inchworm "
                 "-K 25 -L 25 --monitor 1"]
        _, line_dict, _, _, _ = make_data_files.build_datasets(generate(lines))
        self.assertEquals(lines, line_dict['inchworm'])

    def test_build_datasets_adds_bowtie_build_to_line_dict_if_no_path(self):
        lines = ["20140919 11:11:24 44143 user 20 44131 0 R 547552 0 502040 535608 164 260 3244 9 0.01 4.98 99 "
                 "02:27.86 0 960 0 960 0 124 0 0 0 bowtie-build /N/test/trinity_out/jaccard_clip_workdir/iworm.fa "
                 "/N/trinity_out/jaccard_clip_workdir/iworm.fa",
                 "20140919 11:11:29 44143 user 20 44131 0 R 547552 0 502040 535608 164 260 3244 9 0.01 4.98 99 02:27.86"
                 " 0 960 0 960 0 124 0 0 0 bowtie-build /N/trinity_out/jaccard_clip_workdir/iworm.fa "
                 "/N/trinity_out/jaccard_clip_workdir/iworm.fa",
                 "20140919 11:11:34 44143 user 20 44131 0 R 547552 0 502040 535608 164 260 3244 9 0.01 4.98 99 02:27.86"
                 " 0 960 0 960 0 124 0 0 0 some_app /N/trinity_out/bowtie-build/jaccard_clip_workdir/iworm.fa "
                 "/N/trinity_out/jaccard_clip_workdir/iworm.fa"]
        _, line_dict, _, _, _ = make_data_files.build_datasets(generate(lines))
        self.assertEquals(lines[0:2], line_dict['bowtie-build'])
        
    def test_Chrysalis_is_identified_correctly(self):
        lines = ["20140920 11:29:54 49383 befulton 20 43684 0 S 18468 0 1464 300 92 236 3416 8 0.00 0.00 0 00:02.44 "
                 "0 0 0 0 0 0 0 0 0 /trinity/Chrysalis/Chrysalis "
                 "-i both.fa -iworm /jobdir/trinity_out/inchworm.K25.L25.fa.clipped.fa -o "
                 "/jobdir/trinity_out/chrysalis -cpu 16 -min_glue 2 "
                 "-min_iso_ratio 0.05 -glue_factor 0.05 -weldmer_size 48 -min 200 -dist 500 -max_reads 200000 "
                 "-sort_buffer_size 20G -max_mem_reads 1000000 -strand 1 -paired -reads_for_pairs both.fa "
                 "-butterfly /trinity/Butterfly/Butterfly.jar"]
        _, line_dict, _, _, _ = make_data_files.build_datasets(generate(lines))
        self.assertEquals(lines, line_dict['Chrysalis'])

    def test_Butterfly_is_identified_correctly(self):
        lines = ["20140920 17:28:14 30084 befulton 20 62435 0 S 24215176 0 37768 24080384 92 4 14420 30 0.03 0.17 "
                 "4 00:00.20 0 6 519 0 203 1 0 0 1056 java -Xmx20G -Xms1G -jar /trinity/Butterfly/Butterfly.jar "
                 "-N 100000 -L 200 -F 500 -C /jobdir/trinity_out/chrysalis/Component_bins/Cbin5/c6058.graph "
                 "--max_number_of_paths_per_node=10 --path_reinforcement_distance=75 --triplet-lock"]
        _, line_dict, _, _, _ = make_data_files.build_datasets(generate(lines))
        self.assertEquals(lines, line_dict['Butterfly'])

    def test_shell_c_commands_are_dropped(self):
        lines = ["20140919 10:06:14 43812 befulton 20 43684 0 S 108160 0 1228 196 92 848 1832 1 0.00 0.00 0 00:00.00 " \
                "0 0 0 0 0 0 0 0 0 sh -c /trinity/trinity-plugins/fastool/fastool --rev  --illumina-trinity " \
                "--to-fasta /jobdir/reads.left.fq >> left.fa"]
        _, line_dict, _, _, _ = make_data_files.build_datasets(generate(lines))
        self.assertEquals(0, len(line_dict))

    def test_prettyprocess_combines_sort_and_binsort(self):
        line = "20140919 11:08:49 44009 user 20 44008 0 R 10975748 0 10837764 10957692 92 128 3416 9 0.00 4.99 99 " \
               "01:33:35 0 0 0 0 0 0 0 0 1 /bin/sort"
        self.assertEquals("sort", make_data_files.prettyprocess(line))
        line = "20140919 11:08:49 44009 user 20 44008 0 R 10975748 0 10837764 10957692 92 128 3416 9 0.00 4.99 99 " \
               "01:33:35 0 0 0 0 0 0 0 0 1 sort"
        self.assertEquals("sort", make_data_files.prettyprocess(line))

    def test_sort(self):
        line = "20140919 16:06:24 45210 befulton 20 45209 0 R 973936 0 966268 965740 92 88 1824 10 1.96 0.94 58 " \
               "04:26.85 0 188211 178646 188416 82 92 0 0 0 sort -T . -S 2G -k 1,1 -k 3,3 right_fa.sam"
        self.assertEquals("sort", make_data_files.prettyprocess(line))

    def test_GraphFromFasta_identified(self):
        line = "20141008 04:38:10 15672 befulton 20 15671 0 R 91636 0 11216 73412 92 292 3416 1 0.01 5.18 103" \
               " 00:05.19 0 0 1342 0 337 3 0 0 579 " \
               "/N/dc2/scratch/befulton/TrinityMason/trinityrnaseq_r20140717/Chrysalis/GraphFromFasta " \
               "-i /dev/shm/trinity72/trinity_out/inchworm.K25.L25.fa -r both.fa -min_contig_length 200 -min_glue 2 " \
               "-glue_factor 0.05 -min_iso_ratio 0.05 -t 32 -k 24 -kk 48 -strand -scaffolding iworm_scaffolds.txt"

        self.assertEquals("GraphFromFasta", make_data_files.prettyprocess(line))

    def test_Samtools_identified_with_verb(self):
        line = "20141008 01:51:30 14984 befulton 20 14982 0 S 17752 0 820 200 92 384 2664 2 0.00 0.00 0 00:00.00 99" \
               " 0 1 0 2 0 0 0 61 samtools view -F4 -Sb -"
        self.assertEquals("samtools_view", make_data_files.prettyprocess(line))
        line = "20141008 01:51:30 14985 befulton 20 14982 0 S 17748 0 740 196 92 384 2664 16 0.00 0.00 0 00:00.00 26" \
               " 0 1 0 2 0 0 0 59 samtools sort -no - -"
        self.assertEquals("samtools_sort", make_data_files.prettyprocess(line))

    def test_Jellyfish_identified(self):
        line = "20141008 00:03:55 14423 befulton 20 14226 0 S 11316512 0 9153040 11314672 92 1732 0 17 3.14 19.75" \
               " 457 00:22.89 0 0 26047 0 3257 0 0 10 1129 " \
               "/N/dc2/scratch/befulton/TrinityMason/trinityrnaseq_r20140717/trinity-plugins/jellyfish/bin/jellyfish " \
               "count -t 32 -m 25 -s 1380200998 both.fa"
        self.assertEquals("jellyfish", make_data_files.prettyprocess(line))

    def test_Unrecognized_gives_binary_and_unknown(self):
        line = "20141008 00:03:55 14423 befulton 20 14226 0 S 11316512 0 9153040 11314672 92 1732 0 17 3.14 19.75" \
               " 457 00:22.89 0 0 26047 0 3257 0 0 10 1129 " \
               "/N/dc2/scratch/befulton/TrinityMason/trinityrnaseq_r20140717/trinity-plugins/new_app/bin/new_app " \
               "count -t 32 -m 25 -s 1380200998 both.fa"
        self.assertEquals("new_app_unknown", make_data_files.prettyprocess(line))

    def test_Unknown_apps_in_Chrysalis_directory_are_not_Chrysalis(self):
        line = "20141008 00:03:55 14423 befulton 20 14226 0 S 11316512 0 9153040 11314672 92 1732 0 17 3.14 19.75" \
               " 457 00:22.89 0 0 26047 0 3257 0 0 10 1129 " \
               "/N/dc2/scratch/befulton/TrinityMason/trinityrnaseq_r20140717/Chrysalis/new_app " \
               "count -t 32 -m 25 -s 1380200998 both.fa"
        self.assertEquals("new_app_unknown", make_data_files.prettyprocess(line))



