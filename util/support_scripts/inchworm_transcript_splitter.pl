#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long qw(:config no_ignore_case bundling);

use FindBin;

use Cwd;

$ENV{LC_ALL} = 'C';

my $util_dir = "$FindBin::RealBin";


my $usage = <<_EOUSAGE_;


#########################################################################################
#
# Required:
#
#  --iworm <string>                  inchworm assembled contigs
#
#  --left  <string>                  left fragment file
#  --right <string>                  right fragment file
# 
#      or  --single_but_really_paired <string>      single read file containing both pairs.    
#
#  --seqType <string>                fq|fa
#
# Optional (if strand-specific RNA-Seq):
#  
#  --work_dir <string>               directory to perform data processing (default: workdir.\$pid
#
#  --SS_lib_type <string>            RF or FR
#
#  --CPU <int>                       default: 2
#
###########################################################################################

_EOUSAGE_

;


my $inchworm_contigs;
my $left_file;
my $right_file;
my $single_file;
my $seqType;
my $SS_lib_type;
my $work_dir;
my $CPU = 2;

&GetOptions( 'iworm=s' => \$inchworm_contigs,
             'left=s' => \$left_file,
             'right=s' => \$right_file,
             'seqType=s' => \$seqType,
             'SS_lib_type=s' => \$SS_lib_type,
             'work_dir=s' => \$work_dir,
             'CPU=i' => \$CPU,
             'single_but_really_paired=s' => \$single_file,
             );


unless ($inchworm_contigs && ($single_file || ($left_file && $right_file)) && $seqType) {
  die $usage;
}

unless ($work_dir) {
    $work_dir = cwd() . "/jaccard_clip_workdir";
    unless (-d $work_dir) {
        mkdir($work_dir) or die "Error, cannot mkdir $work_dir";
    }
}


my $SYMLINK = ($ENV{NO_SYMLINK}) ? "cp" : "ln -sf";

main: {
	
  my $curr_dir = cwd();
	
  # create full paths to inputs, if not set.
  foreach my $file ($inchworm_contigs, $left_file, $right_file, $single_file) {
      
      unless ($file) { next; }
      
      unless ($file =~ /^\//) {
          $file = "$curr_dir/$file";
      }
  }
  
  my $outdir = $work_dir;
  unless (-d $outdir) {
      mkdir ($outdir) or die "Error, cannot mkdir $outdir";
  }
  
    

  my $target_iworm_fa = "$outdir/iworm.fa";
  &process_cmd("$SYMLINK $inchworm_contigs $target_iworm_fa") unless (-e $target_iworm_fa);
  	
   
  ## run the bowtie alignment pipeline

  my $bowtie_out = "$outdir";
  my $cmd  = "";

  if ($left_file && $right_file) {
          
      $cmd = "$util_dir/bowtie2_wrapper.pl --seqType $seqType --left $left_file --right $right_file --CPU $CPU --target $target_iworm_fa -o $bowtie_out/bowtie2";
  }
  else {
            
      $cmd = "$util_dir/bowtie2_wrapper.pl --seqType $seqType --single $single_file --CPU $CPU --target $target_iworm_fa -o $bowtie_out/bowtie2";
  }
  
  if ($SS_lib_type) {
    $cmd .= " --SS_lib_type $SS_lib_type ";
  }


  chdir $outdir or die "Error, cannot cd to $outdir";
  
  
  &process_cmd($cmd);
    
	
  my $final_bam_file = ($SS_lib_type) ? "$bowtie_out/bowtie2.coordSorted.bam.+.bam" : "$bowtie_out/bowtie2.coordSorted.bam";
	
  my $alignment_file = "bowtie_alignments.for_jaccard.bam";
  &process_cmd("$SYMLINK $final_bam_file $alignment_file");
  
  ## run Jaccard computation:
  my $jaccard_wig_file = "$alignment_file.J100.wig";
  my $frag_coords_file = "$alignment_file.frag_coords";
  &process_cmd("$util_dir/SAM_ordered_pair_jaccard.pl --sam $alignment_file -W 100 > $jaccard_wig_file") unless (-s $jaccard_wig_file && -s $frag_coords_file);
    
  # The above creates a .frag_coords file.  Use this to compute fragment-level coverage
  my $frag_coverage_file = "$frag_coords_file.wig";
  &process_cmd("$util_dir/fragment_coverage_writer.pl $frag_coords_file > $frag_coverage_file") unless (-e $frag_coverage_file);
    
    
  ## define the transcript clip points:
  my $clips_file = "$jaccard_wig_file.clips";
  &process_cmd("$util_dir/jaccard_wig_clipper.pl --jaccard_wig $jaccard_wig_file --coverage_wig $frag_coverage_file  > $clips_file") unless (-e $clips_file) ;

  ## clip the inchworm transcripts:
  &process_cmd("$util_dir/jaccard_fasta_clipper.pl $inchworm_contigs  $clips_file > $inchworm_contigs.clipped.fa") unless (-e "$inchworm_contigs.clipped.fa");
	
	
  exit(0);
	
}


####
sub process_cmd {
  my ($cmd) = @_;
	
  print "CMD: $cmd\n";

  my $ret = system($cmd);
	
  if ($ret) {
    die "Error, cmd: $cmd died with ret $ret";
  }

  return($ret);
}
