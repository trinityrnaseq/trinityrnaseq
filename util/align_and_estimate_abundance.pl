#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;

use Cwd;
use Carp;

use Getopt::Long qw(:config no_ignore_case bundling pass_through);
use Data::Dumper;

my %aligner_params = ( 

    
    ############
    ## Bowtie-1
    ############
    
    
    'bowtie_RSEM' => '--all --best --strata -m 300 --chunkmbs 512',
    # params used by RSEM itself:
    #  -a -m 200
    
    
    'bowtie_eXpress' => '--all --best --strata -m 300 --chunkmbs 512',
    # bowtie -aS -X 800 --offrate 1  (requires: bowtie-build --offrate 1)
    
    
    #############
    ## Bowtie-2
    #############
    
    'bowtie2_RSEM' => '--no-mixed --no-discordant --gbar 1000 --end-to-end -k 200 ',
    
    ## params used by RSEM itself:
    #   --dpad 0 --gbar 99999999 --mp 1,1 --np 1 --score-min L,0,-0.1 -I 1 -X 1000 --no-mixed --no-discordant -k 200
    
        
    'bowtie2_eXpress' => '--no-mixed --no-discordant --gbar 1000 --end-to-end -k 200 ',
    
    
    # recommended eXpress params: http://bio.math.berkeley.edu/eXpress/faq.html
    # -a -X 600 --rdg 6,5 --rfg 6,5 --score-min L,-.6,-.4 --no-discordant --no-mixed
    
        
    'bowtie_none' => '--all --best --strata -m 300 --chunkmbs 512',
    
    'bowtie2_none' => '--no-mixed --no-discordant --gbar 1000 --end-to-end -k 200 ', 
        
    );

my $rsem_add_opts = "";

my $kallisto_add_opts = "";
my $salmon_add_opts= "";
my $salmon_idx_type = 'quasi';
my $salmon_quasi_kmer_length = 31;
my $salmon_fmd_kmer_length = 19;

my $usage = <<__EOUSAGE__;

#########################################################################
#
########################
#  Essential parameters:
########################
#
#  --transcripts <string>           transcript fasta file
#
#  --seqType <string>               fq|fa
# 
#  If Paired-end:
#
#     --left <string>
#     --right <string>
#  
#   or Single-end:
#
#      --single <string>
#
#   or (preferred):
#
#      --samples_file <string>    tab-delimited text file indicating biological replicate relationships.
#                                   ex.
#                                        cond_A    cond_A_rep1    A_rep1_left.fq    A_rep1_right.fq
#                                        cond_A    cond_A_rep2    A_rep2_left.fq    A_rep2_right.fq
#                                        cond_B    cond_B_rep1    B_rep1_left.fq    B_rep1_right.fq
#                                        cond_B    cond_B_rep2    B_rep2_left.fq    B_rep2_right.fq
#
#                      # if single-end instead of paired-end, then leave the 4th column above empty.
#
#
#
#  --est_method <string>           abundance estimation method.
#                                        alignment_based:  RSEM
#                                        alignment_free: kallisto|salmon
#  
###################################
#  Potentially optional parameters:
###################################
#
# --output_dir <string>            write all files to output directory 
#                                  (note, if using --samples_file, output_dir will be set automatically according to replicate name))
#  
#
#  if alignment_based est_method:
#       --aln_method <string>            bowtie|bowtie2 alignment method.  (note: RSEM requires either bowtie or bowtie2)
#                                       
###########
# Optional:
# #########
#
# --SS_lib_type <string>           strand-specific library type:  paired('RF' or 'FR'), single('F' or 'R').
#                                  
# --samples_idx <int>               restricte processing to sample entry (index starts at one)
#
#
# --thread_count                   number of threads to use (default = 4)
#
# --debug                          retain intermediate files
#
#  --gene_trans_map <string>        file containing 'gene(tab)transcript' identifiers per line.
#     or  
#  --trinity_mode                   Setting --trinity_mode will automatically generate the gene_trans_map and use it.
#
#
#  --prep_reference                 prep reference (builds target index)
#
#
########################################
#
#  Parameters for single-end reads:
#
#  --fragment_length <int>         specify RNA-Seq fragment length (default: 200) 
#  --fragment_std <int>            fragment length standard deviation (defalt: 80)
#
########################################
#  
#   bowtie-related parameters: (note, tool-specific settings are further below)
#
#  --max_ins_size <int>             maximum insert size (bowtie -X parameter, default: 800)
#  --coordsort_bam                  provide coord-sorted bam in addition to the default (unsorted) bam.
#
########################################
#  RSEM opts:
#
#  --bowtie_RSEM <string>          if using 'bowtie', default: \"$aligner_params{bowtie_RSEM}\"
#  --bowtie2_RSEM <string>         if using 'bowtie2', default: \"$aligner_params{bowtie2_RSEM}\"
#                                ** if you change the defaults, specify the full set of parameters to use! **
#
#  --include_rsem_bam              provide the RSEM enhanced bam file including posterior probabilities of read assignments.
#  --rsem_add_opts <string>        additional parameters to pass on to rsem-calculate-expression
#
##########################################################################
#  kallisto opts:
#
#  --kallisto_add_opts <string>  default: $kallisto_add_opts  
#
##########################################################################
#
#  salmon opts:
#
#  --salmon_idx_type <string>    quasi|fmd (defalt: $salmon_idx_type)
#  --salmon_add_opts <string>    default: $salmon_add_opts
#
#
#  Example usage
#
#   ## Just prepare the reference for alignment and abundance estimation
#
#    $0 --transcripts Trinity.fasta --est_method RSEM --aln_method bowtie --trinity_mode --prep_reference
#
#   ## Run the alignment and abundance estimation (assumes reference has already been prepped, errors-out if prepped reference not located.)
#
#    $0 --transcripts Trinity.fasta --seqType fq --left reads_1.fq --right reads_2.fq --est_method RSEM --aln_method bowtie --trinity_mode --output_dir rsem_outdir
#
##  ## prep the reference and run the alignment/estimation
#
#    $0 --transcripts Trinity.fasta --seqType fq --left reads_1.fq --right reads_2.fq --est_method RSEM --aln_method bowtie --trinity_mode --prep_reference --output_dir rsem_outdir
#
#   ## Use a samples.txt file:
#
#    $0 --transcripts Trinity.fasta --est_method RSEM --aln_method bowtie2 --prep_reference --trinity_mode --samples_file samples.txt --seqType fq  
#
#########################################################################


__EOUSAGE__

    ;




my $output_dir;
my $help_flag;
my $transcripts;
my $bam_file;
my $DEBUG_flag = 0;
my $SS_lib_type;
my $thread_count = 4;
my $seqType;
my $left;
my $right;
my $single;
my $gene_trans_map_file;
my $max_ins_size = 800;

my $est_method;
my $aln_method = "";

my $retain_sorted_bam_file = 0;

my $fragment_length = 200;
my $fragment_std = 80;

my $output_prefix = "";

# devel opts
my $prep_reference = 0;

my $trinity_mode;

my $include_rsem_bam;
my $coordsort_bam_flag = 0;

my $samples_file = "";
my $samples_idx = 0;

&GetOptions ( 'help|h' => \$help_flag,
              'transcripts=s' => \$transcripts,
              'name_sorted_bam=s' => \$bam_file,
              'debug' => \$DEBUG_flag,
              'SS_lib_type=s' => \$SS_lib_type,

              'thread_count=i' => \$thread_count,
              
              'gene_trans_map=s' => \$gene_trans_map_file,
              'trinity_mode' => \$trinity_mode,
              
              'seqType=s' => \$seqType,
              'left=s' => \$left,
              'right=s' => \$right,
              'single=s' => \$single,
              'max_ins_size=i' => \$max_ins_size,
              'samples_file=s' => \$samples_file,
              'samples_idx=i' => \$samples_idx,
              
              'output_dir=s' => \$output_dir,
      
              'est_method=s' => \$est_method,
              'aln_method=s' => \$aln_method,


              'include_rsem_bam' => \$include_rsem_bam,

              #'output_prefix=s' => \$output_prefix,
              
              ##  devel opts
              'prep_reference' => \$prep_reference,

              # opts for single-end reads
              'fragment_length=i' => \$fragment_length,
              'fragment_std=i' => \$fragment_std,
              
              #
              'bowtie_RSEM=s' => \($aligner_params{'bowtie_RSEM'}),
              'bowtie2_RSEM=s' => \($aligner_params{'bowtie2_RSEM'}),

              
              'rsem_add_opts=s' => \$rsem_add_opts,
              'kallisto_add_opts=s' => \$kallisto_add_opts,
              'salmon_add_opts=s' => \$salmon_add_opts,
    
              'coordsort_bam' => \$coordsort_bam_flag,

             'salmon_idx_type=s' => \$salmon_idx_type,
             'salmon_quasi_kmer_length=i' => \$salmon_quasi_kmer_length,
             'salmon_fmd_kmer_length=i' => \$salmon_fmd_kmer_length,
    
    );



if (@ARGV) {
    die "Error, don't understand arguments: @ARGV ";
}

if ($help_flag) {
    die $usage;
}

unless ($est_method) {
    die $usage;
}

my @EST_METHODS = qw(RSEM kallisto salmon);
my %ALIGNMENT_BASED_EST_METHODS = map { + $_ => 1 } qw (RSEM);
my %ALIGNMENT_FREE_EST_METHODS = map { + $_ => 1 } qw (kallisto salmon);


unless (
    
    ($est_method && $prep_reference && $transcripts && (! ($single||$left||$right||$samples_file)) ) ## just prep reference
    
    || 
    
    ($transcripts && $est_method && $seqType && ($single || ($left && $right) || $samples_file)) # do alignment
    
    ) {
    
    die "Error, missing parameter. See example usage options below.\n" . $usage;
}


if  ($ALIGNMENT_FREE_EST_METHODS{$est_method}) {
    $aln_method = "none";
}
elsif ($aln_method !~ /bowtie2?/) {
    die "Error, --aln_method must be either 'bowtie' or 'bowtie2' ";
}


unless ($est_method =~ /^(RSEM|kallisto|salmon|none)$/i) {
    die "Error, --est_method @EST_METHODS only\n";
}


my @samples_to_process;
if ($samples_file) {
    @samples_to_process = &parse_samples_file($samples_file);
    if ($samples_idx >= 0) {
        my $num_samples = scalar(@samples_to_process);
        if ($samples_idx > $num_samples) {
            die "Error, sample index $samples_idx > $num_samples num samples ";
        }
    }
}
elsif ( ($left && $right) || $single) {
    
    unless ($output_dir) {
        die "Error, must specify output directory name via: --output_dir   ";
    }
    @samples_to_process = &create_sample_definition($output_dir, $left, $right, $single);
    
}



my $PE_mode = 1;

if ($single || (@samples_to_process && $samples_to_process[0]->{single})) {
    
    unless ($fragment_length) {
        die "Error, specify --fragment_length for single-end reads (note, not the length of the read but the mean fragment length)\n\n";
    }

    $PE_mode = 0;
}

    
$transcripts = &create_full_path($transcripts);

$gene_trans_map_file = &create_full_path($gene_trans_map_file) if $gene_trans_map_file;
if ($gene_trans_map_file && ! -s $gene_trans_map_file) {
    die "Error, $gene_trans_map_file doesn't exist or is empty";
}


if ($SS_lib_type) {
    unless ($SS_lib_type =~ /^(RF|FR|R|F)$/) {
        die "Error, do not recognize SS_lib_type: [$SS_lib_type]\n";
    }
    if ($PE_mode && length($SS_lib_type) != 2 ) {
        die "Error, SS_lib_type [$SS_lib_type] is not compatible with paired reads";
    }
}

if ( $thread_count !~ /^\d+$/ ) {
    die "Error, --thread_count value must be an integer";
}


{  # check for required tools in PATH 
    
    my $missing = 0;
    my @tools = ('samtools');
    if ($aln_method eq 'bowtie') {
        push (@tools, 'bowtie-build', 'bowtie');
    }
    elsif ($aln_method eq 'bowtie2') {
        push (@tools, 'bowtie2', 'bowtie2-build');
    }
    
    if ($est_method =~ /^RSEM$/i) {
        push (@tools, 'rsem-calculate-expression');
    }
    elsif ($est_method eq 'kallisto') {
        push (@tools, 'kallisto');
    }
    elsif ($est_method eq 'salmon') {
        push (@tools, 'salmon');
    }
    
        
    foreach my $tool (@tools) {
        my $p = `which $tool`;
        unless ($p =~ /\w/) {
            warn("ERROR, cannot find $tool in PATH setting: $ENV{PATH}\n\n");
            $missing = 1;
        }
    }
    if ($missing) {
        die "Please be sure the utilities @tools are available via your PATH setting.\n";
    }
}



main: {

    if ($trinity_mode && ! $gene_trans_map_file) {
        $gene_trans_map_file = "$transcripts.gene_trans_map";
        my $cmd = "$FindBin::RealBin/support_scripts/get_Trinity_gene_to_trans_map.pl $transcripts > $gene_trans_map_file";
        &process_cmd($cmd) unless (-e $gene_trans_map_file);
    }


    
    if ($ALIGNMENT_BASED_EST_METHODS{$est_method}) {
        
        &run_alignment_BASED_estimation(@samples_to_process);

    }
    else {
        &run_alignment_FREE_estimation(@samples_to_process);
    }

    exit(0);
}



####
sub run_alignment_FREE_estimation {
    my @samples = @_;

    
    if ($est_method eq "kallisto") {
        &run_kallisto(@samples);
    }
    elsif ($est_method eq "salmon") {
        &run_salmon(@samples);
    }
    else {
        die "Error, not recognizing est_method: $est_method";
        # sholdn't get here
    }
}



####
sub run_alignment_BASED_estimation {
    my @samples = @_;

    
    my $db_index_name = "$transcripts.${aln_method}";
    

    ###############################################
    ## Prepare transcript database for alignments
    ###############################################
    
    
    if ($prep_reference) {
        
        my $cmd = "${aln_method}-build $transcripts $db_index_name";
        
        unless (-e "$db_index_name.ok") { 
            
            if (-e "$db_index_name.started") {
                print STDERR "WARNING - looks like the prep for $db_index_name was already started by another process.  Proceeding with caution.\n";
            }
            
            &process_cmd("touch $db_index_name.started");
            
            &process_cmd($cmd);
            
            rename("$db_index_name.started", "$db_index_name.ok");
            
        }
        
        
    }

    if (! -e "$db_index_name.ok") {
        die "Error, index $db_index_name not prepared.  Be sure to include parameter '--prep_reference' to first prepare the reference for alignment.";
    }
        
    
    
    my $rsem_prefix = &create_full_path("$transcripts.RSEM");
    
    if ($est_method eq 'RSEM') {
                
        if ($prep_reference) {
            
            if (-e "$rsem_prefix.rsem.prepped.started") {
                print STDERR "WARNING - appears that another process has started the rsem-prep step... proceeding with caution.\n";
            }

            unless (-e "$rsem_prefix.rsem.prepped.ok") {
                
                &process_cmd("touch $rsem_prefix.rsem.prepped.started");
                
                my $cmd = "rsem-prepare-reference "; #--no-bowtie"; # update for RSEM-2.15
                
                if ($gene_trans_map_file) {
                    $cmd .= " --transcript-to-gene-map $gene_trans_map_file";
                }
                $cmd .= " $transcripts $rsem_prefix";
                
                &process_cmd($cmd);
                
                rename("$rsem_prefix.rsem.prepped.started", "$rsem_prefix.rsem.prepped.ok");
            }


            unless (-e "$rsem_prefix.rsem.prepped.ok") {
                
                die "Error, the RSEM data must first be prepped. Please rerun with '--prep_reference' parameter.\n"; 
                
            }
        }
                
    }
    
    
    unless (@samples) {
        print STDERR "Only prepping reference. Stopping now.\n";
        exit(0);
    }
    
    print STDERR Dumper(\@samples);
    
    my $curr_workdir = cwd();
    foreach my $sample_href (@samples) {
        chdir $curr_workdir or die "Error, cannot cd to $curr_workdir";
        # process below will cd into output dir
        &run_alignment_do_quant($sample_href, $db_index_name, $rsem_prefix);
    }

}

####
sub run_alignment_do_quant {
    my ($sample_href, $db_index_name, $rsem_prefix) = @_;

    my $output_dir = $sample_href->{output_dir};
    
    #####################
    ## Run alignments
    #####################
    
    unless (-d $output_dir) {
        system("mkdir -p $output_dir");
    }
    chdir $output_dir or die "Error, cannot cd to output directory $output_dir";
    
    my $prefix = $output_prefix;
    if ($prefix) {
        $prefix .= "."; # add separator in filename
    }
    my $bam_file = "${prefix}${aln_method}.bam";
    my $bam_file_ok = "$bam_file.ok";
        

    my $read_type = ($seqType eq "fq") ? "-q" : "-f";
        
    ##############
    ## Align reads
    
    my $bowtie_cmd;
    
    if ($aln_method eq 'bowtie') {
        if ($PE_mode) {
            my ($left_file, $right_file) = ($sample_href->{left}, $sample_href->{right});
            ## PE alignment
            $bowtie_cmd = "set -o pipefail && bowtie $read_type " . $aligner_params{"${aln_method}_${est_method}"} . " -X $max_ins_size -S -p $thread_count $db_index_name -1 $left_file -2 $right_file | samtools view -@ $thread_count -F 4 -S -b | samtools sort -@ $thread_count -n -o $bam_file ";
            
        }
        else {
            my $single_file = $sample_href->{single};
            # SE alignment
            $bowtie_cmd = "set -o pipefail && bowtie $read_type " . $aligner_params{"${aln_method}_${est_method}"} . " -S -p $thread_count $db_index_name $single_file | samtools view -@ $thread_count -F 4 -S -b | samtools sort -@ $thread_count -n -o $bam_file ";
        }
    }
    elsif ($aln_method eq 'bowtie2') {
        
        if ($PE_mode) {
            ## PE alignment
            my ($left_file, $right_file) = ($sample_href->{left}, $sample_href->{right});
            $bowtie_cmd = "set -o pipefail && bowtie2 " . $aligner_params{"${aln_method}_${est_method}"} . " $read_type -X $max_ins_size -x $db_index_name -1 $left_file -2 $right_file -p $thread_count | samtools view -@ $thread_count -F 4 -S -b | samtools sort -@ $thread_count -n -o $bam_file ";
        }
        else {
            # SE alignment
            my $single_file = $sample_href->{single};
            $bowtie_cmd = "set -o pipefail && bowtie2 " . $aligner_params{"${aln_method}_${est_method}"} . " $read_type -x $db_index_name -U $single_file -p $thread_count | samtools view -@ $thread_count -F 4 -S -b | samtools sort -@ $thread_count -n -o $bam_file ";
        }
    }
    
    &process_cmd($bowtie_cmd) unless (-s $bam_file && -e $bam_file_ok);
    
    &process_cmd("touch $bam_file_ok") unless (-e $bam_file_ok);
    
        
    if ($est_method eq "RSEM") {
        
        # convert bam file for use with rsem:
        &process_cmd("convert-sam-for-rsem $bam_file $bam_file.for_rsem");
        
        &run_RSEM("$bam_file.for_rsem.bam", $rsem_prefix, $output_prefix);
    }
    elsif ($est_method eq "none") {
        print STDERR "Not running abundance estimation, stopping now after alignment.\n";
    }
    else {
        die "Error, --est_method $est_method is not supported";
    }
    
    if ($coordsort_bam_flag) {
        
        &sort_bam_file($bam_file);
        
    }
    
    return;
    
}


####
sub sort_bam_file {
    my ($bam_file) = @_;
    my $sorted_bam_file = $bam_file;
    $sorted_bam_file =~ s/bam$/csorted/;
    if (! -e "$sorted_bam_file.bam.ok") {
        ## sort the bam file
        
        my $cmd = "samtools sort $bam_file -o $sorted_bam_file.bam";
        &process_cmd($cmd);
        $cmd = "samtools index $sorted_bam_file.bam";
        &process_cmd($cmd);
        
        &process_cmd("touch $sorted_bam_file.bam.ok");
    }

    return;
}


####
sub run_RSEM {
    my ($bam_file, $rsem_prefix, $output_prefix) = @_;
        

    unless ($output_prefix) {
        $output_prefix = "RSEM";
    }
    
    my $keep_intermediate_files_opt = ($DEBUG_flag) ? "--keep-intermediate-files" : "";
    
    my $fraglength_info_txt = "";
    if ($single) {
        $fraglength_info_txt = "--fragment-length-mean $fragment_length --fragment-length-sd $fragment_std";
    }
    
    my $SS_opt = "";
    if ($SS_lib_type) {
        if ($SS_lib_type =~ /^F/) {
            $SS_opt = "--forward-prob 1.0";
        }
        else {
            $SS_opt = "--forward-prob 0";
        }
    }
    
    my $no_qualities_string = "";
    if ($seqType eq 'fa') {
        $no_qualities_string = "--no-qualities";
    }

    my $paired_flag_text = ($PE_mode) ? "--paired-end" : "";

    my $rsem_bam_flag = ($include_rsem_bam) ? "" : "--no-bam-output";


    my $cmd = "rsem-calculate-expression $no_qualities_string "
        . "$paired_flag_text "
        . " $rsem_add_opts "
        . "-p $thread_count "
        . "$fraglength_info_txt "
        . "$keep_intermediate_files_opt "
        . "$SS_opt $rsem_bam_flag "
        . "--bam $bam_file "
        . "$rsem_prefix "
        . "$output_prefix ";
    
    unless (-e "$output_prefix.isoforms.results.ok") {
        &process_cmd($cmd);
    }
    &process_cmd("touch $output_prefix.isoforms.results.ok");

    return;
}


####
sub process_cmd {
    my ($cmd) = @_;

    unless ($cmd) {
        confess "Error, no cmd specified";
    }
    
    print STDERR "CMD: $cmd\n";
    
    my $ret = system("bash", "-o", "pipefail", "-c", $cmd);
    
    if ($ret) {
        die "Error, cmd: $cmd died with ret: $ret";
    }
    
    return;
}

###
sub create_full_path {
    my ($file_list) = shift;
    
    my $cwd = cwd();

    my @files;

    foreach my $file (split(/,/, $file_list)) {
        
        
        if ($file !~ m|^/|) { # must be a relative path
            $file = $cwd . "/$file";
        }
        
        push (@files, $file);
    }

    $file_list = join(",", @files);

    return($file_list);

    
}

####
sub add_zcat_gz {
    my ($file_listing) = @_;

    my @files;

    foreach my $file (split(/,/, $file_listing)) {
        
        if ($file =~ /\.gz$/) {

            $file = "<(gunzip -c $file)";  # used to be zcat
    
            
        }
        push (@files, $file);
    }

    $file_listing = join(",", @files);

    return($file_listing);
}


####
sub run_kallisto {
    my @samples = @_;
    
    my $kallisto_index = "$transcripts.kallisto_idx";
    
    if ( (! $prep_reference) && (! -e $kallisto_index)) {
        confess "Error, no kallisto index file: $kallisto_index, and --prep_reference not set.  Re-run with --prep_reference";
    }
    if ($prep_reference && ! -e $kallisto_index) {
        
        my $cmd = "kallisto index -i $kallisto_index $transcripts";
        &process_cmd($cmd);
    }


    if ($SS_lib_type) {
        # add strand-specific options for kallisto
        my $kallisto_ss_opt = ($SS_lib_type =~ /^R/) ? "--rf-stranded" : "--fr-stranded";
        if ($kallisto_add_opts !~ /$kallisto_ss_opt/) {
            $kallisto_add_opts .= " $kallisto_add_opts";
        }
    }
        
    foreach my $sample_href (@samples) {
     
        my ($output_dir, $left_file, $right_file, $single_file) = ($sample_href->{output_dir},
                                                                   $sample_href->{left},
                                                                   $sample_href->{right},
                                                                   $sample_href->{single});
   
        if ($left_file && $right_file) {
            
            my $cmd = "kallisto quant -i $kallisto_index $kallisto_add_opts -o $output_dir $left_file $right_file";
            &process_cmd($cmd);
        }
        elsif ($single_file) {
            my $cmd = "kallisto quant -l $fragment_length -s $fragment_std -i $kallisto_index -o $output_dir $kallisto_add_opts --single $single_file";
            &process_cmd($cmd);
        }
        
        
        if ($gene_trans_map_file) {
            
            my $cmd = "$FindBin::RealBin/support_scripts/kallisto_trans_to_gene_results.pl $output_dir/abundance.tsv $gene_trans_map_file > $output_dir/abundance.tsv.genes";
            &process_cmd($cmd);
        }
    }

    return;
}



####
sub run_salmon {
    my (@samples) = @_;
    
    my $salmon_index = "$transcripts.salmon_${salmon_idx_type}.idx";
    
    if ( (! $prep_reference) && (! -e $salmon_index)) {
        confess "Error, no salmon index file: $salmon_index, and --prep_reference not set.  Re-run with --prep_reference";
    }
    if ($prep_reference && ! -e $salmon_index) {
        
        ## Prep salmon index
        my $cmd;
        
        if ($salmon_idx_type eq 'quasi') {
            $cmd = "salmon index -t $transcripts --keepDuplicates -i $salmon_index --type quasi -k $salmon_quasi_kmer_length -p $thread_count";
        }
        elsif ($salmon_idx_type eq 'fmd') {
            $cmd = "salmon index -t $transcripts --keepDuplicates -i $salmon_index --type fmd -p $thread_count";
        }
        else {
            die "Error, not recognizing idx type: $salmon_idx_type";
        }
        
        &process_cmd($cmd);
    }

    my $num_failures = 0;
    
    foreach my $sample_href (@samples) {
        
        my ($output_dir, $left_file, $right_file, $single_file) = ($sample_href->{output_dir},
                                                                   $sample_href->{left},
                                                                   $sample_href->{right},
                                                                   $sample_href->{single});
        

        
        my $outdir = $output_dir; #"$output_dir.$salmon_idx_type";

        if (-s "$outdir/quant.sf") {
            print STDERR "-output already exists: $outdir/quant.sf, skipping.\n";
            next;
        }
        

        eval {
        
            if ($left_file && $right_file) {
                ## PE mode
                my $cmd;
                my $libtype = ($SS_lib_type) ? "IS" . substr($SS_lib_type, 0, 1) : "IU";
                
                if ($salmon_idx_type eq 'quasi') {
                    $cmd = "salmon quant -i $salmon_index -l $libtype -1 $left_file -2 $right_file -o $outdir $salmon_add_opts -p $thread_count";
                }
                elsif ($salmon_idx_type eq 'fmd') {
                    $cmd = "salmon quant -i $salmon_index -l $libtype -1 $left_file -2 $right_file -k $salmon_fmd_kmer_length -o $outdir $salmon_add_opts -p $thread_count";
                }
                else {
                    die "Error, not recognizing salmon_idx_type: $salmon_idx_type";
                }
                
                &process_cmd($cmd);
            
            }
            elsif ($single_file) {
                my $libtype = ($SS_lib_type) ? "S" . substr($SS_lib_type, 0, 1) : "U";
                my $cmd;
                
                if ($salmon_idx_type eq 'quasi') {
                    $cmd = "salmon quant -i $salmon_index -l $libtype -r $single_file -o $outdir $salmon_add_opts -p $thread_count";
                }
                elsif ($salmon_idx_type eq 'fmd') {
                    $cmd = "salmon quant -i $salmon_index -l $libtype -r $single_file -k $salmon_fmd_kmer_length -o $outdir $salmon_add_opts -p $thread_count";
                }
                else {
                    die "Error, not recognizing salmon_idx_type: $salmon_idx_type";
                }
                
                &process_cmd($cmd);
                
            }
            
            if ($gene_trans_map_file) {
                
                my $cmd = "$FindBin::RealBin/support_scripts/salmon_trans_to_gene_results.pl $output_dir/quant.sf $gene_trans_map_file > $output_dir/quant.sf.genes";
                &process_cmd($cmd);
            }
        };
        if ($@) {
            $num_failures++;
            print STDERR "Error detected: $@";
        }
    }


    if ($num_failures) {
        die "Error, encountered $num_failures failed salmon jobs. See errors above";
    }
    
    return;
}



####
sub parse_samples_file {
    my ($samples_file) = @_; 

    my @samples_to_process;

    my %seen;
    open (my $fh, $samples_file) or die "Error, cannot open file: [$samples_file]";
    while (<$fh>) {
        chomp;
        if (/^\#/) { next; }
        unless (/\w/) { next; }
        if (/^\-/) { next; } 
        s/^\s+|\s+$//g; # trim trailing ws
        my @x = split(/\s+/);

        my $sample_name = $x[0];
        my $rep_name = $x[1];
        if ($seen{$rep_name}) {
            die "Error, replicate names must be unique.  Found $rep_name listed multiple times";
        }
        $seen{$rep_name}++;
        
        my $output_dir = $rep_name;
        
        my $left_fq = $x[2];
        my $right_fq = $x[3];
        
        if ($left_fq) {
            unless (-s $left_fq) {
                die "Error, cannot locate file: $left_fq as specified in samples file: $samples_file";
            }
            $left_fq = &create_full_path($left_fq);
            if ($left_fq =~ /\.gz$/) {
                $left_fq = &add_zcat_gz($left_fq) if ($aln_method eq "bowtie");
            }
        }
        else {
            die "Error, cannot parse line $_ of samples file: $samples_file . See usage info for samples file formatting requirements.";
            
        }
        if ($right_fq) {
            unless (-s $right_fq) {
                die "Error, cannot locate file $right_fq as specified in samples file: $samples_file";
            }
            $right_fq = &create_full_path($right_fq);
            if ($right_fq =~ /\.gz$/) {
                $right_fq = &add_zcat_gz($right_fq) if ($aln_method eq "bowtie");
            } 
        }

        if ($left_fq && $right_fq) {

            push (@samples_to_process, { left => $left_fq,
                                         right => $right_fq,
                                         output_dir => $output_dir,
                  } );
        }
        else {
            push (@samples_to_process, { single => $left_fq,
                                         output_dir => $output_dir,
                  } );
        }
        
    }
   
    
    return (@samples_to_process);
}


####
sub create_sample_definition {
    my ($output_dir, $left, $right, $single) = @_;

    $left = &create_full_path($left) if $left;
    $right = &create_full_path($right) if $right;
    $single = &create_full_path($single) if $single;

    if ($left && $left =~ /\.gz$/) {
        $left = &add_zcat_gz($left) if ($aln_method eq "bowtie");
    }
    if ($right && $right =~ /\.gz$/) {
        $right = &add_zcat_gz($right) if ($aln_method eq "bowtie");
    }
    if ($single && $single =~ /\.gz$/) {
        $single = &add_zcat_gz($single) if ($aln_method eq "bowtie");
    }
    
    
    if ($left && $right) {
        return ( { left => $left,
                   right => $right,
                   output_dir => $output_dir,
                 } );
    }
    else {
        return( { single => $single,
                  output_dir => $output_dir,
                } );
    }

}
