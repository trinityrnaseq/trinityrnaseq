#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;

use Cwd;

use Getopt::Long qw(:config no_ignore_case bundling pass_through);


my $RSEM_DIR = "$FindBin::Bin/../trinity-plugins/rsem";
$ENV{PATH} = "$RSEM_DIR:$ENV{PATH}"; # be sure to use the included rsem package over other ones installed.


my %aligner_params = ( 'bowtie_RSEM' => '--all --best --strata -m 300 --chunkmbs 512',
                                              
                       'bowtie2_RSEM' => '--no-mixed --no-discordant --gbar 1000 --end-to-end',
                       
                       'bowtie_eXpress' => '--all --best --strata -m 300 --chunkmbs 512',
                       
                       'bowtie2_eXpress' => '--end-to-end',
                       
    );

my $rsem_add_opts = "";
my $eXpress_add_opts = "";


my $usage = <<__EOUSAGE__;

#########################################################################
#
#  --transcripts <string>           transcript fasta file
#  --seqType <string>               fq|fa
# 
#  If Paired-end:
#
#  --left <string>
#  --right <string>
#  
#    or Single-end:
#
#  --single <string>
#
#  --est_method <string>            RSEM|eXpress       abundance estimation method.
#
#  --aln_method <string>            bowtie|bowtie2|(path to bam file) alignment method.  (note: RSEM requires bowtie)
#                                       (if you already have a bam file, you can use it here instead of rerunning bowtie)
#
# Optional:
# 
# --SS_lib_type <string>           strand-specific library type:  paired('RF' or 'FR'), single('F' or 'R').
#
# --thread_count                   number of threads to use (default = 4)
#
# --max_ins_size <int>             maximum insert size (bowtie -X parameter, default: 800)
#
# --debug                          retain intermediate files
#
# --output_dir <string>            write all files to output directory
#  
#
#  --gene_trans_map <string>        file containing 'gene(tab)transcript' identifiers per line.
#     or  
#  --trinity_mode                   Setting --trinity_mode will automatically generate the gene_trans_map and use it.
#
#
#  --prep_reference                 prep reference set for eXpress (builds bowtie index, etc)
#
#  --output_prefix <string>         prefix for output files.  Defaults to --est_method setting.
#
#
#  --coordsort_bam                  provide coord-sorted bam in addition to the default (unsorted) bam.
#
#  --show_full_usage_info           provide more detailed usage info for customizing the alignment or abundance estimation parameters.
#
#############################
#  RSEM opts:
#  --fragment_length <int>         optionally specify fragment length (not seq length, but frag size ie. 300) for SE reads.
#
#  --include_rsem_bam              provide the RSEM enhanced bam file including posterior probabilities of read assignments.
#
#########################################################################
#
#  Example usage:
#
#   ## Just prepare the reference for alignment and abundance estimation
#
#    $0 --transcripts Trinity.fasta --est_method RSEM --aln_method bowtie --trinity_mode --prep_reference
#
#   ## Run the alignment and abundance estimation (assumes reference has already been prepped, errors-out if prepped reference not located.)
#
#    $0 --transcripts Trinity.fasta --seqType fq --left reads_1.fq --right reads_2.fq --est_method RSEM --aln_method bowtie --trinity_mode 
#
##  ## prep the reference and run the alignment/estimation
#
#    $0 --transcripts Trinity.fasta --seqType fq --left reads_1.fq --right reads_2.fq --est_method RSEM --aln_method bowtie --trinity_mode --prep_reference
#
#########################################################################

__EOUSAGE__

    ;



my $advanced_usage_info = <<__EOADVANCEDUSAGE__;


#############################################################
## Customizing alignment and abundance estimation parameters.
#############################################################
#
#  Default alignment parameters are:
#
#  --bowtie_RSEM <string>    default: \"$aligner_params{bowtie_RSEM}\"
#  --bowtie2_RSEM <string>   default: \"$aligner_params{bowtie2_RSEM}\"
#
#  --bowtie_eXpress <string>  default: \"$aligner_params{bowtie_eXpress}\"
#  --bowtie2_eXpress <string> default: \"$aligner_params{bowtie2_eXpress}\"
#
#  Options to pass on to RSEM or eXpress
#
#  --rsem_add_opts <string>    default: "$rsem_add_opts"
#
#  --eXpress_add_opts <string>  default: "$eXpress_add_opts"
#
#  * note, options for handling strand-specific reads are already taken care of internally, so no need to 
#    pass on those parameters.  
#
##############################################################

__EOADVANCEDUSAGE__

    ;

my $show_full_usage_info;

my $output_dir;
my $help_flag;
my $transcripts;
my $bam_file;
my $paired_flag = "";
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
my $aln_method;


my $retain_sorted_bam_file = 0;

my $fragment_length = "";

my $output_prefix = "";

# devel opts
my $prep_reference = 0;

my $trinity_mode;

my $include_rsem_bam;
my $coordsort_bam_flag = 0;

&GetOptions ( 'h' => \$help_flag,
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
              
              'output_dir=s' => \$output_dir,
      
              'est_method=s' => \$est_method,
              'aln_method=s' => \$aln_method,


              'include_rsem_bam' => \$include_rsem_bam,

              'output_prefix=s' => \$output_prefix,
              
              ##  devel opts
              'prep_reference' => \$prep_reference,

              ## rsem opts
              'fragment_length=i' => \$fragment_length,
              
              #
              'show_full_usage_info' => \$show_full_usage_info,
              
   
              'bowtie_RSEM=s' => \($aligner_params{'bowtie_RSEM'}),
              'bowtie2_RSEM=s' => \($aligner_params{'bowtie2_RSEM'}),
              'bowtie_eXpress=s' => \($aligner_params{'bowtie_eXpress'}),
              'bowtie2_eXpress=s' => \($aligner_params{'bowtie2_eXpress'}),
              
              'rsem_add_opts=s' => \$rsem_add_opts,
              'eXpress_add_opts=s' => \$eXpress_add_opts,
                   
              'coordsort_bam' => \$coordsort_bam_flag,
    
              );



if (@ARGV) {
    die "Error, don't understand arguments: @ARGV ";
}

if ($help_flag) {
    die $usage;
}
if ($show_full_usage_info) {
    die "$usage\n\n$advanced_usage_info\n\n";
}

unless (($est_method && $aln_method && $prep_reference && $transcripts) ## just prep reference

        || ($transcripts && $aln_method && $est_method && $seqType && ($single || ($left && $right))) # do alignment

    ) {

    die "Error, missing parameter. See example usage options below.\n" . $usage;
}

my $PROCESSING_EXISTING_BAM_FLAG = 0;
if ($aln_method =~ /\.bam$/) {
    if (-s $aln_method) {
        # no problem, using the bam file
        $PROCESSING_EXISTING_BAM_FLAG = 1;
    }
    else {
        die "Error, file $aln_method does not exist or is empty.  ";
    }
}
elsif ($aln_method !~ /bowtie2?/) {
    die "Error, --aln_method must be either 'bowtie' or 'bowtie2' ";
}

unless ($est_method =~ /^(RSEM|eXpress)$/) {
    die "Error, --est_method  'RSEM' or 'eXpress' only, and capitalization matters. :) \n";
}




if ($output_dir) {
    $left = &create_full_path($left) if $left;
    $right = &create_full_path($right) if $right;
    $single = &create_full_path($single) if $single;
    
    $transcripts = &create_full_path($transcripts);

    if ($gene_trans_map_file) {
        $gene_trans_map_file = &create_full_path($gene_trans_map_file);
    }
}


if ($left && $left =~ /\.gz$/) {
    $left = &add_zcat_gz($left);
}
if ($right && $right =~ /\.gz$/) {
    $right = &add_zcat_gz($right);
}
if ($single && $single =~ /\.gz$/) {
    $single = &add_zcat_gz($single);
}



if ($SS_lib_type) {
    unless ($SS_lib_type =~ /^(RF|FR|R|F)$/) {
        die "Error, do not recognize SS_lib_type: [$SS_lib_type]\n";
    }
    if ($left && $right && length($SS_lib_type) != 2 ) {
        die "Error, SS_lib_type [$SS_lib_type] is not compatible with paired reads";
    }
}

if ( $thread_count !~ /^\d+$/ ) {
    die "Error, --thread_count value must be an integer";
}


{  # check for RSEM installation in PATH 
    
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
    elsif ($est_method =~ /^eXpress$/i) {
        push (@tools, 'express');
    }
    

    foreach my $tool (@tools) {
        my $p = `which $tool`;
        unless ($p =~ /\w/) {
            warn("ERROR, cannot find $tool in PATH setting: $ENV{PATH}\n\n");
            $missing = 1;
        }
    }
    if ($missing) {
        die "Please be sure bowtie and express are installed and the utilities @tools are available via your PATH setting.\n";
    }
}

main: {


    ###############################################
    ## Prepare transcript database for alignments
    ###############################################

    my $db_index_name = "$transcripts.${aln_method}";
    

    if ($trinity_mode && ! $gene_trans_map_file) {
        $gene_trans_map_file = "$transcripts.gene_trans_map";
        my $cmd = "$FindBin::Bin/support_scripts/get_Trinity_gene_to_trans_map.pl $transcripts > $gene_trans_map_file";
        &process_cmd($cmd) unless (-e $gene_trans_map_file);
    }
        
    unless ($PROCESSING_EXISTING_BAM_FLAG) {
            
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
    
    
    unless ( ($left && $right) || $single) {
        print STDERR "Only prepping reference. Stopping now.\n";
        exit(0);
    }


    #####################
    ## Run alignments
    #####################


    if ($output_dir) {
        
        my $init_dir = cwd();
        
        system("mkdir -p $output_dir");
        chdir $output_dir or die "Error, cannot cd to output directory $output_dir";
    
    }

    my $prefix = $output_prefix;
    if ($prefix) {
        $prefix .= "."; # add separator in filename
    }

    my $read_type = ($seqType eq "fq") ? "-q" : "-f";
    
    ## run bowtie
    my $bowtie_cmd;
    my $bam_file = "${prefix}${aln_method}.bam";
    my $bam_file_ok = "$bam_file.ok";
    if ($PROCESSING_EXISTING_BAM_FLAG) {
        $bam_file = $aln_method;
        unless (-e $bam_file_ok) {
            &process_cmd("touch $bam_file_ok");
        }
    }
        
    if ($left && $right) {
        $paired_flag = 1;
    }
    
    if ($aln_method eq 'bowtie') {
        if ($left && $right) {
            $bowtie_cmd = "bowtie $read_type " . $aligner_params{"${aln_method}_${est_method}"} . " -X $max_ins_size -S -p $thread_count $db_index_name -1 $left -2 $right | samtools view -F 4 -S -b -o $bam_file -";
           
        }
        else {
            $bowtie_cmd = "bowtie $read_type " . $aligner_params{"${aln_method}_${est_method}"} . " -S -p $thread_count $db_index_name $single | samtools view -F 4 -S -b -o $bam_file -";
        }
    }
    elsif ($aln_method eq 'bowtie2') {
        
        if ($left && $right) {
                
            $bowtie_cmd = "bowtie2 " . $aligner_params{"${aln_method}_${est_method}"} . " $read_type -X $max_ins_size -x $db_index_name -1 $left -2 $right -p $thread_count | samtools view -F 4 -S -b -o $bam_file -";
        }
        else {

            $bowtie_cmd = "bowtie2 " . $aligner_params{"${aln_method}_${est_method}"} . " $read_type -x $db_index_name -U $single -p $thread_count | samtools view -F 4 -S -b -o $bam_file -";
        }
    }
    
    &process_cmd($bowtie_cmd) unless (-s $bam_file && -e $bam_file_ok);
    
    &process_cmd("touch $bam_file_ok") unless (-e $bam_file_ok);
    
    if ($est_method eq "eXpress") {
        &run_eXpress($bam_file);
    }
    elsif ($est_method eq "RSEM") {
        &run_RSEM($bam_file, $rsem_prefix, $output_prefix);
    }
    else {
        die "Error, --est_method $est_method is not supported";
    }
    
    if ($coordsort_bam_flag) {
        
        my $sorted_bam_file = $bam_file;
        $sorted_bam_file =~ s/bam$/csorted/;
        if (! -e "$sorted_bam_file.bam.ok") {
            ## sort the bam file
            
            my $cmd = "samtools sort $bam_file $sorted_bam_file";
            &process_cmd($cmd);
            $cmd = "samtools index $sorted_bam_file.bam";
            &process_cmd($cmd);
            
            &process_cmd("touch $sorted_bam_file.bam.ok");
        }
    }
    

    
    exit(0);
    
}


####
sub run_eXpress {
    my ($bam_file) = @_;
    
    my $SS_opt = "";
    if ($SS_lib_type) {
        if ($SS_lib_type eq "F") {
            $SS_opt = "--f-stranded";
        }
        elsif ($SS_lib_type eq "R") {
            $SS_opt = "--r-stranded";
        }
        elsif ($SS_lib_type eq "FR") {
            $SS_opt = "--fr-stranded";
        }
        elsif ($SS_lib_type eq "RF") {
            $SS_opt = "--rf-stranded";
        }
    }
    
    ## run eXpress
    my $express_cmd = "express $SS_opt $eXpress_add_opts $transcripts";
        
    my $cmd = "$express_cmd $bam_file";
    &process_cmd($cmd);
    
    if ($gene_trans_map_file) {
        &make_eXpress_genes_file("results.xprs", $gene_trans_map_file);
    }
    
    return;

}

####
sub make_eXpress_genes_file {
    my ($results_xprs_file, $gene_to_trans_map_file) = @_;

    my $cmd = "$FindBin::Bin/support_scripts/eXpress_trans_to_gene_results.pl $results_xprs_file $gene_to_trans_map_file > $results_xprs_file.genes";
    &process_cmd($cmd);
    
    return;
}


####
sub run_RSEM {
    my ($bam_file, $rsem_prefix, $output_prefix) = @_;
        

    unless ($output_prefix) {
        $output_prefix = "RSEM";
    }
    
    my $keep_intermediate_files_opt = ($DEBUG_flag) ? "--keep-intermediate-files" : "";
    
    if ($fragment_length) {
        $fragment_length = "--fragment-length-mean $fragment_length";
    }
    else {
        $fragment_length = "";
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

    my $paired_flag_text = ($paired_flag) ? "--paired-end" : "";

    my $rsem_bam_flag = ($include_rsem_bam) ? "" : "--no-bam-output";


    my $cmd = "rsem-calculate-expression $no_qualities_string "
        . "$paired_flag_text "
        . " $rsem_add_opts "
        . "-p 4 "
        . "$fragment_length "
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

            $file = "<(zcat $file)";
    
            
        }
        push (@files, $file);
    }

    $file_listing = join(",", @files);

    return($file_listing);
}
