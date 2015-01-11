#!/usr/bin/env perl

use strict;
use warnings;
use Cwd;
use Carp;

my $usage = "usage: $0 left_1.fastq [right_2.fastq]\n\n";

my $left_fq = $ARGV[0] or die $usage;
my $right_fq = $ARGV[1];


my $allpathslg_path = $ENV{ALLPATHSLG_BASEDIR} or die "Error, must have ALLPATHSLG_BASEDIR env variable set to AllpathsLG installation directory";

main: {

	my $workdir = cwd();
	
	## make full paths for inputs.
    if ($left_fq !~ /^\//) {
		$left_fq = "$workdir/$left_fq";
	}
	if ($right_fq && $right_fq !~ /^\//) {
		$right_fq = "$workdir/$right_fq";
	}
	
	my $errCorDir = "errorCorrectionDir.$$";
	mkdir ($errCorDir) or die "Error, cannot mkdir $errCorDir";
	
	if ($right_fq) {
		&process_cmd("cat $left_fq $right_fq > $errCorDir/target.fastq");
	}
	else {
		&process_cmd("ln -s $left_fq $errCorDir/target.fastq");
	}
	
	chdir $errCorDir or die "Error, cannot cd to $errCorDir";
	
	## run the AllpathsLG steps

	my $cmd = "$allpathslg_path/bin/FastqToFastbQualb FASTQ=target.fastq OUT_HEAD=target";
	&process_cmd($cmd);

	mkdir("target") or die "Error, cannot mkdir target";
	rename("target.fastb", "target/target.fastb") or die $!;
	rename("target.qualb", "target/target.qualb") or die $!;
	
	$cmd = "$allpathslg_path/bin/FindErrors K=25 IN_HEAD=target OUT_EDIT_HEAD=target_edit OUT_CORR_HEAD=target_corr";


	if (my $num_threads = $ENV{OMP_NUM_THREADS}) {
		$cmd .= " NUM_THREADS=$num_threads";
	}

	&process_cmd($cmd);	

	chdir("target_edit") or die "Error, cannot cd to target_edit/";

	## convert the error-corrected sequences back to fastq format:
	$cmd = "$allpathslg_path/bin/FastbAndQualb2Fastq HEAD=target_edit";
	&process_cmd($cmd);

	## restore the original accession names.
	my $error_corrected_fastq = "target_edit.fastq";
	my $acc_adjusted_file = $error_corrected_fastq . ".accAdj";
	&replace_accessions("../target.fastq", $error_corrected_fastq, $acc_adjusted_file);
	
	## write final outputs, separating pairs again if paired.
	my $new_left_fq = "$left_fq.ErrCor.fq";
	if ($right_fq) {
		
		my $new_right_fq = "$right_fq.ErrCor.fq";
		
		open (my $left_fh, ">$new_left_fq") or die "Error, cannot write to $new_left_fq";
		open (my $right_fh, ">$new_right_fq") or die "Error, cannot write to $new_right_fq";
		
		open (my $fh, $acc_adjusted_file);
		while (my @record = &get_record($fh)) {
			my $acc_line = $record[0];
			chomp $acc_line;
			if ($acc_line =~ /\/1$/) {
				print $left_fh join("", @record);
			}
			elsif ($acc_line =~ /\/2$/) {
				print $right_fh join("", @record);
			}
			else {
				croak "Error, cannot determine if left or right fragment read of pair:\n" . join("", @record);
			}
		}
		close $fh;
		close $left_fh;
		close $right_fh;
		
		print STDERR "Done. Error-corrected files provided as:\n"
			. "$new_left_fq\n"
			. "$new_right_fq\n\n";
		
	}
	else {
		rename($acc_adjusted_file, $new_left_fq) or die "Error, cannot rename $acc_adjusted_file to $new_left_fq";
		print STDERR "Done.  Error-corrected file provided as:\n$new_left_fq\n\n";
		
	}
	
	
	exit(0);
	


}

####
sub replace_accessions {
	my ($input_fq, $corrected_fq, $output_fq) = @_;

	open (my $ofh, ">$output_fq") or die "Error, cannot write to $output_fq";
	
	open (my $in_fq_fh, $input_fq) or die "Error, cannot open file $input_fq";
	open (my $cor_fq_fh, $corrected_fq) or die "Error, cannot open file $corrected_fq";
	
	while (my @original_record = &get_record($in_fq_fh)) {
		my @corrected_record = &get_record($cor_fq_fh);

		$corrected_record[0] = $original_record[0];
		
		print $ofh join("", @corrected_record);
	}

	close $ofh;
	close $in_fq_fh;
	close $cor_fq_fh;

	return;
}

####
sub get_record {
	my ($fh) = @_;

	my @record;
	for (1..4) {
		my $line = <$fh>;
		push (@record, $line) if ($line);
	}

	
    if (@record && $record[0] !~ /^\@/) {
		croak "Error, fastq record doesn't start with an accession: " . join("", @record);
	}

	return(@record);
}


####
sub process_cmd {
	my ($cmd) = @_;

	print "CMD: $cmd\n";

	my $ret = system($cmd);
	if ($ret) {
		die "Error, command $cmd died with ret $ret";
	}

	return;
}

