#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib ("$FindBin::Bin/../../PerlLib");
use SAM_entry;

my $usage = "\n\nusage: $0 nameSorted.sam iworm_fasta_file\n\n";

my $name_sorted_sam_file = $ARGV[0] or die $usage;
my $iworm_fasta = $ARGV[1] or die $usage;


my $SAM_OFH;


main: {
        
    open ($SAM_OFH, ">scaffolding_entries.sam") or die $!;
    
    my %iworm_acc_to_fasta_index;
    {
        
        my $counter = 0;
        open (my $fh, $iworm_fasta) or die "Error, cannot open file $iworm_fasta";
        while (<$fh>) {
            if (/^>(\S+)/) {
                my $acc = $1;
                $iworm_acc_to_fasta_index{$acc} = $counter; # starts at zero
                $counter++;
            }
        }
        close $fh;
    }
    

    my %paired_iworm_contigs;

    my $prev_core_acc = "";
    my %end_to_iworm;
    
    my $num_warnings = 0;

    my $num_unrecognized_iworm_contig_names = 0;

    if ($name_sorted_sam_file =~ /\.bam$/) {
        $name_sorted_sam_file = "samtools view $name_sorted_sam_file | ";
    }
    
    open (my $fh, $name_sorted_sam_file) or die "Error, cannot open file $name_sorted_sam_file";
    while (<$fh>) {
        my $line = $_;
        chomp;
        my @x = split(/\t/);
        my $read_acc = $x[0];
        my $iworm_acc = $x[2];

        unless (defined $iworm_acc) { next; }
        
        if ($iworm_acc eq '*') { 
            # unmapped read
            next;
        }
        
        unless ($iworm_acc =~ /^a\d+;\d+/) {
            $num_unrecognized_iworm_contig_names++;
            if ($num_unrecognized_iworm_contig_names <= 10) {
                print STDERR "warning, inchworm contig ($iworm_acc) isn't recognized. Ignoring entry: [[$line]]\n";
            }
            if ($num_unrecognized_iworm_contig_names == 11) {
                print STDERR "warning, too many unrecognized inchworm contig names.  Will report summary of counts later.\n";
            }
            next;
        }


        my $core_acc;
        my $frag_end;
        

        if ($read_acc =~ /^(\S+)\/([12])$/) { 
            
            $core_acc = $1;
            $frag_end = $2;
        }
        
        # capture long read mappings
        elsif ($read_acc =~ /^(LR\$\|\S+)/) {
            $core_acc = $1;
            $frag_end = "LR";
        }
        
        else {
            # must have mixed in a single read with the pairs...
            if ($num_warnings <= 10) {
                print STDERR "warning, ignoring read: $read_acc since cannot decipher if /1 or /2 of a pair.\n";
            }
            elsif ($num_warnings == 11) {
                print STDERR "number of read warnings exceeded 10.  Turning off warning messages from here out.\n";
            }
            $num_warnings++;
            
            next;
        }
        
        if ($core_acc ne $prev_core_acc) {
            
            &examine_frags(\%end_to_iworm, \%paired_iworm_contigs);
            %end_to_iworm = ();
            
        }

        $end_to_iworm{$frag_end}->{$iworm_acc} = $line;
        
        $prev_core_acc = $core_acc;

    }
    close $fh;

    ## get last one
    &examine_frags(\%end_to_iworm, \%paired_iworm_contigs);
    
    if ($num_warnings) {
        print STDERR "WARNING: note there were $num_warnings reads that could not be deciphered as being /1 or /2 of a PE fragment.  Hopefully, these were SE reads that should have been ignored. Otherwise, please research this further.\n\n";
    }
    if ($num_unrecognized_iworm_contig_names) {
        print STDERR "WARNING: note, there were $num_unrecognized_iworm_contig_names inchworm contig names in the SAM file that were ignored due to the inchworm contig accession not being recognized.\n";
    }
    
    foreach my $pairing (reverse sort {$paired_iworm_contigs{$a}<=>$paired_iworm_contigs{$b}} keys %paired_iworm_contigs) {

        my ($iworm_acc_A, $iworm_acc_B) = split(/\t/, $pairing);
      
        my $index_A = $iworm_acc_to_fasta_index{$iworm_acc_A};
        unless (defined $index_A) {
            print STDERR "WARNING, no index for iworm acc: $iworm_acc_A\n";
            next;
        }
        my $index_B = $iworm_acc_to_fasta_index{$iworm_acc_B};
        unless (defined $index_B) {
            print STDERR "WARNING, no index for iworm acc: $iworm_acc_B\n";
            next;
        }
        
        print join("\t", $iworm_acc_A, $index_A, $iworm_acc_B, $index_B, $paired_iworm_contigs{$pairing}) . "\n";
    
    }
    

    close $SAM_OFH;
    
    exit(0);
}

####
sub examine_frags {
    my ($end_to_iworm_href, $paired_iworm_contigs_href) = @_;

    if (exists ($end_to_iworm_href->{LR})) {
        my @LR_read_names = keys %{$end_to_iworm_href->{LR}};
        
        for (my $i = 0; $i < $#LR_read_names; $i++) {
            for (my $j = $i + 1; $j <= $#LR_read_names; $j++) {
                
                my $pair  = join("\t", sort ($LR_read_names[$i], $LR_read_names[$j]) );
                $paired_iworm_contigs_href->{$pair}++;

                #print STDERR "-got LR pair: $pair\n";
            }
        }
        
        return;
    }
    

    ## handle the paired-end reads
    
    my @ends = keys %$end_to_iworm_href;

    unless (scalar @ends == 2) {
        ## no pairs
        return;
    }
    
    
    my @iworm_left = keys %{$end_to_iworm_href->{1}};
    my @iworm_right = keys %{$end_to_iworm_href->{2}};

    
    if (scalar(@iworm_left) == 1 && scalar(@iworm_right) == 1
        &&
        $iworm_left[0] ne $iworm_right[0]) {
        
        ## simplest case, pairs each aligning to different contigs.
    
        my $i_left = $iworm_left[0];
        my $i_right = $iworm_right[0];

        my $pair = join("\t", sort ($i_left, $i_right));
        #print STDERR "Got pairing: $pair\n";
        $paired_iworm_contigs_href->{$pair}++;
    
    }
    else {
        # examine split reads
        
        &check_for_split_reads($end_to_iworm_href->{1}, $paired_iworm_contigs_href);
        &check_for_split_reads($end_to_iworm_href->{2}, $paired_iworm_contigs_href);
        
        
    }
        
    ## output all alignments for these reads.
    foreach my $sam_line (values %{$end_to_iworm_href->{1}}, values %{$end_to_iworm_href->{2}}) {
        print $SAM_OFH $sam_line;
    }
        
    
    return;
}

        
####
sub check_for_split_reads {
    my ($read_mappings_href, $paired_iworm_contigs_href) = @_;

    my @iworm_names = keys %$read_mappings_href;

    if (scalar(@iworm_names) < 2) {
        # nothing to do here.
        return;
    }
    
    ## see if we have a read where different parts of the read are aligning to different iworm contigs.
    for (my $i = 0; $i < $#iworm_names; $i++) {
        my $iworm_name_i = $iworm_names[$i];
        my $align_i = new SAM_entry($read_mappings_href->{$iworm_name_i});
        
        for (my $j = $i + 1; $j <= $#iworm_names; $j++) {
            my $iworm_name_j = $iworm_names[$j];
            my $align_j = new SAM_entry($read_mappings_href->{$iworm_name_j});

            if (&overlap_not_involving_containment($align_i, $align_j)) {
                ## store pair
                
                my $pair = join("\t", sort($iworm_name_i, $iworm_name_j));
                $paired_iworm_contigs_href->{$pair}++;
                
            }
        }
    }

    return;
}
    
####
sub overlap_not_involving_containment {
    my ($align_i, $align_j) = @_;
        
    my ($read_i_lend, $read_i_rend) = $align_i->get_read_span();
            
    my ($read_j_lend, $read_j_rend) = $align_j->get_read_span();
    
    #print join("\t", "$read_i_lend-$read_i_rend", "$read_j_lend-$read_j_rend") . "\n";
    
    if (  ($read_i_lend <= $read_j_lend && $read_i_rend >= $read_j_rend) # i encapsulates j
          ||
          ($read_j_lend <= $read_i_lend && $read_j_rend >= $read_i_rend)  # j encapsulates i
        ) {
        ## containment
        #print "\t* containment\n";
        return(0);
    }
    elsif ( $read_i_lend < $read_j_rend && $read_i_rend > $read_j_lend) {
        #print "\t OVERLAP\n";
        return(1);
    }
    else {
        #print "\tok\n";
        
        return(1);
    }
        
}
    
