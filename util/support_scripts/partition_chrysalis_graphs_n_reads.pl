#!/usr/bin/env perl

use strict;
use warnings;
use File::Basename;
use Getopt::Long qw(:config no_ignore_case bundling pass_through);

my $usage = <<__EOUSAGE__;

#########################################################################
#
#  --deBruijns <string>           bundled_iworm_contigs.fasta.deBruijn
#  --componentReads <string>      readsToComponents.out.sort
#  
#  -N <int>                       number of graphs per partition
#  -L <int>                       min contig length
#
##########################################################################

__EOUSAGE__


    ;



my $help_flag;
my $deBruijns_file;
my $componentReads_file;
my $num_graphs_per_partition;
my $min_contig_length;
my $components_directory;

&GetOptions ( 'h' => \$help_flag,
              'deBruijns=s' => \$deBruijns_file,
              'componentReads=s' => \$componentReads_file,
              'N=s' => \$num_graphs_per_partition,
              'L=s' => \$min_contig_length,
	      'compdir=s' => \$components_directory
              
    );


if ($help_flag) {
    die $usage;
}

unless ($deBruijns_file && $componentReads_file && $num_graphs_per_partition && $min_contig_length) {
    die $usage;
}

main: {

    print "Partitioning chrysalis graphs and reads\n";
    my $outdir_base = dirname($deBruijns_file);
    my $components_directory = $outdir_base . '/Component_bins';

    my %component_id_to_partition_base;

    my $component_reader = Component_reader->new($deBruijns_file);

    my $num_components = 0;

    unless (-e $components_directory) {
        mkdir $components_directory or die "Error, cannot mkdir $components_directory";
    }
    
    while (my $component = $component_reader->next_component()) {

        my $id = $component->{component_id};
        my $num_kmers = $component->{num_kmers};
        my $graph_text = $component->{graph_text};

        #print join("\t", $id, $num_kmers) . "\n";

        if ($num_kmers + 24 < $min_contig_length) { # assuming the 25-mer kmer size
            next;
        }
        
        $num_components++;
        my $outdir = "$components_directory/Cbin" . int($num_components/$num_graphs_per_partition);
        unless (-d $outdir) {
            mkdir $outdir or die "Error, cannot mkdir $outdir";
        }
        my $component_file = "$outdir/c$id.graph.tmp";
        print STDERR "\r[$num_components] writing to $component_file     " if $num_components % 1000 == 0;
        open (my $ofh, ">$component_file") or die "Error, cannot write to $component_file";
        print $ofh $graph_text;
        close $ofh;

        $component_id_to_partition_base{$id} = "$outdir/c$id";
        
    }
    
    print STDERR "\n\nDone partitioning graphs.\nPartitioning reads...\n\n";
    
    ## Now, partition the reads
    my $prev_comp = -1;
    open (my $fh, "$componentReads_file") or die $!;
    my $outfh;
    while (my $line = <$fh>) {
        chomp $line;
        my ($comp, $acc, $pct, $read) = split(/\t/, $line);
        if ($comp != $prev_comp) {
            ## see if we need to capture these reads
            $outfh = undef;
            if (my $base = $component_id_to_partition_base{$comp}) {
                my $outfile = $base . ".reads.tmp";
                open ($outfh, ">$outfile") or die "Error, cannot write to $outfile";
            }
        }
        if ($outfh) {
            print $outfh "$acc $pct\n$read\n";
        }
        $prev_comp = $comp;
    }
    close $outfh if $outfh;
    close $fh;
        
    print STDERR "Done partitioning reads.\n\n";
    
    ## identify components to target for quantifygraph:
    {
	my $list_output = $outdir_base.'/component_base_listing.txt';
        print STDERR "-writing $list_output\n";
        
        open (my $ofh, ">$list_output") or die $!;
        foreach my $component_id (sort {$a<=>$b} keys %component_id_to_partition_base) {
            
            my $base = $component_id_to_partition_base{$component_id};
            my $tmp_graph_file = $base . ".graph.tmp";
            my $tmp_reads_file = $base . ".reads.tmp";
            if (-s $tmp_graph_file && -s $tmp_reads_file) {

                print $ofh join("\t", $component_id, $base) . "\n";
            }
            else {
                print STDERR "Warning: no reads written for entry: $tmp_graph_file\n";
            }
        }
        close $ofh;
        
        print STDERR "\nDone.\n\n";
    }
    
    exit(0);
}


##########################
##########################
package Component_reader;

use strict;
use warnings;
use Carp;

sub new {
    my ($packagename) = shift;
    my ($filename) = @_;

    my $self = { filename => $filename,
                 fh => undef,
                 prev_line => undef,
    };

    bless ($self, $packagename);

    $self->_init();
    
    return($self);
}

####
sub _init {
    my $self = shift;
    
    open (my $fh, $self->{filename}) or die "Error, cannot open file $self-><{filename}";
    $self->{fh} = $fh;
    $self->{prev_line} = <$fh>;

    unless ($self->{prev_line} =~ /^Component/) {
        confess "Error, did not extract component line from file: $self->{filename}";
    }
    
    return;
    
    
}


####
sub next_component {
    my $self = shift;

    my @lines;
    my $fh = $self->{fh};
    
    my $component_id;
    my $component_line = $self->{prev_line};
    if ($component_line =~ /Component (\d+)/) {
        $component_id = $1;
    }
    

    while (my $line = <$fh>) {
        if ($line =~ /^Component/) {
            $self->{prev_line} = $line;
            last;
        }
        else {
            push (@lines, $line);
        }
    }
    if (@lines) {
        my $num_kmers = scalar(@lines);
        my $graph_text = join("", $component_line, @lines);  # newlines already included
        my $component = Component->new($component_id, $num_kmers, $graph_text);
        return($component);
    }
    else {
        return(undef);
    }
}



####################
####################
package Component;

use strict;
use warnings;
use Carp;

sub new {
    my $packagename = shift;
    my ($component_id, $num_kmers, $graph_text) = @_;
    
    my $self = { component_id => $component_id,
                 num_kmers => $num_kmers,
                 graph_text => $graph_text,
    };

    bless ($self, $packagename);

    return($self);
}


