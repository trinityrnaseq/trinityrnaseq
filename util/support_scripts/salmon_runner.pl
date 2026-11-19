#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib ("$FindBin::Bin/../../PerlLib");
use Process_cmd;

my $usage = "\n\n\tusage: $0 Trinity.fasta reads.fa [threads=1]\n\n";

my $trin_fa = $ARGV[0] or die $usage;
my $reads_fa = $ARGV[1] or die $usage;
my $CPU = $ARGV[2] || 1;


my $ZERO_ASSIGNED_MSG = "Salmon was only able to assign 0 fragments to transcripts";

main: {

    my $salmon_major_version = &get_salmon_major_version();
    
    my $salmon_index = "$trin_fa.salmon.idx";
    my $salmon_stderr = "_salmon.$$.stderr";
    my $cmd = "salmon --no-version-check index -t $trin_fa -i $salmon_index -k 25 -p $CPU  > $salmon_stderr 2>&1";
    &run_cmd_capture_stderr($cmd, $salmon_stderr);

    ## note: selective alignment (formerly --validateMappings) is the default as of salmon 1.0
    $cmd = "salmon --no-version-check quant -i $salmon_index -l U -r $reads_fa -o salmon_outdir -p $CPU";
    if ($salmon_major_version < 2) {
        # --minAssignedFrags was removed in salmon 2.0 (rust rewrite)
        $cmd .= " --minAssignedFrags 1";
    }
    $cmd .= " > $salmon_stderr 2>&1";
    
    eval {
        &run_cmd_capture_stderr($cmd, $salmon_stderr);
    };
    if ($@) {
        my $errmsg = $@;
        # salmon >= 2.0 writes its outputs but exits non-zero when no fragments map.
        # Report it the same way salmon 1.x does so Trinity can recognize it.
        if ($salmon_major_version >= 2 && &no_fragments_mapped("salmon_outdir")) {
            die "$ZERO_ASSIGNED_MSG.\n$errmsg";
        }
        die $errmsg;
    }

    unlink($salmon_stderr);
    
    exit(0);

}


####
sub run_cmd_capture_stderr {
    my ($cmd, $stderr_file) = @_; 
    eval {
        &process_cmd($cmd);
    };
    if ($@) {
        my $errmsg = `cat $stderr_file`;
        die "Error, cmd: $cmd failed with msg: $errmsg $@";
    }
    return;
}

####
sub get_salmon_major_version {
    my $version_info = `salmon --version 2>&1`;
    if ($version_info =~ /(\d+)\.(\d+)\.(\d+)/) {
        return($1);
    }
    die "Error, cannot determine salmon version from: $version_info";
}

####
sub no_fragments_mapped {
    my ($salmon_outdir) = @_;

    my $meta_info_file = "$salmon_outdir/aux_info/meta_info.json";
    if (-s $meta_info_file) {
        open(my $fh, $meta_info_file) or return(0);
        my $json = do { local $/; <$fh> };
        close $fh;
        if ($json =~ /\"num_mapped\"\s*:\s*0\s*[,\}]/) {
            return(1);
        }
    }
    return(0);
}
