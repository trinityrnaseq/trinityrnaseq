#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;

use lib ("$ENV{EUK_MODULES}");
use Gene_obj;


while (<>) {
    chomp;
    my ($scaff, $read_acc, $lend, $rend) = split(/\t/);
    
    my $gene_obj = new Gene_obj();
    $gene_obj->populate_gene_object({}, {$lend => $rend});
    $gene_obj->{asmbl_id} = $scaff;
    
    $gene_obj->{TU_feat_name} = $read_acc;
    $gene_obj->{Model_feat_name} = $read_acc;


    print $gene_obj->to_BED_format();

    
    
}

exit(0);



