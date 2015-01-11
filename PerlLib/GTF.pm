# -*- Mode: Perl; tab-width: 8; perl-indent-level: 2; indent-tabs-mode: nil -*- 
use strict;
###############################################################################
# GTF
###############################################################################
# The GTF object parses and stores all data from a gtf file.
package GTF;
use Carp;
# GTF::new(hash)
#    This is the contructor for GTF objects.  It takes one required argument
#    which is a hash containing the following optional fields:
#        gtf_filename: The filename of the gtf file to be loaded.  If no filename 
#                      is given it just creates an empty object.
#        seq_filename: The filename of the sequence corresponding to this gtf file.
#                      If given additional validity checks are made, additional stats
#                      are kept and the transcript can be output.
#        tx_out_fh:    The filehandle to output the genes transcripts to.  Only works
#                      If seq_filename is given.
#        conseq_filename: The filename of the conservation sequence corresponding to 
#                      this gtf file.  If given additional stats are kept.
#        warning_fh:   The filehandle to output gtf format warnings to.  If none 
#                      is given warnings are disregarded.
sub new {
  my ($info) = @_;
  my $gtf = bless {Genes => [], Transcripts => [], Exons => [], CDS   => [], 
                   Inter => [], Inter_CNS => [], Filename => "", Sequence => "",
                   Comments => "", Modified => 0};
  $gtf->{Total_Conseq} = {0 => -1,
                          1 => -1,
                          2 => -1};
  $gtf->{Total_Seq} = {A => -1,
                       C => -1,
                       G => -1,
                       T => -1,
                       N => -1};
  if(defined($info->{seq_filename})){
    if(-e $info->{seq_filename}){
      $gtf->{Sequence} = $info->{seq_filename};
    }
    else{
      $gtf->{Sequence} = 0;
      print STDERR "GTF: $info->{seq_filename} does not exist.\n"; 
    }
  }
  else{
    $gtf->{Sequence} = 0;
  }
  if(defined($info->{strict})){
    $gtf->{Strict} = $info->{strict};
  }
  else{
    $gtf->{Strict} = 0;
  }
  if(defined($info->{tx_out_fh})){
    if(defined($gtf->{Sequence})){
      $gtf->{Tx} = $info->{tx_out_fh};
    }
    else{
      $gtf->{Tx} = 0;
      print STDERR "GTF: Cannot ouput transcripts without the sequence.\n"; 
    }
  }
  else{
    $gtf->{Tx} = 0;
  }
  if(defined($info->{conseq_filename})){
    if(-e $info->{conseq_filename}){
      $gtf->{Conseq} = $info->{conseq_filename};
    }
    else{
      $info->{Conseq} = 0;
      print STDERR "GTF: $info->{conseq_filename} does not exist.\n"; 
    }
  }
  else{
    $info->{Conseq} = 0;
  }
  if((defined($info->{inframe_stops})) && ($info->{inframe_stops})){
    $gtf->{Inframe_Stops} = 1;
  }
  else{
    $gtf->{Inframe_Stops} = 0;
  }
  if((defined($info->{fix_gtf})) &&
     ($info->{fix_gtf} == 1)){
    $gtf->{Fix_GTF} = 1;
  }
  else{
    $gtf->{Fix_GTF} = 0;
  }
  if(defined($info->{warning_skips})){
    $gtf->{Warning_Skips} = $info->{warning_skips};
  }
  else{
    $gtf->{Warning_Skips} = [];
  }
  if(defined($info->{bad_list})){
    $gtf->{Bad_List} = $info->{bad_list};
  }
  else{
    $gtf->{Bad_List} = [];
  }
  if(defined($info->{detailed_error_count})){
    $gtf->{Detailed_Error_Count} = $info->{detailed_error_count};
  }
  if(defined($info->{mark_ase})) {
    $gtf->{Mark_ASE} = $info->{mark_ase};
  }
  else {
    $gtf->{Mark_ASE} = 0;
  }
  if(defined($info->{gtf_filename})){
    if(-e $info->{gtf_filename}){
      $gtf->{Filename} = $info->{gtf_filename};
      if((defined($info->{no_check})) && ($info->{no_check})){
        $gtf->_load_no_check($info->{gtf_filename});
      }
      else{
        $gtf->_parse_file($info->{warning_fh});
      }
      if($gtf->{Mark_ASE}) {
        $gtf->mark_ase;
      }
    }
    else{
      print STDERR "GTF: $info->{gtf_filename} does not exist.\n"; 
    }		    
  }
  else{
    $gtf->{Filename} = "";
  }
  return $gtf;
}

# GTF::_update
#    This function resorts all lists.  It is called when a list is requested.
#    It does nothing unless the $gtf->{Modified} value is 1.  Sets 
#    $gene->{Modified} to 0;
sub _update{
  my ($gtf) = @_;
  unless($gtf->{Modified}){
    return;
  }
  $gtf->{Modified} = 0;
  my $genes = $gtf->{Genes};
  my $txs = $gtf->{Transcripts};
  my $cds = $gtf->{CDS};    
  my $inter = $gtf->{Inter};
  $genes = [sort {$a->start <=> $b->start} @$genes];
  $txs = [sort {$a->start <=> $b->start} @$txs];
  $cds = [sort {$a->start <=> $b->start} @$cds];
  $inter = [sort {$a->start <=> $b->start} @$inter];
  $gtf->{Genes} = $genes;
  $gtf->{Transcripts} = $txs;
  $gtf->{CDS} = $cds;
  $gtf->{Inter} = $inter;
}

# GTF::genes()
#    This function returns an array of refernces to GTF::Gene objects.
#    The array contains one reference for each gene in the gtf file.
#    The array is sorted by the start position of each gene.
sub genes{
  my ($gtf) = @_;
  $gtf->_update;
  return $gtf->{Genes};
}

sub transcripts{
  my ($gtf) = @_;
  $gtf->_update;
  return $gtf->{Transcripts};
}

# GTF::add_gene(gene_ref)
#    This function take a reference to a GTF::Gene object and adds it to the 
#    list of genes stored by this object.
sub add_gene{
  my ($gtf,$gene) = @_;
  my $genes = $gtf->{Genes};
  my $txs = $gtf->{Transcripts};
  my $cds = $gtf->{CDS};    
  push @$genes, $gene;
  push @$txs, @{$gene->transcripts};
  push @$cds, @{$gene->cds};
  $gtf->{Modified} = 1;
}

# GTF::add_feature(feature_ref,seqname,source,strand,id)
#   This function will take a reference to an intergenic feature (such as 
#   inter or inter_CNS) and add it to the appropriate feature list.
#   Pass a sequence name, the strand ("+" is OK for intergenic)
#   and an ID to identify the feature.
#   Note that the feature must be of the types inter or inter_CNS.
sub add_feature{
  my ($gtf,$feature,$seqname,$source,$strand,$id) = @_;
  die "Only inter or inter_CNS features may be directly added to a GTF.\n"
    if($feature->type ne "inter" && $feature->type ne "inter_CNS");
  my $gene = GTF::Gene::new($id,$seqname,$source,$strand);
  my $tx = GTF::Transcript::new($id);
  $gene->add_transcript($tx);
  $tx->add_feature($feature);
  my $features;
  if($feature->type eq "inter")         { $features = $gtf->{Inter}     }
  elsif($feature->type eq "inter_CNS")  { $features = $gtf->{Inter_CNS} }
  push @$features, $feature;
  $gtf->{Modified} = 1;
}

# GTF::set_genes(gene_list)
#    This function take a reference to a list of GTF::Gene objects and adds each to the 
#    list of genes stored by this object.
sub set_genes{
  my ($gtf,$genes) = @_;
  $gtf->{Genes} = $genes;
  my @txs;
  my @cds;
  foreach my $gene (@$genes){
    push @txs, @{$gene->transcripts};
    push @cds, @{$gene->cds};
  }
  $gtf->{Transcripts} = \@txs;
  $gtf->{CDS} = \@cds;
  $gtf->{Modified} = 1;
}

# GTF::remove_gene(gene_id)
#    This function takes a gene_id and removes the gene with that gene_id
#    from the list of genes stored by this object, if a gene with that gene_id
#    can be found.  Otherwise does nothing.  Returns true if a gene was removed
#    and false if no gene with that gene_id was found.
sub remove_gene{
  my ($gtf,$gene_id) = @_;
  my $genes = $gtf->{Genes};
  my $removed = 0;
  my @new_genes;
  my @new_txs;
  my @new_cds;
  foreach my $gene (@$genes){
    if($gene->gene_id eq $gene_id){
      $removed++;
    }
    else{
      push @new_genes, $gene;
      push @new_txs, @{$gene->transcripts};
      push @new_cds, @{$gene->cds};
    }
  }
  if($removed){
    $gtf->{Genes} = \@new_genes;
    $gtf->{Transcripts} = \@new_txs;
    $gtf->{CDS}   = \@new_cds;
  }
  return $removed;
}

# GTF::remove_feature(feature)
#   This function will remove a feature from the list of features.
#   Note the feature must be of type inter or inter_CNS.
sub remove_feature {
    my($gtf,$feature) = @_;
    my $removed = 0;
    my $feature_key;
    if($feature->type eq "inter") {
        $feature_key = "Inter";
    }
    elsif($feature->type eq "inter_CNS") {
        $feature_key = "Inter_CNS";
    }
    my @new_features;
    foreach my $stored_feature (@{$gtf->{$feature_key}}) {
        if($stored_feature == $feature) {
            $removed++;
        }
        else {
            push @new_features, $stored_feature;
        }
    }
    if($removed) {
        $gtf->{$feature_key} = \@new_features;
    }
    return $removed;
}

# GTF::set_filename(filename)
#    Sets the filename to be returned by the filename function 
sub set_filename{
  my ($gtf,$filename) = @_;
  $gtf->{Filename} = $filename;
}

# GTF::offset(ammount)
#    This function will offset all genes in this gene by the given ammount
sub offset{
  my ($gtf,$offset) = @_;
  unless(defined($offset)){
    print STDERR "Undefined value passed to GTF::offset.\n";
    return;
  }
  my $genes = $gtf->{Genes};
  foreach my $gene (@$genes){
    $gene->offset($offset);
  }    
  $gtf->{Modified} = 1;
}

# GTF::reverse_complement(seq_length)
#    This function takes the length of the sequence this gtf file was based on and 
#    and reverse complements everything in the file.  So the positive strand becomes
#    the negative strand and the files starts at seq_length and moves downward to 0.
sub reverse_complement{
  my ($gtf, $seq_length) = @_;
  unless(defined($seq_length)){
    print STDERR 
        "Undefined value for seq_length passed to GTF::reverse_complement.\n";
  }
  my $genes = $gtf->{Genes};
  foreach my $gene(@$genes){
    $gene->reverse_complement($seq_length);
  }
  $gtf->{Modified} = 1;
}

# GTF::cds()
#    This function returns an array of refernces to GTF::Feature objects. 
#    The array contains one refernece to each CDS feature in the gtf file.
#    The array is sorted by the <start> of each CDS.
sub cds{
  my ($gtf) = @_;
  $gtf->_update;
  return $gtf->{CDS};
}

sub inter {
    my ($gtf) = @_;
    $gtf->_update;
    return $gtf->{Inter};
}

sub inter_cns {
    my ($gtf) = @_;
    $gtf->_update;
    return $gtf->{Inter_CNS};
}

# GTF::mark_ase()
# If you are looking for alternative splicing events in your GTF file,
# call this function to marke them.  This function is separate from the 
# _update methods in order to keep the running time low.

sub mark_ase {
  my ($gtf) = @_;

  #set the alternative splicing event flags for each feature
  #currently, only optional coding exons

  # The main idea is to check all pairs of transcripts in a gene,
  # and find any incident where there are two exons that have the same
  # end coord and two exons which have the same start coord
  # and that an exon in only one of the transcripts lies between them.
  # It is accomplished by sorting all the exons in both transcripts
  # and checking for a supercasette, i.e.
  #
  # ===|||||||||===||||||||======||||||||||====
  # ==||||||||||=================||||||||||||==
  #             <--------------->
  #               supercassette
  #
  # Note that it doesn't matter if the opposite ends of the exon are equal.
  
  # rpz

  my @genes = @{$gtf->genes};
  for my $gene(@{$gtf->genes}) {
    my @txs = (@{$gene->transcripts}); 
    if(scalar(@txs) > 1) {
      for (my $i = 0; $i < scalar(@txs)-1; $i++) {
        next if(scalar(@{$txs[$i]->cds}) == 0);
        for (my $j = $i+1; $j < scalar(@txs); $j++) {
          next if(scalar(@{$txs[$j]->cds}) == 0);

          my @ocds = sort {$a->start <=> $b->start || $a->stop <=> $b->stop}
            (@{$txs[$i]->cds}, @{$txs[$j]->cds});

          for(my $k = 4; $k < scalar(@ocds); $k++) {
            if(     $ocds[$k-4]->stop == $ocds[$k-3]->stop &&
                    $ocds[$k-1]->start == $ocds[$k]->start) {
              if($ocds[$k-2]->length%3 == 0) {
                $ocds[$k-2]->set_ase("InframeOptional");
              }
              else {
                $ocds[$k-2]->set_ase("FrameshiftOptional");
              }
            }
          }
        }
      }
    }
  }

}

# GTF::infer_exons()
# Skips infering if Exons already exist
# Output is undefined if features overlap, except for start_codon
# Currently knows about utr5, start,stop,cds, utr3
sub infer_exons {
  my ($gtf) = @_; 
  $gtf->_update;

  for my $tx (@{$gtf->transcripts}) {
    next if (@{$tx->exons});

    if ($tx->strand eq '+') {    
      for  ( @{$tx->utr5}, @{$tx->cds}, @{$tx->stop_codons}, @{$tx->utr3} ) {

        if (!@{$tx->exons} ||  ${$tx->exons}[-1]->stop + 1 < $_->start ) {
          $tx->add_feature( GTF::Feature::new('exon', $_->start, $_->stop, 0, 0) );

        }  elsif ( ${$tx->exons}[-1]->stop + 1 eq $_->start ) {
          ${$tx->exons}[-1]->set_stop($_->stop);            
        }
      }
    } else {
      for  ( @{$tx->utr3}, @{$tx->stop_codons}, @{$tx->cds}, @{$tx->utr5} ) {

        if (!@{$tx->exons} || ${$tx->exons}[-1]->stop + 1 < $_->start  ) {
          $tx->add_feature(GTF::Feature::new('exon', $_->start, $_->stop, 0, 0));
        } elsif ( ${$tx->exons}[-1]->stop + 1 eq $_->start ) {
          ${$tx->exons}[-1]->set_stop($_->stop);            
        }
      }
    }
  }
}

# GTF::remove_exons()
# remove exons from all transcripts
sub remove_exons {
  my ($gtf) = @_; 
  $gtf->_update;

  $_->{Exons} = [] for (@{$gtf->transcripts});
}

# GTF::filename()
#    This function returns the filename of the gtf file.
sub filename           {shift->{Filename}}

# GTF::comments()
#    This function returns all full line commentsfrom the gtf file.
sub comments           {shift->{Comments}};

# GTF::conservation_count()
#    Returns a hash containing counts for each conservation character (0,1,2)
sub conservation_count {shift->{Total_Conseq}};

# GTF::sequence_count()
#    Returns a hash containing counts for each sequence character (A,C,T,G,N,X),
#    Where X is all other characters
sub sequence_count {shift->{Total_Seq}};

# GTF::output_gtf_file([file_handle]);
#    This function moves through each gene in the gtf file and outputs 
#    all data in valid gtf2 format.  If takes and optional file_handle
#    as the only argument to which it outputs the data.  If no argument is 
#    given it outputs the data to stdout.
sub output_gtf_file{
  my ($gtf,$out_handle) = @_;
  $gtf->_update;
  unless(defined $out_handle){
    $out_handle = \*STDOUT;
  }
  if($gtf->comments){
    foreach my $comment (@{$gtf->comments}){
      print $out_handle "$comment\n";
    }
  }

  for my $gene (@{$gtf->genes}) { $gene->_update }
  my @top_level_features = sort {$a->start <=> $b->start} 
    (@{$gtf->genes}, @{$gtf->inter_cns}, @{$gtf->inter});
  for my $feature (@top_level_features) {
    $feature->output_gtf($out_handle);
  }
}

# GTF::output_gff_file([file_handle]);
#    This function moves through each gene in the gtf file and outputs 
#    all data in valid GFF format.  If takes and optional file_handle
#    as the only argument to which it outputs the data.  If no argument is 
#    given it outputs the data to stdout.
sub output_gff_file{
  my ($gtf,$out_handle) = @_;
  $gtf->_update;
  unless(defined $out_handle){
    $out_handle = \*STDOUT;
  }
  if($gtf->genes){
    foreach my $gene (@{$gtf->genes}){
      $gene->output_gff($out_handle);
    }
  }
}

# GTF::_parse_file()
#    This function is used internally by the GTF constructor GTF::new
#    to read, parse, and store the GTF file;
sub _parse_file{
  my ($gtf,$warnings) = @_;
  unless(defined($warnings)){
    if(open(WARN, ">/dev/null")){
      $warnings = \*WARN;
    }
    else{
      die "Could not open /dev/null for write\n";
    }
  }
  my $fix_gtf = $gtf->{Fix_GTF};
  my $strict = $gtf->{Strict};
  #open file 
  open(GTF, "<$gtf->{Filename}") 
      or die "Could not open $gtf->{Filename}.\n";
  my (%GeneObj,%TxObj);
  my (@all_genes,@all_txs,@all_cds,@all_inter,@all_cns);
  my $all_line_comments = [];
  #prepare error counts
  my $max_error_count = 5;
  if(defined($gtf->{Detailed_Error_Count})){
    $max_error_count = $gtf->{Detailed_Error_Count};
  }
  my @error_msgs = $gtf->_get_error_messages();
  my @errors;
  my @bad_list;
  my $no_warn = $gtf->{Warning_Skips};
  for(my $i = 0;$i <= $#error_msgs;$i++){
    if($$no_warn[$i]){
      $errors[$i] = $max_error_count;
    }
    else{
      $errors[$i] = 0;
    }
  }
  #read the file
  my $line_num = 0;
  while(my $in_line = <GTF>){
    $line_num++;
    chomp $in_line;
    if($in_line !~ /\S/){
      next;
    }
    if($in_line =~ /^\s*\#/){
      push @$all_line_comments, $in_line;
      next;
    }
    else{
      #remove leading whitespace and get comments
      my $comments = "";	  
      while($in_line =~ /^(.*)\#(.*)$/){
        $in_line = $1;
        $comments = $2.$comments;
      }
      if($in_line =~ /^\s+(.*)$/){
        $in_line = $1;
      }
      my @data = split /\s+/,$in_line;
      #verify line is correct length
      if($#data < 8){
        if($errors[0] < $max_error_count){
          print $warnings
              "Not enough fields on line $line_num.\n";
        }
        $errors[0]++;
        next;
      }
      #check for correct whitespace
      my $field = 1;
      while(($field < 8)&&($in_line =~ /\S+(\s+)/g)){
        unless($1 =~ /^\t$/){
          if($errors[1] < $max_error_count){
            print $warnings 
                "Incorrect type of whitespace between fields ",
                "$field and ",$field + 1, " on line $line_num.  ",
                "Should be tab.\n";
          }
          $errors[1]++;
        }
        $field++;
      }
      #verify correct <feature> field
      $data[2] = lc($data[2]);
      unless(($data[2] eq "cds") || ($data[2] eq "exon")
             || ($data[2] eq "5utr") || ($data[2] eq "3utr")
             || ($data[2] eq "start_codon") 
             || ($data[2] eq "stop_codon")
             || ($data[2] eq "sec") 
             || ($data[2] eq "intron_cns")
             || ($data[2] eq "inter_cns")
             || ($data[2] eq "inter")) {
        if($errors[2] < $max_error_count){
          print $warnings
              "Incorrect value for <feature> field, \"$data[2]\", ",
              "on line $line_num.\n";
        }
        $errors[2]++;
      }		
      if (($data[2] eq "cds") || ($data[2] eq "5utr")
          || ($data[2] eq "3utr") || ($data[2] eq "sec")) {
        $data[2] = uc($data[2]);
      }
      if ($data[2] eq "inter_cns")  { $data[2] = "inter_CNS"  }
      if ($data[2] eq "intron_cns") { $data[2] = "intron_CNS" }
      #verify correct <start> format
      if($data[3] !~ /^\d*$/){
        if($errors[3] < $max_error_count){
          print $warnings 
              "Non-numerical value for <start> field, \"$data[3]\"",
              ", on line $line_num.\n";
        }
        $errors[3]++;
      }	  
      #verify correct <stop> format
      if($data[4] !~ /^\d*$/){		
        if($errors[4] < $max_error_count){
          print $warnings 
              "Non-numerical value for <stop> field, \"$data[4]\"",
              ", on line $line_num.\n";
        }
        $errors[4]++;
      }
      #verify start is < stop
      if($data[3] > $data[4]){
        if($errors[38] < $max_error_count){
          print $warnings 
              "Start field is greater than stop field on line $line_num.\n";
        }
        $errors[38]++;
        my $start = $data[3];
        $data[3] = $data[4];
        $data[4] = $start;
      }
      #verify correct <score> format
      if($data[5] !~ /^-?\d*\.?\d*$/){ # should check for float
        if($data[5] ne '.'){
          if($errors[5] < $max_error_count){
            print $warnings
                "Value for <score> field, \"$data[5]\", is ",
                "neither numerical nor \".\" on line $line_num\n";
          }
          $errors[5]++;
        }
      }
      if(!$strict && $data[5] eq '.'){
        $data[5] = 0;
      }
      #verify correct <strand> field
      unless(($data[6] eq "+")||($data[6] eq "-")||
             ($data[6] eq ".")){
        if($errors[6] < $max_error_count){
          print $warnings
              "Value of <strand> field is \"$data[6]\" on line ",
              "$line_num.  Should be \"+\", \"-\", or \".\".\n";
        }
        $errors[6]++;
        if(!$strict){
          if($data[6] eq "1"){
            $data[6] = "+";
          }
          elsif($data[6] eq "-1"){
            $data[6] = "-";
          }
          elsif($data[6] eq "0"){
            $data[6] = ".";
          }
        }
      }
      #verify correct <frame> field
      unless(($data[7] eq "0")||($data[7] eq "1")||
             ($data[7] eq "2")||($data[7] eq ".")){
        if($errors[7] < $max_error_count){
          print $warnings
              "Value of <frame> field is \"$data[7]\" on line ",
              "$line_num.  Should be \"0\", \"1\", \"2\", or \".\"",
              ".\n";
        }
        $errors[7]++;
      }	    
      #prepare <attributes> and <comment> fields
      my $attribute_string = "";
      if($in_line =~ /^\S+\t\S+\t\S+\t\S+\t\S+\t\S+\t\S+\t\S+\t(.*)$/){
        $attribute_string = "$1 ";
      }	       
      #check for equals
      if($attribute_string =~ / = /){
        if($errors[30] < $max_error_count){
          print $warnings
              "Illegal \'=\' character in <attributes> field on ",
              "line $line_num\n";
        }
        $errors[30]++;
        $attribute_string =~ s/=/ /g;
      }
      # check for tab character.
      if($attribute_string =~ /\t/){
        if($errors[9] < $max_error_count){
          print $warnings
              "Illegal tab character in <attributes> field on ",
              "line $line_num.\n";
        }
        $errors[9]++;
        $attribute_string =~ s/\t/ /g;
      }
      # make sure it has semicolon terminator after each attribute
      # otherwise fake it.
      unless($attribute_string =~ /;/){
        if($errors[10] < $max_error_count){
          print $warnings
              "<attributes> field conatins no semicolons",
              " terminating each attribute on line $line_num.\n";
        }
        $errors[10]++;
        my @att = split /\s+/, $attribute_string;
        $attribute_string = "";
        my $att_count = $#att - (($#att+1) % 2);
        for(my $i = 0;$i <= $att_count;$i+=2){
          $attribute_string .= $att[$i]." ".$att[$i+1]."; ";
        }
      }
      # parse <attributes> field
      my $attributes = {};
      while($attribute_string =~ /(\S+)(\s+)(\S+)(\s*)/g){
        my $name = $1;
        my $space1 = $2;
        my $value = $3;
        my $space2 = $4;
        if($value =~ /^\"(\S*)\";$/){
          $value = $1;
        }
        else{
          if($value =~ /^(\S+);$/){
            $value = $1;
          }
          else{
            if($errors[10] < $max_error_count){
              $value =~ /(.)$/;
              print $warnings
                  "Bad terminator, \"$1\", after ",
                  "<attributes> name-value pair, $name $value, on ",
                  "line $line_num. Should be \";\".\n";
            }
            $errors[10]++;
          }
          if($value =~ /^\"(\S+)\"$/){
            $value = $1;
          }
          else{
            if($errors[11] < $max_error_count){
              print $warnings
                  "<attributes> field missing quotes around",
                  " values on line $line_num.\n";
            }
            $errors[11]++;
          }
        }
        $attributes->{$name} = $value;
        unless($space1 =~ /^ $/){
          if($errors[12] < $max_error_count){
            print $warnings
                "Bad separator, $space1, between <attributes> ",
                "name-value pair, $name $value, on line ",
                "$line_num.\n";
          }
          $errors[12]++;
        }
        unless($space2 =~ /^ $/){
          if($errors[12] < $max_error_count){
            print $warnings
                "Bad separator, \"$space2\", between <attributes>",
                " name-value pairs on line $line_num.  Should be",
                " \" \".\n$name $value\n";
          }
          $errors[12]++;
        }		    
      }
      # check for gene_id and transcript_id attributes
      my $gid = "missing";
      if(defined($attributes->{gene_id})){
        $gid = $attributes->{gene_id};
      }
      else{
        if($errors[13] < $max_error_count){
          print $warnings
              "gene_id not set in <attributes> field on line ",
              "$line_num.\n";
        }
        $errors[13]++;
        $attributes->{gene_id} = "missing";		   
      }
      my $tid = "missing";
      if(defined($attributes->{transcript_id})){
        $tid = $attributes->{transcript_id};
      }		
      else{
        if($errors[14] < $max_error_count){
          print $warnings
              "transcript_id not set in <attributes> field on line ",
              "$line_num.\n";
        }
        $errors[14]++;
        $attributes->{transcript_id} = "missing";
      }
      #create objects
      my $feature = 
        GTF::Feature::new($data[2], $data[3], $data[4], $data[5], $data[7]);
      if($feature->type eq "CDS"){
        push @all_cds, $feature;
      }
      elsif($feature->type eq "inter") {
        push @all_inter, $feature;
      }
      elsif($feature->type eq "inter_CNS") {
        push @all_cns, $feature;
      }
      
      my $tx;
      my $gene;
      unless(defined($TxObj{$tid}) && length($tid)){
        $tx = GTF::Transcript::new($tid);
        $TxObj{$tid} = $tx;
        unless($feature->type eq "inter" || $feature->type eq "inter_CNS") {
          push @all_txs, $tx;
        }
        unless(defined($GeneObj{$gid}) && length($gid)){
          $gene = GTF::Gene::new($gid,$data[0],$data[1],$data[6]);
          $GeneObj{$gid} = $gene;
          unless($feature->type eq "inter" || $feature->type eq "inter_CNS") {
            push @all_genes, $gene;
          }
        }
        $gene = $GeneObj{$gid};
        $gene->add_transcript($tx);
      }
      else{
        $gene = $GeneObj{$gid};
      }
      $tx = $TxObj{$tid};
      #check that strand value is consistent with the strand of the tx 
      unless($data[6] eq $tx->strand){
        if($errors[20] < $max_error_count){
          print $warnings
            "Inconsistent <strand> value across gene_id = ".
                $tx->gene_id.".\n";
      }
      $errors[20]++;
      if (_check_errors_against_badlist($gtf,20)) {
        push @bad_list, $tx->id;
      }
      if($strict){
        die;
      }
        #last;???
      }
      $tx->add_feature($feature);
    }
  }
  close(GTF);
  #check each gene object
  foreach my $tx (@all_txs){
    my $exons = $tx->exons;
    my $cds = $tx->cds;
    my $utr3 = $tx->utr3;
    my $utr5 = $tx->utr5;
    my $sec = $tx->sec;
    my $inter = $tx->inter;
    my $inter_cns = $tx->inter_cns;
    if($#$exons >= 0 && $#$utr3 == -1 && $#$utr5 == -1){
      $tx->create_utr_objects_from_exons;
    }
    my $starts = $tx->start_codons;
    my $stops = $tx->stop_codons;
    my $strand = $tx->strand;
    #make sure this tx contains at least one CDS, inter_CNS, or inter feature
    if($#$cds == -1 && $#$inter_cns == -1 && $#$inter == -1){
      if($errors[15] < $max_error_count){
        print $warnings
          "Transcript ".$tx->id." contains no CDS, inter or inter_CNS features.\n";
      }
      $errors[15]++;
      if (_check_errors_against_badlist($gtf,15)) {
        push @bad_list, $tx->id;
      }
      next;
    }
    #verify start/stop codon exists and has length 3
    my $start_codon_start = -1;
    my $start_codon_stop = -1;
    if($#$starts == -1){
      if($errors[17] < $max_error_count){
        print $warnings
            "Missing start_codon for transcript \"".$tx->id."\".\n";
      }
      $errors[17]++;
      if (_check_errors_against_badlist($gtf,17)) {
        push @bad_list, $tx->id;
      }
    }
    else{	    
      my $length = 0;
      foreach my $start (@$starts){
        $length += $start->stop - $start->start + 1;
        if(($start_codon_start < 0) || ($start_codon_start > $start->start)){
          $start_codon_start = $start->start;
        }
        if(($start_codon_stop < 0) || ($start_codon_stop < $start->stop)){
          $start_codon_stop = $start->stop;
        }
      }
      if($length != 3){
        if($errors[16] < $max_error_count){
          print $warnings
              "Start Codon length is not three in transcript \"".$tx->id."\".\n";
        }
        $errors[16]++;
        if (_check_errors_against_badlist($gtf,16)) {
          push @bad_list, $tx->id;
        }
      }
    }
    my $stop_codon_start = -1;
    my $stop_codon_stop = -1;
    if($#$stops == -1){
      if($errors[19] < $max_error_count){
        print $warnings
            "Missing stop_codon for transcript \"".$tx->id."\".\n";
      }
      $errors[19]++;
      if (_check_errors_against_badlist($gtf,19)) {
        push @bad_list, $tx->id;
      }
    }
    else{
      my $length = 0;
      foreach my $stop (@$stops){
        $length += $stop->stop - $stop->start + 1;
        if(($stop_codon_start < 0) || ($stop_codon_start > $stop->stop)){
          $stop_codon_start = $stop->start;
        }
        if(($stop_codon_stop < 0) || ($stop_codon_stop < $stop->stop)){
          $stop_codon_stop = $stop->stop;
        }
      }
      if($length != 3){
        if($errors[18] < $max_error_count){
          print $warnings
              "Stop Codon length is not three in gene \"".$tx->id."\".\n";
        }
        $errors[18]++;
        if (_check_errors_against_badlist($gtf,18)) {
          push @bad_list, $tx->id;
        }
      }
    }

    #check that no cds or utr features overlap each other
    # also counts introns that are too short
    # introns shorter than 4 bp are in all likelihood 
    # biologically impossible (at least 2 bp for each splice signal)
    my $p_exons;
    my $overlap = 0;
    my $zero_intron_count = 0;
    my $short_intron_count = 0;
    foreach $p_exons ($utr5, $starts, $cds, $stops, $utr3)
    {
      my $last_stop;
      if ($#$p_exons > 0)
      {
        $last_stop = $$p_exons[0]->stop;
      }

      for(my $i = 1;$i <= $#$p_exons;$i++){
        my $intron_length = $$p_exons[$i]->start - $last_stop - 1;
        if($intron_length < 0){
          $overlap++;
        }
        elsif($intron_length < 4){
          $zero_intron_count++;
        }
        elsif($intron_length < 20){
          $short_intron_count++;
        }
        $last_stop = $$p_exons[$i]->stop;
      }
    }

    # check CDS - UTR junctions for overlaps
    # here, introns with length 0 are allowed because some
    # exons could be in part UTR and in part CDS
    if ($tx->strand eq "+")
    {
      my $ilen;
      if (($#$utr5 > -1) && ($#$starts > -1)){
        $ilen = $$starts[0]->start - $$utr5[$#$utr5]->stop - 1;
        $overlap++ if ($ilen < 0);
        $zero_intron_count++ if (($ilen > 0) && ($ilen < 4));
      }
      if (($#$cds > -1) && ($#$stops > -1)){
        # assumes stop codons are not included in CDS
        $ilen = $$stops[0]->start - $$cds[$#$cds]->stop - 1;
        $overlap++ if ($ilen < 0);
        $zero_intron_count++ if (($ilen > 0) && ($ilen < 4));
      }
      if (($#$stops > -1) && ($#$utr3 > -1)){
        $ilen = $$utr3[0]->start - $$stops[$#$stops]->stop - 1;
        $overlap++ if ($ilen < 0);
        $zero_intron_count++ if (($ilen > 0) && ($ilen < 4));
      }
    }
    else{
      my $ilen;
      if (($#$utr3 > -1) && ($#$stops > -1)){
        $ilen = $$stops[0]->start - $$utr3[$#$utr3]->stop - 1;
        $overlap++ if ($ilen < 0);
        $zero_intron_count++ if (($ilen > 0) && ($ilen < 4));
      }
      if (($#$stops > -1) && ($#$cds > -1)){
        # assumes stop codon is not included in CDS
        $ilen = $$cds[0]->start - $$stops[$#$stops]->stop - 1;
        $overlap++ if ($ilen < 0);
        $zero_intron_count++ if (($ilen > 0) && ($ilen < 4));
      }
      if (($#$starts > -1) && ($#$utr5 > -1)){
        $ilen = $$utr5[0]->start - $$starts[$#$starts]->stop - 1;
        $overlap++ if ($ilen < 0);
        $zero_intron_count++ if (($ilen > 0) && ($ilen < 4));
      }
    }

    if($overlap > 0){
      if($errors[45] < $max_error_count){
        print $warnings
            "Overlapping CDS or UTR features in transcript \"".$tx->id."\".\n";
      }
      $errors[45]++;
      if (_check_errors_against_badlist($gtf,45)) {
        push @bad_list, $tx->id;
      }
    }
    if($zero_intron_count > 0){
      if($errors[48] < $max_error_count){
        print $warnings "$zero_intron_count intron(s) shorter than 4 bp found on in transcript \"".$tx->id."\".\n" ;
      }
      $errors[48]++;
      if (_check_errors_against_badlist($gtf,48)) {
        push @bad_list, $tx->id;
      }
    }
    if($short_intron_count > 0){
      if($errors[49] < $max_error_count){
        print $warnings "$short_intron_count short (<20bp) intron(s) found on in transcript \"".$tx->id."\".\n" ;
      }
      $errors[49]++;
      if (_check_errors_against_badlist($gtf,49)) {
        push @bad_list, $tx->id;
      }
    }
    #check that each cds is between the start and stop and that the initial and 
    #terminal exons are placed correctly with respect to the start/stop codons
    if($start_codon_start >= 0){
      if($strand eq '-'){
        for(my $i = $#$cds;$i >= 0;$i--){
          if($$cds[$i]->end > $start_codon_stop){
            if($errors[21] < $max_error_count){
              print $warnings
                  "CDS before start codon in transcript_id = \"".
                      "".$tx->id."\".\n";
            }
            $errors[21]++;
            if (_check_errors_against_badlist($gtf,21)) {
              push @bad_list, $tx->id;
            }
          }
        }
        my $cds = $tx->initial_exon;
        if($cds->end != $start_codon_stop){
          if($errors[33] < $max_error_count){
            print $warnings
                "Start codon location not consistent with initial exon ".
                    "locaiton in transcript_id = \"".$tx->id."\".\n";
          }
          $errors[33]++;
          if (_check_errors_against_badlist($gtf,33)) {
            push @bad_list, $tx->id;
          }
        }
        if($cds->length < 3){
          my $cds = $tx->cds;
          my $error = 0;
          for(my $i = 0;$i <= $#$starts;$i++){
            if($i == $#$starts){
              unless($$cds[$#$cds-$i]->start <= $$starts[$#$starts-$i]->start){
                $error = 1;
              }
            }
            else{
              unless($$cds[$#$cds-$i]->start == $$starts[$#$starts-$i]->start){
                $error = 1;
              }
            }
          }
          if($error){
            if($errors[47] < $max_error_count){
              print $warnings "Start codon annotated in intron region on ".
                  "transcript_id = \"".$tx->id."\".\n";
            }
            $errors[47]++;
            if (_check_errors_against_badlist($gtf,47)) {
              push @bad_list, $tx->id;
            }
            my $pos = 0;
            my $start_pos = 0;
            my $cds_pos = $#$cds;
            #fix the start codon by splitting it across the cds
            while($pos < 3){
              my $end;
              if($$cds[$cds_pos]->length < 3-$pos){
                $end = $$cds[$cds_pos]->start;
                $pos += $$cds[$cds_pos]->length;
              }
              else{
                $end = $$cds[$cds_pos]->end-2+$pos;
                $pos = 3;
              }
              if($start_pos <= $#$starts){
                $$starts[$start_pos]->set_start($end);
                $$starts[$start_pos]->set_stop($$cds[$cds_pos]->end);
              }
              else{
                my $sc = GTF::Feature::new("start_codon", $end, $$cds[$cds_pos]->end, ".", ".");
                $tx->add_feature($sc);
              }
              $cds_pos--;
              $start_pos++;
            }
            $starts = [sort {$a->start <=> $b->start} @$starts];
          }
        }
      }
      else{
        for(my $i = 0;$i <= $#$cds;$i++){
          if($$cds[$i]->start < $start_codon_start){
            if($errors[21] < $max_error_count){
              print $warnings
                  "CDS before start codon in transcript_id = \"".
                      "".$tx->id."\".\n";
            }
            $errors[21]++;
            if (_check_errors_against_badlist($gtf,21)) {
              push @bad_list, $tx->id;
            }
          }
        }
        my $cds = $tx->initial_exon;
        if($cds->start != $start_codon_start){
          if($errors[33] < $max_error_count){
            print $warnings
                "Start codon location not consistent with initial exon ".
                    "locaiton in transcript_id = \"".$tx->id."\".\n";
          }
          $errors[33]++;
          if (_check_errors_against_badlist($gtf,33)) {
            push @bad_list, $tx->id;
          }
        }
        if($cds->length < 3){
          my $cds = $tx->cds;
          my $error = 0;
          for(my $i = 0; $i <= $#$starts;$i++){
            if($i == $#$starts){
              unless($$cds[$i]->stop >= $$starts[$i]->stop){
                $error = 1;
              }
            }
            else{
              unless($$cds[$i]->stop == $$starts[$i]->stop){
                $error = 1;
              }
            }
          }
          if($error){
            if($errors[47] < $max_error_count){
              print $warnings "Start codon annotated in intron region on ".
                  "transcript_id = \"".$tx->id."\".\n";
            }
            $errors[47]++;
            if (_check_errors_against_badlist($gtf,47)) {
              push @bad_list, $tx->id;
            }
            my $pos = 0;
            my $start_pos = 0;
            my $cds_pos = 0;
            #fix the start codon by splitting it across the cds
            while($pos < 3){
              my $end;
              if($$cds[$cds_pos]->length < 3-$pos){
                $end = $$cds[$cds_pos]->stop;
                $pos += $$cds[$cds_pos]->length;
              }
              else{
                $end = $$cds[$cds_pos]->start+2-$pos;
                $pos = 3;
              }
              if($start_pos <= $#$starts){
                $$starts[$start_pos]->set_start($$cds[$cds_pos]->start);
                $$starts[$start_pos]->set_stop($end);
              }
              else{
                my $sc = GTF::Feature::new("start_codon", $$cds[$cds_pos]->start, $end, ".", ".");
                $tx->add_feature($sc);
              }
              $cds_pos++;
              $start_pos++;
            }
            $starts = [sort {$a->start <=> $b->start} @$starts];
          }
        }
      }
    }
    if($stop_codon_start >= 0){
      if($strand =~ /\-/){
        for(my $i = 0;$i <= $#$cds;$i++){
          if($$cds[$i]->start < $stop_codon_start){
            if($errors[22] < $max_error_count){
              print $warnings
                  "CDS after stop codon in transcript_id = \"".
                      "".$tx->id."\".\n";
            }
            $errors[22]++;
            if (_check_errors_against_badlist($gtf,22)) {
              push @bad_list, $tx->id;
            }
          }
        }
        my $cds = $tx->terminal_exon;
        if($cds->start == $stop_codon_start){
          unless($cds->stop < $stop_codon_stop + 1){
            $cds->set_start($stop_codon_stop + 1);
          }
          if($errors[46] < $max_error_count){
            print $warnings
                "Stop codon included in terminal CDS transcript_id = ".
                    "\"".$tx->id."\".\n";
          }
          $errors[46]++;
          if (_check_errors_against_badlist($gtf,46)) {
            push @bad_list, $tx->id;
          }
        }
        if($cds->start != $stop_codon_stop + 1){
          if($errors[34] < $max_error_count){
            print $warnings
                "Stop codon location not consistent with terminal exon ".
                    "locaiton in transcript_id = \"".$tx->id."\".\n";
          }
          $errors[34]++;
          if (_check_errors_against_badlist($gtf,34)) {
            push @bad_list, $tx->id;
          }
          
        }
      }
      else{
        for(my $i = $#$cds;$i >= 0;$i--){
          if($$cds[$i]->end > $stop_codon_stop){
            if($errors[22] < $max_error_count){
              print $warnings
                  "CDS after stop codon in transcript_id = \"".$tx->id."\".\n";
            }
            $errors[22]++;
            if (_check_errors_against_badlist($gtf,22)) {
              push @bad_list, $tx->id;
            }
          }		    
        }
        my $cds = $tx->terminal_exon;
        if($cds->stop == $stop_codon_stop){
          unless($cds->start > $stop_codon_start -1){
            $cds->set_stop($stop_codon_start - 1);
          }
          if($errors[46] < $max_error_count){
            print $warnings
                "Stop codon included in terminal CDS transcript_id = ".
                    "\"".$tx->id."\".\n";
          }
          $errors[46]++;
          if (_check_errors_against_badlist($gtf,46)) {
            push @bad_list, $tx->id;
          }
        }
        if($cds->stop != $stop_codon_start - 1){
          if($errors[34] < $max_error_count){
            print $warnings
                "Stop codon location not consistent with terminal exon ".
                    "locaiton in transcript_id = \"".$tx->id."\".\n";
          }
          $errors[34]++;
          if (_check_errors_against_badlist($gtf,34)) {
            push @bad_list, $tx->id;
          }
        }
      }
    }
    #check that complete transcripts have even number of codons
    if(($start_codon_start >= 0) && ($stop_codon_stop >= 0)){
      my $len = 0;
      for(my $i = $#$cds;$i >= 0;$i--){
        $len += $$cds[$i]->length;
      }
      unless(($len % 3) == 0){
        if($errors[41] < $max_error_count){
          print $warnings
              "Complete transcript length not a multiple of three".
                  " on transcript_id = \"".$tx->id."\".\n";
        }
        $errors[41]++;
        if (_check_errors_against_badlist($gtf,41)) {
          push @bad_list, $tx->id;
        }
      }
    }
    #check for consistent phase
    my $phase = -1;
    if($start_codon_start >= 0){
      $phase = 0;
    }
    elsif($start_codon_stop >= 0){
      my $len = 0;
      for(my $i = $#$cds;$i >= 0;$i--){
        $len += $$cds[$i]->length;
      }
      $phase = $len % 3;
    }
    if($strand eq "-"){
      if($phase == -1){
        $phase = $$cds[$#$cds]->frame;
        if($phase eq "."){
          my $pos = -1;
          for(my $i = $#$cds; $i >= 0;$i--){
            unless($$cds[$i]->frame eq "."){
              $pos = $i;
              $phase = $$cds[$i]->frame;
              last;
            }
          }
          if($pos >= 0){
            for(my $i = $pos + 1; $i <= $#$cds;$i++){
              my $len = $$cds[$i]->length + $phase;
              $phase = $len % 3;
            }
          }
        }
        
      }
      unless($phase eq "."){
        my $bad_frame = 0;
        for(my $i = $#$cds;$i >= 0;$i--){
          unless($phase eq $$cds[$i]->frame){
            unless($$cds[$i]->frame eq "."){
              $bad_frame = 1;
            }
            $$cds[$i]->set_frame($phase);
          }
          $phase = ($$cds[$i]->length - $phase) % 3;
          if($phase == 1){
            $phase = 2;
          }
          elsif($phase == 2){
            $phase = 1;
          }
        }
        if($bad_frame){
          if($errors[42] < $max_error_count){
            print $warnings
                "Inconsistent frame accross transcript_id = ".
                    "\"".$tx->id."\".\n";
          }
          $errors[42]++;
          if (_check_errors_against_badlist($gtf,42)) {
            push @bad_list, $tx->id;
          }
        }
      }
    }
    else{
      if($phase == -1){
        $phase = $$cds[0]->frame;
        if($phase eq "."){
          my $pos = -1;
          for(my $i = 0; $i <= $#$cds;$i++){
            unless($$cds[$i]->frame eq "."){
              $pos = $i;
              $phase = $$cds[$i]->frame;
              last;
            }
          }
          if($pos >= 0){
            for(my $i = $pos - 1; $i >= 0;$i--){
              my $len = $$cds[$i]->length + $phase;
              $phase = $len % 3;
            }
          }
        }
      }
      unless($phase eq "."){
        my $bad_frame = 0;
        for(my $i = 0; $i <= $#$cds;$i++){
          unless($phase eq $$cds[$i]->frame){
            unless($$cds[$i]->frame eq "."){
              $bad_frame = 1;
            }
            $$cds[$i]->set_frame($phase);
          }
          $phase = ($$cds[$i]->length - $phase) % 3;
          if($phase == 1){
            $phase = 2;
          }
          elsif($phase == 2){
            $phase = 1;
          }
        }
        if($bad_frame){
          if($errors[42] < $max_error_count){
            print $warnings
                "Inconsistent frame accross transcript_id = ".
                    "\"".$tx->id."\".\n";
          }
          $errors[42]++;
          if (_check_errors_against_badlist($gtf,42)) {
            push @bad_list, $tx->id;
          }
        }
      }
    }
    # check selenocysteine (SEC) features here.  must be inside a 
    # CDS and have frame 0
    if($#$sec >= 0){
      foreach my $s (@$sec){
        if($s->frame ne "0" && $s->frame ne "."){
          #bad SEC frame
          if($errors[51] < $max_error_count){
            print $warnings
                "Bad frame, ".$s->frame."on SEC feature on tx ".$tx->id.".  Must of 0 or \".\"\n";
          }
          $errors[51]++;
          if (_check_errors_against_badlist($gtf,51)) {
            push @bad_list, $tx->id;
          }
        }
        my $good = 0;
        foreach my $c (@$cds){
          if($c->stop < $s->start){
            next;
          }
          if($c->start <= $s->start &&
             $c->stop >= $s->stop){
            $good = 1;
            last;
          }
          if($c->start > $s->stop){
            last;
          }
        }
        if($good == 0){
          # SEC outside CDS
          if($errors[52] < $max_error_count){
            print $warnings
                "SEC feature outside of CDS region on tx ".$tx->id.".\n";
          }
          $errors[52]++;
          if (_check_errors_against_badlist($gtf,52)) {
            push @bad_list, $tx->id;
          }
          
        }
      }
    }
  }

  #if the conservation sequence is available store it's info;
  my $conseqname = $gtf->{Conseq};
  my %total_conserv = (0 => 0,1 => 0,2 => 0);
  if($conseqname){
    open(Conseq, $conseqname)
        or die "Could not open conservation sequence file: \"$conseqname\".\n";
    my @conseq_info = stat(Conseq);
    my @comp;
    foreach my $tx (@all_txs){
      my $cds = $tx->cds;
      my $exons = $tx->exons;
      my $starts = $tx->start_codons;
      my $stops = $tx->stop_codons;	    
      foreach my $start (@$starts){
        push @comp,{feature => $start, "0" => 0, "1" => 0, "2" => 0};
      }
      foreach my $stop (@$stops){
        push @comp,{feature => $stop, "0" => 0, "1" => 0, "2" => 0};
      }
      foreach my $exon (@$cds){
        push @comp,{feature => $exon, "0" => 0, "1" => 0, "2" => 0};
      }
      foreach my $exon (@$exons){
        push @comp,{feature => $exon, "0" => 0, "1" => 0, "2" => 0};
      }
    }
    @comp = sort {$a->{feature}->start <=> $b->{feature}->start} @comp;
    my $base = "X";
    my $nuc_pos = 0;
    my $file_pos = 0;
    my $line = "X";
    my $next_comp = 0;
    my $max_comps = $#comp;
    my %checking;
    while(($next_comp <= $max_comps) and ($file_pos < $conseq_info[7])){
      $base = getc(Conseq);	    
      unless(defined($base)){
        print $warnings "Features beyond end of conservation sequence.\n";
        last;
      }
      $file_pos++;		   
      if($base =~ /\d/){
        $nuc_pos++;
        if(defined($total_conserv{$base})){
          $total_conserv{$base}++;
        }
        else{
          $total_conserv{X}++;
          if($errors[43] < $max_error_count){
            print $warnings "Bad conservation sequence character, $base.  ".
                "Should be \'0\',\'1\', or \'2\'.\n";
          }
          $errors[43]++;
        }
        #compile composition stats
        while(($next_comp <= $max_comps) and
              ($nuc_pos == $comp[$next_comp]->{feature}->start)){
          $checking{$comp[$next_comp]} = $comp[$next_comp];
          $next_comp++;
        }
        foreach my $id (keys %checking){
          $checking{$id}->{$base}++;
          if($checking{$id}->{feature}->stop == $nuc_pos){
            delete $checking{$id};
          }
        }
      }
    }
    while($file_pos < $conseq_info[7]){
      $base = getc(Conseq);	    
      $file_pos++;
      if($base =~ /\d/){
        $nuc_pos++;
        if(defined($total_conserv{$base})){
          $total_conserv{$base}++;
        }
        else{
          $total_conserv{X}++;
          if($errors[43] < $max_error_count){
            print $warnings "Bad conservation sequence character, $base.  ".
                "Should be \'0\',\'1\', or \'2\'.\n";
          }
          $errors[43]++;
        }
      }
    }
    close(Conseq);
    foreach my $info (@comp){
      $info->{feature}->set_conseq($info->{0},$info->{1},$info->{2});
    }
  }
  $gtf->{Total_Conseq} = \%total_conserv;
  #if the sequence is available use it to fix any problems
  my @inframe_stops;
  my $seqname = $gtf->{Sequence};
  my %tx_phase;
  foreach my $tx (@all_txs){
    $tx_phase{$tx->id} = -1;
  }
  my %total_seq = ('A' => 0,'C' => 0,'G' =>0,'T' => 0,'N' => 0,'X' => 0);
  #EVAN need to rewrite all this sequence loading stuff because it is fugly
  if($seqname){
    open(SEQ, "<$seqname")
        or die "Could not open sequence file: \"$seqname\".\n";
    my @seq_info = stat(SEQ);
    #try to fix things with the sequence info here
    #check start codons, stop codons, splice site, internal stop
    my %sequence;
    my %seq_pos;
    my @checks;
    foreach my $tx (@all_txs){
      $sequence{$tx->id} = "";
      $seq_pos{$tx->id} = 0;
      my $starts = $tx->start_codons;
      my $stops = $tx->stop_codons;
      my $utr3 = $tx->utr3;
      my $utr5 = $tx->utr5;
      my $cds = $tx->cds;
      if($#$cds == -1){
        next;
      }
      my $exons = $tx->exons;
      foreach my $start (@$starts){
        push @checks, {pos => $start->start,
                       type => "start_sequence",
                       tx => $tx,
                       object => $start,
                       stop => $start->stop};
      }
      foreach my $stop (@$stops){
        push @checks, {pos => $stop->start,
                       type => "stop_sequence",
                       tx => $tx,
                       object => $stop,
                       stop => $stop->stop};
      }
      foreach my $exon (@$exons){
        push @checks, {pos => $exon->start,
                       type => "exon_sequence",
                       tx => $tx,
                       object => $exon,
                       stop => $exon->stop};
      }
      foreach my $coding (@$cds){
        push @checks, {pos => $coding->start,
                       type => "coding_sequence",
                       tx => $tx,
                       object => $coding,
                       stop => $coding->stop};
      }
    
      my($p_exons, $splice_type1, $splice_type2);
      if ($tx->strand eq "+"){
        $splice_type1 = "donor_splice";
        $splice_type2 = "acceptor_splice";
      }
      else{
        $splice_type1 = "acceptor_splice";
        $splice_type2 = "donor_splice";
      }

      foreach $p_exons ($utr5, $starts, $cds, $stops, $utr3){
        next if ($#$p_exons < 1);

        my($i, $subtype);
        $subtype = ($$p_exons[0]->type =~ /UTR/) ? "utr" : "cds";
        for ($i = 1; $i <= $#$p_exons; $i++){
          push @checks, {pos  => $$p_exons[$i-1]->stop+2, 
                         type => $subtype."_".$splice_type1, 
                         tx => $tx,
                         stop => 0};
          push @checks, {pos  => $$p_exons[$i]->start-1, 
                         type => $subtype."_".$splice_type2, 
                         tx => $tx,
                         stop => 0};
        }
      }

      # this assumes that the innermost start codon is included in CDS
      # therefore, no checks are run on start_codon - CDS junction
      # also, this allows zero-length introns at UTR5 - start_codon 
      # CDS - stop_codon and UTR3 - stop_codon junctions 
      if ($tx->strand eq "+"){
        if (($#$utr5 > -1) && ($#$starts > -1)){
          my $ilen = $$starts[0]->start - $$utr5[$#$utr5]->stop - 1;
          if ($ilen != 0){
            push @checks, {pos  => $$utr5[$#$utr5]->stop+2, 
                           type => "utr_donor_splice", 
                           tx => $tx,
                           stop => 0};
            push @checks, {pos  => $$starts[0]->start-1, 
                           type => "utr_acceptor_splice", 
                           tx => $tx,
                           stop => 0};
          }
        }
        if (($#$cds > -1) && ($#$stops > -1)){
          my $ilen = $$stops[0]->start - $$cds[$#$cds]->stop - 1;
          if ($ilen != 0){
            push @checks, {pos  => $$cds[$#$cds]->stop+2, 
                           type => "cds_donor_splice", 
                           tx => $tx,
                           stop => 0};
            push @checks, {pos  => $$stops[0]->start-1, 
                           type => "cds_acceptor_splice", 
                           tx => $tx,
                           stop => 0};
          }
        }  
        if (($#$stops > -1) && ($#$utr3 > -1)){
          my $ilen = $$utr3[0]->start - $$stops[$#$stops]->stop - 1;
          if ($ilen != 0){
            push @checks, {pos  => $$stops[$#$stops]->stop+2, 
                           type => "utr_donor_splice", 
                           tx => $tx,
                           stop => 0};
            push @checks, {pos  => $$utr3[0]->start-1, 
                           type => "utr_acceptor_splice", 
                           tx => $tx,
                           stop => 0};
          }
        }
        if ($#{$stops} < 0){
          push @checks, {pos => $$cds[$#{$cds}]->stop+1,
                         type => "stop_sequence?",
                         tx => $tx,
                         stop => $$cds[$#{$cds}]->stop+3};
        }
        if ($#{$starts} < 0){
          if ($$cds[0]->start-3 > 0){
            push @checks, {pos => $$cds[0]->start-3,
                           type => "start_sequence?",
                           tx => $tx,
                           stop => $$cds[0]->start-1};
          }
          else{
            $sequence{$tx->id} = "NNN";
          }
        }
      }	    
      else{
        if (($#$utr3 > -1) && ($#$stops > -1)){
          my $ilen = $$stops[0]->start - $$utr3[$#$utr3]->stop - 1;
          if ($ilen != 0){
            push @checks, {pos  => $$utr3[$#$utr3]->stop+2, 
                           type => "utr_acceptor_splice", 
                           tx => $tx,
                           stop => 0};
            push @checks, {pos  => $$stops[0]->start-1, 
                           type => "utr_donor_splice", 
                           tx => $tx,
                           stop => 0};
          }
        }
        if (($#$stops > -1) && ($#$cds > -1)){
          my $ilen = $$cds[0]->start - $$stops[$#$stops]->stop - 1;
          if ($ilen != 0){
            push @checks, {pos  => $$stops[$#$stops]->stop+2, 
                           type => "cds_acceptor_splice", 
                           tx => $tx,
                           stop => 0};
            push @checks, {pos  => $$cds[0]->start-1, 
                           type => "cds_donor_splice", 
                           tx => $tx,
                           stop => 0};
          }
        }  
        if (($#$starts > -1) && ($#$utr5 > -1)){
          my $ilen = $$utr5[0]->start - $$starts[$#$starts]->stop - 1;
          if ($ilen != 0){
            push @checks, {pos  => $$starts[$#$starts]->stop+2, 
                           type => "utr_acceptor_splice", 
                           tx => $tx,
                           stop => 0};
            push @checks, {pos  => $$utr5[0]->start-1, 
                           type => "utr_donor_splice", 
                           tx => $tx,
                           stop => 0};
          }
        }
        if ($#{$stops} < 0){
          if ($$cds[0]->start-3 > 0){
            push @checks, {pos => $$cds[0]->start-3,
                           type => "stop_sequence?",
                           tx => $tx,
                           stop => $$cds[0]->start-1};
          }
        }
        if ($#{$starts} < 0){
          push @checks, {pos => $$cds[$#{$cds}]->stop+1,
                         type => "start_sequence?",
                         tx => $tx,
                         stop => $$cds[$#{$cds}]->stop+3};
        }
      }	    
    }
    
    @checks  = sort{$a->{pos} <=> $b->{pos} || $a->{stop} <=> $b->{stop}} @checks;
    my $file_pos = length(<SEQ>);
    my $base = "X";
    my $nuc_pos = 0;
    my $line = "X";
    my $next = 0;
    my $prev_base = "X";
    my $max_checks = $#checks;
    my %checking;
    my $splice_type;
    
    while((($next <= $max_checks) || ((keys %checking) >= 0)) 
          && ($file_pos < $seq_info[7])){
      $base = getc(SEQ);
      unless(defined($base)){
        print $warnings "Features beyond end of sequence.\n";
        last;
      }
      $file_pos++;		   
      if($base =~ /\S/){
        $base = uc($base);
        if((defined($total_seq{$base})) && ($base ne 'X')){
          $total_seq{$base}++;
        }
        else{
          $total_seq{X}++;
          if($errors[44] < $max_error_count){
            print $warnings "Bad fasta sequence character, $base. Should be".
                " \'A\',\'C\',\'T\',\'G\',or \'N\'.\n";
          }
          $errors[44]++;
        }
        $nuc_pos++;
        while(($next <= $max_checks) && ($checks[$next]{pos} == $nuc_pos)){
          if($checks[$next]{type} =~ /sequence/){
            $checking{$checks[$next]{tx}->id . "\t" . 
                          $checks[$next]{type}} = {
                            check => $checks[$next],
                            A => 0,
                            C => 0,
                            G => 0,
                            T => 0,
                            N => 0};
          }
          elsif(($checks[$next]{type} =~ /^(\w+)_acceptor_splice$/)
                && ($splice_type = $1)) {
            if($checks[$next]{tx}->strand eq "-"){
              unless(($prev_base eq "C") && ($base eq "T")){
                if (($prev_base eq "G") && ($base eq "T")) {
                  if ($errors[50] < $max_error_count) {
                    print $warnings
                        "Possible non-canonical splice site sequence, ".
                            "$prev_base$base, on negative strand transcript \"".
                                $checks[$next]{tx}->id.
                                    "\".\n";
                  }
                  $errors[50]++;
                  if ((_check_errors_against_badlist($gtf,50)) 
                      && ($splice_type ne "utr")) {
                    push @bad_list, $checks[$next]{tx}->id;
                  }
                }
                else {
                  if($errors[26] < $max_error_count){
                    print $warnings
                        "Bad acceptor splice site sequence, $prev_base".
                            "$base, on negative strand transcript \"".
                                $checks[$next]{tx}->id.
                                    "\".  Should be CT or GT.\n";
                  }
                  $errors[26]++;
                  if (_check_errors_against_badlist($gtf,26)) {
                    push @bad_list, $checks[$next]{tx}->id;
                  }
                }
              }
              
            }
            else{
              unless(($prev_base eq "A") && ($base eq "G")){
                if (($prev_base eq "A") && ($base eq "C")) {
                  if ($errors[50] < $max_error_count) {
                    print $warnings
                        "Possible non-canonical splice site sequence, ".
                            "$prev_base$base, on transcript \"".
                                $checks[$next]{tx}->id.
                                    "\".\n";
                  }
                  $errors[50]++;
                  if ((_check_errors_against_badlist($gtf,50)) 
                      && ($splice_type ne "utr")) {
                    push @bad_list, $checks[$next]{tx}->id;
                  }
                }
                else {
                  if($errors[26] < $max_error_count){
                    print $warnings
                        "Bad acceptor splice site sequence, $prev_base".
                            "$base, on transcript \"".
                                $checks[$next]{tx}->id.
                                    "\".  Should be AG or AC.\n";
                  }
                  $errors[26]++;
                  if (_check_errors_against_badlist($gtf,26)) {
                    push @bad_list, $checks[$next]{tx}->id;
                  }
                }
              }
            }
          }
          elsif(($checks[$next]{type} =~ /^(\w+)_donor_splice$/)
                && ($splice_type = $1)) {
            if($checks[$next]{tx}->strand eq "-"){
              unless(($prev_base eq "A") && ($base eq "C")){
                if (($prev_base eq "A") && ($base eq "T")) {
                  if ($errors[50] < $max_error_count) {
                    print $warnings
                        "Possible non-canonical splice site sequence, ".
                            "$prev_base$base, on negative strand transcript \"".
                                $checks[$next]{tx}->id.
                                    "\".\n";
                  }
                  $errors[50]++;
                  if ((_check_errors_against_badlist($gtf,50)) 
                      && ($splice_type ne "utr")) {
                    push @bad_list, $checks[$next]{tx}->id;
                  }
                }
                elsif (($prev_base eq "G") && ($base eq "C")) {
                  if ($errors[50] < $max_error_count) {
                    print $warnings
                        "Possible non-canonical splice site sequence, ".
                            "$prev_base$base, on negative strand transcript \"".
                                $checks[$next]{tx}->id.
                                    "\".\n";
                  }
                  $errors[50]++;
                  if ((_check_errors_against_badlist($gtf,50)) 
                      && ($splice_type ne "utr")) {
                    push @bad_list, $checks[$next]{tx}->id;
                  }
                }					
                else {
                  if ($errors[27] < $max_error_count){
                    print $warnings
                        "Bad donor splice site sequence, $prev_base".
                            "$base, on negative strand transcript \"".
                                $checks[$next]{tx}->id.
                                    "\".  Should be AC, GC, or AT.\n";
                  }
                  $errors[27]++;
                  if (_check_errors_against_badlist($gtf,27)) {
                    push @bad_list, $checks[$next]{tx}->id;
                  }
                }
              }
            }
            else{
              unless(($prev_base eq "G") && ($base eq "T")){
                if (($prev_base eq "G") && ($base eq "C")) {
                  if ($errors[50] < $max_error_count) {
                    print $warnings
                        "Possible non-canonical splice site sequence, ".
                            "$prev_base$base, on transcript \"".
                                $checks[$next]{tx}->id.
                                    "\".\n";
                  }
                  $errors[50]++;
                  if ((_check_errors_against_badlist($gtf,50)) 
                      && ($splice_type ne "utr")) {
                    push @bad_list, $checks[$next]{tx}->id;
                  }
                }
                elsif (($prev_base eq "A") && ($base eq "T")) {
                  if ($errors[50] < $max_error_count) {
                    print $warnings
                        "Possible non-canonical splice site sequence, ".
                            "$prev_base$base, on transcript \"".
                                $checks[$next]{tx}->id.
                                    "\".\n";
                  }
                  $errors[50]++;
                  if ((_check_errors_against_badlist($gtf,50)) 
                      && ($splice_type ne "utr")) {
                    push @bad_list, $checks[$next]{tx}->id;
                  }
                }					
                else {
                  if ($errors[27] < $max_error_count){
                    print $warnings
                        "Bad donor splice site sequence, $prev_base".
                            "$base, on transcript \"".
                                $checks[$next]{tx}->id.
                                    "\".  Should be GT, GC, or AT.\n";
                  }
                  $errors[27]++;
                  if (_check_errors_against_badlist($gtf,27)) {
                    push @bad_list, $checks[$next]{tx}->id;
                  }
                }
              }
            }
          }
          $next++;
        }
        foreach my $check (keys %checking){		
          $check =~ /^(.+)\t.+$/;
          my $tx_id = $1;
          unless(($checking{$check}{check}{type} =~ /exon/) || 
                 ($seq_pos{$tx_id} >= $nuc_pos)){
            $sequence{$tx_id}.= uc($base);
            $seq_pos{$tx_id} = $nuc_pos;
          }
          $checking{$check}{uc($base)}++;
          if($checking{$check}{check}{stop} eq $nuc_pos){			
            my $object = $checking{$check}{check}{object};
            if(defined($object)){
              if($object->strand eq "-"){
                $object->set_bases($checking{$check}{T},
                                   $checking{$check}{G},
                                   $checking{$check}{C},
                                   $checking{$check}{A},
                                   $checking{$check}{N});
                
              }
              else{
                $object->set_bases($checking{$check}{A},
                                   $checking{$check}{C},
                                   $checking{$check}{G},
                                   $checking{$check}{T},
                                   $checking{$check}{N});
              }
            }
            delete $checking{$check};
          }
        }	
        $prev_base = $base;
      }
    }
    foreach my $tx (@all_txs){
      my $cds = $tx->cds;
      unless($#$cds >= 0){
        next;
      }
      my $sec = $tx->sec;
      my $starts = $tx->start_codons;
      my $stops = $tx->stop_codons;
      my $strand = $tx->strand;
      my $tx_id = $tx->id;
      my $phase = $$cds[0]->frame;
      my $old_phase = $phase;
      if($strand eq "-"){
        $phase = $$cds[$#$cds]->frame;
      }
      my $seq = $sequence{$tx->id};	    
      if($seq =~ /^$/){
        print STDERR "No sequence for $tx_id\n";
        next;
      }
      if($tx->strand eq "-"){
        $seq = _rev_comp($seq);
      }
      my $start = 0;
      my $stop = 0;
      my $found_start = -1;
      my $found_stop = -1;	    
      my @seq = split //, $seq;	    
      #check some SEC stuff here;
      foreach my $s (@$sec){
        my $pos = 0;
        if($tx->strand eq '-'){
          foreach my $c (reverse @$cds){
            if($c->start > $s->stop){
              $pos += $c->stop-$c->start+1;
            }
            else{
              $pos += $c->stop-$s->stop;
              last;
            }
          }
          if($#$stops == -1){
            $pos += 3;
          }
        }
        else{
          foreach my $c (@$cds){
            if($c->stop < $s->start){
              $pos += $c->stop-$c->start+1;
            }
            else{
              $pos += $s->start-$c->start;
              last;
            }
          }
          if($#$starts == -1){
            $pos += 3;
          }
        }
        my $codon = substr($seq,$pos,3);
        unless($codon eq "TGA"){
          #bad SEC codon
          $error_msgs[53] = "SEC feature had bad sequence.  Must be \"TGA\"";
          if($errors[53] < $max_error_count){
            print $warnings "SEC feature had bad sequence, $codon, on tx, ".$tx->id.
                ".  Must be \"TGA\"\n";
          }
          $errors[53]++;
          if (_check_errors_against_badlist($gtf,53)) {
            push @bad_list, $tx->id;
          }
        }
        substr($seq,$pos,3) = "NNN";
      }
      if($#$starts >= 0){
        my $sc = substr($seq,0,3);
        unless($sc eq "ATG"){
          if($errors[24] < $max_error_count){
            print $warnings
                "Bad start codon sequence, $sc, on tx ".
                    $tx->id.".  Should be ATG.\n";
          }
          $errors[24]++;
          if (_check_errors_against_badlist($gtf,24)) {
            push @bad_list, $tx->id;
          }
        }
      }
      else{
        my $sc1 = substr($seq,0,3);
        my $sc2 = substr($seq,3,3);
        if($sc2 eq "ATG"){
          $found_start = 3;
          $start = 3;
        }
        elsif($sc1 eq "ATG"){
          $found_start = 0;
        }
        else{
          $start = 3;
        }
      }
      if($#$stops >= 0){
        my $sc = substr($seq,$#seq-2,3);
        unless(($sc eq "TAA") ||
               ($sc eq "TAG") ||
               ($sc eq "TGA")){
          if($errors[25] < $max_error_count){
            print $warnings
                "Bad stop codon sequence, $sc, on tx ".
                    $tx->id.".  Should be TAA, TAG, or TGA.\n";
          }
          $errors[25]++;
          if (_check_errors_against_badlist($gtf,25)) {
            push @bad_list, $tx->id;
          }
        }		       
      }
      else{
        my $sc1 = substr($seq,$#seq - 2,3);
        my $sc2 = substr($seq,$#seq - 5,3);
        if(($sc2 eq "TAA") ||
           ($sc2 eq "TAG") ||
           ($sc2 eq "TGA")){
          $found_stop = $#seq - 5;
          $stop = 3;
        }
        elsif(($sc1 eq "TAA") ||
              ($sc1 eq "TAG") ||
              ($sc1 eq "TGA")){
          $found_stop = $#seq - 2;
        }
        else{
          $stop = 3;
        }
      }
      my @sc = (0,0,0);
      my $scphase = $phase;
      $scphase = 0 if($scphase eq '.');
      for(my $i = $start;$i <= $#seq - $stop - 5;$i++){
        my $codon = substr($seq,$i,3);
        if(($codon eq "TAA") ||
           ($codon eq "TAG") ||
           ($codon eq "TGA")){
          $sc[($i+$scphase) % 3] = 1;
        }
      }
      if(($phase ne ".") &&
         ($sc[0])){
        if($errors[29] < $max_error_count){
          print $warnings "In frame stop codon(s) found on transcript ".
              $tx->id."\n";
        }
        $errors[29]++;
        if (_check_errors_against_badlist($gtf,29)) {
          push @bad_list, $tx->id;
        }
        push @inframe_stops, $tx->id;
      }
      else{
        if(($sc[0] == 1) &&
           ($sc[1] == 1) &&
           ($sc[2] == 1)){
          push @inframe_stops, $tx->id;
          if($errors[37] < $max_error_count){
            print $warnings "Stop Codons in all frames on transcript ".
                $tx->id."\n";
          }
          $errors[37]++;
          if (_check_errors_against_badlist($gtf,37)) {
            push @bad_list, $tx->id;
          }
        }
      }
      # add start or stop if possible
      if($found_stop != -1){
        if(($phase eq ".") ||
           ($phase == (($#seq + 1) % 3))){
          if($sc[($#seq +1) % 3] == 0){
            $phase = ($#seq + 1) % 3;
            if($errors[32] < $max_error_count){
              print $warnings
                  "Unannotated stop codon found on transcript ".
                      $tx->id."\n";
            }
            $errors[32]++;
            if (_check_errors_against_badlist($gtf,32)) {
              push @bad_list, $tx->id;
            }
            if($fix_gtf){
              if($strand eq "+"){
                my $terminal_exon = $$cds[$#$cds];
                my $sc_start = $terminal_exon->stop + 1 - $stop;
                my $stop_codon = 
                  GTF::Feature::new("stop_codon",$sc_start,
                                    $sc_start + 2,0,0);
                $terminal_exon->set_stop($sc_start - 1);
                $tx->add_feature($stop_codon);
              }
              else{
                my $terminal_exon = $$cds[0];
                my $sc_start = $terminal_exon->start - 3 + $stop;
                my $stop_codon = 
                  GTF::Feature::new("stop_codon",$sc_start,
                                    $sc_start + 2,0,0);
                $terminal_exon->set_start($sc_start + 3);
                $tx->add_feature($stop_codon);
              }
            }
          }
        }
      }
      if($found_start != -1){
        if(($phase eq ".") ||
           ($phase == 0)){
          if($sc[0] == 0){
            if($errors[31] < $max_error_count){
              print $warnings
                  "Unannotated start codon found on transcript ".
                      $tx->id."\n";
            }
            $errors[31]++;
            if (_check_errors_against_badlist($gtf,31)) {
              push @bad_list, $tx->id;
            }
            if($fix_gtf){
              if($strand eq "+"){
                my $initial_exon = $$cds[0];
                my $sc_start = $initial_exon->start - 3 + $start;
                my $start_codon = 
                  GTF::Feature::new("start_codon",$sc_start,
                                    $sc_start + 2,0,0);
                $initial_exon->set_start($sc_start);
                $tx->add_feature($start_codon);
              }
              else{
                my $initial_exon = $$cds[$#$cds];
                my $sc_start = $initial_exon->stop + 1 - $start;
                my $start_codon = 
                  GTF::Feature::new("start_codon",$sc_start,
                                    $sc_start + 2,0,0);
                $initial_exon->set_stop($sc_start + 2);
                $tx->add_feature($start_codon);
              }
            }
          }
        }
      }
      # set phase if unset here;
      if(($phase ne ".") && ($old_phase eq ".") && ($fix_gtf)){
        if($strand eq "+"){
          for(my $i = 0;$i <= $#$cds;$i++){
            $$cds[$i]->set_frame($phase);
            $phase = ($$cds[$i]->length - $phase) % 3;
            if($phase == 1){
              $phase = 2;
            }
            elsif($phase == 2){
              $phase = 1;
            }
          }
        }
        else{
          for(my $i = $#$cds;$i >= 0;$i--){
            $$cds[$i]->set_frame($phase);
            $phase = ($$cds[$i]->length - $phase) % 3;
            if($phase == 1){
              $phase = 2;
            }
            elsif($phase == 2){
              $phase = 1;
            }
          }
        }
      }
      if($gtf->{Tx}){
        my $tx_fh = $gtf->{Tx};
        print $tx_fh ">".$tx->id ."\t".$tx->start."\t".$tx->stop."\n".
            substr($seq,$start,$#seq - $stop - $start + 1)."\n";
      }
    }
  }
  $gtf->{Total_Seq} = \%total_seq;
  for(my $i = 0;$i <= $#error_msgs;$i++){
    if(($errors[$i] > 0) && !($$no_warn[$i])){
      print $warnings 
          "\nWarnings encountered:\n",
          "Count\tDescription\n";
      for(my $j = $i;$j <= $#error_msgs;$j++){
        if(($errors[$j] > 0) && !($$no_warn[$j])){
          print $warnings
              $errors[$j],"\t$error_msgs[$j]\n";
        }
      }
      last;
    }
  }
  
  #sort genes, transcripts, cds and utr by start position
  @all_genes  = sort {$a->start <=> $b->start} @all_genes;
  @all_txs  = sort {$a->start <=> $b->start} @all_txs;
  @all_cds  = sort {$a->start <=> $b->start} @all_cds;
  @all_inter = sort {$a->start <=> $b->start} @all_inter;
  @all_cns = sort {$a->start <=> $b->start} @all_cns;

  print $warnings "\nStatistics:\n";
  print $warnings 
      "\t", $#all_genes + 1, " genes with ", $#all_txs + 1,
      " transcripts containing ", $#all_cds + 1 , " cds.\n";
  if($gtf->{Inframe_Stops}){
    print $warnings "Inframe Stop Genes:\n";
    foreach my $stop (@inframe_stops){
      print $warnings "$stop\n";
    }
  }
  if ($gtf->{Bad_List}) {
    print $warnings "Bad Genes:\n";
    my %bad_out;
    foreach my $bad (@bad_list) {
      $bad_out{$bad}++;
    }
    my @out = keys %bad_out; 
    while (@out) {
      print $warnings pop(@out),"\n";
    }
  }
  close(WARN);
  close(GTF);
  $gtf->{Genes}        = \@all_genes;
  $gtf->{Transcripts}  = \@all_txs;
  $gtf->{CDS}          = \@all_cds;
  $gtf->{Inter}        = \@all_inter;
  $gtf->{Inter_CNS}    = \@all_cns;
  $gtf->{Comments}     = $all_line_comments;
  $gtf->{Parse_Errors} = \@errors;
}

# GTF::_get_error_messages()
#    This function returns an array of the error messages associated with
#    each error number.
sub _get_error_messages{
  my @error_msgs;
  $error_msgs[0] = "Not enough fields on line.";
  $error_msgs[1] = "Illegal type of whitespace between fields.  Sould be tab.";
  $error_msgs[2] = "Illegal value for <feature> field.  Should be ".
      "\'CDS\', \'exon\', \'start_codon\', or \'stop_codon\'.";
  $error_msgs[3] = "Illegal value for <start> field.  Should be numerical.";
  $error_msgs[4] = "Illegal value for <stop> field.  Should be numerical.";
  $error_msgs[5] = "Illegal value for <score> field.  Should be numerical or \'.\'.";
  $error_msgs[6] = "Illegal value for <strand> field.  ".
      "Should be \'+\', \'-\', or \'.\'.";
  $error_msgs[7] = "Illegal value for <frame> field.  ".
      "Should be \'0\', \'1\', \'2\', or \'.\'.";
  $error_msgs[8] = "Illegal \'\#\' character in attribute name.";
  $error_msgs[9] = "Illegal tab character in <attributes> field.";
  $error_msgs[10] = "Missing \';\' terminator after attribute in <attributes> field.";
  $error_msgs[11] = "Missing \'\"\'s around attribute values in  <attributes> field.";
  $error_msgs[12] = "Incorrect separator between attribute name-value pairs in ".
      "<attributes> field.  Should be \' \'.";
  $error_msgs[13] = "Missing \'gene_id\' attribute in <attributes> field.";
  $error_msgs[14] = "Missing \'transcript_id\' attribute in <attributes> field.";
  $error_msgs[15] = "Transcript contains no CDS, inter or inter_CNS features.";
  $error_msgs[16] = "Start Codon length is not three.";
  $error_msgs[17] = "Transcript has no start codon.";
  $error_msgs[18] = "Stop Codon length is not three.";
  $error_msgs[19] = "Transcript has no stop codon.";	    
  $error_msgs[20] = "Inconsistent <strand> value across gene_id.";
  $error_msgs[21] = "CDS before start codon.";
  $error_msgs[22] = "CDS after stop codon.";
  $error_msgs[23] = "Multiple genes using the same transcript_id.";
  $error_msgs[24] = "Bad start codon sequence.  Should be ATG.";
  $error_msgs[25] = "Bad stop codon sequence.  Should be TAA, TAG, or TGA.";
  $error_msgs[26] = "Bad acceptor splice site sequence.";
  $error_msgs[27] = "Bad donor splice site sequence.";
  $error_msgs[28] = "In frame stop codon found.";
  $error_msgs[29] = "In frame stop codon in transcript.";
  $error_msgs[30] = "Illegal \'=\' character in <attributes> field.";
  $error_msgs[31] = "Unannotated start codon found.";
  $error_msgs[32] = "Unannotated stop codon found.";
  $error_msgs[33] = "Start codon location inconsistent with initial exon location.";
  $error_msgs[34] = "Stop codon location inconsistent with terminal exon location.";
  $error_msgs[35] = "Start codon length is wrong.";
  $error_msgs[36] = "Start codon length is wrong.";
  $error_msgs[37] = "Transcript does not translate in any frame.";
  $error_msgs[38] = "<start> field is greater than <stop> field.";
  $error_msgs[39] = "Wrong phase value.";
  $error_msgs[40] = "Missing phase value.";
  $error_msgs[41] = "Complete transcript length not multiple of 3.";
  $error_msgs[42] = "Inconsistant phase value accross transcript.";
  $error_msgs[43] = "Bad conservation sequence character.  Should be \'0\',\'1\', ".
      "or \'2\'.";
  $error_msgs[44] = "Bad fasta sequence character. Should be \'A\',\'C\',\'T\',\'G\',".
      "or \'N\'.";
  $error_msgs[45] = "Overlapping CDS features in transcript.";
  $error_msgs[46] = "Stop codon included in terminal CDS.";
  $error_msgs[47] = "Start codon annotated in intron region.";
  $error_msgs[48] = "Transcript contains zero length intron(s).";
  $error_msgs[49] = "Transcript contains short (<20bp) intron(s).";
  $error_msgs[50] = "Possible non-canonical splice site sequence.";
  $error_msgs[51] = "Bad frame on SEC feature.  Must be 0 or \".\"";
  $error_msgs[52] = "SEC feature outside of CDS feature.";
  $error_msgs[53] = "SEC feature had bad sequence.  Must be \"TGA\"";
  return @error_msgs;
}

#this checks to see if the error is one of the errors that
#we want to return the transcript name and remove
sub _check_errors_against_badlist{
  my ($gtf,$error_num) = @_;
  my $bad = $gtf->{Bad_List};
  if($$bad[$error_num]) {
    return 1;
  }
  else {
    return 0;
  }
}


#this loads the gtf file without checking for errors
sub _load_no_check{
  my ($gtf,$warnings) = @_;
  my (@all_genes,@all_txs,@all_cds,@all_inter,@all_cns);
  my %TxObj;
  my %GeneObj;
  my $strict = $gtf->{Strict};
  #open file 
  open(GTF, "<$gtf->{Filename}") 
      or die "Could not open $gtf->{Filename}.\n";
  while(my $in_line = <GTF>){
    chomp $in_line;
    if($in_line =~ /^(.*)\#/){
      $in_line = $1;
    }
    if($in_line =~ /^\s+(.*)$/){
      $in_line = $1;
    }
    if($in_line !~ /\S/){
      next;
    }
    my @data = split /\t/,$in_line;
    if($#data != 8){
      die "Wrong number of fields on line:\n$in_line\n";
    }
    if(!$strict && $data[5] eq "."){
      $data[5] = 0;
    }
    my %attributes;
    while($data[8] =~ /^\s*(\S+) \"(\S*)\";(.*)$/){
      my $key = $1;
      my $val = $2;
      $data[8] = $3;
      $attributes{$key} = $val;
    }
    my $gid;
    if(defined($attributes{gene_id})){
      $gid = $attributes{gene_id};
    }
    else{
      die "Missing gene_id\n";
    }
    my $tid;
    if(defined($attributes{transcript_id})){
      $tid = $attributes{transcript_id};
    }
    else{
      die "Missing transcript_id\n";
    }
    #create objects
    my $feature = 
      GTF::Feature::new($data[2], $data[3], $data[4], $data[5], $data[7]);
    if($feature->type eq "CDS"){
      push @all_cds, $feature;
    }
    elsif($feature->type eq "inter") {
      push @all_inter, $feature;
    }
    elsif($feature->type eq "inter_CNS") {
      push @all_cns, $feature;
    }

    my $tx;
    my $gene;
    unless(defined($TxObj{$tid}) && length($tid)){
      $tx = GTF::Transcript::new($tid);
      $TxObj{$tid} = $tx;
      unless($feature->type eq "inter" || $feature->type eq "inter_CNS") {
        push @all_txs, $tx;
      }
      unless(defined($GeneObj{$gid}) && length($gid)){
        $gene = GTF::Gene::new($gid,$data[0],$data[1],$data[6]);
        $GeneObj{$gid} = $gene;
        unless($feature->type eq "inter" || $feature->type eq "inter_CNS") {
          push @all_genes, $gene;
        }
      }
      $gene = $GeneObj{$gid};
      $gene->add_transcript($tx);
    }
    else{
      $gene = $GeneObj{$gid};
    }
    $tx = $TxObj{$tid};
    if($tx->strand ne $data[6]){
      die "Features on transcript ".$tx->id." have different strands\n";
    }
    $tx->add_feature($feature);
  }
  close(GTF);
  #sort genes, exons, and cds by start position
  @all_genes  = sort {$a->start <=> $b->start} @all_genes;
  @all_txs  = sort {$a->start <=> $b->start} @all_txs;
  @all_cds  = sort {$a->start <=> $b->start} @all_cds;
  @all_inter = sort {$a->start <=> $b->start} @all_inter;
  @all_cns = sort {$a->start <=> $b->start} @all_cns;
  $gtf->{Genes}        = \@all_genes;
  $gtf->{Transcripts}  = \@all_txs;
  $gtf->{CDS}          = \@all_cds;
  $gtf->{Inter_CNS}    = \@all_cns;
  $gtf->{Inter}        = \@all_inter;
}

sub _rev_comp{
  my ($seq) = @_;
  $seq = uc $seq;
  my @bases = split //, $seq;
  my $out = "";
  foreach my $base (@bases){
    if($base eq 'A'){
      $base = 'T';
    }
    elsif($base eq 'C'){
      $base = 'G';
    }
    elsif($base eq 'G'){
      $base = 'C';
    }
    elsif($base eq 'T'){
      $base = 'A';
    }			
    $out = $base.$out;
  }
  return $out;
  
}

###############################################################################
# GTF::Gene
###############################################################################
# The GTF::Gene object stores all gtf data for one gene.  A gene contains one
# or more transcripts
package GTF::Gene;
# GTF::Gene::new(gene_id,seqname,source,strand)
#    This is the contrcutor for a GTF::Gene object.  It has 4 required 
#    parameters.  The are: 
#      gene_id    - The gene_id for this gene
#      seqname   - The <seqname> field of the line.  Should be a string 
#                   containing no whitespace.
#      source     - The <source> field of the line.  Should be a string 
#                   containing no whitespace.
#      strand     - The <strand> field of the line.  Should be either
#                   '+', '-', or '.'.
sub new {
  my $gene = bless {};
  ($gene->{Id},$gene->{Seqname}, $gene->{Source}, $gene->{Strand}) = @_;
  $gene->{Transcripts} = [];
  $gene->{CDS} = [];
  $gene->{Modified} = 0;
  $gene->{Start} = -1;
  $gene->{Stop} = -1;
  return $gene;
}

sub add_transcript{
  my ($gene,$tx) = @_;
  my $txs = $gene->{Transcripts};
  push @$txs, $tx;
  $tx->_set_gene($gene);
  $tx->_update;
  $gene->{Modified} = 1;
}

sub id {shift->{Id}}

sub set_id{
  my ($gene,$id) = @_;
  $gene->{Id} = $id;
}

# GTF::Feature::seqname()
#    This function returns the sequnce name, <seqname> field, that this 
#    feature came from as a string.
sub seqname    {shift->{Seqname}}

# GTF::Feature::set_source()
#    This function sets the sequence name, <seqname> field, of this line
sub set_seqname{
  my ($feature,$seqname) = @_;
  $feature->{Seqname} = $seqname;
}

# GTF::Feature::source()
#    This function returns the source, <source> field, of this line as a 
#    string.
sub source      {shift->{Source}}

# GTF::Feature::set_source()
#    This function sets the source, <source> field, of this line
sub set_source{
  my ($feature,$source) = @_;
  $feature->{Source} = $source;
}

# GTF::Feature::strand()
#    This function returns the strand, <strand> field, of this line as a 
#   string.  Should be '+', '-', or '.'.
sub strand      {shift->{Strand}}

# GTF::Gene::id()
#    This function returns this transcript's unique transcript_id
sub gene_id{
  my ($gene) = @_;
  return $gene->{Id};
}

# GTF::Gene::transcript_ids()
#    This function returns an array of all trnacsripts in this gene
sub transcripts{
  my ($gene) = @_;
  $gene->_update;
  return $gene->{Transcripts};# [sort {$a->start <=> $b->start || $a->stop <=> $b->stop} @{$gene->{Transcripts}}];
}

# GTF::Gene::remove_transcript()
sub remove_transcript{
  my($gene, $tx_id) = @_;
  my $txs = $gene->transcripts;
  my @new_txs = ();

  @new_txs = grep($_->id ne $tx_id, @$txs);
  $gene->{Transcripts} = \@new_txs if ($#new_txs < $#$txs);
}

sub cds{
  my ($gene) = @_;

  my $txs = $gene->transcripts;
  my @cds;
  foreach my $tx (@$txs){
    push @cds, @{$tx->cds}; 
  }

  # employ a slow sort to avoid perl sort problems (see GTF::Gene::_update) rpz
  my $tmp;
  for(my $a = 0; $a <= $#cds; $a++){
    for(my $b = $a+1; $b <= $#cds; $b++){
      if(($cds[$a]->start > $cds[$b]->start) ||
         (($cds[$a]->start == $cds[$b]->start) &&
          ($cds[$a]->stop > $cds[$b]->stop))){
        $tmp = $cds[$a];
        $cds[$a] = $cds[$b];
        $cds[$b] = $tmp;
      }
    }
  }

  return \@cds;
}

# GTF::Gene::copy()
#    This function returns a new object which is a copy of this object.  Also
#    creates copies of all transcripts of this gene
sub copy{
  my ($gene) = @_;
  my $copy = GTF::Gene::new($gene->id,$gene->seqname,$gene->source,$gene->strand);
  my $txs = $gene->transcripts;
  my @new_txs;
  my $next = 0;
  foreach my $tx (@$txs){	
    $next++;
    my $ntx = $tx->copy;
    $copy->add_transcript($ntx);
  }
  return($copy);
}

# GTF::Gene::offset(ammount)
#    This function will offset all locations in this gene by the given ammount
sub offset{
  my ($gene,$offset) = @_;
  my $txs = $gene->transcripts;
  foreach my $tx (@$txs){
    $tx->offset($offset);
  }
}

# GTF::Gene::output_gtf([file_handle])
#    This function outputs the information for this gene in standard 
#    gtf2 format.  If takes and optional file_handle as the only argument 
#    to which it outputs the data.  If no argument is given it outputs 
#    the data to stdout.
sub output_gtf{
  my ($gene,$out_handle) = @_;
  unless(defined($out_handle)){
    $out_handle = \*STDOUT;
  }
  my $txs = $gene->transcripts;
  foreach my $tx (@$txs){
    $tx->output_gtf($out_handle);
    print $out_handle "\n";
  }
}

# GTF::Gene::output_gff([file_handle])
#    This function outputs the information for this gene in standard 
#    GFF format.  If takes and optional file_handle as the only argument 
#    to which it outputs the data.  If no argument is given it outputs 
#    the data to stdout.
sub output_gff{
  my ($gene,$out_handle) = @_;
  unless(defined($out_handle)){
    $out_handle = \*STDOUT;
  }
  my $txs = $gene->transcripts;
  foreach my $tx (@$txs){
    $tx->output_gff($out_handle);
  }
}

# GTF::Gene::reverse_complement(seq_length)
#    This function takes the length of the sequence this gene came from and 
#    reverse complements everything in the gene.  
sub reverse_complement{
  my ($gene, $seq_length) = @_;
  my $txs = $gene->transcripts;
  foreach my $tx (@$txs){
    $tx->reverse_complement($seq_length);
  }
  $gene->{Modified} = 1;
}

# GTF::Gene::start()
#    This function return the lowest coord of any feature in any tx of this gene
sub start{
  my ($gene) = @_;
  $gene->_update;
  return $gene->{Start};
}

# GTF::Gene::stop()
#    This function returns the highest coord of any feature in any tx of this gene
sub stop{
  my ($gene) = @_;
  $gene->_update;
  return $gene->{Stop};
}

sub equals{
  my ($gene,$compare) = @_;
  $gene->_update;
  unless($gene && $compare){
    return 0;
  }
  unless(($gene->start == $compare->start) &&
         ($gene->stop == $compare->stop) &&
         ($gene->strand eq $compare->strand)){
    return 0;
  }
  my $gene_txs = $gene->transcripts;
  my $compare_txs = $compare->transcripts;
  unless($#$gene_txs == $#$compare_txs){
    return 0;
  }
  for(my $i = 0;$i <= $#$gene_txs;$i++){
    unless($$gene_txs[$i]->equals($$compare_txs[$i])){
      return 0;
    }
  }
  return 1;
}

sub length{
  my ($gene) = @_;
  my $txs = $gene->transcripts;
  if($#$txs == -1){
    return 0;
  }
  my $start = $$txs[0]->start;
  my $stop = $$txs[0]->stop;
  my $total = 0;
  foreach my $tx (@$txs){
    if($tx->start < $start){
      $start = $tx->start;
    }
    if($tx->stop > $stop){
      $stop = $tx->stop;
    }
  }
  return($stop - $start + 1);
}

sub gc_percentage{
  my ($gene) = @_;
  my $txs = $gene->transcripts;
  if($#$txs == -1){
    return 0;
  }
  my $total = 0;
  foreach my $tx (@$txs){
    $total += $tx->gc_percentage;
  }
  return($total/($#$txs+1));
}

sub match_percentage{
  my ($gene) = @_;
  my $txs = $gene->transcripts;
  if($#$txs == -1){
    return 0;
  }
  my $total = 0;
  foreach my $tx (@$txs){
    $total += $tx->match_percentage;
  }
  return($total/($#$txs+1));
}

sub mismatch_percentage{
  my ($gene) = @_;
  my $txs = $gene->transcripts;
  if($#$txs == -1){
    return 0;
  }
  my $total = 0;
  foreach my $tx (@$txs){
    $total += $tx->mismatch_percentage;
  }
  return($total/($#$txs+1));
}

sub unaligned_percentage{
  my ($gene) = @_;
  my $txs = $gene->transcripts;
  if($#$txs == -1){
    return 0;
  }
  my $total = 0;
  foreach my $tx (@$txs){
    $total += $tx->unaligned_percentage;
  }
  return($total/($#$txs+1));
}

sub tag{
  my ($gene) = @_;
  return $gene->{Tag};
}

sub set_tag{
  my ($gene,$tag) = @_;
  $gene->{Tag} = $tag;
}

sub _update{
  my ($gene) = @_;
  unless($gene->{Modified}){
    return;
  }
  $gene->{Modified} = 0;
  my $txs = $gene->{Transcripts};


  # this slow sort is being done because the commented out sort below
  # causes crashes sometimes and I have no idea why since the uncommented
  # sort has no problems
  # -- the reason for the problems is that perl's sort routine does not
  #    afford you the capability of modifying the $a or $b variable in 
  #    the usersub callback, since they are aliases.  rpz
  $gene->{Start} = -1;
  $gene->{Stop} = -1;
  my $tmp;
  for(my $a = 0; $a <= $#$txs; $a++){
    for(my $b = $a+1; $b <= $#$txs; $b++){
      if(($$txs[$a]->start > $$txs[$b]->start) ||
         (($$txs[$a]->start == $$txs[$b]->start) &&
          ($$txs[$a]->stop > $$txs[$b]->stop))){
        $tmp = $$txs[$a];
        $$txs[$a] = $$txs[$b];
        $$txs[$b] = $tmp;
      }
    }
    if($$txs[$a]->stop > $gene->{Stop}){
      # highest stop may not be stop of last tx so we check them all
      $gene->{Stop} = $$txs[$a]->stop;
    }
  }
  if($#$txs != -1){
    # lowest start is start of first tx  
    $gene->{Start} = $$txs[0]->start;
  }
}

###############################################################################
# GTF::Transcript
###############################################################################
# The GTF::Transcript object stores all data for a single transcript of a gtf
# file.  
package GTF::Transcript;
# GTF::Transcript::new(transcript_id)
#    This is the constructor for GTF::Feature objects.  It takes 1 arguments.
#      transciprt_id - the id for this transcript
sub new{
  my $tx = bless {};
  ($tx->{Id}) = @_;
  $tx->{Exons} = [];
  $tx->{CDS} = [];
  $tx->{UTR3} = [];
  $tx->{UTR5} = [];
  $tx->{Introns} = [];
  $tx->{Starts} = [];
  $tx->{Stops} = [];
  $tx->{SEC} = [];
  $tx->{Intron_CNS} = [];
  $tx->{Inter_CNS} = [];
  $tx->{Inter} = [];
  $tx->{Coding_Start} = -1;
  $tx->{Coding_Stop} = -1;
  $tx->{Coding_Length} = -1;
  $tx->{Gene} = -1;
  $tx->{Start} = -1;
  $tx->{Stop} = -1;
  $tx->{Modified} = 1;
  return $tx;
}

sub _set_gene{
  my ($tx,$gene) = @_;
  $tx->{Gene} = $gene;
}

sub add_feature{
  my ($tx,$feature) = @_;
  my $type = $feature->type;
  if($type eq 'CDS'){
    push @{$tx->{CDS}}, $feature; 
    $tx->{CDS} = [sort {$a->start <=> $b->start} @{$tx->{CDS}}];
  }
  elsif($type eq '3UTR'){
    push @{$tx->{UTR3}}, $feature; 
    $tx->{UTR3} = [sort {$a->start <=> $b->start} @{$tx->{UTR3}}];
  }
  elsif($type eq '5UTR'){
    push @{$tx->{UTR5}}, $feature; 
    $tx->{UTR5} = [sort {$a->start <=> $b->start} @{$tx->{UTR5}}];
  }
  elsif($type eq 'start_codon'){
    push @{$tx->{Starts}}, $feature; 
    $tx->{Starts} = [sort {$a->start <=> $b->start} @{$tx->{Starts}}];
  }
  elsif($type eq 'stop_codon'){
    push @{$tx->{Stops}}, $feature; 
    $tx->{Stops} = [sort {$a->start <=> $b->start} @{$tx->{Stops}}];
  }
  elsif($type eq 'exon'){
    push @{$tx->{Exons}}, $feature; 
    $tx->{Exons} = [sort {$a->start <=> $b->start} @{$tx->{Exons}}];
  }
  elsif($type eq 'SEC'){
    push @{$tx->{SEC}}, $feature; 
    $tx->{SEC} = [sort {$a->start <=> $b->start} @{$tx->{SEC}}];
  }
  elsif($type eq 'intron_CNS'){
    push @{$tx->{Intron_CNS}}, $feature;
    $tx->{Intron_CNS} = [sort {$a->start <=> $b->start} @{$tx->{Intron_CNS}}];
  }
  elsif($type eq 'inter_CNS'){
    push @{$tx->{Inter_CNS}}, $feature;
    $tx->{Inter_CNS} = [sort {$a->start <=> $b->start} @{$tx->{Inter_CNS}}];
  }
  elsif($type eq 'inter'){
    push @{$tx->{Inter}}, $feature;
    $tx->{Inter} = [sort {$a->start <=> $b->start} @{$tx->{Inter}}];
  }
  $feature->set_transcript($tx);
  $tx->{Modified} = 1;
  if($tx->{Gene} != -1){ $tx->{Gene}->{Modified} = 1;}
}

sub transfer_missing_utr{
  my($dest, $src, $mode) = @_;
  my $dest_utr3 = $dest->utr3;
  my $dest_utr5 = $dest->utr5;
  my $src_utr3 = $src->utr3;
  my $src_utr5 = $src->utr5;
  
  $mode = "Full" if ($#_ == 1); # default

  if (($mode eq "3UTR") || ($mode eq "Full")) {
    if (($#$dest_utr3 == -1) && ($#$src_utr3 > -1)){
      print "Transfering 3UTR from ".$src->id." to ".$dest->id."\n";
      foreach my $src_exon (@$src_utr3)
      {
        my $utr_exon = $src_exon->copy;
        $utr_exon->set_transcript($dest);
        $dest->add_feature($utr_exon);
      }
    }
  }
  if (($mode eq "5UTR") || ($mode eq "Full")) {  
    if (($#$dest_utr5 == -1) && ($#$src_utr5 > -1)){
      print "Transfering 5UTR from ".$src->id." to ".$dest->id."\n";
      foreach my $src_exon (@$src_utr5)
      {
        my $utr_exon = $src_exon->copy;
        $utr_exon->set_transcript($dest);
        $dest->add_feature($utr_exon);
      }
    }
  }
}

# copies the exon obejcts info into utr5 and utr3 objects
# WARNING: don't call this if you already have utr5 and utr3 info
sub create_utr_objects_from_exons{
  my ($tx) = @_;
  my $exons = $tx->exons;
  my $start = $tx->coding_start;
  my $stop = $tx->coding_stop;
  if($start != -1){
    my $type = "5UTR";
    if($tx->strand eq "-"){
      $type = "3UTR";
    }
    for(my $i = 0; $i <= $#$exons; $i++){
      if($$exons[$i]->start < $start){
        my $utr = $$exons[$i]->copy;
        $utr->set_type($type);
        if($start < $utr->stop){
          $utr->set_stop($start);
        }
        $tx->add_feature($utr);
      }
      else{
        last;
      }
    }
  }
  if($stop != -1){
    my $type = "3UTR";
    if($tx->strand eq "-"){
      $type = "5UTR";
    }
    for(my $i = $#$exons; $i >= 0; $i--){
      if($$exons[$i]->stop > $stop){
        my $utr = $$exons[$i]->copy;
        $utr->set_type($type);
        if($stop > $utr->start){
          $utr->set_start($stop);
        }
        $tx->add_feature($utr);
      }
      else{
        last;
      }
    }    
  }
}

sub exons{
  my ($tx) = @_;
  $tx->_update;
  return $tx->{Exons};
}

sub cds{
  my ($tx) = @_;
  $tx->_update;
  return $tx->{CDS};
}

sub utr3{
  my ($tx) = @_;
  $tx->_update;
  return $tx->{UTR3};
}

sub utr5{
  my ($tx) = @_;
  $tx->_update;
  return $tx->{UTR5};
}

sub introns{
  my ($tx) = @_;
  $tx->_update;
  return $tx->{Introns};
}

sub start_codons{
  my ($tx) = @_;
  $tx->_update;
  return $tx->{Starts};
}

sub stop_codons{
  my ($tx) = @_;
  $tx->_update;
  return $tx->{Stops};
}

sub sec{
  my ($tx) = @_;
  $tx->_update;
  return $tx->{SEC};
}

sub intron_cns{
  my ($tx) = @_;
  $tx->_update;
  return $tx->{Intron_CNS};
}

sub inter_cns{
  my ($tx) = @_;
  $tx->_update;
  return $tx->{Inter_CNS};
}

sub inter{
  my ($tx) = @_;
  $tx->_update;
  return $tx->{Inter};
}

sub seqname{
  my ($tx) = @_;
  my $gene = $tx->gene;
  return $gene->seqname;
}

sub source{
  my ($tx) = @_;
  my $gene = $tx->gene;
  return $gene->source;
}

sub strand{
  my ($tx) = @_;
  my $gene = $tx->gene;
  return $gene->strand;
}

sub set_id{
  my ($tx,$id) = @_;
  $tx->{Id} = $id;
}

sub id{
  my ($tx) = @_;
  return $tx->{Id};
}

sub gene_id{
  my ($tx) = @_;
  my $gene = $tx->gene;
  return $gene->{Id};
}

sub gene{
  my ($tx) = @_;
  return $tx->{Gene};
}

# GTF::Transcript::remove_exon(exon start pos)
#    Allows you to remove an exon from the list
#    Removes the all exons at the specified position
sub remove_exon{
  my ($tx, $pos) = @_;
  my $new_exons = [];
  my $removed = 0;
  foreach my $exon (@{$tx->{Exons}}){
    if($exon->start != $pos){
      push @$new_exons, $exon;
    }
    else{
      $removed++;
    }
  }
  $tx->{Exons} = $new_exons;
  return $removed;
}

# GTF::Transcript::remove_cds(cds start pos)
#    Allows you to remove a cds from the list
#    Removes the all cdss at the specified position
sub remove_cds{
  my ($tx, $pos) = @_;
  my $new_cds = [];
  my $removed = 0;
  foreach my $cds (@{$tx->{CDS}}){
    if($cds->start != $pos){
      push @$new_cds, $cds;
    }
    else{
      $removed++;
    }
  }
  $tx->{CDS} = $new_cds;
  return $removed;
}

# GTF::Transcript::length()
#    This function return the length of the tx from the start to the stop.
sub length{
  my ($tx) = @_;
  return($tx->stop - $tx->start + 1);
}

# GTF::Transcript::start()
#    This function returns the start of this tx (lowest coord).  
sub start{    
  my ($tx) = @_;
  $tx->_update;
  return $tx->{Start};
}

# GTF::Transcript::stop()
#    This function return the end of the tx (highest coord). 
sub stop{
  my ($tx) = @_;
  $tx->_update;
  return $tx->{Stop};
}

# GTF::Transcript::coding_start()
#    This function returns the start position (lowest coord) of coding region of 
#    this tx.  Same as the start codon's <start> field.
sub coding_start{
  my ($tx) = @_;
  $tx->_update;
  return $tx->{Coding_Start};
}

# GTF::Transcript::coding_stop()
#    This function return the end position (highest coord) of the coding region of 
#    the tx.  Same as the stop codon's <stop> field.
sub coding_stop{
  my ($tx) = @_;
  $tx->_update;
  return $tx->{Coding_Stop};
}

# GTF::Transcript::coding_length()
#    This function return the coding_length of the tx from the coding_start 
#    to the coding_stop.
sub coding_length{
  my ($tx) = @_;
  $tx->_update;
  return $tx->{Coding_Length};
}

sub score{
  my ($tx) = @_;
  my $features = $tx->_all_features;
  my $score = 0;
  foreach my $feature (@$features){
    if($feature->score ne "."){
      $score += $feature->score;
    }
  }
  return $score;
}

# GTF::Transcript::initial_exon()
#    This function returns a CDS::Feature object representing the  
#    initial_exon (cds) of this tx
#    NOTE:  Isn't exon the wrong word for this?
sub initial_exon{
  my ($tx) = @_;
  my $starts = $tx->start_codons;
  if($#$starts == -1){
    return 0;
  }
  $tx->_update;
  my $exons = $tx->{CDS};
  if($tx->strand eq "-"){
    return $$exons[$#{$exons}];
  }
  else{
    return $$exons[0];
  }
}

# GTF::Transcript::terminal_exon()
#    This function returns a CDS::Feature object representing the  
#    terminal exon (cds) of this tx
#    NOTE:  Isn't exon the wrong word for this?
sub terminal_exon{
  my ($tx) = @_;
  my $stops = $tx->stop_codons;
  if($#$stops == -1){
    return 0;
  }
  $tx->_update;
  my $exons = $tx->cds;
  if($tx->strand eq "-"){
    return $$exons[0];
  }
  else{
    return $$exons[$#{$exons}];
  }
}

sub equals{
  my ($tx,$compare, $mode) = @_;
  $mode = "Full" if ($#_ == 1); # default

  unless($tx && $compare){
    return 0;
  }
  unless(($tx->start == $compare->start) &&
         ($tx->stop == $compare->stop) &&
         ($tx->strand eq $compare->strand)){
    return 0;
  }

  if (($mode eq "CDS") || ($mode eq "Full")) {
    my $tx_cds = $tx->cds;
    my $comp_cds = $compare->cds;
    unless($#$tx_cds == $#$comp_cds){
      return 0;
    }
    for(my $i = 0;$i <= $#$tx_cds;$i++){
      unless($$tx_cds[$i]->equals($$comp_cds[$i])){
        return 0;
      }
    }
    my $tx_starts = $tx->start_codons;
    my $comp_starts = $compare->start_codons;
    unless($#$tx_starts == $#$comp_starts){
      return 0;
    }
    for(my $i = 0;$i <= $#$tx_starts;$i++){
      unless($$tx_starts[$i]->equals($$comp_starts[$i])){
        return 0;
      }
    }
    my $tx_stops = $tx->stop_codons;
    my $comp_stops = $compare->stop_codons;
    unless($#$tx_stops == $#$comp_stops){
      return 0;
    }
    for(my $i = 0;$i <= $#$tx_stops;$i++){
      unless($$tx_stops[$i]->equals($$comp_stops[$i])){
        return 0;
      }
    }
  }

  if (($mode eq "3UTR") || ($mode eq "Full")) {
    my $tx_utr3 = $tx->utr3;
    my $comp_utr3 = $compare->utr3;
    unless($#$tx_utr3 == $#$comp_utr3){
      return 0;
    }
    for(my $i = 0;$i <= $#$tx_utr3;$i++){
      unless($$tx_utr3[$i]->equals($$comp_utr3[$i])){
        return 0;
      }
    }
  }

  if (($mode eq "5UTR") || ($mode eq "Full")) {
    my $tx_utr5 = $tx->utr5;
    my $comp_utr5 = $compare->utr5;
    unless($#$tx_utr5 == $#$comp_utr5){
      return 0;
    }
    for(my $i = 0;$i <= $#$tx_utr5;$i++){
      unless($$tx_utr5[$i]->equals($$comp_utr5[$i])){
        return 0;
      }
    }
  }

  my $tx_exon = $tx->exons;
  my $comp_exon = $compare->exons;
  unless($#$tx_exon == $#$comp_exon){
    return 0;
  }
  for(my $i = 0;$i <= $#$tx_exon;$i++){
    unless($$tx_exon[$i]->equals($$comp_exon[$i])){
      return 0;
    }
  }

  my $tx_sec = $tx->sec;
  my $comp_sec = $compare->sec;
  unless($#$tx_sec == $#$comp_sec){
    return 0;
  }
  for(my $i = 0;$i <= $#$tx_sec;$i++){
    unless($$tx_sec[$i]->equals($$comp_sec[$i])){
      return 0;
    }
  }
  return 1;
}

# GTF::Transcript::offset(ammount)
#    This function will offset all locations in this tx by the given ammount
sub offset{
  my ($tx,$offset) = @_;
  unless(defined($offset)){
    print STDERR "Undefined value passed to GTF::Gene::offset.\n";
    return;
  }
  my $features = $tx->_all_features;
  foreach my $feature (@$features){
    if(defined($feature)){
      $feature->offset($offset);
    }
  }
  $tx->{Modified} = 1;
  if($tx->{Gene} != -1){ $tx->{Gene}->{Modified} = 1;}
}

# GTF::Transcript::output_gtf([file_handle])
#    This function outputs the information for this tx in standard 
#    gtf2 format.  If takes and optional file_handle as the only argument 
#    to which it outputs the data.  If no argument is given it outputs 
#    the data to stdout.
sub output_gtf{
  my ($tx,$out_handle) = @_;
  unless(defined $out_handle){
    $out_handle = \*STDOUT;
  }
  my @features = 
      sort {$a->start <=> $b->start || $a->stop <=> $b->stop} @{$tx->_all_features};
  foreach my $feature (@features){
    $feature->output_gtf($out_handle);
  }
}

# GTF::Transcript::output_gff([file_handle])
#    This function outputs the information for this tx in standard 
#    GFF format.  If takes and optional file_handle as the only argument 
#    to which it outputs the data.  If no argument is given it outputs 
#    the data to stdout.
sub output_gff{
  my ($tx,$out_handle) = @_;
  unless(defined $out_handle){
    $out_handle = \*STDOUT;
  }
  my @features = 
      sort {$a->start <=> $b->start || $a->stop <=> $b->stop} @{$tx->_all_features};
  foreach my $feature (@features){
    $feature->output_gff($out_handle);
  }
}

# GTF::Transcript::reverse_complement(seq_length)
#    This function takes the length of the sequence this tx came from and 
#    reverse complements everything in the gene.  
sub reverse_complement{
  my ($tx, $seq_length) = @_;
  unless(defined($seq_length)){
    print STDERR 
        "Undefined value for seq_length passed to GTF::Gene::reverse_complement.\n";
  }
  my $features = $tx->_all_features;
  foreach my $feature (@$features){
    $feature->reverse_complement($seq_length);
  }
}

sub _all_features{
  my ($tx) = @_;
  my @features;
  my $starts = $tx->start_codons;    
  push @features, @$starts;
  my $stops = $tx->stop_codons;
  push @features, @$stops;
  my $exons = $tx->exons;
  push @features, @$exons;
  my $cds = $tx->cds;
  push @features, @$cds;
  my $utr3 = $tx->utr3;
  push @features, @$utr3;
  my $utr5 = $tx->utr5;
  push @features, @$utr5;
  my $sec = $tx->sec;
  push @features, @$sec;
  my $cns = $tx->intron_cns;
  push @features, @$cns;
  return \@features;
}

sub copy{
  my ($tx) = @_;
  my $copy = GTF::Transcript::new($tx->id);
  my $starts = $tx->start_codons;
  foreach my $old (@$starts){
    my $new = $old->copy;
    $copy->add_feature($new);
  }
  my $stops = $tx->stop_codons;
  foreach my $old (@$stops){
    my $new = $old->copy;
    $copy->add_feature($new);
  }
  my $exons = $tx->exons;
  foreach my $old (@$exons){
    my $new = $old->copy;
    $copy->add_feature($new);
  }
  my $cds = $tx->cds;
  foreach my $old (@$cds){
    my $new = $old->copy;
    $copy->add_feature($new);
  }
  my $utr3 = $tx->utr3;
  foreach my $old (@$utr3){
    my $new = $old->copy;
    $copy->add_feature($new);
  }
  my $utr5 = $tx->utr5;
  foreach my $old (@$utr5){
    my $new = $old->copy;
    $copy->add_feature($new);
  }
  my $sec = $tx->sec;
  foreach my $old (@$sec){
    my $new = $old->copy;
    $copy->add_feature($new);
  }
  return($copy);
}

sub gc_percentage{
  my ($tx) = @_;
  my $features = $tx->_all_features;
  my $total = 0;
  my $percent = 0;
  foreach my $feature (@$features){
    $percent += $feature->length * $feature->gc_percentage;
    $total += $feature->length;
  }
  $percent /= $total;
  return $percent;
}

sub match_percentage{
  my ($tx) = @_;
  my $features = $tx->_all_features;
  my $total = 0;
  my $percent = 0;
  foreach my $feature (@$features){
    $percent += $feature->length * $feature->match_percentage;
    $total += $feature->length;
  }
  $percent /= $total;
  return $percent;
}

sub mismatch_percentage{
  my ($tx) = @_;
  my $features = $tx->_all_features;
  my $total = 0;
  my $percent = 0;
  foreach my $feature (@$features){
    $percent += $feature->length * $feature->mismatch_percentage;
    $total += $feature->length;
  }
  $percent /= $total;
  return $percent;
}

sub unaligned_percentage{
  my ($tx) = @_;
  my $features = $tx->_all_features;
  my $total = 0;
  my $percent = 0;
  foreach my $feature (@$features){
    $percent += $feature->length * $feature->unaligned_percentage;
    $total += $feature->length;
  }
  $percent /= $total;
  return $percent;
}

sub tag{
  my ($tx) = @_;
  return $tx->{Tag};
}

sub set_tag{
  my ($tx,$tag) = @_;
  $tx->{Tag} = $tag;
}

sub _update{
  my ($tx) = @_;
  unless($tx->{Modified}){
    return;
  }
  $tx->{Modified} = 0;
  my $starts = $tx->start_codons;    
  my $stops = $tx->stop_codons;
  my $cds = $tx->cds;
  my $utr5 = $tx->utr5;
  my $utr3 = $tx->utr3;
  my $exons = $tx->exons;
  foreach my $utr (@$utr5){
    $utr->set_subtype("UTR5");
  }
  foreach my $utr (@$utr3){
    $utr->set_subtype("UTR3");
  }
  #set subtypes
  if($#$cds == 0){
    if(($#$starts >= 0) && ($#$stops >= 0)){
      $$cds[0]->set_subtype("Single");
    }
    elsif($#$starts >= 0){
      $$cds[0]->set_subtype("Initial");
    }
    elsif($#$stops >= 0){
      $$cds[0]->set_subtype("Terminal");
    }
    else{
      $$cds[0]->set_subtype("Internal");	    
    }
  }
  elsif($#$cds >= 0){
    my $strand = $tx->strand;
    my $starts = $tx->start_codons;
    my $stops = $tx->stop_codons;
    my $start = 0;
    my $stop = $#$cds;
    if($strand eq '-'){
      if($#$starts >= 0){
        $$cds[$#$cds]->set_subtype("Initial");
        $stop--;
      }
      if($#$stops >= 0){
        $$cds[0]->set_subtype("Terminal");
        $start++;
      }
    }
    else{
      if($#$starts >= 0){
        $$cds[0]->set_subtype("Initial");
        $start++;
      }
      if($#$stops >= 0){
        $$cds[$#$cds]->set_subtype("Terminal");
        $stop--;
      }
    }
    for(my $i = $start;$i <= $stop;$i++){
      $$cds[$i]->set_subtype("Internal");
    }
  }
  #get introns and coding length
  my $clen = 0;
  my @introns;
  for(my $i = 0;$i < $#$cds;$i++){
    $clen += $$cds[$i]->stop - $$cds[$i]->start +1;
    my $start = $$cds[$i]->stop + 1;
    my $stop = $$cds[$i+1]->start - 1;
    if($start <= $stop){
      my $intron = GTF::Feature::new("Intron",$start,$stop,0,'.');
      $intron->set_subtype("Intron");
      $intron->set_transcript($tx);
      push @introns, $intron;
    }
  }
  if($#$cds >= 0){
    $clen += $$cds[$#$cds]->stop - $$cds[$#$cds]->start +1;
  }
  $tx->{Coding_Length} = $clen;
  $tx->{Introns} = \@introns;
  #get start/coding start
  my $start = -1;
  if(($#{$cds} >= 0) &&
     (($start == -1) || ($start > $$cds[0]->start))){
    $start = $$cds[0]->start;
  }
  $tx->{Coding_Start} = $start;
  if(($#{$starts} >= 0) &&
     (($start == -1) || ($start > $$starts[0]->start))){
    $start = $$starts[0]->start;
  }
  if(($#{$stops} >= 0) &&
     (($start == -1) || ($start > $$stops[0]->start))){
    $start = $$stops[0]->start;
  }
  if(($#{$exons} >= 0) &&
     (($start == -1) || ($start > $$exons[0]->start))){
    $start = $$exons[0]->start;
  }
  if(($#{$utr5} >= 0) &&
     (($start == -1) || ($start > $$utr5[0]->start))){
    $start = $$utr5[0]->start;
  }
  if(($#{$utr3} >= 0) &&
     (($start == -1) || ($start > $$utr3[0]->start))){
    $start = $$utr3[0]->start;
  }
  $tx->{Start} = $start;
  #get stop/coding stop
  my $stop = -1;
  if(($#{$cds} >= 0) &&
     (($stop == -1) || ($stop < $$cds[$#$cds]->stop))){
    $stop = $$cds[$#$cds]->stop;
  }
  $tx->{Coding_Stop} = $stop;
  if(($#{$starts} >= 0) && 
     (($stop == -1) || ($stop < $$starts[$#$starts]->stop))){
    $stop = $$starts[$#$starts]->stop;
  }
  if(($#{$stops} >= 0) &&
     (($stop == -1) || ($stop < $$stops[$#$stops]->stop))){
    $stop = $$stops[$#$stops]->stop;
  }
  if(($#{$exons} >= 0) &&
     (($stop == -1) || ($stop < $$exons[$#$exons]->stop))){
    $stop = $$exons[$#$exons]->stop;
  }
  if(($#{$utr5} >= 0) &&
     (($stop == -1) || ($stop < $$utr5[$#$utr5]->stop))){
    $stop = $$utr5[$#$utr5]->stop;
  }
  if(($#{$utr3} >= 0) &&
     (($stop == -1) || ($stop < $$utr3[$#$utr3]->stop))){
    $stop = $$utr3[$#$utr3]->stop;
  }
  $tx->{Stop} = $stop;
}



###############################################################################
# GTF::Feature
###############################################################################
# The GTF::Feature object stores all data for a single feature of a gtf
# file.  This object basically stores all data for a single non-comment
# line of a gtf file.
package GTF::Feature;
# GTF::Feature::new(feature, start, end, score, frame) 
#    This is the constructor for GTF::Feature objects.  It takes 5 required
#    arguments.  They are:
#      type       - The <feature> field of the line.  Should be either
#                   'cds', 'exon', 'start_codon', or 'stop_codon'.
#      start      - The <start> field of the line.  Should be a number.
#      end        - The <end> field of the line.  Should be a number.
#      score      - The <score> field of the line.  Should be a number or '.'.
#      frame      - The <frame> field of the line.  Should be either
#                   '0', '1', '2', or '.'.
sub new {
  my $feature = bless {};
  ($feature->{Type}, $feature->{Start}, $feature->{End}, $feature->{Score}, 
   $feature->{Frame}) = @_;    
  my $start = $feature->{Start};
  my $stop = $feature->{Stop};
  if(defined($start) && defined($stop) && ($start > $stop)){
    $feature->{Start} = $stop;
    $feature->{Stop} = $start;
  }
  $feature->{Transcript} = 0;
  $feature->{Seq} = 0;
  $feature->{Conseq} = 0;
  $feature->{Subtype} = 0;
  $feature->{ASE} = 0;          # rpz - Alternative Splicing Event
  $feature->{Match} = -1;
  $feature->{Mismatch} = -1;
  $feature->{Unaligned} = -1;
  $feature->{GC} = -1;
  $feature->{A} = -1;
  $feature->{C} = -1;
  $feature->{G} = -1;
  $feature->{T} = -1;
  $feature->{N} = -1;
  $feature->{0} = -1;
  $feature->{1} = -1;
  $feature->{2} = -1;
  return $feature;
}

# GTF::Feature::copy()
#    This function returns a new object which is a copy of this object
sub copy{
  my ($feature) = @_;
  my $copy = GTF::Feature::new($feature->type,$feature->start,$feature->stop,
                               $feature->score,$feature->frame);
  $copy->set_bases($feature->get_a_count,$feature->get_c_count,$feature->get_g_count,
                   $feature->get_t_count,$feature->get_n_count);
  $copy->set_conseq($feature->get_match_count,$feature->get_mismatch_count,
                    $feature->get_unaligned_count);
  return($copy);
}


# GTF::Feature::offset(ammount)
#    This function will offset all locations in this Feature by the given ammount
sub offset{
  my ($feature,$offset) = @_;
  unless(defined($offset)){
    print STDERR "Undefined value passed to GTF::Feature::offset.\n";
    return;
  }
  $feature->{Start} += $offset;
  $feature->{End} += $offset;
}

# GTF::Feature::reverse_complement(seq_length)
#    This function takes the length of the sequence this feature came from and 
#    reverse complements it.  
sub reverse_complement{
  my ($feature, $seq_length) = @_;
  unless(defined($seq_length)){
    print STDERR 
        "Undefined value for seq_length passed to ",
        "GTF::Feature::reverse_complement.\n";
  }
  my $start = $seq_length - $feature->{Start} + 1;
  my $stop = $seq_length - $feature->{End} + 1;
  if($start < $stop){
    $feature->{Start} = $start;
    $feature->{End} = $stop;
  }
  else{
    $feature->{Start} = $stop;
    $feature->{End} = $start;
  }
}   

# GTF::Feature::type()
#    This function returns the feature type, <feature> field of this line as
#    a string.  Should be either 'CDS', 'exon', 'start_codon', or 'stop_codon'.
sub type     {shift->{Type}}

sub set_type{
  my ($feature,$type) = @_;
  $feature->{Type} = $type;
}

# GTF::Feature::type()
#    This function returns the feature sub type, currently inly implemented for 
#    cds, should be 'Initial', 'Internal', 'Terminal', or 'Single'
sub subtype     {shift->{Subtype}}

sub set_subtype{
  my ($feature,$subtype) = @_;
  $feature->{Subtype} = $subtype;
}

# GTF::Feature::ase()
#   This function returns the feature's alternative splicing event, 
#   if there is one.  Currently the only implemented ASE is 'Optional'
# rpz
sub ase         {shift->{ASE}}

sub set_ase {
    my ($feature, $ASE) = @_;
    $feature->{ASE} = $ASE;
}

# GTF::Feature::length()
#    This function returns the length of this feature in base pairs
#    as a number.
sub length{
  my ($feature) = @_;
  return($feature->stop - $feature->start + 1);
}

# GTF::Feature::start()
#    This function returns the start location ,<start> field, of this feature
#    as a number.
sub start       {shift->{Start}}

# GTF::Feature::set_start()
#    This function sets the start location ,<start> field, of this feature
#    as a number.
sub set_start{
  my ($feature,$start) = @_;
  $feature->{Start} = $start;
}

# GTF::Feature::end()
#    This function returns the end location, <end> field, of this feature
#    as a number.
sub end         {shift->{End}}

# GTF::Feature::stop()
#    This function returns the end location, <end> field, of this feature
#    as a number.
#    Same as GTF::Feature::end();
sub stop         {shift->{End}}

# GTF::Feature::set_stop()
#    This function sets the stop location ,<stop> field, of this feature
#    as a number.
sub set_stop{
  my ($feature,$stop) = @_;
  $feature->{End} = $stop;
}

# GTF::Feature::score()
#    This function returns the score, <score> field, of this line.  This will 
#    be either a number or '.'.
sub score       {shift->{Score}}


# GTF::Feature::frame()
#    This function returns the frame, <frame> field, of this line.  Should be
#    '0', '1', '2', or '.'.
sub frame       {shift->{Frame}}

# GTF::Feature::set_frame()
#    This function sets the frame, <frame> field, of this line.  Should be
#    '0', '1', '2', or '.'.
sub set_frame{
  my ($feature,$frame) = @_;
  if(($frame eq '0') || ($frame eq '1') ||
     ($frame eq '2') || ($frame eq '.')){
    $feature->{Frame} = $frame;
  }
  else{
    print STDERR "Bad Frame Value \'$frame\' should be \'0\', \'1\',".
        "\'2\', or \'.\'.\n";
  }
}

# GTF::Feature::transcript_id()
#    This function returns the transcript_id of this exon's transcript
sub transcript_id{
  my ($feature) = @_;
  my $tx = $feature->{Transcript};
  return $tx->id;
}

# GTF::Feature::transcript()
#    This function returns the transcript of this exon's transcript
sub transcript{
  my ($feature) = @_;
  return $feature->{Transcript};
}

# GTF::Feature::transcript()
#    This function returns the transcript of this exon's transcript
sub set_transcript{
  my ($feature,$tx) = @_;
  $feature->{Transcript} = $tx;
}

# GTF::Feature::transcript_id()
#    This function returns the transcript_id of this exon's transcript
sub gene_id{
  my ($feature) = @_;
  my $tx = $feature->transcript;
  my $gene = $tx->gene;
  return $gene->id;
}

# GTF::Feature::transcript()
#    This function returns the transcript of this exon's transcript
sub gene{
  my ($feature) = @_;
  my $tx = $feature->{Transcript};
  return $tx->gene;
}


# GTF::Feature::seqname()
#    This function returns the sequnce name, <seqname> field, that this 
#    feature came from as a string.
sub seqname{
  my ($feature) = @_;
  my $tx = $feature->{Transcript};
  return $tx->seqname;
}    

# GTF::Feature::source()
#    This function returns the source, <source> field, of this line as a 
#    string.
sub source{
  my ($feature) = @_;
  my $tx = $feature->{Transcript};
  return $tx->source;
}

# GTF::Feature::strand()
#    This function returns the strand, <strand> field, of this line as a 
#   string.  Should be '+', '-', or '.'.
sub strand{
  my ($feature) = @_;
  my $tx = $feature->{Transcript};
  return $tx->strand;
}

sub equals{
  my ($feature,$compare) = @_;
  unless($feature && $compare){
    return 0;
  }
  if(($feature->start == $compare->start) &&
     ($feature->stop == $compare->stop) &&
     ($feature->strand eq $compare->strand) &&
     ($feature->type eq $compare->type)){
    return 1;
  }
  return 0;
}

# GTF::Feature::output_gtf([file_handle])
#    This function outputs the information for this feature in standard 
#    gtf2 format.  If takes and optional file_handle as the only argument 
#    to which it outputs the data.  If no argument is given it outputs 
#    the data to stdout.
sub output_gtf{
  my ($feature,$out_handle) = @_;
  unless(defined $out_handle){
    $out_handle = \*STDOUT;
  }
  print $out_handle
      "".$feature->seqname."\t".$feature->source."\t$feature->{Type}\t",
      "$feature->{Start}\t$feature->{End}\t$feature->{Score}\t".$feature->strand."\t",
      "$feature->{Frame}\tgene_id \"".$feature->gene_id."\"; transcript_id \"",
      $feature->transcript_id."\";\n";
}   

# GTF::Feature::output_gff([file_handle])
#    This function outputs the information for this feature in standard 
#    GFF format.  If takes and optional file_handle as the only argument 
#    to which it outputs the data.  If no argument is given it outputs 
#    the data to stdout.
sub output_gff{
  my ($feature,$out_handle) = @_;
  unless(defined $out_handle){
    $out_handle = \*STDOUT;
  }
  
  print $out_handle
      "".$feature->seqname."\t".$feature->source."\t$feature->{Type}\t",
      "$feature->{Start}\t$feature->{End}\t$feature->{Score}\t".$feature->strand."\t",
      "$feature->{Frame}\t".$feature->transcript_id."\n";
}   

# GTF::Feature::set_bases(A,C,G,T,N)
#    This funtion allows you to set the number of A,C,G,T nucleotides 
#    in this feature
sub set_bases{    
  my ($feature,$a,$c,$g,$t,$n) = @_;
  $feature->{A} = $a;
  $feature->{C} = $c;
  $feature->{G} = $g;
  $feature->{T} = $t;    
  $feature->{N} = $n;    
  $feature->{Seq} = 1;
}

# GTF::Feature::set_conseq(n0,n1,n2)
#    This funtion allows you to set the number of 0,1,2 in the  
#    conservation sequnce for this feature
sub set_conseq{
  my ($feature,$n0,$n1,$n2) = @_;
  $feature->{0} = $n0;
  $feature->{1} = $n1;
  $feature->{2} = $n2;
  $feature->{Conseq} = 1;
}

# GTF::Feature::gc_percentage()
#    This funciton returns the gc_percentage for this feature if the 
#    sequence was loaded otherwise it returns -1.
sub gc_percentage{
  my ($feature) = @_;
  if($feature->{GC} != -1){
    return $feature->{GC};
  }
  unless($feature->{Seq}){
    return -1;
  }
  my $a = $feature->{A};
  my $c = $feature->{C};
  my $g = $feature->{G};
  my $t = $feature->{T};
  my $n = $feature->{N};
  my $total = $a + $c + $g + $t + $n;
  my $gc = $g + $c;
  my $percent = $gc/$total;
  $feature->{GC} = $percent;
  return $percent;
}

# GTF::Feature::match_percentage()
#    This funciton returns the percentage of matches in the conservation sequence
#    for this feature if it was loaded otherwise it returns -1.
sub match_percentage{
  my ($feature) = @_;
  if($feature->{Match} != -1){
    return $feature->{Match};
  }
  unless($feature->{Conseq}){
    return -1;
  }
  my $match = $feature->{1};
  my $mismatch = $feature->{0};
  my $unaligned = $feature->{2};
  my $total = $match + $mismatch + $unaligned;
  my $percent = $match/$total;
  $feature->{Match} = $percent;
  return $percent;
}

# GTF::Feature::mismatch_percentage()
#    This funciton returns the percentage of mismatches in the conservation sequence
#    for this feature if it was loaded otherwise it returns -1.
sub mismatch_percentage{
  my ($feature) = @_;
  if($feature->{Mismatch} != -1){
    return $feature->{Mismatch};
  }
  unless($feature->{Conseq}){
    return -1;
  }
  my $match = $feature->{1};
  my $mismatch = $feature->{0};
  my $unaligned = $feature->{2};
  my $total = $match + $mismatch + $unaligned;
  my $percent = $mismatch/$total;
  $feature->{Mismatch} = $percent;
  return $percent;
}

# GTF::Feature::unaligned_percentage()
#    This funciton returns the percentage of unaligned nucs in the conservation sequence
#    for this feature if it was loaded otherwise it returns -1.
sub unaligned_percentage{
  my ($feature) = @_;
  if($feature->{Unaligned} != -1){
    return $feature->{Unaligned};
  }
  unless($feature->{Conseq}){
    return -1;
  }
  my $match = $feature->{1};
  my $mismatch = $feature->{0};
  my $unaligned = $feature->{2};
  my $total = $match + $mismatch + $unaligned;
  my $percent = $unaligned/$total;
  $feature->{Unaligned} = $percent;
  return $percent;
}

# GTF::Feature::get_a_count()
#    This funciton returns the number of A nucleotides in this features sequence
#    if the sequence was loaded, otherwise it returns -1.
sub get_a_count{   
  my ($feature) = @_;
  unless($feature->{Seq}){
    return -1;
  }
  return $feature->{A};
}

# GTF::Feature::get_c_count()
#    This funciton returns the number of C nucleotides in this features sequence
#    if the sequence was loaded, otherwise it returns -1.
sub get_c_count{
  my ($feature) = @_;
  unless($feature->{Seq}){
    return -1;
  }
  return $feature->{C};
}

# GTF::Feature::get_g_count()
#    This funciton returns the number of G nucleotides in this features sequence
#    if the sequence was loaded, otherwise it returns -1.
sub get_g_count{
  my ($feature) = @_;
  unless($feature->{Seq}){
    return -1;
  }
  return $feature->{G};
}

# GTF::Feature::get_t_count()
#    This funciton returns the number of T nucleotides in this features sequence
#    if the sequence was loaded, otherwise it returns -1.
sub get_t_count{
  my ($feature) = @_;
  unless($feature->{Seq}){
    return -1;
  }
  return $feature->{T};
}

# GTF::Feature::get_n_count()
#    This funciton returns the number of N nucleotides in this features sequence
#    if the sequence was loaded, otherwise it returns -1.
sub get_n_count{
  my ($feature) = @_;
  unless($feature->{Seq}){
    return -1;
  }
  return $feature->{N};
}

# GTF::Feature::get_match_count()
#    This funciton returns the number of matched nucleotides in this features 
#    conservation sequence if it was loaded, otherwise it returns -1.
sub get_match_count{
  my ($feature) = @_;
  unless($feature->{Conseq}){
    return -1;
  }
  return $feature->{1};
}

# GTF::Feature::get_mismatch_count()
#    This funciton returns the number of mismatched nucleotides in this features 
#    conservation sequence if it was loaded, otherwise it returns -1.
sub get_mismatch_count{
  my ($feature) = @_;
  unless($feature->{Conseq}){
    return -1;
  }
  return $feature->{0};
}

# GTF::Feature::get_unaligned_count()
#    This funciton returns the number of unaligned nucleotides in this features 
#    conservation sequence if it was loaded, otherwise it returns -1.
sub get_unaligned_count{
  my ($feature) = @_;
  unless($feature->{Conseq}){
    return -1;
  }
  return $feature->{2};
}

sub tag{
  my ($feature) = @_;
  return $feature->{Tag};
}

sub set_tag{
  my ($feature,$tag) = @_;
  $feature->{Tag} = $tag;
}

1;
__END__
