#!/usr/bin/env perl

use strict;
use warnings;

use CGI;
use CGI::Carp qw(fatalsToBrowser);
use FindBin;
use File::Basename;
use POSIX;


$|++;


#my $BASE_URL = "ftp://ftp.broad.mit.edu/pub/users/bhaas/trinity_sample";
my $BASE_URL = "http://www.broadinstitute.org/~bhaas";

main: {
        
    my $cgi = new CGI();
    print $cgi->header();

    print $cgi->start_html(-title => "RNA-Seq Sample Report",
                           -style => { -src => ['CSS/common.css',
                                                ]
                                                },
                           -script => [ 
                                        { -type => 'text/javascript',
                                          -src => 'http://www.java.com/js/deployJava.js',
                                      },
                                        { -type => 'text/javascript',
                                          -src => 'http://genomeview.org/start/genomeview.js',
                                      },
                                        { -type => 'text/javascript',
                                          -src => 'http://ajax.googleapis.com/ajax/libs/jquery/1.7.1/jquery.min.js',
                                      }
                                        ]
                           
                           );

    print "<script>  <!-- needed for genomeview communication -->\n"
        . " instanceID = 'ALL'\n"
        . "</script>\n";
    
    

    my $genomeview_js = <<__EOJS__;


       <script type="text/javascript"> 
           var instanceID="ALL";
    var gv_url= "$BASE_URL/Trinity.fasta";
    var gv_config= "$BASE_URL/gv_config.txt";
    var gv_location= null;
    var gv_extra    = __ANNOT_TRACKS__;
              
    var wsUrl = "http://genomeview.org/start/launch.jnlp? --config " + gv_config + " --url " + gv_url + " " + gv_extra;
    

    function launchGV () {
                   
                  // launch as applet, but tied to existing page
                      //startGV(gv_url,gv_location,gv_config,gv_extra,500,500);
                   

                  // launch using webstart, remains active across multiple pages.
                      jQuery('#frame1').attr("src", wsUrl);

              }


         </script>

__EOJS__

;

    my @rest_tracks = ("$BASE_URL/best_candidates.gff3", "$BASE_URL/bowtie_out.coordSorted.sam.plus.bam");
    
    my $annot_tracks_text = "\"" . join("\"\n +\" ", @rest_tracks) . "\"";
    
    $genomeview_js =~ s/__ANNOT_TRACKS__/$annot_tracks_text/;
    
    print $genomeview_js;
    print "<iframe id=\"frame1\" style=\"display:none\"></iframe>\n";
    

    ## toolbar buttons
    print "<div id='toolbar' style=\"border:1px solid black; margin:auto; width:60%; background-color:#004242\">\n"
        . "    <button id='GenomeViewBttn'>GenomeView</button>\n"
        . "</div>\n\n";
        
    ## add action initialization code

    print "<script>\n"
        . "jQuery(function() { \n"
        . "    jQuery('#GenomeViewBttn').click( function() { launchGV(); });\n"
        . "});\n"
        . "</script>\n";




    ## add some example links.
    my @accs = qw( comp65_c1_seq6
                   comp28_c0_seq1
                   comp33_c0_seq1 );


    print "<ul>Example links:\n";
    foreach my $acc (@accs){
        print "<li><a href='javascript:void(0);' onClick=\"javascript:positionGV(\'$acc:1:100\');\" >$acc</a>\n";
    }
    print "</ul>\n";
    

    
    print $cgi->end_html();


    exit(0);

}
