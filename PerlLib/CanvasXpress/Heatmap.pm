package CanvasXpress::Heatmap;

use strict;
use warnings;


sub draw {
    my %inputs = @_;

    # structure of input hash:
    #
    #   %inputs = ( samples => [ 'sampleA', 'sampleB', 'sampleC', ...],
    #               value_matrix => [ ['featureA', valA, valB, valC, ...],
    #                                 ['featureB', valA, valB, valC, ...],  ... ,                                 
    #                                ],
    #                    ## and optionally:
    #               feature_tree => "", # string containing the newick formatted tree
    #               sample_tree => "",  # ditto above
    #               feature_descriptions => ['name of feature A', 'name of feature B', ...],
    #               cluster_features => 0|1,   ## clustering done in browser, exclusive with feature_tree option
    #               cluster_samples => 0|1,
    #             );
    #
    
    
    


    my $html = "<!--[if IE]><script type=\"text/javascript\" src=\"canvasXpress-SF/js/excanvas.js\"></script><![endif]-->\n";
    $html .= "<script type=\"text/javascript\" src=\"http://canvasxpress.org/js/canvasXpress.min.js\"></script>\n";
    
    
    $html .= "<script>\n";
    
    $html .= "    var make_heatmap = function() {\n";
    $html .= "    var cx = new CanvasXpress(\"canvas\", {\n";
    

#$html .= "          \"x\": {\n";
#$html .= "          \"Desc\": [\n";
#$html .= "               \"Sample-1\",\n";
#$html .= "               \"Sample-2\",\n";
#$html .= "               \"Sample-3\"\n";
#$html .= "           ] },\n";


    ## reorganize some of the data.
    my @gene_ids;
    my @values;
    my $value_matrix_aref = $inputs{value_matrix};
    foreach my $row (@$value_matrix_aref) {
        my @data = @$row;
        my $gene_id = shift @data;
        push (@gene_ids, $gene_id);
        push (@values, [@data]);
    }
    

    $html .= "            \"y\": {\n";
    $html .= "            \"vars\": [ \n";

    for (my $i = 0; $i <= $#gene_ids; $i++) {
        
        $html .= "             \"$gene_ids[$i]\"";
        if ($i != $#gene_ids) {
            $html .= ",";
        }
        $html .= "\n";
    }
    $html .= "                ],\n";

    my @sample_ids = @{$inputs{samples}};

    $html .= "                \"smps\": [\n";
    
    $html .= "            \"" . join("\",\"", @sample_ids) . "\"\n";
    $html .= "                          ],\n";
    
    
    $html .= "            \"data\": [\n";
    
    for (my $i = 0; $i <= $#values; $i++) {
        $html .= "[ " . join(",", @{$values[$i]}) . "]";
        if ($i != $#values) {
            $html .= ",";
        }
        $html .= "\n";
    }
    $html .= "                   ]\n";
    $html .= "      }\n";
        
    if (my $feature_descriptions_aref = $inputs{feature_descriptions}) {
        $html .= "     ,\n";
        $html .= "   \"z\": { \n";
        $html .= "            \"Desc\": [\n";
#$html .= "                       \"Expression\"\n";
        for (my $i = 0; $i <= $#gene_ids; $i++) {
            
            $html .= "             \"Info on $gene_ids[$i]\"";
            if ($i != $#gene_ids) {
                $html .= ",";
            }
            $html .= "\n";
        }
        $html .= "                   ] },\n";
    };

    if (my $feature_tree = $inputs{feature_tree}) {
    
        $html .= "    ,\n";
        ## gene tree
        $html .= "  \"t\" : {\n";
        $html .= "     \"vars\" : \"$feature_tree\"\n";
        $html .= "    }\n"; 
        
    }
    
    $html .= " },\n";
    
    $html .= " {\n";
    $html .= "     \"graphType\": \"Heatmap\"\n";
    $html .= "     ,\"zoomSamplesDisable\": true\n";
    $html .= "     ,\"smpLabelScaleFontFactor\" : 2\n";

#$html .= "      ,\"smpLabelDescription\": \"Desc\"\n";
    if ($inputs{feature_descriptions}) {
        $html .= "      ,\"varLabelDescription\": \"Desc\"\n";
    }
    if ($inputs{feature_tree}) {
        $html .= "     ,\"showVarDendrogram\": true,\n";
    }
#$html .= "     \"showSmpDendrogram\": true,\n";
    $html .= " });\n";
    
    if ($inputs{cluster_features} ) {
        $html .= "cx.clusterVariables();\n";
    }
    if ($inputs{cluster_samples}) {
        $html .= "cx.clusterSamples();\n";
    }
    
    $html .= " }\n";
    
    $html .= <<__EOJS__;
    
    </script>
        
    <body onload="make_heatmap();">
      <div>
        <canvas id="canvas" width="613" height="500"></canvas> 
      </div>
   </body>
        
__EOJS__

    ;

    return($html);
}
 

1; #EOM

   
