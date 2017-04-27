
test:
# --min_expr_any
	../../util/filter_low_expr_transcripts.pl -m kallisto-trans.TMM.EXPR.matrix --transcripts ../test_DATA/Trinity.fasta --min_expr_any 1.0 > filt.min_expr_any.1.0.out 2>filt.min_expr_any.1.0.err
# --min_pct_dom_iso
	../../util/filter_low_expr_transcripts.pl -t ../test_DATA/Trinity.fasta -m kallisto-trans.TMM.EXPR.matrix --trinity_mode --min_pct_dom_iso 70 > filt.min_pct_dom_iso.70.out 2>filt.min_pct_dom_iso.70.err
# --highest_iso_only
	../../util/filter_low_expr_transcripts.pl -t ../test_DATA/Trinity.fasta -m kallisto-trans.TMM.EXPR.matrix --trinity_mode --highest_iso_only > filt.highest_iso_only.out 2>filt.highest_iso_only.err
	touch test


clean:
	rm -f ./filt.*
	rm -f ./test
