/*
 * $Id: AlignCommandLine.java,v 1.7 2005/04/18 05:37:52 ahmed Exp $
 * 
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 * 
 * (bhaas: modified to provide a simple interface to the Jaligner library) 04-02-2011
 * 
 */

import jaligner.Alignment;
import jaligner.NeedlemanWunschGotohBanded;
import jaligner.Sequence;
import jaligner.SmithWatermanGotoh;
import jaligner.NeedlemanWunschGotoh;
import jaligner.formats.Pair;
import jaligner.matrix.Matrix;
import jaligner.matrix.MatrixLoader;
import jaligner.matrix.MatrixGenerator;
import jaligner.util.Commons;
import jaligner.util.SequenceParser;

import java.io.File;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Command line interface for JAligner.
 * 
 * @author Ahmed Moustafa (ahmed@users.sf.net)
 */

public class NWalign {
	
	private static final Logger logger = Logger.getLogger(NWalign.class.getName());
	
	/**
	 * @param args The command line arguments 
	 */
	public static void main(String[] args) {
	    logger.info( Commons.getJAlignerInfo() ); 
				
        if (args.length < 3) {
			logger.severe( "Invalid number of arguments: " + args.length );
			printUsage();
			System.exit(1);
        
        } else {
			
			try {
				
				String f1 = args[0];					// file name of sequence #1
				String f2 = args[1];					// file name of sequence #1
				String mode = args[2]; // S or N  (smith-waterman or needleman-wunsch)
				
				
				
				//String m = args[2];						// scoring matrix id or file name user-defined scoring matrix
				//float o = Float.parseFloat(args[3]);	// open gap penalty
				//float e = Float.parseFloat(args[4]);	// extend gap penalty
				
				Sequence s1 = SequenceParser.parse(new File(f1));
				Sequence s2 = SequenceParser.parse(new File(f2));
				
				int match_score = 4;
				int mismatch_penalty = -5;
				int open_penalty = 10;
				int extend_penalty = 1;
				
				Matrix matrix = MatrixGenerator.generate(match_score, mismatch_penalty);
				
				//Matrix matrix = MatrixLoader.load(m);
				Alignment alignment;
				if (mode.equals("S")) {
					alignment = SmithWatermanGotoh.align (s1, s2, matrix, open_penalty, extend_penalty);
				}
				else if (mode.equals("N")){
					alignment = NeedlemanWunschGotoh.align(s1, s2, matrix, open_penalty, extend_penalty);
				}
				else if (mode.equals("NB")) {
					int bandwidth = Integer.parseInt(args[3]);
					alignment = NeedlemanWunschGotohBanded.align(s1, s2, matrix, open_penalty, extend_penalty, bandwidth);
				}
				else {
					printUsage();
					throw new RuntimeException("error");
				}
				
				System.out.println (alignment.getSummary());
				System.out.println (new Pair().format(alignment));
				
				AlignmentStats stats = get_alignment_stats(alignment);
				
				int alignment_length = stats.alignment_length;
				int matches = stats.matches;
				int mismatches = stats.mismatches;
				int gaps = stats.gaps;
				
				float percent_A_in_alignment = (float) stats.get_count_of_bases_in_aligned_region(s1.getId()) / (s1.length()) * 100;
				float percent_B_in_alignment = (float) stats.get_count_of_bases_in_aligned_region(s2.getId()) / (s2.length()) * 100;
					
				float max_percent_aligned = Math.max(percent_A_in_alignment, percent_B_in_alignment);
					
				float percent_identity = (float)matches/(matches+mismatches) * 100;
				float percent_gapped = (float)gaps/alignment_length * 100;
				
				System.out.println("Matches: " + matches + ", Mismatches: " + mismatches + ", gaps: " + gaps + ", align_len: " + alignment_length);
				System.out.println("percent_identity: " + percent_identity + ", percent_gapped: " + percent_gapped);
				System.out.println("max_percent_aligned: " + max_percent_aligned + "\n");
				
			} catch (Exception e) {
				logger.log(Level.SEVERE, "Failed processing the command line: " + e.getMessage(), e); 
				System.exit(1);
			}
		} 
		
	}
	
	/**
	 * Prints the syntax for using JAligner 
	 */
	private static void printUsage( ) {
		StringBuffer buffer = new StringBuffer();
		buffer.append ( "\n" );
		buffer.append ( "Usage:\n" );
		buffer.append ( "------\n" );
		buffer.append ( "java -jar jaligner.jar <s1> <s2> (S, N, or NB) [bandwidth=10]\n" );
		buffer.append("\tS: Smith-Waterman\n");
		buffer.append("\tN: Needleman-Wunsch\n");
		buffer.append ( "\n" ) ;
		logger.info(buffer.toString());
	}
	
	// utility methods:
	
	public static Alignment run_SW_alignment(String name1, String s1, String name2, String s2, int match, int mismatch, int open, int extend) {
		
		Matrix matrix = MatrixGenerator.generate(match, mismatch);
		
		Sequence seq1 = new Sequence(name1, s1);
		Sequence seq2 = new Sequence(name2, s2);
		Alignment alignment = SmithWatermanGotoh.align(seq1, seq2, matrix, open, extend);
		
		return(alignment);
	}
	
	public static Alignment run_NW_alignment(String name1, String s1, String name2, String s2, 
				int match, int mismatch, int open, int extend) {
		
		Matrix matrix = MatrixGenerator.generate(match, mismatch);
		
		Sequence seq1 = new Sequence(name1, s1);
		Sequence seq2 = new Sequence(name2, s2);
		Alignment alignment = NeedlemanWunschGotoh.align(seq1, seq2, matrix, open, extend);
		
		return(alignment);
	}
	

	public static Alignment run_NW_banded_alignment(String name1, String s1, String name2, String s2, 
			int match, int mismatch, int open, int extend, int bandwidth) {

		Matrix matrix = MatrixGenerator.generate(match, mismatch);

		Sequence seq1 = new Sequence(name1, s1);
		Sequence seq2 = new Sequence(name2, s2);
		Alignment alignment = NeedlemanWunschGotohBanded.align(seq1, seq2, matrix, open, extend, bandwidth);

		return(alignment);
	}

	
	
		
	public static AlignmentStats get_alignment_stats(Alignment a) {
		AlignmentStats stats = new AlignmentStats(a);
		return(stats);
	}
	
	
}
	
	
	

