/*
 * $Id: SmithWatermanGotohExample.java,v 1.1 2006/02/03 03:46:18 ahmed Exp $
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
 */

package jaligner.example;

import jaligner.Alignment;
import jaligner.Sequence;
import jaligner.SmithWatermanGotoh;
import jaligner.formats.Pair;
import jaligner.matrix.MatrixLoader;
import jaligner.util.SequenceParser;

import java.io.IOException;
import java.io.InputStream;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Example of using JAligner API to align P53 human aganist
 * P53 mouse using Smith-Waterman-Gotoh algorithm.
 *
 * @author Ahmed Moustafa (ahmed@users.sf.net)
 */

public class SmithWatermanGotohExample {
	
	/**
	 * 
	 */
	private static final String SAMPLE_SEQUENCE_P35_HUMAN = "jaligner/example/sequences/p53_human.fasta";
	
	/**
	 * 
	 */
	private static final String SAMPLE_SEQUENCE_P35_MOUSE = "jaligner/example/sequences/p53_mouse.fasta";
	
	/**
	 * Logger
	 */
	private static final Logger logger = Logger.getLogger(SmithWatermanGotohExample.class.getName());
	
	/**
	 * 
	 * @param args
	 */
	public static void main(String[] args) {
        try {
        	logger.info("Running example...");
        	
			Sequence s1 = SequenceParser.parse(loadP53Human());  
			Sequence s2 = SequenceParser.parse(loadP53Mouse());
	        
	        Alignment alignment = SmithWatermanGotoh.align(s1, s2, MatrixLoader.load("BLOSUM62"), 10f, 0.5f);
	        
	        System.out.println ( alignment.getSummary() );
	        System.out.println ( new Pair().format(alignment) );
	        
	        logger.info("Finished running example");
        } catch (Exception e) {
        	logger.log(Level.SEVERE, "Failed running example: " + e.getMessage(), e);
        }
    }
	
	/**
	 * 
	 * @param path location of the sequence
	 * @return sequence string
	 * @throws IOException
	 */
	private static String loadSampleSequence(String path) throws IOException {
		InputStream inputStream = NeedlemanWunschExample.class.getClassLoader().getResourceAsStream(path);
        StringBuffer buffer = new StringBuffer();
        int ch;
        while ((ch = inputStream.read()) != -1) {
            buffer.append((char)ch);
        }
        return buffer.toString();
	}
	
	/**
	 * 
	 * @return sequence string
	 * @throws IOException
	 */
	public static String loadP53Human( ) throws IOException {
		return loadSampleSequence(SAMPLE_SEQUENCE_P35_HUMAN);
	}

	/**
	 * 
	 * @return sequence string
	 * @throws IOException
	 */
	public static String loadP53Mouse( ) throws IOException {
		return loadSampleSequence(SAMPLE_SEQUENCE_P35_MOUSE);
	}
}