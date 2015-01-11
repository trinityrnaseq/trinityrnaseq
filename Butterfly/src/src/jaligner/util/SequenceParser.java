/*
 * @author Ahmed Moustafa (ahmed at users.sourceforge.net)
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

package jaligner.util;


import jaligner.Sequence;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.InputStreamReader;
import java.util.StringTokenizer;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * SequenceParser to sequences from different formats.
 * <br>
 * Currently the supported formats are:
 * <ul>
 * <li>Plain sequence</li>
 * <li><a href="http://www.ncbi.nlm.nih.gov/BLAST/fasta.html">Sequence</a></li>
 * </ul>
 *
 * @author Ahmed Moustafa (ahmed@users.sf.net)
 */

public class SequenceParser {
    
    /**
     * Logger
     */
    private static final Logger logger = Logger.getLogger(SequenceParser.class.getName());
	
	/**
	 * Returns a parsed Sequence from a string.
	 * @param sequence string to parse
	 * @return parsed sequence
	 * @throws SequenceParserException
	 * @see Sequence
	 */
	public static Sequence parse (String sequence) throws SequenceParserException {
		if (sequence == null) {
			throw new SequenceParserException ( "Null sequence" );
		}

		if (sequence.trim().length() == 0) {
			throw new SequenceParserException ( "Empty sequence" );
		}
		
		sequence = sequence.replaceAll("\r\n", "\n");
		
		String sequenceName = null;
		String sequenceDescription = null;
		
		if (sequence.startsWith(">")) {
			// FASTA format
			int index = sequence.indexOf("\n");
				
			if (index == -1) {
				throw new SequenceParserException ( "Invalid sequence format" );
			}
				
			String first = sequence.substring(1, index);
			sequence = sequence.substring(index);
				
			index = 0;
			for (int i = 0; i < first.length() && first.charAt(i) != ' ' && first.charAt(i) != '\t'; i++, index++) {
			    // Skip white spaces
			}
			sequenceName = first.substring(0, index);
			StringTokenizer stringTokenizer = new StringTokenizer(sequenceName, "|");
			while (stringTokenizer.hasMoreTokens()) {
				sequenceName = stringTokenizer.nextToken();
			}
			sequenceDescription = index + 1 > first.length() ? "" : first.substring(index + 1);
		} else {
			// Plain format ... nothing to do here
		}
		
    	Sequence s = new Sequence(prepare(sequence), sequenceName, sequenceDescription, Sequence.PROTEIN);
    	
    	return s;
	}

	/**
	 * Returns a Sequence parsed and loaded from a file
	 * @param file to parse
	 * @return parsed sequence 
	 * @throws SequenceParserException
	 * @see Sequence
	 */
	public static Sequence parse (File file) throws SequenceParserException {
	    String sequenceName = null;
		String sequenceDescription = null;
		BufferedReader reader = null;
		try {
			reader = new BufferedReader(new InputStreamReader(new FileInputStream(file)));
			StringBuffer buffer = new StringBuffer();
			
			// Read & parse the first line
			String line = reader.readLine();
				
			if (line.startsWith(">")) {
				// FASTA sequence
				
				line = line.substring(1).trim();
				int index = 0;
				for (int i = 0; i < line.length() && line.charAt(i) != ' ' && line.charAt(i) != '\t'; i++, index++) {
				    // Skip white spaces
				}
				
				sequenceName = line.substring(0, index);
				StringTokenizer stringTokenizer = new StringTokenizer(sequenceName, "|");
				while (stringTokenizer.hasMoreTokens()) {
					sequenceName = stringTokenizer.nextToken();
				}
				sequenceDescription = index + 1 > line.length() ? "" : line.substring(index + 1);
			} else {
				// Plain sequence
				buffer.append(prepare(line));
			}
			
			// Read the remaining the file (the actual sequence)
			while ((line = reader.readLine()) != null) {
				buffer.append(prepare(line));
			}
			reader.close();
	
	    	Sequence s = new Sequence(buffer.toString(), sequenceName, sequenceDescription, Sequence.PROTEIN);
	    	return s;
		} catch (Exception e) {
		    throw new SequenceParserException(e.getMessage());
		} finally {
		    if (reader != null) {
		        try {
		            reader.close();
		        } catch (Exception silent) {
		            logger.log(Level.WARNING, "Failed closing reader: " + silent.getMessage(), silent);
		        }
		    }
		}
    	
	}

	/**
	 * Removes whitespaces from a sequence and validates other characters.
	 * @param sequence sequence to be prepared
	 * @return prepared array of characters
	 * @throws SequenceParserException
	 */
	private static String prepare (String sequence) throws SequenceParserException {
		StringBuffer buffer = new StringBuffer();
		String copy = sequence.trim().toUpperCase();
		
		for (int i = 0, n = copy.length(); i < n; i++) {
			switch ( copy.charAt(i) ) {
				// skip whitespaces
				case 9:
				case 10:
				case 13:
				case 32: break;
				
				// add a valid character
				case 'A':
				case 'B':
				case 'C':
				case 'D':
				case 'E':
				case 'F':
				case 'G':
				case 'H':
				case 'I':
				case 'K':
				case 'L':
				case 'M':
				case 'N':
				case 'P':
				case 'Q':
				case 'R':
				case 'S':
				case 'T':
				case 'U':
				case 'V':
				case 'W':
				case 'Y':
				case 'Z':
				case 'X':
				
				case '-':
				case '*': buffer.append(copy.charAt(i)); break;
							
				// throw an exception for anything else
				default: throw new SequenceParserException( "Invalid sequence character: '" + copy.charAt(i) + "'"); 
			}
		}
		return buffer.toString();
	}
}