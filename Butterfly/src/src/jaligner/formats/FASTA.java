/*
 * $Id: FASTA.java,v 1.1 2005/05/25 19:56:30 ahmed Exp $
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

package jaligner.formats;

import jaligner.Alignment;
import jaligner.Sequence;

/**
 * <a href="http://www.ncbi.nlm.nih.gov/BLAST/fasta.html">FASTA</a> format.
 * 
 * @author Ahmed Moustafa (ahmed@users.sf.net)
 */

public class FASTA extends Format {
	
	/**
	 * Number of characters per line
	 */
    private static final int LINE_WIDTH = 60;
	
    /**
	 * Constructor for FASTA.
	 */
	public FASTA() {
		super( );
		setId("FASTA");
	}
	
	/**
	 * Returns the name, description and sequence combined in one string.
	 * The length of each line in the sequence is FASTA.LINE_LENGTH
	 * 
	 * @return String
	 */
	public String format (Sequence sequence) {
		StringBuffer buffer = new StringBuffer (">");
		buffer.append(sequence.getId() == null ? "" : sequence.getId());
		buffer.append("\n");
        for (int i = 0, n = sequence.length(); i * LINE_WIDTH < n; i++) {
        	for (int j = i * LINE_WIDTH, m = (i + 1) * LINE_WIDTH < n ? (i + 1) * LINE_WIDTH: n; j < m; j++) {
        		buffer.append(sequence.subsequence(j, 1));
        	}
        	buffer.append("\n");
        }
		return buffer.toString( );
	}

	/**
	 * 
	 * @param alignment
	 * @return FASTA format of the input alignment
	 */
	public String format (Alignment alignment) {
		StringBuffer buffer = new StringBuffer();
		StringBuffer s1 = new StringBuffer();
		StringBuffer s2 = new StringBuffer();
		s1.append(alignment.getSequence1());
		s2.append(alignment.getSequence2());
		buffer.append(format(new Sequence(s1.toString(), alignment.getName1(), "", Sequence.PROTEIN)));
		buffer.append(format(new Sequence(s2.toString(), alignment.getName2(), "", Sequence.PROTEIN)));
		return buffer.toString();
	}
}