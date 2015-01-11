/*
 * $Id: CLUSTAL.java,v 1.1 2005/05/25 19:56:30 ahmed Exp $
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

/**
 * CLUSTAL format.
 * Example:
 * <pre>
 * <small>
 * CLUSTAL_FORMAT W(1.60) multiple sequence alignment
 *
 *
 * JC2395          NVSDVNLNK---YIWRTAEKMK---ICDAKKFARQHKIPESKIDEIEHNSPQDAAE----
 * KPEL_DROME      MAIRLLPLPVRAQLCAHLDAL-----DVWQQLATAVKLYPDQVEQISSQKQRGRS-----
 * FASA_MOUSE      NASNLSLSK---YIPRIAEDMT---IQEAKKFARENNIKEGKIDEIMHDSIQDTAE----
 *
 *
 * JC2395          -------------------------QKIQLLQCWYQSHGKT--GACQALIQGLRKANRCD
 * KPEL_DROME      -------------------------ASNEFLNIWGGQYN----HTVQTLFALFKKLKLHN
 * FASA_MOUSE      -------------------------QKVQLLLCWYQSHGKS--DAYQDLIKGLKKAECRR
 *
 *
 * JC2395          IAEEIQAM
 * KPEL_DROME      AMRLIKDY
 * FASA_MOUSE      TLDKFQDM
 * </small>
 * </pre>
 *
 * @author Ahmed Moustafa (ahmed@users.sf.net)
 */

public class CLUSTAL extends Format {
	
	/**
	 * Name width 
	 */
    private static final int NAME_WIDTH = 36;
	
    /**
	 * Sequence width
	 */
    private static final int SEQUENCE_WIDTH = 50;
	
    /**
	 * CLUSTAL header
	 */
    private static final String HEADER = "CLUSTAL_FORMAT W(1.60) multiple sequence alignment\n\n";
	
	/**
	 * Constructor
	 */
	public CLUSTAL() {
		super( );
		setId("CLUSTAL");
	}
	
	/**
	 * Returns CLUSTAL format
	 * @param names array of the names of the sequences.
	 * @param sequences	array of the sequences 
	 */
	public String format(String[] names, String[] sequences) {
		StringBuffer buffer = new StringBuffer (HEADER);
		int maxSequenceLength = 0;
		for (int i = 0; i < sequences.length; i++) {
			if (sequences[i].length() > maxSequenceLength) {
				maxSequenceLength = sequences[i].length(); 
			}
		}
		
		for (int i = 0; i * SEQUENCE_WIDTH < maxSequenceLength; i++) {
			
			for (int j = 0; j < sequences.length; j++) {
				buffer.append(  NAME_WIDTH <= names[j].length() ? names[j].substring(0, NAME_WIDTH - 1) : names[j] ); 
				for (int k = names[j].length(); k < NAME_WIDTH; k++) {
					buffer.append(" ");
				}
				if (names[j].length() >= NAME_WIDTH) {
					buffer.append(" ");
				}
				buffer.append(sequences[j].substring(i * SEQUENCE_WIDTH, ((i + 1) * SEQUENCE_WIDTH) < sequences[j].length() ? (i + 1) * SEQUENCE_WIDTH: sequences[j].length()));
				if (j < sequences.length) {
					buffer.append("\n");
				}
			}
			if ((i + 1) * SEQUENCE_WIDTH < maxSequenceLength) {
				buffer.append("\n\n");
			}
		}
		return buffer.toString();
	}

	/**
	 * Returns CLUSTAL format of the alignment
	 * @param alignment ({@link Alignment})
	 * @return CLUSTAL format of the alignment 
	 */
	public String format (Alignment alignment) {
		String[] sequences = {new String(alignment.getSequence1()), new String(alignment.getSequence2())};
		String[] names = {alignment.getName1(), alignment.getName2()};
		return format (names, sequences);
	}
}