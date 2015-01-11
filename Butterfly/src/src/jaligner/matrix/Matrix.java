/*
 * $Id: Matrix.java,v 1.2 2006/01/18 20:16:37 ahmed Exp $
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

package jaligner.matrix;

import java.io.Serializable;

/**
 * Scoring matrix.
 * 
 * @author Ahmed Moustafa (ahmed@users.sf.net)
 */

public class Matrix implements Serializable {
	/**
	 * 
	 */
	private static final long serialVersionUID = 3833742170619524400L;
	
	/**
	 * The size of the scoring matrix. It is the number of the characters in the
	 * ASCII table. It is more than the 20 amino acids just to save the
	 * processing time of the mapping.
	 */
	public static final int SIZE = 127;

	/**
	 * Matrix id (or name)
	 */
	private final String id;

	/**
	 * Scores
	 */
	private final float[][] scores;

	public Matrix(String id, float[][] scores) {
		this.id = id;
		this.scores = scores;
	}

	/**
	 * @return Returns the id.
	 */
	public String getId() {
		return this.id;
	}

	/**
	 * @return Returns the scores.
	 */
	public float[][] getScores() {
		return this.scores;
	}

	/**
	 * 
	 * @param a
	 * @param b
	 * @return score
	 */
	public float getScore(char a, char b) {
		return this.scores[a][b];
	}
}