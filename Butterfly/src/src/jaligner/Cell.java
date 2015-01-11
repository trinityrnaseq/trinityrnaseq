/*
 * $Id: Cell.java,v 1.1 2005/05/25 19:56:30 ahmed Exp $
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

package jaligner;

/**
 * A cell in a similarity matrix, to hold row, column and score.
 * 
 * @author Ahmed Moustafa (ahmed@users.sf.net)
 */

public class Cell {
	/**
	 * Row of the cell
	 */
	private int row;
	/**
	 * Column of the cell
	 */
	private int col;
	/**
	 * Alignment score at this cell
	 */
	private float score;
	
	/**
	 * Constructor
	 */
	public Cell() {
		super();
		this.row = 0;
		this.col = 0;
		this.score = Float.NEGATIVE_INFINITY;
	}
	/**
	 * @return Returns the col.
	 */
	public int getCol() {
		return this.col;
	}
	/**
	 * @param col The col to set.
	 */
	public void setCol(int col) {
		this.col = col;
	}
	/**
	 * @return Returns the row.
	 */
	public int getRow() {
		return this.row;
	}
	/**
	 * @param row The row to set.
	 */
	public void setRow(int row) {
		this.row = row;
	}
	/**
	 * @return Returns the score.
	 */
	public float getScore() {
		return this.score;
	}
	/**
	 * @param score The score to set.
	 */
	public void setScore(float score) {
		this.score = score;
	}
	
	/**
	 * Sets the row, column and score of the cell.
	 * @param row The row to set.
	 * @param col The col to set.
	 * @param score The score to set.
	 */
	public void set(int row, int col, float score) {
		this.row = row;
		this.col = col;
		this.score = score;
	}
}