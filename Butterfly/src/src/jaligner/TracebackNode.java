/*
 * $Id: TracebackNode.java,v 1.1 2006/01/18 20:12:04 ahmed Exp $
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
 * Traceback node
 * 
 * @author Bram Minnaert
 * @author Ahmed Moustafa
 */

public final class TracebackNode {

	/**
	 * Direction
	 */
	private byte direction;

	/**
	 * Length of gap (for up and left directions)
	 */
	private int length;

	/**
	 * Constructor
	 */
	public TracebackNode() {
		super();
	}

	/**
	 * Sets the direction to diagonal
	 * 
	 */
	public void setDiagonal() {
		direction = Directions.DIAGONAL;
		length = 0;
	}

	/**
	 * Sets the direction to up
	 * 
	 * @param length
	 *            length of gap
	 */
	public void setUp(int length) {
		direction = Directions.UP;
		this.length = length;
	}

	/**
	 * Sets the direction to left
	 * 
	 * @param length
	 *            length of gap
	 */
	public void setLeft(int length) {
		direction = Directions.LEFT;
		this.length = length;
	}

	/**
	 * Sets the direction to stop
	 */
	public void setStop() {
		direction = Directions.STOP;
		length = 0;
	}

	/**
	 * @return Returns the direction.
	 */
	public byte getDirection() {
		return direction;
	}

	/**
	 * @return Returns the length.
	 */
	public int getLength() {
		return length;
	}
}