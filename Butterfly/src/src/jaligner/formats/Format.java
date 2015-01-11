/*
 * $Id: Format.java,v 1.1 2005/05/25 19:56:30 ahmed Exp $
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
 * Abstract format
 * 
 * @author Ahmed Moustafa (ahmed@users.sf.net)
 */

public abstract class Format {
	
	/**
	 * Format id
	 */
	private String id = null;
	
	/**
	 * Formats alignment
	 * @param alignment
	 * @return formatted alignment
	 * @see Alignment
	 */
	public abstract String format(Alignment alignment);
	
	/**
	 * Sets format id
	 * @param id to set
	 */
	public void setId (String id) {
		this.id =  id;
	}
	
	/**
	 * Returns format id
	 * @return id
	 */
	public String getId ( ) {
		return this.id == null ? this.getClass().getName() : this.id;
	}
}