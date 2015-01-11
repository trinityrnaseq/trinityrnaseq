/*
 * $Id: FormatFactory.java,v 1.1 2005/05/25 19:56:30 ahmed Exp $
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

import java.util.Collection;
import java.util.HashMap;

/**
 * Formats factory.
 * 
 * @author Ahmed Moustafa (ahmed@users.sf.net)
 */

public class FormatFactory {
    /**
     * Instance of {@link FormatFactory}
     */
	private static FormatFactory instance = null;
    
    /**
     * {@link HashMap} of {@link Format}
     */
	private HashMap formats = new HashMap();
	
	/**
	 * Hidden constructor
	 *
	 */
	private FormatFactory ( ) {
		super();
	}
	
	/**
     * Returns an instance for {@link FormatFactory}.
     * @return {@link FormatFactory}
     */
    public static FormatFactory getInstance() {
    	if (instance == null) {
    		instance = new FormatFactory();
    	}
    	return instance;
    }
    
    /**
     * Registers format.
     * @param format instance of format
     */
    public void registerFormat(Format format) {
    	formats.put(format.getId(), format);
    }
    
    /**
     * Returns an instance of {@link Format}.
     * @param id format id
     * @return {@link Format} or null if id not found
     */
    public Format getFormat(String id) {
    	return (Format) formats.get(id);
    }
    
    /**
     * Returns a list of registered formats
     * @return {@link Collection}
     */
    public Collection getFormats( ) {
    	return formats.keySet();
    }
}