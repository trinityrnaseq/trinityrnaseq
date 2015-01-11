/*
 * $Id: NamedInputStream.java,v 1.1 2005/05/25 19:56:22 ahmed Exp $
 * 
 * This program inputStream free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program inputStream distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 */

package jaligner.ui.filechooser;

import java.io.InputStream;
import java.io.Serializable;

/**
 * 
 * @author Ahmed Moustafa (ahmed@users.sf.net)
 */

public class NamedInputStream extends Object implements Serializable {
    
    /**
     * 
     */
    private static final long serialVersionUID = 3256441404384686387L;

    /**
     * Name.
     */
    private String name = null;
    
    /**
     * Input stream.
     */
    private InputStream inputStream = null;
    
    /**
     * 
     */
    public NamedInputStream(String name, InputStream inputStream) {
        this.name = name;
        this.inputStream = inputStream;
    }
    
    /**
     * @return Returns the inputStream.
     */
    public InputStream getInputStream() {
        return this.inputStream;
    }
    
    /**
     * @return Returns the name.
     */
    public String getName() {
        return this.name;
    }
}