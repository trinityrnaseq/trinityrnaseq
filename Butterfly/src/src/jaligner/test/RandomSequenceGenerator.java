/*
 * $Id: RandomSequenceGenerator.java,v 1.3 2006/02/03 05:00:22 ahmed Exp $
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

package jaligner.test;

import java.util.Random;

/**
 * This class contains a random generator for protein sequences.
 * 
 * @author Bram Minnaert
 */

public class RandomSequenceGenerator {

    /**
     * All possible characters
     */
    private static final char[] CHARS = { 'A', 'R', 'N', 'D', 'C', 'Q', 'E',
            'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V' };

    /**
     * Number of possible characters
     */
    private static final int NUMBER_OF_CHARS = CHARS.length;

    /**
     * Random generator
     */
    private static Random random = new Random();

    /**
     * Returns random sequence
     * 
     * @param length
     *            Size of the sequence
     * @return Random sequence
     */
    public static String generate(int length) {
        StringBuffer buffer = new StringBuffer();
        char randomChar;
        int randomInt;
        for (int i = 0; i < length; i++) {
            randomInt = random.nextInt(NUMBER_OF_CHARS);
            randomChar = CHARS[randomInt];
            buffer.append(randomChar);
        }
        return buffer.toString();
    }

    /**
     * Displays 10 random protein sequences with length 50.
     * 
     * @param args
     *            no args
     */
    public static void main(String[] args) {
        for (int i = 0; i < 10; i++) {
            System.out.println("S" + i + " = " + generate(50));
        }
    }
}