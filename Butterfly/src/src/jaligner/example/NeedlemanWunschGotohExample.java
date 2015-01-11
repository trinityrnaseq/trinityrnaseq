/*
 * $Id: NeedlemanWunschGotohExample.java,v 1.1 2006/02/03 05:00:23 ahmed Exp $
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

package jaligner.example;

import jaligner.Alignment;
import jaligner.NeedlemanWunschGotoh;
import jaligner.Sequence;
import jaligner.formats.Pair;
import jaligner.matrix.Matrix;
import jaligner.matrix.MatrixGenerator;
import jaligner.util.SequenceParser;

import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Example of using JAligner API to perform global pairwise sequence alignment
 * with {@link jaligner.NeedlemanWunschGotoh}.
 * 
 * @author Ahmed Moustafa (ahmed@users.sf.net)
 */

public class NeedlemanWunschGotohExample {

    /**
     * Logger
     */
    private static final Logger logger = Logger
            .getLogger(NeedlemanWunschGotohExample.class.getName());

    /**
     * 
     * @param args
     */
    public static void main(String[] args) {
        try {
            logger.info("Running example...");

            Sequence s1 = SequenceParser.parse("GAATTCAGTTA");
            Sequence s2 = SequenceParser.parse("GGATCGA");

            float match = 2;
            float mismatch = -1;
            Matrix matrix = MatrixGenerator.generate(match, mismatch);
            float open = 2;
            float extend = 2;

            Alignment alignment = NeedlemanWunschGotoh.align(s1, s2, matrix,
                    open, extend);

            System.out.println(alignment.getSummary());
            System.out.println(new Pair().format(alignment));

            logger.info("Finished running example");
        } catch (Exception e) {
            logger.log(Level.SEVERE, "Failed running example: "
                    + e.getMessage(), e);
        }
    }
}