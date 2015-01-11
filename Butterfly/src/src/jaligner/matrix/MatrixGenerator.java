/*
 * $Id: MatrixGenerator.java,v 1.3 2006/02/26 20:26:41 ahmed Exp $
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

/**
 * A utility class to generate a scoring matrix from a fixed match and mismatch
 * scoring scheme.
 * 
 * @author Ahmed Moustafa (ahmed@users.sf.net)
 */

public final class MatrixGenerator {

    /**
     * Returns scoring matrix from a fixed match/mismatch scoring scheme.
     * 
     * @param match
     *            match score
     * @param mismatch
     *            mistmatch score
     * 
     * @return scoring matrix
     */
    public static Matrix generate(float match, float mismatch) {
        float[][] scores = new float[Matrix.SIZE][Matrix.SIZE];

        // Fill the matrix with the scores
        for (int i = 0; i < Matrix.SIZE; i++) {
            for (int j = 0; j < Matrix.SIZE; j++) {
                if (i == j || i == 'N' || j == 'N') {
                    scores[i][j] = match;
                } else {
                    scores[i][j] = mismatch;
                }
            }
        }

        // Generate some id for the matrix (hopefully to be somehow unique)
        String id = new Long(System.currentTimeMillis()).toString();

        Matrix matrix = new Matrix(id, scores);

        return matrix;
    }
}
