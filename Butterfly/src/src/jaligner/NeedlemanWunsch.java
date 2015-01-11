/*
 * $Id: NeedlemanWunsch.java,v 1.4 2006/02/09 12:23:24 ahmed Exp $
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

import jaligner.matrix.Matrix;

/**
 * An implementation of the Needleman-Wunsch algorithm for biological global
 * pairwise sequence alignment.
 * 
 * <br>
 * Reference: <a
 * href="http://www.sbc.su.se/~per/molbioinfo2001/dynprog/adv_dynamic.html">Advanced
 * Dynamic Programming Tutorial</a>.
 * 
 * 
 * @author <a href="ahmed@users.sf.net">Ahmed Moustafa</a>
 */

public final class NeedlemanWunsch {
    /**
     * Hidden constructor
     */
    private NeedlemanWunsch() {
        super();
    }

    /**
     * Aligns two sequences by Needleman-Wunsch (global)
     * 
     * @param s1
     *            sequene #1 ({@link Sequence})
     * @param s2
     *            sequene #2 ({@link Sequence})
     * @param matrix
     *            scoring matrix ({@link Matrix})
     * @param gap
     *            open gap penalty
     * @return alignment object contains the two aligned sequences, the
     *         alignment score and alignment statistics
     * @see Sequence
     * @see Matrix
     */
    public static Alignment align(Sequence s1, Sequence s2, Matrix matrix,
            float gap) {
        float[][] scores = matrix.getScores();

        int m = s1.length() + 1;
        int n = s2.length() + 1;

        byte[][] pointers = new byte[m][n];

        // Initializes the element (0,0) of the traceback matrix to STOP.
        pointers[0][0] = Directions.STOP;

        // Initializes the boundaries of the traceback matrix.
        for (int i = 1; i < m; i++) {
            pointers[i][0] = Directions.UP;
        }
        for (int j = 1; j < n; j++) {
            pointers[0][j] = Directions.LEFT;
        }

        Cell cell = construct(s1, s2, scores, gap, pointers);
        Alignment alignment = traceback(s1, s2, matrix, pointers, cell);

        alignment.setOriginalSequence1(s1);
        alignment.setOriginalSequence2(s2);
        alignment.setMatrix(matrix);
        alignment.setOpen(gap);
        alignment.setExtend(gap);

        if (s1.getId() != null) {
            alignment.setName1(s1.getId());
        }
        if (s2.getId() != null) {
            alignment.setName2(s2.getId());
        }

        return alignment;
    }

    /**
     * Constructs directions matrix for the traceback.
     * 
     * @param s1
     *            sequence #1
     * @param s2
     *            sequence #2
     * @param matrix
     *            scoring matrix
     * @param gap
     *            gap penalty
     * @param pointers
     *            traceback matrix
     * 
     * @return The cell where the traceback starts.
     */
    private static Cell construct(Sequence s1, Sequence s2, float[][] matrix,
            float gap, byte[][] pointers) {

        Cell cell = new Cell();

        char[] a1 = s1.toArray();
        char[] a2 = s2.toArray();

        int m = s1.length() + 1; // number of rows in similarity matrix
        int n = s2.length() + 1; // number of columns in similarity matrix

        // float[] v = new float[n]; // optimal alignment
        float[] v = new float[n];

        // Initialization of v
        for (int j = 1; j < n; j++) {
            v[j] = j * -gap;
        }
        v[0] = 0;

        float x, y, z;
        float vOld = 0;

        // Fill the matrices
        for (int i = 1; i < m; i++) { // for all rows
            v[0] = i * -gap;

            for (int j = 1; j < n; j++) { // for all columns

                x = v[j] - gap;
                y = v[j - 1] - gap;
                z = vOld + matrix[a1[i - 1]][a2[j - 1]];

                vOld = v[j];
                v[j] = maximum(x, y, z);

                // Determine the traceback direction
                if (v[j] == x) {
                    pointers[i][j] = Directions.UP;
                } else if (v[j] == y) {
                    pointers[i][j] = Directions.LEFT;
                } else {
                    pointers[i][j] = Directions.DIAGONAL;
                }

            } // loop columns

            vOld = i * -gap;
        } // loop rows

        // cell contains the row number, the column number
        // and the score of the cell with the maximum score

        // Set the traceback start at the last cell m, n
        // because we are doing global alignment
        cell.set(m - 1, n - 1, v[n - 1]);
        return cell;
    }

    /**
     * Returns the alignment of two sequences based on the passed array of
     * pointers
     * 
     * @param s1
     *            sequence #1
     * @param s2
     *            sequence #2
     * @param m
     *            scoring matrix
     * @param pointers
     *            traceback matrix
     * @param cell
     *            The cell where the traceback starts.
     * @return {@link Alignment} with the two aligned sequences and alignment
     *         score.
     * @see Cell
     * @see Alignment
     */
    private static Alignment traceback(Sequence s1, Sequence s2, Matrix m,
            byte[][] pointers, Cell cell) {
        char[] array1 = s1.toArray();
        char[] array2 = s2.toArray();
        float[][] scores = m.getScores();

        Alignment alignment = new Alignment();
        alignment.setScore(cell.getScore());

        // maximum length after the aligned sequences
        int maxlen = s1.length() + s2.length();

        char[] reversed1 = new char[maxlen]; // reversed sequence #1
        char[] reversed2 = new char[maxlen]; // reversed sequence #2
        char[] reversed3 = new char[maxlen]; // reversed markup

        int len1 = 0; // length of sequence #1 after alignment
        int len2 = 0; // length of sequence #2 after alignment
        int len3 = 0; // length of the markup line

        int identity = 0; // count of identitcal pairs
        int similarity = 0; // count of similar pairs
        int gaps = 0; // count of gaps

        char c1, c2;

        int i = cell.getRow(); // traceback start row
        int j = cell.getCol(); // traceback start col

        // Traceback flag, where true => continue and false => stop
        boolean stillGoing = true;
        while (stillGoing) {
            switch (pointers[i][j]) {
            case Directions.UP:
                reversed1[len1++] = array1[--i];
                reversed2[len2++] = Alignment.GAP;
                reversed3[len3++] = Markups.GAP;
                gaps++;
                break;
            case Directions.DIAGONAL:
                c1 = array1[--i];
                c2 = array2[--j];
                reversed1[len1++] = c1;
                reversed2[len2++] = c2;
                if (c1 == c2) {
                    reversed3[len3++] = Markups.IDENTITY;
                    identity++;
                    similarity++;
                } else if (scores[c1][c2] > 0) {
                    reversed3[len3++] = Markups.SIMILARITY;
                    similarity++;
                } else {
                    reversed3[len3++] = Markups.MISMATCH;
                }
                break;
            case Directions.LEFT:
                reversed1[len1++] = Alignment.GAP;
                reversed2[len2++] = array2[--j];
                reversed3[len3++] = Markups.GAP;
                gaps++;
                break;
            case Directions.STOP:
                stillGoing = false;
            }
        }

        alignment.setSequence1(reverse(reversed1, len1));
        alignment.setStart1(i);
        alignment.setSequence2(reverse(reversed2, len2));
        alignment.setStart2(j);
        alignment.setMarkupLine(reverse(reversed3, len3));
        alignment.setIdentity(identity);
        alignment.setGaps(gaps);
        alignment.setSimilarity(similarity);

        return alignment;
    }

    /**
     * Returns the maximum of two float numbers.
     * 
     * @param a
     *            float #1
     * @param b
     *            float #2
     * @param c
     *            float #3
     * @return the maximum of a and b
     */
    private static float maximum(float a, float b, float c) {
        if (a > b) {
            return a > c ? a : c;
        } else {
            return b > c ? b : c;
        }
    }

    /**
     * Reverses an array of chars
     * 
     * @param a
     * @param len
     * @return the input array of char reserved
     */
    private static char[] reverse(char[] a, int len) {
        char[] b = new char[len];
        for (int i = len - 1, j = 0; i >= 0; i--, j++) {
            b[j] = a[i];
        }
        return b;
    }
}