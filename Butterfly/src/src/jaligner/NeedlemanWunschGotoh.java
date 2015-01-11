/*
 * $Id: NeedlemanWunschGotoh.java,v 1.12 2006/04/05 00:18:30 ahmed Exp $
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

import java.util.logging.Logger;

/**
 * An implementation of the Needleman-Wunsch algorithm with Gotoh's improvement
 * for biological global pairwise sequence alignment.
 * 
 * Modified to use int[[] for lengths array for deal with long sequences by Josh Bowden, CSIRO
 * 
 * @author Ahmed Moustafa
 * @author Bram Minnaert
 * @version $Revision: 1.12 $
 */

public final class NeedlemanWunschGotoh {

    /**
     * Logger
     */
    private static final Logger logger = Logger
            .getLogger(NeedlemanWunschGotoh.class.getSimpleName());

    /**
     * Hidden constructor
     */
    private NeedlemanWunschGotoh() {
        super();
    }

    /**
     * Aligns two sequences by Needleman-Wunsch (global)
     * 
     * @param s1
     *            sequene #1 ({@link Read})
     * @param s2
     *            sequene #2 ({@link Read})
     * @param matrix
     *            scoring matrix ({@link Matrix})
     * @param o
     *            open gap penalty
     * @param e
     *            extend gap penalty
     * @return alignment object contains the two aligned sequences, the
     *         alignment score and alignment statistics
     * @see Read
     * @see Matrix
     */
    public static Alignment align(Sequence s1, Sequence s2, Matrix matrix,
            float o, float e) {

        float[][] scores = matrix.getScores();

        Sequence _s1;
        Sequence _s2;

        if (s1.length() < s2.length()) {
            _s1 = s2;
            _s2 = s1;
        } else {
            _s1 = s1;
            _s2 = s2;
        }

        int m = _s1.length() + 1;
        int n = _s2.length() + 1;

        byte[] pointers = new byte[m * n];

        int[] lengths = new int[m * n];  // was short

        // Initializes the element (0,0) of the traceback matrix to STOP.
        pointers[0] = Directions.STOP;

        // Initializes the boundaries of the traceback matrix.
        for (int i = 1, k = n; i < m; i++, k += n) {
            pointers[k] = Directions.UP;
            lengths[k] = (int) i; // was short
        }
        for (int j = 1; j < n; j++) {
            pointers[j] = Directions.LEFT;
            lengths[j] = (int) j;  // was short
        }

        Cell cell = construct(_s1, _s2, scores, o, e, pointers, lengths);

        Alignment alignment = traceback(_s1, _s2, matrix, pointers, cell,
                lengths);

        alignment.setMatrix(matrix);
        alignment.setOpen(o);
        alignment.setExtend(e);
        alignment.setName1(_s1.getId());
        alignment.setName2(_s2.getId());
        alignment.setOriginalSequence1(_s1);
        alignment.setOriginalSequence2(_s2);

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
     * @param o
     *            open gap penalty
     * @param e
     *            extend gap penalty
     * @param pointers
     *            traceback matrix
     * 
     * @return The cell where the traceback starts.
     */
    private static Cell construct(Sequence s1, Sequence s2, float[][] matrix,
            float o, float e, byte[] pointers, int[] lengths) {  // was short[] lengths

        //logger.info("Started...");

        char[] a1 = s1.toArray();
        char[] a2 = s2.toArray();

        int m = s1.length() + 1; // number of rows in similarity matrix
        int n = s2.length() + 1; // number of columns in similarity matrix

        float[] v = new float[n];
        float vDiagonal = 0;// Float.NEGATIVE_INFINITY; // best score in cell
        float f = Float.NEGATIVE_INFINITY; // score from diagonal
        float h = Float.NEGATIVE_INFINITY; // best score ending with gap from
        // left
        float[] g = new float[n]; // best score ending with gap from above

        // Initialization of v and g
        g[0] = Float.NEGATIVE_INFINITY;
        for (int j = 1; j < n; j++) {
            v[j] = 0;// -o - (j - 1) * e;
            g[j] = Float.NEGATIVE_INFINITY;
        }

        int lengthOfHorizontalGap = 0;
        int[] lengthOfVerticalGap = new int[n];

        float similarityScore;
        float maximumScore = Float.NEGATIVE_INFINITY;
        int maxi = 0;
        int maxj = 0;

        // Fill the matrices
        for (int i = 1, k = n; i < m; i++, k += n) { // for all rows
            v[0] = -o - (i - 1) * e;
            for (int j = 1, l = k + 1; j < n; j++, l++) { // for all columns

                similarityScore = matrix[a1[i - 1]][a2[j - 1]];

                f = vDiagonal + similarityScore;// from diagonal

                // Which cell from the left?
                if (h - e >= v[j - 1] - o) {
                    h -= e;
                    lengthOfHorizontalGap++;
                } else {
                    h = v[j - 1] - o;
                    lengthOfHorizontalGap = 1;
                }

                // Which cell from above?
                if (g[j] - e >= v[j] - o) {
                    g[j] = g[j] - e;
                    lengthOfVerticalGap[j] = lengthOfVerticalGap[j] + 1;
                } else {
                    g[j] = v[j] - o;
                    lengthOfVerticalGap[j] = 1;
                }

                vDiagonal = v[j];
                v[j] = maximum(f, g[j], h); // best one
                if (v[j] > maximumScore) {
                    maximumScore = v[j];
                    maxi = i;
                    maxj = j;
                }

                // Determine the traceback direction
                if (v[j] == f) {
                    pointers[l] = Directions.DIAGONAL;
                } else if (v[j] == g[j]) {
                    pointers[l] = Directions.UP;
                    lengths[l] = (int) lengthOfVerticalGap[j]; // was short
                } else if (v[j] == h) {
                    pointers[l] = Directions.LEFT;
                    lengths[l] = (int) lengthOfHorizontalGap;  // was short
                }

            } // loop columns

            // Reset
            h = Float.NEGATIVE_INFINITY;
            vDiagonal = 0;// -o - (i - 1) * e;

            lengthOfHorizontalGap = 0;

        } // loop rows

        Cell cell = new Cell();
        cell.set(maxi, maxj, v[n - 1]);

        //logger.info("Finished.");

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
            byte[] pointers, Cell cell, int[] lengths) { // was short[] lengths
        //logger.info("Started...");

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
        int n = s2.length() + 1;
        int row = i * n;

        int a = s1.length() - 1;
        int b = s2.length() - 1;
        if (a - i > b - j) {
            for (; a - i > b - j; a--) {
                reversed1[len1++] = array1[a];
                reversed2[len2++] = Alignment.GAP;
                reversed3[len3++] = Markups.GAP;
                gaps++;
            }
            for (; b > j - 1; a--, b--) {
                c1 = array1[a];
                c2 = array2[b];

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
            }
        } else {
            for (; b - j > a - i; b--) {
                reversed1[len1++] = Alignment.GAP;
                reversed2[len2++] = array2[b];
                reversed3[len3++] = Markups.GAP;
                gaps++;
            }
            for (; a > i - 1; a--, b--) {
                c1 = array1[a];
                c2 = array2[b];

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
            }
        }

        // Traceback flag, where true => continue and false => stop
        boolean stillGoing = true;
        while (stillGoing) {
            int l = row + j;
            switch (pointers[l]) {
            case Directions.UP:
                for (int k = 0, len = lengths[l]; k < len; k++) {
                    reversed1[len1++] = array1[--i];
                    reversed2[len2++] = Alignment.GAP;
                    reversed3[len3++] = Markups.GAP;
                    row -= n;
                    gaps++;
                }
                break;
            case Directions.DIAGONAL:
                c1 = array1[--i];
                c2 = array2[--j];
                reversed1[len1++] = c1;
                reversed2[len2++] = c2;
                row -= n;
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
                for (int k = 0, len = lengths[l]; k < len; k++) {
                    reversed1[len1++] = Alignment.GAP;
                    reversed2[len2++] = array2[--j];
                    reversed3[len3++] = Markups.GAP;
                    gaps++;
                }
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

        //logger.info("Finished.");

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