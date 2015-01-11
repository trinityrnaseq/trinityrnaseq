/*
 * $Id: SmithWatermanGotoh.java,v 1.10 2006/02/09 13:27:36 ahmed Exp $
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
 * An implementation of the Smith-Waterman algorithm with Gotoh's improvement
 * for biological local pairwise sequence alignment.
 * 
 * <strong>Recursive definition:</strong>
 * <ul>
 * <li>
 * <strong>Base conditions:</strong>
 * <ul>
 * <li><code>V(0, 0) = 0</code></li>
 * <li><code>V(i, 0) = E(i, 0) = W<sub>g</sub> + iW<sub>s</sub></code></li>
 * <li><code>V(0, j) = F(0, j) = W<sub>g</sub> + jW<sub>s</sub></code></li>
 * </ul>
 * </li>
 * <li>
 * <strong>Recurrence relation:</strong>
 * <ul>
 * <li><code>V(i, j) = max{E(i, j), F(i, j), G(i, j)}</code>, where:</li>
 * <li><code>G(i, j) = V(i - 1, j - 1) + similarity(S<sub>i</sub>, T<sub>j</sub>)</code></li>
 * <li><code>E(i, j) = max{E(i, j - 1) + W<sub>s</sub>, V(i, j - 1) + W<sub>g</sub> + W<sub>s</sub>}</code></li>
 * <li><code>F(i, j) = max{F(i - 1, j) + W<sub>s</sub>, V(i - 1, j) + W<sub>g</sub> + W<sub>s</sub>}</code></li>
 * </ul>
 * </ul> 
 * 
 * @author Ahmed Moustafa (ahmed@users.sf.net)
 */

public final class SmithWatermanGotoh {
	/**
	 * Hidden constructor
	 */
	private SmithWatermanGotoh() {
		super();
	}

	/**
	 * Logger
	 */
	private static final Logger logger = Logger
			.getLogger(SmithWatermanGotoh.class.getName());

	/**
	 * Aligns two sequences by Smith-Waterman (local)
	 * 
	 * @param s1
	 *            sequene #1 ({@link Sequence})
	 * @param s2
	 *            sequene #2 ({@link Sequence})
	 * @param matrix
	 *            scoring matrix ({@link Matrix})
	 * @param o
	 *            open gap penalty
	 * @param e
	 *            extend gap penalty
	 * @return alignment object contains the two aligned sequences, the
	 *         alignment score and alignment statistics
	 * @see Sequence
	 * @see Matrix
	 */
	public static Alignment align(Sequence s1, Sequence s2, Matrix matrix,
			float o, float e) {
		//logger.info("Started...");
		long start = System.currentTimeMillis();
		float[][] scores = matrix.getScores();

		int m = s1.length() + 1;
		int n = s2.length() + 1;

		byte[] pointers = new byte[m * n];

		// Initializes the boundaries of the traceback matrix to STOP.
		for (int i = 0, k = 0; i < m; i++, k += n) {
			pointers[k] = Directions.STOP;
		}
		for (int j = 1; j < n; j++) {
			pointers[j] = Directions.STOP;
		}

		short[] sizesOfVerticalGaps = new short[m * n];
		short[] sizesOfHorizontalGaps = new short[m * n];
		for (int i = 0, k = 0; i < m; i++, k += n) {
			for (int j = 0; j < n; j++) {
				sizesOfVerticalGaps[k + j] = sizesOfHorizontalGaps[k + j] = 1;
			}
		}

		Cell cell = SmithWatermanGotoh.construct(s1, s2, scores, o, e,
				pointers, sizesOfVerticalGaps, sizesOfHorizontalGaps);
		Alignment alignment = SmithWatermanGotoh.traceback(s1, s2, matrix,
				pointers, cell, sizesOfVerticalGaps, sizesOfHorizontalGaps);
		alignment.setOriginalSequence1(s1);
		alignment.setOriginalSequence2(s2);
		alignment.setMatrix(matrix);
		alignment.setOpen(o);
		alignment.setExtend(e);
		if (s1.getId() != null) {
			alignment.setName1(s1.getId());
		}
		if (s2.getId() != null) {
			alignment.setName2(s2.getId());
		}
		/*logger.info("Finished in " + (System.currentTimeMillis() - start)
				+ " milliseconds");
		*/
		
		return alignment;
	}

	/**
	 * Constructs directions matrix for the traceback
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
	 * @return The cell where the traceback starts.
	 */
	private static Cell construct(Sequence s1, Sequence s2, float[][] matrix,
			float o, float e, byte[] pointers, short[] sizesOfVerticalGaps,
			short[] sizesOfHorizontalGaps) {
		//logger.info("Started...");
		long start = System.currentTimeMillis();

		char[] a1 = s1.toArray();
		char[] a2 = s2.toArray();

		int m = s1.length() + 1;
		int n = s2.length() + 1;

		float f; // score of alignment x1...xi to y1...yi if xi aligns to yi
		float[] g = new float[n]; // score if xi aligns to a gap after yi
		float h; // score if yi aligns to a gap after xi
		float[] v = new float[n]; // best score of alignment x1...xi to
		// y1...yi
		float vDiagonal;

		g[0] = Float.NEGATIVE_INFINITY;
		h = Float.NEGATIVE_INFINITY;
		v[0] = 0;

		for (int j = 1; j < n; j++) {
			g[j] = Float.NEGATIVE_INFINITY;
			v[j] = 0;
		}

		float similarityScore, g1, g2, h1, h2;

		Cell cell = new Cell();

		for (int i = 1, k = n; i < m; i++, k += n) {
			h = Float.NEGATIVE_INFINITY;
			vDiagonal = v[0];
			for (int j = 1, l = k + 1; j < n; j++, l++) {
				similarityScore = matrix[a1[i - 1]][a2[j - 1]];

				// Fill the matrices
				f = vDiagonal + similarityScore;

				g1 = g[j] - e;
				g2 = v[j] - o;
				if (g1 > g2) {
					g[j] = g1;
					sizesOfVerticalGaps[l] = (short) (sizesOfVerticalGaps[l - n] + 1);
				} else {
					g[j] = g2;
				}

				h1 = h - e;
				h2 = v[j - 1] - o;
				if (h1 > h2) {
					h = h1;
					sizesOfHorizontalGaps[l] = (short) (sizesOfHorizontalGaps[l - 1] + 1);
				} else {
					h = h2;
				}

				vDiagonal = v[j];
				v[j] = maximum(f, g[j], h, 0);

				// Determine the traceback direction
				if (v[j] == 0) {
					pointers[l] = Directions.STOP;
				} else if (v[j] == f) {
					pointers[l] = Directions.DIAGONAL;
				} else if (v[j] == g[j]) {
					pointers[l] = Directions.UP;
				} else {
					pointers[l] = Directions.LEFT;
				}

				// Set the traceback start at the current cell i, j and score
				if (v[j] > cell.getScore()) {
					cell.set(i, j, v[j]);
				}
			}
		}
		/*logger.info("Finished in " + (System.currentTimeMillis() - start)
				+ " milliseconds");
		 */	
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
	 * @param cell
	 *            The cell where the traceback starts.
	 * @return {@link Alignment}with the two aligned sequences and alignment
	 *         score.
	 * @see Cell
	 * @see Alignment
	 */
	private static Alignment traceback(Sequence s1, Sequence s2, Matrix m,
			byte[] pointers, Cell cell, short[] sizesOfVerticalGaps,
			short[] sizesOfHorizontalGaps) {
		//logger.info("Started...");
		long start = System.currentTimeMillis();

		char[] a1 = s1.toArray();
		char[] a2 = s2.toArray();

		float[][] scores = m.getScores();

		int n = s2.length() + 1;

		Alignment alignment = new Alignment();
		alignment.setScore(cell.getScore());

		int maxlen = s1.length() + s2.length(); // maximum length after the
		// aligned sequences

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
		int k = i * n;

		boolean stillGoing = true; // traceback flag: true -> continue & false
		// -> stop

		while (stillGoing) {
			switch (pointers[k + j]) {
			case Directions.UP:
				for (int l = 0, len = sizesOfVerticalGaps[k + j]; l < len; l++) {
					reversed1[len1++] = a1[--i];
					reversed2[len2++] = Alignment.GAP;
					reversed3[len3++] = Markups.GAP;
					k -= n;
					gaps++;
				}
				break;
			case Directions.DIAGONAL:
				c1 = a1[--i];
				c2 = a2[--j];
				k -= n;
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
				for (int l = 0, len = sizesOfHorizontalGaps[k + j]; l < len; l++) {
					reversed1[len1++] = Alignment.GAP;
					reversed2[len2++] = a2[--j];
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

		/*logger.info("Finished in " + (System.currentTimeMillis() - start)
				+ " milliseconds");
		*/
		
		return alignment;
	}

	/**
	 * Returns the maximum of 4 float numbers.
	 * 
	 * @param a
	 *            float #1
	 * @param b
	 *            float #2
	 * @param c
	 *            float #3
	 * @param d
	 *            float #4
	 * @return The maximum of a, b, c and d.
	 */
	private static float maximum(float a, float b, float c, float d) {
		if (a > b) {
			if (a > c) {
				return a > d ? a : d;
			} else {
				return c > d ? c : d;
			}
		} else if (b > c) {
			return b > d ? b : d;
		} else {
			return c > d ? c : d;
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