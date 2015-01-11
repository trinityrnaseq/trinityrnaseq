/*
 * $Id: Alignment.java,v 1.11 2006/08/19 00:42:35 ahmed Exp $
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
import jaligner.util.Commons;

import java.text.DecimalFormat;

/**
 * Holds the output of a pairwise sequences alignment.
 * 
 * @author Ahmed Moustafa (ahmed@users.sf.net)
 */

public final class Alignment {

	/**
	 * Gap character
	 */
	public static final char GAP = '-';

	/**
	 * Default name for sequence #1
	 */
	private static final String SEQUENCE1 = "jaligner_1";

	/**
	 * Default name for sequence #2
	 */
	private static final String SEQUENCE2 = "jaligner_2";

	/**
	 * Scoring matrix
	 */
	private Matrix matrix;

	/**
	 * Gap open cost
	 */
	private float open;

	/**
	 * Gap extend cost
	 */
	private float extend;

	/**
	 * Alignment score
	 */
	private float score;

	/**
	 * Aligned sequence #1
	 */
	private char[] sequence1;

	/**
	 * Name of sequence #1
	 */
	private String name1;

	/**
	 * Alignment start location in sequence #1
	 */
	private int start1;

	/**
	 * Aligned sequence #2
	 */
	private char[] sequence2;

	/**
	 * Name of sequence #2
	 */
	private String name2;

	/**
	 * Alignment start location in sequence #2
	 */
	private int start2;

	/**
	 * Markup line
	 */
	private char[] markupLine;

	/**
	 * Count of identical locations
	 */
	private int identity;

	/**
	 * Count of similar locations
	 */
	private int similarity;

	/**
	 * Count of gap locations
	 */
	private int gaps;

	private Sequence originalSequence1;

	private Sequence originalSequence2;

	/**
	 * Constructor for Alignment
	 */

	public Alignment() {
		super();
	}

	/**
	 * @return Returns the extend.
	 */
	public float getExtend() {
		return extend;
	}

	/**
	 * @param extend
	 *            The extend to set.
	 */
	public void setExtend(float extend) {
		this.extend = extend;
	}

	/**
	 * @return Returns the matrix.
	 */
	public Matrix getMatrix() {
		return matrix;
	}

	/**
	 * @param matrix
	 *            The matrix to set.
	 */
	public void setMatrix(Matrix matrix) {
		this.matrix = matrix;
	}

	/**
	 * @return Returns the name1.
	 */
	public String getName1() {
		return name1 == null || name1.trim().length() == 0 ? SEQUENCE1 : name1;
	}

	/**
	 * @param name1
	 *            The name1 to set.
	 */
	public void setName1(String name1) {
		this.name1 = name1;
	}

	/**
	 * @return Returns the name2.
	 */
	public String getName2() {
		return name2 == null || name2.trim().length() == 0 ? SEQUENCE2 : name2;
	}

	/**
	 * @param name2
	 *            The name2 to set.
	 */
	public void setName2(String name2) {
		this.name2 = name2;
	}

	/**
	 * @return Returns the open.
	 */
	public float getOpen() {
		return open;
	}

	/**
	 * @param open
	 *            The open to set.
	 */
	public void setOpen(float open) {
		this.open = open;
	}

	/**
	 * @return Returns the score.
	 */
	public float getScore() {
		return score;
	}

	/**
	 * @param score
	 *            The score to set.
	 */
	public void setScore(float score) {
		this.score = score;
	}

	/**
	 * Returns the length of the alignment
	 * 
	 * @return Alignment length
	 */
	public int getLength() {
		return this.sequence1.length;
	}

	/**
	 * @return Returns the sequence1.
	 */
	public char[] getSequence1() {
		return sequence1;
	}

	/**
	 * @param sequence1
	 *            The sequence1 to set.
	 */
	public void setSequence1(char[] sequence1) {
		this.sequence1 = sequence1;
	}

	/**
	 * @return Returns the sequence2.
	 */
	public char[] getSequence2() {
		return sequence2;
	}

	/**
	 * @param sequence2
	 *            The sequence2 to set.
	 */
	public void setSequence2(char[] sequence2) {
		this.sequence2 = sequence2;
	}

	/**
	 * @return Returns the start1.
	 */
	public int getStart1() {
		return start1;
	}

	/**
	 * @param start1
	 *            The start1 to set.
	 */
	public void setStart1(int start1) {
		this.start1 = start1;
	}

	/**
	 * @return Returns the start2.
	 */
	public int getStart2() {
		return start2;
	}

	/**
	 * @param start2
	 *            The start2 to set.
	 */
	public void setStart2(int start2) {
		this.start2 = start2;
	}

	/**
	 * @return Returns the gaps.
	 */
	public int getGaps() {
		return gaps;
	}

	/**
	 * @param gaps
	 *            The gaps to set.
	 */
	public void setGaps(int gaps) {
		this.gaps = gaps;
	}

	/**
	 * @return Returns the identity.
	 */
	public int getIdentity() {
		return identity;
	}

	/**
	 * @param identity
	 *            The identity to set.
	 */
	public void setIdentity(int identity) {
		this.identity = identity;
	}

	/**
	 * @return Returns the markupLine.
	 */
	public char[] getMarkupLine() {
		return markupLine;
	}

	/**
	 * @param markupLine
	 *            The markupLine to set.
	 */
	public void setMarkupLine(char[] markupLine) {
		this.markupLine = markupLine;
	}

	/**
	 * @return Returns the similarity.
	 */
	public int getSimilarity() {
		return similarity;
	}

	/**
	 * @param similarity
	 *            The similarity to set.
	 */
	public void setSimilarity(int similarity) {
		this.similarity = similarity;
	}

	/**
	 * Returns a summary for alignment
	 * 
	 * @return {@link String} alignment summary
	 */
	public String getSummary() {
		StringBuffer buffer = new StringBuffer();
		DecimalFormat f1 = new DecimalFormat("0.00");
		DecimalFormat f2 = new DecimalFormat("0.00%");

		int length = getSequence1().length;

		buffer.append("Sequence #1: " + getName1());
		buffer.append(Commons.getLineSeparator());
		buffer.append("Sequence #2: " + getName2());
		buffer.append(Commons.getLineSeparator());
		buffer.append("Length #1: " + getOriginalSequence1().length());
		buffer.append(Commons.getLineSeparator());
		buffer.append("Length #2: " + getOriginalSequence2().length());
		buffer.append(Commons.getLineSeparator());
		buffer.append("Matrix: "
				+ (matrix.getId() == null ? "" : matrix.getId()));
		buffer.append(Commons.getLineSeparator());
		buffer.append("Gap open: " + open);
		buffer.append(Commons.getLineSeparator());
		buffer.append("Gap extend: " + extend);
		buffer.append(Commons.getLineSeparator());
		buffer.append("Length: " + length);
		buffer.append(Commons.getLineSeparator());
		buffer.append("Identity: " + identity + "/" + length + " ("
				+ f2.format(identity / (float) length) + ")");
		buffer.append(Commons.getLineSeparator());
		buffer.append("Similarity: " + similarity + "/" + length + " ("
				+ f2.format(similarity / (float) length) + ")");
		buffer.append(Commons.getLineSeparator());
		buffer.append("Gaps: " + gaps + "/" + length + " ("
				+ f2.format(gaps / (float) length) + ")");
		buffer.append(Commons.getLineSeparator());
		buffer.append("Score: " + f1.format(score));
		buffer.append(Commons.getLineSeparator());

		return buffer.toString();
	}

	/**
	 * Calculate the score of the alignment, not using the score field (the
	 * function only uses sequence1, sequence2, matrix and gap penalties).
	 * 
	 * @return the calculated score (By: Bram Minnaert)
	 */
	public float calculateScore() {
		// The calculated score
		float calcScore = 0;

		// In the previous step there was a gap in the first sequence
		boolean previous1wasGap = false;

		// In the previous step there was a gap in the second sequence
		boolean previous2wasGap = false;

		int start = 0;
		int end = sequence1.length - 1;

		char c1, c2; // the next character
		for (int i = start; i <= end; i++) {
			c1 = sequence1[i];
			c2 = sequence2[i];
			// the next character in the first sequence is a gap
			if (c1 == GAP) {
				if (previous1wasGap) {
					calcScore -= extend;
				} else {
					calcScore -= open;
				}
				previous1wasGap = true;
				previous2wasGap = false;
			}
			// the next character in the second sequence is a gap
			else if (c2 == GAP) {
				if (previous2wasGap) {
					calcScore -= extend;
				} else {
					calcScore -= open;
				}
				previous1wasGap = false;
				previous2wasGap = true;
			}
			// the next characters in boths sequences are not gaps
			else {
				calcScore += matrix.getScore(c1, c2);
				previous1wasGap = false;
				previous2wasGap = false;
			}
		}
		return calcScore;
	}

	/**
	 * Calculate the score of the alignment without the terminal gaps.
	 * 
	 */
	public float getScoreWithNoTerminalGaps() {
		// The calculated score
		float calcScore = 0;

		// In the previous step there was a gap in the first sequence
		boolean previous1wasGap = false;

		// In the previous step there was a gap in the second sequence
		boolean previous2wasGap = false;

		int start = 0;
		int end = sequence1.length - 1;

		if (sequence1[start] == GAP) {
			while (sequence1[start] == GAP) {
				start++;
			}
		} else if (sequence2[start] == GAP) {
			while (sequence2[start] == GAP) {
				start++;
			}
		}

		if (sequence1[end] == GAP) {
			while (sequence1[end] == GAP) {
				end--;
			}
		} else if (sequence2[end] == GAP) {
			while (sequence2[end] == GAP) {
				end--;
			}
		}

		char c1, c2; // the next character
		for (int i = start; i <= end; i++) {
			c1 = sequence1[i];
			c2 = sequence2[i];
			// the next character in the first sequence is a gap
			if (c1 == GAP) {
				if (previous1wasGap) {
					calcScore -= extend;
				} else {
					calcScore -= open;
				}
				previous1wasGap = true;
				previous2wasGap = false;
			}
			// the next character in the second sequence is a gap
			else if (c2 == GAP) {
				if (previous2wasGap) {
					calcScore -= extend;
				} else {
					calcScore -= open;
				}
				previous1wasGap = false;
				previous2wasGap = true;
			}
			// the next characters in boths sequences are not gaps
			else {
				calcScore += matrix.getScore(c1, c2);
				previous1wasGap = false;
				previous2wasGap = false;
			}
		}
		return calcScore;
	}

	/**
	 * Check if the calculated score matches the field score.
	 * 
	 * @return true if equal, else false. (By: Bram Minnaert)
	 */
	public boolean checkScore() {
		if (calculateScore() == score) {
			return true;
		} else {
			return false;
		}
	}

	/**
	 * Returns original {@link Sequence} #1
	 * 
	 * @return original {@link Sequence} #1
	 */
	public Sequence getOriginalSequence1() {
		return originalSequence1;
	}

	/**
	 * 
	 * @param originalSequence1
	 */
	public void setOriginalSequence1(Sequence originalSequence1) {
		this.originalSequence1 = originalSequence1;
	}

	/**
	 * Returns original {@link Sequence} #2
	 * 
	 * @return original {@link Sequence} #2
	 */
	public Sequence getOriginalSequence2() {
		return originalSequence2;
	}

	/**
	 * 
	 * @param originalSequence2
	 */
	public void setOriginalSequence2(Sequence originalSequence2) {
		this.originalSequence2 = originalSequence2;
	}

	/**
	 * Returns the number of gaps of the aligned sequence #1
	 * 
	 * @return the number of gaps of the aligned sequence #1
	 */
	public int getGaps1() {
		int count = 0;
		for (int i = 0, n = sequence1.length; i < n; i++) {
			if (sequence1[i] == Alignment.GAP) {
				count++;
			}
		}
		return count;
	}

	/**
	 * Returns the number of gaps of the aligned sequence #2
	 * 
	 * @return the number of gaps of the aligned sequence #2
	 */
	public int getGaps2() {
		int count = 0;
		for (int i = 0, n = sequence2.length; i < n; i++) {
			if (sequence2[i] == Alignment.GAP) {
				count++;
			}
		}
		return count;
	}

	/**
	 * Adds this alignment to another alignment, the order is important.
	 * 
	 * @param a1
	 *            The 1st alignment
	 * @param a2
	 *            The 2nd alignment
	 * @return the sum of two alignment
	 */
	public static Alignment add(Alignment a1, Alignment a2) {

		if (a1 == null) {
			if (a2 == null) {
				return null;
			} else {
				return copy(a2);
			}
		} else {
			if (a2 == null) {
				return copy(a1);
			} else {
				Alignment sum = new Alignment();

				StringBuffer buffer = new StringBuffer();

				buffer.append(a1.getOriginalSequence1());
				buffer.append(a2.getOriginalSequence1());

				sum.setOriginalSequence1(new Sequence(buffer.toString()));

				buffer = new StringBuffer();

				buffer.append(a1.getOriginalSequence2());
				buffer.append(a2.getOriginalSequence2());

				sum.setOriginalSequence2(new Sequence(buffer.toString()));

				buffer = new StringBuffer();

				buffer.append(a1.getSequence1());
				buffer.append(a2.getSequence1());

				sum.setSequence1(buffer.toString().toCharArray());

				buffer = new StringBuffer();

				buffer.append(a1.getSequence2());
				buffer.append(a2.getSequence2());

				sum.setSequence2(buffer.toString().toCharArray());

				buffer = new StringBuffer();

				buffer.append(a1.getMarkupLine());
				buffer.append(a2.getMarkupLine());

				sum.setMarkupLine(buffer.toString().toCharArray());

				sum.setScore(a1.getScore() + a2.getScore());
				sum.setGaps(a1.getGaps() + a2.getGaps());

				sum.setStart1(a1.getStart1());
				sum.setStart2(a1.getStart2());
				sum.setExtend(a1.getExtend());
				sum.setOpen(a1.getOpen());
				sum.setMatrix(a1.getMatrix());

				return sum;
			}
		}
	}

	/**
	 * Copies an alignment to another alignment.
	 * 
	 * @param alignment
	 *            Alignment
	 */
	public static Alignment copy(Alignment alignment) {

		Alignment copy = new Alignment();

		copy.setSequence1(alignment.getSequence1());
		copy.setSequence2(alignment.getSequence2());
		copy.setMarkupLine(alignment.getMarkupLine());
		copy.setScore(alignment.getScore());
		copy.setGaps(alignment.getGaps());
		copy.setStart1(alignment.getStart1());
		copy.setStart2(alignment.getStart2());
		copy.setExtend(alignment.getExtend());
		copy.setOpen(alignment.getOpen());
		copy.setMatrix(alignment.getMatrix());
		copy.setOriginalSequence1(alignment.getOriginalSequence1());
		copy.setOriginalSequence2(alignment.getOriginalSequence2());

		return copy;
	}
}