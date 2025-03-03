import java.util.ArrayList;
import java.util.List;

public class SequenceAlignment {

    private static final String SEQ1 = "CCATGGTACATTTGGCTAGGTTTTATAGCTGGCTTGATTGCCATAGTAATGGTGACAATTATGCTTTGCTGTATGACCAGTTGCTGTAGCGC";
    private static final String SEQ2 = "CCGTATGTTTGTTAGCAAAACAAGTATCTGTAGATGCTATGTCACGAGTGACACCACCATCAATAGCCTTGTATCCTATGATTTCACTTGACC";

    public static void main(String[] args) {
        // Needleman-Wunsch (Global Alignment)
        System.out.println("\nNeedleman-Wunsch (Global Alignment):");
        AlignmentResult nwResult = needlemanWunsch(SEQ1, SEQ2);
        System.out.println("Score: " + nwResult.score);
        displayAlignment(nwResult.alignment1, nwResult.alignment2);
        System.out.println();

        // Smith-Waterman (Local Alignment)
        System.out.println("Smith-Waterman (Local Alignment):");
        AlignmentResult swResult = smithWaterman(SEQ1, SEQ2);
        System.out.println("Score: " + swResult.score);
        displayAlignment(swResult.alignment1, swResult.alignment2);
        System.out.println();

        // BLAST (Pairwise Heuristic Local Alignment)
        System.out.println("BLAST (Pairwise Heuristic Local Alignment):");
        blastAlignment(SEQ1, SEQ2);
    }

    // Needleman-Wunsch algorithm for Global Alignment
    public static AlignmentResult needlemanWunsch(String seq1, String seq2) {
        int len1 = seq1.length();
        int len2 = seq2.length();
        int[][] scoreMatrix = new int[len1 + 1][len2 + 1];
        int match = 2;
        int mismatch = -1;
        int gap = -2;

        // Initialize score matrix
        for (int i = 0; i <= len1; i++) {
            scoreMatrix[i][0] = i * gap;
        }
        for (int j = 0; j <= len2; j++) {
            scoreMatrix[0][j] = j * gap;
        }

        // Fill the score matrix
        for (int i = 1; i <= len1; i++) {
            for (int j = 1; j <= len2; j++) {
                int diagScore = scoreMatrix[i - 1][j - 1] + (seq1.charAt(i - 1) == seq2.charAt(j - 1) ? match : mismatch);
                int leftScore = scoreMatrix[i][j - 1] + gap;
                int upScore = scoreMatrix[i - 1][j] + gap;
                scoreMatrix[i][j] = Math.max(Math.max(diagScore, leftScore), upScore);
            }
        }

        // Traceback to find alignment
        StringBuilder alignedSeq1 = new StringBuilder();
        StringBuilder alignedSeq2 = new StringBuilder();
        int i = len1;
        int j = len2;
        while (i > 0 || j > 0) {
            if (i > 0 && j > 0 && scoreMatrix[i][j] == scoreMatrix[i - 1][j - 1] + (seq1.charAt(i - 1) == seq2.charAt(j - 1) ? match : mismatch)) {
                alignedSeq1.insert(0, seq1.charAt(i - 1));
                alignedSeq2.insert(0, seq2.charAt(j - 1));
                i--;
                j--;
            } else if (i > 0 && scoreMatrix[i][j] == scoreMatrix[i - 1][j] + gap) {
                alignedSeq1.insert(0, seq1.charAt(i - 1));
                alignedSeq2.insert(0, '-');
                i--;
            } else {
                alignedSeq1.insert(0, '-');
                alignedSeq2.insert(0, seq2.charAt(j - 1));
                j--;
            }
        }

        return new AlignmentResult(scoreMatrix[len1][len2], alignedSeq1.toString(), alignedSeq2.toString());
    }

    // Smith-Waterman algorithm for Local Alignment
    public static AlignmentResult smithWaterman(String seq1, String seq2) {
        int len1 = seq1.length();
        int len2 = seq2.length();
        int[][] scoreMatrix = new int[len1 + 1][len2 + 1];
        int match = 2;
        int mismatch = -1;
        int gap = -2;

        // Initialize score matrix
        for (int i = 0; i <= len1; i++) {
            scoreMatrix[i][0] = 0;
        }
        for (int j = 0; j <= len2; j++) {
            scoreMatrix[0][j] = 0;
        }

        // Fill the score matrix
        int maxScore = 0;
        int maxI = 0;
        int maxJ = 0;
        for (int i = 1; i <= len1; i++) {
            for (int j = 1; j <= len2; j++) {
                int diagScore = scoreMatrix[i - 1][j - 1] + (seq1.charAt(i - 1) == seq2.charAt(j - 1) ? match : mismatch);
                int leftScore = scoreMatrix[i][j - 1] + gap;
                int upScore = scoreMatrix[i - 1][j] + gap;
                scoreMatrix[i][j] = Math.max(0, Math.max(Math.max(diagScore, leftScore), upScore));
                if (scoreMatrix[i][j] > maxScore) {
                    maxScore = scoreMatrix[i][j];
                    maxI = i;
                    maxJ = j;
                }
            }
        }

        // Traceback to find alignment
        StringBuilder alignedSeq1 = new StringBuilder();
        StringBuilder alignedSeq2 = new StringBuilder();
        int i = maxI;
        int j = maxJ;
        while (i > 0 && j > 0 && scoreMatrix[i][j] != 0) {
            if (scoreMatrix[i][j] == scoreMatrix[i - 1][j - 1] + (seq1.charAt(i - 1) == seq2.charAt(j - 1) ? match : mismatch)) {
                alignedSeq1.insert(0, seq1.charAt(i - 1));
                alignedSeq2.insert(0, seq2.charAt(j - 1));
                i--;
                j--;
            } else if (scoreMatrix[i][j] == scoreMatrix[i - 1][j] + gap) {
                alignedSeq1.insert(0, seq1.charAt(i - 1));
                alignedSeq2.insert(0, '-');
                i--;
            } else {
                alignedSeq1.insert(0, '-');
                alignedSeq2.insert(0, seq2.charAt(j - 1));
                j--;
            }
        }

        return new AlignmentResult(maxScore, alignedSeq1.toString(), alignedSeq2.toString());
    }

    // BLAST (Pairwise Heuristic Local Alignment)
    public static void blastAlignment(String seq1, String seq2) {
        // Perform a simplified BLAST-like heuristic local alignment
        int matchScore = 2;
        int mismatchScore = -1;
        int gapPenalty = -2;
        int threshold = 15; // Set a score threshold for reporting alignments

        List<AlignmentResult> alignments = new ArrayList<>();

        // Perform local alignments using sliding window approach
        for (int start = 0; start <= seq2.length() - 20; start += 10) {
            String subSeq2 = seq2.substring(start, start + 20);
            AlignmentResult result = smithWaterman(seq1, subSeq2);
            if (result.score >= threshold) {
                alignments.add(result);
            }
        }

        // Display all significant alignments found
        for (AlignmentResult alignment : alignments) {
            System.out.println("Alignment:");
            System.out.println("Score: " + alignment.score);
            displayAlignment(alignment.alignment1, alignment.alignment2);
            System.out.println();
        }
    }

    // Display alignment
    public static void displayAlignment(String align1, String align2) {
        int len = align1.length();
        for (int i = 0; i < len; i += 75) {
            System.out.println(align1.substring(i, Math.min(len, i + 75)));
            System.out.println(getAlignmentString(align1, align2, i, Math.min(len, i + 75)));
            System.out.println(align2.substring(i, Math.min(len, i + 75)));
            System.out.println();
        }
    }

    // Helper method to get alignment line with symbols
    private static String getAlignmentString(String align1, String align2, int start, int end) {
        StringBuilder alignmentString = new StringBuilder();
        for (int i = start; i < end; i++) {
            if (align1.charAt(i) == align2.charAt(i)) {
                alignmentString.append('|');
            } else if (align1.charAt(i) == '-' || align2.charAt(i) == '-') {
                alignmentString.append(' ');
            } else {
                alignmentString.append(':');
            }
        }
        return alignmentString.toString();
    }

    // Alignment result class
    private static class AlignmentResult {
        int score;
        String alignment1;
        String alignment2;

        AlignmentResult(int score, String alignment1, String alignment2) {
            this.score = score;
            this.alignment1 = alignment1;
            this.alignment2 = alignment2;
        }
    }
}
