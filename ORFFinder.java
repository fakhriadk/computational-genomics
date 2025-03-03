import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class ORFFinder {

    public static void main(String[] args) {
        String fastaFile = "/Users/macbook/IdeaProjects/Komputasi Genomik/src/FastaFiles/examplefile.fasta";
        try {
            orfFinder(fastaFile);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static void orfFinder(String fastaFile) throws IOException {
        BufferedReader reader = new BufferedReader(new FileReader(fastaFile));
        String line;
        List<String> headers = new ArrayList<>();
        List<String> sequences = new ArrayList<>();

        while ((line = reader.readLine()) != null) {
            if (line.startsWith(">")) {
                headers.add(line.substring(1).trim());
                sequences.add("");
            } else {
                int lastIndex = sequences.size() - 1;
                sequences.set(lastIndex, sequences.get(lastIndex) + line.trim());
            }
        }
        reader.close();

        for (int k = 0; k < sequences.size(); k++) {
            String header = headers.get(k);
            String sequence = sequences.get(k);

            System.out.printf("\nSequence file: %s\n%s\n", header, sequence);
            printCompleteCDSAndInvComplement(sequence);

            int seqLength = sequence.length();
            System.out.printf("\nSequence length: %d\n", seqLength);

            // Sequence statistics for the entire sequence
            sequenceStatistics(sequence);

            // Finding ORFs for the entire sequence
            List<ORF> orfs = findORFs(sequence);

            System.out.printf("ORFs Found: %d\n", orfs.size());
            printORFs(orfs);
        }
    }

    public static List<ORF> findORFs(String sequence) {
        int minORFLength = 100;
        String startCodon = "ATG";
        String[] stopCodons = {"TAA", "TAG", "TGA"};
        List<ORF> orfs = new ArrayList<>();

        for (int strand : new int[]{1, -1}) {
            String seq = (strand == -1) ? reverseComplement(sequence) : sequence;
            for (int frame = 0; frame < 3; frame++) {
                int length = 3 * ((seq.length() - frame) / 3);
                for (int i = frame; i < length; i += 3) {
                    String codon = seq.substring(i, i + 3);
                    if (codon.equals(startCodon)) {
                        for (int j = i + 3; j < length; j += 3) {
                            String stopCodon = seq.substring(j, j + 3);
                            if (contains(stopCodons, stopCodon)) {
                                if (j + 3 - i >= minORFLength) {
                                    String orfSeq = seq.substring(i, j + 3);
                                    orfs.add(new ORF(strand, frame + 1, i, j + 2, orfSeq));
                                }
                                break;
                            }
                        }
                    }
                }
            }
        }
        return orfs;
    }

    public static void printORFs(List<ORF> orfs) {
        for (int index = 0; index < orfs.size(); index++) {
            ORF orf = orfs.get(index);
            System.out.printf("\n[Frame Type: %s%d][Start: %d Stop: %d] [Length: %d]\n",
                    (orf.strand == 1) ? "+" : "-", orf.frame, orf.start, orf.end, orf.end - orf.start + 1);
            System.out.printf("ORF[%d]: %s\n", index + 1, orf.sequence);
            System.out.println("Complete ORF Data");
            printCodons(orf.sequence);
        }
        System.out.println();
    }

    public static void printCodons(String sequence) {
        Map<String, String> geneticCode = getGeneticCode();
        for (int i = 0; i <= sequence.length() - 3; i += 3) {
            String codon = sequence.substring(i, i + 3);
            String aminoAcid = geneticCode.get(codon);
            String label = "";
            if (codon.equals("ATG")) {
                label = " (Start)************";
            } else if (codon.equals("TAA") || codon.equals("TAG") || codon.equals("TGA")) {
                label = " (Stop Codon) *****************";
            }
            System.out.printf("Codon[%d,%d]: %s %s%s\n", (i / 3) + 1, i / 3 + 1, codon, aminoAcid, label);
        }
        System.out.println();
    }

    public static void sequenceStatistics(String sequence) {
        Map<Character, Integer> baseCounts = new HashMap<>();
        int totalBases = sequence.length();

        for (char base : sequence.toCharArray()) {
            baseCounts.put(base, baseCounts.getOrDefault(base, 0) + 1);
        }

        System.out.println("\nSequence Statistics:");
        for (char base : "ATGC".toCharArray()) {
            System.out.printf("total %c: %d\n", base, baseCounts.getOrDefault(base, 0));
        }
        System.out.printf("\nTotal Base Pair: %d\n", totalBases);

        for (char base : "ATGC".toCharArray()) {
            double freq = (double) baseCounts.getOrDefault(base, 0) / totalBases;
            System.out.printf("%c frequency: %.4f\n", base, freq);
        }
        System.out.println();
    }

    public static void printCompleteCDSAndInvComplement(String sequence) {
        System.out.println("\nComplete CDS:");
        System.out.println(sequence);
        System.out.println("\nInverse Complement:");
        System.out.println(reverseComplement(sequence));
    }

    public static String reverseComplement(String sequence) {
        StringBuilder complement = new StringBuilder();
        for (int i = sequence.length() - 1; i >= 0; i--) {
            char base = sequence.charAt(i);
            switch (base) {
                case 'A':
                    complement.append('T');
                    break;
                case 'T':
                    complement.append('A');
                    break;
                case 'G':
                    complement.append('C');
                    break;
                case 'C':
                    complement.append('G');
                    break;
            }
        }
        return complement.toString();
    }

    public static boolean contains(String[] array, String target) {
        for (String s : array) {
            if (s.equals(target)) {
                return true;
            }
        }
        return false;
    }

    public static Map<String, String> getGeneticCode() {
        Map<String, String> geneticCode = new HashMap<>();
        geneticCode.put("TTT", "F");
        geneticCode.put("TTC", "F");
        geneticCode.put("TTA", "L");
        geneticCode.put("TTG", "L");
        geneticCode.put("CTT", "L");
        geneticCode.put("CTC", "L");
        geneticCode.put("CTA", "L");
        geneticCode.put("CTG", "L");
        geneticCode.put("ATT", "I");
        geneticCode.put("ATC", "I");
        geneticCode.put("ATA", "I");
        geneticCode.put("ATG", "M");
        geneticCode.put("GTT", "V");
        geneticCode.put("GTC", "V");
        geneticCode.put("GTA", "V");
        geneticCode.put("GTG", "V");
        geneticCode.put("TCT", "S");
        geneticCode.put("TCC", "S");
        geneticCode.put("TCA", "S");
        geneticCode.put("TCG", "S");
        geneticCode.put("CCT", "P");
        geneticCode.put("CCC", "P");
        geneticCode.put("CCA", "P");
        geneticCode.put("CCG", "P");
        geneticCode.put("ACT", "T");
        geneticCode.put("ACC", "T");
        geneticCode.put("ACA", "T");
        geneticCode.put("ACG", "T");
        geneticCode.put("GCT", "A");
        geneticCode.put("GCC", "A");
        geneticCode.put("GCA", "A");
        geneticCode.put("GCG", "A");
        geneticCode.put("TAT", "Y");
        geneticCode.put("TAC", "Y");
        geneticCode.put("TAA", "*");
        geneticCode.put("TAG", "*");
        geneticCode.put("CAT", "H");
        geneticCode.put("CAC", "H");
        geneticCode.put("CAA", "Q");
        geneticCode.put("CAG", "Q");
        geneticCode.put("AAT", "N");
        geneticCode.put("AAC", "N");
        geneticCode.put("AAA", "K");
        geneticCode.put("AAG", "K");
        geneticCode.put("GAT", "D");
        geneticCode.put("GAC", "D");
        geneticCode.put("GAA", "E");
        geneticCode.put("GAG", "E");
        geneticCode.put("TGT", "C");
        geneticCode.put("TGC", "C");
        geneticCode.put("TGA", "*");
        geneticCode.put("TGG", "W");
        geneticCode.put("CGT", "R");
        geneticCode.put("CGC", "R");
        geneticCode.put("CGA", "R");
        geneticCode.put("CGG", "R");
        geneticCode.put("AGT", "S");
        geneticCode.put("AGC", "S");
        geneticCode.put("AGA", "R");
        geneticCode.put("AGG", "R");
        geneticCode.put("GGT", "G");
        geneticCode.put("GGC", "G");
        geneticCode.put("GGA", "G");
        geneticCode.put("GGG", "G");
        return geneticCode;
    }

    public static class ORF {
        int strand;
        int frame;
        int start;
        int end;
        String sequence;

        public ORF(int strand, int frame, int start, int end, String sequence) {
            this.strand = strand;
            this.frame = frame;
            this.start = start;
            this.end = end;
            this.sequence = sequence;
        }
    }
}
