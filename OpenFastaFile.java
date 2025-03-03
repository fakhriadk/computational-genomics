import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;

public class OpenFastaFile {

    public static void main(String[] args) {
        String path = "/Users/macbook/IdeaProjects/Komputasi Genomik/src/FastaFiles/examplefile.fasta";
        readAndAnalyzeFasta(path);
    }

    public static void readAndAnalyzeFasta(String fastaFile) {
        HashMap<String, String> sequences = readFasta(fastaFile);

        for (String sequenceId : sequences.keySet()) {
            String sequence = sequences.get(sequenceId);
            System.out.println("\nIndex: " + sequenceId);
            System.out.println("\nGenome Data: " + sequence);
            calculateNucleotideFrequency(sequence);
        }
    }

    public static HashMap<String, String> readFasta(String fastaFile) {
        HashMap<String, String> sequences = new HashMap<>();
        try (BufferedReader br = new BufferedReader(new FileReader(fastaFile))) {
            String line;
            String sequenceId = "";
            StringBuilder sequence = new StringBuilder();

            while ((line = br.readLine()) != null) {
                if (line.startsWith(">")) {
                    if (!sequenceId.isEmpty()) {
                        sequences.put(sequenceId, sequence.toString());
                    }
                    sequenceId = line.trim();
                    sequence = new StringBuilder();
                } else {
                    sequence.append(line.trim());
                }
            }
            if (!sequenceId.isEmpty()) {
                sequences.put(sequenceId, sequence.toString());
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        return sequences;
    }

    public static void calculateNucleotideFrequency(String sequence) {
        char[] nucleotides = {'A', 'C', 'G', 'T'};
        int[] nucleotideCounts = new int[4];

        for (char nucleotide : sequence.toCharArray()) {
            switch (nucleotide) {
                case 'A':
                    nucleotideCounts[0]++;
                    break;
                case 'C':
                    nucleotideCounts[1]++;
                    break;
                case 'G':
                    nucleotideCounts[2]++;
                    break;
                case 'T':
                    nucleotideCounts[3]++;
                    break;
            }
        }

        int totalNucleotides = sequence.length();
        if (totalNucleotides > 0) {
            System.out.println("\nFrequency of nucleotides:");
            for (int i = 0; i < nucleotides.length; i++) {
                double frequency = (double) nucleotideCounts[i] / totalNucleotides;
                System.out.printf("%c: %.4f%n", nucleotides[i], frequency);
            }

            int cgCount = nucleotideCounts[1] + nucleotideCounts[2]; // C + G
            int atCount = nucleotideCounts[0] + nucleotideCounts[3]; // A + T
            int totalBp = cgCount + atCount;

            if (totalBp > 0) {
                double cgRatio = (double) cgCount / totalBp;
                double atRatio = (double) atCount / totalBp;
                System.out.println("\nBase pair ratios:");
                System.out.printf("CG ratio: %.4f%n", cgRatio);
                System.out.printf("AT ratio: %.4f%n", atRatio);
            } else {
                System.out.println("No base pairs found.");
            }
        } else {
            System.out.println("No nucleotides found.");
        }
    }
}
