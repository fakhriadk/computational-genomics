import java.io.*;
import java.util.*;

public class CodonMap {

    private static Map<String, String> codonMap = new HashMap<>();

    static {
        codonMap.put("AUG", "M");
        codonMap.put("UUU", "F");
        codonMap.put("UUC", "F");
        codonMap.put("UUA", "L");
        codonMap.put("UUG", "L");
        codonMap.put("CUU", "L");
        codonMap.put("CUC", "L");
        codonMap.put("CUA", "L");
        codonMap.put("CUG", "L");
        codonMap.put("AUU", "I");
        codonMap.put("AUC", "I");
        codonMap.put("AUA", "I");
        codonMap.put("GUU", "V");
        codonMap.put("GUC", "V");
        codonMap.put("GUA", "V");
        codonMap.put("GUG", "V");
        codonMap.put("UCU", "S");
        codonMap.put("UCC", "S");
        codonMap.put("UCA", "S");
        codonMap.put("UCG", "S");
        codonMap.put("AGU", "S");
        codonMap.put("AGC", "S");
        codonMap.put("CCU", "P");
        codonMap.put("CCC", "P");
        codonMap.put("CCA", "P");
        codonMap.put("CCG", "P");
        codonMap.put("ACU", "T");
        codonMap.put("ACC", "T");
        codonMap.put("ACA", "T");
        codonMap.put("ACG", "T");
        codonMap.put("GCU", "A");
        codonMap.put("GCC", "A");
        codonMap.put("GCA", "A");
        codonMap.put("GCG", "A");
        codonMap.put("UAU", "Y");
        codonMap.put("UAC", "Y");
        codonMap.put("CAU", "H");
        codonMap.put("CAC", "H");
        codonMap.put("CAA", "Q");
        codonMap.put("CAG", "Q");
        codonMap.put("AAU", "N");
        codonMap.put("AAC", "N");
        codonMap.put("AAA", "K");
        codonMap.put("AAG", "K");
        codonMap.put("GAU", "D");
        codonMap.put("GAC", "D");
        codonMap.put("GAA", "E");
        codonMap.put("GAG", "E");
        codonMap.put("UGU", "C");
        codonMap.put("UGC", "C");
        codonMap.put("UGG", "W");
        codonMap.put("CGU", "R");
        codonMap.put("CGC", "R");
        codonMap.put("CGA", "R");
        codonMap.put("CGG", "R");
        codonMap.put("AGA", "R");
        codonMap.put("AGG", "R");
        codonMap.put("GGU", "G");
        codonMap.put("GGC", "G");
        codonMap.put("GGA", "G");
        codonMap.put("GGG", "G");
    }

    public static void main(String[] args) {
        String path = "/Users/macbook/IdeaProjects/Komputasi Genomik/src/FastaFiles/examplefile.fasta";
        Scanner scanner = new Scanner(System.in);
        System.out.print("Enter the protein sequence to encode: ");
        String proteinSequence = scanner.nextLine();
        Map<String, String> sequences = readFasta(path);

        for (String header : sequences.keySet()) {
            String sequence = sequences.get(header);
            System.out.println("\nIndex: >" + header);
            System.out.println("\nGenome Data: " + sequence);

            List<String> rnaSequences = findSequences(proteinSequence, "RNA");
            System.out.println("\nPossible RNA sequences encoding the protein sequence '" + proteinSequence + "':");
            for (String rnaSeq : rnaSequences) {
                System.out.println(rnaSeq);
            }

            List<String> dnaSequences = findSequences(proteinSequence, "DNA");
            System.out.println("\nPossible DNA sequences encoding the protein sequence '" + proteinSequence + "':");
            for (String dnaSeq : dnaSequences) {
                System.out.println(dnaSeq);
            }
        }
        scanner.close();
    }

    public static Map<String, String> readFasta(String filePath) {
        Map<String, String> sequences = new HashMap<>();
        try (BufferedReader br = new BufferedReader(new FileReader(filePath))) {
            String line;
            StringBuilder sequence = new StringBuilder();
            String header = null;
            while ((line = br.readLine()) != null) {
                line = line.trim();
                if (line.startsWith(">")) {
                    if (header != null) {
                        sequences.put(header, sequence.toString());
                        sequence = new StringBuilder();
                    }
                    header = line.substring(1);
                } else {
                    sequence.append(line);
                }
            }
            if (header != null) {
                sequences.put(header, sequence.toString());
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        return sequences;
    }

    public static List<String> findSequences(String protein, String nucleicAcid) {
        Map<String, List<String>> aminoToCodons = new HashMap<>();
        for (String codon : codonMap.keySet()) {
            String amino = codonMap.get(codon);
            aminoToCodons.computeIfAbsent(amino, k -> new ArrayList<>()).add(codon);
        }

        List<String> possibleSequences = new ArrayList<>();
        possibleSequences.add("");

        for (int i = 0; i < protein.length(); i++) {
            String amino = String.valueOf(protein.charAt(i));
            if (aminoToCodons.containsKey(amino)) {
                List<String> aminoCodons = aminoToCodons.get(amino);
                List<String> newSequences = new ArrayList<>();
                for (String seq : possibleSequences) {
                    for (String codon : aminoCodons) {
                        newSequences.add(seq + codon);
                    }
                }
                possibleSequences = newSequences;
            }
        }

        if (nucleicAcid.equals("DNA")){
            for (int i = 0; i < possibleSequences.size(); i++) {
                possibleSequences.set(i, possibleSequences.get(i).replace('U', 'T'));
            }
        }

        return possibleSequences;
    }
}
