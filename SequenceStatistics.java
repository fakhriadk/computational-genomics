import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

public class SequenceStatistics {

    public static void main(String[] args) {
        String path = "/Users/macbook/IdeaProjects/Komputasi Genomik/src/FastaFiles/examplefile.fasta";
        String[] fastaContent = readFasta(path);
        String annotation = fastaContent[0];
        String genomeData = fastaContent[1];

        System.out.println("\nIndex:\n" + annotation);
        System.out.println("Genome Data:\n" + genomeData);

        Map<String, double[]> baseStatistics = calculateBaseStatistics(genomeData);
        System.out.println("\nBase Statistics:");
        for (Map.Entry<String, double[]> entry : baseStatistics.entrySet()) {
            System.out.printf("%s: %.2f%% (%.2f)\n", entry.getKey(), entry.getValue()[0], entry.getValue()[1]);
        }

        int windowSize = 3;  // Adjust window size as needed
        Map<String, Map<String, Double>> windowStatistics = calculateWindowStatistics(genomeData, windowSize);
        System.out.println("\nFrequency of each base nitrogen:");
        for (Map.Entry<String, Map<String, Double>> entry : windowStatistics.entrySet()) {
            String index = entry.getKey();
            Map<String, Double> frequencies = entry.getValue();
            System.out.printf("%s: A: %.2f%%, C: %.2f%%, G: %.2f%%, T: %.2f%%\n", index, frequencies.get("A"), frequencies.get("C"), frequencies.get("G"), frequencies.get("T"));
        }
    }

    private static String[] readFasta(String fastaFile) {
        StringBuilder annotation = new StringBuilder();
        StringBuilder genomeData = new StringBuilder();
        try (BufferedReader br = new BufferedReader(new FileReader(fastaFile))) {
            String line;
            while ((line = br.readLine()) != null) {
                if (line.startsWith(">")) {
                    annotation.append(line).append("\n");
                } else {
                    genomeData.append(line.trim());
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        return new String[]{annotation.toString(), genomeData.toString()};
    }

    private static Map<String, double[]> calculateBaseStatistics(String sequence) {
        Map<String, double[]> baseStatistics = new HashMap<>();
        baseStatistics.put("A", new double[]{0, 0});
        baseStatistics.put("T", new double[]{0, 0});
        baseStatistics.put("C", new double[]{0, 0});
        baseStatistics.put("G", new double[]{0, 0});

        int totalBases = sequence.length();
        for (int i = 0; i < totalBases; i++) {
            char base = sequence.charAt(i);
            if (baseStatistics.containsKey(String.valueOf(base))) {
                baseStatistics.get(String.valueOf(base))[1]++;
            }
        }

        for (Map.Entry<String, double[]> entry : baseStatistics.entrySet()) {
            double count = entry.getValue()[1];
            double percentage = count / totalBases * 100;
            entry.setValue(new double[]{percentage, count / totalBases});
        }
        return baseStatistics;
    }

    private static Map<String, Map<String, Double>> calculateWindowStatistics(String sequence, int windowSize) {
        Map<String, Map<String, Double>> windowStatistics = new HashMap<>();
        for (int i = 0; i <= sequence.length() - windowSize; i++) {
            String window = sequence.substring(i, i + windowSize);
            Map<String, Integer> nucleotideCounts = countBases(window);
            double totalNucleotides = nucleotideCounts.values().stream().mapToInt(Integer::intValue).sum();
            Map<String, Double> windowFrequency = new HashMap<>();
            windowFrequency.put("A", 0.0);
            windowFrequency.put("C", 0.0);
            windowFrequency.put("G", 0.0);
            windowFrequency.put("T", 0.0);

            for (Map.Entry<String, Double> entry : windowFrequency.entrySet()) {
                String nucleotide = entry.getKey();
                if (totalNucleotides > 0) {
                    double frequency = nucleotideCounts.get(nucleotide) / totalNucleotides * 100;
                    windowFrequency.put(nucleotide, frequency);
                } else {
                    windowFrequency.put(nucleotide, 0.0);
                }
            }
            windowStatistics.put("W" + (i + 1), windowFrequency);
        }
        return windowStatistics;
    }

    private static Map<String, Integer> countBases(String sequence) {
        Map<String, Integer> nucleotideCounts = new HashMap<>();
        nucleotideCounts.put("A", 0);
        nucleotideCounts.put("C", 0);
        nucleotideCounts.put("G", 0);
        nucleotideCounts.put("T", 0);

        for (int i = 0; i < sequence.length(); i++) {
            char base = sequence.charAt(i);
            if (nucleotideCounts.containsKey(String.valueOf(base))) {
                nucleotideCounts.put(String.valueOf(base), nucleotideCounts.get(String.valueOf(base)) + 1);
            }
        }
        return nucleotideCounts;
    }
}
