import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;
import java.util.List;

public class FinalProject extends JFrame {

    private JTextArea outputArea;
    private File fastaFile;

    public FinalProject() {
        setTitle("INTEGRATED SEQUENCE ANALYSIS PROGRAM");
        setSize(800, 600);
        setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        setLocationRelativeTo(null);

        JButton uploadButton = new JButton("Upload FASTA File");
        uploadButton.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                uploadFastaFile();
            }
        });

        JButton sequenceStatsButton = new JButton("Calculate Sequence Statistics");
        sequenceStatsButton.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                if (fastaFile != null) {
                    calculateSequenceStatistics();
                } else {
                    JOptionPane.showMessageDialog(null, "Please upload a FASTA file first.");
                }
            }
        });

        JButton proteinMappingButton = new JButton("Protein Mapping");
        proteinMappingButton.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                if (fastaFile != null) {
                    performProteinMapping();
                } else {
                    JOptionPane.showMessageDialog(null, "Please upload a FASTA file first.");
                }
            }
        });

        JButton orfFinderButton = new JButton("ORF Finder");
        orfFinderButton.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                if (fastaFile != null) {
                    try {
                        orfFinder(fastaFile.getPath());
                    } catch (IOException ex) {
                        ex.printStackTrace();
                        JOptionPane.showMessageDialog(null, "Error reading FASTA file", "Error", JOptionPane.ERROR_MESSAGE);
                    }
                } else {
                    JOptionPane.showMessageDialog(null, "Please upload a FASTA file first.");
                }
            }
        });

        // New button for GC Content Analysis
        JButton gcContentButton = new JButton("GC Content Analysis");
        gcContentButton.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                if (fastaFile != null) {
                    calculateGCContent();
                } else {
                    JOptionPane.showMessageDialog(null, "Please upload a FASTA file first.");
                }
            }
        });

        outputArea = new JTextArea();
        outputArea.setEditable(false);

        JPanel buttonPanel = new JPanel();
        buttonPanel.setLayout(new GridLayout(1, 5)); // Ensure the layout can accommodate all buttons
        buttonPanel.add(uploadButton);
        buttonPanel.add(sequenceStatsButton);
        buttonPanel.add(proteinMappingButton);
        buttonPanel.add(orfFinderButton);
        buttonPanel.add(gcContentButton); // Add the new button to the panel

        JPanel panel = new JPanel();
        panel.setLayout(new BorderLayout());
        panel.add(buttonPanel, BorderLayout.NORTH);
        panel.add(new JScrollPane(outputArea), BorderLayout.CENTER);

        add(panel);
    }

    private void uploadFastaFile() {
        JFileChooser fileChooser = new JFileChooser();
        int returnValue = fileChooser.showOpenDialog(this);
        if (returnValue == JFileChooser.APPROVE_OPTION) {
            fastaFile = fileChooser.getSelectedFile();
            outputArea.append("FASTA file uploaded: " + fastaFile.getName() + "\n");
        }
    }

    private void calculateSequenceStatistics() {
        try {
            String[] fastaContent = readFasta(fastaFile.getPath());
            String annotation = fastaContent[0];
            String genomeData = fastaContent[1];

            outputArea.append("\nIndex:\n" + annotation);
            outputArea.append("\nGenome Data:\n" + genomeData);

            Map<String, double[]> baseStatistics = calculateBaseStatistics(genomeData);
            outputArea.append("\n\nBase Statistics:\n");
            for (Map.Entry<String, double[]> entry : baseStatistics.entrySet()) {
                outputArea.append(String.format("%s: %.2f%% (%.2f)\n", entry.getKey(), entry.getValue()[0], entry.getValue()[1]));
            }

            int windowSize = 3;  // Adjust window size as needed
            Map<String, Map<String, Double>> windowStatistics = calculateWindowStatistics(genomeData, windowSize);
            outputArea.append("\nFrequency of each base nitrogen:\n");
            for (Map.Entry<String, Map<String, Double>> entry : windowStatistics.entrySet()) {
                String index = entry.getKey();
                Map<String, Double> frequencies = entry.getValue();
                outputArea.append(String.format("%s: A: %.2f%%, C: %.2f%%, G: %.2f%%, T: %.2f%%\n", index, frequencies.get("A"), frequencies.get("C"), frequencies.get("G"), frequencies.get("T")));
            }
        } catch (IOException e) {
            e.printStackTrace();
            JOptionPane.showMessageDialog(this, "Error reading FASTA file", "Error", JOptionPane.ERROR_MESSAGE);
        }
    }

    private void performProteinMapping() {
        try {
            String proteinSequence = JOptionPane.showInputDialog(this, "Enter the protein sequence to encode:");
            if (proteinSequence != null && !proteinSequence.trim().isEmpty()) {
                Map<String, String> sequences = readFastaMap(fastaFile.getPath());

                for (String header : sequences.keySet()) {
                    String sequence = sequences.get(header);

                    java.util.List<String> rnaSequences = findSequences(proteinSequence, "RNA");
                    outputArea.append("\n\nPossible RNA sequences encoding the protein sequence '" + proteinSequence + "':\n");
                    for (String rnaSeq : rnaSequences) {
                        outputArea.append(rnaSeq + "\n");
                    }

                    java.util.List<String> dnaSequences = findSequences(proteinSequence, "DNA");
                    outputArea.append("\n\nPossible DNA sequences encoding the protein sequence '" + proteinSequence + "':\n");
                    for (String dnaSeq : dnaSequences) {
                        outputArea.append(dnaSeq + "\n");
                    }
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
            JOptionPane.showMessageDialog(this, "Error reading FASTA file", "Error", JOptionPane.ERROR_MESSAGE);
        }
    }

    private void orfFinder(String fastaFile) throws IOException {
        BufferedReader reader = new BufferedReader(new FileReader(fastaFile));
        String line;
        java.util.List<String> headers = new ArrayList<>();
        java.util.List<String> sequences = new ArrayList<>();

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

            outputArea.append("\nSequence file: " + header + "\n" + sequence);
            printCompleteCDSAndInvComplement(sequence);

            int seqLength = sequence.length();
            outputArea.append("\nSequence length: " + seqLength + "\n");

            // Sequence statistics for the entire sequence
            sequenceStatistics(sequence);

            // Finding ORFs for the entire sequence
            java.util.List<ORF> orfs = findORFs(sequence);

            outputArea.append("\nORFs Found: " + orfs.size() + "\n");
            printORFs(orfs);
        }
    }

    private String[] readFasta(String fastaFile) throws IOException {
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
        }
        return new String[]{annotation.toString(), genomeData.toString()};
    }

    private Map<String, String> readFastaMap(String filePath) throws IOException {
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
        }
        return sequences;
    }

    private Map<String, double[]> calculateBaseStatistics(String sequence) {
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
            double[] values = entry.getValue();
            values[0] = (values[1] / totalBases) * 100;
        }

        return baseStatistics;
    }

    private Map<String, Map<String, Double>> calculateWindowStatistics(String sequence, int windowSize) {
        Map<String, Map<String, Double>> windowStatistics = new HashMap<>();
        int sequenceLength = sequence.length();

        for (int i = 0; i <= sequenceLength - windowSize; i++) {
            String window = sequence.substring(i, i + windowSize);
            Map<String, Double> baseCounts = new HashMap<>();
            baseCounts.put("A", 0.0);
            baseCounts.put("C", 0.0);
            baseCounts.put("G", 0.0);
            baseCounts.put("T", 0.0);

            for (char base : window.toCharArray()) {
                baseCounts.put(String.valueOf(base), baseCounts.get(String.valueOf(base)) + 1);
            }

            for (Map.Entry<String, Double> entry : baseCounts.entrySet()) {
                baseCounts.put(entry.getKey(), (entry.getValue() / windowSize) * 100);
            }

            windowStatistics.put("Window " + (i + 1), baseCounts);
        }

        return windowStatistics;
    }

    private java.util.List<String> findSequences(String proteinSequence, String type) {
        List<String> sequences = new ArrayList<>();
        for (char aminoAcid : proteinSequence.toCharArray()) {
            List<String> codons = getCodons(aminoAcid, type);
            if (sequences.isEmpty()) {
                sequences.addAll(codons);
            } else {
                List<String> newSequences = new ArrayList<>();
                for (String sequence : sequences) {
                    for (String codon : codons) {
                        newSequences.add(sequence + codon);
                    }
                }
                sequences = newSequences;
            }
        }
        return sequences;
    }

    private List<String> getCodons(char aminoAcid, String type) {
        Map<Character, String[]> rnaCodonTable = new HashMap<>();
        rnaCodonTable.put('A', new String[]{"GCU", "GCC", "GCA", "GCG"});
        rnaCodonTable.put('R', new String[]{"CGU", "CGC", "CGA", "CGG", "AGA", "AGG"});
        rnaCodonTable.put('N', new String[]{"AAU", "AAC"});
        rnaCodonTable.put('D', new String[]{"GAU", "GAC"});
        rnaCodonTable.put('C', new String[]{"UGU", "UGC"});
        rnaCodonTable.put('Q', new String[]{"CAA", "CAG"});
        rnaCodonTable.put('E', new String[]{"GAA", "GAG"});
        rnaCodonTable.put('G', new String[]{"GGU", "GGC", "GGA", "GGG"});
        rnaCodonTable.put('H', new String[]{"CAU", "CAC"});
        rnaCodonTable.put('I', new String[]{"AUU", "AUC", "AUA"});
        rnaCodonTable.put('L', new String[]{"UUA", "UUG", "CUU", "CUC", "CUA", "CUG"});
        rnaCodonTable.put('K', new String[]{"AAA", "AAG"});
        rnaCodonTable.put('M', new String[]{"AUG"});
        rnaCodonTable.put('F', new String[]{"UUU", "UUC"});
        rnaCodonTable.put('P', new String[]{"CCU", "CCC", "CCA", "CCG"});
        rnaCodonTable.put('S', new String[]{"UCU", "UCC", "UCA", "UCG", "AGU", "AGC"});
        rnaCodonTable.put('T', new String[]{"ACU", "ACC", "ACA", "ACG"});
        rnaCodonTable.put('W', new String[]{"UGG"});
        rnaCodonTable.put('Y', new String[]{"UAU", "UAC"});
        rnaCodonTable.put('V', new String[]{"GUU", "GUC", "GUA", "GUG"});
        rnaCodonTable.put('*', new String[]{"UAA", "UAG", "UGA"});

        String[] rnaCodons = rnaCodonTable.get(aminoAcid);
        if (rnaCodons == null) {
            throw new IllegalArgumentException("Invalid amino acid: " + aminoAcid);
        }

        List<String> codons = new ArrayList<>(Arrays.asList(rnaCodons));
        if (type.equalsIgnoreCase("DNA")) {
            List<String> dnaCodons = new ArrayList<>();
            for (String rnaCodon : codons) {
                dnaCodons.add(rnaCodon.replace('U', 'T'));
            }
            return dnaCodons;
        }
        return codons;
    }

    private void printCompleteCDSAndInvComplement(String sequence) {
        // Complete CDS
        String completeCDS = sequence;
        outputArea.append("\n\nComplete CDS: " + completeCDS + "\n");

        // Inverted complementary chain
        String invertedComplement = new StringBuilder(sequence).reverse().toString();
        invertedComplement = invertedComplement.replace('A', 't').replace('T', 'a').replace('C', 'g').replace('G', 'c').toUpperCase();
        outputArea.append("\nInverted Complementary Chain: " + invertedComplement + "\n");
    }

    private void sequenceStatistics(String sequence) {
        // Frequency of each nucleotide in the sequence
        Map<Character, Integer> frequencyMap = new HashMap<>();
        for (char nucleotide : sequence.toCharArray()) {
            frequencyMap.put(nucleotide, frequencyMap.getOrDefault(nucleotide, 0) + 1);
        }
        outputArea.append("\nNucleotide Frequency:\n");
        for (Map.Entry<Character, Integer> entry : frequencyMap.entrySet()) {
            outputArea.append(entry.getKey() + ": " + entry.getValue() + "\n");
        }
    }

    private java.util.List<ORF> findORFs(String sequence) {
        List<ORF> orfs = new ArrayList<>();
        String[] startCodons = {"ATG"};
        String[] stopCodons = {"TAA", "TAG", "TGA"};

        for (int frame = 0; frame < 3; frame++) {
            for (int i = frame; i < sequence.length() - 2; i += 3) {
                String codon = sequence.substring(i, i + 3);
                if (Arrays.asList(startCodons).contains(codon)) {
                    for (int j = i + 3; j < sequence.length() - 2; j += 3) {
                        String stopCodon = sequence.substring(j, j + 3);
                        if (Arrays.asList(stopCodons).contains(stopCodon)) {
                            orfs.add(new ORF(i, j + 3, sequence.substring(i, j + 3)));
                            break;
                        }
                    }
                }
            }
        }

        return orfs;
    }

    private void printORFs(java.util.List<ORF> orfs) {
        for (ORF orf : orfs) {
            outputArea.append("ORF: Start = " + orf.start + ", End = " + orf.end + ", Sequence = " + orf.sequence + "\n");
        }
    }

    private static class ORF {
        int start;
        int end;
        String sequence;

        ORF(int start, int end, String sequence) {
            this.start = start;
            this.end = end;
            this.sequence = sequence;
        }
    }

    private void calculateGCContent() {
        try {
            String[] fastaContent = readFasta(fastaFile.getPath());
            String genomeData = fastaContent[1];

            int gcCount = 0;
            for (char base : genomeData.toCharArray()) {
                if (base == 'G' || base == 'C') {
                    gcCount++;
                }
            }

            double gcContent = ((double) gcCount / genomeData.length()) * 100;
            outputArea.append("\nGC Content: " + String.format("%.2f", gcContent) + "%\n");
        } catch (IOException e) {
            e.printStackTrace();
            JOptionPane.showMessageDialog(this, "Error reading FASTA file", "Error", JOptionPane.ERROR_MESSAGE);
        }
    }

    public static void main(String[] args) {
        SwingUtilities.invokeLater(new Runnable() {
            @Override
            public void run() {
                new FinalProject().setVisible(true);
            }
        });
    }
}