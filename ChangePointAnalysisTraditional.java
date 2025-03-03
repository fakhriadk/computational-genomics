import javax.swing.*;
import java.awt.*;
import java.awt.event.*;
import java.io.*;
import java.util.HashMap;
import java.util.Map;

public class ChangePointAnalysisTraditional extends JFrame {
    private static final String path = "/Users/macbook/IdeaProjects/Komputasi Genomik/src/FastaFiles/examplefile.fasta";

    private String recSeq = "";
    private Map<String, String> metadata = new HashMap<>();
    private int windowSize = 10;
    private double[][] baseFrequencies;
    private double[] atDensity;
    private double[] gcDensity;

    public ChangePointAnalysisTraditional() {
        setTitle("DNA Sequence Analysis");
        setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);

        JPanel mainPanel = new JPanel();
        mainPanel.setLayout(new BoxLayout(mainPanel, BoxLayout.Y_AXIS));

        JPanel controlPanel = new JPanel();
        controlPanel.add(new JLabel("Window Size:"));

        JSlider windowSizeSlider = new JSlider(JSlider.HORIZONTAL, 10, 500, windowSize);
        windowSizeSlider.setPreferredSize(new Dimension(150, 50));
        windowSizeSlider.addChangeListener(e -> {
            windowSize = windowSizeSlider.getValue();
            updatePlot();
        });
        controlPanel.add(windowSizeSlider);
        mainPanel.add(controlPanel);

        JPanel plotPanel = new JPanel();
        plotPanel.setLayout(new BoxLayout(plotPanel, BoxLayout.Y_AXIS));
        plotPanel.setPreferredSize(new Dimension(800, 800));
        mainPanel.add(plotPanel);

        JPanel nucleotidePanel = new JPanel() {
            @Override
            protected void paintComponent(Graphics g) {
                super.paintComponent(g);
                plotBaseFrequencies(g);
                drawNucleotideLegend(g);
                drawTitle(g, "Nucleotide density with sliding window", 20);
                drawAxes(g, 50, 50, 700, 350, "Position", "Density");
            }
        };
        nucleotidePanel.setPreferredSize(new Dimension(800, 400));
        plotPanel.add(nucleotidePanel);

        JPanel densityPanel = new JPanel() {
            @Override
            protected void paintComponent(Graphics g) {
                super.paintComponent(g);
                plotAtGcDensity(g);
                drawDensityLegend(g);
                drawTitle(g, "Base pair density with sliding window", 20);
                drawAxes(g, 50, 50, 500, 350, "Position", "Density");
            }
        };
        densityPanel.setPreferredSize(new Dimension(800, 400));
        plotPanel.add(densityPanel);

        loadFastaFile();

        getContentPane().add(mainPanel);
        pack();
        setLocationRelativeTo(null);
        setVisible(true);
    }

    private void loadFastaFile() {
        try (BufferedReader br = new BufferedReader(new FileReader(path))) {
            String line;
            StringBuilder seqBuilder = new StringBuilder();
            while ((line = br.readLine()) != null) {
                if (line.startsWith(">")) {
                    String[] headerParts = line.substring(1).split("\\|");
                    metadata.put("id", headerParts[0]);
                    metadata.put("name", headerParts.length > 1 ? headerParts[1] : "N/A");
                    metadata.put("organism", headerParts.length > 2 ? headerParts[2] : "N/A");
                } else {
                    seqBuilder.append(line.trim());
                }
            }
            recSeq = seqBuilder.toString();
            updatePlot();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private void updatePlot() {
        baseFrequencies = calculateBaseFrequencies(recSeq, windowSize);
        double[][] densities = calculateAtGcDensity(recSeq, windowSize);
        atDensity = densities[0];
        gcDensity = densities[1];

        // Repaint the panels
        repaint();
    }

    private double[][] calculateBaseFrequencies(String seq, int windowSize) {
        double[][] baseFreq = new double[4][seq.length() - windowSize + 1];

        for (int i = 0; i <= seq.length() - windowSize; i++) {
            String window = seq.substring(i, i + windowSize).toUpperCase();
            for (int j = 0; j < windowSize; j++) {
                char base = window.charAt(j);
                switch (base) {
                    case 'A':
                        baseFreq[0][i]++;
                        break;
                    case 'C':
                        baseFreq[1][i]++;
                        break;
                    case 'G':
                        baseFreq[2][i]++;
                        break;
                    case 'T':
                        baseFreq[3][i]++;
                        break;
                }
            }
        }

        // Normalize frequencies
        for (int row = 0; row < 4; row++) {
            for (int col = 0; col < seq.length() - windowSize + 1; col++) {
                baseFreq[row][col] /= windowSize;
            }
        }

        return baseFreq;
    }

    private double[][] calculateAtGcDensity(String seq, int windowSize) {
        double[] atDensity = new double[seq.length() - windowSize + 1];
        double[] gcDensity = new double[seq.length() - windowSize + 1];

        for (int i = 0; i <= seq.length() - windowSize; i++) {
            String window = seq.substring(i, i + windowSize).toUpperCase();
            int atCount = 0;
            int gcCount = 0;
            for (int j = 0; j < windowSize; j++) {
                char base = window.charAt(j);
                if (base == 'A' || base == 'T') {
                    atCount++;
                } else if (base == 'G' || base == 'C') {
                    gcCount++;
                }
            }
            atDensity[i] = (double) atCount / windowSize;
            gcDensity[i] = (double) gcCount / windowSize;
        }

        return new double[][]{atDensity, gcDensity};
    }

    private void plotBaseFrequencies(Graphics g) {
        if (recSeq.isEmpty()) return;

        int plotWidth = 800;
        int plotHeight = 400;
        int margin = 50;
        int plotAreaWidth = plotWidth - 2 * margin;
        int plotAreaHeight = plotHeight - 2 * margin;

        // Draw axes
        g.drawLine(margin, plotHeight - margin, plotWidth - margin, plotHeight - margin);
        g.drawLine(margin, plotHeight - margin, margin, margin);

        // Draw base frequencies
        drawGraph(g, baseFrequencies[0], plotAreaWidth, plotAreaHeight, margin, Color.BLUE); // A
        drawGraph(g, baseFrequencies[1], plotAreaWidth, plotAreaHeight, margin, Color.GREEN); // C
        drawGraph(g, baseFrequencies[2], plotAreaWidth, plotAreaHeight, margin, Color.RED); // G
        drawGraph(g, baseFrequencies[3], plotAreaWidth, plotAreaHeight, margin, Color.MAGENTA); // T
    }

    private void plotAtGcDensity(Graphics g) {
        if (recSeq.isEmpty()) return;

        int plotWidth = 800;
        int plotHeight = 400;
        int margin = 50;
        int plotAreaWidth = plotWidth - 2 * margin;
        int plotAreaHeight = plotHeight - 2 * margin;

        // Draw axes
        g.drawLine(margin, plotHeight - margin, plotWidth - margin, plotHeight - margin);
        g.drawLine(margin, plotHeight - margin, margin, margin);

        // Draw AT-GC density
        drawGraph(g, atDensity, plotAreaWidth, plotAreaHeight, margin, Color.ORANGE); // AT density
        drawGraph(g, gcDensity, plotAreaWidth, plotAreaHeight, margin, Color.YELLOW); // GC density
    }

    private void drawGraph(Graphics g, double[] data, int plotAreaWidth, int plotAreaHeight, int margin, Color color) {
        int numPoints = data.length;
        int dx = plotAreaWidth / numPoints;
        int x0 = margin;
        int y0 = plotAreaHeight + margin;

        // Find max value in data to normalize y-axis
        double maxValue = 0;
        for (double value : data) {
            if (value > maxValue) {
                maxValue = value;
            }
        }

        g.setColor(color);
        for (int i = 0; i < numPoints - 1; i++) {
            int x1 = x0 + i * dx;
            int y1 = y0 - (int) (data[i] * plotAreaHeight / maxValue);
            int x2 = x0 + (i + 1) * dx;
            int y2 = y0 - (int) (data[i + 1] * plotAreaHeight / maxValue);
            g.drawLine(x1, y1, x2, y2);
        }
    }

    private void drawNucleotideLegend(Graphics g) {
        int legendX = 700;
        int legendY = 50;
        int legendWidth = 80;
        int legendHeight = 90;
        g.setColor(Color.WHITE);
        g.fillRect(legendX, legendY, legendWidth, legendHeight);
        g.setColor(Color.BLACK);
        g.drawRect(legendX, legendY, legendWidth, legendHeight);

        g.setColor(Color.BLUE);
        g.fillRect(legendX + 10, legendY + 10, 10, 10);
        g.setColor(Color.BLACK);
        g.drawString("A", legendX + 30, legendY + 20);

        g.setColor(Color.GREEN);
        g.fillRect(legendX + 10, legendY + 30, 10, 10);
        g.setColor(Color.BLACK);
        g.drawString("C", legendX + 30, legendY + 40);

        g.setColor(Color.RED);
        g.fillRect(legendX + 10, legendY + 50, 10, 10);
        g.setColor(Color.BLACK);
        g.drawString("G", legendX + 30, legendY + 60);

        g.setColor(Color.MAGENTA);
        g.fillRect(legendX + 10, legendY + 70, 10, 10);
        g.setColor(Color.BLACK);
        g.drawString("T", legendX + 30, legendY + 80);
    }

    private void drawDensityLegend(Graphics g) {
        int legendX = 700;
        int legendY = 50;
        int legendWidth = 80;
        int legendHeight = 50;
        g.setColor(Color.WHITE);
        g.fillRect(legendX, legendY, legendWidth, legendHeight);
        g.setColor(Color.BLACK);
        g.drawRect(legendX, legendY, legendWidth, legendHeight);

        g.setColor(Color.ORANGE);
        g.fillRect(legendX + 10, legendY + 10, 10, 10);
        g.setColor(Color.BLACK);
        g.drawString("A-T", legendX + 30, legendY + 20);

        g.setColor(Color.YELLOW);
        g.fillRect(legendX + 10, legendY + 30, 10, 10);
        g.setColor(Color.BLACK);
        g.drawString("C-G", legendX + 30, legendY + 40);
    }

    private void drawTitle(Graphics g, String title, int yPosition) {
        g.setColor(Color.BLACK);
        g.setFont(new Font("Arial", Font.BOLD, 20));
        FontMetrics metrics = g.getFontMetrics();
        int titleWidth = metrics.stringWidth(title);
        int xPosition = (800 - titleWidth) / 2; // Centered horizontally
        g.drawString(title, xPosition, yPosition);
    }

    private void drawAxes(Graphics g, int x0, int y0, int width, int height, String xLabel, String yLabel) {
        g.setColor(Color.BLACK);
        g.setFont(new Font("Arial", Font.PLAIN, 12));

        // Draw x-axis numbers
        int numDivisions = 5;
        int dx = width / numDivisions;
        for (int i = 0; i <= numDivisions; i++) {
            int x = x0 + i * dx;
            g.drawString(String.valueOf(i * (baseFrequencies[0].length / numDivisions)), x - 10, y0 + height + 20);
        }

        // Draw y-axis numbers
        int dy = height / numDivisions;
        for (int i = 0; i <= numDivisions; i++) {
            int y = y0 + height - i * dy;
            g.drawString(String.format("%.1f", i * (1.0 / numDivisions)), x0 - 30, y + 5);
        }

        // Draw labels
        g.drawString(xLabel, x0 + width / 2 - 20, y0 + height + 40);
        Graphics2D g2d = (Graphics2D) g;
        g2d.rotate(-Math.PI / 2);
        g.drawString(yLabel, -(y0 + height / 2 + 20), x0 - 30);
        g2d.rotate(Math.PI / 2);
    }

    public static void main(String[] args) {
        SwingUtilities.invokeLater(new Runnable() {
            public void run() {
                new ChangePointAnalysisTraditional();
            }
        });
    }
}
