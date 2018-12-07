package uk.ac.ebi.pride.spectracluster.clusteringfileconverter.util;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.HashMap;
import java.util.Map;

/**
 * Created by jg on 21.05.15.
 */
public class FastaFile {
    private final File inputFile;
    private final Map<String, String> proteinSequences;

    public FastaFile(File inputFile) throws Exception {
        this.inputFile = inputFile;
        this.proteinSequences = loadFastaFile(inputFile);
    }

    protected Map<String, String> loadFastaFile(File file) throws Exception {
        FileReader fileReader = new FileReader(file);
        BufferedReader reader = new BufferedReader(fileReader);

        String line;
        String currentAnnotation = "UNKNOWN";
        StringBuilder currentSequence = new StringBuilder();
        Map<String, String> sequences = new HashMap<String, String>();

        while ((line = reader.readLine()) != null) {
            line = line.trim();

            // ignore empty lines
            if (line.length() < 1)
                continue;

            // ignore any comment line and store the previous protein sequence
            if (line.startsWith(">")) {
                if (currentSequence.length() > 0) {
                    sequences.put(currentAnnotation, currentSequence.toString());
                    currentSequence = new StringBuilder();
                }

                currentAnnotation = line.substring(1); // just ignore the ">"

                continue;
            }

            currentSequence.append(line);
        }

        if (currentSequence.length() > 0) {
            sequences.put(currentAnnotation, currentSequence.toString());
        }

        return sequences;
    }

    /**
     * Retruns the protein annotation for the given peptide or null in case
     * the peptide cannot be mapped to any protein.
     * @param peptideSequence
     * @return
     */
    public String getProteinAnnotation(String peptideSequence) {
        for (String annotation : proteinSequences.keySet()) {
            String proteinSequence = proteinSequences.get(annotation);

            if (proteinSequence.contains(peptideSequence))
                return annotation;
        }

        return null;
    }
}
