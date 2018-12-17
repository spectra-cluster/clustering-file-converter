package uk.ac.ebi.pride.spectracluster.clusteringfileconverter.util;

import uk.ac.ebi.pride.spectracluster.clusteringfilereader.objects.ICluster;
import uk.ac.ebi.pride.spectracluster.clusteringfilereader.objects.IPeptideSpectrumMatch;
import uk.ac.ebi.pride.spectracluster.clusteringfilereader.objects.ISpectrumReference;

import java.util.*;

/**
 * Created by jg on 14.07.14.
 */
public class ClusterUtilities {
    public static final String UNIDENTIFIED_SEQUENCE="UNLDENTLFLED";

    private ICluster currentCluster;

    private String maxSequence;
    private float maxILAngosticRatio;
    /**
     * This ratio is only based on the sequences
     * observed (ie. if a spectrum was identified
     * as 2 sequences, it will be counted twice)
     */
    private float maxILAngosticSequenceRatio;
    private int maxSequenceCount;
    private int nProjects;
    private int nAssays;
    private double mzRange;
    private Set<String> species;
    private Map<String, Integer> sequenceCounts;
    private int charge = 0;
    private IPeptideSpectrumMatch mostCommonPsm;

    private String secondMaxSequence;
    private int secondMaxSequenceCount;

    private String thirdMaxSequence;
    private int thirdMaxSequenceCount;

    private boolean isAverageCharge;

    public ClusterUtilities() {

    }

    public ClusterUtilities(ICluster cluster) {
        processCluster(cluster);
    }

    /**
     * Processed the passed cluster. This call overwrites
     * any results that may have been created before.
     * @param cluster
     */
    public void processCluster(ICluster cluster) {
        // this is only saved for potential future reference
        currentCluster = cluster;

        updateNumberOfProjects(cluster);

        // update the sequences
        updateSequences(cluster);
        updatePrecursorMzRange(cluster);
        updateSpecies(cluster);
        try {
            charge = SpectrumAnnotator.estimateCharge(this.getMostCommonPsm(), cluster.getAvPrecursorMz());
            this.isAverageCharge = false;
        }
        catch (Exception e) {
            charge = this.calculateAverageCharge(cluster);
            this.isAverageCharge = true;
        }
    }

    private void updateSequences(ICluster cluster) {
        if (cluster.getIdentifiedSpecCount() < 1)
            return;

        // update the most common sequence
        List<Object> maxSequenceProperties = getMaxSequence(cluster, Collections.EMPTY_SET);

        this.maxSequence = (String) maxSequenceProperties.get(0);
        this.maxSequenceCount = (Integer) maxSequenceProperties.get(1);
        this.maxILAngosticRatio = (float) this.maxSequenceCount / cluster.getIdentifiedSpecCount();
        this.sequenceCounts = createSequenceCounts(cluster);

        // get the second most common sequence
        String maxIlAgnosticSequence = this.maxSequence.replaceAll("I", "L");
        Set<String> knownSequence = new HashSet<String>();
        knownSequence.add(maxIlAgnosticSequence);

        maxSequenceProperties = getMaxSequence(cluster, knownSequence);

        this.secondMaxSequence = (String) maxSequenceProperties.get(0);
        this.secondMaxSequenceCount = (Integer) maxSequenceProperties.get(1);

        // get the third max sequence
        if (secondMaxSequence != null) {
            knownSequence.add(secondMaxSequence.replace("I", "L"));
            maxSequenceProperties = getMaxSequence(cluster, knownSequence);

            this.thirdMaxSequence = (String) maxSequenceProperties.get(0);
            this.thirdMaxSequenceCount = (Integer) maxSequenceProperties.get(1);
        }
        else {
            this.thirdMaxSequence = null;
            this.thirdMaxSequenceCount = 0;
        }

        // create the max i/l agnostic sequence ratio
        Map<String, Integer> ilSequenceCounts = new HashMap<String, Integer>();
        Map<String, Integer> completeSequenceCounts = createCompleteSequenceCounts(cluster);
        int totalPsms = 0;

        for (String sequence : completeSequenceCounts.keySet()) {
            String ilSequence = sequence.replaceAll("I", "L");

            if (ilSequenceCounts.containsKey(ilSequence)) {
                ilSequenceCounts.put(ilSequence,
                        ilSequenceCounts.get(ilSequence)
                                + completeSequenceCounts.get(sequence));
            }
            else {
                ilSequenceCounts.put(ilSequence, completeSequenceCounts.get(sequence));
            }

            totalPsms += completeSequenceCounts.get(sequence);
        }

        String maxIlSequence = maxSequence.replaceAll("I", "L");
        this.maxILAngosticSequenceRatio = (float) ilSequenceCounts.get(maxIlSequence) /
                (float) totalPsms;
    }

    /**
     * Calculate the average charge based on all spectra
     * in the cluster.
     * @param cluster
     * @return
     */
    private int calculateAverageCharge(ICluster cluster) {
        // calculate the average charge
        int sumCharge = 0;

        int nValidSpecRefs = 0;

        for (ISpectrumReference specRef : cluster.getSpectrumReferences()) {
            if (specRef.getCharge() < 1)
                continue;

            nValidSpecRefs++;
            sumCharge += specRef.getCharge();
        }

        float avCharge = (float) sumCharge / (float) nValidSpecRefs;
        int avChargeRounded = (int) (avCharge + 0.5);

        return avChargeRounded;
    }

    private void updateSpecies(ICluster cluster) {
        species = new HashSet<String>();

        for (ISpectrumReference specRef : cluster.getSpectrumReferences()) {
            species.add(specRef.getSpecies());
        }
    }

    public static String cleanSequence(String sequence) {
        if (sequence == null) {
            return null;
        }

        return sequence.toUpperCase().replaceAll("[^A-Z]", "");
    }

    private Map<String, Integer> createSequenceCounts(ICluster cluster) {
        Map<String, Integer> sequenceCounts = new HashMap<String, Integer>();

        for (ISpectrumReference spectrumReference : cluster.getSpectrumReferences()) {
            if (!spectrumReference.isIdentified())
                continue;

            IPeptideSpectrumMatch peptideSpectrumMatch = spectrumReference.getMostCommonPSM();
            String cleanSequence = cleanSequence(peptideSpectrumMatch.getSequence());

            if (!sequenceCounts.containsKey(cleanSequence))
                sequenceCounts.put(cleanSequence, 1);
            else
                sequenceCounts.put(cleanSequence, sequenceCounts.get(cleanSequence) + 1);
        }

        return sequenceCounts;
    }

    /**
     * Create the complete sequence counts taking all PSMs into consideration.
     * @param cluster
     * @return
     */
    private Map<String, Integer> createCompleteSequenceCounts(ICluster cluster) {
        Map<String, Integer> sequenceCounts = new HashMap<String, Integer>();

        for (ISpectrumReference spectrumReference : cluster.getSpectrumReferences()) {
            if (!spectrumReference.isIdentified())
                continue;

            for (IPeptideSpectrumMatch peptideSpectrumMatch : spectrumReference.getPSMs()) {
                String cleanSequence = cleanSequence(peptideSpectrumMatch.getSequence());

                if (!sequenceCounts.containsKey(cleanSequence))
                    sequenceCounts.put(cleanSequence, 1);
                else
                    sequenceCounts.put(cleanSequence, sequenceCounts.get(cleanSequence) + 1);
            }
        }

        return sequenceCounts;
    }

    /**
     * Retruns the sequence and the max count when ignoring the sequences set in knownMaxSequence
     * @param cluster
     * @param knownMaxSequence
     * @return
     */
    private List<Object> getMaxSequence(ICluster cluster, Set<String> knownMaxSequence) {
        if (cluster.getIdentifiedSpecCount() < 1)
            return null;

        Map<String, Integer> sequenceCounts = new HashMap<String, Integer>();
        Map<String, String> ilCorrectedToOriginalSequence = new HashMap<String, String>();

        for (ISpectrumReference spectrumReference : cluster.getSpectrumReferences()) {
            // only works on identified specRefs
            if (!spectrumReference.isIdentified())
                continue;

            IPeptideSpectrumMatch peptideSpectrumMatch = spectrumReference.getMostCommonPSM();
            String ilAgnosticSequence = peptideSpectrumMatch.getSequence().replaceAll("I", "L");
            ilAgnosticSequence = cleanSequence(ilAgnosticSequence);

            // ignore previously processed sequences
            if (knownMaxSequence.contains(ilAgnosticSequence))
                continue;

            if (!ilCorrectedToOriginalSequence.containsKey(ilAgnosticSequence)) {
                ilCorrectedToOriginalSequence.put(ilAgnosticSequence, peptideSpectrumMatch.getSequence());
            }

            if (!sequenceCounts.containsKey(ilAgnosticSequence)) {
                sequenceCounts.put(ilAgnosticSequence, 0);
            }

            sequenceCounts.put(ilAgnosticSequence, sequenceCounts.get(ilAgnosticSequence) + 1);
        }

        // get the max count
        int maxCount = 0;
        String maxSequence = null;

        for (String ilAgnosticSequence : sequenceCounts.keySet()) {
            int count = sequenceCounts.get(ilAgnosticSequence);

            if (count > maxCount) {
                maxSequence = ilCorrectedToOriginalSequence.get(ilAgnosticSequence);
                maxCount = count;
            }
        }

        List<Object> result = new ArrayList<Object>(2);
        result.add(cleanSequence(maxSequence));
        result.add(new Integer(maxCount));

        return result;
    }

    private void updateNumberOfProjects(ICluster cluster) {
        Set<String> projects = new HashSet<String>();
        Set<String> assays = new HashSet<String>();

        for (ISpectrumReference specRef : cluster.getSpectrumReferences()) {
            String id = specRef.getSpectrumId();

            String[] fields = id.split(";");

            // PXD000807;index=2376,splib_sequence=KYLYEIAR,score=0.639,peptideR2=,scoreR2=
            projects.add(fields[0]);
            assays.add((fields.length >= 3) ? fields[1] : "Unknown");
        }

        this.nProjects = projects.size();
        this.nAssays = assays.size();
    }

    private void updatePrecursorMzRange(ICluster cluster) {
        double minMZ = Float.MAX_VALUE, maxMz = 0;

        for (ISpectrumReference specRef : cluster.getSpectrumReferences()) {
            double mz = specRef.getPrecursorMz();

            if (mz < minMZ) {
                minMZ = mz;
            }
            if (mz > maxMz) {
                maxMz = mz;
            }
        }

        this.mzRange = maxMz - minMZ;
    }

    /**
     * Returns one PSM object matching the most common sequence
     * choosing the PSM object with the lowest delta mass.
     * @return
     */
    public IPeptideSpectrumMatch getMostCommonPsm() throws Exception {
        if (mostCommonPsm != null)
            return mostCommonPsm;

        IPeptideSpectrumMatch psm = null;
        double minDelta = Double.MAX_VALUE;

        for (ISpectrumReference specRef : currentCluster.getSpectrumReferences()) {
            for (IPeptideSpectrumMatch currentPsm : specRef.getPSMs()) {
                String psmSequence = cleanSequence(currentPsm.getSequence());

                if (psmSequence.equals(maxSequence)) {
                    double delta = Math.abs(SpectrumAnnotator.getDeltaMass(currentPsm, currentCluster.getAvPrecursorMz()));

                    if (delta < minDelta) {
                        psm = currentPsm;
                        minDelta = delta;
                    }
                }
            }
        }

        mostCommonPsm = psm;
        return mostCommonPsm;
    }

    public ICluster getCurrentCluster() {
        return currentCluster;
    }

    public String getMaxSequence() {
        return maxSequence;
    }

    public float getMaxILAngosticRatio() {
        return maxILAngosticRatio;
    }

    public int getMaxSequenceCount() {
        return maxSequenceCount;
    }

    public int getnProjects() {
        return nProjects;
    }

    public int getnAssays() {
        return nAssays;
    }

    public double getMzRange() {
        return mzRange;
    }

    public Set<String> getSpecies() {
        return Collections.unmodifiableSet(species);
    }

    public String getSecondMaxSequence() {
        return secondMaxSequence;
    }

    public int getSecondMaxSequenceCount() {
        return secondMaxSequenceCount;
    }

    public Map<String, Integer> getSequenceCounts() {
        return Collections.unmodifiableMap(sequenceCounts);
    }

    public String getThirdMaxSequence() {
        return thirdMaxSequence;
    }

    public int getThirdMaxSequenceCount() {
        return thirdMaxSequenceCount;
    }

    public boolean isStable() {
        if (currentCluster.getIdentifiedSpecCount() >= 10 & maxILAngosticRatio > 0.7)
            return true;

        return false;
    }

    public int getCharge() {
        return charge;
    }

    /**
     * Indicates whether the returned charge is the average charge
     * of all spectra within the cluster or whether it was calculated
     * based on the most common PSM and its modifications.
     * @return True if the charge is based on the average charge.
     */
    public boolean isAverageCharge() {
        return isAverageCharge;
    }

    /**
     * Returns the maximum sequence counting every PSM once (ie. if a
     * spectrum was identified as two peptides, each will be counted
     * separately).
     * @return The maximum ratio based on the PSMs observed.
     */
    public float getMaxILAngosticSequenceRatio() {
        return maxILAngosticSequenceRatio;
    }
}
