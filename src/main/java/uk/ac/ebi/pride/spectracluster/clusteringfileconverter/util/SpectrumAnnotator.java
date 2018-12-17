package uk.ac.ebi.pride.spectracluster.clusteringfileconverter.util;

import uk.ac.ebi.pride.spectracluster.clusteringfilereader.objects.ICluster;
import uk.ac.ebi.pride.spectracluster.clusteringfilereader.objects.IModification;
import uk.ac.ebi.pride.spectracluster.clusteringfilereader.objects.IPeptideSpectrumMatch;
import uk.ac.ebi.pride.utilities.exception.IllegalAminoAcidSequenceException;
import uk.ac.ebi.pride.utilities.iongen.impl.DefaultPrecursorIon;
import uk.ac.ebi.pride.utilities.iongen.model.*;
import uk.ac.ebi.pride.utilities.mol.Group;
import uk.ac.ebi.pride.utilities.mol.MoleculeUtilities;
import uk.ac.ebi.pride.utilities.mol.PTModification;
import uk.ac.ebi.pride.utilities.mol.Peptide;
import uk.ac.ebi.pridemod.ModReader;
import uk.ac.ebi.pridemod.model.PTM;

import java.util.*;

/**
 * Created by jg on 21.10.15.
 */
public class SpectrumAnnotator {
    private SpectrumAnnotator() {

    }

    public static float getAnnotatedTic(ICluster cluster, float fragmentIonTolerance) throws Exception {
        ClusterUtilities clusterUtilities = new ClusterUtilities(cluster);

        IPeptideSpectrumMatch maxPsm = clusterUtilities.getMostCommonPsm();

        if (maxPsm == null)
            return 0;

        int estimatedCharge = clusterUtilities.getCharge();

        ProductIonSet productIonSet = getProductIonSet(cluster, maxPsm, estimatedCharge);

        if (productIonSet == null)
            return 0;

        List<Integer> annotatedPeakIndices = getMatchedPeakIndices(cluster, productIonSet, fragmentIonTolerance);

        double totalTIC = 0;
        double annotatedTIC = 0;

        for (int peakIndex = 0; peakIndex < cluster.getConsensusIntensValues().size(); peakIndex++) {
            float intensity = cluster.getConsensusIntensValues().get(peakIndex);

            totalTIC += intensity;

            if (annotatedPeakIndices.contains(peakIndex)) {
                annotatedTIC += intensity;
            }
        }

        return (float) annotatedTIC / (float) totalTIC;
    }

    /**
     * Returns the peaks' annotations in the form of ProductIonS. The key
     * is the 0-based peak's index.
     * @param cluster
     * @param tolerance
     * @param clusterUtilities
     * @return
     * @throws Exception
     */
    public static Map<Integer, ProductIon> getPeakProductIons(ICluster cluster, float tolerance, ClusterUtilities clusterUtilities) throws Exception {
        ProductIonSet productIonSet = getProductIonSet(cluster, clusterUtilities.getMostCommonPsm(), clusterUtilities.getCharge());

        Map<Integer, ProductIon> peakProductIons = new HashMap<Integer, ProductIon>();
        List<Float> sortedMz = new ArrayList<Float>(cluster.getConsensusMzValues());
        Collections.sort(sortedMz);

        int lastIndex = 0;

        // iterate over all possible product ions and try to match them against the peaks
        // - product ions are sorted according to m/z
        for (ProductIon productIon : productIonSet) {
            double productMz = productIon.getMassOverCharge();

            for (int i = lastIndex; i < sortedMz.size(); i++) {
                float peakMz = sortedMz.get(i);

                // if it's too low, go on
                if (peakMz < productMz - tolerance)
                    continue;

                // if the peak is outside the range, it's no match
                if (peakMz > productMz + tolerance)
                    break;

                // it's a match
                peakProductIons.put(i, productIon);
                lastIndex = i + 1;
                break;
            }
        }

        return peakProductIons;
    }

    /**
     * Returns the indices of all peaks that are explained by a
     * b- or y-ion.
     * @param cluster
     * @param productIons
     * @param tolerance
     * @return
     */
    private static List<Integer> getMatchedPeakIndices(ICluster cluster, ProductIonSet productIons, float tolerance) {
        List<Integer> matchedPeaksIndices = new ArrayList<Integer>();
        List<Float> sortedMz = new ArrayList<Float>(cluster.getConsensusMzValues());
        Collections.sort(sortedMz);

        int lastIndex = 0;

        // iterate over all possible product ions and try to match them against the peaks
        // - product ions are sorted according to m/z
        for (ProductIon productIon : productIons) {
            // only match b- and y- ions
            if (!productIon.getType().getName().equals("b") && !productIon.getType().getName().equals("y"))
                continue;

            double productMz = productIon.getMassOverCharge();

            for (int i = lastIndex; i < sortedMz.size(); i++) {
                float peakMz = sortedMz.get(i);

                // if it's too low, go on
                if (peakMz < productMz - tolerance)
                    continue;

                // if the peak is outside the range, it's no match
                if (peakMz > productMz + tolerance)
                    break;

                // it's a match
                matchedPeaksIndices.add(i);
                lastIndex = i + 1;
                break;
            }
        }

        return matchedPeaksIndices;
    }

    /**
     * Estimates the charge by calculating the theoretical mass (including PTMs)
     * and deviding it by the observed precursor mass.
     * @param psm The IPeptideSpectrumMatch representing the identification
     * @param precursorMz The measured precursor m/z
     * @return The calculated charge.
     * @throws Exception
     */
    public static int estimateCharge(IPeptideSpectrumMatch psm, double precursorMz) throws Exception {
        double calculatedMass = calculateTheoreticalMass(psm);
        return Math.round((float) calculatedMass / (float) precursorMz);
    }

    /**
     * Estimate the cluster's charge by dividing the theoretical
     * mass by the precursor m/z.
     * @param sequence
     * @param precursorMz
     * @return
     */
    @Deprecated
    public static int estimateCharge(String sequence, double precursorMz) throws Exception {
        sequence = ClusterUtilities.cleanSequence(sequence);

        try {
            double estimatedMass = MoleculeUtilities.calculateTheoreticalMass(sequence);
            return Math.round((float) (estimatedMass / precursorMz));
        }
        catch(IllegalAminoAcidSequenceException e) {
            throw new Exception(e);
        }
    }

    private static double[] getPtmMasses(IPeptideSpectrumMatch psm) throws Exception {
        // add all modification delta masses to a mod array
        ModReader modReader = ModReader.getInstance();

        double[] ptmMasses = new double[psm.getModifications().size()];
        for (int i = 0; i < psm.getModifications().size(); i++) {
            String modAccession = psm.getModifications().get(i).getAccession();

            // fix a typo in the annotated data
            if ("MOD:010900".equals(modAccession))
                modAccession = "MOD:01090";

            PTM ptm = modReader.getPTMbyAccession(modAccession);

            if (ptm == null)
                throw new Exception("Unknown PTM encountered: " + modAccession);
            if (ptm.getMonoDeltaMass() == null)
                throw new Exception("No mass specified for " + modAccession);

            ptmMasses[i] = ptm.getMonoDeltaMass();
        }

        return ptmMasses;
    }

    public static Double calculateTheoreticalMass(IPeptideSpectrumMatch psm) throws Exception {
        double[] ptmMasses = getPtmMasses(psm);

        double estimatedMass = MoleculeUtilities.calculateTheoreticalMass(psm.getSequence(), ptmMasses);

        return estimatedMass;
    }

    public static Double getDeltaMass(IPeptideSpectrumMatch psm, double precursorMz) throws Exception {
        try {
            // first estimate the charge
            int charge = estimateCharge(psm, precursorMz);
            if (charge == 0)
                throw new Exception("Failed to estimate charge");

            double[] ptmMasses = getPtmMasses(psm);
            List<Double> ptmMassList = new ArrayList<Double>(ptmMasses.length);
            for (double mass : ptmMasses)
                ptmMassList.add(mass);

            return MoleculeUtilities.calculateDeltaMz(psm.getSequence(), precursorMz, charge, ptmMassList);
        }
        catch (Exception e) {
            throw e;
        }
    }

    /**
     * Creates all possible product ions for the given peptide.
     * @param cluster
     * @param maxPsm
     * @param charge
     * @return
     */
    public static ProductIonSet getProductIonSet(ICluster cluster, IPeptideSpectrumMatch maxPsm, int charge) throws Exception {
        // build the peak set
        double mzValues[] = new double[cluster.getConsensusMzValues().size()];
        double intensValues[] = new double[cluster.getConsensusMzValues().size()];

        // convert into Doubles
        for (int i = 0; i < cluster.getConsensusMzValues().size(); i++) {
            mzValues[i] = new Double(cluster.getConsensusMzValues().get(i));
            intensValues[i] = new Double(cluster.getConsensusIntensValues().get(i));
        }

        PeakSet peakSet = PeakSet.getInstance(mzValues, intensValues);

        // build the modification list
        Map<Integer, PTModification> modifications = buildModificationMap(maxPsm);

        try {
            Peptide peptide = new Peptide(maxPsm.getSequence(), modifications);
            PrecursorIon precursorIon = new DefaultPrecursorIon(peptide, charge);
            PeptideScore peptideScore = new PeptideScore(precursorIon, peakSet);

            return peptideScore.getProductIonSet();
        }
        catch (IllegalArgumentException e) {
            throw new Exception(e);
        }
    }

    /**
     * Extracts the defined PTMs from the IPeptideSpectrumMatch and
     * returns a map with the AA index as key and the modification as
     * value.
     * @param psm
     * @return
     */
    private static Map<Integer, PTModification> buildModificationMap(IPeptideSpectrumMatch psm) throws Exception {
        ModReader modReader = ModReader.getInstance();
        Map<Integer, PTModification> modificationMap = new HashMap<Integer, PTModification>();
        Group nTerm = null, cTerm = null;

        for (IModification mod : psm.getModifications()) {
            int nPosition = mod.getPosition() - 1; // 0-based positions

            // put n- / c-term mods on the first or last AA
            if (mod.getPosition() == 0)
                nPosition = 0;
            if (mod.getPosition() > psm.getSequence().length())
                nPosition = nPosition - 1;

            PTM ptm = modReader.getPTMbyAccession(mod.getAccession());

            if (ptm == null || ptm.getMonoDeltaMass() == null)
                throw new Exception("Failed to resolve modification " + mod.getAccession());

            List<Double> avgMass = new ArrayList<Double>(1);
            avgMass.add(ptm.getAveDeltaMass());
            List<Double> monoMass = new ArrayList<Double>(1);
            monoMass.add(ptm.getMonoDeltaMass());
            PTModification ptModification = new PTModification(ptm.getName(), ptm.getName(), ptm.getDescription(), monoMass, avgMass);

            modificationMap.put(nPosition, ptModification);
        }

        return modificationMap;
    }
}
