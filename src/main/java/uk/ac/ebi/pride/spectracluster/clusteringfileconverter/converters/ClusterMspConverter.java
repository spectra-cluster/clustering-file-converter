package uk.ac.ebi.pride.spectracluster.clusteringfileconverter.converters;

import uk.ac.ebi.pride.spectracluster.clusteringfileconverter.util.ClusterUtilities;
import uk.ac.ebi.pride.spectracluster.clusteringfileconverter.util.ModificationPositionComparator;
import uk.ac.ebi.pride.spectracluster.clusteringfileconverter.util.SpectrumAnnotator;
import uk.ac.ebi.pride.spectracluster.clusteringfilereader.objects.ICluster;
import uk.ac.ebi.pride.spectracluster.clusteringfilereader.objects.IModification;
import uk.ac.ebi.pride.spectracluster.clusteringfilereader.objects.IPeptideSpectrumMatch;
import uk.ac.ebi.pride.utilities.iongen.model.ProductIon;
import uk.ac.ebi.pride.utilities.iongen.model.ProductIonSet;
import uk.ac.ebi.pridemod.ModReader;
import uk.ac.ebi.pridemod.model.PTM;

import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.util.*;

/**
 * Created by jg on 01.08.14.
 */
public class ClusterMspConverter extends AbstractClusterConverter {
    public final static String FILE_EXTENSION = "msp";
    public final static Double MOD_TOLERANCE = 0.001;
    public final static float ANNOTATION_TOLERANCE = 0.5F;
    public final int BASE_PEAK_INTENSITY = 10000; // intensity used for base peak normalization

    private Map<String, Double> modToDeltaMap;
    private boolean normalizeSpectra = false;
    private boolean addAnnotationString = false;
    private ProductIon lastProductIon = null;

    // logging
    private Set<String> unresolvedAccessions = new HashSet<String>();
    private Set<String> missingModMassAccessions = new HashSet<String>();
    private Set<String> unmatchedModificationAccessions = new HashSet<String>();

    public ClusterMspConverter() {
        modToDeltaMap = loadMspModMap();
    }

        @Override
    public String getFileHeader() {
        return "";
    }

    @Override
    public String getFiletypeExtension() {
        return FILE_EXTENSION;
    }

    @Override
    public String convertCluster(ICluster cluster) {
        try {
            ClusterUtilities clusterUtilities = new ClusterUtilities(cluster);

            // TODO: filter clusters where the most common PSM was identified with two different sequences

            if (clusterUtilities.isAverageCharge()) {
                throw new Exception("Failed to calculate charge for cluster " + cluster.getId());
            }

            StringBuilder mspString = new StringBuilder();

            mspString.append("Name: ").append(generateClusterName(clusterUtilities.getMaxSequence(), clusterUtilities.getCharge())).append("\n");

            double molecularWeight = (cluster.getAvPrecursorMz() * clusterUtilities.getCharge()) - (clusterUtilities.getCharge() * 1.008);
            mspString.append("MW: ").append(molecularWeight).append("\n");

            String commentString = generateComments(cluster, clusterUtilities);
            mspString.append("Comment: ").append(commentString).append("\n");

            // normalize intensity values
            List<Float> normalizedIntensities = null;
            if (normalizeSpectra)
                normalizedIntensities = normalizeIntensities(cluster.getConsensusIntensValues(), BASE_PEAK_INTENSITY);

            mspString.append("Num peaks: ").append(getPeakCount(cluster, normalizedIntensities)).append("\n");

            String annotatedPeakList = createPeakList(cluster, clusterUtilities, normalizedIntensities);

            mspString.append(annotatedPeakList);

            return mspString.toString();
        }
        catch (Exception e) {
            e.printStackTrace();
            return "";
        }
    }

    private String createPeakList(ICluster cluster, ClusterUtilities clusterUtilities, List<Float> normalizedIntensities) throws Exception {
        StringBuilder peakListString = new StringBuilder();
        ProductIonSet productIonSet = SpectrumAnnotator.getProductIonSet(cluster, clusterUtilities.getMostCommonPsm(), clusterUtilities.getCharge());

        // create the peak annotation set
        for (int i = 0; i < cluster.getConsensusMzValues().size(); i++) {
            // ignore peaks with 0 m/z and 0 intensity
            if (cluster.getConsensusMzValues().get(i) == 0)
                continue;
            boolean emptyIntensity = (normalizeSpectra) ? normalizedIntensities.get(i) == 0 : cluster.getConsensusIntensValues().get(i) == 0;
            if (emptyIntensity)
                continue;

            String intensityString = (normalizeSpectra) ? normalizedIntensities.get(i).toString() : cluster.getConsensusIntensValues().get(i).toString();

            // get the ion set
            String annotation = "";
            if (addAnnotationString) {
                annotation = "\"" + getBestPeakAnnotation(cluster.getConsensusMzValues().get(i), productIonSet.iterator(), cluster, clusterUtilities) + "\"";
            }

            peakListString.append(
                    cluster.getConsensusMzValues().get(i))
                    .append(" ")
                    .append(intensityString)
                    .append(" ")
                    .append(annotation)
                    .append("\n");
        }

        return peakListString.toString();
    }

    private int getPeakCount(ICluster cluster, List<Float> normalizedIntensities) {
        // count the peaks ignoring any peaks with 0 m/z or intensity
        int nPeaks = 0;
        for (int i = 0; i < cluster.getConsensusMzValues().size(); i++) {
            // ignore peaks with 0 m/z and 0 intensity
            if (cluster.getConsensusMzValues().get(i) == 0) {
                continue;
            }

            if (normalizeSpectra) {
                if (normalizedIntensities.get(i) == 0) {
                    continue;
                }
            } else {
                if (cluster.getConsensusIntensValues().get(i) == 0)
                    continue;
            }

            nPeaks++;
        }

        return nPeaks;
    }

    private String getBestPeakAnnotation(Float mz, Iterator<ProductIon> iterator, ICluster cluster, ClusterUtilities clusterUtilities) {
        // get all possible annotations
        List<String> annotations = getPeakIonAnnotations(mz, iterator, clusterUtilities, cluster);

        if (annotations.size() < 1)
            throw new IllegalStateException("No annotation returned for peak");

        // first parent
        for (String annotation : annotations) {
            if (annotation.contains("p/"))
                return annotation;
        }

        // then y
        for (String annotation : annotations) {
            if (annotation.contains("y") & !annotation.contains("^"))
                return annotation;
        }

        // then b
        for (String annotation : annotations) {
            if (annotation.contains("b") & !annotation.contains("^"))
                return annotation;
        }

        // then y with charge
        for (String annotation : annotations) {
            if (annotation.contains("y"))
                return annotation;
        }

        // then b with charge
        for (String annotation : annotations) {
            if (annotation.contains("b"))
                return annotation;
        }

        return annotations.get(0);
    }

    private List<String> getPeakIonAnnotations(Float mz, Iterator<ProductIon> iterator, ClusterUtilities clusterUtilities, ICluster cluster) {
        List<String> annotations = new ArrayList<String>();

        // TODO: optimize peak calculation...
        int minPeaks = Math.round(cluster.getSpecCount() * 0.7F);
        int correctSpecs = clusterUtilities.getMaxSequenceCount();
        String peakCount = correctSpecs + "/" + minPeaks;

        // TODO: calculate median deviation in m/z * 100
        String deviation = "1.0";

        // check if it's the precursor
        // TODO: add neutral losses to precursor identification
        if (mz < cluster.getAvPrecursorMz() + ANNOTATION_TOLERANCE && mz > cluster.getAvPrecursorMz() - ANNOTATION_TOLERANCE) {
            String annotation =  "p/" + " " + peakCount + " " + String.format(".2f", cluster.getAvPrecursorMz() - mz);
            annotations.add(annotation);
        }

        // TODO: write while clause
        while(lastProductIon != null || iterator.hasNext()) {
            if (lastProductIon == null)
                lastProductIon = iterator.next();

            float lowerMz = mz - ANNOTATION_TOLERANCE;
            float upperMz = mz + ANNOTATION_TOLERANCE;

            if (lastProductIon.getMassOverCharge() < lowerMz) {
                lastProductIon = null; // jump to next ion by setting it to null
                continue;
            }

            // no more to be found
            if (lastProductIon.getMassOverCharge() > upperMz) {
                break;
            }

            // within range, add string to annotation string list
            String annotationString = lastProductIon.getType().toString().toLowerCase() +
                    lastProductIon.getPosition();
            if (lastProductIon.getCharge() > 1)  {
                annotationString += "^" + lastProductIon.getCharge();
            }

            annotationString += "/" + String.format("%.2f", lastProductIon.getMassOverCharge() - mz) +
                        " " + peakCount + " " + deviation;

            annotations.add(annotationString);
            lastProductIon = null;
        }

        // if there is no annotation, mark peak as unknown
        if (annotations.size() < 1) {
            annotations.add("? " + peakCount + " 0.5");
        }

        return annotations;
    }

    private List<Float> normalizeIntensities(List<Float> consensusIntensValues, int fixedHighestPeakIntensity) {
        // get the highest intensity
        float maxIntensity = Collections.max(consensusIntensValues);

        // calculate the factor
        float factor = fixedHighestPeakIntensity / maxIntensity;

        List<Float> normalizedIntensities = new ArrayList<Float>(consensusIntensValues.size());

        for (Float intensity : consensusIntensValues) {
            normalizedIntensities.add(intensity * factor);
        }

        return normalizedIntensities;
    }

    private String generateComments(ICluster cluster, ClusterUtilities clusterUtilities) throws Exception {
        StringBuilder commentString = new StringBuilder();

        String modString = generateModString(cluster, clusterUtilities);

        commentString.append("Spec=Consensus");
        commentString.append(" Mods=" + modString);
        commentString.append(" Parent=").append(String.format("%.3f", cluster.getAvPrecursorMz()));
        commentString.append(" Nreps=").append(cluster.getSpecCount());
        commentString.append(" Naa=").append(clusterUtilities.getMaxSequence().length());
        commentString.append(" MaxRatio=").append(String.format("%.3f", clusterUtilities.getMaxILAngosticRatio()));
        commentString.append(" PrecursorMzRange=").append(String.format("%.4f", cluster.getSpectrumPrecursorMzRange()));
        commentString.append(" DeltaMass=").append(String.format("%.2f", SpectrumAnnotator.getDeltaMass(clusterUtilities.getMostCommonPsm(), cluster.getAvPrecursorMz())));
        if (cluster.getId() != null)
            commentString.append(" ClusterId=").append(cluster.getId());

        if (currentAnnotation != null)
            commentString.append(" Protein=" + currentAnnotation);

        return commentString.toString();
    }

    private String generateModString(ICluster cluster, ClusterUtilities clusterUtilities) throws Exception {
        // get the PSM for the most common sequence, use the one with most modifications annotated
        IPeptideSpectrumMatch psm = clusterUtilities.getMostCommonPsm();

        // this case should never happen
        if (psm == null) {
            throw new IllegalStateException("Failed to retrieve PSM object for most common sequence.");
        }

        // build the modification string based on the psm's modifications
        if (psm.getModifications().size() < 1) {
            return "0";
        }

        ModReader modReader = ModReader.getInstance();

        // sort by position and only use unique modifications
        List<IModification> modifications = new ArrayList<IModification>( new HashSet<IModification>(psm.getModifications()) );
        Collections.sort(modifications, new ModificationPositionComparator());

        StringBuilder modificationString = new StringBuilder();
        modificationString.append(modifications.size());

        for (IModification modification : modifications) {
            String modMspName;

            // first try to get the modification's delta mass
            PTM ptmObject = modReader.getPTMbyAccession(modification.getAccession());

            // if we can't get the delta, we'll use the accession as a name
            if (ptmObject == null) {
                unresolvedAccessions.add(modification.getAccession());
                modMspName = modification.getAccession();
            } else {
                Double delta = ptmObject.getMonoDeltaMass();

                // if there is no delta for the modification (incorrect modification supplied), use the accession
                if (delta == null) {
                    missingModMassAccessions.add(modification.getAccession());
                    modMspName = modification.getAccession();
                } else {
                    // using the delta we try to match it to the MSP name
                    modMspName = getMspNameForDetla(delta);

                    // if the delta cannot be matched, we'll use the delta
                    if (modMspName == null) {
                        unmatchedModificationAccessions.add(modification.getAccession() + " > " + delta);
                        modMspName = delta.toString();
                    }
                }
            }

            // extract the sequence position
            Character sequenceChar;
            String psmSequence = ClusterUtilities.cleanSequence(psm.getSequence());

            if (modification.getPosition() < 1) {
                sequenceChar = '^';
            } else if (modification.getPosition() >= psmSequence.length()) {
                sequenceChar = '$';
            } else {
                sequenceChar = psmSequence.charAt(modification.getPosition() - 1);
            }

            int zeroBasedPosition = (modification.getPosition() > 0) ? modification.getPosition() - 1 : 0;

            // add the modification to the string
            modificationString.append("/" + zeroBasedPosition + "," + sequenceChar + "," + modMspName);
        }

        return modificationString.toString();
    }

    private String getMspNameForDetla(Double delta) {
        for (String mspName : modToDeltaMap.keySet()) {
            // if the delta is within the tolerance, use this name
            Double modDelta = modToDeltaMap.get(mspName);

            if (modDelta > delta - MOD_TOLERANCE && modDelta < delta + MOD_TOLERANCE) {
                return mspName;
            }
        }

        return null;
    }

    private String generateClusterName(String sequence, int charge) {
        return String.format("%s/%d", sequence, charge);
    }

    @Override
    public void onNewClusterRead(ICluster newCluster) {
        if (!shouldClusterBeExported(newCluster))
            return;

        String mspString = convertCluster(newCluster);

        writeStringToFile(mspString + "\n", getFileHeader());
    }

    public boolean isNormalizeSpectra() {
        return normalizeSpectra;
    }

    public void setNormalizeSpectra(boolean normalizeSpectra) {
        this.normalizeSpectra = normalizeSpectra;
    }

    public boolean isAddAnnotationString() {
        return addAnnotationString;
    }

    public void setAddAnnotationString(boolean addAnnotationString) {
        this.addAnnotationString = addAnnotationString;
    }

    private Map<String, Double> loadMspModMap() {
        BufferedReader reader = new BufferedReader(new InputStreamReader(
                ClusterMspConverter.class.getClassLoader().getResourceAsStream("msp_mod_mapping.tsv")));

        Map<String, Double> mspNameToDeltaMap = new HashMap<String, Double>();
        try {
            String line;

            while ((line = reader.readLine()) != null) {
                // remove comments
                int index = line.indexOf('#');
                if (index >= 0) {
                    line = line.substring(0, index);
                }

                // ignore empty lines
                if (line.length() < 1)
                    continue;

                // split based on fields
                String[] fields = line.split("\t");

                if (fields.length < 2)
                    throw new IllegalStateException("Invalid line encountered in msp_mod_mapping.tsv");

                String mspName = fields[0];
                Double deltaMass = Double.parseDouble(fields[1]);

                mspNameToDeltaMap.put(mspName, deltaMass);
            }

            return mspNameToDeltaMap;
        }
        catch (Exception e) {
            throw new IllegalStateException("Failed to load MSP mod name map.");
        }
    }
}
