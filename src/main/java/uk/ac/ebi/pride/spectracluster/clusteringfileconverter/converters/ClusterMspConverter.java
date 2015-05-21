package uk.ac.ebi.pride.spectracluster.clusteringfileconverter.converters;

import uk.ac.ebi.pride.spectracluster.analysis.util.ClusterUtilities;
import uk.ac.ebi.pride.spectracluster.clusteringfileconverter.util.ModificationPositionComparator;
import uk.ac.ebi.pride.spectracluster.clusteringfilereader.objects.*;
import uk.ac.ebi.pride.utilities.iongen.impl.DefaultPrecursorIon;
import uk.ac.ebi.pride.utilities.iongen.ion.FragmentIonUtilities;
import uk.ac.ebi.pride.utilities.iongen.model.*;
import uk.ac.ebi.pride.utilities.mol.Element;
import uk.ac.ebi.pride.utilities.mol.Group;
import uk.ac.ebi.pride.utilities.mol.PTModification;
import uk.ac.ebi.pride.utilities.mol.Peptide;
import uk.ac.ebi.pridemod.ModReader;
import uk.ac.ebi.pridemod.model.PTM;

import java.util.*;

/**
 * Created by jg on 01.08.14.
 */
public class ClusterMspConverter extends AbstractClusterConverter {
    public final static String FILE_EXTENSION = "msp";
    public final static Double MOD_TOLERANCE = 0.01;
    public final static float ANNOTATION_TOLERANCE = 0.8F; // as per the MSP file specification
    public final int BASE_PEAK_INTENSITY = 10000; // intensity used for base peak normalization

    private Map<String, Double> modToDeltaMap = new HashMap<String, Double>();
    private ModReader modReader = ModReader.getInstance();
    private boolean normalizeSpectra = false;
    private boolean addAnnotationString = false;
    private ProductIon lastProductIon = null;

    public ClusterMspConverter() {
        loadModToDeltaMap();
    }

    /**
     * This code was extracted from the SpectraST source code. I was
     * unable to find these modifications in the MSP format
     * specification.
     */
    private void loadModToDeltaMap() {
        modToDeltaMap.put("ICAT_light", 227.126991);
        modToDeltaMap.put("ICAT-C", 227.126991); // PSI new name

        modToDeltaMap.put("ICAT_heavy", 236.157185);
        modToDeltaMap.put("ICAT-C:13C(9)", 236.157185); // PSI new name

        modToDeltaMap.put("AB_old_ICATd0", 442.224991);
        modToDeltaMap.put("ICAT-D", 442.224991); // PSI new name

        modToDeltaMap.put("AB_old_ICATd8", 450.275205);
        modToDeltaMap.put("ICAT-D:2H(8)", 450.275205); // PSI new name


        modToDeltaMap.put("Carbamidomethyl", 57.021464);

        modToDeltaMap.put("Carboxymethyl", 58.005479);

        modToDeltaMap.put("Propionamide", 71.037114); // alkylation of acrylamide to cysteines
        modToDeltaMap.put("Propionamide:2H(3)", 74.055944); // alkylation of heavy acrylamide to cysteines
        modToDeltaMap.put("Propionamide:13C(3)", 74.047178); // alkylation of heavy acrylamide to cysteines

        modToDeltaMap.put("Oxidation", 15.99491);

        modToDeltaMap.put("Acetyl", 42.010565); // acetylation of N terminus

        modToDeltaMap.put("Deamidation", 0.984016);
        modToDeltaMap.put("Deamidated", 0.984016); // PSI new name

        modToDeltaMap.put("Pyro-cmC", 39.994915); // cyclicization of N-terminal CAM-cysteine (FIXED value 01/27/07)
        modToDeltaMap.put("Pyro-carbamidomethyl", 39.994915); // PSI new name

        modToDeltaMap.put("Pyro-glu", -17.026549); // loss of NH3 from glutamine
        modToDeltaMap.put("Gln->pyro-Glu", -17.026549); // PSI new name

        modToDeltaMap.put("Pyro_glu", -18.010565); // loss of H2O from glutamic acid
        modToDeltaMap.put("Glu->pyro-Glu", -18.010565); // PSI new name

        modToDeltaMap.put("Amide", -0.984016); // amidation of C terminus
        modToDeltaMap.put("Amidated", -0.984016); // PSI new name

        modToDeltaMap.put("Phospho", 79.966331); // phosphorylation

        modToDeltaMap.put("Thiophospho", 95.943487); // phosphorylation

        modToDeltaMap.put("Sulfo", 79.956815); // O-sulfonation

        modToDeltaMap.put("Methyl", 14.015650); // methylation

        //  modToDeltaMap.put("Deimination", 0.984016); // deamidation on R

        modToDeltaMap.put("Carbamyl", 43.005814); // carbamylation of N terminus or lysines

        modToDeltaMap.put("iTRAQ4plex", 144.102063); // iTRAQ 4-plex

        modToDeltaMap.put("iTRAQ4plexAcetyl", 186.112628); // iTRAQ 4-plex

        modToDeltaMap.put("iTRAQ8plex:13C(6)15N(2)", 304.19904); // iTRAQ on N terminus or K
        modToDeltaMap.put("iTRAQ8plex", 304.20536); // iTRAQ on N terminus or K


        modToDeltaMap.put("TMT6plex", 229.162932); // TMT 6-plex

        modToDeltaMap.put("PEO-Iodoacetyl-LC-Biotin", 414.52); // Hui Zhang's PEO alkylation agent on cysteines

        modToDeltaMap.put("Label:2H(3)", 3.018830); // SILAC heavy leucine (+3)

        modToDeltaMap.put("Label:2H(4)", 4.025107); // Lys4 label (+4)

        modToDeltaMap.put("Label:13C(6)", 6.020129); // SILAC heavy lysine and arginine (+6)
        modToDeltaMap.put("Label:13C(6)15N(1)", 7.017165);
        modToDeltaMap.put("Label:13C(6)15N(2)", 8.014199); // SILAC heavy lysine (+8)
        modToDeltaMap.put("Label:13C(6)15N(3)", 9.011235);
        modToDeltaMap.put("Label:13C(6)15N(4)", 10.008269); // SILAC heavy arginine (+10)

        modToDeltaMap.put("Methylthio", 45.987721); // methylthiolated cysteine (cys blocking by MMTS)

        modToDeltaMap.put("Leucyl", 113.08406); // leucine added to N-term or K
        modToDeltaMap.put("Leucyl:13C(6)15N(1)", 120.101224); // heavy leucine added to N-term or K


        modToDeltaMap.put("Nitro", 44.985078);
        modToDeltaMap.put("Dimethyl", 28.031300);
        modToDeltaMap.put("Trimethyl", 42.046950);

        modToDeltaMap.put("Bromo", 77.910511);

        // Ubl chains
        modToDeltaMap.put("SUMO_1", 2135.920495); // SUMO-1 Tryptic/LysC tail
        modToDeltaMap.put("SUMO_2_3_Tryp", 3549.536567); // SUMO-2/3 Tryptic tail
        modToDeltaMap.put("Smt3_R93A_Tryp", 3812.747563); // Smt3_R93A Tryptic tail
        modToDeltaMap.put("Smt3_R93A_LysC", 4544.074787); // Smt3_R93A LysC tail
        modToDeltaMap.put("NEDD8_LysC", 1555.956231); // NEDD8 LysC tail
        modToDeltaMap.put("Rub1_LysC", 2454.341699); // Rub1 LysC tail
        modToDeltaMap.put("Ub_LysC", 1431.831075); // Ubiquitin LysC tail
        modToDeltaMap.put("GlyGly", 114.042927); // Ubiquitin/NEDD8 Tryptic tail (2 glycines)

        // added based on PSI-MOD entries
        modToDeltaMap.put("Formyl", 27.994915);
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
        ClusterUtilities clusterUtilities = new ClusterUtilities(cluster);

        StringBuilder mspString = new StringBuilder();

        mspString.append("Name: ").append(generateClusterName(clusterUtilities)).append("\n");
        // TODO: calculate molecular weight based on (m/z * charge) - (charge * 1.008) but this line is optional
        // mspString.append("MW: ").append(cluster.getAvPrecursorMz()).append("\n");

        String commentString = generateComments(cluster, clusterUtilities);
        mspString.append("Comment: ").append(commentString).append("\n");

        // normalize intensity values
        List<Integer> normalizedIntensities = null;
        if (normalizeSpectra)
            normalizedIntensities = normalizeIntensities(cluster.getConsensusIntensValues(), BASE_PEAK_INTENSITY);

        // get the product ion set
        ProductIonSet productIonSet = getProductIonSet(cluster, clusterUtilities);

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
            }
            else {
                if (cluster.getConsensusIntensValues().get(i) == 0)
                    continue;
            }

            nPeaks++;
        }

        mspString.append("Num peaks: ").append(nPeaks).append("\n");

        Iterator<ProductIon> productIonIterator = productIonSet.iterator();
        lastProductIon = null;

        // create the peak annotation set
        for (int i = 0; i < cluster.getConsensusMzValues().size(); i++) {
            // ignore peaks with 0 m/z and 0 intensity
            if (cluster.getConsensusMzValues().get(i) == 0)
                continue;
            boolean emptyIntensity = (normalizeSpectra) ? normalizedIntensities.get(i) == 0 : cluster.getConsensusIntensValues().get(i) == 0;
            if (emptyIntensity)
                continue;

            String intensityString = (normalizeSpectra) ? normalizedIntensities.get(i).toString() : cluster.getConsensusIntensValues().get(i).toString();

            // TODO: add peak annotation field (in quotes, fields: [ion = ?/y/b...],

            // get the ion set
            String annotation = "";
            if (addAnnotationString) {
                annotation = "\"" + getBestPeakAnnotation(cluster.getConsensusMzValues().get(i), productIonIterator, cluster, clusterUtilities) + "\"";
            }

            mspString.append(
                    cluster.getConsensusMzValues().get(i))
                    .append(" ")
                    .append(intensityString)
                    .append(" ")
                    .append(annotation)
                    .append("\n");
        }

        return mspString.toString();
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

    private ProductIonSet getProductIonSet(ICluster cluster, ClusterUtilities clusterUtilities) {
        // build the peak set
        double mzValues[] = new double[cluster.getConsensusMzValues().size()];
        double intensValues[] = new double[cluster.getConsensusMzValues().size()];

        for (int i = 0; i < cluster.getConsensusMzValues().size(); i++) {
            mzValues[i] = new Double(cluster.getConsensusMzValues().get(i));
            intensValues[i] = new Double(cluster.getConsensusIntensValues().get(i));
        }

        PeakSet peakSet = PeakSet.getInstance(mzValues, intensValues);

        // get the max psm ref
        IPeptideSpectrumMatch maxPsm = getMaxSpecRef(cluster.getSpectrumReferences(), clusterUtilities.getMaxSequence());

        // build the modification list
        Map<Integer, PTModification> modifications = new HashMap<Integer, PTModification>();
        Group nTerm = null, cTerm = null;

        // TODO: find out how to report n- / c-term mods
        for (IModification mod : maxPsm.getModifications()) {
            if (mod.getPosition() == 0 || mod.getPosition() > maxPsm.getSequence().length())
                continue;

            int nPosition = mod.getPosition() - 1; // 0-based positions
            PTM ptm = modReader.getPTMbyAccession(mod.getAccession());

            if (ptm == null || ptm.getMonoDeltaMass() == null)
                continue;

            List<Double> avgMass = new ArrayList<Double>(1);
            avgMass.add(ptm.getAveDeltaMass());
            List<Double> monoMass = new ArrayList<Double>(1);
            monoMass.add(ptm.getMonoDeltaMass());
            PTModification ptModification = new PTModification(ptm.getName(), ptm.getName(), ptm.getDescription(), monoMass, avgMass);

            modifications.put(nPosition, ptModification);
        }

        // build the ion set
        int charge = clusterUtilities.getCharge();
        if (charge < 1)
            charge = 1;

        Peptide peptide = new Peptide(clusterUtilities.getMaxSequence(), modifications); // TODO: add n / c term mods
        PrecursorIon precursorIon = new DefaultPrecursorIon(peptide, charge);
        PeptideScore peptideScore = new PeptideScore(precursorIon, peakSet);

        return peptideScore.getProductIonSet();
    }

    private IPeptideSpectrumMatch getMaxSpecRef(List<ISpectrumReference> spectrumReferences, String maxSequence) {
        // TODO: optimize this function - there might be incorrectly reported PTMs
        int maxPtms = -1;
        IPeptideSpectrumMatch maxPSM = null;

        for (ISpectrumReference specRef : spectrumReferences) {
            for (IPeptideSpectrumMatch psm : specRef.getPSMs()) {
                String sequence = ClusterUtilities.cleanSequence(psm.getSequence());
                if (!sequence.equals(maxSequence))
                    continue;

                if (psm.getModifications().size() > maxPtms) {
                    maxPtms = psm.getModifications().size();
                    maxPSM = psm;
                }
            }
        }

        return maxPSM;
    }

    private List<Integer> normalizeIntensities(List<Float> consensusIntensValues, int fixedHighestPeakIntensity) {
        // get the highest intensity
        float maxIntensity = Collections.max(consensusIntensValues);

        // calculate the factor
        float factor = fixedHighestPeakIntensity / maxIntensity;

        List<Integer> normalizedIntensities = new ArrayList<Integer>(consensusIntensValues.size());

        for (Float intensity : consensusIntensValues) {
            // TODO: warning: roundings destroy several spectra (if the precursor wasn't removed f.e.)
            normalizedIntensities.add(Math.round(intensity * factor));
        }

        return normalizedIntensities;
    }

    private String generateComments(ICluster cluster, ClusterUtilities clusterUtilities) {
        StringBuilder commentString = new StringBuilder();

        String modString = generateModString(cluster, clusterUtilities);

        commentString.append("Spec=Consensus");
        commentString.append(" Mods=" + modString);
        commentString.append(" Parent=").append(String.format("%.3f", cluster.getAvPrecursorMz()));
        commentString.append(" Nreps=").append(cluster.getSpecCount());
        commentString.append(" Naa=").append(clusterUtilities.getMaxSequence().length());
        commentString.append(" MaxRatio=").append(String.format("%.3f", clusterUtilities.getMaxILAngosticRatio()));
        commentString.append(" PrecursorMzRange=").append(String.format("%.4f", cluster.getSpectrumPrecursorMzRange()));
        if (currentAnnotation != null)
            commentString.append(" Protein=" + currentAnnotation);
        // TODO: add PRIDE Cluster specific fields

        return commentString.toString();
    }

    private String generateModString(ICluster cluster, ClusterUtilities clusterUtilities) {
        // get the PSM for the most common sequence, use the one with most modifications annotated
        String mostCommonSequence = clusterUtilities.getMaxSequence();
        IPeptideSpectrumMatch psm = null;
        int maxModNum = Integer.MIN_VALUE;

        for (ISpectrumReference specRef : cluster.getSpectrumReferences()) {
            for (IPeptideSpectrumMatch currentPsm : specRef.getPSMs()) {
                String psmSequence = ClusterUtilities.cleanSequence(currentPsm.getSequence());

                if (psmSequence.equals(mostCommonSequence) && currentPsm.getModifications().size() > maxModNum) {
                    psm = currentPsm;
                    maxModNum = currentPsm.getModifications().size();
                }
            }
        }

        // this case should never happen
        if (psm == null) {
            throw new IllegalStateException("Failed to retrieve PSM object for most common sequence.");
        }

        // build the modification string based on the psm's modifications
        if (psm.getModifications().size() < 1) {
            return "0";
        }

        StringBuilder modificationString = new StringBuilder();
        modificationString.append(psm.getModifications().size());

        // sort by position and only use unique modifications
        List<IModification> modifications = new ArrayList<IModification>( new HashSet<IModification>(psm.getModifications()) );
        Collections.sort(modifications, new ModificationPositionComparator());


        for (IModification modification : modifications) {
            String modMspName;

            // first try to get the modification's delta mass
            PTM ptmObject = modReader.getPTMbyAccession(modification.getAccession());

            // if we can't get the delta, we'll use the accession as a name
            if (ptmObject == null) {
                System.out.println("Warning: Failed to resolve modification " + modification.getAccession());
                modMspName = modification.getAccession();
            } else {
                Double delta = ptmObject.getMonoDeltaMass();

                // if there is no delta for the modification (incorrect modification supplied), use the accession
                if (delta == null) {
                    System.out.println("Warning: Modification accession supplied that is not associated with a mass delta: " + modification.getAccession());
                    modMspName = modification.getAccession();
                } else {
                    // using the delta we try to match it to the MSP name
                    modMspName = getMspNameForDetla(delta);

                    // if the delta cannot be matched, we'll use the delta
                    if (modMspName == null) {
                        System.out.println("Warning: Failed to match modification delta (" + modification.getAccession() + " > " + delta + ") to MSP name.");
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

    private String generateClusterName(ClusterUtilities clusterUtilities) {
        return String.format("%s/%d", clusterUtilities.getMaxSequence(), clusterUtilities.getCharge());
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
}
