package uk.ac.ebi.pride.spectracluster.clusteringfileconverter;

import org.junit.Test;
import uk.ac.ebi.pride.utilities.iongen.impl.DefaultPeptideIon;
import uk.ac.ebi.pride.utilities.iongen.impl.DefaultPrecursorIon;
import uk.ac.ebi.pride.utilities.iongen.model.*;
import uk.ac.ebi.pride.utilities.mol.Peptide;

import java.io.File;
import java.util.Iterator;

/**
 * Created by jg on 08.01.15.
 */
public class FragmentIonTest {
    private final static File mgfTestFile = new File("/home/jg/Projects/ebi-pride/cluster-result-assessor/src/test/resources/spectra-size_100-ratio_0.5/spectra-size_100-ratio_0.541f624a0-28fe-4811-8d9e-d2aabf54f10f.mgf");

    @Test
    public void testFragmentIonGeneration() {
        Peptide peptide = new Peptide("EDKTLQTPR");
        DefaultPeptideIon defaultPeptideIon = new DefaultPeptideIon(peptide, 2);

        double[] mzArray = new double[] {203.976, 173.07, 220.623, 260.248, 216.118, 301.019, 237.7, 283.045, 233.076, 160.99, 271.32, 286.199, 245.364, 215.105, 258.232, 199.221, 287.183, 259.059, 247.027, 244.707, 295.223, 265.188, 284.162, 159.012, 223.139, 171.096, 126.822, 123.109, 189.18, 116.397, 255.153, 246.263, 225.941, 256.875, 230.967, 272.975, 241.083, 110.002, 193.144, 238.813, 235.268, 240.014, 105.045, 213.205, 253.002, 167.017, 243.237, 209.064, 206.078, 228.166};
        double[] intensArray = new double[] {2.917, 6.853, 12.344, 3.567, 8.23, 0, 4.251, 74.218, 7.26, 7.609, 4.348, 7.654, 27.731, 36.632, 25.16, 9.381, 6.415, 24.669, 3.44, 14.666, 8.063, 7.837, 7.367, 7.073, 6.159, 7.397, 5.768, 3.161, 15.661, 2.78, 47.395, 18.168, 3.034, 29.669, 10.494, 30.227, 11.456, 5.48, 2.286, 3.31, 4.087, 6.44, 8.767, 3.83, 3.986, 4.441, 16.507, 9.97, 4.692, 7.048};

        PeakSet peakSet = PeakSet.getInstance(mzArray, intensArray);

        System.out.println(defaultPeptideIon.toString());

        PrecursorIon precursorIon = new DefaultPrecursorIon(peptide, 2);
        PeptideScore peptideScore = new PeptideScore(precursorIon, peakSet);

        ProductIonSet productIonSet = peptideScore.getProductIonSet();

        Iterator<ProductIon> productIonIterator = productIonSet.iterator();

        System.out.println(productIonSet);

        while (productIonIterator.hasNext()) {
            ProductIon productIon = productIonIterator.next();


            System.out.println(productIon + " " +
                    productIon.getType().toString() +
                    productIon.getCharge() +
                    productIon.getPosition() + " " +
                    productIon.getMassOverCharge()
            );
        }
    }
}
