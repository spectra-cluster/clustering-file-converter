package uk.ac.ebi.pride.spectracluster.clusteringfileconverter;

import junit.framework.Assert;
import org.junit.Before;
import org.junit.Test;
import uk.ac.ebi.pride.spectracluster.clusteringfileconverter.converters.ClusterMspConverter;
import uk.ac.ebi.pride.spectracluster.clusteringfilereader.io.ClusteringFileReader;
import uk.ac.ebi.pride.spectracluster.clusteringfilereader.objects.ICluster;

import java.io.File;
import java.net.URI;
import java.util.List;

/**
 * Created by jg on 02.01.15.
 */
public class MissingAAModTest {
    ClusteringFileReader reader;
    List<ICluster> clusters;
    ClusterMspConverter converter;

    @Before
    public void setUp() throws Exception {
        URI testFileUri = ClusterMspConverterTest.class.getClassLoader().getResource("missing_aa_mod.clustering").toURI();
        reader = new ClusteringFileReader(new File(testFileUri));
        clusters = reader.readAllClusters();
        converter = new ClusterMspConverter();
        converter.setMinSize(10);
        converter.setMinRatio(0.7F);
    }

    @Test
    public void testConversion() throws Exception {
        int nCluster = 0;

        for (ICluster cluster : clusters) {
            nCluster++;

            String mspString = converter.convertCluster(cluster);

            if (mspString.length() < 1)
                continue;

            String[] lines = mspString.split("\n");
            String sequence = lines[0].substring(6, lines[0].length() - 2);
            String commentLine = lines[2].substring(31);

            if (!commentLine.startsWith("P")) {
                String modPosString = commentLine.substring(0, commentLine.indexOf(','));
                if (!modPosString.contains("/")) {
                    int firstModPos = Integer.parseInt(modPosString);
                    Assert.assertTrue("Impossible mod position (" + cluster.getId() + " - " + nCluster + ")", firstModPos < sequence.length());
                }
            }
        }
    }
}
