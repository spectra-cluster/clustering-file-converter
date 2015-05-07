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
public class ZeroPositionTest {
    ClusteringFileReader reader;
    List<ICluster> clusters;
    ClusterMspConverter converter;

    @Before
    public void setUp() throws Exception {
        URI testFileUri = ClusterMspConverterTest.class.getClassLoader().getResource("position_0_mod.clustering").toURI();
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

            String[] lines = mspString.split("\n");
            String commentLine = lines[1].substring(31);

            if (nCluster == 101 || nCluster == 103)
                continue;

            Assert.assertFalse("0 position mod encountered: " + commentLine + " (" + cluster.getId() + " - " + nCluster + ")", commentLine.contains("/0,") && commentLine.charAt(0) != '0');
        }
    }
}
