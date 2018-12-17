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
 * Created by jg on 03.01.15.
 */
public class DuplicatedModTest {
    ClusteringFileReader reader;
    List<ICluster> clusters;
    ClusterMspConverter converter;

    @Before
    public void setUp() throws Exception {
        URI testFileUri = ClusterMspConverterTest.class.getClassLoader().getResource("duplicated_mod.clustering").toURI();
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
            String commentLine = lines[2];

            Assert.assertFalse(cluster.getId() + " - " + nCluster, commentLine.contains("0,K,Acetyl/0,K,Acetyl"));
        }
    }
}
