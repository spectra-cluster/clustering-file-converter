package uk.ac.ebi.pride.spectracluster.clusteringfileconverter;

import junit.framework.Assert;
import org.junit.Before;
import org.junit.Test;
import uk.ac.ebi.pride.spectracluster.clusteringfileconverter.converters.ClusterMspConverter;
import uk.ac.ebi.pride.spectracluster.clusteringfilereader.io.ClusteringFileReader;
import uk.ac.ebi.pride.spectracluster.clusteringfilereader.io.IClusterSourceListener;
import uk.ac.ebi.pride.spectracluster.clusteringfilereader.objects.ICluster;

import java.io.File;
import java.net.URI;
import java.util.ArrayList;
import java.util.List;

/**
 * Created by jg on 09.08.14.
 */
public class ClusterMspConverterTest {
    ClusteringFileReader reader;
    List<ICluster> clusters;
    ClusterMspConverter converter;

    @Before
    public void setUp() throws Exception {
        URI testFileUri = ClusterMspConverterTest.class.getClassLoader().getResource("testfile.clustering").toURI();
        reader = new ClusteringFileReader(new File(testFileUri));
        clusters = reader.readAllClusters();
        converter = new ClusterMspConverter();
    }

    @Test
    public void testConversion() {
        converter.setAddAnnotationString(true);
        converter.setNormalizeSpectra(true);
        String mspCluster10 = converter.convertCluster(clusters.get(10));

        Assert.assertEquals("Name: KNYGK/0\n" +
                "Comment: Spec=Consensus Mods=0 Parent=305.010 Nreps=1 Naa=5 MaxRatio=1.000 PrecursorMzRange=0.0000\n" +
                "Num peaks: 33\n" +
                "93.084 800 \"? 1/1 0.5\"\n" +
                "95.164 725 \"? 1/1 0.5\"\n" +
                "110.259 5043 \"? 1/1 0.5\"\n" +
                "129.109 2071 \"b1/-0.01 1/1 1.0\"\n" +
                "171.111 1752 \"? 1/1 0.5\"\n" +
                "175.276 2674 \"? 1/1 0.5\"\n" +
                "191.088 2614 \"? 1/1 0.5\"\n" +
                "259.323 4537 \"? 1/1 0.5\"\n" +
                "262.339 3797 \"? 1/1 0.5\"\n" +
                "263.435 3425 \"? 1/1 0.5\"\n" +
                "269.336 3398 \"? 1/1 0.5\"\n" +
                "270.288 3059 \"? 1/1 0.5\"\n" +
                "319.355 3158 \"? 1/1 0.5\"\n" +
                "320.215 3173 \"? 1/1 0.5\"\n" +
                "359.311 3278 \"? 1/1 0.5\"\n" +
                "370.553 8321 \"? 1/1 0.5\"\n" +
                "391.858 10000 \"? 1/1 0.5\"\n" +
                "401.437 1715 \"? 1/1 0.5\"\n" +
                "406.148 9367 \"b3/0.06 1/1 1.0\"\n" +
                "406.774 1861 \"? 1/1 0.5\"\n" +
                "421.269 3029 \"? 1/1 0.5\"\n" +
                "425.916 1886 \"? 1/1 0.5\"\n" +
                "502.375 387 \"? 1/1 0.5\"\n" +
                "508.51 347 \"? 1/1 0.5\"\n" +
                "525.282 409 \"? 1/1 0.5\"\n" +
                "554.575 535 \"? 1/1 0.5\"\n" +
                "589.398 357 \"? 1/1 0.5\"\n" +
                "601.327 271 \"? 1/1 0.5\"\n" +
                "602.224 258 \"? 1/1 0.5\"\n" +
                "603.364 1033 \"? 1/1 0.5\"\n" +
                "611.135 206 \"? 1/1 0.5\"\n" +
                "630.366 360 \"? 1/1 0.5\"\n" +
                "737.475 196 \"? 1/1 0.5\"\n", mspCluster10);
    }

    @Test
    public void testCompleteConversion() throws Exception {
        File tmpFile = File.createTempFile("ConversionTest", ".msp");

        converter.setOutputPath(tmpFile.getPath());

        List<IClusterSourceListener> listeners = new ArrayList<IClusterSourceListener>(1);
        listeners.add(converter);

        reader.readClustersIteratively(listeners);
        converter.close();
    }

    @Test
    public void testExternalFileConversion() throws Exception {
        File externalFile = new File("/tmp/ClusteringBin0300.clustering");

        if (!externalFile.exists()  || !externalFile.isFile())
            return;

        ClusteringFileReader externalReader = new ClusteringFileReader(externalFile);


        converter.setOutputPath(externalFile.getPath() + ".msp");

        List<IClusterSourceListener> listeners = new ArrayList<IClusterSourceListener>(1);
        listeners.add(converter);

        externalReader.readClustersIteratively(listeners);
        converter.close();
    }
}
