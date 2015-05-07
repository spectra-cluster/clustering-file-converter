package uk.ac.ebi.pride.spectracluster.clusteringfileconverter.converters;

import uk.ac.ebi.pride.spectracluster.analysis.util.ClusterUtilities;
import uk.ac.ebi.pride.spectracluster.clusteringfilereader.objects.ICluster;

import java.io.BufferedWriter;
import java.io.FileWriter;

/**
 * Created by jg on 09.08.14.
 */
public class ClusterMgfConverter extends AbstractClusterConverter {
    public static String FILE_EXTENSION = "mgf";
    private int clusterCounter = 0;

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

        StringBuilder stringBuilder = new StringBuilder("BEGIN IONS\n");

        stringBuilder.append(String.format("TITLE=%s\n", (cluster.getId() != null) ? cluster.getId() : clusterUtilities.getMaxSequence(), clusterCounter));
        int charge = clusterUtilities.getCharge();
        stringBuilder.append(String.format("PEPMASS=%.3f\n", cluster.getAvPrecursorMz()));
        stringBuilder.append(String.format("CHARGE=%d%c\n", Math.abs(charge), (charge > 0 ? '+' : '-')));

        // add the peak list
        for (int i = 0; i < cluster.getConsensusMzValues().size(); i++) {
            // ignore peaks with 0 m/z and 0 intensity
            if (cluster.getConsensusMzValues().get(i) == 0) {
                System.out.println("Warning: Cluster " + cluster.getId() + " contains empty peak (m/z).");
                continue;
            }
            if (cluster.getConsensusIntensValues().get(i) == 0) {
                System.out.println("Warning: Cluster " + cluster.getId() + " contains empty peak (intensity).");
                continue;
            }

            stringBuilder.append(cluster.getConsensusMzValues().get(i)).append(" ").append(cluster.getConsensusIntensValues().get(i)).append("\n");
        }

        stringBuilder.append("END IONS\n\n");

        return stringBuilder.toString();
    }

    @Override
    public void close() throws Exception {
        if (writer != null) {
            writer.close();
            writer = null;
            if (!append)
                clusterCounter = 0;
        }
    }

    @Override
    public void onNewClusterRead(ICluster newCluster) {
        if (!shouldClusterBeExported(newCluster))
            return;

        writeStringToFile(convertCluster(newCluster), getFileHeader());
        clusterCounter++;
    }
}
