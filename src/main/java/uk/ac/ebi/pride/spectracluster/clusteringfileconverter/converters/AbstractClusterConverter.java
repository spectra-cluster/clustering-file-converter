package uk.ac.ebi.pride.spectracluster.clusteringfileconverter.converters;

import uk.ac.ebi.pride.spectracluster.clusteringfilereader.objects.ICluster;
import uk.ac.ebi.pride.spectracluster.clusteringfilereader.objects.IPeptideSpectrumMatch;
import uk.ac.ebi.pride.spectracluster.clusteringfilereader.objects.ISpectrumReference;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.Set;

/**
 * Created by jg on 10.08.14.
 */
public abstract class AbstractClusterConverter implements IClusterConverter {
    protected String outputPath;
    protected BufferedWriter writer;
    protected boolean append = false;

    protected int minSize = 0;
    protected int maxSize = Integer.MAX_VALUE;
    protected float minRatio = 0;
    protected float maxRatio = 1;
    protected Set<String> species = null;

    @Override
    public void setOutputPath(String outputPath) {
        this.outputPath = outputPath;
    }

    @Override
    public String getOuputPath() {
        return outputPath;
    }

    @Override
    public void setAppend(boolean append) {
        this.append = append;
    }

    @Override
    public abstract String getFileHeader();

    @Override
    public abstract String getFiletypeExtension();

    @Override
    public abstract String convertCluster(ICluster cluster);

    @Override
    public void close() throws Exception {
        if (writer != null) {
            writer.close();
            writer = null;
        }
    }

    @Override
    public void setMinSize(int minSize) {
        this.minSize = minSize;
    }

    @Override
    public void setMaxSize(int maxSize) {
        this.maxSize = maxSize;
    }

    @Override
    public void setMinRatio(float minRatio) {
        this.minRatio = minRatio;
    }

    @Override
    public void setMaxRatio(float maxRatio) {
        this.maxRatio = maxRatio;
    }

    @Override
    public int getMinSize() {
        return minSize;
    }

    @Override
    public int getMaxSize() {
        return maxSize;
    }

    @Override
    public float getMinRatio() {
        return minRatio;
    }

    @Override
    public float getMaxRatio() {
        return maxRatio;
    }

    @Override
    public abstract void onNewClusterRead(ICluster newCluster);

    /**
     * Checks whether the cluster should be exported based
     * on the set minSize, maxSize, minRatio, maxRatio, and taxonomyId.
     * @param cluster
     * @return
     */
    protected boolean shouldClusterBeExported(ICluster cluster) {
        if (cluster.getSpecCount() < minSize)
            return false;
        if (cluster.getSpecCount() > maxSize)
            return false;
        if (cluster.getMaxRatio() < minRatio)
            return false;
        if (cluster.getMaxRatio() > maxRatio)
            return false;

        // check if the taxonomy needs to be taken into consideration
        boolean containsSpecies = false;

        if (species != null && species.size() > 0) {
            for (ISpectrumReference specRef : cluster.getSpectrumReferences()) {
                // ignore if no species is available
                if (specRef.getSpecies() == null)
                    continue;

                String[] specRefSpecies = specRef.getSpecies().split(",");

                for (String s : specRefSpecies) {
                    if (species.contains(s)) {
                        containsSpecies = true;
                        break;
                    }
                }
            }
        }
        else {
            // disable the species test in case no species was set
            containsSpecies = true;
        }

        if (!containsSpecies)
            return false;

        return true;
    }

    /**
     * Writes the passed string to the defined output file. If
     * necessary, the BufferedWriter object is created by this
     * function.
     * @param string
     * @param fileHeader The header string written to a new file.
     */
    protected void writeStringToFile(String string, String fileHeader) {
        if (outputPath == null)
            throw new IllegalStateException("OutputPath must be set before clusters can be written.");

        try {
            if (writer == null) {
                File outputFile = new File(outputPath);
                boolean fileExists = outputFile.exists();

                writer = new BufferedWriter(new FileWriter(outputFile, append));

                // write the file header in case a new file was created or the current
                // file is being overwritten.
                if (!fileExists || !append)
                    writer.write(fileHeader);
            }

            writer.write(string);
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }

    @Override
    public Set<String> getSpecies() {
        return species;
    }

    public void setSpecies(Set<String> species) {
        this.species = species;
    }
}
