package uk.ac.ebi.pride.spectracluster.clusteringfileconverter.converters;

import uk.ac.ebi.pride.spectracluster.clusteringfileconverter.util.FastaFile;
import uk.ac.ebi.pride.spectracluster.clusteringfilereader.io.IClusterSourceListener;
import uk.ac.ebi.pride.spectracluster.clusteringfilereader.objects.ICluster;

import java.util.Set;

/**
 * Created by jg on 01.08.14.
 */
public interface IClusterConverter extends IClusterSourceListener {
    /**
     * Sets the path of the result file.
     * @param outputPath
     */
    public void setOutputPath(String outputPath);

    /**
     * Retrieves the current path of the result file.
     * @return
     */
    public String getOuputPath();

    /**
     * Defines whether the processed spectra should be appended
     * to the output file if it exists. Otherwise, the output file
     * will be overwritten.
     * @param append
     */
    public void setAppend(boolean append);

    /**
     * Retruns the header of the output file if required by the
     * format.
     * @return
     */
    public String getFileHeader();

    /**
     * Returns the default filetype extension (without the ".") of the
     * file format.
     * @return
     */
    public String getFiletypeExtension();

    /**
     * Convert a cluster to the corresponding representation of the file
     * format.
     * @param cluster
     * @return
     */
    public String convertCluster(ICluster cluster);

    /**
     * Close the handle to the output file. This function must be
     * called after writing is complete.
     * @throws Exception
     */
    public void close() throws Exception;

    public void setMinSize(int minSize);
    public void setMaxSize(int maxSize);
    public void setMinRatio(float minRatio);
    public void setMaxRatio(float maxRatio);
    public void setSpecies(Set<String> taxonomyIds);
    public void setFastaFile(FastaFile fastFile);


    public int getMinSize();
    public int getMaxSize();
    public float getMinRatio();
    public float getMaxRatio();
    public Set<String> getSpecies();
    public FastaFile getFastaFile();
}
