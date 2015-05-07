package uk.ac.ebi.pride.spectracluster.clusteringfileconverter.util;

import uk.ac.ebi.pride.spectracluster.clusteringfilereader.objects.IModification;

import java.util.Comparator;

/**
 * Created by jg on 02.01.15.
 */
public class ModificationPositionComparator implements Comparator<IModification> {
    @Override
    public int compare(IModification o1, IModification o2) {
        return Integer.compare(o1.getPosition(), o2.getPosition());
    }
}
