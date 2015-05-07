package uk.ac.ebi.pride.spectracluster.clusteringfileconverter.converters;

/**
 * Created by jg on 09.08.14.
 */
public class ConverterFactory {
    private ConverterFactory() {}

    public enum CONVERTER {
        MSP_CONVERTER("msp"),
        MGF_CONVERTER("mgf");

        private String name;
        private CONVERTER(String name) {
            this.name = name;
        }

        public String getName() {
            return name;
        }

        public static String[] getAllNames() {
            String[] allNames = new String[CONVERTER.values().length];

            for (int i = 0; i < CONVERTER.values().length; i++) {
                allNames[i] = CONVERTER.values()[i].getName();
            }

            return allNames;
        }

        public static CONVERTER getConverterByName(String name) {
            for (CONVERTER c : values()) {
                if (c.getName().equals(name))
                    return c;
            }

            return null;
        }
    }

    public static IClusterConverter getConverter(String converterName) throws Exception {
        CONVERTER converter = CONVERTER.getConverterByName(converterName);

        if (converter == null)
            throw new Exception("Unknown converter '" + converterName + "'");

        return getConverter(converter);
    }

    public static IClusterConverter getConverter(CONVERTER converter) {
        switch(converter) {
            case MSP_CONVERTER:
                return new ClusterMspConverter();
            case MGF_CONVERTER:
                return new ClusterMgfConverter();
            default:
                throw new IllegalStateException("Unsupported converter type passed.");
        }
    }
}
