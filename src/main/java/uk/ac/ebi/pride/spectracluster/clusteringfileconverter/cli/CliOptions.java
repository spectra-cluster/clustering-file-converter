package uk.ac.ebi.pride.spectracluster.clusteringfileconverter.cli;

import org.apache.commons.cli.Option;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;

@SuppressWarnings("static-access")
public class CliOptions {

	public enum OPTIONS {
		HELP("help"),
        MAX_RATIO("max_ratio"),
        MIN_RATIO("min_ratio"),
        MIN_SIZE("min_size"),
        MAX_SIZE("max_size"),
        SPECIES("species"),
        COMBINE("combine"),
        FORMAT("format"),
        SPEC_LIB_ANNOTATION("spec_lib_add_annotation"),
        SPEC_LIB_NORMALIZE("spec_lib_normalize"),
        OUTPUT_PATH("output_path");

		private String value;

		OPTIONS(String value) {
			this.value = value;
		}
		
		public String getValue() {
			return value;
		}

		@Override
		public String toString() {
			return value;
		}
	}

	private static final Options options = new Options();

	static {
        Option minRatio = OptionBuilder
                .withDescription("limits the minimum RATIO a cluster may have to be included.")
                .withArgName("RATIO")
                .withType(Float.class)
                .hasArg()
                .create(OPTIONS.MIN_RATIO.getValue());
        options.addOption(minRatio);

        Option maxRatio = OptionBuilder
                .withDescription("limits the maximum RATIO a cluster may have to be included.")
                .hasArg()
                .withArgName("RATIO")
                .withType(Float.class)
                .create(OPTIONS.MAX_RATIO.getValue());
        options.addOption(maxRatio);

        Option minSize = OptionBuilder
                .withDescription("limits the minimum SIZE a cluster may have to be included.")
                .hasArg()
                .withArgName("SIZE")
                .withType(Integer.class)
                .create(OPTIONS.MIN_SIZE.getValue());
        options.addOption(minSize);

        Option maxSize = OptionBuilder
                .withDescription("limits the maximum SIZE a cluster may have to be included.")
                .hasArg()
                .withArgName("SIZE")
                .withType(Integer.class)
                .create(OPTIONS.MAX_SIZE.getValue());
        options.addOption(maxSize);

        Option species = OptionBuilder
                .withDescription("only exports cluster that contain at least on spectrum from the specified species.")
                .hasArg()
                .withArgName("TAXONOMY ID")
                .withType(Integer.class)
                .create(OPTIONS.SPECIES.getValue());
        options.addOption(species);

        Option format = OptionBuilder
                .withDescription("defines a file format of the output. Multiple formats may be specified simultaneously.")
                .hasArg()
                .withArgName("FORMAT")
                .create(OPTIONS.FORMAT.getValue());
        options.addOption(format);

        Option combine = OptionBuilder
                .withDescription("if set all passed .clustering files will be combined in a single output file.")
                .create(OPTIONS.COMBINE.getValue());
        options.addOption(combine);

        Option outputPath = OptionBuilder
                .withDescription("path to the output file. For each format the format specific extension will be appended to this path")
                .hasArg()
                .create(OPTIONS.OUTPUT_PATH.getValue());
        options.addOption(outputPath);

        Option specLibAnnotation = OptionBuilder
                .withDescription("if set peak annotations are added to the MSP file")
                .create(OPTIONS.SPEC_LIB_ANNOTATION.getValue());
        options.addOption(specLibAnnotation);

        Option specLibNormalize = OptionBuilder
                .withDescription("if set spectra are normalized according to the MSP specification (highest peak = 10,000), intensity reported as integer")
                .create(OPTIONS.SPEC_LIB_NORMALIZE.getValue());
        options.addOption(specLibNormalize);

		Option help = OptionBuilder
                .withDescription("print this help.")
                .create(OPTIONS.HELP.getValue());
		options.addOption(help);
	}

	public static Options getOptions() {
		return options;
	}
}
