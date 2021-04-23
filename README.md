# clustering-file-converter

## Introduction

This Java tool is used to convert clustering results into a 
different file format. It's most important feature is that 
consensus spectra can be converted into the MSP file format. 
Thereby, a spectral library can quickly be created based on 
the clustering results.

## Installation

1) Download the latest ZIP file from the 
   [release section](https://github.com/spectra-cluster/clustering-file-converter/releases)
2) Extract the ZIP file
3) Launch the tool using `java -jar clustering-file-converter-[VERSION].jar` where 
   `[VERSION]` will vary depending on the release.

## Usage

The complete list of options is shown using the `-help` parameter.

The basic conversion to an MSP file looks like this:

```bash
java -jar clustering-file-converter-1.0-SNAPSHOT.jar \
  -format msp \
  -output_path my_library.msp \
  -min_ratio 0.7 \
  -min_size 3 \
  -spec_lib_add_annotation \
  my_cluster_result.clustering  
```

This will convert the `my_cluster_resul.clustering` clustering result
file into one MSP file.

Only clusters where 70% of the spectra were identified as the same
peptide (`-min_ration 0.7`) with at least 3 spetra (`-min_size 3`)
will be included in the MSP file.

The addition of annotations (`-spec_lib_add_annotation`) is optional
but is required by some downstream tools.
