# MLST to PhyloViz

Script to convert output from github.com/tseemann/mlst to PhyloViz/eBurst format.

## Why?

I found myself manually converting MLST summary tables to a format that works for PhyloViz. 
Additionally, I often had to do this by trial and error because I kept forgetting certain samples to remove, headers to add, etc.
This script performs a couple of standard preprocessing steps and checks a couple of things, before outputting a PhyloViz-ready table.

The script only supports classical MLST. cg/wgMLST is not supported and if the scheme uses any other number than seven genes, this must be made explicit on the command line.

## Processing steps

The script performs a couple of steps:

1. Read MLST output data from separate files (preferred) or a summary table obtained by running `cat mlst_out/* > mlst_summary.tsv`. 
  * Preparing a table using `cat` is faster, but some checks are not supported for that input.
  * When using separate files as input, the script checks whether the file is not empty, has the expected number of loci or a single row.
2. Create a sensible header based on the first sample.
3. Select the appropriate scheme. This can be a user-defined scheme (`--scheme` flag, corresponding to scheme names as used by tseemann/mlst), or the script can select the most common scheme in the dataset.
  * The script outputs a warning if less than 75% of the dataset is typed according to the user-defined scheme or most common scheme.
4. Samples typed according to the selected scheme are kept.
  * Samples typed according to another scheme can be written to a separate file to trace which samples might be contaminated (`--other-scheme-out`).
5. Samples without complete typing information are filtered. This includes samples with novel/incomplete alleles (see [the mlst github](https://github.com/tseemann/mlst#missing-data)) and samples with complete alleles, but novel sequence types.
  * The excluded samples can also be written to a separate file (`--failed-out`). This is where you'll find new STs to submit to databases.
6. Finally, the script writes the allelic profiles to a file.
  * File names and scheme can be included by running the script with `--include-filename`, but this renders the output PhyloViz-incompatible.

## Dependencies

- Python (https://www.python.org/)
- pandas (https://pandas.pydata.org/)

## Usage

```
usage: convert_mlst_to_phyloviz.py [-h] [--mlst-summary MLST_OUTPUT_SUMMARY_FILE] [--mlst-files MLST_OUTPUT_FILES [MLST_OUTPUT_FILES ...]]
                                   [--scheme MLST_SCHEME] [--out OUTPUT_FILE] [--other-scheme-out OUTPUT_FILE_FOR_MISMATCHING_SCHEME]
                                   [--failed-out OUTPUT_FILE_FOR_MISSING_DATA] [--loci NUMBER_OF_LOCI] [--include-filename]

Convert MLST output data produced using github.com/tseemann/mlst for PhyloViz/eBurst analysis. To ensure correct parsing of MLST files, use separate
MLST files as input (--mlst-files) and provide number of loci in scheme (--loci).

optional arguments:
  -h, --help            show this help message and exit
  --mlst-summary MLST_OUTPUT_SUMMARY_FILE
                        Summary file of mlst files, obtained through cat
  --mlst-files MLST_OUTPUT_FILES [MLST_OUTPUT_FILES ...]
                        Output files from tseemann/mlst, can be globbed
  --scheme MLST_SCHEME  Scheme to parse results for (default: select most common scheme)
  --out OUTPUT_FILE     Output data for PhyloViz/eBurst (default: allelic_profiles.tsv)
  --other-scheme-out OUTPUT_FILE_FOR_MISMATCHING_SCHEME
                        Output file of samples which did not match selected scheme (default: don't output)
  --failed-out OUTPUT_FILE_FOR_MISSING_DATA
                        Output file of samples which had failed allele calls or novel STs (default: don't output)
  --loci NUMBER_OF_LOCI
                        Number of loci in MLST scheme. Used to parse MLST output files correctly (default: 7)
  --include-filename    Include file name and scheme in output. NB: output file will not be directly compatible with PhyloViz anymore.
```

## Examples

#### Minimal

```
# Call STs using mlst
mkdir mlst_out

for file in assemblies/*.fasta
do
  NAME=$(basename $file .fasta)
  mlst "$file" > mlst_out/"$NAME".tsv
done

# Call mlst_to_phyloviz
convert_mlst_to_phyloviz.py --mlst-files mlst_out/*.tsv

head -n 5 allelic_profiles.tsv

>>ST      gapA    infB    mdh     pgi     phoE    rpoB    tonB
>>23      2       1       1       1       9       4       12
>>38      2       1       2       1       2       2       2
>>67      2       1       9       1       15      5       28
>>134     3       1       2       1       1       1       4
```

#### Advanced

```
# Use parallel to speed up for larger analyses
parallel --jobs 16 NAME'=$('basename {} .fasta');' mlst {} '>' mlst_out'/$'NAME.tsv ::: assemblies/*.fasta

# Also output additional files, force scheme and use custom output file name
convert_mlst_to_phyloviz.py --mlst-files mlst_out/*.tsv --scheme ssuis --other-scheme-out samples_other_scheme.tsv --failed-out samples_failed.tsv --out ssuis_allelic_profiles.tsv
```
