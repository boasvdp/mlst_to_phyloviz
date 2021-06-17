#!/usr/bin/env python3

import pandas as pd
import argparse
import os
import sys

def read_mlst_files(mlst_files, nr_loci):
  '''
  Read separate mlst files.

  Reads and checks MLST files produced using github.com/tseemann/mlst. Checks whether files are not empty and have expected number of loci.

  Parameters
  ----------
  mlst_files: list
    List of file names to read.
  nr_loci: int
    Number of loci expected in the scheme. Defaults to 7 through argparse.

  Returns
  -------
  df : pd.DataFrame
    Contains MLST profiles that were read and passed quality control.

  '''
  df_list = []
  nr_cols = nr_loci + 3
  for file in mlst_files:
    if os.stat(file).st_size > 0:
      tmp_df = pd.read_table(file, header=None)
      if tmp_df.shape == (1, nr_cols):
        df_list.append(tmp_df)
      else:
        print(str(file) + ' did not have expected dimensions. Skipping...')
    else:
      print(str(file) + ' is an empty file. Skipping...')
  nr_dfs = len(df_list)
  df = pd.concat(df_list)
  print('Combined ' + str(nr_dfs) + ' files')
  return df

def read_mlst_summary(mlst_summary):
  '''
  Read MLST summary file.

  Expects a checked summary file of tseemann/mlst output (can be produced by Unix `cat`).

  Parameters
  ----------
  mlst_summary : str
    Filename of MLST summary file.

  Returns
  -------
  df : pandas.DataFrame
    The parsed MLST summary file.

  Notes
  -----
  This functions employs no checks aside from standard Pandas reading.
  Ensure the MLST summary file is of good quality before reading using this function.
  Alternatively supply separate MLST files which will be checked.

  '''
  df = pd.read_csv(mlst_summary, sep = '\t', header=None)
  return df

def fix_header(dataframe):
  '''
  Add a meaningful header to MLST output DataFrame. This is parsed from the first 

  Paramaters
  ----------
  dataframe : pandas.DataFrame
    DataFrame read from file(s), with missing header

  Returns
  -------
  dataframe : pandas.DataFrame
    DataFrame with updated header
  '''
  header = ['file', 'scheme', 'ST']
  for allele in dataframe.iloc[0,3:]:
    gene_name = allele.split('(')[0]
    header.append(gene_name)
  dataframe.columns = header
  return dataframe

def select_scheme(dataframe, scheme):
  '''
  Select the correct MLST scheme.

  Scheme selection can be provided through argparse. Otherwise, the most common scheme is selected from the data.

  Parameters
  ----------
  dataframe : pandas.DataFrame
    MLST DataFrame.
  scheme : str
    Scheme to use and check.

  Returns
  -------
  selected_scheme : str
    Selected scheme.

  Notes
  -----
  Outputs a warning if the selected scheme is less than 75% of the dataset.

  '''
  valuecounts = dataframe['scheme'].value_counts()
  total = sum(valuecounts)
  if scheme == None:
    selected_scheme = valuecounts.index[0]
  else:
    selected_scheme = scheme
  selected_pct = (int(valuecounts[valuecounts.index == selected_scheme]) / total) * 100
  if selected_pct < 75:
    print_pct = "{:.2f}".format(selected_pct)
    print(selected_pct)
    print(type(selected_pct))
    print('WARNING: the selected scheme ' + str(selected_scheme) + ' applies to ' + str(print_pct) + '% of the dataset.')
  return selected_scheme

def clean_dataframe_scheme(dataframe, scheme):
  '''
  Clean DataFrame based on a selected MLST scheme.

  Parameters
  ----------
  dataframe : pandas.DataFrame
    DataFrame of MLST results with meaningful header.
  scheme : str
    Selected MLST scheme.

  Returns
  -------
  dataframe_clean : pandas.DataFrame
    DataFrame of samples typed according to selected scheme.
  dataframe_not_scheme : pandas.DataFrame
    DataFrame of samples NOT typed according to selected scheme.

  '''
  filter_scheme = dataframe['scheme'] == scheme
  dataframe_clean = dataframe[filter_scheme]
  dataframe_not_scheme = dataframe[~filter_scheme]
  dataframe_not_scheme = dataframe_not_scheme[dataframe_not_scheme['scheme'] != '-']
  return dataframe_clean, dataframe_not_scheme

def clean_dataframe_alleles(dataframe):
  '''
  Clean DataFrame based on failed typing.

  Parameters
  ----------
  dataframe : pandas.DataFrame
    DataFrame filtered on scheme.

  Returns
  -------
  dataframe_clean : pandas.DataFrame
    DataFrame with complete typing information.
  dataframe_failed_alleles : pandas.DataFrame
    DataFrame where typing could not identify all alleles.

  Notes
  -----
  Failed allele typing might indicate contamination or low quality of sequencing, but might also indicate novel alleles and MLSTs.
  Failed typings can be written to a file for submission to relevant MLST databases for novel assignments.
  '''
  dataframe_clean = dataframe.copy()
  for col in dataframe_clean.columns[3:]:
    stripped_series = dataframe_clean[col].str.replace('[a-zA-Z()]', '')
    dataframe_clean[col] =  pd.to_numeric(stripped_series, errors='coerce')
  filter_nan = dataframe_clean.isna().iloc[:,3:].any(axis=1)
  filter_no_ST = dataframe_clean['ST'] == '-'
  dataframe_clean = dataframe_clean[~filter_nan & ~filter_no_ST]
  dataframe_failed_alleles = dataframe[filter_nan | filter_no_ST]
  gene_dtypes = {}
  for gene in dataframe_clean.columns[3:]:
    gene_dtypes[gene] = 'int32'
  dataframe_clean = dataframe_clean.astype(gene_dtypes)
  dataframe_failed_alleles = dataframe_failed_alleles.astype(gene_dtypes, errors='ignore')
  return dataframe_clean, dataframe_failed_alleles

def main(args):
  if args.mlst_files != None:
    df = read_mlst_files(args.mlst_files, args.nr_loci)
  if args.mlst_summary != None:
    df = read_mlst_summary(args.mlst_summary)
  df = fix_header(df)
  scheme = select_scheme(df, args.scheme)
  df, df_not_scheme = clean_dataframe_scheme(df, scheme)
  df, df_failed_alleles = clean_dataframe_alleles(df)
  if args.include_filename == False:
    df = df.drop(['file', 'scheme'], axis=1)
  df.to_csv(args.out, sep = '\t', index=None)
  if args.not_scheme != None:
    df_not_scheme.to_csv(args.not_scheme, sep='\t', index=None)
  if args.failed_alleles != None:
    df_failed_alleles.to_csv(args.failed_alleles, sep='\t', index=None)


if __name__ == '__main__':
  # Parse arguments if executed as script
  parser = argparse.ArgumentParser(description='Convert MLST output data produced using github.com/tseemann/mlst for PhyloViz/eBurst analysis. To ensure correct parsing of MLST files, use separate MLST files as input (--mlst-files) and provide number of loci in scheme (--loci).')

  parser.add_argument('--mlst-summary', dest="mlst_summary", help="Summary file of mlst files, obtained through cat", type=str, metavar='MLST_OUTPUT_SUMMARY_FILE')
  parser.add_argument('--mlst-files', dest="mlst_files", nargs='+', help="Output files from tseemann/mlst, can be globbed", type=str, metavar='MLST_OUTPUT_FILES')
  parser.add_argument('--scheme', dest="scheme", help="Scheme to parse results for (default: select most common scheme)", type=str, metavar='MLST_SCHEME')
  parser.add_argument('--out', dest="out", help="Output data for PhyloViz/eBurst (default: allelic_profiles.tsv)", type=str, default='allelic_profiles.tsv', metavar='OUTPUT_FILE')
  parser.add_argument('--other-scheme-out', dest="not_scheme", help="Output file of samples which did not match selected scheme (default: don't output)", type=str, metavar='OUTPUT_FILE_FOR_MISMATCHING_SCHEME')
  parser.add_argument('--failed-out', dest="failed_alleles", help="Output file of samples which had failed allele calls or novel STs (default: don't output)", type=str, metavar='OUTPUT_FILE_FOR_MISSING_DATA')
  parser.add_argument('--loci', dest='nr_loci', help="Number of loci in MLST scheme. Used to parse MLST output files correctly (default: 7)", type=int, default=7, metavar='NUMBER_OF_LOCI')
  parser.add_argument('--include-filename', dest='include_filename', help='Include file name and scheme in output. NB: output file will not be directly compatible with PhyloViz anymore.', action='store_true')

  args = parser.parse_args()

  # exit if none or both of mlst_files and mlst_summary have been set
  if (args.mlst_files == None) == (args.mlst_summary == None):
    sys.exit('Need to set exactly one of either --mlst-dir or --mlst-summary')

  main(args)
