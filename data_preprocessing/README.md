## Data Pre-processing 


This document describes how to acquire the different datasets used and how (if at all) they were pre-processed before being used in our code.


### STCRDab 
The raw data files can be acquired as TSV files from http://opig.stats.ox.ac.uk/webapps/stcrdab/.

Run the script `get_cdrs_stcr.py` with the command `python get_cdrs_stcr.py <cdr a3/b3> <data directory> <IMGT structures directory>` to extract the CDRs.


### Twins 

The raw data files for the 6 pairs of twins can be acquired at http://labcfg.ibch.ru/tcr.html#MZTwins -> Clonesets (txt.tar.gz).

We chose to use the top 1000 sequences for analysis. The raw data files sort the CDR sequences by most frequent to least frequent read count. Therefore, we parse the first 1000 CDR sequences.

The script `parse_twins_data.py` parses the sequences into text files that only contain the CDR sequences.

Run the script with the command `python parse_twins_data.py <full path to raw data file> <output directory> <individual: TwA1/TwA2/TwC1/TwC2/TwD1/TwD2> <A/B> <number of sequences to parse>` to parse each raw data file into a file of just the CDR sequences ready for use with our code.

Example: 

```
python parse_twins_data.py datasets/twins_data/MZTwins_txt/a1alpha.txt datasets/twins_data/ TwA1 A 1000
```

### Glanville 
The raw data files can be acquired in the supplementary material (Table 1) of the paper. The excel file needs to be converted to CSV format to be read. No further pre-processing is needed.

### Dash

To acquire the raw data as a TSV file, go to https://vdjdb.cdr3.net/search and filter the data with PMID 28636592. Export the data as a TSV file with paired genes.

The script `parse_dash_data.py` parses the sequences into several text files which just contain the CDR3 sequences, one per alpha/beta epitope-specific repertoire, in the specified output directory. 

Run the script with the command `python parse_dash_data.py <full path to TSV data file> <output directory>`.

### McPAS

Acquire the data file at http://friedmanlab.weizmann.ac.il/McPAS-TCR/ -> download the complete database. This CSV file is the raw data file passed in to generate figure 6. We did not pre-process this data.