__author__ = 'neerjathakkar'

import csv, sys

# args
# 1: raw data tsv file -- ...datasets/dash_data/dash_raw_data.tsv
# 2: output dir for parsed files -- ...datasets/dash_data/

# to get data file: https://vdjdb.cdr3.net/search, filter by PMID 28636592, export as TSV with paired gene

data_file = sys.argv[1]
output_dir = sys.argv[2]
if output_dir[-1] != '/': output_dir += '/'

dash_data = open(data_file, "rb")
dash_data_reader = csv.reader(dash_data, delimiter="\t")

gene_ind = 1
cdr_ind = 2
epitope_ind = 10

beta_epitopes = {}
alpha_epitopes = {}

i = 0

for row in dash_data_reader:
    if i == 0:
        i += 1
        continue

    cdr = row[cdr_ind]
    epitope = row[epitope_ind]

    if row[gene_ind] == "TRB":
        if epitope in beta_epitopes:
            beta_epitopes[epitope].append(cdr)
        else:
            beta_epitopes[epitope] = [cdr]
    else:
        if epitope in alpha_epitopes:
            alpha_epitopes[epitope].append(cdr)
        else:
            alpha_epitopes[epitope] = [cdr]

for key in alpha_epitopes:
    f = open(output_dir + key + "_alpha.txt", "w")
    for seq in alpha_epitopes[key]:
        f.write(seq)
        f.write("\n")
    f.close()

for key in beta_epitopes:
    f = open(output_dir + key + "_beta.txt", "w")
    for seq in beta_epitopes[key]:
        f.write(seq)
        f.write("\n")
    f.close()

