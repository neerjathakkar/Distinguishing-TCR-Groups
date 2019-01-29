__author__ = 'neerjathakkar'

import sys

# args
# 1: raw data txt file -- ...datasets/twins_data/MZTwins_txt/a1alpha.txt
# 2: output dir for parsed files -- ...datasets/twins_data/
# 3: individual: TwA1/TwA2/TwC1/TwC2/TwD1/TwD2 -- TwA1
# 4: A or B cdrs -- A
# 4: number of sequences to parse -- 1000

# to get data files:  http://labcfg.ibch.ru/tcr.html#MZTwins -> Clonesets (txt.tar.gz).

data_file = sys.argv[1]
output_dir = sys.argv[2]
if output_dir[-1] != '/': output_dir += '/'
individual = sys.argv[3]
cdr = sys.argv[4]
seqs = int(sys.argv[5])

with open(data_file) as f:
    lines = f.readlines()
    lines = lines[1:seqs+1]

cdr_i = 6
out_list = []
for line in lines:
    arr = line.split()
    out_list.append(arr[cdr_i])

out_filename = output_dir + individual + "_" + cdr + "_top_" + str(seqs) + ".txt"

outfile = open(out_filename, "w")
for seq in out_list:
    outfile.write(seq)
    outfile.write("\n")
