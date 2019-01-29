# By Neerja Thakkar for "Balancing sensitivity and specificity in
# distinguishing TCR groups by CDR sequence similarity"
# See README for license information

import sys
from nearest_neighbor_threshold_analysis import run_threshold_analysis

# args
# 1: root of output directory
# 2: output file title
# 3: list of repertoire files (full path)

# inputs (from command line): list of cluster files, directory, title
if len(sys.argv) < 3:
    print("need the following arguments: output directory, file title, repertoire files")
else:
    output_directory = sys.argv[1]
    title = sys.argv[2]
    repertoire_file_list = sys.argv[3].split(',')

run_threshold_analysis(output_directory, title, repertoire_file_list)
