# By Neerja Thakkar for "Balancing sensitivity and specificity in
# distinguishing TCR groups by CDR sequence similarity"
# See README for license information

__author__ = 'neerjathakkar'

import sys
import csv

from sw_scoring import getDistanceSW
from parse_sequences import txttoseqlist

# inputs (from command line): list of cluster files, directory, title
if len(sys.argv) < 3:
    print("need the following arguments: output directory, a3 or b3, repertoire files")
    sys.exit()
else:
    output_directory = sys.argv[1]
    cdr_type = sys.argv[2]
    repertoire_file_list = sys.argv[3].split(',')


all_seqs = []

# read in the clusters
seq_to_rep = {}
rep_to_seq = {}
i = 0
for file in repertoire_file_list:
    i += 1
    content = txttoseqlist(file, trim_first=False)
    print(file)
    print(len(content))
    rep_to_seq[i] = content
    all_seqs += content
    for cdr in content:
        if cdr in seq_to_rep:
            found = False
            for nums in seq_to_rep[cdr]:
                if nums == i:
                    found = True
            if not found:
                seq_to_rep[cdr].append(i)
        else:
            seq_to_rep[cdr] = [i]

# closest distances
# [cdr1, cdr1 cluster, cdr2, cdr2 cluster, distance]
closest_distances = []

for curr_i in rep_to_seq:
    # extract cluster
    cluster = rep_to_seq[curr_i]

    # find smallest within-cluster distance
    for cdr1 in cluster:
        min_distance = 2.0
        min_distance_cdr = None
        for cdr2 in cluster:
            if cdr1 != cdr2:
                dist = getDistanceSW(cdr1, cdr2)
                if dist < min_distance:
                    min_distance = dist
                    min_distance_cdr = cdr2
        closest_distances.append([cdr1, curr_i, min_distance_cdr, curr_i, min_distance])


    # find smallest between-cluster distance
    for cdr1 in cluster:
        min_distance = 2.0
        min_distance_cdr = None
        min_distance_cluster = None
        for structural_cluster_i in rep_to_seq:
            if structural_cluster_i != curr_i:
                for cdr2 in rep_to_seq[structural_cluster_i]:
                    dist = getDistanceSW(cdr1, cdr2)
                    if dist < min_distance:
                        min_distance = dist
                        min_distance_cdr = cdr2
                        min_distance_cluster = structural_cluster_i
        closest_distances.append([cdr1, curr_i, min_distance_cdr, min_distance_cluster, min_distance])

file = output_directory + str(cdr_type) + "_stcr_pairs.csv"

with open(file, "w") as output_full:
    writer = csv.writer(output_full, lineterminator='\n')
    writer.writerows(closest_distances)


