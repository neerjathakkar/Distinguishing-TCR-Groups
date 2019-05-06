# By Neerja Thakkar for "Balancing sensitivity and specificity in
# distinguishing TCR groups by CDR sequence similarity"
# See README for license information

from parse_sequences import txttoseqlist
from cluster_CDRs import cluster_progression, cluster_data
from cluster_distance_scores import score_repertoire, create_motif
import sys


# inputs (from command line): list of cluster files, directory, title
if len(sys.argv) < 5:
    print("need the following arguments: output directory, threshold to cluster to, "
          "human/mice, which epitope to get the motifs of, data directory")
    print("human epitopes: BMLF1, p65, M1")
    print("mice epitopes: PB1, PB1-F2, NP, PA, M38, M45, m139")
    sys.exit()

output_dir = sys.argv[1]
distance = float(sys.argv[2])
human = str(sys.argv[3])
epitope = sys.argv[4]
data_dir = sys.argv[5]          # datasets/dash_data/

epitope_filename = data_dir + epitope + "_beta.txt"

print(epitope_filename)
Xlabels = txttoseqlist(epitope_filename, trim_first=True, remove_special=True)


human_files = [data_dir + "BMLF1_beta.txt", data_dir + "p65_beta.txt",
               data_dir + "M1_beta.txt"]
mice_files = [data_dir + "PB1_beta.txt", data_dir + "PB1-F2_beta.txt",
              data_dir + "NP_beta.txt", data_dir + "PA_beta.txt",
              data_dir + "M38_beta.txt", data_dir + "M45_beta.txt",
              data_dir + "m139_beta.txt"]


other_files = []
if human.lower() == 'true':
    for file in human_files:
        if file != epitope_filename:
            other_files.append(file)
else:
    for file in mice_files:
        if file != epitope_filename:
            other_files.append(file)

print(other_files)

clusters2 = cluster_data(Xlabels, max_distance=0.2)
clusters3 = cluster_data(Xlabels, max_distance=0.3)
clusters4 = cluster_data(Xlabels, max_distance=0.4)

# convert clusters to list of clusters
clusters = []
for key in clusters2:
    clusters.append(clusters2[key])

other_clusters = []
for file in other_files:
    XlabelsOther = txttoseqlist(file, trim_first=True, remove_special=True)
    clusters2other = cluster_data(XlabelsOther, max_distance=0.3)
    clusters_other_list = []
    for key in clusters2other:
        clusters_other_list.append(clusters2other[key])
    other_clusters.append(clusters_other_list)

scores = score_repertoire(clusters,  other_clusters)
for score_data in scores:
    print(score_data)

log_file = open(output_dir + "output" + "_" + str(epitope) + ".txt", mode = 'w')
log_file.write("\nrepertoire: " + str(epitope_filename))

log_file.write("\nlowest scoring and therefore least specific to this repertoire:")
lowest = scores[0]
log_file.write(str(lowest[1]))
log_file.write("\nits pair cluster: ")
log_file.write(str(lowest[3]))
log_file.write("\nfrom repertoire: ")
log_file.write(str(other_files[lowest[2] - 1]))
log_file.write("\nthe score: ")
log_file.write(str(lowest[0]))
log_file.write("\ncreating motif...")
create_motif(lowest[1], "least_specific" + "_" + str(epitope), output_dir)
log_file.write("\ngetting cluster progression of motif...")
clustersprog = cluster_progression(output_dir + epitope + "_least_specific_", clusters2, clusters3, clusters4, clust = lowest[1])
log_file.write("\ncreating pair motif...")
create_motif(lowest[3], "least_specific_pair" + "_" + str(epitope), output_dir)

score_progression = score_repertoire(clustersprog, other_clusters)
for data in score_progression:
    log_file.write(str(data))

log_file.write("\nhighest scoring and therefore most specific to this repertoire:")
highest = scores[-1]
# score cannot be 1.0
if float(highest[0]) == 1.0:
    i = -2
    while float(highest[0]) == 1.0:
        highest = scores[i]
        i -=1

log_file.write(str(highest[1]))
log_file.write("\nits pair cluster:")
log_file.write(str(highest[3]))
log_file.write("\nfrom repertoire")
log_file.write(str(other_files[highest[2] - 1]))
log_file.write("\nthe score: ")
log_file.write(str(highest[0]))
log_file.write("\ncreating motif...")
create_motif(highest[1], "most_specific" + "_" + str(epitope), output_dir)
log_file.write("\ngetting cluster progression of motif...")
clustersprog = cluster_progression(output_dir + epitope + "_most_specific_", clusters2, clusters3, clusters4, clust = highest[1])
log_file.write("\ncreating pair motif...")
create_motif(highest[3], "most_specific_pair" + "_" + str(epitope), output_dir)

score_progression = score_repertoire(clustersprog, other_clusters)
for data in score_progression:
    log_file.write(str(data))

log_file.close()

